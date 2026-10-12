! ***********************************************************************
!
!   Copyright (C) 2025-2026  Henrik Knierim
!
!   This file is part of the element_sedimentation extension for MESA.
!   It is distributed under the MIT License; see the LICENSE file in the
!   root of the element_sedimentation repository for the full text.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
! ***********************************************************************

module sedimentation_phase_diagram

   ! Immiscibility phase diagram of a two-component subsystem (e.g. H-He).
   !
   ! The diagram is read from a CSV file with one header line followed by rows
   !
   !    y_lo, y_hi, T [K], P [Mbar]
   !
   ! where y_lo ("He-poor") and y_hi ("He-rich") are the number-fraction limits
   ! of the sedimenting species in the two-component subsystem. Rows are grouped
   ! in blocks of constant P; T is strictly increasing inside a block and P is
   ! strictly increasing between blocks. P is stored in dyn/cm^2.
   !
   ! phase_diagram_limits interpolates bilinearly: linear in T inside the two
   ! bracketing P blocks, then linear in P. It is pure (no module state, no I/O)
   ! and therefore safe to call from OpenMP loops.

   use const_def, only: dp
   use utils_lib, only: is_bad

   implicit none

   private
   public :: phase_diagram
   public :: load_phase_diagram, free_phase_diagram, check_phase_diagram, print_phase_diagram
   public :: phase_diagram_limits
   public :: phase_diagram_error_message

   ! error codes (ierr = 0 means OK, as everywhere in MESA)
   integer, parameter, public :: phase_diagram_ok = 0
   integer, parameter, public :: phase_diagram_err_not_loaded = 1    ! phase_diagram_limits on a diagram that is not loaded
   integer, parameter, public :: phase_diagram_err_bad_input = 2     ! NaN or Inf in P or T
   integer, parameter, public :: phase_diagram_err_out_of_range = 3  ! interpolated limit outside [0,1] (corrupt table)
   integer, parameter, public :: phase_diagram_err_file = 4          ! file cannot be opened, read or parsed
   integer, parameter, public :: phase_diagram_err_structure = 5     ! table violates the block structure or value ranges

   type phase_diagram
      character(len=1024) :: filename = ''
      logical :: loaded = .false.
      logical :: extrapolate_in_T = .false.
      real(dp) :: T_shift = 0d0
      integer :: n_points = 0
      real(dp), allocatable :: y_lo(:), y_hi(:), T(:), P(:)  ! (n_points); T includes T_shift, P in dyn/cm^2
      integer :: num_unique_P = 0
      integer, allocatable :: unique_P_index(:), num_T_at_unique_P(:)  ! (num_unique_P); first row, number of rows of each P block
      real(dp), allocatable :: unique_P(:)  ! (num_unique_P)
   end type phase_diagram

   real(dp), parameter :: Mbar = 1d12  ! dyn/cm^2
   integer, parameter :: num_columns = 4
   integer, parameter :: max_line_len = 1024
   integer, parameter :: min_points_per_block = 2  ! needed for linear interpolation in T
   integer, parameter :: min_num_blocks = 2        ! needed for linear interpolation in P
   character(len=*), parameter :: numeric_chars = '0123456789+-.eEdD'

contains

   subroutine load_phase_diagram(pd, filename, T_shift, extrapolate_in_T, ierr)
      ! read the CSV table, apply the temperature shift, build the P-block
      ! structure and validate it. On any error pd is left unloaded, ierr is
      ! nonzero and a message naming the file (and line where possible) is
      ! written to stdout.
      type(phase_diagram), intent(inout) :: pd
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: T_shift
      logical, intent(in) :: extrapolate_in_T
      integer, intent(out) :: ierr
      integer, allocatable :: line_no(:)
      integer :: iounit, ios

      ierr = 0
      call free_phase_diagram(pd)

      pd% filename = filename
      pd% T_shift = T_shift
      pd% extrapolate_in_T = extrapolate_in_T

      open(newunit=iounit, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(a)') 'load_phase_diagram: unable to open file ' // trim(filename)
         ierr = phase_diagram_err_file
         call free_phase_diagram(pd)
         return
      end if

      call read_csv_table(iounit, filename, pd% n_points, pd% y_lo, pd% y_hi, pd% T, pd% P, line_no, ierr)
      close(iounit)
      if (ierr /= 0) then
         call free_phase_diagram(pd)
         return
      end if

      pd% T(1:pd% n_points) = pd% T(1:pd% n_points) + T_shift
      pd% P(1:pd% n_points) = pd% P(1:pd% n_points)*Mbar

      call build_P_blocks(pd)

      call validate_structure(pd, ierr, line_no)
      if (ierr /= 0) then
         write(*,'(a)') 'load_phase_diagram: invalid phase diagram in file ' // trim(filename)
         call free_phase_diagram(pd)
         return
      end if

      pd% loaded = .true.

   end subroutine load_phase_diagram


   subroutine free_phase_diagram(pd)
      type(phase_diagram), intent(inout) :: pd
      if (allocated(pd% y_lo)) deallocate(pd% y_lo)
      if (allocated(pd% y_hi)) deallocate(pd% y_hi)
      if (allocated(pd% T)) deallocate(pd% T)
      if (allocated(pd% P)) deallocate(pd% P)
      if (allocated(pd% unique_P_index)) deallocate(pd% unique_P_index)
      if (allocated(pd% num_T_at_unique_P)) deallocate(pd% num_T_at_unique_P)
      if (allocated(pd% unique_P)) deallocate(pd% unique_P)
      pd% filename = ''
      pd% loaded = .false.
      pd% extrapolate_in_T = .false.
      pd% T_shift = 0d0
      pd% n_points = 0
      pd% num_unique_P = 0
   end subroutine free_phase_diagram


   subroutine check_phase_diagram(pd, ierr)
      ! structure checks of an already loaded diagram; messages use row numbers
      ! (data rows, header not counted) since file line numbers are not stored.
      type(phase_diagram), intent(in) :: pd
      integer, intent(out) :: ierr
      ierr = 0
      if (.not. pd% loaded) then
         write(*,'(a)') 'check_phase_diagram: phase diagram is not loaded'
         ierr = phase_diagram_err_not_loaded
         return
      end if
      call validate_structure(pd, ierr)
   end subroutine check_phase_diagram


   subroutine print_phase_diagram(pd, iounit)
      type(phase_diagram), intent(in) :: pd
      integer, intent(in) :: iounit
      integer :: i, j
      write(iounit,'(a)') '------------------------------------------------------------'
      write(iounit,'(a)') 'phase diagram: ' // trim(pd% filename)
      write(iounit,'(a,l2)') '   loaded = ', pd% loaded
      write(iounit,'(a,l2)') '   extrapolate_in_T = ', pd% extrapolate_in_T
      write(iounit,'(a,1pd26.16)') '   T_shift [K] = ', pd% T_shift
      write(iounit,'(a,i0)') '   n_points = ', pd% n_points
      write(iounit,'(a,i0)') '   num_unique_P = ', pd% num_unique_P
      if (.not. pd% loaded) then
         write(iounit,'(a)') '------------------------------------------------------------'
         return
      end if
      do j = 1, pd% num_unique_P
         write(iounit,'(a,i0,a,1pd26.16,a,i0,a)') '   block ', j, ': P [dyn/cm^2] = ', pd% unique_P(j), &
            ' with ', pd% num_T_at_unique_P(j), ' temperature points'
      end do
      write(iounit,'(a6,4a26)') 'row', 'P [dyn/cm^2]', 'T [K]', 'y_lo', 'y_hi'
      do i = 1, pd% n_points
         write(iounit,'(i6,4(1pd26.16))') i, pd% P(i), pd% T(i), pd% y_lo(i), pd% y_hi(i)
      end do
      write(iounit,'(a)') '------------------------------------------------------------'
   end subroutine print_phase_diagram


   pure subroutine phase_diagram_limits(pd, P, T, y_lo, y_hi, has_gap, ierr)
      ! Number-fraction limits of the immiscibility gap at (P, T) [cgs].
      ! Bilinear interpolation: linear in T inside the two bracketing P blocks,
      ! then linear in P between the blocks.
      !
      ! has_gap = .false. (and y_lo = y_hi = -1d0) when P is outside the P range
      ! of the table, T is below the lowest T of either bracketing block, T is
      ! above the highest T of either bracketing block and extrapolate_in_T is
      ! off, or the interpolated limits satisfy y_lo >= y_hi (closed gap).
      ! With extrapolate_in_T, T above a block's range is extrapolated linearly
      ! from the block's last two points and both limits are clamped to [0,1].
      !
      ! ierr /= 0 only for an unloaded diagram, NaN/Inf input, or an interpolated
      ! limit outside [0,1] without extrapolation (impossible for a table that
      ! passed check_phase_diagram). No I/O, no module state: thread safe.
      type(phase_diagram), intent(in) :: pd
      real(dp), intent(in) :: P, T
      real(dp), intent(out) :: y_lo, y_hi
      logical, intent(out) :: has_gap
      integer, intent(out) :: ierr
      integer :: j, j0, j1
      real(dp) :: fP, y_lo_block(2), y_hi_block(2)
      logical :: ok

      ierr = 0
      has_gap = .false.
      y_lo = -1d0
      y_hi = -1d0

      if (.not. pd% loaded) then
         ierr = phase_diagram_err_not_loaded
         return
      end if
      if (is_bad(P) .or. is_bad(T)) then
         ierr = phase_diagram_err_bad_input
         return
      end if

      if (P < pd% unique_P(1) .or. P > pd% unique_P(pd% num_unique_P)) return

      ! bracketing P blocks j0 <= j1 = j0 + 1 with unique_P(j0) <= P <= unique_P(j1)
      j0 = locate_interval(P, pd% num_unique_P, pd% unique_P)
      j1 = j0 + 1

      do j = 1, 2
         call limits_in_block(pd, j0 + j - 1, T, y_lo_block(j), y_hi_block(j), ok)
         if (.not. ok) return  ! T outside the block's range: no gap
      end do

      fP = (P - pd% unique_P(j0))/(pd% unique_P(j1) - pd% unique_P(j0))
      y_lo = (1d0 - fP)*y_lo_block(1) + fP*y_lo_block(2)
      y_hi = (1d0 - fP)*y_hi_block(1) + fP*y_hi_block(2)

      if (pd% extrapolate_in_T) then
         ! extrapolated values may leave [0,1]
         y_lo = min(1d0, max(0d0, y_lo))
         y_hi = min(1d0, max(0d0, y_hi))
      else
         ! interpolation of a valid table cannot leave [0,1]; guard against corrupt data
         if (y_lo < 0d0 .or. y_lo > 1d0 .or. y_hi < 0d0 .or. y_hi > 1d0) then
            ierr = phase_diagram_err_out_of_range
            y_lo = -1d0
            y_hi = -1d0
            return
         end if
      end if

      if (y_lo >= y_hi) then  ! closed gap
         y_lo = -1d0
         y_hi = -1d0
         return
      end if

      has_gap = .true.

   end subroutine phase_diagram_limits


   function phase_diagram_error_message(ierr) result(msg)
      integer, intent(in) :: ierr
      character(len=64) :: msg
      select case (ierr)
      case (phase_diagram_ok)
         msg = 'ok'
      case (phase_diagram_err_not_loaded)
         msg = 'phase diagram is not loaded'
      case (phase_diagram_err_bad_input)
         msg = 'NaN or Inf in P or T'
      case (phase_diagram_err_out_of_range)
         msg = 'interpolated limit outside [0,1]'
      case (phase_diagram_err_file)
         msg = 'cannot open, read or parse the phase diagram file'
      case (phase_diagram_err_structure)
         msg = 'phase diagram table has an invalid structure'
      case default
         msg = 'unknown phase diagram error code'
      end select
   end function phase_diagram_error_message


   ! ----------------------------------------------------------------------
   ! private helpers
   ! ----------------------------------------------------------------------

   pure subroutine limits_in_block(pd, j, T, y_lo, y_hi, ok)
      ! linear interpolation in T inside P block j. ok = .false. if T is below
      ! the block's range, or above it without extrapolate_in_T.
      type(phase_diagram), intent(in) :: pd
      integer, intent(in) :: j
      real(dp), intent(in) :: T
      real(dp), intent(out) :: y_lo, y_hi
      logical, intent(out) :: ok
      integer :: i0, i1, i, n
      real(dp) :: fT
      logical :: extrapolating

      y_lo = -1d0
      y_hi = -1d0
      ok = .false.

      i0 = pd% unique_P_index(j)
      n = pd% num_T_at_unique_P(j)
      i1 = i0 + n - 1

      if (T < pd% T(i0)) return
      extrapolating = (T > pd% T(i1))
      if (extrapolating .and. .not. pd% extrapolate_in_T) return

      if (extrapolating) then
         i = i1 - 1  ! last two points
      else
         i = i0 - 1 + locate_interval(T, n, pd% T(i0:i1))
      end if

      fT = (T - pd% T(i))/(pd% T(i+1) - pd% T(i))
      y_lo = (1d0 - fT)*pd% y_lo(i) + fT*pd% y_lo(i+1)
      y_hi = (1d0 - fT)*pd% y_hi(i) + fT*pd% y_hi(i+1)

      if (.not. extrapolating) then
         ! keep round-off from pushing an interpolated value past its end points
         y_lo = min(max(y_lo, min(pd% y_lo(i), pd% y_lo(i+1))), max(pd% y_lo(i), pd% y_lo(i+1)))
         y_hi = min(max(y_hi, min(pd% y_hi(i), pd% y_hi(i+1))), max(pd% y_hi(i), pd% y_hi(i+1)))
      end if

      ok = .true.

   end subroutine limits_in_block


   pure integer function locate_interval(x, n, xs) result(i)
      ! bisection: largest i in [1, n-1] with xs(i) <= x, for strictly increasing
      ! xs with xs(1) <= x <= xs(n). Returns n-1 for x == xs(n).
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp), intent(in) :: xs(n)
      integer :: lo, hi, mid
      lo = 1
      hi = n
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (xs(mid) <= x) then
            lo = mid
         else
            hi = mid
         end if
      end do
      i = lo
   end function locate_interval


   subroutine build_P_blocks(pd)
      ! group consecutive rows of equal P into blocks
      type(phase_diagram), intent(inout) :: pd
      integer :: i, nb
      if (allocated(pd% unique_P_index)) deallocate(pd% unique_P_index)
      if (allocated(pd% num_T_at_unique_P)) deallocate(pd% num_T_at_unique_P)
      if (allocated(pd% unique_P)) deallocate(pd% unique_P)
      pd% num_unique_P = 0
      if (pd% n_points <= 0) return
      nb = 1
      do i = 2, pd% n_points
         if (pd% P(i) /= pd% P(i-1)) nb = nb + 1
      end do
      allocate(pd% unique_P_index(nb), pd% num_T_at_unique_P(nb), pd% unique_P(nb))
      pd% num_unique_P = nb
      nb = 1
      pd% unique_P_index(1) = 1
      pd% unique_P(1) = pd% P(1)
      do i = 2, pd% n_points
         if (pd% P(i) /= pd% P(i-1)) then
            nb = nb + 1
            pd% unique_P_index(nb) = i
            pd% unique_P(nb) = pd% P(i)
         end if
      end do
      do i = 1, nb - 1
         pd% num_T_at_unique_P(i) = pd% unique_P_index(i+1) - pd% unique_P_index(i)
      end do
      pd% num_T_at_unique_P(nb) = pd% n_points - pd% unique_P_index(nb) + 1
   end subroutine build_P_blocks


   subroutine validate_structure(pd, ierr, line_no)
      ! all structure and value checks. line_no (optional, size n_points) gives
      ! the file line of each row for the messages.
      type(phase_diagram), intent(in) :: pd
      integer, intent(out) :: ierr
      integer, intent(in), optional :: line_no(:)
      integer :: i, j, i0, i1, n
      character(len=64) :: where_str

      ierr = phase_diagram_err_structure

      if (pd% n_points < 1) then
         write(*,'(a)') 'check_phase_diagram: the table contains no data rows'
         return
      end if
      if (.not. (allocated(pd% y_lo) .and. allocated(pd% y_hi) .and. allocated(pd% T) .and. allocated(pd% P))) then
         write(*,'(a)') 'check_phase_diagram: data arrays are not allocated'
         return
      end if
      if (size(pd% y_lo) < pd% n_points .or. size(pd% y_hi) < pd% n_points &
            .or. size(pd% T) < pd% n_points .or. size(pd% P) < pd% n_points) then
         write(*,'(a)') 'check_phase_diagram: data arrays are smaller than n_points'
         return
      end if
      if (.not. (allocated(pd% unique_P_index) .and. allocated(pd% num_T_at_unique_P) .and. allocated(pd% unique_P))) then
         write(*,'(a)') 'check_phase_diagram: pressure block arrays are not allocated'
         return
      end if
      if (pd% num_unique_P < 1 .or. size(pd% unique_P_index) < pd% num_unique_P &
            .or. size(pd% num_T_at_unique_P) < pd% num_unique_P .or. size(pd% unique_P) < pd% num_unique_P) then
         write(*,'(a)') 'check_phase_diagram: pressure block arrays are inconsistent with num_unique_P'
         return
      end if
      if (sum(pd% num_T_at_unique_P(1:pd% num_unique_P)) /= pd% n_points) then
         write(*,'(a)') 'check_phase_diagram: pressure blocks do not cover all rows'
         return
      end if

      ! row-wise value checks
      do i = 1, pd% n_points
         where_str = row_label(i, line_no)
         if (is_bad(pd% y_lo(i)) .or. is_bad(pd% y_hi(i)) .or. is_bad(pd% T(i)) .or. is_bad(pd% P(i))) then
            write(*,'(a)') 'check_phase_diagram: NaN or Inf in ' // trim(where_str)
            return
         end if
         if (pd% y_lo(i) < 0d0 .or. pd% y_lo(i) > 1d0) then
            write(*,'(a,1pd26.16)') 'check_phase_diagram: y_lo outside [0,1] in ' // trim(where_str) // ': y_lo = ', pd% y_lo(i)
            return
         end if
         if (pd% y_hi(i) < 0d0 .or. pd% y_hi(i) > 1d0) then
            write(*,'(a,1pd26.16)') 'check_phase_diagram: y_hi outside [0,1] in ' // trim(where_str) // ': y_hi = ', pd% y_hi(i)
            return
         end if
         if (pd% y_lo(i) >= pd% y_hi(i)) then
            write(*,'(a,2(1pd26.16))') 'check_phase_diagram: y_lo >= y_hi in ' // trim(where_str) // ': y_lo, y_hi = ', &
               pd% y_lo(i), pd% y_hi(i)
            return
         end if
         if (pd% P(i) <= 0d0) then
            write(*,'(a,1pd26.16)') 'check_phase_diagram: P <= 0 in ' // trim(where_str) // ': P [dyn/cm^2] = ', pd% P(i)
            return
         end if
         if (pd% T(i) <= 0d0) then
            write(*,'(a,2(1pd26.16))') 'check_phase_diagram: T <= 0 in ' // trim(where_str) // &
               ' (T includes T_shift): T [K], T_shift [K] = ', pd% T(i), pd% T_shift
            return
         end if
      end do

      ! block structure
      if (pd% num_unique_P < min_num_blocks) then
         write(*,'(a,i0,a,i0)') 'check_phase_diagram: need at least ', min_num_blocks, &
            ' blocks of constant P for interpolation in P, found ', pd% num_unique_P
         return
      end if
      do j = 1, pd% num_unique_P
         i0 = pd% unique_P_index(j)
         n = pd% num_T_at_unique_P(j)
         i1 = i0 + n - 1
         if (i0 < 1 .or. i1 > pd% n_points .or. n < 1) then
            write(*,'(a,i0)') 'check_phase_diagram: inconsistent row range for P block ', j
            return
         end if
         if (pd% unique_P(j) /= pd% P(i0)) then
            write(*,'(a,i0)') 'check_phase_diagram: unique_P does not match the first row of P block ', j
            return
         end if
         if (n < min_points_per_block) then
            write(*,'(a,i0,a,1pd26.16,a,i0,a)') 'check_phase_diagram: P block ', j, ' (P [dyn/cm^2] = ', pd% unique_P(j), &
               ') starting in ' // trim(row_label(i0, line_no)) // ' has only ', n, ' temperature point(s); need at least 2'
            return
         end if
         do i = i0, i1
            if (pd% P(i) /= pd% unique_P(j)) then
               write(*,'(a,i0)') 'check_phase_diagram: P is not constant inside P block ', j
               return
            end if
            if (i > i0) then
               if (pd% T(i) <= pd% T(i-1)) then
                  write(*,'(a,2(1pd26.16))') 'check_phase_diagram: T is not strictly increasing within a P block at ' // &
                     trim(row_label(i, line_no)) // ': T(previous), T = ', pd% T(i-1), pd% T(i)
                  return
               end if
            end if
         end do
         if (j > 1) then
            if (pd% unique_P(j) <= pd% unique_P(j-1)) then
               write(*,'(a,2(1pd26.16))') 'check_phase_diagram: P is not strictly increasing between P blocks at ' // &
                  trim(row_label(i0, line_no)) // ': P(previous block), P [dyn/cm^2] = ', pd% unique_P(j-1), pd% unique_P(j)
               return
            end if
         end if
      end do

      ierr = 0

   end subroutine validate_structure


   function row_label(i, line_no) result(str)
      ! 'data row i' or 'data row i (file line L)'
      integer, intent(in) :: i
      integer, intent(in), optional :: line_no(:)
      character(len=64) :: str
      if (present(line_no)) then
         if (i >= 1 .and. i <= size(line_no)) then
            write(str,'(a,i0,a,i0,a)') 'data row ', i, ' (file line ', line_no(i), ')'
            return
         end if
      end if
      write(str,'(a,i0)') 'data row ', i
   end function row_label


   subroutine read_csv_table(iounit, filename, n, y_lo, y_hi, T, P, line_no, ierr)
      ! read all data rows of the CSV file on iounit. Blank lines are skipped,
      ! the first non-blank line is the header. Each data row must contain
      ! exactly num_columns comma separated numbers. P is returned as read
      ! (Mbar), T without shift. line_no(i) is the file line of data row i.
      integer, intent(in) :: iounit
      character(len=*), intent(in) :: filename
      integer, intent(out) :: n
      real(dp), allocatable, intent(out) :: y_lo(:), y_hi(:), T(:), P(:)
      integer, allocatable, intent(out) :: line_no(:)
      integer, intent(out) :: ierr
      character(len=max_line_len) :: line
      character(len=max_line_len) :: fields(num_columns)
      real(dp) :: vals(num_columns)
      real(dp), allocatable :: buf(:,:)
      integer, allocatable :: buf_lines(:)
      integer :: ios, nline, nfields, c, capacity
      logical :: header_seen

      ierr = 0
      n = 0
      nline = 0
      header_seen = .false.
      capacity = 64
      allocate(buf(num_columns, capacity), buf_lines(capacity))

      do
         read(iounit, '(a)', iostat=ios) line
         if (ios < 0) exit  ! end of file
         nline = nline + 1
         if (ios > 0) then
            write(*,'(a,i0,a)') 'load_phase_diagram: read error at line ', nline, ' of file ' // trim(filename)
            ierr = phase_diagram_err_file
            return
         end if
         call clean_line(line)
         if (len_trim(line) == 0) cycle

         call split_fields(line, nfields, fields)

         if (.not. header_seen) then
            header_seen = .true.
            ! a first line consisting of numbers only means the header is missing
            if (nfields == num_columns) then
               if (all_numeric(fields)) then
                  write(*,'(a,i0,a)') 'load_phase_diagram: line ', nline, ' of file ' // trim(filename) // &
                     ' contains numbers but a header line is required (e.g. "He-poor","He-rich","T","P")'
                  ierr = phase_diagram_err_file
                  return
               end if
            end if
            cycle
         end if

         if (nfields /= num_columns) then
            write(*,'(a,i0,a,i0,a,i0,a)') 'load_phase_diagram: line ', nline, ' of file ' // trim(filename) // &
               ': expected ', num_columns, ' comma separated values (y_lo, y_hi, T, P), found ', nfields, &
               ': ' // trim(line)
            ierr = phase_diagram_err_file
            return
         end if
         do c = 1, num_columns
            if (.not. parse_real(fields(c), vals(c))) then
               write(*,'(a,i0,a,i0,a)') 'load_phase_diagram: line ', nline, ' of file ' // trim(filename) // &
                  ': cannot parse column ', c, ' as a number: "' // trim(fields(c)) // '"'
               ierr = phase_diagram_err_file
               return
            end if
         end do

         if (n == capacity) call grow_buffers(buf, buf_lines, capacity)
         n = n + 1
         buf(:, n) = vals
         buf_lines(n) = nline
      end do

      if (.not. header_seen) then
         write(*,'(a)') 'load_phase_diagram: file ' // trim(filename) // ' is empty'
         ierr = phase_diagram_err_file
         return
      end if
      if (n == 0) then
         write(*,'(a)') 'load_phase_diagram: file ' // trim(filename) // ' contains a header but no data rows'
         ierr = phase_diagram_err_file
         return
      end if

      allocate(y_lo(n), y_hi(n), T(n), P(n), line_no(n))
      y_lo = buf(1, 1:n)
      y_hi = buf(2, 1:n)
      T = buf(3, 1:n)
      P = buf(4, 1:n)
      line_no = buf_lines(1:n)

   end subroutine read_csv_table


   subroutine grow_buffers(buf, buf_lines, capacity)
      real(dp), allocatable, intent(inout) :: buf(:,:)
      integer, allocatable, intent(inout) :: buf_lines(:)
      integer, intent(inout) :: capacity
      real(dp), allocatable :: new_buf(:,:)
      integer, allocatable :: new_lines(:)
      integer :: new_capacity
      new_capacity = 2*capacity
      allocate(new_buf(num_columns, new_capacity), new_lines(new_capacity))
      new_buf(:, 1:capacity) = buf(:, 1:capacity)
      new_lines(1:capacity) = buf_lines(1:capacity)
      call move_alloc(new_buf, buf)
      call move_alloc(new_lines, buf_lines)
      capacity = new_capacity
   end subroutine grow_buffers


   subroutine clean_line(line)
      ! replace carriage returns and tabs by blanks (CRLF files, tab separated whitespace)
      character(len=*), intent(inout) :: line
      integer :: i
      do i = 1, len(line)
         if (line(i:i) == achar(13) .or. line(i:i) == achar(9)) line(i:i) = ' '
      end do
   end subroutine clean_line


   subroutine split_fields(line, nfields, fields)
      ! split a line at commas; returns at most num_columns fields (blank padded,
      ! leading blanks removed) and the true number of fields in the line.
      character(len=*), intent(in) :: line
      integer, intent(out) :: nfields
      character(len=*), intent(out) :: fields(:)
      integer :: start, pos, l
      fields = ''
      nfields = 0
      l = len_trim(line)
      start = 1
      do
         nfields = nfields + 1
         pos = index(line(start:l), ',')
         if (pos == 0) then
            if (nfields <= size(fields)) fields(nfields) = adjustl(line(start:l))
            exit
         end if
         if (nfields <= size(fields)) fields(nfields) = adjustl(line(start:start+pos-2))
         start = start + pos
         if (start > l) then  ! trailing comma: one more (empty) field
            nfields = nfields + 1
            exit
         end if
      end do
   end subroutine split_fields


   logical function parse_real(field, val)
      ! strict conversion of a single numeric token
      character(len=*), intent(in) :: field
      real(dp), intent(out) :: val
      integer :: ios
      val = 0d0
      parse_real = .false.
      if (len_trim(field) == 0) return
      if (verify(trim(field), numeric_chars) /= 0) return
      read(field, *, iostat=ios) val
      if (ios /= 0) return
      if (is_bad(val)) return
      parse_real = .true.
   end function parse_real


   logical function all_numeric(fields)
      character(len=*), intent(in) :: fields(:)
      real(dp) :: dummy
      integer :: c
      all_numeric = .true.
      do c = 1, size(fields)
         if (.not. parse_real(fields(c), dummy)) then
            all_numeric = .false.
            return
         end if
      end do
   end function all_numeric

end module sedimentation_phase_diagram
