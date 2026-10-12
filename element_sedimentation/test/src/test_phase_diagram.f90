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

module test_phase_diagram

   ! unit tests for sedimentation_phase_diagram. Run from the test directory
   ! (the data file is addressed relative to it).

   use const_def, only: dp
   use sedimentation_phase_diagram
   use test_support, only: assert_true, assert_close

   implicit none

   private
   public :: run_phase_diagram_tests

   character(len=*), parameter :: data_file = 'data/Schoettler+2018.csv'
   character(len=*), parameter :: tmp_file = 'test_phase_diagram_tmp_table.csv'
   real(dp), parameter :: Mbar = 1d12
   real(dp), parameter :: rtol = 1d-12, atol = 1d-14
   real(dp), parameter :: tight_rtol = 1d-14, tight_atol = 1d-15

   ! expected structure of Schoettler+2018.csv
   integer, parameter :: expected_n_points = 61
   integer, parameter :: expected_num_blocks = 8
   real(dp), parameter :: expected_P_Mbar(expected_num_blocks) = &
      [0.5d0, 1d0, 1.2d0, 1.5d0, 2d0, 4d0, 10d0, 24d0]
   integer, parameter :: expected_num_T(expected_num_blocks) = [2, 3, 6, 8, 9, 11, 11, 11]

contains

   subroutine run_phase_diagram_tests(nfail)
      integer, intent(inout) :: nfail
      type(phase_diagram) :: pd
      integer :: ierr

      write(*,'(a)') 'test_phase_diagram'

      call load_phase_diagram(pd, data_file, 0d0, .false., ierr)
      call assert_true('load ' // data_file, ierr == 0 .and. pd% loaded, nfail)
      if (ierr /= 0) then
         write(*,'(a)') '   cannot continue: the tester must run from the test directory'
         return
      end if

      call test_structure(pd, nfail)
      call test_nodes(pd, nfail)
      call test_midpoints(pd, nfail)
      call test_out_of_range(pd, nfail)
      call test_extrapolation(pd, nfail)
      call test_temperature_shift(pd, nfail)
      call test_check_and_print(pd, nfail)
      call test_reload(nfail)
      call test_openmp(pd, nfail)
      call test_tolerant_formats(nfail)
      call test_malformed_files(nfail)

      call free_phase_diagram(pd)
      call assert_true('free_phase_diagram resets the state', &
         .not. pd% loaded .and. pd% n_points == 0 .and. pd% num_unique_P == 0 &
         .and. .not. allocated(pd% T) .and. .not. allocated(pd% unique_P), nfail)

   end subroutine run_phase_diagram_tests


   subroutine test_structure(pd, nfail)
      type(phase_diagram), intent(in) :: pd
      integer, intent(inout) :: nfail
      integer :: j, expected_index(expected_num_blocks)
      real(dp) :: max_err

      call assert_true('n_points', pd% n_points == expected_n_points, nfail)
      call assert_true('num_unique_P', pd% num_unique_P == expected_num_blocks, nfail)
      if (pd% num_unique_P /= expected_num_blocks) return

      max_err = 0d0
      do j = 1, expected_num_blocks
         max_err = max(max_err, abs(pd% unique_P(j)/Mbar - expected_P_Mbar(j)))
      end do
      call assert_close('unique_P in Mbar', max_err, 0d0, 0d0, 1d-13, nfail)
      call assert_true('num_T_at_unique_P', all(pd% num_T_at_unique_P == expected_num_T), nfail)

      expected_index(1) = 1
      do j = 2, expected_num_blocks
         expected_index(j) = expected_index(j-1) + expected_num_T(j-1)
      end do
      call assert_true('unique_P_index', all(pd% unique_P_index == expected_index), nfail)
      call assert_true('P stored in cgs', pd% P(1) == 0.5d0*Mbar, nfail)
      call assert_true('T_shift stored', pd% T_shift == 0d0, nfail)
      call assert_true('extrapolate_in_T stored', .not. pd% extrapolate_in_T, nfail)
      call assert_true('filename stored', trim(pd% filename) == data_file, nfail)
      call assert_true('all rows have y_lo < y_hi', all(pd% y_lo < pd% y_hi), nfail)
   end subroutine test_structure


   subroutine test_nodes(pd, nfail)
      ! limits at the table nodes reproduce the table
      type(phase_diagram), intent(in) :: pd
      integer, intent(inout) :: nfail
      integer :: i, ierr, num_bad_ierr, num_no_gap
      real(dp) :: y_lo, y_hi, max_err_lo, max_err_hi
      logical :: has_gap

      num_bad_ierr = 0
      num_no_gap = 0
      max_err_lo = 0d0
      max_err_hi = 0d0
      do i = 1, pd% n_points
         call phase_diagram_limits(pd, pd% P(i), pd% T(i), y_lo, y_hi, has_gap, ierr)
         if (ierr /= 0) num_bad_ierr = num_bad_ierr + 1
         if (.not. has_gap) num_no_gap = num_no_gap + 1
         max_err_lo = max(max_err_lo, abs(y_lo - pd% y_lo(i)))
         max_err_hi = max(max_err_hi, abs(y_hi - pd% y_hi(i)))
      end do
      call assert_true('nodes: no error', num_bad_ierr == 0, nfail)
      call assert_true('nodes: all have a gap', num_no_gap == 0, nfail)
      call assert_close('nodes: y_lo reproduces the table', max_err_lo, 0d0, 0d0, tight_atol, nfail)
      call assert_close('nodes: y_hi reproduces the table', max_err_hi, 0d0, 0d0, tight_atol, nfail)
   end subroutine test_nodes


   subroutine test_midpoints(pd, nfail)
      type(phase_diagram), intent(in) :: pd
      integer, intent(inout) :: nfail
      integer :: i1, i2, i3, i4, ierr
      real(dp) :: y_lo, y_hi, e_lo, e_hi, fT, fP
      logical :: has_gap

      ! midpoint in T inside the 1.5 Mbar block
      i1 = find_row(pd, 1.5d0, 3000d0)
      i2 = find_row(pd, 1.5d0, 4000d0)
      call assert_true('midpoint T: rows found', i1 > 0 .and. i2 == i1 + 1, nfail)
      call phase_diagram_limits(pd, 1.5d0*Mbar, 3500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('midpoint T: ok and gap', ierr == 0 .and. has_gap, nfail)
      call assert_close('midpoint T: y_lo', y_lo, 0.5d0*(pd% y_lo(i1) + pd% y_lo(i2)), rtol, atol, nfail)
      call assert_close('midpoint T: y_hi', y_hi, 0.5d0*(pd% y_hi(i1) + pd% y_hi(i2)), rtol, atol, nfail)

      ! midpoint in P between the 1.5 and 2 Mbar blocks at a common T node
      i1 = find_row(pd, 1.5d0, 5000d0)
      i2 = find_row(pd, 2d0, 5000d0)
      call assert_true('midpoint P: rows found', i1 > 0 .and. i2 > 0, nfail)
      call phase_diagram_limits(pd, 1.75d0*Mbar, 5000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('midpoint P: ok and gap', ierr == 0 .and. has_gap, nfail)
      call assert_close('midpoint P: y_lo', y_lo, 0.5d0*(pd% y_lo(i1) + pd% y_lo(i2)), rtol, atol, nfail)
      call assert_close('midpoint P: y_hi', y_hi, 0.5d0*(pd% y_hi(i1) + pd% y_hi(i2)), rtol, atol, nfail)

      ! bilinear midpoint: average of four nodes
      i1 = find_row(pd, 1.5d0, 5000d0)
      i2 = find_row(pd, 1.5d0, 6000d0)
      i3 = find_row(pd, 2d0, 5000d0)
      i4 = find_row(pd, 2d0, 6000d0)
      call phase_diagram_limits(pd, 1.75d0*Mbar, 5500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('bilinear midpoint: ok and gap', ierr == 0 .and. has_gap, nfail)
      call assert_close('bilinear midpoint: y_lo', y_lo, &
         0.25d0*(pd% y_lo(i1) + pd% y_lo(i2) + pd% y_lo(i3) + pd% y_lo(i4)), rtol, atol, nfail)
      call assert_close('bilinear midpoint: y_hi', y_hi, &
         0.25d0*(pd% y_hi(i1) + pd% y_hi(i2) + pd% y_hi(i3) + pd% y_hi(i4)), rtol, atol, nfail)

      ! general bilinear point: fT = 0.25, fP = 0.2
      fT = 0.25d0
      fP = 0.2d0
      call phase_diagram_limits(pd, (1.5d0 + fP*0.5d0)*Mbar, 5000d0 + fT*1000d0, y_lo, y_hi, has_gap, ierr)
      e_lo = (1d0 - fP)*((1d0 - fT)*pd% y_lo(i1) + fT*pd% y_lo(i2)) + fP*((1d0 - fT)*pd% y_lo(i3) + fT*pd% y_lo(i4))
      e_hi = (1d0 - fP)*((1d0 - fT)*pd% y_hi(i1) + fT*pd% y_hi(i2)) + fP*((1d0 - fT)*pd% y_hi(i3) + fT*pd% y_hi(i4))
      call assert_true('bilinear general: ok and gap', ierr == 0 .and. has_gap, nfail)
      call assert_close('bilinear general: y_lo', y_lo, e_lo, rtol, atol, nfail)
      call assert_close('bilinear general: y_hi', y_hi, e_hi, rtol, atol, nfail)

      ! inside the last P interval (10 .. 24 Mbar)
      i1 = find_row(pd, 10d0, 8000d0)
      i2 = find_row(pd, 24d0, 8000d0)
      call phase_diagram_limits(pd, 17d0*Mbar, 8000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('last P interval: ok and gap', ierr == 0 .and. has_gap, nfail)
      call assert_close('last P interval: y_lo', y_lo, 0.5d0*(pd% y_lo(i1) + pd% y_lo(i2)), rtol, atol, nfail)
      call assert_close('last P interval: y_hi', y_hi, 0.5d0*(pd% y_hi(i1) + pd% y_hi(i2)), rtol, atol, nfail)
   end subroutine test_midpoints


   subroutine test_out_of_range(pd, nfail)
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      type(phase_diagram), intent(in) :: pd
      integer, intent(inout) :: nfail
      type(phase_diagram) :: empty
      integer :: ierr
      real(dp) :: y_lo, y_hi, nan
      logical :: has_gap

      call phase_diagram_limits(pd, 0.4d0*Mbar, 2000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('P below range: no gap', ierr == 0 .and. .not. has_gap .and. y_lo == -1d0 .and. y_hi == -1d0, nfail)
      call phase_diagram_limits(pd, -1d0, 2000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('negative P: no gap', ierr == 0 .and. .not. has_gap, nfail)
      call phase_diagram_limits(pd, 24.5d0*Mbar, 5000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('P above range: no gap', ierr == 0 .and. .not. has_gap .and. y_lo == -1d0 .and. y_hi == -1d0, nfail)
      call phase_diagram_limits(pd, 0.5d0*Mbar, 1500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('P at lower edge: gap', ierr == 0 .and. has_gap, nfail)
      call phase_diagram_limits(pd, 24d0*Mbar, 10500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('P at upper edge: gap', ierr == 0 .and. has_gap, nfail)

      call phase_diagram_limits(pd, 1.5d0*Mbar, 999.999d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T below range: no gap', ierr == 0 .and. .not. has_gap .and. y_lo == -1d0, nfail)
      call phase_diagram_limits(pd, 1.5d0*Mbar, 1000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T at lower edge: gap', ierr == 0 .and. has_gap, nfail)

      ! 0.5 Mbar block ends at 2000 K: T above either bracketing block gives no gap
      call phase_diagram_limits(pd, 0.75d0*Mbar, 2500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T above lower block: no gap', ierr == 0 .and. .not. has_gap .and. y_hi == -1d0, nfail)
      ! at P = 1 Mbar exactly the 1 Mbar block (max 3000 K) still limits the range
      call phase_diagram_limits(pd, 1d0*Mbar, 3500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T above block at node P: no gap', ierr == 0 .and. .not. has_gap, nfail)
      call phase_diagram_limits(pd, 1.2d0*Mbar, 6000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T at upper edge of block: gap', ierr == 0 .and. has_gap, nfail)
      call phase_diagram_limits(pd, 1.2d0*Mbar, 6000.001d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T just above block: no gap', ierr == 0 .and. .not. has_gap, nfail)
      ! T above the upper bracketing block only
      call phase_diagram_limits(pd, 3d0*Mbar, 9500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T above upper block only: no gap', ierr == 0 .and. .not. has_gap, nfail)

      nan = ieee_value(1d0, ieee_quiet_nan)
      call phase_diagram_limits(pd, 1.5d0*Mbar, nan, y_lo, y_hi, has_gap, ierr)
      call assert_true('NaN T: error', ierr == phase_diagram_err_bad_input .and. .not. has_gap, nfail)
      call phase_diagram_limits(pd, nan, 3000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('NaN P: error', ierr == phase_diagram_err_bad_input .and. .not. has_gap, nfail)

      call phase_diagram_limits(empty, 1.5d0*Mbar, 3000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('unloaded diagram: error', ierr == phase_diagram_err_not_loaded .and. .not. has_gap, nfail)
      call assert_true('error message text', &
         len_trim(phase_diagram_error_message(phase_diagram_err_not_loaded)) > 0 &
         .and. trim(phase_diagram_error_message(0)) == 'ok', nfail)
   end subroutine test_out_of_range


   subroutine test_extrapolation(pd, nfail)
      type(phase_diagram), intent(in) :: pd
      integer, intent(inout) :: nfail
      type(phase_diagram) :: pdx
      integer :: ierr, i1, i2
      real(dp) :: y_lo, y_hi, fT
      logical :: has_gap

      call load_phase_diagram(pdx, data_file, 0d0, .true., ierr)
      call assert_true('extrapolation: load', ierr == 0 .and. pdx% extrapolate_in_T, nfail)
      if (ierr /= 0) return

      ! within range the result is unchanged
      call phase_diagram_limits(pdx, 1.75d0*Mbar, 5500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('extrapolation: in-range gap', ierr == 0 .and. has_gap, nfail)
      call check_same(pd, pdx, 1.75d0*Mbar, 5500d0, 'extrapolation: in-range value unchanged', nfail)

      ! P = 1 Mbar, T = 3500 K: 1 Mbar block ends at 3000 K; extrapolate from 2000 and 3000 K
      i1 = find_row(pdx, 1d0, 2000d0)
      i2 = find_row(pdx, 1d0, 3000d0)
      fT = 1.5d0
      call phase_diagram_limits(pdx, 1d0*Mbar, 3500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('extrapolation: gap above block', ierr == 0 .and. has_gap, nfail)
      call assert_close('extrapolation: y_lo', y_lo, (1d0 - fT)*pdx% y_lo(i1) + fT*pdx% y_lo(i2), rtol, atol, nfail)
      call assert_close('extrapolation: y_hi', y_hi, (1d0 - fT)*pdx% y_hi(i1) + fT*pdx% y_hi(i2), rtol, atol, nfail)

      ! far above the 0.5 Mbar block the extrapolated gap closes (y_lo >= y_hi): no gap
      call phase_diagram_limits(pdx, 0.5d0*Mbar, 8000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('extrapolation: closed gap', ierr == 0 .and. .not. has_gap .and. y_lo == -1d0 .and. y_hi == -1d0, nfail)

      ! below the range there is still no gap
      call phase_diagram_limits(pdx, 1.5d0*Mbar, 500d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('extrapolation: T below range still no gap', ierr == 0 .and. .not. has_gap, nfail)
      ! P outside the range is still no gap
      call phase_diagram_limits(pdx, 30d0*Mbar, 5000d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('extrapolation: P above range still no gap', ierr == 0 .and. .not. has_gap, nfail)

      call free_phase_diagram(pdx)

      ! clamping to [0,1] with a synthetic table where both limits leave [0,1]
      call write_tmp_file([character(len=64) :: &
         '"He-poor","He-rich","T","P"', &
         '0.3,0.5,1000,1', &
         '0.1,0.9,2000,1', &
         '0.3,0.5,1000,2', &
         '0.1,0.9,2000,2'])
      call load_phase_diagram(pdx, tmp_file, 0d0, .true., ierr)
      call assert_true('clamping: load synthetic table', ierr == 0, nfail)
      if (ierr == 0) then
         call phase_diagram_limits(pdx, 1.5d0*Mbar, 3000d0, y_lo, y_hi, has_gap, ierr)
         call assert_true('clamping: gap', ierr == 0 .and. has_gap, nfail)
         call assert_true('clamping: y_lo clamped to 0', y_lo == 0d0, nfail)
         call assert_true('clamping: y_hi clamped to 1', y_hi == 1d0, nfail)
         call phase_diagram_limits(pdx, 1.5d0*Mbar, 2500d0, y_lo, y_hi, has_gap, ierr)
         call assert_true('clamping: partial', ierr == 0 .and. has_gap, nfail)
         call assert_close('clamping: y_lo at 2500 K', y_lo, 0d0, 0d0, atol, nfail)
         call assert_close('clamping: y_hi at 2500 K', y_hi, 1d0, 0d0, atol, nfail)
         call free_phase_diagram(pdx)
      end if
      ! same table without extrapolation: no gap above the range
      call load_phase_diagram(pdx, tmp_file, 0d0, .false., ierr)
      call assert_true('clamping: reload without extrapolation', ierr == 0, nfail)
      if (ierr == 0) then
         call phase_diagram_limits(pdx, 1.5d0*Mbar, 3000d0, y_lo, y_hi, has_gap, ierr)
         call assert_true('no extrapolation: no gap above range', ierr == 0 .and. .not. has_gap, nfail)
         call free_phase_diagram(pdx)
      end if
      call delete_tmp_file()
   end subroutine test_extrapolation


   subroutine test_temperature_shift(pd, nfail)
      type(phase_diagram), intent(in) :: pd
      integer, intent(inout) :: nfail
      type(phase_diagram) :: pds
      integer :: ierr
      real(dp) :: y_lo, y_hi, y_lo_ref, y_hi_ref
      logical :: has_gap, has_gap_ref

      call load_phase_diagram(pds, data_file, 500d0, .false., ierr)
      call assert_true('T shift: load', ierr == 0 .and. pds% T_shift == 500d0, nfail)
      if (ierr /= 0) return
      call assert_true('T shift: table T shifted', pds% T(1) == 1500d0 .and. all(pds% T == pd% T + 500d0), nfail)
      call assert_true('T shift: P unchanged', all(pds% P == pd% P), nfail)

      ! limits at T + 500 equal unshifted limits at T
      call phase_diagram_limits(pd, 1.6d0*Mbar, 5250d0, y_lo_ref, y_hi_ref, has_gap_ref, ierr)
      call phase_diagram_limits(pds, 1.6d0*Mbar, 5750d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T shift: gap consistent', ierr == 0 .and. has_gap .and. has_gap_ref, nfail)
      call assert_close('T shift: y_lo consistent', y_lo, y_lo_ref, rtol, atol, nfail)
      call assert_close('T shift: y_hi consistent', y_hi, y_hi_ref, rtol, atol, nfail)

      ! the shifted range starts at 1500 K
      call phase_diagram_limits(pds, 1.5d0*Mbar, 1200d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T shift: below shifted range', ierr == 0 .and. .not. has_gap, nfail)
      call phase_diagram_limits(pd, 1.5d0*Mbar, 1200d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('T shift: unshifted has gap there', ierr == 0 .and. has_gap, nfail)
      call free_phase_diagram(pds)

      ! negative shift
      call load_phase_diagram(pds, data_file, -200d0, .false., ierr)
      call assert_true('negative T shift: load', ierr == 0, nfail)
      if (ierr /= 0) return
      call phase_diagram_limits(pds, 1.5d0*Mbar, 900d0, y_lo, y_hi, has_gap, ierr)
      call assert_true('negative T shift: gap at 900 K', ierr == 0 .and. has_gap, nfail)
      call phase_diagram_limits(pd, 1.5d0*Mbar, 1100d0, y_lo_ref, y_hi_ref, has_gap_ref, ierr)
      call assert_close('negative T shift: y_lo consistent', y_lo, y_lo_ref, rtol, atol, nfail)
      call assert_close('negative T shift: y_hi consistent', y_hi, y_hi_ref, rtol, atol, nfail)
      call free_phase_diagram(pds)
   end subroutine test_temperature_shift


   subroutine test_check_and_print(pd, nfail)
      type(phase_diagram), intent(in) :: pd
      integer, intent(inout) :: nfail
      type(phase_diagram) :: empty
      integer :: ierr, iounit, ios, nlines
      character(len=256) :: line

      call check_phase_diagram(pd, ierr)
      call assert_true('check_phase_diagram on valid table', ierr == 0, nfail)

      write(*,'(a)') '   (an error message about an unloaded diagram is expected next)'
      call check_phase_diagram(empty, ierr)
      call assert_true('check_phase_diagram on unloaded table', ierr == phase_diagram_err_not_loaded, nfail)

      open(newunit=iounit, status='scratch', action='readwrite', form='formatted')
      call print_phase_diagram(pd, iounit)
      rewind(iounit)
      nlines = 0
      do
         read(iounit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         nlines = nlines + 1
      end do
      close(iounit)
      ! 7 summary lines + one per P block + table header + rows + closing line
      call assert_true('print_phase_diagram line count', nlines == 7 + pd% num_unique_P + 1 + pd% n_points + 1, nfail)

      open(newunit=iounit, status='scratch', action='readwrite', form='formatted')
      call print_phase_diagram(empty, iounit)
      rewind(iounit)
      nlines = 0
      do
         read(iounit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         nlines = nlines + 1
      end do
      close(iounit)
      call assert_true('print_phase_diagram of unloaded table', nlines == 8, nfail)
   end subroutine test_check_and_print


   subroutine test_reload(nfail)
      integer, intent(inout) :: nfail
      type(phase_diagram) :: pd
      integer :: ierr
      call load_phase_diagram(pd, data_file, 0d0, .false., ierr)
      call assert_true('reload: first load', ierr == 0, nfail)
      call load_phase_diagram(pd, data_file, 100d0, .true., ierr)
      call assert_true('reload: second load into the same variable', ierr == 0 .and. pd% loaded &
         .and. pd% n_points == expected_n_points .and. pd% T_shift == 100d0 .and. pd% extrapolate_in_T, nfail)
      call free_phase_diagram(pd)
   end subroutine test_reload


   subroutine test_openmp(pd, nfail)
      ! phase_diagram_limits from an OpenMP loop gives the serial result
      type(phase_diagram), intent(in) :: pd
      integer, intent(inout) :: nfail
      integer, parameter :: nP = 120, nT = 130
      real(dp), allocatable :: Ps(:), Ts(:), y_lo_s(:,:), y_hi_s(:,:), y_lo_p(:,:), y_hi_p(:,:)
      logical, allocatable :: gap_s(:,:), gap_p(:,:)
      integer, allocatable :: ierr_p(:,:)
      integer :: i, j, ierr, num_bad
      logical :: same

      allocate(Ps(nP), Ts(nT), y_lo_s(nT,nP), y_hi_s(nT,nP), y_lo_p(nT,nP), y_hi_p(nT,nP), &
         gap_s(nT,nP), gap_p(nT,nP), ierr_p(nT,nP))
      do i = 1, nP
         Ps(i) = (0.3d0 + (26d0 - 0.3d0)*dble(i-1)/dble(nP-1))*Mbar
      end do
      do j = 1, nT
         Ts(j) = 500d0 + (13000d0 - 500d0)*dble(j-1)/dble(nT-1)
      end do

      num_bad = 0
      do i = 1, nP
         do j = 1, nT
            call phase_diagram_limits(pd, Ps(i), Ts(j), y_lo_s(j,i), y_hi_s(j,i), gap_s(j,i), ierr)
            if (ierr /= 0) num_bad = num_bad + 1
         end do
      end do
      call assert_true('openmp: serial reference has no errors', num_bad == 0, nfail)
      call assert_true('openmp: grid has gap and no-gap points', &
         count(gap_s) > 0 .and. count(.not. gap_s) > 0, nfail)

      !$OMP PARALLEL DO PRIVATE(i, j) SCHEDULE(dynamic, 4)
      do i = 1, nP
         do j = 1, nT
            call phase_diagram_limits(pd, Ps(i), Ts(j), y_lo_p(j,i), y_hi_p(j,i), gap_p(j,i), ierr_p(j,i))
         end do
      end do
      !$OMP END PARALLEL DO

      same = all(ierr_p == 0) .and. all(gap_p .eqv. gap_s) .and. all(y_lo_p == y_lo_s) .and. all(y_hi_p == y_hi_s)
      call assert_true('openmp: parallel equals serial', same, nfail)
   end subroutine test_openmp


   subroutine test_tolerant_formats(nfail)
      ! CRLF endings, blank lines, blanks and tabs around numbers, exponents
      integer, intent(inout) :: nfail
      type(phase_diagram) :: pd
      integer :: ierr
      real(dp) :: y_lo, y_hi
      logical :: has_gap
      character(len=1) :: cr, tab

      cr = achar(13)
      tab = achar(9)
      call write_tmp_file([character(len=64) :: &
         '',  &
         ' "He-poor" , "He-rich" , "T" , "P" ' // cr, &
         '0.1, 0.5, 1d3, 1' // cr, &
         '', &
         '  0.2 ,' // tab // '0.6, +2.0E3 ,1.0' // cr, &
         '1.0e-1,5e-1,1000,2.' // cr, &
         '.2,.6,2000,2' // cr, &
         '', &
         '   '])
      call load_phase_diagram(pd, tmp_file, 0d0, .false., ierr)
      call assert_true('tolerant formats: load', ierr == 0 .and. pd% n_points == 4 .and. pd% num_unique_P == 2, nfail)
      if (ierr == 0) then
         call assert_true('tolerant formats: values', &
            pd% y_lo(1) == 0.1d0 .and. pd% y_hi(2) == 0.6d0 .and. pd% T(2) == 2000d0 .and. pd% P(3) == 2d0*Mbar &
            .and. pd% T(1) == 1000d0 .and. pd% y_lo(3) == 0.1d0, nfail)
         call phase_diagram_limits(pd, 1.5d0*Mbar, 1500d0, y_lo, y_hi, has_gap, ierr)
         call assert_true('tolerant formats: limits', ierr == 0 .and. has_gap, nfail)
         call assert_close('tolerant formats: y_lo', y_lo, 0.15d0, rtol, atol, nfail)
         call assert_close('tolerant formats: y_hi', y_hi, 0.55d0, rtol, atol, nfail)
      end if
      call free_phase_diagram(pd)
      call delete_tmp_file()
   end subroutine test_tolerant_formats


   subroutine test_malformed_files(nfail)
      integer, intent(inout) :: nfail
      type(phase_diagram) :: pd
      integer :: ierr
      character(len=*), parameter :: header = '"He-poor","He-rich","T","P"'

      write(*,'(a)') '   (error messages from load_phase_diagram are expected below)'

      call expect_failure('missing column', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2,0.6,2000', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_file, nfail)
      call expect_failure('extra column', [character(len=64) :: header, &
         '0.1,0.5,1000,1,7', '0.2,0.6,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_file, nfail)
      call expect_failure('trailing comma', [character(len=64) :: header, &
         '0.1,0.5,1000,1,', '0.2,0.6,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_file, nfail)
      call expect_failure('non-numeric field', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2,abc,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_file, nfail)
      call expect_failure('two numbers in one field', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2 0.3,0.6,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_file, nfail)
      call expect_failure('empty field', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2,,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_file, nfail)
      call expect_failure('missing header', [character(len=64) :: &
         '0.1,0.5,1000,1', '0.2,0.6,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_file, nfail)
      call expect_failure('header only', [character(len=64) :: header], phase_diagram_err_file, nfail)
      call expect_failure('empty file', [character(len=64) :: ''], phase_diagram_err_file, nfail)
      call expect_failure('non-monotonic T in block', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2,0.6,3000,1', '0.3,0.7,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], &
         phase_diagram_err_structure, nfail)
      call expect_failure('repeated T in block', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2,0.6,2000,1', '0.3,0.7,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], &
         phase_diagram_err_structure, nfail)
      call expect_failure('P decreasing between blocks', [character(len=64) :: header, &
         '0.1,0.5,1000,2', '0.2,0.6,2000,2', '0.1,0.5,1000,1', '0.2,0.6,2000,1'], &
         phase_diagram_err_structure, nfail)
      call expect_failure('P block repeated later', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2,0.6,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2', '0.1,0.5,1000,1', '0.2,0.6,2000,1'], &
         phase_diagram_err_structure, nfail)
      call expect_failure('y_lo >= y_hi', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.6,0.6,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_structure, nfail)
      call expect_failure('y_hi > 1', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2,1.6,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_structure, nfail)
      call expect_failure('y_lo < 0', [character(len=64) :: header, &
         '-0.1,0.5,1000,1', '0.2,0.6,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_structure, nfail)
      call expect_failure('single P block', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2,0.6,2000,1'], phase_diagram_err_structure, nfail)
      call expect_failure('block with one T point', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_structure, nfail)
      call expect_failure('non-positive P', [character(len=64) :: header, &
         '0.1,0.5,1000,0', '0.2,0.6,2000,0', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_structure, nfail)
      call expect_failure('T not positive after shift', [character(len=64) :: header, &
         '0.1,0.5,1000,1', '0.2,0.6,2000,1', '0.1,0.5,1000,2', '0.2,0.6,2000,2'], phase_diagram_err_structure, nfail, &
         T_shift=-1000d0)

      call load_phase_diagram(pd, 'test_phase_diagram_does_not_exist.csv', 0d0, .false., ierr)
      call assert_true('malformed: nonexistent file', ierr == phase_diagram_err_file .and. .not. pd% loaded, nfail)
   end subroutine test_malformed_files


   subroutine expect_failure(name, lines, expected_ierr, nfail, T_shift)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: expected_ierr
      integer, intent(inout) :: nfail
      real(dp), intent(in), optional :: T_shift
      type(phase_diagram) :: pd
      real(dp) :: shift
      integer :: ierr
      shift = 0d0
      if (present(T_shift)) shift = T_shift
      call write_tmp_file(lines)
      call load_phase_diagram(pd, tmp_file, shift, .false., ierr)
      call assert_true('malformed: ' // name, ierr == expected_ierr .and. .not. pd% loaded &
         .and. pd% n_points == 0 .and. .not. allocated(pd% T), nfail)
      call delete_tmp_file()
   end subroutine expect_failure


   subroutine write_tmp_file(lines)
      character(len=*), intent(in) :: lines(:)
      integer :: iounit, i
      open(newunit=iounit, file=tmp_file, status='replace', action='write')
      do i = 1, size(lines)
         write(iounit, '(a)') trim(lines(i))
      end do
      close(iounit)
   end subroutine write_tmp_file


   subroutine delete_tmp_file()
      integer :: iounit, ios
      open(newunit=iounit, file=tmp_file, status='old', iostat=ios)
      if (ios == 0) close(iounit, status='delete')
   end subroutine delete_tmp_file


   integer function find_row(pd, P_Mbar, T) result(i)
      ! index of the table row with the given P [Mbar] and T [K]; 0 if absent
      type(phase_diagram), intent(in) :: pd
      real(dp), intent(in) :: P_Mbar, T
      integer :: k
      i = 0
      do k = 1, pd% n_points
         if (abs(pd% P(k) - P_Mbar*Mbar) <= 1d-6*Mbar .and. abs(pd% T(k) - T) <= 1d-6) then
            i = k
            return
         end if
      end do
   end function find_row


   subroutine check_same(pd1, pd2, P, T, name, nfail)
      type(phase_diagram), intent(in) :: pd1, pd2
      real(dp), intent(in) :: P, T
      character(len=*), intent(in) :: name
      integer, intent(inout) :: nfail
      real(dp) :: y_lo1, y_hi1, y_lo2, y_hi2
      logical :: g1, g2
      integer :: ierr1, ierr2
      call phase_diagram_limits(pd1, P, T, y_lo1, y_hi1, g1, ierr1)
      call phase_diagram_limits(pd2, P, T, y_lo2, y_hi2, g2, ierr2)
      call assert_true(name, ierr1 == 0 .and. ierr2 == 0 .and. (g1 .eqv. g2) &
         .and. abs(y_lo1 - y_lo2) <= tight_atol .and. abs(y_hi1 - y_hi2) <= tight_atol, nfail)
   end subroutine check_same

end module test_phase_diagram
