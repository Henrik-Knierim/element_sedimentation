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

module test_miscibility

   ! unit tests for sedimentation_miscibility. Run from the test directory.

   use const_def, only: dp
   use sedimentation_phase_diagram
   use sedimentation_miscibility
   use test_support, only: assert_true, assert_close

   implicit none

   private
   public :: run_miscibility_tests

   character(len=*), parameter :: data_file = 'data/Schoettler+2018.csv'
   real(dp), parameter :: Mbar = 1d12
   real(dp), parameter :: A_h1 = 1.00794d0, A_he4 = 4.002602d0, A_c12 = 12.0107d0, A_o16 = 15.9994d0
   real(dp), parameter :: rtol = 1d-13, atol = 1d-15

contains

   subroutine run_miscibility_tests(nfail)
      integer, intent(inout) :: nfail
      write(*,'(a)') 'test_miscibility'
      call test_pair_number_fraction(nfail)
      call test_two_species(nfail)
      call test_four_species(nfail)
      call test_round_trip(nfail)
      call test_old_algorithm_equivalence(nfail)
      call test_degenerate_cases(nfail)
      call test_profile(nfail)
      call test_degenerate_cells_not_in_gap(nfail)
      call test_errors(nfail)
      call test_large_profile(nfail)
   end subroutine run_miscibility_tests


   subroutine test_pair_number_fraction(nfail)
      integer, intent(inout) :: nfail
      real(dp) :: xa2(2), A2(2), xa4(4), A4(4), expected

      A2 = [A_h1, A_he4]
      xa2 = [0.75d0, 0.25d0]
      expected = (0.25d0/A_he4)/(0.25d0/A_he4 + 0.75d0/A_h1)
      call assert_close('pair_number_fraction: H-He', pair_number_fraction(xa2, A2, 2, 1), expected, rtol, atol, nfail)
      call assert_close('pair_number_fraction: complement', pair_number_fraction(xa2, A2, 1, 2), 1d0 - expected, rtol, atol, nfail)
      xa2 = [0d0, 0d0]
      call assert_true('pair_number_fraction: both zero gives 0', pair_number_fraction(xa2, A2, 2, 1) == 0d0, nfail)
      xa2 = [0d0, 1d0]
      call assert_true('pair_number_fraction: pure He gives 1', pair_number_fraction(xa2, A2, 2, 1) == 1d0, nfail)
      xa2 = [1d0, 0d0]
      call assert_true('pair_number_fraction: pure H gives 0', pair_number_fraction(xa2, A2, 2, 1) == 0d0, nfail)

      ! metals do not enter the pair fraction
      A4 = [A_h1, A_he4, A_c12, A_o16]
      xa4 = [0.7d0, 0.25d0, 0.03d0, 0.02d0]
      expected = (0.25d0/A_he4)/(0.25d0/A_he4 + 0.7d0/A_h1)
      call assert_close('pair_number_fraction: 4 species', pair_number_fraction(xa4, A4, 2, 1), expected, rtol, atol, nfail)
      ! unnormalized input gives the same fraction
      call assert_close('pair_number_fraction: scale invariant', pair_number_fraction(3d0*xa4, A4, 2, 1), &
         expected, rtol, atol, nfail)
   end subroutine test_pair_number_fraction


   subroutine test_two_species(nfail)
      ! with only H and He: X_He = y A_He / (y A_He + (1 - y) A_H), independent of the start
      integer, intent(inout) :: nfail
      real(dp) :: xa(2), A(2), y, x, expected, y_cur
      integer :: i
      real(dp), parameter :: ys(5) = [0d0, 0.1d0, 0.5d0, 0.9d0, 1d0]
      character(len=64) :: name

      A = [A_h1, A_he4]
      xa = [0.75d0, 0.25d0]
      do i = 1, size(ys)
         y = ys(i)
         x = mass_fraction_for_pair_number_fraction(xa, A, 2, 1, y)
         expected = y*A_he4/(y*A_he4 + (1d0 - y)*A_h1)
         write(name, '(a,1pd10.3)') 'two species: X_He at y = ', y
         call assert_close(name, x, expected, rtol, atol, nfail)
         ! same answer from a different starting composition
         call assert_close(trim(name) // ' (other start)', &
            mass_fraction_for_pair_number_fraction([0.1d0, 0.9d0], A, 2, 1, y), expected, rtol, atol, nfail)
      end do
      ! requesting the current fraction returns the current mass fraction
      y_cur = pair_number_fraction(xa, A, 2, 1)
      call assert_close('two species: identity', mass_fraction_for_pair_number_fraction(xa, A, 2, 1, y_cur), &
         0.25d0, rtol, atol, nfail)
      call assert_true('two species: X(0) = 0 and X(1) = 1', &
         mass_fraction_for_pair_number_fraction(xa, A, 2, 1, 0d0) == 0d0 &
         .and. mass_fraction_for_pair_number_fraction(xa, A, 2, 1, 1d0) == 1d0, nfail)
   end subroutine test_two_species


   subroutine test_four_species(nfail)
      ! the recipe of the specification evaluated explicitly
      integer, intent(inout) :: nfail
      real(dp) :: xa(4), A(4), n(4), n_new(4), y_t, y_cur, f, expected, x, mass_new(4)
      integer :: j

      A = [A_h1, A_he4, A_c12, A_o16]
      xa = [0.70d0, 0.26d0, 0.03d0, 0.01d0]
      y_t = 0.15d0
      n = xa/A
      y_cur = n(2)/(n(1) + n(2))
      f = (1d0 - y_t)/(1d0 - y_cur)
      do j = 1, 4
         n_new(j) = f*n(j)
      end do
      n_new(2) = y_t*(n_new(1))/(1d0 - y_t)   ! n_i/(n_i + f n_p) = y_t
      mass_new = n_new*A
      expected = mass_new(2)/sum(mass_new)

      x = mass_fraction_for_pair_number_fraction(xa, A, 2, 1, y_t)
      call assert_close('four species: X_He for y = 0.15', x, expected, rtol, atol, nfail)
      call assert_true('four species: result in [0,1]', x >= 0d0 .and. x <= 1d0, nfail)
      ! number ratios of C and O to H are unchanged in the new composition by construction
      call assert_close('four species: C/H ratio kept', n_new(3)/n_new(1), n(3)/n(1), rtol, atol, nfail)
      call assert_close('four species: O/H ratio kept', n_new(4)/n_new(1), n(4)/n(1), rtol, atol, nfail)
      ! the new composition has the requested pair fraction
      call assert_close('four species: new pair fraction', pair_number_fraction(mass_new, A, 2, 1), y_t, rtol, atol, nfail)
      ! monotonic in y
      call assert_true('four species: monotonic in y', &
         mass_fraction_for_pair_number_fraction(xa, A, 2, 1, 0.1d0) < x .and. &
         x < mass_fraction_for_pair_number_fraction(xa, A, 2, 1, 0.2d0), nfail)
   end subroutine test_four_species


   subroutine test_round_trip(nfail)
      ! X -> composition with others scaled -> pair_number_fraction == y_target
      integer, intent(inout) :: nfail
      real(dp) :: A(4), xa(4), xa_new(4), x, y_back, others, max_err
      real(dp), parameter :: ys(7) = [0d0, 0.01d0, 0.05d0, 0.3d0, 0.5d0, 0.95d0, 1d0]
      integer :: i, c, j
      real(dp) :: xa_cases(4,4)

      A = [A_h1, A_he4, A_c12, A_o16]
      xa_cases(:,1) = [0.70d0, 0.26d0, 0.03d0, 0.01d0]
      xa_cases(:,2) = [0.98d0, 0.01d0, 0.005d0, 0.005d0]
      xa_cases(:,3) = [0.02d0, 0.97d0, 0.005d0, 0.005d0]
      xa_cases(:,4) = [0.50d0, 0.0d0, 0.25d0, 0.25d0]

      max_err = 0d0
      do c = 1, 4
         xa = xa_cases(:,c)
         do i = 1, size(ys)
            x = mass_fraction_for_pair_number_fraction(xa, A, 2, 1, ys(i))
            others = 0d0
            do j = 1, 4
               if (j /= 2) others = others + xa(j)
            end do
            xa_new = xa*(1d0 - x)/others
            xa_new(2) = x
            y_back = pair_number_fraction(xa_new, A, 2, 1)
            max_err = max(max_err, abs(y_back - ys(i)))
            if (abs(sum(xa_new) - 1d0) > 1d-14) max_err = 1d0
         end do
      end do
      call assert_close('round trip: max |y_back - y_target| over 28 cases', max_err, 0d0, 0d0, 1d-14, nfail)
      ! round trip starting from the pair fraction
      xa = xa_cases(:,1)
      call assert_close('round trip: X(y(xa)) == xa(i)', &
         mass_fraction_for_pair_number_fraction(xa, A, 2, 1, pair_number_fraction(xa, A, 2, 1)), xa(2), rtol, atol, nfail)
   end subroutine test_round_trip


   subroutine test_old_algorithm_equivalence(nfail)
      ! the new closed form reproduces the old from_scaled_number_fraction_to_mass_fraction
      integer, intent(inout) :: nfail
      real(dp) :: A(4), xa(4), max_err
      real(dp), parameter :: ys(5) = [0.02d0, 0.05d0, 0.2d0, 0.5d0, 0.8d0]
      integer :: i

      A = [A_h1, A_he4, A_c12, A_o16]
      xa = [0.70d0, 0.26d0, 0.03d0, 0.01d0]
      max_err = 0d0
      do i = 1, size(ys)
         max_err = max(max_err, abs(mass_fraction_for_pair_number_fraction(xa, A, 2, 1, ys(i)) &
            - old_algorithm(xa, A, 2, 1, ys(i))))
      end do
      call assert_close('old algorithm equivalence: 4 species', max_err, 0d0, 0d0, 1d-14, nfail)
      xa = [0.10d0, 0.80d0, 0.05d0, 0.05d0]
      max_err = 0d0
      do i = 1, size(ys)
         max_err = max(max_err, abs(mass_fraction_for_pair_number_fraction(xa, A, 2, 1, ys(i)) &
            - old_algorithm(xa, A, 2, 1, ys(i))))
      end do
      call assert_close('old algorithm equivalence: He rich', max_err, 0d0, 0d0, 1d-14, nfail)
      ! also with the sedimenting species first and the partner last
      max_err = 0d0
      do i = 1, size(ys)
         max_err = max(max_err, abs(mass_fraction_for_pair_number_fraction([0.26d0, 0.03d0, 0.01d0, 0.70d0], &
            [A_he4, A_c12, A_o16, A_h1], 1, 4, ys(i)) - old_algorithm([0.26d0, 0.03d0, 0.01d0, 0.70d0], &
            [A_he4, A_c12, A_o16, A_h1], 1, 4, ys(i))))
      end do
      call assert_close('old algorithm equivalence: permuted indices', max_err, 0d0, 0d0, 1d-14, nfail)
   end subroutine test_old_algorithm_equivalence


   real(dp) function old_algorithm(xa, A, i_sed, i_partner, y_m) result(x)
      ! transcription of the old from_scaled_number_fraction_to_mass_fraction
      real(dp), intent(in) :: xa(:), A(:), y_m
      integer, intent(in) :: i_sed, i_partner
      real(dp), allocatable :: yi(:)
      real(dp) :: scaling, dyi, f, mmw
      allocate(yi(size(xa)))
      yi = xa/A
      scaling = yi(i_sed) + yi(i_partner)
      yi = yi/scaling
      dyi = y_m - yi(i_sed)
      f = 1d0 - dyi/(1d0 - yi(i_sed))
      yi = f*yi
      yi(i_sed) = y_m
      yi = yi/sum(yi)
      mmw = sum(yi*A)
      x = yi(i_sed)*A(i_sed)/mmw
   end function old_algorithm


   subroutine test_degenerate_cases(nfail)
      integer, intent(inout) :: nfail
      real(dp) :: A(4)
      A = [A_h1, A_he4, A_c12, A_o16]
      ! no partner: a pair fraction below 1 needs X_He = 0
      call assert_true('degenerate: no partner, y < 1', &
         mass_fraction_for_pair_number_fraction([0d0, 0.9d0, 0.05d0, 0.05d0], A, 2, 1, 0.3d0) == 0d0, nfail)
      ! no partner and y = 1: current fraction returned
      call assert_close('degenerate: no partner, y = 1', &
         mass_fraction_for_pair_number_fraction([0d0, 0.9d0, 0.05d0, 0.05d0], A, 2, 1, 1d0), 0.9d0, rtol, atol, nfail)
      ! pure He: y = 0 gives 0; any other y is unreachable and returns the current fraction 1
      call assert_true('degenerate: pure He', &
         mass_fraction_for_pair_number_fraction([0d0, 1d0, 0d0, 0d0], A, 2, 1, 1d0) == 1d0 .and. &
         mass_fraction_for_pair_number_fraction([0d0, 1d0, 0d0, 0d0], A, 2, 1, 0.5d0) == 1d0 .and. &
         mass_fraction_for_pair_number_fraction([0d0, 1d0, 0d0, 0d0], A, 2, 1, 0d0) == 0d0, nfail)
      ! all zero
      call assert_true('degenerate: all zero', &
         mass_fraction_for_pair_number_fraction([0d0, 0d0, 0d0, 0d0], A, 2, 1, 0.5d0) == 0d0, nfail)
      ! unnormalized input is normalized in the answer
      call assert_close('degenerate: unnormalized input', &
         mass_fraction_for_pair_number_fraction(2d0*[0.70d0, 0.26d0, 0.03d0, 0.01d0], A, 2, 1, 0.15d0), &
         mass_fraction_for_pair_number_fraction([0.70d0, 0.26d0, 0.03d0, 0.01d0], A, 2, 1, 0.15d0), rtol, atol, nfail)
   end subroutine test_degenerate_cases


   subroutine test_profile(nfail)
      ! synthetic 12-cell profile: cells 1-3 outside the P range, 4-8 in the gap,
      ! 9-10 He-poor side of the gap, 11 He-rich side, 12 outside the P range
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 12, species = 4
      real(dp) :: xa(species, nz), A(species), P(nz), T(nz), x_lo(nz), x_hi(nz), y_lo, y_hi, xa_new(species), others
      real(dp) :: xa_perm(species, nz), A_perm(species), x_lo_perm(nz), x_hi_perm(nz)
      logical :: has_gap(nz), in_gap(nz), expected_has_gap(nz), expected_in_gap(nz), has_gap_perm(nz), in_gap_perm(nz), gap_k
      type(phase_diagram) :: pd
      integer :: ierr, k, j
      real(dp), parameter :: P_Mbar(nz) = [0.1d0, 0.3d0, 0.49d0, 1.5d0, 1.75d0, 2d0, 2.5d0, 3d0, 4d0, 10d0, 10d0, 30d0]
      real(dp), parameter :: T_K(nz) = &
         [500d0, 800d0, 1500d0, 3000d0, 3500d0, 4000d0, 4500d0, 5000d0, 9000d0, 10000d0, 6000d0, 12000d0]

      call load_phase_diagram(pd, data_file, 0d0, .false., ierr)
      call assert_true('profile: load phase diagram', ierr == 0, nfail)
      if (ierr /= 0) return

      A = [A_h1, A_he4, A_c12, A_o16]
      do k = 1, nz
         xa(:,k) = [0.71d0, 0.27d0, 0.012d0, 0.008d0]
         P(k) = P_Mbar(k)*Mbar
         T(k) = T_K(k)
      end do
      xa(:,11) = [0.008d0, 0.98d0, 0.007d0, 0.005d0]  ! He rich
      expected_has_gap = [.false., .false., .false., .true., .true., .true., .true., .true., .true., .true., .true., .false.]
      expected_in_gap = [.false., .false., .false., .true., .true., .true., .true., .true., .false., .false., .false., .false.]

      call evaluate_miscibility(pd, nz, species, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('profile: no error', ierr == 0, nfail)
      call assert_true('profile: has_gap pattern', all(has_gap .eqv. expected_has_gap), nfail)
      call assert_true('profile: in_gap pattern', all(in_gap .eqv. expected_in_gap), nfail)
      call assert_true('profile: in_gap exactly where x_lo < xa < x_hi', &
         all(in_gap .eqv. (has_gap .and. x_lo < xa(2,:) .and. xa(2,:) < x_hi)), nfail)
      call assert_true('profile: limits are -1 without gap', &
         all(x_lo == -1d0 .or. has_gap) .and. all(x_hi == -1d0 .or. has_gap), nfail)
      call assert_true('profile: 0 <= x_lo < x_hi <= 1 with gap', &
         all((x_lo >= 0d0 .and. x_lo < x_hi .and. x_hi <= 1d0) .or. .not. has_gap), nfail)
      call assert_true('profile: He-poor cells below x_lo', xa(2,9) <= x_lo(9) .and. xa(2,10) <= x_lo(10), nfail)
      call assert_true('profile: He-rich cell above x_hi', xa(2,11) >= x_hi(11), nfail)

      ! the mass-fraction limits map back to the number-fraction limits of the diagram
      do k = 4, 8
         call phase_diagram_limits(pd, P(k), T(k), y_lo, y_hi, gap_k, ierr)
         others = 0d0
         do j = 1, species
            if (j /= 2) others = others + xa(j,k)
         end do
         xa_new = xa(:,k)*(1d0 - x_lo(k))/others
         xa_new(2) = x_lo(k)
         call assert_close('profile: x_lo maps to y_lo', pair_number_fraction(xa_new, A, 2, 1), y_lo, 1d-12, 1d-14, nfail)
         xa_new = xa(:,k)*(1d0 - x_hi(k))/others
         xa_new(2) = x_hi(k)
         call assert_close('profile: x_hi maps to y_hi', pair_number_fraction(xa_new, A, 2, 1), y_hi, 1d-12, 1d-14, nfail)
      end do

      ! permuted species order (c12, he4, o16, h1) gives the same result
      A_perm = [A_c12, A_he4, A_o16, A_h1]
      do k = 1, nz
         xa_perm(:,k) = [xa(3,k), xa(2,k), xa(4,k), xa(1,k)]
      end do
      call evaluate_miscibility(pd, nz, species, xa_perm, A_perm, 2, 4, P, T, has_gap_perm, in_gap_perm, &
         x_lo_perm, x_hi_perm, ierr)
      call assert_true('profile: permuted species no error', ierr == 0, nfail)
      call assert_true('profile: permuted species same flags', &
         all(has_gap_perm .eqv. has_gap) .and. all(in_gap_perm .eqv. in_gap), nfail)
      call assert_close('profile: permuted species same x_lo', maxval(abs(x_lo_perm - x_lo)), 0d0, 0d0, 1d-14, nfail)
      call assert_close('profile: permuted species same x_hi', maxval(abs(x_hi_perm - x_hi)), 0d0, 0d0, 1d-14, nfail)

      ! oversized arrays (as with s% xa) are fine
      call evaluate_miscibility(pd, 8, species, xa, A, 2, 1, P, T, has_gap_perm, in_gap_perm, x_lo_perm, x_hi_perm, ierr)
      call assert_true('profile: nz smaller than the arrays', ierr == 0 .and. all(in_gap_perm(1:8) .eqv. in_gap(1:8)), nfail)

      call free_phase_diagram(pd)
   end subroutine test_profile


   subroutine test_degenerate_cells_not_in_gap(nfail)
      ! Cells without the partner species (pair number fraction exactly 1) or without either
      ! pair species (y = 0) must never be inside the gap, also when sum(xa) differs from 1
      ! by round-off and the He-rich limit of the table is exactly 1 (P = 2 Mbar, T = 1000 K).
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 4, species = 4
      type(phase_diagram) :: pd
      real(dp) :: xa(species, nz), A(species), P(nz), T(nz), x_lo(nz), x_hi(nz)
      logical :: has_gap(nz), in_gap(nz)
      integer :: ierr

      call load_phase_diagram(pd, data_file, 0d0, .false., ierr)
      call assert_true('degenerate cells: phase diagram loaded', ierr == 0, nfail)
      if (ierr /= 0) return
      A = [A_h1, A_he4, A_c12, A_o16]
      P = 2d0*Mbar
      T = 1000d0
      xa(:,1) = [0d0, nearest(0.98d0, -1d0), 0.012d0, 0.008d0]   ! no H, sum = 1 - ulp
      xa(:,2) = [0d0, nearest(1d0, -1d0), 0d0, 0d0]              ! pure He, sum = 1 - ulp
      xa(:,3) = [0d0, 0d0, 0.6d0, 0.4d0]                         ! no H and no He
      xa(:,4) = [0.70d0, 0.27d0, 0.02d0, 0.01d0]                 ! normal cell for comparison
      call evaluate_miscibility(pd, nz, species, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('degenerate cells: no error', ierr == 0, nfail)
      call assert_true('degenerate cells: gap exists at (P, T)', all(has_gap), nfail)
      call assert_true('degenerate cells: no partner, sum = 1 - ulp, not in gap', .not. in_gap(1), nfail)
      call assert_true('degenerate cells: pure He, sum = 1 - ulp, not in gap', .not. in_gap(2), nfail)
      call assert_true('degenerate cells: no pair species, not in gap', .not. in_gap(3), nfail)
      call assert_true('degenerate cells: normal cell in gap', in_gap(4), nfail)
      call free_phase_diagram(pd)
   end subroutine test_degenerate_cells_not_in_gap


   subroutine test_errors(nfail)
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 6, species = 4
      real(dp) :: xa(species, nz), A(species), P(nz), T(nz), x_lo(nz), x_hi(nz), T_bad(nz), A_bad(species)
      logical :: has_gap(nz), in_gap(nz), short(nz-1)
      type(phase_diagram) :: pd, empty
      integer :: ierr, k

      write(*,'(a)') '   (error messages from evaluate_miscibility are expected below)'
      call load_phase_diagram(pd, data_file, 0d0, .false., ierr)
      if (ierr /= 0) return
      A = [A_h1, A_he4, A_c12, A_o16]
      do k = 1, nz
         xa(:,k) = [0.71d0, 0.27d0, 0.012d0, 0.008d0]
         P(k) = (1d0 + 0.5d0*k)*Mbar
         T(k) = 3000d0 + 500d0*k
      end do

      call evaluate_miscibility(pd, nz, species, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: baseline ok', ierr == 0 .and. all(in_gap), nfail)

      call evaluate_miscibility(pd, nz, species, xa, A, 2, 2, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: i_sed == i_partner', ierr == miscibility_err_bad_input, nfail)
      call evaluate_miscibility(pd, nz, species, xa, A, 0, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: i_sed out of range', ierr == miscibility_err_bad_input, nfail)
      call evaluate_miscibility(pd, nz, species, xa, A, 2, 5, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: i_partner out of range', ierr == miscibility_err_bad_input, nfail)
      call evaluate_miscibility(pd, nz, 1, xa, A, 1, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: single species', ierr == miscibility_err_bad_input, nfail)
      call evaluate_miscibility(pd, nz, species + 1, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: species larger than xa', ierr == miscibility_err_bad_input, nfail)
      call evaluate_miscibility(pd, nz + 1, species, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: nz larger than arrays', ierr == miscibility_err_bad_input, nfail)
      call evaluate_miscibility(pd, nz, species, xa, A, 2, 1, P, T, short, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: output array too short', ierr == miscibility_err_bad_input, nfail)
      A_bad = A
      A_bad(3) = 0d0
      call evaluate_miscibility(pd, nz, species, xa, A_bad, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: zero atomic weight', ierr == miscibility_err_bad_input, nfail)
      call evaluate_miscibility(empty, nz, species, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: unloaded phase diagram', ierr == phase_diagram_err_not_loaded, nfail)

      xa(3,3) = -1d-3
      call evaluate_miscibility(pd, nz, species, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: negative mass fraction', ierr == miscibility_err_composition, nfail)
      xa(3,3) = 0.012d0
      xa(1,2) = ieee_value(1d0, ieee_quiet_nan)
      call evaluate_miscibility(pd, nz, species, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: NaN mass fraction', ierr == miscibility_err_composition, nfail)
      xa(1,2) = 0.71d0

      T_bad = T
      T_bad(5) = ieee_value(1d0, ieee_quiet_nan)
      call evaluate_miscibility(pd, nz, species, xa, A, 2, 1, P, T_bad, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: NaN temperature', ierr == phase_diagram_err_bad_input, nfail)
      call assert_true('errors: other cells still evaluated', in_gap(1) .and. in_gap(6) .and. .not. has_gap(5), nfail)

      call evaluate_miscibility(pd, 0, species, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('errors: nz = 0 is trivially ok', ierr == 0, nfail)

      call assert_true('errors: message text', &
         trim(miscibility_error_message(miscibility_err_composition)) /= trim(miscibility_error_message(0)) &
         .and. trim(miscibility_error_message(phase_diagram_err_not_loaded)) == &
         trim(phase_diagram_error_message(phase_diagram_err_not_loaded)), nfail)

      call free_phase_diagram(pd)
   end subroutine test_errors


   subroutine test_large_profile(nfail)
      ! the OpenMP loop of evaluate_miscibility agrees with a serial evaluation
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 5000, species = 4
      real(dp), allocatable :: xa(:,:), P(:), T(:), x_lo(:), x_hi(:)
      logical, allocatable :: has_gap(:), in_gap(:)
      real(dp) :: A(species), y_lo, y_hi, x_lo_s, x_hi_s, Y, Z, max_err
      logical :: gap_s, ok
      type(phase_diagram) :: pd
      integer :: ierr, k, num_mismatch

      call load_phase_diagram(pd, data_file, 0d0, .false., ierr)
      if (ierr /= 0) return
      allocate(xa(species, nz), P(nz), T(nz), x_lo(nz), x_hi(nz), has_gap(nz), in_gap(nz))
      A = [A_h1, A_he4, A_c12, A_o16]
      do k = 1, nz
         P(k) = (0.3d0 + 25.7d0*dble(k-1)/dble(nz-1))*Mbar
         T(k) = 500d0 + 12500d0*(0.5d0 + 0.5d0*sin(0.37d0*dble(k)))
         Y = 0.05d0 + 0.92d0*(0.5d0 + 0.5d0*cos(0.11d0*dble(k)))  ! Y <= 0.97 so that X = 1 - Y - Z > 0
         Z = 0.02d0*(0.5d0 + 0.5d0*sin(0.05d0*dble(k)))
         xa(:,k) = [1d0 - Y - Z, Y, 0.6d0*Z, 0.4d0*Z]
      end do

      call evaluate_miscibility(pd, nz, species, xa, A, 2, 1, P, T, has_gap, in_gap, x_lo, x_hi, ierr)
      call assert_true('large profile: no error', ierr == 0, nfail)
      call assert_true('large profile: nontrivial', count(in_gap) > 0 .and. count(has_gap .and. .not. in_gap) > 0 &
         .and. count(.not. has_gap) > 0, nfail)

      num_mismatch = 0
      max_err = 0d0
      do k = 1, nz
         call phase_diagram_limits(pd, P(k), T(k), y_lo, y_hi, gap_s, ierr)
         if (ierr /= 0) num_mismatch = num_mismatch + 1
         if (gap_s) then
            x_lo_s = mass_fraction_for_pair_number_fraction(xa(:,k), A, 2, 1, y_lo)
            x_hi_s = mass_fraction_for_pair_number_fraction(xa(:,k), A, 2, 1, y_hi)
            ok = has_gap(k) .and. (in_gap(k) .eqv. (x_lo_s < xa(2,k) .and. xa(2,k) < x_hi_s))
            max_err = max(max_err, abs(x_lo_s - x_lo(k)), abs(x_hi_s - x_hi(k)))
         else
            ok = .not. has_gap(k) .and. .not. in_gap(k) .and. x_lo(k) == -1d0 .and. x_hi(k) == -1d0
         end if
         if (.not. ok) num_mismatch = num_mismatch + 1
      end do
      call assert_true('large profile: flags agree with serial evaluation', num_mismatch == 0, nfail)
      call assert_close('large profile: limits agree with serial evaluation', max_err, 0d0, 0d0, 1d-15, nfail)
      call free_phase_diagram(pd)
   end subroutine test_large_profile

end module test_miscibility
