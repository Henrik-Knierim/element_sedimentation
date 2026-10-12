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

! Unit tests for module advection_diffusion.

module test_advection_diffusion

   use const_def, only: dp, i8
   use advection_diffusion
   use test_support, only: section, assert_true, assert_close, assert_equal_int, next_random

   implicit none

   private

   public :: run_advection_diffusion_tests

contains

   subroutine run_advection_diffusion_tests(nfail)
      integer, intent(inout) :: nfail
      call section('advection_diffusion')
      call test_input_validation(nfail)
      call test_timescales_and_messages(nfail)
      call test_pure_diffusion(nfail)
      call test_upwind_advection(nfail)
      call test_upward_advection(nfail)
      call test_floor(nfail)
      call test_cap(nfail)
      call test_positivity(nfail)
      call test_flux_consistency(nfail)
      call test_passive_step(nfail)
      call test_passive_counterflow(nfail)
      call test_time_accurate(nfail)
      call test_change_limit(nfail)
      call test_failure_handling(nfail)
      call test_stiff_problems(nfail)
      call test_time_accurate_long(nfail)
      call test_capped_receiver_counterflow(nfail)
      call test_no_false_steady_state(nfail)
      call test_passive_species_in_steady_test(nfail)
   end subroutine run_advection_diffusion_tests


   real(dp) function total_mass(p, x)
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: x(:)
      integer :: k
      total_mass = 0d0
      do k = 1, p% nz
         total_mass = total_mass + p% dm(k)*x(k)
      end do
   end function total_mass


   real(dp) function quiet_nan()
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      quiet_nan = ieee_value(1d0, ieee_quiet_nan)
   end function quiet_nan


   subroutine make_random_problem(p, nz, seed, with_floor_cap, x)
      ! random non-uniform dm, sig, w (either sign) and, optionally, floors and caps
      type(adv_diff_problem), intent(inout) :: p
      integer, intent(in) :: nz
      integer(i8), intent(inout) :: seed
      logical, intent(in) :: with_floor_cap
      real(dp), intent(out) :: x(:)
      integer :: k, ierr
      call alloc_adv_diff_problem(p, nz, ierr)
      if (ierr /= 0) return
      do k = 1, nz
         p% dm(k) = 0.5d0 + next_random(seed)
         x(k) = next_random(seed)
      end do
      do k = 2, nz
         p% sig(k) = next_random(seed)
         p% w(k) = 2d0*next_random(seed) - 1d0
      end do
      if (with_floor_cap) then
         p% cap_width = 1d-2
         do k = 1, nz
            if (next_random(seed) < 0.5d0) then
               p% x_floor(k) = 0.3d0*next_random(seed)
            else
               p% x_floor(k) = adv_diff_huge
            end if
            if (next_random(seed) < 0.5d0) then
               p% x_cap(k) = 0.7d0 + 0.3d0*next_random(seed)
            else
               p% x_cap(k) = 1d0   ! a mass fraction can never exceed 1 (cap = 1 where there is no gap)
            end if
         end do
      end if
   end subroutine make_random_problem


   ! ---------------------------------------------------------------------------
   ! 1. input validation
   ! ---------------------------------------------------------------------------

   subroutine test_input_validation(nfail)
      integer, intent(inout) :: nfail
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(3), x4(4), xp(3,1), xp_new(3,1)
      character(len=256) :: msg
      integer :: ierr

      call alloc_adv_diff_problem(p, 0, ierr)
      call assert_true('alloc with nz = 0 rejected', ierr == adv_diff_err_bad_input, nfail)

      call alloc_adv_diff_problem(p, 3, ierr)
      call assert_true('alloc ok', ierr == 0 .and. p% nz == 3 .and. size(p% sig) == 4, nfail)
      call assert_true('alloc defaults', all(p% x_cap == adv_diff_huge) .and. all(p% scale == 1d0) &
         .and. all(p% x_floor == 0d0) .and. all(p% sig == 0d0) .and. all(p% w == 0d0), nfail)
      p% dm = 1d0
      x = 0.5d0
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('valid problem passes check', ierr == 0 .and. len_trim(msg) == 0, nfail)

      call adv_diff_integrate(p, x4, 0d0, opts, stats, ierr)
      call assert_true('wrong size of x -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      call assert_true('wrong size of x -> message', len_trim(stats% message) > 0, nfail)

      p% dm(2) = -1d0
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('negative dm -> bad_input (check)', ierr == adv_diff_err_bad_input .and. index(msg, 'dm') > 0, nfail)
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('negative dm -> bad_input (integrate)', ierr == adv_diff_err_bad_input, nfail)
      p% dm(2) = 1d0

      p% sig(2) = quiet_nan()
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('NaN sig -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      p% sig(2) = -1d0
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('negative sig -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      p% sig(2) = 0d0

      p% w(3) = quiet_nan()
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('NaN w -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      p% w(3) = 0d0

      p% w(1) = 1d0
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('nonzero flux at face 1 -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      p% w(1) = 0d0
      p% sig(4) = 1d0
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('nonzero flux at face nz+1 -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      p% sig(4) = 0d0

      p% cap_width = 0d0
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('cap_width = 0 -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      p% cap_width = 1d-3

      p% scale(1) = -1d0
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('negative scale -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      p% scale(1) = 1d0

      p% x_floor(2) = quiet_nan()
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('NaN x_floor -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      p% x_floor(2) = 0d0

      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('problem valid again', ierr == 0, nfail)

      p% sig(2) = 1d0
      x(2) = quiet_nan()
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('NaN in x -> bad_input', ierr == adv_diff_err_bad_input .and. stats% k_fail == 2, nfail)
      x(2) = -1d-3
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('x < -1d-12 -> bounds error', ierr == adv_diff_err_bounds .and. stats% k_fail == 2, nfail)
      x(2) = 0.5d0

      opts% tol_steady = 0d0
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('invalid option -> bad_input', ierr == adv_diff_err_bad_input &
         .and. index(stats% message, 'tol_steady') > 0, nfail)
      opts% tol_steady = 1d-8
      opts% tol_residual = -1d0
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('invalid tol_residual -> bad_input', ierr == adv_diff_err_bad_input &
         .and. index(stats% message, 'tol_residual') > 0, nfail)
      opts% tol_residual = 1d-8

      call adv_diff_integrate(p, x, quiet_nan(), opts, stats, ierr)
      call assert_true('NaN total_time -> bad_input', ierr == adv_diff_err_bad_input, nfail)

      call adv_diff_implicit_step(p, 0d0, x, x4(1:3), opts, stats, ierr)
      call assert_true('implicit step with dt = 0 -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      call adv_diff_implicit_step(p, 1d0, x4, x4(1:3), opts, stats, ierr)
      call assert_true('implicit step with wrong size -> bad_input', ierr == adv_diff_err_bad_input, nfail)

      xp(:,1) = 0.5d0
      call adv_diff_passive_step(p, 1d0, 0, xp, xp_new, ierr)
      call assert_true('passive step with nvar = 0 -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      call adv_diff_passive_step(p, -1d0, 1, xp, xp_new, ierr)
      call assert_true('passive step with dt < 0 -> bad_input', ierr == adv_diff_err_bad_input, nfail)
      call adv_diff_passive_step(p, 1d0, 1, xp, xp_new, ierr)
      call assert_true('passive step valid -> ok', ierr == 0, nfail)

      call free_adv_diff_problem(p)
      call assert_true('free resets nz', p% nz == 0 .and. .not. allocated(p% dm), nfail)
      call check_adv_diff_problem(p, ierr, msg)
      call assert_true('unallocated problem fails check', ierr == adv_diff_err_bad_input, nfail)
   end subroutine test_input_validation


   ! ---------------------------------------------------------------------------
   ! timescales, fluxes and error messages
   ! ---------------------------------------------------------------------------

   subroutine test_timescales_and_messages(nfail)
      integer, intent(inout) :: nfail
      type(adv_diff_problem) :: p
      real(dp) :: tau_min, tau_max, x(3), fa(4), fd(4)
      integer :: ierr, n, code
      logical :: ok

      call alloc_adv_diff_problem(p, 3, ierr)
      p% dm = [1d0, 2d0, 4d0]
      p% scale = [2d0, 1d0, 0.5d0]
      p% sig(2) = 0.5d0
      p% w(2) = 0.25d0      ! donor cell 1: rate = 0.5 + 0.25*2 = 1, tau = min(1,2)/1 = 1
      p% w(3) = -1d0        ! donor cell 3: rate = 0 + 1*0.5 = 0.5, tau = min(2,4)/0.5 = 4
      call adv_diff_timescales(p, tau_min, tau_max, n)
      call assert_close('timescales: tau_min', tau_min, 1d0, 1d-15, 0d0, nfail)
      call assert_close('timescales: tau_max', tau_max, 4d0, 1d-15, 0d0, nfail)
      call assert_equal_int('timescales: active faces', n, 2, nfail)

      ! fluxes: cell 1 above its floor, receiver 2 inside the throttle ramp, w < 0 at face 3
      p% x_floor = [0.1d0, 0d0, adv_diff_huge]
      p% x_cap = [adv_diff_huge, 0.5d0, adv_diff_huge]
      p% cap_width = 0.1d0
      x = [0.4d0, 0.45d0, 0.2d0]
      call adv_diff_fluxes(p, x, fa, fd)
      call assert_close('fluxes: diffusive face 2', fd(2), 0.5d0*(0.4d0 - 0.45d0), 1d-15, 0d0, nfail)
      ! w(2) scale(1) (x1 - floor1) throttle(x2) = 0.25*2*0.3*((0.5-0.45)/0.1)
      call assert_close('fluxes: advective face 2 (floor and throttle)', fa(2), 0.25d0*2d0*0.3d0*0.5d0, 1d-14, 0d0, nfail)
      call assert_true('fluxes: donor with huge floor does not advect', fa(3) == 0d0, nfail)
      call assert_true('fluxes: boundary faces are zero', fa(1) == 0d0 .and. fd(1) == 0d0 .and. fa(4) == 0d0 .and. fd(4) == 0d0, nfail)
      p% x_floor(3) = 0.1d0
      call adv_diff_fluxes(p, x, fa, fd)
      ! w(3) scale(3) (x3 - floor3) throttle(x2) = -1*0.5*0.1*0.5, receiver 2 is in its ramp
      call assert_close('fluxes: upward advection (w < 0) is negative', fa(3), -1d0*0.5d0*0.1d0*0.5d0, 1d-14, 0d0, nfail)
      x(2) = 0.6d0
      call adv_diff_fluxes(p, x, fa, fd)
      call assert_true('fluxes: receiver above cap gets nothing', fa(2) == 0d0, nfail)
      p% w = 0d0
      p% sig = 0d0
      call adv_diff_timescales(p, tau_min, tau_max, n)
      call assert_true('timescales: no active faces', n == 0 .and. tau_min == 0d0 .and. tau_max == 0d0, nfail)
      call free_adv_diff_problem(p)

      ok = .true.
      do code = 0, 8
         ok = ok .and. len_trim(adv_diff_error_message(code)) > 0
      end do
      call assert_true('error messages: non-empty for all codes', ok, nfail)
      call assert_true('error messages: unknown code', index(adv_diff_error_message(99), 'unknown') > 0, nfail)
      call assert_true('error messages: ok', trim(adv_diff_error_message(adv_diff_ok)) == 'ok', nfail)
   end subroutine test_timescales_and_messages


   ! ---------------------------------------------------------------------------
   ! 2. pure diffusion with zero-flux boundaries on a non-uniform grid
   ! ---------------------------------------------------------------------------

   subroutine test_pure_diffusion(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 20
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), x0(nz), xbar, mass0
      integer(i8) :: seed
      integer :: ierr, k

      seed = 4242_i8
      call make_random_problem(p, nz, seed, .false., x)
      p% w = 0d0
      do k = 2, nz
         p% sig(k) = 0.1d0 + 2d0*next_random(seed)
      end do
      x0 = x
      mass0 = total_mass(p, x0)
      xbar = mass0/sum(p% dm)

      opts% tol_steady = 1d-14
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('diffusion: no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_equal_int('diffusion: stop reason steady state', stats% stop_reason, adv_diff_stop_steady_state, nfail)
      call assert_close('diffusion: max |x - mean|', maxval(abs(x - xbar)), 0d0, 0d0, 1d-10, nfail)
      call assert_close('diffusion: mass conserved', total_mass(p, x), mass0, 1d-13, 0d0, nfail)
      call assert_true('diffusion: few steps', stats% steps > 0 .and. stats% steps <= 60, nfail)
      ! a linear step needs one solve plus at most one refinement solve when dt sig >> dm
      call assert_true('diffusion: at most two Newton solves per linear step', &
         stats% newton_iters <= 2*stats% steps, nfail)
      call assert_true('diffusion: stats consistent', stats% residual <= opts% tol_steady &
         .and. stats% tau_min > 0d0 .and. stats% tau_max >= stats% tau_min .and. stats% rejected_steps == 0 &
         .and. stats% time > 0d0 .and. stats% conservation_error <= opts% conservation_rtol, nfail)
      write(*,'(a,i6,a,1pd12.4,a,1pd12.4)') '          steps = ', stats% steps, '  residual = ', stats% residual, &
         '  conservation error = ', stats% conservation_error
      call free_adv_diff_problem(p)
   end subroutine test_pure_diffusion


   ! ---------------------------------------------------------------------------
   ! 3. pure upwind advection: single step exactness and geometric steady state
   ! ---------------------------------------------------------------------------

   subroutine test_upwind_advection(nfail)
      integer, intent(inout) :: nfail
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x_old(3), x_new(3), x_hand(3), dt, w2, w3
      real(dp) :: x(6), mass0, ratio
      integer :: ierr, k
      logical :: ok

      ! single step, nz = 3, hand-computed tridiagonal solution
      call alloc_adv_diff_problem(p, 3, ierr)
      p% dm = [1d0, 2d0, 3d0]
      w2 = 0.7d0
      w3 = 1.3d0
      p% w(2) = w2
      p% w(3) = w3
      x_old = [0.5d0, 0.2d0, 0.1d0]
      dt = 0.9d0
      x_hand(1) = p% dm(1)*x_old(1)/(p% dm(1) + dt*w2)
      x_hand(2) = (p% dm(2)*x_old(2) + dt*w2*x_hand(1))/(p% dm(2) + dt*w3)
      x_hand(3) = x_old(3) + dt*w3*x_hand(2)/p% dm(3)
      call adv_diff_implicit_step(p, dt, x_old, x_new, opts, stats, ierr)
      call assert_true('advection step: no error', ierr == 0, nfail)
      call assert_close('advection step: x(1)', x_new(1), x_hand(1), 1d-14, 0d0, nfail)
      call assert_close('advection step: x(2)', x_new(2), x_hand(2), 1d-14, 0d0, nfail)
      call assert_close('advection step: x(3)', x_new(3), x_hand(3), 1d-14, 0d0, nfail)
      call assert_close('advection step: mass conserved', total_mass(p, x_new), total_mass(p, x_old), 1d-15, 0d0, nfail)
      call assert_equal_int('advection step: one Newton solve for a linear step', stats% newton_iters, 1, nfail)
      call free_adv_diff_problem(p)

      ! steady state for constant sig and w on a uniform grid: x_k/x_{k-1} = 1 + w/sig
      call alloc_adv_diff_problem(p, 6, ierr)
      p% dm = 1d0
      p% sig(2:6) = 1d0
      p% w(2:6) = 2d0
      x = 0.1d0
      mass0 = total_mass(p, x)
      opts% tol_steady = 1d-14
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('advection steady: no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      ok = .true.
      do k = 2, 6
         ratio = x(k)/x(k-1)
         ok = ok .and. abs(ratio - 3d0) <= 1d-10*3d0
      end do
      call assert_true('advection steady: geometric profile x_k/x_{k-1} = 1 + w/sig', ok, nfail)
      call assert_close('advection steady: x(1) analytic', x(1), mass0/364d0, 1d-10, 0d0, nfail)
      call assert_close('advection steady: mass conserved', total_mass(p, x), mass0, 1d-13, 0d0, nfail)
      call assert_equal_int('advection steady: stop reason', stats% stop_reason, adv_diff_stop_steady_state, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_upwind_advection


   subroutine test_upward_advection(nfail)
      ! same as the geometric steady state but with w < 0 (advection toward smaller k)
      integer, intent(inout) :: nfail
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(6), mass0, ratio
      integer :: ierr, k
      logical :: ok

      call alloc_adv_diff_problem(p, 6, ierr)
      p% dm = 1d0
      p% sig(2:6) = 1d0
      p% w(2:6) = -2d0
      x = 0.1d0
      mass0 = total_mass(p, x)
      opts% tol_steady = 1d-14
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('upward advection: no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      ok = .true.
      do k = 2, 6
         ratio = x(k-1)/x(k)
         ok = ok .and. abs(ratio - 3d0) <= 1d-10*3d0
      end do
      call assert_true('upward advection: geometric profile x_{k-1}/x_k = 1 + |w|/sig', ok, nfail)
      call assert_close('upward advection: x(6) analytic', x(6), mass0/364d0, 1d-10, 0d0, nfail)
      call assert_close('upward advection: mass conserved', total_mass(p, x), mass0, 1d-13, 0d0, nfail)
      call assert_true('upward advection: all positive', all(x > 0d0), nfail)
      call free_adv_diff_problem(p)
   end subroutine test_upward_advection


   ! ---------------------------------------------------------------------------
   ! 4. floor: envelope (huge sig) / rain region (floor) / deposit region (sig, huge floor)
   ! ---------------------------------------------------------------------------

   subroutine test_floor(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 12, k_rain_lo = 5, k_rain_hi = 8
      real(dp), parameter :: floor = 0.2d0
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), mass0, x_dep, m_above
      integer :: ierr, k

      call alloc_adv_diff_problem(p, nz, ierr)
      do k = 1, nz
         p% dm(k) = 1d0 + 0.1d0*k
      end do
      p% x_floor = adv_diff_huge
      ! envelope cells 1..4: huge diffusion, connected by diffusion to the top rain cell
      p% sig(2:k_rain_lo) = 1d6
      ! rain region cells 5..8: advection only, with a floor; rain out of cell 8 into the deposit
      p% x_floor(k_rain_lo:k_rain_hi) = floor
      p% w(k_rain_lo+1:k_rain_hi+1) = 1d0
      ! deposit cells 9..12: diffusion, never rain
      p% sig(k_rain_hi+2:nz) = 1d0
      x(1:k_rain_hi) = 0.5d0
      x(k_rain_hi+1:nz) = 0.1d0
      mass0 = total_mass(p, x)
      m_above = sum(p% dm(1:k_rain_hi))
      x_dep = (mass0 - floor*m_above)/sum(p% dm(k_rain_hi+1:nz))

      opts% tol_steady = 1d-14
      opts% x_upper_bound = 1d0
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('floor: no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_equal_int('floor: stop reason steady state', stats% stop_reason, adv_diff_stop_steady_state, nfail)
      call assert_close('floor: rain cells drained to the floor', maxval(x(k_rain_lo:k_rain_hi)) - floor, 0d0, 0d0, 1d-8, nfail)
      call assert_true('floor: rain cells never below the floor', all(x(k_rain_lo:k_rain_hi) >= floor - 1d-12), nfail)
      call assert_close('floor: envelope depleted to the top rain cell', maxval(abs(x(1:k_rain_lo-1) - x(k_rain_lo))), &
         0d0, 0d0, 1d-10, nfail)
      call assert_close('floor: deposit homogeneous', maxval(x(k_rain_hi+1:nz)) - minval(x(k_rain_hi+1:nz)), 0d0, 0d0, 1d-10, nfail)
      call assert_close('floor: deposit value from mass budget', x(nz), x_dep, 1d-8, 0d0, nfail)
      call assert_close('floor: mass conserved', total_mass(p, x), mass0, 1d-13, 0d0, nfail)
      call assert_true('floor: all x in [0,1]', all(x >= 0d0) .and. all(x <= 1d0), nfail)
      write(*,'(a,i6,a,i6,a,i6)') '          steps = ', stats% steps, '  newton_iters = ', stats% newton_iters, &
         '  rejected = ', stats% rejected_steps
      call free_adv_diff_problem(p)
   end subroutine test_floor


   ! ---------------------------------------------------------------------------
   ! 5. cap: pure advection into a bottom cell with cap 0.9
   ! ---------------------------------------------------------------------------

   subroutine test_cap(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 5
      real(dp), parameter :: cap = 0.9d0
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), mass0
      integer :: ierr

      call alloc_adv_diff_problem(p, nz, ierr)
      p% dm = 1d0
      p% w(2:nz) = 1d0
      p% x_cap(nz) = cap
      x = 0.3d0
      mass0 = total_mass(p, x)
      opts% tol_steady = 1d-14
      opts% x_upper_bound = 1d0
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('cap: no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_true('cap: bottom cell does not exceed the cap', x(nz) <= cap + 1d-12, nfail)
      call assert_close('cap: bottom cell saturated', x(nz), cap, 0d0, 1d-6, nfail)
      call assert_close('cap: cell above fills with the rest', x(nz-1), mass0 - cap, 1d-6, 0d0, nfail)
      call assert_close('cap: upper cells drained', maxval(x(1:nz-2)), 0d0, 0d0, 1d-8, nfail)
      call assert_close('cap: mass conserved', total_mass(p, x), mass0, 1d-13, 0d0, nfail)
      call assert_true('cap: no x > 1 and no x < 0', all(x <= 1d0) .and. all(x >= 0d0), nfail)
      call assert_equal_int('cap: stop reason steady state', stats% stop_reason, adv_diff_stop_steady_state, nfail)
      write(*,'(a,i6,a,i6,a,i6)') '          steps = ', stats% steps, '  newton_iters = ', stats% newton_iters, &
         '  rejected = ', stats% rejected_steps

      ! a second cap above the first: the layer grows upward cell by cell
      x = 0.3d0
      p% x_cap(nz-1) = 0.5d0
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('cap: two capped cells, no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_true('cap: two capped cells respect their caps', x(nz) <= cap + 1d-12 .and. x(nz-1) <= 0.5d0 + 1d-12, nfail)
      call assert_close('cap: two capped cells, third fills', x(nz-2), mass0 - cap - 0.5d0, 1d-6, 0d0, nfail)
      call assert_close('cap: two capped cells, mass conserved', total_mass(p, x), mass0, 1d-13, 0d0, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_cap


   ! ---------------------------------------------------------------------------
   ! 6. positivity and bounds with random coefficients, floors and caps
   ! ---------------------------------------------------------------------------

   subroutine test_positivity(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 30
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), x0(nz), mass0, tau_min, tau_max
      integer(i8) :: seed
      integer :: ierr, n

      seed = 777_i8
      call make_random_problem(p, nz, seed, .true., x)
      x0 = x
      mass0 = total_mass(p, x0)
      opts% x_upper_bound = 1d0
      opts% dt_growth_max = 1.3d0
      opts% max_dt_factor = 1d3

      ! many time-accurate steps
      call adv_diff_timescales(p, tau_min, tau_max, n)
      call adv_diff_integrate(p, x, 200d0*tau_max, opts, stats, ierr)
      call assert_true('positivity (time-accurate): no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_true('positivity (time-accurate): many steps', stats% steps >= 20, nfail)
      call assert_true('positivity (time-accurate): x in [0,1]', all(x >= 0d0) .and. all(x <= 1d0), nfail)
      call assert_close('positivity (time-accurate): mass conserved', total_mass(p, x), mass0, 1d-12, 0d0, nfail)
      call assert_equal_int('positivity (time-accurate): stop reason', stats% stop_reason, adv_diff_stop_total_time, nfail)
      write(*,'(a,i6,a,i6,a,i6,a,i6)') '          steps = ', stats% steps, '  newton_iters = ', stats% newton_iters, &
         '  rejected = ', stats% rejected_steps, '  clipped = ', stats% num_clipped

      ! steady state of the same problem
      x = x0
      opts% dt_growth_max = 10d0
      opts% max_dt_factor = 1d4
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('positivity (steady): no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_true('positivity (steady): x in [0,1]', all(x >= 0d0) .and. all(x <= 1d0), nfail)
      call assert_close('positivity (steady): mass conserved', total_mass(p, x), mass0, 1d-12, 0d0, nfail)
      call assert_equal_int('positivity (steady): stop reason', stats% stop_reason, adv_diff_stop_steady_state, nfail)
      write(*,'(a,i6,a,i6,a,i6,a,i6)') '          steps = ', stats% steps, '  newton_iters = ', stats% newton_iters, &
         '  rejected = ', stats% rejected_steps, '  clipped = ', stats% num_clipped
      call free_adv_diff_problem(p)
   end subroutine test_positivity


   ! ---------------------------------------------------------------------------
   ! fluxes after an implicit step are consistent with (x_new - x_old) dm / dt
   ! ---------------------------------------------------------------------------

   subroutine test_flux_consistency(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 25
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x_old(nz), x_new(nz), fa(nz+1), fd(nz+1), tau_min, tau_max, dt, lhs, rhs, worst, scale
      integer(i8) :: seed
      integer :: ierr, k, n

      seed = 99_i8
      call make_random_problem(p, nz, seed, .true., x_old)
      call adv_diff_timescales(p, tau_min, tau_max, n)
      dt = 2d0*tau_max
      call adv_diff_implicit_step(p, dt, x_old, x_new, opts, stats, ierr)
      call assert_true('flux consistency: step ok', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call adv_diff_fluxes(p, x_new, fa, fd)
      worst = 0d0
      do k = 1, nz
         lhs = (x_new(k) - x_old(k))*p% dm(k)/dt
         rhs = (fa(k) + fd(k)) - (fa(k+1) + fd(k+1))
         scale = max(abs(lhs), abs(fa(k)) + abs(fd(k)) + abs(fa(k+1)) + abs(fd(k+1)), p% dm(k)/dt)
         worst = max(worst, abs(lhs - rhs)/scale)
      end do
      call assert_close('flux consistency: (x_new - x_old) dm/dt = F_k - F_k+1', worst, 0d0, 0d0, 1d-12, nfail)
      call assert_true('flux consistency: Newton iterations counted', stats% newton_iters >= 1, nfail)
      call assert_close('flux consistency: mass conserved', total_mass(p, x_new), total_mass(p, x_old), 1d-13, 0d0, nfail)
      write(*,'(a,i6,a,1pd12.4)') '          newton_iters = ', stats% newton_iters, '  worst relative mismatch = ', worst
      call free_adv_diff_problem(p)
   end subroutine test_flux_consistency


   ! ---------------------------------------------------------------------------
   ! passive step: linearity, positivity, several right-hand sides
   ! ---------------------------------------------------------------------------

   subroutine test_passive_step(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 15, nvar = 3
      type(adv_diff_problem) :: p
      real(dp) :: x_old(nz, nvar), x_new(nz, nvar), one_old(nz, 1), one_new(nz, 1), xs(nz)
      real(dp) :: tau_min, tau_max, dt, worst
      integer(i8) :: seed
      integer :: ierr, j, k, n
      logical :: ok

      seed = 31_i8
      call make_random_problem(p, nz, seed, .true., xs)
      do k = 1, nz
         p% scale(k) = 0.5d0 + next_random(seed)
      end do
      do j = 1, nvar
         do k = 1, nz
            x_old(k, j) = next_random(seed)
         end do
      end do
      call adv_diff_timescales(p, tau_min, tau_max, n)
      dt = 5d0*tau_max
      call adv_diff_passive_step(p, dt, nvar, x_old, x_new, ierr)
      call assert_true('passive step: no error', ierr == 0, nfail)
      call assert_true('passive step: positivity', all(x_new >= 0d0), nfail)
      ok = .true.
      do j = 1, nvar
         ok = ok .and. abs(total_mass(p, x_new(:, j)) - total_mass(p, x_old(:, j))) <= 1d-13*total_mass(p, x_old(:, j))
      end do
      call assert_true('passive step: mass of every variable conserved', ok, nfail)
      ! several right-hand sides give the same answer as separate solves
      worst = 0d0
      do j = 1, nvar
         one_old(:, 1) = x_old(:, j)
         call adv_diff_passive_step(p, dt, 1, one_old, one_new, ierr)
         worst = max(worst, maxval(abs(one_new(:, 1) - x_new(:, j))))
      end do
      call assert_close('passive step: multi-rhs equals single-rhs solves', worst, 0d0, 0d0, 1d-15, nfail)
      ! linearity: the solve of the sum is the sum of the solves
      one_old(:, 1) = x_old(:, 1) + x_old(:, 2)
      call adv_diff_passive_step(p, dt, 1, one_old, one_new, ierr)
      call assert_close('passive step: linear', maxval(abs(one_new(:, 1) - x_new(:, 1) - x_new(:, 2))), 0d0, 0d0, 1d-14, nfail)
      ! dt = 0 is the identity
      call adv_diff_passive_step(p, 0d0, nvar, x_old, x_new, ierr)
      call assert_true('passive step: dt = 0 is the identity', ierr == 0 .and. maxval(abs(x_new - x_old)) <= 1d-15, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_passive_step


   ! ---------------------------------------------------------------------------
   ! 7. passive counter-flow: sum of all species stays 1, every species conserved
   ! ---------------------------------------------------------------------------

   subroutine test_passive_counterflow(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 10, npas = 3
      type(adv_diff_problem) :: p, p_pas
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), xp(nz, npas), xp_new(nz, npas), mass_pas0(npas), mass_act0, r1, r2, r3
      real(dp) :: worst_sum, worst_mass, min_pas, tau_min, tau_max
      integer(i8) :: seed
      integer :: ierr, k, j, n, num_calls

      seed = 2024_i8
      call make_random_problem(p, nz, seed, .false., x)
      p% cap_width = 1d-2
      p% w(2:nz) = 0.5d0
      p% x_floor = adv_diff_huge
      p% x_floor(2:5) = 0.1d0
      p% x_cap(8:nz) = 0.8d0
      do k = 1, nz
         x(k) = 0.2d0 + 0.3d0*next_random(seed)
         r1 = next_random(seed)
         r2 = next_random(seed)
         r3 = next_random(seed)
         xp(k, 1) = (1d0 - x(k))*r1/(r1 + r2 + r3)
         xp(k, 2) = (1d0 - x(k))*r2/(r1 + r2 + r3)
         xp(k, 3) = (1d0 - x(k))*r3/(r1 + r2 + r3)
      end do
      mass_act0 = total_mass(p, x)
      do j = 1, npas
         mass_pas0(j) = total_mass(p, xp(:, j))
      end do

      call alloc_adv_diff_problem(p_pas, nz, ierr)
      p_pas% dm = p% dm
      p_pas% sig = p% sig

      worst_sum = 0d0
      worst_mass = 0d0
      min_pas = huge(1d0)
      num_calls = 0
      call adv_diff_timescales(p, tau_min, tau_max, n)
      opts% x_upper_bound = 1d0
      opts% dt_growth_max = 2d0
      opts% max_dt_factor = 10d0
      call adv_diff_integrate(p, x, 50d0*tau_max, opts, stats, ierr, after_step_passive)
      call assert_true('counter-flow: no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_true('counter-flow: callback called once per accepted step', num_calls == stats% steps .and. num_calls > 5, nfail)
      call assert_close('counter-flow: max |sum_j x_j - 1| over all steps and cells', worst_sum, 0d0, 0d0, 1d-13, nfail)
      call assert_close('counter-flow: max relative mass change of any species', worst_mass, 0d0, 0d0, 1d-12, nfail)
      call assert_true('counter-flow: passives stay >= 0', min_pas >= 0d0, nfail)
      call assert_close('counter-flow: active species mass conserved', total_mass(p, x), mass_act0, 1d-13, 0d0, nfail)
      call assert_true('counter-flow: active species moved', stats% max_change_total > 1d-3, nfail)
      write(*,'(a,i6,a,1pd12.4,a,1pd12.4)') '          steps = ', stats% steps, '  worst |sum - 1| = ', worst_sum, &
         '  worst mass change = ', worst_mass
      call free_adv_diff_problem(p)
      call free_adv_diff_problem(p_pas)

      contains

      subroutine after_step_passive(dt, x_new, flux_adv, flux_diff, change_other, ierr)
         ! the counter-flow of the passive species, as element_sedimentation does it
         real(dp), intent(in) :: dt
         real(dp), intent(in) :: x_new(:), flux_adv(:), flux_diff(:)
         real(dp), intent(out) :: change_other
         integer, intent(out) :: ierr
         integer :: kk, jj
         real(dp) :: s
         ierr = 0
         change_other = 0d0
         num_calls = num_calls + 1
         p_pas% w(1:nz+1) = -flux_adv(1:nz+1)
         do kk = 1, nz
            p_pas% scale(kk) = 1d0/max(1d0 - x_new(kk), 1d-30)
         end do
         call adv_diff_passive_step(p_pas, dt, npas, xp, xp_new, ierr)
         if (ierr /= 0) return
         change_other = maxval(abs(xp_new - xp))
         xp = xp_new
         do kk = 1, nz
            s = x_new(kk)
            do jj = 1, npas
               s = s + xp(kk, jj)
               min_pas = min(min_pas, xp(kk, jj))
            end do
            worst_sum = max(worst_sum, abs(s - 1d0))
         end do
         do jj = 1, npas
            worst_mass = max(worst_mass, abs(total_mass(p_pas, xp(:, jj)) - mass_pas0(jj))/mass_pas0(jj))
         end do
         if (size(flux_diff) /= nz+1) ierr = -1
      end subroutine after_step_passive

   end subroutine test_passive_counterflow


   ! ---------------------------------------------------------------------------
   ! 8. time-accurate mode against the analytic two-cell exchange
   ! ---------------------------------------------------------------------------

   subroutine test_time_accurate(nfail)
      integer, intent(inout) :: nfail
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(2), xbar, tau, t_end, x1_exact, err1, err2, dt
      integer :: ierr

      call alloc_adv_diff_problem(p, 2, ierr)
      p% dm = [1d0, 3d0]
      p% sig(2) = 0.5d0
      tau = p% dm(1)*p% dm(2)/((p% dm(1) + p% dm(2))*p% sig(2))
      x = [1d0, 0d0]
      xbar = total_mass(p, x)/sum(p% dm)
      t_end = 2d0
      x1_exact = xbar + (1d0 - xbar)*exp(-t_end/tau)

      opts% dt_growth_max = 1d0        ! constant dt
      opts% dt_initial_factor = 0.01d0  ! dt = 0.01 tau_min = 0.02
      opts% max_dt_factor = 1d6
      call adv_diff_integrate(p, x, t_end, opts, stats, ierr)
      dt = opts% dt_initial_factor*stats% tau_min
      err1 = abs(x(1) - x1_exact)
      call assert_true('time-accurate: no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_equal_int('time-accurate: stop reason total_time', stats% stop_reason, adv_diff_stop_total_time, nfail)
      call assert_true('time-accurate: stats% time == total_time exactly', stats% time == t_end, nfail)
      call assert_equal_int('time-accurate: number of steps', stats% steps, nint(t_end/dt), nfail)
      call assert_true('time-accurate: backward Euler within first-order error bound', err1 <= (1d0 - xbar)*dt/tau, nfail)
      call assert_true('time-accurate: error is nonzero (not a trivial pass)', err1 > 1d-8, nfail)
      call assert_close('time-accurate: mass conserved', total_mass(p, x), xbar*sum(p% dm), 1d-14, 0d0, nfail)

      ! halving dt roughly halves the error (first order)
      x = [1d0, 0d0]
      opts% dt_initial_factor = 0.005d0
      call adv_diff_integrate(p, x, t_end, opts, stats, ierr)
      err2 = abs(x(1) - x1_exact)
      call assert_true('time-accurate: first-order convergence in dt', ierr == 0 .and. err2 <= 0.6d0*err1 .and. err2 >= 0.4d0*err1, nfail)
      write(*,'(a,1pd12.4,a,1pd12.4,a,1pd12.4)') '          error(dt) = ', err1, '  error(dt/2) = ', err2, '  ratio = ', err2/err1

      ! the end time is hit exactly also when the last step is shortened
      x = [1d0, 0d0]
      opts% dt_initial_factor = 0.3d0
      opts% dt_growth_max = 10d0
      call adv_diff_integrate(p, x, t_end, opts, stats, ierr)
      call assert_true('time-accurate (growing dt): ends exactly at total_time', ierr == 0 &
         .and. stats% time == t_end .and. stats% stop_reason == adv_diff_stop_total_time, nfail)
      call assert_true('time-accurate (growing dt): last dt not larger than total_time', stats% dt_last <= t_end, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_time_accurate


   ! ---------------------------------------------------------------------------
   ! 9. max_abs_change stops the integration with the bound honoured
   ! ---------------------------------------------------------------------------

   subroutine test_change_limit(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 20
      real(dp), parameter :: limit = 0.05d0
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), x0(nz), change
      integer(i8) :: seed
      integer :: ierr

      seed = 4242_i8
      call make_random_problem(p, nz, seed, .false., x)
      x0 = x
      opts% max_abs_change = limit

      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      change = maxval(abs(x - x0))
      call assert_true('change limit (steady mode): no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_equal_int('change limit (steady mode): stop reason', stats% stop_reason, adv_diff_stop_change_limit, nfail)
      call assert_true('change limit (steady mode): bound honoured', change <= limit*(1d0 + 1d-12), nfail)
      call assert_true('change limit (steady mode): limit actually reached', change >= 0.9d0*limit, nfail)
      call assert_close('change limit (steady mode): stats% max_change_total', stats% max_change_total, change, 1d-15, 0d0, nfail)
      call assert_close('change limit (steady mode): mass conserved', total_mass(p, x), total_mass(p, x0), 1d-13, 0d0, nfail)

      x = x0
      call adv_diff_integrate(p, x, 1d3*stats% tau_max, opts, stats, ierr)
      change = maxval(abs(x - x0))
      call assert_true('change limit (time mode): no error and stop reason', ierr == 0 &
         .and. stats% stop_reason == adv_diff_stop_change_limit, nfail)
      call assert_true('change limit (time mode): bound honoured and reached', change <= limit*(1d0 + 1d-12) &
         .and. change >= 0.9d0*limit, nfail)

      x = x0
      opts% max_abs_change = 10d0
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('change limit never reached: steady state instead', ierr == 0 &
         .and. stats% stop_reason == adv_diff_stop_steady_state, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_change_limit


   ! ---------------------------------------------------------------------------
   ! 10. failure handling: clean error codes and messages, no crash
   ! ---------------------------------------------------------------------------

   subroutine test_failure_handling(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 5
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), x_new(nz), x0(nz)
      integer :: ierr

      ! nonlinear cap problem with a single Newton iteration allowed
      call alloc_adv_diff_problem(p, nz, ierr)
      p% dm = 1d0
      p% w(2:nz) = 1d0
      p% x_cap(nz) = 0.9d0
      x = 0.3d0
      x0 = x
      opts% max_newton_iters = 1
      call adv_diff_implicit_step(p, 1d3, x, x_new, opts, stats, ierr)
      call assert_true('newton limit (step): converged or clean Newton error', ierr == 0 .or. ierr == adv_diff_err_newton, nfail)
      call assert_true('newton limit (step): message when failing', ierr == 0 .or. len_trim(stats% message) > 0, nfail)
      write(*,'(a,i3,a,a)') '          ierr = ', ierr, '  ', trim(stats% message)
      ! with a single Newton iteration per step only tiny (effectively linear) steps succeed, so the
      ! continuation either crawls to max_steps or gives up on dt; all of these must be clean errors
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('newton limit (integrate): converged or clean error', ierr == 0 .or. ierr == adv_diff_err_newton &
         .or. ierr == adv_diff_err_dt_too_small .or. ierr == adv_diff_err_no_steady_state, nfail)
      call assert_true('newton limit (integrate): message when failing', ierr == 0 .or. len_trim(stats% message) > 0, nfail)
      call assert_true('newton limit (integrate): rejected steps counted when failing', ierr == 0 .or. stats% rejected_steps > 0, nfail)
      call assert_true('newton limit (integrate): state stays in bounds', all(x >= 0d0) .and. all(x <= 1d0), nfail)
      call assert_close('newton limit (integrate): mass conserved', total_mass(p, x), total_mass(p, x0), 1d-13, 0d0, nfail)
      write(*,'(a,i3,a,a)') '          ierr = ', ierr, '  ', trim(stats% message)
      opts% max_newton_iters = 50

      ! steady state not reached within max_steps
      x = x0
      opts% max_steps = 1
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_equal_int('max_steps = 1: no_steady_state error', ierr, adv_diff_err_no_steady_state, nfail)
      call assert_true('max_steps = 1: message', len_trim(stats% message) > 0 .and. index(stats% message, 'max_steps') > 0, nfail)
      call assert_true('max_steps = 1: one step taken, no stop reason', stats% steps == 1 .and. stats% stop_reason == 0, nfail)
      write(*,'(a,a)') '          ', trim(stats% message)
      x = x0
      call adv_diff_integrate(p, x, 1d3, opts, stats, ierr)
      call assert_true('max_steps = 1 (time mode): clean error', ierr == adv_diff_err_no_steady_state &
         .and. index(stats% message, 'total_time') > 0, nfail)
      opts% max_steps = 1000

      ! callback failure
      x = x0
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr, failing_callback)
      call assert_equal_int('failing callback: error code', ierr, adv_diff_err_callback, nfail)
      call assert_true('failing callback: message and one step', len_trim(stats% message) > 0 .and. stats% steps == 1, nfail)

      ! nothing to do
      p% w = 0d0
      x = x0
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('no active face: nothing_to_do', ierr == 0 .and. stats% stop_reason == adv_diff_stop_nothing_to_do &
         .and. stats% steps == 0 .and. all(x == x0), nfail)
      call adv_diff_integrate(p, x, 5d0, opts, stats, ierr)
      call assert_true('no active face (time mode): time set', ierr == 0 .and. stats% time == 5d0, nfail)
      call free_adv_diff_problem(p)

      contains

      subroutine failing_callback(dt, x_new, flux_adv, flux_diff, change_other, ierr)
         real(dp), intent(in) :: dt
         real(dp), intent(in) :: x_new(:), flux_adv(:), flux_diff(:)
         real(dp), intent(out) :: change_other
         integer, intent(out) :: ierr
         change_other = 0d0
         ierr = 7
         if (dt < 0d0 .or. size(x_new) < 1 .or. size(flux_adv) < 1 .or. size(flux_diff) < 1) ierr = 8
      end subroutine failing_callback

   end subroutine test_failure_handling



   ! ---------------------------------------------------------------------------
   ! 11. very stiff problems: dt*sig/dm far beyond 1/epsilon (pivots must stay exact)
   ! ---------------------------------------------------------------------------

   subroutine test_stiff_problems(nfail)
      integer, intent(inout) :: nfail
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(3), x_new(3), x0(3), mass0
      integer :: ierr, i
      real(dp), parameter :: weak(3) = [1d-13, 1d-30, 1d-60]

      ! one implicit step with dt*sig/dm = 1d17: the classic Thomas recurrence cancels to a
      ! zero pivot here; the column-sum recurrence must give the equilibrated state
      call alloc_adv_diff_problem(p, 3, ierr)
      p% dm = [1d0, 1d0, 2d0]
      p% sig(2) = 1d0
      p% sig(3) = 1d0
      x0 = [1d0, 0d0, 0.5d0]
      mass0 = total_mass(p, x0)
      call adv_diff_implicit_step(p, 1d17, x0, x_new, opts, stats, ierr)
      call assert_true('stiff step (dt*sig/dm = 1d17): no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_close('stiff step: equilibrated', maxval(x_new) - minval(x_new), 0d0, 0d0, 1d-14, nfail)
      call assert_close('stiff step: mass conserved', total_mass(p, x_new), mass0, 1d-15, 0d0, nfail)

      ! steady state with one fast and one very slow face (tau_max/tau_min up to 1d60)
      do i = 1, size(weak)
         p% sig(3) = weak(i)
         x = [1d0, 0d0, 0d0]   ! mass has to cross the weak face
         mass0 = total_mass(p, x)
         call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
         call assert_true('stiff steady state: no error', ierr == 0, nfail)
         if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
         call assert_close('stiff steady state: uniform', maxval(x) - minval(x), 0d0, 0d0, 1d-6, nfail)
         ! dt must reach the slow face's timescale; for sig(3) <= 1d-30 that means dt*sig/dm >> 1d17 at
         ! the fast face, where the cancelling Thomas recurrence would give a zero pivot
         call assert_true('stiff steady state: dt grew past the slowest timescale', &
            stats% dt_last >= stats% tau_max, nfail)
         call assert_close('stiff steady state: mass conserved', total_mass(p, x), mass0, 1d-13, 0d0, nfail)
         write(*,'(a,1pd9.2,a,i5,a,i4)') '          sig(3) = ', weak(i), '  steps = ', stats% steps, &
            '  rejected = ', stats% rejected_steps
      end do
      call free_adv_diff_problem(p)
   end subroutine test_stiff_problems


   ! ---------------------------------------------------------------------------
   ! 12. time-accurate mode over many face timescales (a MESA time step) within max_steps
   ! ---------------------------------------------------------------------------

   subroutine test_time_accurate_long(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 200
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), mass0, t_end
      integer :: ierr, k

      call alloc_adv_diff_problem(p, nz, ierr)
      p% dm = 1d27
      do k = 2, nz
         p% sig(k) = 1d24
         p% w(k) = 3d22
      end do
      p% x_cap = 1d0
      p% x_floor = 0.2d0
      x = 0.27d0
      mass0 = total_mass(p, x)
      t_end = 1d13   ! ~1d10 face timescales
      call adv_diff_integrate(p, x, t_end, opts, stats, ierr)
      call assert_true('time-accurate (long): no error with default options', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_true('time-accurate (long): steady before total_time, time = total_time', stats% time == t_end &
         .and. stats% stop_reason == adv_diff_stop_steady_state, nfail)
      call assert_true('time-accurate (long): well within max_steps', stats% steps < opts% max_steps/10, nfail)
      call assert_true('time-accurate (long): dt capped at max_dt_factor*tau_max', &
         stats% dt_last <= opts% max_dt_factor*stats% tau_max*(1d0 + 1d-12), nfail)
      call assert_close('time-accurate (long): mass conserved', total_mass(p, x), mass0, 1d-13, 0d0, nfail)
      call assert_true('time-accurate (long): in bounds', all(x >= 0d0) .and. all(x <= 1d0), nfail)
      write(*,'(a,i5,a,1pd12.4,a,1pd12.4)') '          steps = ', stats% steps, '  tau_max = ', stats% tau_max, &
         '  last dt = ', stats% dt_last
      call free_adv_diff_problem(p)
   end subroutine test_time_accurate_long


   ! ---------------------------------------------------------------------------
   ! 13. counter-flow with a capped receiver and no diffusion around it: the accepted
   !     Newton iterate must satisfy the cell equation to round-off, otherwise the passive
   !     species no longer sum to 1 - x
   ! ---------------------------------------------------------------------------

   subroutine test_capped_receiver_counterflow(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 200, npas = 3
      real(dp), parameter :: widths(3) = [1d-1, 1d-3, 1d-4]
      type(adv_diff_problem) :: p, p_pas
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), xp(nz, npas), xp_new(nz, npas), mass_act0, mass_pas0(npas), worst_sum, worst_mass
      integer :: ierr, k, j, iw

      do iw = 1, size(widths)
         call alloc_adv_diff_problem(p, nz, ierr)
         call alloc_adv_diff_problem(p_pas, nz, ierr)
         p% dm = 1d27
         p% dm(nz/2:nz) = 1d27/33d0
         do k = 2, nz
            if (k <= nz/2) p% sig(k) = 1d24   ! convective envelope; no diffusion below
            p% w(k) = 3d22
         end do
         p% cap_width = widths(iw)
         p% x_cap = 1d0
         p% x_floor(1:nz/2+4) = adv_diff_huge
         do k = nz/2+5, nz-20
            p% x_floor(k) = 0.2d0 + 0.05d0*dble(k - nz/2)/dble(nz/2)
            p% x_cap(k) = 0.9d0
         end do
         p% x_floor(nz-19:nz) = adv_diff_huge
         p_pas% dm = p% dm
         p_pas% sig = p% sig
         x = 0.27d0
         do k = 1, nz
            xp(k, 1) = 0.70d0 - 1d-3*dble(k)/dble(nz)
            xp(k, 2) = 0.02d0 + 1d-3*dble(k)/dble(nz)
            xp(k, 3) = 1d0 - x(k) - xp(k, 1) - xp(k, 2)
         end do
         mass_act0 = total_mass(p, x)
         do j = 1, npas
            mass_pas0(j) = total_mass(p, xp(:, j))
         end do
         worst_sum = 0d0
         worst_mass = 0d0
         opts% x_upper_bound = 1d0
         call adv_diff_integrate(p, x, 0d0, opts, stats, ierr, after_step_passive)
         call assert_true('capped receiver: no error', ierr == 0, nfail)
         if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
         call assert_true('capped receiver: some cells reach the cap', count(x >= 0.9d0 - widths(iw)) > 0, nfail)
         ! raw defect of the identity sum_j x_j = 1 - x (the driver rescales the passive species to
         ! 1 - x after every step and only fails above 1d-6); before the keep-best-iterate rule in
         ! the Newton loop this reached ~1d-8
         call assert_close('capped receiver: max |sum_j x_j - 1| (raw, no rescaling)', &
            worst_sum, 0d0, 0d0, 1d-9, nfail)
         call assert_close('capped receiver: species masses', worst_mass, 0d0, 0d0, 1d-10, nfail)
         call assert_close('capped receiver: active mass', total_mass(p, x), mass_act0, 1d-13, 0d0, nfail)
         write(*,'(a,1pd9.2,a,i4,a,i4,a,1pd10.3)') '          cap_width = ', widths(iw), '  steps = ', stats% steps, &
            '  rejected = ', stats% rejected_steps, '  worst |sum - 1| = ', worst_sum
         call free_adv_diff_problem(p)
         call free_adv_diff_problem(p_pas)
      end do

      contains

      subroutine after_step_passive(dt, x_new, flux_adv, flux_diff, change_other, ierr)
         real(dp), intent(in) :: dt
         real(dp), intent(in) :: x_new(:), flux_adv(:), flux_diff(:)
         real(dp), intent(out) :: change_other
         integer, intent(out) :: ierr
         integer :: kk, jj
         ierr = 0
         change_other = 0d0
         p_pas% w(1:nz+1) = -flux_adv(1:nz+1)
         do kk = 1, nz
            p_pas% scale(kk) = 1d0/max(1d0 - x_new(kk), 1d-30)
         end do
         call adv_diff_passive_step(p_pas, dt, npas, xp, xp_new, ierr)
         if (ierr /= 0) return
         change_other = maxval(abs(xp_new - xp))
         xp = xp_new
         do kk = 1, nz
            worst_sum = max(worst_sum, abs(x_new(kk) + sum(xp(kk, 1:npas)) - 1d0))
         end do
         do jj = 1, npas
            worst_mass = max(worst_mass, abs(total_mass(p_pas, xp(:, jj)) - mass_pas0(jj))/mass_pas0(jj))
         end do
         if (size(flux_diff) /= nz+1) ierr = -1
      end subroutine after_step_passive

   end subroutine test_capped_receiver_counterflow



   ! ---------------------------------------------------------------------------
   ! 14. no false steady state when a slow process feeds a donor sitting at its floor
   ! ---------------------------------------------------------------------------

   subroutine test_no_false_steady_state(nfail)
      ! cell 1 diffuses slowly into cell 2 (tau ~ 3d12); cell 2 sits exactly at its floor and
      ! is the donor of a fast advective face to cell 3 (tau = 0.1). The first steps change
      ! nothing at round-off level; the solver must not stop there.
      integer, intent(inout) :: nfail
      type(adv_diff_problem) :: p
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(3), x0(3)
      integer :: ierr

      call alloc_adv_diff_problem(p, 3, ierr)
      p% dm = 1d0
      p% sig(2) = 1d-13/0.3d0
      p% w(3) = 10d0
      p% x_floor = adv_diff_huge
      p% x_floor(2) = 0.2d0
      x0 = [0.5d0, 0.2d0, 0d0]

      x = x0
      call adv_diff_integrate(p, x, 0d0, opts, stats, ierr)
      call assert_true('slow feed (steady): no error', ierr == 0, nfail)
      if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
      call assert_close('slow feed (steady): x(1) drained to the floor', x(1), 0.2d0, 0d0, 1d-6, nfail)
      call assert_close('slow feed (steady): x(3) received the rain', x(3), 0.3d0, 0d0, 1d-6, nfail)
      call assert_close('slow feed (steady): mass conserved', total_mass(p, x), total_mass(p, x0), 1d-13, 0d0, nfail)

      x = x0
      call adv_diff_integrate(p, x, 1d14, opts, stats, ierr)
      call assert_true('slow feed (time 1d14 >> tau_max): no error', ierr == 0, nfail)
      call assert_close('slow feed (time 1d14): x(3) received the rain', x(3), 0.3d0, 0d0, 1d-4, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_no_false_steady_state


   ! ---------------------------------------------------------------------------
   ! 15. the change reported by the callback enters the steady-state test
   ! ---------------------------------------------------------------------------

   subroutine test_passive_species_in_steady_test(nfail)
      ! the active species is already steady (uniform, every donor at its floor), the passive
      ! species have a gradient: both modes must keep going until the passive species are mixed
      integer, parameter :: nz = 20, npas = 2
      integer, intent(inout) :: nfail
      type(adv_diff_problem) :: p, p_pas
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      real(dp) :: x(nz), xp(nz, npas), xp_new(nz, npas), xp0(nz, npas)
      real(dp) :: total_time
      integer :: ierr, k, imode

      call alloc_adv_diff_problem(p, nz, ierr)
      call alloc_adv_diff_problem(p_pas, nz, ierr)
      p% dm = 1d0
      p% sig(2:nz) = 1d0
      p% w(2:nz) = 0.5d0
      p% x_floor = 0.3d0
      p_pas% dm = p% dm
      p_pas% sig = p% sig
      do k = 1, nz
         xp0(k, 1) = 0.7d0*dble(k - 1)/dble(nz - 1)
         xp0(k, 2) = 0.7d0 - xp0(k, 1)
      end do
      do imode = 1, 2
         x = 0.3d0
         xp = xp0
         total_time = 0d0
         if (imode == 2) total_time = 1d4
         call adv_diff_integrate(p, x, total_time, opts, stats, ierr, after_step_passive)
         call assert_true('passive in steady test: no error', ierr == 0, nfail)
         if (ierr /= 0) write(*,'(a,a)') '      ', trim(stats% message)
         call assert_close('passive in steady test: passive species mixed', maxval(xp(:, 1)) - minval(xp(:, 1)), &
            0d0, 0d0, 1d-6, nfail)
         call assert_close('passive in steady test: mean kept', sum(xp(:, 1))/dble(nz), 0.35d0, 1d-12, 0d0, nfail)
      end do
      call free_adv_diff_problem(p)
      call free_adv_diff_problem(p_pas)

      contains

      subroutine after_step_passive(dt, x_new, flux_adv, flux_diff, change_other, ierr)
         real(dp), intent(in) :: dt
         real(dp), intent(in) :: x_new(:), flux_adv(:), flux_diff(:)
         real(dp), intent(out) :: change_other
         integer, intent(out) :: ierr
         integer :: kk
         ierr = 0
         p_pas% w(1:nz+1) = -flux_adv(1:nz+1)
         do kk = 1, nz
            p_pas% scale(kk) = 1d0/max(1d0 - x_new(kk), 1d-30)
         end do
         call adv_diff_passive_step(p_pas, dt, npas, xp, xp_new, ierr)
         change_other = maxval(abs(xp_new - xp))
         xp = xp_new
         if (size(flux_diff) /= nz+1) ierr = -1
      end subroutine after_step_passive

   end subroutine test_passive_species_in_steady_test

end module test_advection_diffusion
