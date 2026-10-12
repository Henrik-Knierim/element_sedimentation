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

! Generic conservative 1-D advection-diffusion solver on a Lagrangian (mass) grid.
!
! Cells k = 1..nz (k = 1 outermost), faces k = 1..nz+1 where face k is the outer
! face of cell k.  Both boundary faces carry zero flux.  For every cell
!
!    dm_k dx_k/dt = F_k - F_{k+1},   F_k = flux through face k, positive toward larger k
!
!    F_k = sig_k (x_{k-1} - x_k)                                       (diffusion)
!        + w_k scale_d max(x_d - x_floor_d, 0) throttle_r              (advection)
!
! with donor cell d = k-1 and receiver r = k for w_k > 0 (d = k, r = k-1 for
! w_k < 0) and throttle_r = min(1, max(0, (x_cap_r - x_r)/cap_width)).
!
! Time integration is backward Euler with first-order upwinding for the
! advective term.  The nonlinear (piecewise linear) step is solved by a
! semismooth Newton iteration with backtracking on a tridiagonal Jacobian.
! Steady states are found by pseudo-transient continuation.  Total mass of x is
! conserved to round-off by construction (the cell equations telescope).
!
! The module has no state of its own and is re-entrant.

module advection_diffusion

   use const_def, only: dp
   use utils_lib, only: mesa_error, is_bad

   implicit none

   private

   public :: adv_diff_problem, adv_diff_options, adv_diff_stats
   public :: alloc_adv_diff_problem, free_adv_diff_problem, check_adv_diff_problem
   public :: adv_diff_fluxes, adv_diff_timescales
   public :: adv_diff_implicit_step, adv_diff_passive_step, adv_diff_integrate
   public :: adv_diff_error_message, adv_diff_after_step_interface

   integer, parameter, public :: adv_diff_ok = 0
   integer, parameter, public :: adv_diff_err_bad_input = 1       ! sizes, NaN, negative dm/sig, nonzero boundary flux, cap_width <= 0, ...
   integer, parameter, public :: adv_diff_err_singular = 2        ! zero/negative pivot in tridiagonal solve
   integer, parameter, public :: adv_diff_err_newton = 3          ! Newton did not converge within max_newton_iters
   integer, parameter, public :: adv_diff_err_no_steady_state = 4 ! max_steps reached before tol_steady (or before total_time)
   integer, parameter, public :: adv_diff_err_conservation = 5    ! mass conservation violated beyond conservation_rtol
   integer, parameter, public :: adv_diff_err_bounds = 6          ! x < -1d-12 or (if enabled) x > x_upper_bound + 1d-12
   integer, parameter, public :: adv_diff_err_dt_too_small = 7    ! dt shrank below dt_min_factor * tau_min
   integer, parameter, public :: adv_diff_err_callback = 8        ! after_step callback returned ierr /= 0

   integer, parameter, public :: adv_diff_stop_steady_state = 1
   integer, parameter, public :: adv_diff_stop_total_time = 2
   integer, parameter, public :: adv_diff_stop_change_limit = 3
   integer, parameter, public :: adv_diff_stop_nothing_to_do = 4  ! no active face

   real(dp), parameter, public :: adv_diff_huge = 1d99             ! "no floor"/"no cap" sentinel

   ! private numerical constants
   real(dp), parameter :: bounds_tol = 1d-12                 ! tolerated excursion below 0 / above x_upper_bound (clipped)
   real(dp), parameter :: roundoff_factor = 1d2              ! residual <= roundoff_factor*epsilon*scale counts as zero
   real(dp), parameter :: change_limit_fraction = 0.9d0      ! change limit counts as reached above this fraction of it
   real(dp), parameter :: change_limit_overshoot = 1d-12     ! relative overshoot of the change limit that is tolerated
   integer, parameter :: max_backtracks = 30                 ! halvings of the Newton step in the line search

   type adv_diff_problem
      integer :: nz = 0
      real(dp), allocatable :: dm(:)        ! (nz)   cell masses [g], > 0
      real(dp), allocatable :: sig(:)       ! (nz+1) diffusive conductance across face k [g/s], >= 0; sig(1) = sig(nz+1) = 0
      real(dp), allocatable :: w(:)         ! (nz+1) advective carrier mass flow rate across face k [g/s], > 0 toward larger k; w(1) = w(nz+1) = 0
      real(dp), allocatable :: x_floor(:)   ! (nz)   only (x - x_floor)_+ is advected out of a cell; default 0; adv_diff_huge means "never advects"
      real(dp), allocatable :: x_cap(:)     ! (nz)   advective inflow throttled to zero as x -> x_cap; default adv_diff_huge
      real(dp), allocatable :: scale(:)     ! (nz)   multiplies the advected quantity of the donor cell; default 1
      real(dp) :: cap_width = 1d-3          ! absolute width of the throttle ramp below x_cap
   end type adv_diff_problem

   type adv_diff_options
      real(dp) :: tol_steady = 1d-8
      real(dp) :: tol_newton = 1d-12
      real(dp) :: tol_residual = 1d-8       ! Newton converged only if every cell-equation defect |r_k|/dm_k (in units of x)
                                            ! is <= this or at the round-off floor of its Jacobian row
      integer  :: max_newton_iters = 50
      integer  :: max_steps = 1000
      real(dp) :: dt_initial_factor = 0.5d0
      real(dp) :: dt_growth_max = 10d0
      real(dp) :: dt_growth_min = 2d0       ! dt grows by at least this factor per accepted step (> 1 prevents a stall when the residual plateaus)
      real(dp) :: dt_shrink = 0.25d0
      real(dp) :: dt_min_factor = 1d-12     ! give up when dt < dt_min_factor * tau_min
      real(dp) :: max_dt_factor = 1d4       ! dt <= max_dt_factor * tau_max
      real(dp) :: max_abs_change = -1d0     ! > 0: stop once max |x - x_start| reaches this
      real(dp) :: conservation_rtol = 1d-10
      real(dp) :: x_upper_bound = -1d0      ! > 0: error if x exceeds it by more than 1d-12 (use 1 for mass fractions)
      logical  :: report = .false.          ! print one line per accepted step (explicit formats)
   end type adv_diff_options

   type adv_diff_stats
      integer  :: steps = 0                 ! accepted steps
      integer  :: rejected_steps = 0
      integer  :: newton_iters = 0          ! total
      integer  :: stop_reason = 0
      integer  :: k_fail = 0
      integer  :: num_clipped = 0           ! tiny negatives set to zero
      real(dp) :: time = 0d0                ! pseudo- or physical time integrated
      real(dp) :: dt_last = 0d0
      real(dp) :: tau_min = 0d0, tau_max = 0d0
      real(dp) :: residual = 0d0            ! last dimensionless steady-state residual
      real(dp) :: max_change_last = 0d0     ! max |x_new - x_old| of last step
      real(dp) :: max_change_total = 0d0    ! max |x - x_start|
      real(dp) :: conservation_error = 0d0  ! max relative conservation error seen
      character(len=256) :: message = ''
   end type adv_diff_stats

   abstract interface
      subroutine adv_diff_after_step_interface(dt, x_new, flux_adv, flux_diff, change_other, ierr)
         import :: dp
         real(dp), intent(in) :: dt
         real(dp), intent(in) :: x_new(:)       ! (nz)
         real(dp), intent(in) :: flux_adv(:)    ! (nz+1) advective flux of the active species through each face, + toward larger k
         real(dp), intent(in) :: flux_diff(:)   ! (nz+1)
         real(dp), intent(out) :: change_other  ! largest change in this step of whatever the callback advances (e.g. the
                                                ! passive species); enters the steady-state test like the change of x
         integer, intent(out) :: ierr
      end subroutine adv_diff_after_step_interface
   end interface

contains

   ! ---------------------------------------------------------------------------
   ! problem setup
   ! ---------------------------------------------------------------------------

   subroutine alloc_adv_diff_problem(p, nz, ierr)
      ! allocate the arrays of p for nz cells and fill in the defaults
      ! (sig = 0, w = 0, x_floor = 0, x_cap = adv_diff_huge, scale = 1, dm = 0).
      ! dm must be set by the caller.
      type(adv_diff_problem), intent(inout) :: p
      integer, intent(in) :: nz
      integer, intent(out) :: ierr
      integer :: istat
      ierr = adv_diff_ok
      call free_adv_diff_problem(p)
      if (nz < 1) then
         ierr = adv_diff_err_bad_input
         return
      end if
      allocate(p% dm(nz), p% x_floor(nz), p% x_cap(nz), p% scale(nz), &
         p% sig(nz+1), p% w(nz+1), stat=istat)
      if (istat /= 0) call mesa_error(__FILE__, __LINE__, 'alloc_adv_diff_problem: allocate failed')
      p% nz = nz
      p% dm = 0d0
      p% sig = 0d0
      p% w = 0d0
      p% x_floor = 0d0
      p% x_cap = adv_diff_huge
      p% scale = 1d0
      p% cap_width = 1d-3
   end subroutine alloc_adv_diff_problem


   subroutine free_adv_diff_problem(p)
      type(adv_diff_problem), intent(inout) :: p
      if (allocated(p% dm)) deallocate(p% dm)
      if (allocated(p% sig)) deallocate(p% sig)
      if (allocated(p% w)) deallocate(p% w)
      if (allocated(p% x_floor)) deallocate(p% x_floor)
      if (allocated(p% x_cap)) deallocate(p% x_cap)
      if (allocated(p% scale)) deallocate(p% scale)
      p% nz = 0
   end subroutine free_adv_diff_problem


   logical function problem_is_allocated(p)
      type(adv_diff_problem), intent(in) :: p
      problem_is_allocated = p% nz >= 1 .and. allocated(p% dm) .and. allocated(p% sig) &
         .and. allocated(p% w) .and. allocated(p% x_floor) .and. allocated(p% x_cap) &
         .and. allocated(p% scale)
   end function problem_is_allocated


   subroutine check_adv_diff_problem(p, ierr, message)
      ! validate p; ierr = adv_diff_err_bad_input and a description in message on failure
      type(adv_diff_problem), intent(in) :: p
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer :: k, nz
      ierr = adv_diff_ok
      message = ''
      nz = p% nz
      if (.not. problem_is_allocated(p)) then
         call fail('adv_diff_problem is not allocated (call alloc_adv_diff_problem first)')
         return
      end if
      if (size(p% dm) /= nz .or. size(p% x_floor) /= nz .or. size(p% x_cap) /= nz &
            .or. size(p% scale) /= nz .or. size(p% sig) /= nz+1 .or. size(p% w) /= nz+1) then
         call fail('adv_diff_problem array sizes are inconsistent with nz')
         return
      end if
      if (is_bad(p% cap_width) .or. p% cap_width <= 0d0) then
         call fail('cap_width must be > 0')
         return
      end if
      do k = 1, nz
         if (is_bad(p% dm(k)) .or. p% dm(k) <= 0d0) then
            call fail_k('dm must be > 0 and finite at cell', k)
            return
         end if
         if (is_bad(p% x_floor(k))) then
            call fail_k('x_floor is NaN or Inf at cell', k)
            return
         end if
         if (is_bad(p% x_cap(k))) then
            call fail_k('x_cap is NaN or Inf at cell', k)
            return
         end if
         if (is_bad(p% scale(k)) .or. p% scale(k) < 0d0) then
            call fail_k('scale must be >= 0 and finite at cell', k)
            return
         end if
      end do
      do k = 1, nz+1
         if (is_bad(p% sig(k)) .or. p% sig(k) < 0d0) then
            call fail_k('sig must be >= 0 and finite at face', k)
            return
         end if
         if (is_bad(p% w(k))) then
            call fail_k('w is NaN or Inf at face', k)
            return
         end if
      end do
      if (p% sig(1) /= 0d0 .or. p% w(1) /= 0d0) then
         call fail('boundary face 1 must have sig = w = 0')
         return
      end if
      if (p% sig(nz+1) /= 0d0 .or. p% w(nz+1) /= 0d0) then
         call fail_k('boundary face nz+1 must have sig = w = 0, nz+1 =', nz+1)
         return
      end if

      contains

      subroutine fail(str)
         character(len=*), intent(in) :: str
         ierr = adv_diff_err_bad_input
         message = str
      end subroutine fail

      subroutine fail_k(str, kk)
         character(len=*), intent(in) :: str
         integer, intent(in) :: kk
         ierr = adv_diff_err_bad_input
         write(message, '(a,1x,i0)') str, kk
      end subroutine fail_k

   end subroutine check_adv_diff_problem


   subroutine check_adv_diff_options(opts, ierr, message)
      type(adv_diff_options), intent(in) :: opts
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      ierr = adv_diff_ok
      message = ''
      if (is_bad(opts% tol_steady) .or. opts% tol_steady <= 0d0) then
         message = 'options: tol_steady must be > 0'
      else if (is_bad(opts% tol_newton) .or. opts% tol_newton <= 0d0) then
         message = 'options: tol_newton must be > 0'
      else if (is_bad(opts% tol_residual) .or. opts% tol_residual <= 0d0) then
         message = 'options: tol_residual must be > 0'
      else if (opts% max_newton_iters < 1) then
         message = 'options: max_newton_iters must be >= 1'
      else if (opts% max_steps < 1) then
         message = 'options: max_steps must be >= 1'
      else if (is_bad(opts% dt_initial_factor) .or. opts% dt_initial_factor <= 0d0) then
         message = 'options: dt_initial_factor must be > 0'
      else if (is_bad(opts% dt_growth_max) .or. opts% dt_growth_max < 1d0) then
         message = 'options: dt_growth_max must be >= 1'
      else if (is_bad(opts% dt_growth_min) .or. opts% dt_growth_min < 1d0) then
         message = 'options: dt_growth_min must be >= 1'
      else if (is_bad(opts% dt_shrink) .or. opts% dt_shrink <= 0d0 .or. opts% dt_shrink >= 1d0) then
         message = 'options: dt_shrink must be in (0,1)'
      else if (is_bad(opts% dt_min_factor) .or. opts% dt_min_factor <= 0d0) then
         message = 'options: dt_min_factor must be > 0'
      else if (is_bad(opts% max_dt_factor) .or. opts% max_dt_factor <= 0d0) then
         message = 'options: max_dt_factor must be > 0'
      else if (is_bad(opts% max_abs_change)) then
         message = 'options: max_abs_change is NaN or Inf'
      else if (is_bad(opts% conservation_rtol) .or. opts% conservation_rtol < 0d0) then
         message = 'options: conservation_rtol must be >= 0'
      else if (is_bad(opts% x_upper_bound)) then
         message = 'options: x_upper_bound is NaN or Inf'
      end if
      if (len_trim(message) > 0) ierr = adv_diff_err_bad_input
   end subroutine check_adv_diff_options


   ! ---------------------------------------------------------------------------
   ! fluxes and timescales
   ! ---------------------------------------------------------------------------

   subroutine adv_diff_timescales(p, tau_min, tau_max, num_active_faces)
      ! tau_k = min(dm_{k-1}, dm_k) / (sig_k + |w_k| scale_d) for every interior face with
      ! a nonzero coefficient; tau_min/tau_max are the extrema over these faces (0 if none).
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(out) :: tau_min, tau_max
      integer, intent(out) :: num_active_faces
      integer :: k
      real(dp) :: rate, tau
      tau_min = 0d0
      tau_max = 0d0
      num_active_faces = 0
      do k = 2, p% nz
         rate = p% sig(k)
         if (p% w(k) > 0d0) then
            rate = rate + p% w(k)*p% scale(k-1)
         else if (p% w(k) < 0d0) then
            rate = rate - p% w(k)*p% scale(k)
         end if
         if (rate <= 0d0) cycle
         tau = min(p% dm(k-1), p% dm(k))/rate
         if (num_active_faces == 0) then
            tau_min = tau
            tau_max = tau
         else
            tau_min = min(tau_min, tau)
            tau_max = max(tau_max, tau)
         end if
         num_active_faces = num_active_faces + 1
      end do
   end subroutine adv_diff_timescales


   subroutine adv_diff_fluxes(p, x, flux_adv, flux_diff)
      ! fluxes through all faces for state x; boundary faces carry zero flux
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: x(:)             ! (nz)
      real(dp), intent(out) :: flux_adv(:)     ! (nz+1)
      real(dp), intent(out) :: flux_diff(:)    ! (nz+1)
      call face_fluxes(p, x, flux_adv, flux_diff)
   end subroutine adv_diff_fluxes


   subroutine face_fluxes(p, x, flux_adv, flux_diff, dadv_dxm, dadv_dxp, kink_slope)
      ! Fluxes through all faces and, optionally, the generalized derivatives of the
      ! advective flux with respect to x(k-1) (dadv_dxm) and x(k) (dadv_dxp).
      ! At an exact kink of the piecewise linear flux (x_d == x_floor_d, x_r == x_cap_r or
      ! x_r == x_cap_r - cap_width) the derivative of the inactive side (zero) is used unless
      ! kink_slope is present and true, in which case the sloped side is used.
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: flux_adv(:), flux_diff(:)
      real(dp), intent(out), optional :: dadv_dxm(:), dadv_dxp(:)
      logical, intent(in), optional :: kink_slope
      integer :: k, d, r, nz
      real(dp) :: ws, y, z, q, dq, thr, dthr, dfd, dfr
      logical :: want_derivs, use_slope

      nz = p% nz
      want_derivs = present(dadv_dxm) .and. present(dadv_dxp)
      use_slope = .false.
      if (present(kink_slope)) use_slope = kink_slope

      flux_adv(1) = 0d0
      flux_diff(1) = 0d0
      flux_adv(nz+1) = 0d0
      flux_diff(nz+1) = 0d0
      if (want_derivs) then
         dadv_dxm(1) = 0d0
         dadv_dxp(1) = 0d0
         dadv_dxm(nz+1) = 0d0
         dadv_dxp(nz+1) = 0d0
      end if

      do k = 2, nz
         flux_diff(k) = p% sig(k)*(x(k-1) - x(k))
         flux_adv(k) = 0d0
         if (want_derivs) then
            dadv_dxm(k) = 0d0
            dadv_dxp(k) = 0d0
         end if
         if (p% w(k) == 0d0) cycle
         if (p% w(k) > 0d0) then
            d = k-1
            r = k
         else
            d = k
            r = k-1
         end if
         ws = p% w(k)*p% scale(d)
         if (ws == 0d0) cycle
         ! donor: only the part above the floor is advected
         q = 0d0
         dq = 0d0
         if (p% x_floor(d) < adv_diff_huge) then
            y = x(d) - p% x_floor(d)
            if (y > 0d0) then
               q = y
               dq = 1d0
            else if (y == 0d0 .and. use_slope) then
               dq = 1d0
            end if
         end if
         ! receiver: inflow throttled to zero as x_r approaches the cap
         thr = 1d0
         dthr = 0d0
         if (p% x_cap(r) < adv_diff_huge) then
            z = (p% x_cap(r) - x(r))/p% cap_width
            if (z <= 0d0) then
               thr = 0d0
               if (z == 0d0 .and. use_slope) dthr = -1d0/p% cap_width
            else if (z < 1d0) then
               thr = z
               dthr = -1d0/p% cap_width
            else if (z == 1d0 .and. use_slope) then
               dthr = -1d0/p% cap_width
            end if
         end if
         flux_adv(k) = ws*q*thr
         if (want_derivs) then
            dfd = ws*dq*thr    ! d flux / d x_donor
            dfr = ws*q*dthr    ! d flux / d x_receiver
            if (d == k-1) then
               dadv_dxm(k) = dfd
               dadv_dxp(k) = dfr
            else
               dadv_dxm(k) = dfr
               dadv_dxp(k) = dfd
            end if
         end if
      end do
   end subroutine face_fluxes


   real(dp) function steady_residual(p, flux_adv, flux_diff, tau_max) result(res)
      ! dimensionless steady-state residual max_k |F_k - F_{k+1}|/dm_k * tau_max
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: flux_adv(:), flux_diff(:), tau_max
      integer :: k
      real(dp) :: net
      res = 0d0
      do k = 1, p% nz
         net = (flux_adv(k) + flux_diff(k)) - (flux_adv(k+1) + flux_diff(k+1))
         res = max(res, abs(net)/p% dm(k))
      end do
      res = res*tau_max
   end function steady_residual


   ! ---------------------------------------------------------------------------
   ! tridiagonal solver (Thomas algorithm, no pivoting)
   ! ---------------------------------------------------------------------------

   subroutine factor_tridiagonal(n, dm, sub, sup, piv, cp, k_fail, ierr)
      ! LU factorization of the tridiagonal matrix A with A(k,k-1) = sub(k) <= 0,
      ! A(k,k+1) = sup(k) <= 0 and column sums equal to dm(k) > 0, i.e. with the diagonal
      ! diag(k) = dm(k) + |sup(k-1)| + |sub(k+1)| (assemble_jacobian builds exactly such matrices).
      ! piv(k) are the pivots, cp(k) = sup(k)/piv(k).
      !
      ! The usual recurrence piv(k) = diag(k) - sub(k)*cp(k-1) subtracts two numbers of size
      ! dt*sig and loses all digits once dt*sig/dm exceeds ~1/epsilon. With the column-sum
      ! identity the same pivots follow without any subtraction:
      !    u(1) = dm(1),  u(k) = dm(k) + |sup(k-1)| u(k-1)/piv(k-1),  piv(k) = u(k) + |sub(k+1)|
      ! (u(k) is the column sum of the eliminated matrix), so piv(k) >= dm(k) > 0 for any
      ! stiffness. A non-positive or NaN pivot can therefore only come from bad input.
      integer, intent(in) :: n
      real(dp), intent(in) :: dm(:), sub(:), sup(:)
      real(dp), intent(out) :: piv(:), cp(:)
      integer, intent(out) :: k_fail, ierr
      integer :: k
      real(dp) :: u
      ierr = adv_diff_ok
      k_fail = 0
      u = dm(1)
      do k = 1, n
         if (k > 1) u = dm(k) + abs(sup(k-1))*u/piv(k-1)
         if (k < n) then
            piv(k) = u + abs(sub(k+1))
         else
            piv(k) = u
         end if
         if (.not. (piv(k) > 0d0) .or. is_bad(piv(k))) then
            ierr = adv_diff_err_singular
            k_fail = k
            return
         end if
         cp(k) = sup(k)/piv(k)
      end do
   end subroutine factor_tridiagonal


   subroutine solve_factored_tridiagonal(n, sub, piv, cp, rhs, x)
      integer, intent(in) :: n
      real(dp), intent(in) :: sub(:), piv(:), cp(:), rhs(:)
      real(dp), intent(out) :: x(:)
      integer :: k
      x(1) = rhs(1)/piv(1)
      do k = 2, n
         x(k) = (rhs(k) - sub(k)*x(k-1))/piv(k)
      end do
      do k = n-1, 1, -1
         x(k) = x(k) - cp(k)*x(k+1)
      end do
   end subroutine solve_factored_tridiagonal


   subroutine solve_tridiagonal(n, dm, sub, sup, rhs, x, k_fail, ierr)
      ! solve the tridiagonal system (see factor_tridiagonal) for one right-hand side
      integer, intent(in) :: n
      real(dp), intent(in) :: dm(:), sub(:), sup(:), rhs(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: k_fail, ierr
      real(dp) :: piv(n), cp(n)
      call factor_tridiagonal(n, dm, sub, sup, piv, cp, k_fail, ierr)
      if (ierr /= 0) return
      call solve_factored_tridiagonal(n, sub, piv, cp, rhs, x)
   end subroutine solve_tridiagonal


   subroutine assemble_jacobian(p, dt, dF_dxm, dF_dxp, sub, diag, sup)
      ! Jacobian of r_k = dm_k (x_k - x_old_k) - dt (F_k - F_{k+1}) given the face derivatives
      ! dF_dxm(k) = dF_k/dx_{k-1} and dF_dxp(k) = dF_k/dx_k (both including diffusion).
      ! Column sums equal dm_k, off-diagonals are <= 0: a column diagonally dominant M-matrix.
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: dF_dxm(:), dF_dxp(:)
      real(dp), intent(out) :: sub(:), diag(:), sup(:)
      integer :: k
      do k = 1, p% nz
         diag(k) = p% dm(k) - dt*(dF_dxp(k) - dF_dxm(k+1))
         sub(k) = -dt*dF_dxm(k)
         sup(k) = dt*dF_dxp(k+1)
      end do
   end subroutine assemble_jacobian


   ! ---------------------------------------------------------------------------
   ! one backward Euler step of the nonlinear problem
   ! ---------------------------------------------------------------------------

   subroutine adv_diff_implicit_step(p, dt, x_old, x_new, opts, stats, ierr)
      ! Solve dm_k (x_new_k - x_old_k) = dt (F_k(x_new) - F_{k+1}(x_new)) by semismooth Newton
      ! with backtracking. stats% newton_iters counts the linear solves; stats% k_fail and
      ! stats% message are set on failure. After the solve the step is checked for mass
      ! conservation and bounds (tiny negatives are clipped and counted in stats% num_clipped).
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: x_old(:)
      real(dp), intent(out) :: x_new(:)
      type(adv_diff_options), intent(in) :: opts
      type(adv_diff_stats), intent(inout) :: stats
      integer, intent(out) :: ierr

      integer :: nz, k, iter, ib, k_fail, num_small
      real(dp) :: merit, merit_trial, rel_res, rel_res_trial, dmax, lambda, eps
      logical :: converged, accepted, kink_slope, small_improves
      real(dp), allocatable :: flux_adv(:), flux_diff(:), dadv_dxm(:), dadv_dxp(:), &
         dF_dxm(:), dF_dxp(:), r(:), r_trial(:), sub(:), diag(:), sup(:), delta(:), &
         x_trial(:), fa_trial(:), fd_trial(:)

      ierr = adv_diff_ok
      if (.not. problem_is_allocated(p)) &
         call mesa_error(__FILE__, __LINE__, 'adv_diff_implicit_step: problem not allocated')
      nz = p% nz
      if (size(x_old) /= nz .or. size(x_new) /= nz) then
         ierr = adv_diff_err_bad_input
         stats% message = 'adv_diff_implicit_step: size of x_old or x_new differs from nz'
         return
      end if
      if (is_bad(dt) .or. dt <= 0d0) then
         ierr = adv_diff_err_bad_input
         stats% message = 'adv_diff_implicit_step: dt must be > 0 and finite'
         return
      end if
      do k = 1, nz
         if (is_bad(x_old(k))) then
            ierr = adv_diff_err_bad_input
            stats% k_fail = k
            write(stats% message, '(a,i0)') 'adv_diff_implicit_step: x_old is NaN or Inf at cell ', k
            return
         end if
      end do

      allocate(flux_adv(nz+1), flux_diff(nz+1), dadv_dxm(nz+1), dadv_dxp(nz+1), &
         dF_dxm(nz+1), dF_dxp(nz+1), fa_trial(nz+1), fd_trial(nz+1), &
         r(nz), r_trial(nz), sub(nz), diag(nz), sup(nz), delta(nz), x_trial(nz))

      eps = epsilon(1d0)
      kink_slope = .false.
      converged = .false.
      num_small = 0
      x_new(1:nz) = x_old(1:nz)
      call face_fluxes(p, x_new, flux_adv, flux_diff, dadv_dxm, dadv_dxp, kink_slope)
      call step_residual(p, dt, x_old, x_new, flux_adv, flux_diff, r, merit, rel_res)

      newton_loop: do iter = 1, opts% max_newton_iters

         ! residual at round-off: nothing left to do
         if (rel_res <= roundoff_factor*eps) then
            converged = .true.
            exit newton_loop
         end if

         ! generalized Jacobian and Newton correction
         stats% newton_iters = stats% newton_iters + 1
         do k = 1, nz+1
            dF_dxm(k) = p% sig(k) + dadv_dxm(k)
            dF_dxp(k) = -p% sig(k) + dadv_dxp(k)
         end do
         call assemble_jacobian(p, dt, dF_dxm, dF_dxp, sub, diag, sup)
         r_trial(1:nz) = -r(1:nz)
         call solve_tridiagonal(nz, p% dm, sub, sup, r_trial, delta, k_fail, ierr)
         if (ierr /= 0) then
            stats% k_fail = k_fail
            write(stats% message, '(a,i0,a,i0,a,1pd16.8,a)') &
               'adv_diff_implicit_step: singular pivot in tridiagonal solve at cell ', k_fail, &
               ' (Newton iteration ', iter, ', dt = ', dt, ')'
            return
         end if
         dmax = 0d0
         do k = 1, nz
            dmax = max(dmax, abs(delta(k)))
         end do
         if (is_bad(dmax)) then
            ierr = adv_diff_err_newton
            write(stats% message, '(a,i0,a,1pd16.8)') &
               'adv_diff_implicit_step: Newton correction is NaN or Inf at iteration ', iter, ', dt = ', dt
            return
         end if

         ! converged when the correction is negligible (twice in a row, or once the
         ! residual of the updated state is at round-off, which is checked at the top).
         ! A negligible correction is only applied directly if it does not increase the
         ! residual: next to a kink of the piecewise linear flux (e.g. a receiver cell at its
         ! cap, where the Jacobian carries dt*w/cap_width) a correction of 1d-12 in x can mean
         ! a large residual, and keeping the better iterate is what keeps sum(xa) = 1 to
         ! round-off in the passive counter-flow. If it does increase the residual, the line
         ! search below decides; only if that fails too is the iterate at its round-off floor.
         small_improves = .false.
         if (dmax <= opts% tol_newton) then
            x_trial(1:nz) = x_new(1:nz) + delta(1:nz)
            call face_fluxes(p, x_trial, fa_trial, fd_trial, dadv_dxm, dadv_dxp, kink_slope)
            call step_residual(p, dt, x_old, x_trial, fa_trial, fd_trial, r_trial, merit_trial, rel_res_trial)
            small_improves = merit_trial <= merit
         end if
         if (small_improves) then
            x_new(1:nz) = x_trial(1:nz)
            flux_adv(1:nz+1) = fa_trial(1:nz+1)
            flux_diff(1:nz+1) = fd_trial(1:nz+1)
            r(1:nz) = r_trial(1:nz)
            merit = merit_trial
            rel_res = rel_res_trial
            num_small = num_small + 1
            if (num_small >= 2 .and. residual_resolved(p, x_new, r, sub, diag, sup, opts% tol_residual)) then
               converged = .true.
               exit newton_loop
            end if
            cycle newton_loop
         end if
         num_small = 0

         ! backtracking line search on the max-norm of the residual (in units of x)
         accepted = .false.
         lambda = 1d0
         do ib = 1, max_backtracks
            x_trial(1:nz) = x_new(1:nz) + lambda*delta(1:nz)
            call face_fluxes(p, x_trial, fa_trial, fd_trial, dadv_dxm, dadv_dxp, kink_slope)
            call step_residual(p, dt, x_old, x_trial, fa_trial, fd_trial, r_trial, merit_trial, rel_res_trial)
            if (merit_trial < merit .or. rel_res_trial <= roundoff_factor*eps) then
               accepted = .true.
               exit
            end if
            lambda = 0.5d0*lambda
         end do
         if (accepted) then
            x_new(1:nz) = x_trial(1:nz)
            flux_adv(1:nz+1) = fa_trial(1:nz+1)
            flux_diff(1:nz+1) = fd_trial(1:nz+1)
            r(1:nz) = r_trial(1:nz)
            merit = merit_trial
            rel_res = rel_res_trial
         else if (dmax <= opts% tol_newton .and. &
               residual_resolved(p, x_new, r, sub, diag, sup, opts% tol_residual)) then
            ! the correction is negligible, no fraction of it reduces the residual, and the
            ! residual is small: the iterate is at its round-off floor. Keep it if it conserves
            ! mass; otherwise take the full correction, which conserves mass by construction
            ! (the Jacobian's column sums are dm), so the refinement that restores conservation
            ! is not skipped.
            if (relative_conservation_error(p, x_old, x_new) > 0.1d0*opts% conservation_rtol) then
               x_new(1:nz) = x_new(1:nz) + delta(1:nz)
               call face_fluxes(p, x_new, flux_adv, flux_diff, dadv_dxm, dadv_dxp, kink_slope)
               call step_residual(p, dt, x_old, x_new, flux_adv, flux_diff, r, merit, rel_res)
            end if
            converged = .true.
            exit newton_loop
         else if (.not. kink_slope) then
            ! the state may sit exactly on a kink of the piecewise linear flux with the
            ! wrong one-sided derivative chosen; retry once with the other choice
            kink_slope = .true.
            call face_fluxes(p, x_new, flux_adv, flux_diff, dadv_dxm, dadv_dxp, kink_slope)
         else
            ierr = adv_diff_err_newton
            stats% k_fail = maxloc(abs(r(1:nz))/p% dm(1:nz), dim=1)
            write(stats% message, '(a,i0,a,1pd16.8,a,1pd16.8,a,1pd16.8)') &
               'adv_diff_implicit_step: line search failed at Newton iteration ', iter, &
               ', dt = ', dt, ', residual = ', merit, ', max|delta| = ', dmax
            return
         end if

      end do newton_loop

      ! the last accepted correction may already have brought the residual to round-off
      if (.not. converged) converged = rel_res <= roundoff_factor*eps

      if (.not. converged) then
         ierr = adv_diff_err_newton
         stats% k_fail = maxloc(abs(r(1:nz))/p% dm(1:nz), dim=1)
         write(stats% message, '(a,i0,a,1pd16.8,a,1pd16.8)') &
            'adv_diff_implicit_step: Newton did not converge within max_newton_iters = ', &
            opts% max_newton_iters, ', dt = ', dt, ', residual = ', merit
         return
      end if

      ! conservation and bounds of the accepted solution
      call check_conservation(p, x_old, x_new, opts, stats, ierr)
      if (ierr /= 0) return
      call check_bounds(p, x_new, opts, stats, ierr)
      if (ierr /= 0) return

   end subroutine adv_diff_implicit_step


   subroutine step_residual(p, dt, x_old, x, flux_adv, flux_diff, r, merit, rel_res)
      ! r_k = dm_k (x_k - x_old_k) - dt (F_k - F_{k+1}); merit = max_k |r_k|/dm_k;
      ! rel_res = max_k |r_k| / (round-off scale of row k).
      ! The residual is evaluated in flux form (face differences of x first, then times sig),
      ! so its round-off scale is set by the fluxes themselves and not by |A||x|. A solution
      ! from the tridiagonal solve has a residual of order epsilon*|A||x| instead, so for
      ! stiff problems (dt sig >> dm) the Newton loop performs one or two extra solves that
      ! act as iterative refinement; this is what keeps mass conservation at round-off.
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: x_old(:), x(:), flux_adv(:), flux_diff(:)
      real(dp), intent(out) :: r(:)
      real(dp), intent(out) :: merit, rel_res
      integer :: k
      real(dp) :: net, rscale
      merit = 0d0
      rel_res = 0d0
      do k = 1, p% nz
         net = (flux_adv(k) + flux_diff(k)) - (flux_adv(k+1) + flux_diff(k+1))
         r(k) = p% dm(k)*(x(k) - x_old(k)) - dt*net
         merit = max(merit, abs(r(k))/p% dm(k))
         rscale = p% dm(k)*(abs(x(k)) + abs(x_old(k))) &
            + dt*(abs(flux_adv(k)) + abs(flux_diff(k)) + abs(flux_adv(k+1)) + abs(flux_diff(k+1)))
         rel_res = max(rel_res, abs(r(k))/max(rscale, tiny(1d0)))
      end do
   end subroutine step_residual


   logical function residual_resolved(p, x, r, sub, diag, sup, tol) result(ok)
      ! True if every cell-equation defect |r_k|/dm_k is <= tol or within roundoff_factor*epsilon
      ! of what the Jacobian row can resolve for this state (|J_k,k-1 x_k-1| + |J_kk x_k| + |J_k,k+1 x_k+1|).
      ! For very stiff rows (dt*sig/dm >> 1) the defect of an exact solution is dominated by
      ! round-off in the fluxes, so an absolute tolerance alone cannot be met.
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: x(:), r(:), sub(:), diag(:), sup(:), tol
      integer :: k, nz
      real(dp) :: jscale, allow
      nz = p% nz
      ok = .true.
      do k = 1, nz
         jscale = abs(diag(k)*x(k))
         if (k > 1) jscale = jscale + abs(sub(k)*x(k-1))
         if (k < nz) jscale = jscale + abs(sup(k)*x(k+1))
         allow = max(tol, roundoff_factor*epsilon(1d0)*jscale/p% dm(k))
         if (.not. (abs(r(k))/p% dm(k) <= allow)) then
            ok = .false.
            return
         end if
      end do
   end function residual_resolved


   real(dp) function relative_conservation_error(p, x_old, x_new) result(rel)
      ! |sum dm (x_new - x_old)| / max(sum dm |x_old|, 1d-12 sum dm)
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: x_old(:), x_new(:)
      integer :: k
      real(dp) :: mass_old, mass_tot, err
      mass_old = 0d0
      mass_tot = 0d0
      err = 0d0
      do k = 1, p% nz
         mass_old = mass_old + p% dm(k)*abs(x_old(k))
         mass_tot = mass_tot + p% dm(k)
         err = err + p% dm(k)*(x_new(k) - x_old(k))
      end do
      rel = abs(err)/max(mass_old, 1d-12*mass_tot)
      if (is_bad(rel)) rel = huge(1d0)
   end function relative_conservation_error


   subroutine check_conservation(p, x_old, x_new, opts, stats, ierr)
      ! |sum dm (x_new - x_old)| <= conservation_rtol * max(sum dm |x_old|, 1d-12 sum dm)
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: x_old(:), x_new(:)
      type(adv_diff_options), intent(in) :: opts
      type(adv_diff_stats), intent(inout) :: stats
      integer, intent(out) :: ierr
      integer :: k
      real(dp) :: mass_old, mass_tot, err, rel
      ierr = adv_diff_ok
      mass_old = 0d0
      mass_tot = 0d0
      err = 0d0
      do k = 1, p% nz
         mass_old = mass_old + p% dm(k)*abs(x_old(k))
         mass_tot = mass_tot + p% dm(k)
         err = err + p% dm(k)*(x_new(k) - x_old(k))
      end do
      rel = abs(err)/max(mass_old, 1d-12*mass_tot)
      if (is_bad(rel)) rel = huge(1d0)
      stats% conservation_error = max(stats% conservation_error, rel)
      if (rel > opts% conservation_rtol) then
         ierr = adv_diff_err_conservation
         write(stats% message, '(a,1pd16.8,a,1pd16.8,a,1pd16.8)') &
            'adv_diff_implicit_step: mass conservation violated: relative error ', rel, &
            ' > conservation_rtol = ', opts% conservation_rtol, ', absolute error = ', err
      end if
   end subroutine check_conservation


   subroutine check_bounds(p, x, opts, stats, ierr)
      ! x >= 0 (and x <= x_upper_bound if enabled); excursions up to bounds_tol are clipped and
      ! counted in stats% num_clipped, larger ones give adv_diff_err_bounds with k_fail
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(inout) :: x(:)
      type(adv_diff_options), intent(in) :: opts
      type(adv_diff_stats), intent(inout) :: stats
      integer, intent(out) :: ierr
      integer :: k
      ierr = adv_diff_ok
      do k = 1, p% nz
         if (is_bad(x(k))) then
            ierr = adv_diff_err_bounds
            stats% k_fail = k
            write(stats% message, '(a,i0)') 'adv_diff: x is NaN or Inf at cell ', k
            return
         end if
         if (x(k) < 0d0) then
            if (x(k) >= -bounds_tol) then
               x(k) = 0d0
               stats% num_clipped = stats% num_clipped + 1
            else
               ierr = adv_diff_err_bounds
               stats% k_fail = k
               write(stats% message, '(a,i0,a,1pd16.8)') 'adv_diff: x < 0 beyond tolerance at cell ', k, ': x = ', x(k)
               return
            end if
         end if
         if (opts% x_upper_bound > 0d0) then
            if (x(k) > opts% x_upper_bound) then
               if (x(k) <= opts% x_upper_bound + bounds_tol) then
                  x(k) = opts% x_upper_bound
                  stats% num_clipped = stats% num_clipped + 1
               else
                  ierr = adv_diff_err_bounds
                  stats% k_fail = k
                  write(stats% message, '(a,i0,a,1pd16.8,a,1pd16.8)') &
                     'adv_diff: x > x_upper_bound beyond tolerance at cell ', k, ': x = ', x(k), &
                     ', x_upper_bound = ', opts% x_upper_bound
                  return
               end if
            end if
         end if
      end do
   end subroutine check_bounds


   ! ---------------------------------------------------------------------------
   ! linear step for passive species (several right-hand sides, one factorization)
   ! ---------------------------------------------------------------------------

   subroutine adv_diff_passive_step(p, dt, nvar, x_old, x_new, ierr)
      ! One backward Euler step of the linear problem (x_floor and x_cap of p are ignored,
      ! i.e. x_floor = 0 and x_cap = huge) for nvar variables sharing the same coefficients.
      ! x_old and x_new are (nz, nvar). The matrix is an M-matrix, so x_new >= 0 if x_old >= 0.
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: dt
      integer, intent(in) :: nvar
      real(dp), intent(in) :: x_old(:,:)
      real(dp), intent(out) :: x_new(:,:)
      integer, intent(out) :: ierr

      integer :: nz, k, j, k_fail
      character(len=256) :: message
      real(dp), allocatable :: dF_dxm(:), dF_dxp(:), sub(:), diag(:), sup(:), piv(:), cp(:), rhs(:)

      ierr = adv_diff_ok
      if (.not. problem_is_allocated(p)) &
         call mesa_error(__FILE__, __LINE__, 'adv_diff_passive_step: problem not allocated')
      call check_adv_diff_problem(p, ierr, message)
      if (ierr /= 0) return
      nz = p% nz
      if (nvar < 1 .or. size(x_old, dim=1) /= nz .or. size(x_new, dim=1) /= nz &
            .or. size(x_old, dim=2) < nvar .or. size(x_new, dim=2) < nvar) then
         ierr = adv_diff_err_bad_input
         return
      end if
      if (is_bad(dt) .or. dt < 0d0) then
         ierr = adv_diff_err_bad_input
         return
      end if
      do j = 1, nvar
         do k = 1, nz
            if (is_bad(x_old(k,j))) then
               ierr = adv_diff_err_bad_input
               return
            end if
         end do
      end do

      allocate(dF_dxm(nz+1), dF_dxp(nz+1), sub(nz), diag(nz), sup(nz), piv(nz), cp(nz), rhs(nz))

      dF_dxm(1) = 0d0
      dF_dxp(1) = 0d0
      dF_dxm(nz+1) = 0d0
      dF_dxp(nz+1) = 0d0
      do k = 2, nz
         ! F_k = sig_k (x_{k-1} - x_k) + w_k scale_d x_d with upwind donor d
         dF_dxm(k) = p% sig(k) + max(p% w(k), 0d0)*p% scale(k-1)
         dF_dxp(k) = -p% sig(k) + min(p% w(k), 0d0)*p% scale(k)
      end do
      call assemble_jacobian(p, dt, dF_dxm, dF_dxp, sub, diag, sup)
      call factor_tridiagonal(nz, p% dm, sub, sup, piv, cp, k_fail, ierr)
      if (ierr /= 0) return

      do j = 1, nvar
         do k = 1, nz
            rhs(k) = p% dm(k)*x_old(k,j)
         end do
         call solve_factored_tridiagonal(nz, sub, piv, cp, rhs, x_new(1:nz,j))
         do k = 1, nz
            if (is_bad(x_new(k,j))) then
               ierr = adv_diff_err_singular
               return
            end if
         end do
      end do

   end subroutine adv_diff_passive_step


   ! ---------------------------------------------------------------------------
   ! driver: pseudo-transient continuation / time-accurate integration
   ! ---------------------------------------------------------------------------

   subroutine adv_diff_integrate(p, x, total_time, opts, stats, ierr, after_step)
      ! total_time <= 0: pseudo-transient continuation to the steady state (stop when the
      !    dimensionless residual max_k |dx_k|/dt * tau_max <= tol_steady).
      ! total_time > 0: integrate to exactly total_time; if the state becomes steady (same
      !    criterion) before that, stop with stop_reason = steady_state and time = total_time.
      ! opts% max_abs_change > 0 stops the integration once max |x - x_start| reaches it.
      ! after_step (optional) is called after every accepted step with the fluxes evaluated
      ! at the new state; a nonzero ierr from it aborts with adv_diff_err_callback.
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: total_time
      type(adv_diff_options), intent(in) :: opts
      type(adv_diff_stats), intent(out) :: stats
      integer, intent(out) :: ierr
      procedure(adv_diff_after_step_interface), optional :: after_step

      integer :: nz, k, num_active, ierr_cb
      real(dp) :: tau_min, tau_max, dt, dt_min, dt_max, time, remaining, res, res_prev, &
         growth, change_last, change_total, frac, change_other
      logical :: steady, last_step
      real(dp), allocatable :: x_start(:), x_old(:), x_new(:), flux_adv(:), flux_diff(:)

      stats = adv_diff_stats()
      ierr = adv_diff_ok

      if (.not. problem_is_allocated(p)) &
         call mesa_error(__FILE__, __LINE__, 'adv_diff_integrate: problem not allocated')
      call check_adv_diff_problem(p, ierr, stats% message)
      if (ierr /= 0) return
      call check_adv_diff_options(opts, ierr, stats% message)
      if (ierr /= 0) return
      nz = p% nz
      if (size(x) /= nz) then
         ierr = adv_diff_err_bad_input
         write(stats% message, '(a,i0,a,i0)') 'adv_diff_integrate: size(x) = ', size(x), ' differs from nz = ', nz
         return
      end if
      do k = 1, nz
         if (is_bad(x(k))) then
            ierr = adv_diff_err_bad_input
            stats% k_fail = k
            write(stats% message, '(a,i0)') 'adv_diff_integrate: x is NaN or Inf at cell ', k
            return
         end if
      end do
      if (is_bad(total_time)) then
         ierr = adv_diff_err_bad_input
         stats% message = 'adv_diff_integrate: total_time is NaN or Inf'
         return
      end if

      ! the initial state must satisfy the bounds (tiny excursions are clipped)
      call check_bounds(p, x, opts, stats, ierr)
      if (ierr /= 0) return

      call adv_diff_timescales(p, tau_min, tau_max, num_active)
      stats% tau_min = tau_min
      stats% tau_max = tau_max
      steady = total_time <= 0d0

      if (num_active == 0) then
         stats% stop_reason = adv_diff_stop_nothing_to_do
         if (.not. steady) stats% time = total_time
         return
      end if

      allocate(x_start(nz), x_old(nz), x_new(nz), flux_adv(nz+1), flux_diff(nz+1))
      x_start(1:nz) = x(1:nz)

      dt_min = opts% dt_min_factor*tau_min
      ! cap of the time step: beyond max_dt_factor*tau_max every mode is relaxed, and larger
      ! steps only raise the round-off floor of the cell equations (the Jacobian grows like
      ! dt/tau). In time-accurate mode the integration ends early once the state is steady
      ! (see below), so the cap does not limit how long a time can be covered.
      dt_max = opts% max_dt_factor*tau_max
      if (.not. steady) dt_max = min(dt_max, total_time)
      dt = min(opts% dt_initial_factor*tau_min, dt_max)
      time = 0d0

      call face_fluxes(p, x, flux_adv, flux_diff)
      res_prev = steady_residual(p, flux_adv, flux_diff, tau_max)

      step_loop: do

         if (stats% steps >= opts% max_steps) then
            ierr = adv_diff_err_no_steady_state
            if (steady) then
               write(stats% message, '(a,i0,a,1pd16.8,a,1pd16.8)') &
                  'adv_diff_integrate: no steady state after max_steps = ', opts% max_steps, &
                  ' steps; residual = ', stats% residual, ' > tol_steady = ', opts% tol_steady
            else
               write(stats% message, '(a,i0,a,1pd16.8,a,1pd16.8)') &
                  'adv_diff_integrate: max_steps = ', opts% max_steps, &
                  ' steps reached at time = ', time, ' < total_time = ', total_time
            end if
            return
         end if

         last_step = .false.
         if (.not. steady) then
            remaining = total_time - time
            if (dt >= remaining*(1d0 - 1d-12)) then
               dt = remaining
               last_step = .true.
            end if
         end if

         x_old(1:nz) = x(1:nz)
         call adv_diff_implicit_step(p, dt, x_old, x_new, opts, stats, ierr)
         if (ierr == adv_diff_err_newton .or. ierr == adv_diff_err_singular) then
            ! retry with a smaller step
            stats% rejected_steps = stats% rejected_steps + 1
            if (opts% report) write(*,'(a,i0,a,1pd16.8,a,a)') &
               ' adv_diff: step rejected (ierr = ', ierr, ', dt = ', dt, ') ', trim(stats% message)
            dt = dt*opts% dt_shrink
            if (dt < dt_min) then
               ierr = adv_diff_err_dt_too_small
               stats% dt_last = dt
               stats% message = 'adv_diff_integrate: dt below dt_min_factor*tau_min after repeated Newton failures; last: ' &
                  // trim(stats% message)
               return
            end if
            ierr = adv_diff_ok
            cycle step_loop
         else if (ierr /= 0) then
            stats% dt_last = dt
            return
         end if

         change_last = 0d0
         change_total = 0d0
         do k = 1, nz
            change_last = max(change_last, abs(x_new(k) - x_old(k)))
            change_total = max(change_total, abs(x_new(k) - x_start(k)))
         end do

         if (opts% max_abs_change > 0d0) then
            if (change_total > opts% max_abs_change*(1d0 + change_limit_overshoot)) then
               ! overshoot of the change limit: shrink dt toward the step that just reaches it
               stats% rejected_steps = stats% rejected_steps + 1
               frac = (opts% max_abs_change - stats% max_change_total)/(change_total - stats% max_change_total)
               dt = dt*max(0.1d0, min(0.9d0, 0.9d0*frac))
               if (opts% report) write(*,'(a,1pd16.8,a,1pd16.8,a,1pd16.8)') &
                  ' adv_diff: step rejected, max change ', change_total, ' > max_abs_change ', &
                  opts% max_abs_change, ', new dt = ', dt
               if (dt < dt_min) then
                  ierr = adv_diff_err_dt_too_small
                  stats% dt_last = dt
                  write(stats% message, '(a,1pd16.8)') &
                     'adv_diff_integrate: dt below dt_min_factor*tau_min while enforcing max_abs_change = ', &
                     opts% max_abs_change
                  return
               end if
               cycle step_loop
            end if
         end if

         ! accept the step
         x(1:nz) = x_new(1:nz)
         stats% steps = stats% steps + 1
         stats% dt_last = dt
         if (last_step) then
            time = total_time
         else
            time = time + dt
         end if
         stats% time = time
         stats% max_change_last = change_last
         stats% max_change_total = change_total

         call face_fluxes(p, x, flux_adv, flux_diff)

         change_other = 0d0
         if (present(after_step)) then
            ierr_cb = 0
            call after_step(dt, x, flux_adv, flux_diff, change_other, ierr_cb)
            if (ierr_cb /= 0) then
               ierr = adv_diff_err_callback
               write(stats% message, '(a,i0,a,i0)') &
                  'adv_diff_integrate: after_step callback returned ierr = ', ierr_cb, ' at step ', stats% steps
               return
            end if
            if (is_bad(change_other) .or. change_other < 0d0) change_other = 0d0
         end if

         ! dimensionless residual: largest change per unit time, in units of the slowest
         ! face timescale, of x and of whatever the callback advances
         res = max(change_last, change_other)/dt*tau_max
         stats% residual = res

         if (opts% report) write(*,'(a,i8,a,i8,a,1pd16.8,a,1pd16.8,a,1pd16.8,a,1pd16.8)') &
            ' adv_diff: step ', stats% steps, '  newton_iters ', stats% newton_iters, &
            '  time ', time, '  dt ', dt, '  residual ', res, '  max_change ', change_last

         ! stopping criteria. The steady test is only trusted once dt has grown past the
         ! slowest face timescale: a tiny early step can show no change at round-off level
         ! although slow modes are still far from equilibrium.
         if (steady) then
            if (res <= opts% tol_steady .and. dt >= tau_max) then
               stats% stop_reason = adv_diff_stop_steady_state
               return
            end if
         else if (last_step) then
            stats% stop_reason = adv_diff_stop_total_time
            return
         else if (res <= opts% tol_steady .and. dt >= tau_max) then
            ! time-accurate mode, steady state (of x and of what the callback advances)
            ! reached before total_time: the coefficients are frozen during the integration,
            ! so the rest of the interval changes nothing
            stats% stop_reason = adv_diff_stop_steady_state
            stats% time = total_time
            return
         end if
         if (opts% max_abs_change > 0d0) then
            if (change_total >= change_limit_fraction*opts% max_abs_change) then
               stats% stop_reason = adv_diff_stop_change_limit
               return
            end if
         end if

         ! switched evolution relaxation: grow dt by the residual ratio, but by at least
         ! dt_growth_min (a residual plateau, typical while a slow mode relaxes at constant
         ! rate after the fast modes have equilibrated, must not stall the continuation)
         if (res > 0d0) then
            growth = min(opts% dt_growth_max, max(opts% dt_growth_min, res_prev/res))
         else
            growth = opts% dt_growth_max
         end if
         res_prev = res
         dt = min(dt*growth, dt_max)

      end do step_loop

   end subroutine adv_diff_integrate


   function adv_diff_error_message(ierr) result(msg)
      integer, intent(in) :: ierr
      character(len=64) :: msg
      select case (ierr)
      case (adv_diff_ok)
         msg = 'ok'
      case (adv_diff_err_bad_input)
         msg = 'bad input (sizes, NaN, signs, boundary flux, options)'
      case (adv_diff_err_singular)
         msg = 'singular pivot in tridiagonal solve'
      case (adv_diff_err_newton)
         msg = 'Newton iteration did not converge'
      case (adv_diff_err_no_steady_state)
         msg = 'max_steps reached before steady state / total_time'
      case (adv_diff_err_conservation)
         msg = 'mass conservation violated'
      case (adv_diff_err_bounds)
         msg = 'solution out of bounds'
      case (adv_diff_err_dt_too_small)
         msg = 'time step shrank below dt_min_factor*tau_min'
      case (adv_diff_err_callback)
         msg = 'after_step callback failed'
      case default
         msg = 'unknown advection_diffusion error code'
      end select
   end function adv_diff_error_message

end module advection_diffusion
