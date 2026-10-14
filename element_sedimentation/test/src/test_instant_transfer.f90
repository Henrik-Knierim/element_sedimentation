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

module test_instant_transfer

   ! unit tests for adv_diff_instant_transfer: the steady state of the
   ! floor/cap transport for infinitely fast advection (the 'instant_rain'
   ! method), including the growth of a layer at the cap (helium ocean)

   use const_def, only: dp
   use advection_diffusion
   use test_support, only: assert_true, assert_close, assert_equal_int

   implicit none

   private
   public :: run_instant_transfer_tests

   real(dp), parameter :: rtol = 1d-12, atol = 1d-300

contains

   subroutine run_instant_transfer_tests(nfail)
      integer, intent(inout) :: nfail
      write(*,'(a)') 'test_instant_transfer'
      call test_pure_mixing(nfail)
      call test_envelope_gap_sink(nfail)
      call test_undersaturated_envelope(nfail)
      call test_deposit_block_with_floors(nfail)
      call test_ocean_growth(nfail)
      call test_no_room(nfail)
      call test_bad_input(nfail)
      call test_fluxes_match_changes(nfail)
      call test_floors_without_rain_faces(nfail)
      call test_exact_fill_to_cap(nfail)
   end subroutine run_instant_transfer_tests


   subroutine test_floors_without_rain_faces(nfail)
      ! floors act only through outgoing advection: with no rain face at all, a mixed
      ! region with floors below its mean is homogenized, not clamped
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 3
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), x0(nz), flux(nz+1), dm(nz)
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr
      dm = 1d0
      mixed = .true.; mixed(1) = .false.; mixed(nz+1) = .false.
      rain = .false.
      call setup(p, nz, dm, mixed, rain, [0.20d0, 0.22d0, 0.24d0], [0.9d0, 0.9d0, 0.9d0])
      x0 = 0.30d0
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 10, st, ierr)
      call assert_equal_int('no rain faces: ierr', ierr, 0, nfail)
      call assert_true('no rain faces: uniform at the mean', all(abs(x - 0.30d0) < 1d-14), nfail)
      call assert_true('no rain faces: zero fluxes', all(abs(flux) < 1d-14), nfail)
      call free_adv_diff_problem(p)
   end subroutine test_floors_without_rain_faces


   subroutine test_exact_fill_to_cap(nfail)
      ! a cell filled to its cap is exactly at the cap (cap - x_old not representable)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 2
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), flux(nz+1), dm(nz)
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr
      dm = [1d0, 0.3d0]
      mixed = .false.
      rain = .false.; rain(2) = .true.
      call setup(p, nz, dm, mixed, rain, [0.05d0, adv_diff_huge], [1d0, 0.45d0])
      x = [0.5d0, 0.1d0]
      call adv_diff_instant_transfer(p, x, flux, 10, st, ierr)
      call assert_equal_int('exact cap: ierr', ierr, 0, nfail)
      call assert_true('exact cap: sink exactly at its cap', x(2) == 0.45d0, nfail)
      call assert_close('exact cap: rest back-filled into the source', x(1), 0.5d0 - 0.3d0*(0.45d0 - 0.1d0), rtol, atol, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_exact_fill_to_cap


   subroutine setup(p, nz, dm, mixed, rain, floors, caps)
      ! problem with cell masses dm, mixed faces (sig = 1 across face k if mixed(k)),
      ! rain faces (w = 1 across face k if rain(k)), floors and caps per cell
      type(adv_diff_problem), intent(out) :: p
      integer, intent(in) :: nz
      real(dp), intent(in) :: dm(:), floors(:), caps(:)
      logical, intent(in) :: mixed(:), rain(:)   ! (nz+1)
      integer :: k, ierr
      call alloc_adv_diff_problem(p, nz, ierr)
      p% dm(1:nz) = dm(1:nz)
      do k = 2, nz
         if (mixed(k)) p% sig(k) = 1d0
         if (rain(k)) p% w(k) = 1d0
      end do
      p% x_floor(1:nz) = floors(1:nz)
      p% x_cap(1:nz) = caps(1:nz)
   end subroutine setup


   subroutine check_conservation(name, p, x0, x, flux, nfail)
      character(len=*), intent(in) :: name
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: x0(:), x(:), flux(:)
      integer, intent(inout) :: nfail
      integer :: k, nz
      nz = p% nz
      call assert_close(name // ': mass conserved', dot_product(p% dm(1:nz), x(1:nz)), &
         dot_product(p% dm(1:nz), x0(1:nz)), rtol, atol, nfail)
      call assert_true(name // ': boundary fluxes zero', flux(1) == 0d0 .and. flux(nz+1) == 0d0, nfail)
      do k = 1, nz
         call assert_close(name // ': dm dx = flux in - flux out', p% dm(k)*(x(k) - x0(k)), flux(k) - flux(k+1), &
            1d-12, 1d-14*p% dm(k), nfail)
      end do
      call assert_true(name // ': non-negative', all(x(1:nz) >= 0d0), nfail)
   end subroutine check_conservation


   subroutine test_pure_mixing(nfail)
      ! no rain face: runs of mixed cells are homogenized, unmixed cells untouched
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 6
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), x0(nz), flux(nz+1), dm(nz)
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr
      dm = [1d0, 2d0, 1d0, 1d0, 3d0, 1d0]
      mixed = [.false., .true., .true., .false., .false., .true., .false.]   ! faces 2,3 and 6 mixed
      rain = .false.
      call setup(p, nz, dm, mixed, rain, [(adv_diff_huge, ierr=1,nz)], [(1d0, ierr=1,nz)])
      x0 = [0.3d0, 0.2d0, 0.1d0, 0.5d0, 0.2d0, 0.4d0]
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 10, st, ierr)
      call assert_equal_int('mixing: ierr', ierr, 0, nfail)
      call assert_close('mixing: cells 1-3 mean', x(1), (0.3d0 + 0.4d0 + 0.1d0)/4d0, rtol, atol, nfail)
      call assert_true('mixing: cells 1-3 uniform', x(1) == x(2) .and. x(2) == x(3), nfail)
      call assert_close('mixing: cell 4 untouched', x(4), 0.5d0, rtol, atol, nfail)
      call assert_close('mixing: cells 5-6 mean', x(5), (0.6d0 + 0.4d0)/4d0, rtol, atol, nfail)
      call assert_true('mixing: cells 5-6 uniform', x(5) == x(6), nfail)
      call check_conservation('mixing', p, x0, x, flux, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_pure_mixing


   subroutine test_envelope_gap_sink(nfail)
      ! mixed envelope (no gap) above three gap cells with rising floors, above a mixed
      ! sink region: the envelope drains to the top floor, gap cells sit at their floors,
      ! the sink region takes the rest and homogenizes
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 8
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), x0(nz), flux(nz+1), dm(nz), floors(nz), caps(nz), m_sink, he_sink
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr, k
      dm = [2d0, 2d0, 1d0, 1d0, 1d0, 3d0, 3d0, 4d0]
      mixed = .true.            ! everything convective
      mixed(1) = .false.; mixed(nz+1) = .false.
      rain = .false.
      rain(4) = .true.; rain(5) = .true.; rain(6) = .true.   ! outgoing faces of the gap cells 3, 4, 5
      floors = adv_diff_huge
      floors(3:5) = [0.20d0, 0.22d0, 0.24d0]
      caps = 1d0
      caps(3:5) = 0.9d0
      call setup(p, nz, dm, mixed, rain, floors, caps)
      x0 = 0.28d0
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 10, st, ierr)
      call assert_equal_int('envelope: ierr', ierr, 0, nfail)
      call assert_close('envelope: cell 1 at top floor', x(1), 0.20d0, rtol, atol, nfail)
      call assert_close('envelope: cell 2 at top floor', x(2), 0.20d0, rtol, atol, nfail)
      do k = 3, 5
         call assert_close('envelope: gap cell at its floor', x(k), floors(k), rtol, atol, nfail)
      end do
      m_sink = sum(dm(6:8))
      he_sink = dot_product(dm, x0) - 0.20d0*(dm(1) + dm(2)) - dot_product(dm(3:5), floors(3:5))
      do k = 6, 8
         call assert_close('envelope: sink uniform', x(k), he_sink/m_sink, rtol, atol, nfail)
      end do
      call assert_true('envelope: fluxes toward the centre through the gap (positive)', all(flux(4:6) > 0d0), nfail)
      call assert_true('envelope: few passes', st% steps <= 3, nfail)
      call check_conservation('envelope', p, x0, x, flux, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_envelope_gap_sink


   subroutine test_undersaturated_envelope(nfail)
      ! gap cells above their floors drain; whether the envelope above them takes part
      ! depends on diffusion: with convection the whole region settles below the floors
      ! (helium can diffuse out of a cell that sits at its floor), without convection the
      ! envelope keeps what it has and the sink takes the drained helium
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 5
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), x0(nz), flux(nz+1), dm(nz), floors(nz), caps(nz)
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr
      dm = 1d0
      mixed = .true.; mixed(1) = .false.; mixed(nz+1) = .false.
      rain = .false.; rain(3) = .true.; rain(4) = .true.
      floors = adv_diff_huge; floors(2:3) = [0.25d0, 0.27d0]
      caps = 1d0
      call setup(p, nz, dm, mixed, rain, floors, caps)
      x0 = [0.20d0, 0.26d0, 0.30d0, 0.20d0, 0.20d0]
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 10, st, ierr)
      call assert_equal_int('under: ierr', ierr, 0, nfail)
      ! 1.16 over five cells = 0.232, below both floors: uniform, nothing is held
      call assert_true('under (mixed): uniform', all(abs(x - 0.232d0) < 1d-14), nfail)
      call check_conservation('under (mixed)', p, x0, x, flux, nfail)
      ! the same without diffusion
      mixed = .false.
      call free_adv_diff_problem(p)
      call setup(p, nz, dm, mixed, rain, floors, caps)
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 10, st, ierr)
      call assert_equal_int('under (unmixed): ierr', ierr, 0, nfail)
      call assert_close('under (unmixed): cell 2 at floor', x(2), 0.25d0, rtol, atol, nfail)
      call assert_close('under (unmixed): cell 3 at floor', x(3), 0.27d0, rtol, atol, nfail)
      call assert_close('under (unmixed): envelope untouched', x(1), 0.20d0, rtol, atol, nfail)
      call assert_close('under (unmixed): first sink cell takes it', x(4), 0.24d0, rtol, atol, nfail)
      call assert_close('under (unmixed): last cell untouched', x(5), 0.20d0, rtol, atol, nfail)
      call check_conservation('under (unmixed)', p, x0, x, flux, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_undersaturated_envelope


   subroutine test_deposit_block_with_floors(nfail)
      ! a mixed region of cells that have a gap but are below their floors receives more
      ! than it can hold uniformly: its top cells end at their floors, the rest uniform
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 6
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), x0(nz), flux(nz+1), dm(nz), floors(nz), caps(nz)
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr
      dm = 1d0
      mixed = .true.; mixed(1) = .false.; mixed(nz+1) = .false.
      rain = .false.; rain(2:6) = .true.    ! every cell but the last advects
      floors = [0.10d0, 0.30d0, 0.32d0, 0.34d0, 0.36d0, adv_diff_huge]
      caps = 1d0
      call setup(p, nz, dm, mixed, rain, floors, caps)
      ! cell 1 holds a lot above its floor; cells 2-5 are below their floors; cell 6 is a sink
      x0 = [0.60d0, 0.28d0, 0.28d0, 0.28d0, 0.28d0, 0.28d0]
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 20, st, ierr)
      call assert_equal_int('deposit: ierr', ierr, 0, nfail)
      call assert_close('deposit: source at floor', x(1), 0.10d0, rtol, atol, nfail)
      ! 0.5 of helium mass rains out of cell 1; the five cells below hold 5*0.28 + 0.5 = 1.9,
      ! i.e. 0.38 on average, above every floor: cells 2-5 end at their floors (1.32), the
      ! sink gets the rest (0.58)
      call assert_close('deposit: cell 2 at floor', x(2), 0.30d0, rtol, atol, nfail)
      call assert_close('deposit: cell 3 at floor', x(3), 0.32d0, rtol, atol, nfail)
      call assert_close('deposit: cell 4 at floor', x(4), 0.34d0, rtol, atol, nfail)
      call assert_close('deposit: cell 5 at floor', x(5), 0.36d0, rtol, atol, nfail)
      call assert_close('deposit: sink holds the rest', x(6), 1.9d0 - 1.32d0, rtol, atol, nfail)
      call check_conservation('deposit', p, x0, x, flux, nfail)
      ! less helium: only the top cells of the block reach their floors, the rest stays uniform
      x0 = [0.30d0, 0.28d0, 0.28d0, 0.28d0, 0.28d0, 0.28d0]
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 20, st, ierr)
      call assert_equal_int('deposit2: ierr', ierr, 0, nfail)
      ! 0.2 rains out; 1.6 over five cells = 0.32 > floor of cell 2 (0.30): cell 2 at its
      ! floor; 1.3 over four = 0.325 > floor of cell 3 (0.32): cell 3 at its floor; 0.98 over
      ! three = 0.3267 < floor of cell 4 (0.34): cells 4-6 stay uniform
      call assert_close('deposit2: cell 2 at floor', x(2), 0.30d0, rtol, atol, nfail)
      call assert_close('deposit2: cell 3 at floor', x(3), 0.32d0, rtol, atol, nfail)
      call assert_close('deposit2: cells 4-6 uniform (4)', x(4), 0.98d0/3d0, 1d-12, atol, nfail)
      call assert_close('deposit2: cells 4-6 uniform (5)', x(5), 0.98d0/3d0, 1d-12, atol, nfail)
      call assert_close('deposit2: cells 4-6 uniform (6)', x(6), 0.98d0/3d0, 1d-12, atol, nfail)
      call check_conservation('deposit2', p, x0, x, flux, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_deposit_block_with_floors


   subroutine test_ocean_growth(nfail)
      ! the gap reaches the centre: the innermost cell fills to its cap, then the cells
      ! above it, one after the other; the rest of the column sits at its floors
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 6
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), x0(nz), flux(nz+1), dm(nz), floors(nz), caps(nz), excess
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr
      dm = 1d0
      mixed = .false.                       ! no convection: isolates the drain logic
      rain = .false.; rain(2:6) = .true.    ! cells 1-5 advect; cell 6 (centre) is the sink
      floors = [0.20d0, 0.22d0, 0.24d0, 0.26d0, 0.28d0, 0.30d0]
      caps = [0.80d0, 0.82d0, 0.84d0, 0.86d0, 0.88d0, 0.90d0]
      call setup(p, nz, dm, mixed, rain, floors, caps)
      ! every cell well inside the gap with a lot of helium: 6 x 0.6 = 3.6 of helium. The
      ! floors of cells 1-3 hold 0.66, cells 6, 5 and 4 fill to their caps (2.64), and the
      ! remaining 0.30 stays in cell 3 on top of its floor: the top of the ocean
      x0 = 0.60d0
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 20, st, ierr)
      call assert_equal_int('ocean: ierr', ierr, 0, nfail)
      call assert_close('ocean: centre at cap', x(6), 0.90d0, rtol, atol, nfail)
      call assert_close('ocean: cell 5 at cap', x(5), 0.88d0, rtol, atol, nfail)
      call assert_close('ocean: cell 4 at cap', x(4), 0.86d0, rtol, atol, nfail)
      excess = 3.6d0 - 0.90d0 - 0.88d0 - 0.86d0 - (0.20d0 + 0.22d0)
      call assert_close('ocean: cell 3 partly filled', x(3), excess, rtol, atol, nfail)
      call assert_true('ocean: cell 3 between floor and cap', x(3) > 0.24d0 .and. x(3) < 0.84d0, nfail)
      call assert_close('ocean: cell 1 at floor', x(1), 0.20d0, rtol, atol, nfail)
      call assert_close('ocean: cell 2 at floor', x(2), 0.22d0, rtol, atol, nfail)
      call check_conservation('ocean', p, x0, x, flux, nfail)
      ! a second call on the result changes nothing (cell 3 is inside its gap but is the
      ! top of the ocean: it drains into the full cells below and is refilled, net zero)
      x0 = x
      call adv_diff_instant_transfer(p, x, flux, 20, st, ierr)
      call assert_equal_int('ocean: second call ierr', ierr, 0, nfail)
      call assert_true('ocean: second call is a no-op', all(abs(x - x0) < 1d-14), nfail)
      call assert_true('ocean: second call zero fluxes', all(abs(flux) < 1d-14), nfail)
      call free_adv_diff_problem(p)
   end subroutine test_ocean_growth


   subroutine test_no_room(nfail)
      ! a drain that is blocked by full cells below returns to the same state
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 3
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), x0(nz), flux(nz+1), dm(nz)
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr
      dm = 1d0
      mixed = .false.
      rain = .false.; rain(2:3) = .true.
      call setup(p, nz, dm, mixed, rain, [0.1d0, 0.1d0, 0.1d0], [0.5d0, 0.5d0, 0.5d0])
      x0 = [0.49d0, 0.5d0, 0.5d0]
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 20, st, ierr)
      call assert_equal_int('blocked: ierr', ierr, 0, nfail)
      call assert_true('blocked: state unchanged', all(abs(x - x0) < 1d-14), nfail)
      call assert_true('blocked: zero fluxes', all(abs(flux) < 1d-14), nfail)
      call free_adv_diff_problem(p)
      ! (adv_diff_err_no_room cannot occur when every cap is at least the floor: a draining
      ! cell always has room for what it drained, and a cell at or above its cap is fixed;
      ! the code is a guard against inconsistent limits)
   end subroutine test_no_room


   subroutine test_bad_input(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 3
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), flux(nz+1), dm(nz)
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr
      dm = 1d0
      mixed = .false.
      rain = .false.; rain(2) = .true.
      call setup(p, nz, dm, mixed, rain, [0.1d0, 0.1d0, adv_diff_huge], [1d0, 1d0, 1d0])
      p% w(2) = -1d0   ! advection toward the surface is not supported
      x = 0.3d0
      call adv_diff_instant_transfer(p, x, flux, 20, st, ierr)
      call assert_equal_int('bad input: negative w rejected', ierr, adv_diff_err_bad_input, nfail)
      p% w(2) = 1d0
      x = [0.3d0, -0.1d0, 0.3d0]
      call adv_diff_instant_transfer(p, x, flux, 20, st, ierr)
      call assert_equal_int('bad input: negative x rejected', ierr, adv_diff_err_bad_input, nfail)
      x = 0.3d0
      call adv_diff_instant_transfer(p, x, flux, 0, st, ierr)
      call assert_equal_int('bad input: max_passes < 1 rejected', ierr, adv_diff_err_bad_input, nfail)
      call free_adv_diff_problem(p)
   end subroutine test_bad_input


   subroutine test_fluxes_match_changes(nfail)
      ! random problem: the returned face fluxes reproduce the change of every cell
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 12
      type(adv_diff_problem) :: p
      type(adv_diff_stats) :: st
      real(dp) :: x(nz), x0(nz), flux(nz+1), dm(nz), floors(nz), caps(nz)
      logical :: mixed(nz+1), rain(nz+1)
      integer :: ierr, k
      do k = 1, nz
         dm(k) = 1d0 + 0.1d0*k
         floors(k) = 0.2d0 + 0.01d0*k
         caps(k) = 0.7d0 + 0.01d0*k
      end do
      floors(1:2) = adv_diff_huge
      floors(9:12) = adv_diff_huge
      caps(1:2) = 1d0
      caps(9:12) = 1d0
      mixed = .false.
      mixed(2) = .true.; mixed(10:12) = .true.
      rain = .false.
      rain(4:9) = .true.
      call setup(p, nz, dm, mixed, rain, floors, caps)
      do k = 1, nz
         x0(k) = 0.25d0 + 0.02d0*mod(3*k, 7)
      end do
      x = x0
      call adv_diff_instant_transfer(p, x, flux, 20, st, ierr)
      call assert_equal_int('fluxes: ierr', ierr, 0, nfail)
      call check_conservation('fluxes', p, x0, x, flux, nfail)
      ! through the unmixed rain faces only advection acts: nothing moves outward there
      call assert_true('fluxes: nothing moves outward through the rain faces', all(flux(4:9) >= 0d0), nfail)
      call free_adv_diff_problem(p)
   end subroutine test_fluxes_match_changes

end module test_instant_transfer
