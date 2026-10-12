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

module element_sedimentation_legacy

   ! Legacy element sedimentation methods 'instant_rain' and 'bottom_up_rain'.
   !
   ! Both methods move the sedimenting species (index i_Xi in s% xa) out of
   ! the cells where it is immiscible and deposit it in the region below in
   ! one mass-conserving transfer per call; no transport equation is solved.
   !
   !   instant_rain    one global transfer over the cells 1:nz
   !   bottom_up_rain  the transfer is repeated for the ranges k:nz with
   !                   k = nz, ..., 1, so that regions which start to rain
   !                   only because of the material deposited below them are
   !                   handled within the same call
   !
   ! The miscibility limits come from sedimentation_miscibility
   ! (evaluate_miscibility). The arrays used by the algorithms are those of
   ! the original implementation:
   !
   !   is_miscible(k)     = .not. in_gap(k)
   !   Xi_miscible_min(k) = x_lo(k) where the phase diagram has a usable gap, -1d0 elsewhere
   !   Xi_miscible_max(k) = x_hi(k) where the phase diagram has a usable gap, -1d0 elsewhere
   !
   ! A gap is usable when x_lo(k) < x_hi(k); a cell without hydrogen and helium
   ! (e.g. a pure-metal core) has x_lo = x_hi = 0 and is treated as having no
   ! gap (in the old code such cells produced NaN limits, which were inert).
   !
   ! The algorithms themselves are unchanged. This module is a mechanical
   ! move of the old monolithic element_sedimentation module; the
   ! advection-diffusion method lives in advection_diffusion.f90 and is
   ! driven by element_sedimentation.f90. Results are not bit-for-bit those of
   ! the old code because the miscibility limits are now evaluated by
   ! sedimentation_miscibility: with the total EOS pressure s% Peos (old:
   ! exp(lnPgas)), with atomic weights chem_isos% W (old: integer Z+N, a
   ! ~4d-4 relative shift of the limits), and with the He-rich limit computed
   ! from the original composition (the old code rescaled the composition in
   ! place for the He-poor limit first, which biased the He-rich limit high by
   ! ~6d-4 relative at Z = 0.1).
   !
   ! Mixing: 'instant_rain' always sets D_mix = 0 in the immiscible cells and
   ! switches off do_conv_premix (old behaviour); 'bottom_up_rain' never
   ! changes the mixing. Neither method uses the control
   ! element_sedimentation_suppress_mixing_in_immiscible_regions.
   !
   ! Failures (mass not conserved, mass fraction outside [0,1]) return ierr /= 0
   ! with a message, so that MESA retries the step (the old code stopped with
   ! mesa_error).
   !
   ! Index conventions: arrays named *_monotonic and Xi_temp are allocated
   ! over the working range (klo:khi) and are passed as allocatable dummies
   ! so that their bounds are preserved. check_abundance_array and friends
   ! receive the slice of cells klo:khi (element 1 is cell klo).

   use star_private_def
   use const_def, only: dp, no_mixing
   use utils_lib, only: mesa_error, is_bad
   use sedimentation_phase_diagram, only: phase_diagram, free_phase_diagram
   use sedimentation_miscibility, only: evaluate_miscibility, miscibility_error_message

   implicit none

   private
   public :: init_element_sedimentation_legacy
   public :: shutdown_element_sedimentation_legacy
   public :: do_instant_rain
   public :: do_bottom_up_rain

   ! state set by init_element_sedimentation_legacy
   type(phase_diagram), save :: phase_diag    ! private copy of the immiscibility phase diagram
   integer :: i_Xi = -1                       ! index of the sedimenting species in s% xa
   integer :: i_Xp = -1                       ! index of its partner species (other component of the phase diagram)
   real(dp), allocatable :: atomic_weight(:)  ! (species)
   logical :: initialized = .false.

   ! flags refreshed from the controls at every call
   logical :: verbose = .false.             ! follows s% element_sedimentation_report
   logical :: damp_sedimentation = .false.  ! s% element_sedimentation_maximum_dXi_per_cell > 0

   logical, parameter :: debug = .false.
   logical, parameter :: debug_profile = .false.  ! dump the miscibility arrays to files (print_miscibility_info)
   logical, parameter :: homogenize_regions_after_transfer = .false.  ! was disabled ('if (.false.)') in the old code

   real(dp), parameter :: small = 1d-3  ! relative tolerance of the conservation and monotonicity checks

   ! directions for enforce_monotonicity
   integer, parameter :: enforce_monotonicity_bottom_up = -1
   integer, parameter :: enforce_monotonicity_top_down = 1

   ! monotonicity methods for check_abundance_array
   integer, parameter :: monotonic_increase = 1
   integer, parameter :: no_monotonicity = 0
   integer, parameter :: monotonic_decrease = -1

contains

   ! ----------------------------------------------------------------------
   ! initialization
   ! ----------------------------------------------------------------------

   subroutine init_element_sedimentation_legacy(pd, i_sed, i_partner, A, ierr)
      ! stores a private copy of the phase diagram, the species indices and the atomic weights
      type(phase_diagram), intent(in) :: pd
      integer, intent(in) :: i_sed      ! index of the sedimenting species in s% xa
      integer, intent(in) :: i_partner  ! index of the partner species (other component of the phase diagram)
      real(dp), intent(in) :: A(:)      ! (species) atomic weights
      integer, intent(out) :: ierr
      include 'formats'

      ierr = 0

      if (.not. pd% loaded) then
         write(*,'(a)') 'init_element_sedimentation_legacy: phase diagram is not loaded'
         ierr = 1
         return
      end if
      if (i_sed < 1 .or. i_sed > size(A) .or. i_partner < 1 .or. i_partner > size(A) .or. i_sed == i_partner) then
         write(*,'(a)') 'init_element_sedimentation_legacy: invalid species indices'
         write(*,11) 'i_sed, i_partner, size(A)', i_sed, i_partner, size(A)
         ierr = 1
         return
      end if
      if (any(is_bad(A)) .or. any(A <= 0d0)) then
         write(*,'(a)') 'init_element_sedimentation_legacy: atomic weights must be finite and > 0'
         write(*,1) 'A', A
         ierr = 1
         return
      end if

      phase_diag = pd  ! deep copy (allocatable components)
      i_Xi = i_sed
      i_Xp = i_partner
      if (allocated(atomic_weight)) deallocate(atomic_weight)
      allocate(atomic_weight(size(A)))
      atomic_weight(:) = A(:)
      initialized = .true.

   end subroutine init_element_sedimentation_legacy


   subroutine shutdown_element_sedimentation_legacy()
      ! frees the private copy of the phase diagram
      call free_phase_diagram(phase_diag)
      if (allocated(atomic_weight)) deallocate(atomic_weight)
      i_Xi = -1
      i_Xp = -1
      initialized = .false.
   end subroutine shutdown_element_sedimentation_legacy


   subroutine check_initialized(caller, ierr)
      character(len=*), intent(in) :: caller
      integer, intent(out) :: ierr
      ierr = 0
      if (.not. initialized) then
         write(*,'(a)') trim(caller) // ': init_element_sedimentation_legacy has not been called'
         ierr = 1
      end if
   end subroutine check_initialized


   subroutine set_call_time_flags(s)
      ! module flags that follow the controls
      type(star_info), pointer :: s
      verbose = s% element_sedimentation_report
      damp_sedimentation = (s% element_sedimentation_maximum_dXi_per_cell > 0d0)
   end subroutine set_call_time_flags


   ! ----------------------------------------------------------------------
   ! the two legacy methods
   ! ----------------------------------------------------------------------

   subroutine do_instant_rain(s, ierr)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr

      real(dp), allocatable :: xa_pre(:,:)
      real(dp), allocatable :: dXi(:)
      logical, allocatable :: is_miscible(:)
      real(dp), allocatable :: Xi_miscible_min(:), Xi_miscible_max(:)

      ierr = 0
      call check_initialized('do_instant_rain', ierr)
      if (ierr /= 0) return
      call set_call_time_flags(s)

      ! 1. snapshot of the composition before the rain
      allocate(xa_pre(size(s% xa, dim=1), size(s% xa, dim=2)))
      xa_pre = s% xa

      ! 2. global equilibrium change of the mass fraction of species i_Xi
      if (debug) write(*,'(a)') 'do_instant_rain: calling get_dXi'
      call get_dXi(s, 1, s% nz, dXi, is_miscible, Xi_miscible_min, Xi_miscible_max, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'do_instant_rain: error in get_dXi'
         return
      end if

      if (.not. allocated(dXi)) then
         if (debug) write(*,'(a)') 'do_instant_rain: no precipitating region found.'
         return
      end if

      ! 3. apply the change to s% xa and let the other species counter-flow
      if (debug) write(*,'(a)') 'do_instant_rain: calling adjust_mass_fractions'
      call adjust_mass_fractions(s, dXi, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'do_instant_rain: error in adjust_mass_fractions'
         return
      end if

      ! 4. switch off mixing in the immiscible cells so that MESA does not undo the rain immediately
      if (debug) write(*,'(a)') 'do_instant_rain: calling adjust_mixing_in_sedimenting_region'
      call adjust_mixing_in_sedimenting_region(s, is_miscible, ierr)
      if (ierr /= 0) return

      ! 5. verification
      if (debug) write(*,'(a)') 'do_instant_rain: calling check_after_rain'
      call check_after_rain(s, 1, s% nz, xa_pre, ierr)

   end subroutine do_instant_rain


   subroutine do_bottom_up_rain(s, ierr)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr

      integer :: k_update, k_debug
      real(dp), allocatable :: xa_pre(:,:)
      real(dp), allocatable :: dXi(:)
      logical, allocatable :: is_miscible(:)
      real(dp), allocatable :: Xi_miscible_min(:), Xi_miscible_max(:)
      include 'formats'

      ierr = 0
      call check_initialized('do_bottom_up_rain', ierr)
      if (ierr /= 0) return
      call set_call_time_flags(s)

      allocate(xa_pre(size(s% xa, dim=1), size(s% xa, dim=2)))
      xa_pre = s% xa

      ! A single global transfer can make a region rain out that did not
      ! rain before (because of the material deposited below it). Hence the
      ! transfer is repeated for the ranges k_update:nz, from the centre
      ! outwards, until the whole model has been covered.
      do k_update = s% nz, 1, -1

         if (debug) write(*,'(a,i0)') 'do_bottom_up_rain: k_update = ', k_update

         call get_dXi(s, k_update, s% nz, dXi, is_miscible, Xi_miscible_min, Xi_miscible_max, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'do_bottom_up_rain: error in get_dXi'
            return
         end if

         if (.not. allocated(dXi)) cycle

         if (debug) then
            k_debug = maxloc(abs(dXi), dim=1)
            write(*,3) 'do_bottom_up_rain: k_update, k_max, dXi_max', k_update, k_debug, dXi(k_debug)
         end if

         call adjust_mass_fractions(s, dXi, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'do_bottom_up_rain: error in adjust_mass_fractions'
            return
         end if

         call check_after_rain(s, 1, s% nz, xa_pre, ierr)
         if (ierr /= 0) return

         deallocate(dXi, is_miscible, Xi_miscible_min, Xi_miscible_max)

      end do

   end subroutine do_bottom_up_rain


   ! ----------------------------------------------------------------------
   ! change of the mass fraction of the sedimenting species
   ! ----------------------------------------------------------------------

   subroutine get_dXi(s, klo, khi, dXi, is_miscible, Xi_miscible_min, Xi_miscible_max, ierr)
      ! dXi stays unallocated when fewer than
      ! s% element_sedimentation_cell_threshold_for_precipitation cells are immiscible
      type(star_info), pointer :: s
      integer, intent(in) :: klo, khi  ! range of cells to consider
      real(dp), allocatable, intent(out) :: dXi(:)  ! (1:nz) change of the mass fraction of species i_Xi; zero outside klo:khi
      logical, allocatable, intent(out) :: is_miscible(:)  ! (1:nz)
      real(dp), allocatable, intent(out) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! (1:nz)
      integer, intent(out) :: ierr

      integer :: k, cell_threshold

      if (debug) write(*,'(a)') 'get_dXi: calling get_raw_miscibility'
      call get_raw_miscibility(s, klo, khi, is_miscible, Xi_miscible_min, Xi_miscible_max, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'get_dXi: error in get_raw_miscibility'
         return
      end if

      if (debug_profile) call print_miscibility_info(s, is_miscible, Xi_miscible_min, Xi_miscible_max, &
         .false., 'raw_miscibility_info.dat')

      ! if (nearly) everything is miscible, there is nothing to do
      cell_threshold = s% element_sedimentation_cell_threshold_for_precipitation
      if (count(.not. is_miscible) < cell_threshold) then
         if (verbose) write(*,'(a)') 'get_dXi: no precipitating region found, returning.'
         deallocate(Xi_miscible_min, Xi_miscible_max, is_miscible)
         return
      end if

      if (debug) then
         write(*,'(a)') 'get_dXi: immiscible cells:'
         write(*,'(20i8)') pack([(k, k=1, s% nz)], .not. is_miscible)
      end if

      ! difference between the miscible and the actual mass fraction in the precipitating region
      if (debug) write(*,'(a)') 'get_dXi: calling set_dXi'
      call set_dXi(s, is_miscible, Xi_miscible_min, Xi_miscible_max, klo, khi, dXi, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'get_dXi: error in set_dXi, bypassing mass conservation check'
         return
      end if

      if (debug) write(*,'(a)') 'get_dXi: checking mass conservation'
      if (.not. is_mass_conserved(dXi(klo:khi), s% dm(klo:khi))) then
         write(*,'(a)') 'get_dXi: mass is not conserved after set_dXi. Triggering timestep retry.'
         ierr = 1
         return
      end if

      if (debug) write(*,'(a)') 'get_dXi: finished.'

   end subroutine get_dXi


   subroutine get_raw_miscibility(s, klo, khi, is_miscible, Xi_miscible_min, Xi_miscible_max, ierr)
      ! miscibility of the cells klo:khi without any post-processing.
      ! The arrays are allocated 1:nz; cells outside klo:khi are marked miscible with limits -1d0.
      type(star_info), pointer :: s
      integer, intent(in) :: klo, khi
      logical, allocatable, intent(out) :: is_miscible(:)
      real(dp), allocatable, intent(out) :: Xi_miscible_min(:), Xi_miscible_max(:)
      integer, intent(out) :: ierr

      logical, allocatable :: has_gap(:), in_gap(:), usable(:)
      real(dp), allocatable :: x_lo(:), x_hi(:)
      integer :: n
      include 'formats'

      ierr = 0

      if (klo < 1 .or. khi > s% nz .or. klo > khi) then
         write(*,'(a)') 'get_raw_miscibility: invalid cell range'
         write(*,11) 'klo, khi, nz', klo, khi, s% nz
         ierr = 1
         return
      end if
      if (size(atomic_weight) /= s% species) then
         write(*,'(a)') 'get_raw_miscibility: number of species differs from the one given to init_element_sedimentation_legacy'
         write(*,11) 'size(atomic_weight), s% species', size(atomic_weight), s% species
         ierr = 1
         return
      end if

      allocate(Xi_miscible_min(1:s% nz), Xi_miscible_max(1:s% nz), source=-1d0)
      allocate(is_miscible(1:s% nz), source=.true.)

      n = khi - klo + 1
      allocate(has_gap(n), in_gap(n), usable(n), x_lo(n), x_hi(n))

      call evaluate_miscibility(phase_diag, n, s% species, s% xa(:,klo:khi), atomic_weight, i_Xi, i_Xp, &
         s% Peos(klo:khi), s% T(klo:khi), has_gap, in_gap, x_lo, x_hi, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'get_raw_miscibility: evaluate_miscibility failed: ' // trim(miscibility_error_message(ierr))
         write(*,11) 'klo, khi', klo, khi
         return
      end if

      ! a gap is only usable if x_lo < x_hi: a cell without H and He (e.g. a pure-metal core)
      ! has x_lo = x_hi = 0, and a minimum of 0 would make enforce_monotonicity clamp the
      ! He-poor limit of every cell above it to 0 (complete rain-out of the envelope)
      usable(1:n) = has_gap(1:n) .and. (x_lo(1:n) < x_hi(1:n))
      is_miscible(klo:khi) = .not. in_gap(1:n)
      Xi_miscible_min(klo:khi) = merge(x_lo(1:n), -1d0, usable(1:n))
      Xi_miscible_max(klo:khi) = merge(x_hi(1:n), -1d0, usable(1:n))

   end subroutine get_raw_miscibility


   subroutine set_dXi(s, is_miscible, Xi_miscible_min, Xi_miscible_max, klo, khi, dXi, ierr)
      type(star_info), pointer :: s
      logical, allocatable, intent(in) :: is_miscible(:)  ! (1:nz)
      real(dp), allocatable, intent(in) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! (1:nz)
      integer, intent(in) :: klo, khi
      real(dp), allocatable, intent(out) :: dXi(:)  ! (1:nz) change of the mass fraction of species i_Xi; zero outside klo:khi
      integer, intent(out) :: ierr

      real(dp), allocatable :: Xi_temp(:)  ! (klo:khi) working copy of the mass fractions
      real(dp), allocatable :: Xi_miscible_min_monotonic(:)  ! (klo:khi) miscibility minimum, monotonically increasing inwards
      real(dp), allocatable :: Xi_miscible_max_monotonic(:)  ! (klo:khi) miscibility maximum, monotonically increasing inwards
      integer :: k
      integer :: k_RCB        ! radiative-convective boundary (between atmosphere and convective envelope)
      integer :: k_rain       ! bottom of the precipitating region
      integer :: k_deposite   ! last cell in which material from the precipitating region is deposited
      integer :: k_rain_top   ! top of the precipitating region

      ierr = 0

      allocate(Xi_temp(klo:khi), source=s% xa(i_Xi, klo:khi))

      call init_set_dXi(s, klo, khi, is_miscible, Xi_miscible_min, Xi_miscible_max, &
         k_RCB, k_rain, k_deposite, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'set_dXi: error in init_set_dXi'
         return
      end if

      if (debug_profile) call print_miscibility_info(s, is_miscible, &
         Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, .false., 'monotonic_miscibility_info.dat')

      if (.not. allocated(Xi_miscible_min_monotonic) .or. .not. allocated(Xi_miscible_max_monotonic)) then
         write(*,'(a)') 'set_dXi: monotonic miscibility limits are not allocated'
         ierr = 1
         return
      end if

      ! basic logic: every cell from the RCB down to k_rain is reduced to its
      ! (monotonic) miscibility minimum and the excess is spread uniformly
      ! over the deposit region k_rain+1:k_deposite
      if (debug) write(*,'(a,3i8)') 'set_dXi: transfer Xi; k_RCB, k_rain, k_deposite = ', k_RCB, k_rain, k_deposite
      do k = k_RCB, k_rain
         if (debug) write(*,'(a,i0)') 'set_dXi: transferring Xi at k = ', k
         call transfer_Xi(s, Xi_temp, Xi_miscible_min_monotonic, k, k_rain+1, k_deposite, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'set_dXi: error in transfer_Xi'
            return
         end if
      end do
      if (debug) write(*,'(a)') 'set_dXi: basic transfer complete.'

      if (homogenize_regions_after_transfer) then
         ! homogenize the deposit
         call homogenize_region(s, Xi_temp, k_rain+1, k_deposite, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'set_dXi: error in homogenize_region (deposit)'
            return
         end if
         ! homogenize the convective envelope above the precipitating region
         k_rain_top = findloc(is_miscible, .false., dim=1)
         call homogenize_region(s, Xi_temp, k_RCB, k_rain_top-1, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'set_dXi: error in homogenize_region (envelope)'
            return
         end if
      end if

      ! difference between the new and the actual mass fraction in the region klo:khi
      allocate(dXi(1:s% nz), source=0d0)
      dXi(klo:khi) = Xi_temp(klo:khi) - s% xa(i_Xi, klo:khi)

      if (debug) then
         write(*,'(a)') ''
         write(*,'(a)') 'After set_dXi:'
         write(*,'(a,2i8)') 'klo, khi: ', klo, khi
         write(*,'(a)') 'k, is_miscible, mixing_type, Xi, Xi_temp, dXi, Xi_miscible_min, Xi_miscible_min_monotonic, ' // &
            'Xi_miscible_max, Xi_miscible_max_monotonic:'
         do k = klo, khi
            write(*,'(i12,l5,i12,99(1pd26.16))') k, is_miscible(k), s% mixing_type(k), s% xa(i_Xi, k), Xi_temp(k), dXi(k), &
               Xi_miscible_min(k), Xi_miscible_min_monotonic(k), Xi_miscible_max(k), Xi_miscible_max_monotonic(k)
         end do
         write(*,'(a)') ''
      end if

   end subroutine set_dXi


   subroutine init_set_dXi(s, klo, khi, is_miscible, Xi_miscible_min, Xi_miscible_max, &
         k_RCB, k_rain, k_deposite, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
      ! finds the characteristic cells and the monotonic miscibility limits for set_dXi
      type(star_info), pointer :: s
      integer, intent(in) :: klo, khi  ! region in which we transfer the element
      logical, allocatable, intent(in) :: is_miscible(:)  ! (1:nz)
      real(dp), allocatable, intent(in) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! (1:nz)
      integer, intent(out) :: k_RCB       ! radiative-convective boundary
      integer, intent(out) :: k_rain      ! bottom of the precipitating region
      integer, intent(out) :: k_deposite  ! bottom of the deposit region
      real(dp), allocatable, intent(out) :: Xi_miscible_min_monotonic(:)  ! (klo:khi)
      real(dp), allocatable, intent(out) :: Xi_miscible_max_monotonic(:)  ! (klo:khi)
      integer, intent(out) :: ierr

      call find_RCB(s, klo, khi, k_RCB, ierr)
      if (ierr /= 0) return

      ! first guess of the deposit region with k_rain = k_RCB, needed by adjust_miscibility_gap
      call find_k_deposite(s, klo, khi, k_RCB, k_deposite, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'init_set_dXi: error in find_k_deposite'
         return
      end if

      call adjust_miscibility_gap(s, klo, khi, k_RCB, k_deposite, is_miscible, Xi_miscible_min, Xi_miscible_max, &
         Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'init_set_dXi: error in adjust_miscibility_gap'
         return
      end if

      call find_k_rain(s, Xi_miscible_max_monotonic, is_miscible, k_RCB, khi, k_rain, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'init_set_dXi: error in find_k_rain'
         return
      end if

      ! now the deposit region with the correct k_rain
      call find_k_deposite(s, klo, khi, k_rain, k_deposite, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'init_set_dXi: error in find_k_deposite'
         return
      end if

      if (verbose) then
         write(*,'(a)') 'init_set_dXi: finished'
         write(*,'(a40,i12)') 'klo', klo
         write(*,'(a40,i12)') 'khi', khi
         write(*,'(a40,i12)') 'k_RCB', k_RCB
         write(*,'(a40,i12)') 'k_rain', k_rain
         write(*,'(a40,i12)') 'k_deposite', k_deposite
      end if

   end subroutine init_set_dXi


   subroutine find_RCB(s, klo, khi, k_RCB, ierr)
      ! first mixed cell in klo:khi (klo if there is none)
      type(star_info), pointer :: s
      integer, intent(in) :: klo, khi
      integer, intent(out) :: k_RCB
      integer, intent(out) :: ierr

      integer :: k
      include 'formats'

      ierr = 0
      k_RCB = -1
      do k = klo, khi
         if (s% mixing_type(k) /= no_mixing) then
            k_RCB = k
            exit
         end if
      end do

      if (k_RCB == -1) k_RCB = klo

      if (k_RCB < klo .or. k_RCB > khi) then
         write(*,'(a)') 'find_RCB: RCB out of bounds'
         write(*,11) 'klo, khi, k_RCB', klo, khi, k_RCB
         ierr = 1
         return
      end if

      if (debug) write(*,'(a,i0)') 'find_RCB: found RCB at k = ', k_RCB

   end subroutine find_RCB


   subroutine find_k_deposite(s, klo, khi, k_rain, k_deposite, ierr)
      ! bottom of the deposit region: the cell above the first unmixed cell
      ! that is separated from k_rain by more than a threshold fraction of
      ! the stellar mass (khi if there is none)
      type(star_info), pointer :: s
      integer, intent(in) :: klo, khi
      integer, intent(in) :: k_rain  ! last immiscible cell outside the core
      integer, intent(out) :: k_deposite
      integer, intent(out) :: ierr

      ! minimum mass (as a fraction of the stellar mass) between k_rain and the
      ! deposit region; smaller values make the detection more sensitive
      real(dp), parameter :: mass_separation_to_k_rain_threshold = 0.05d0
      real(dp) :: mass_separation_to_k_rain
      integer :: k
      include 'formats'

      ierr = 0
      k_deposite = -1

      do k = k_rain + 2, khi
         mass_separation_to_k_rain = sum(s% dm(k_rain+1:k-1))/s% mstar
         if (s% mixing_type(k) == no_mixing .and. mass_separation_to_k_rain > mass_separation_to_k_rain_threshold) then
            k_deposite = k - 1
            exit
         end if
      end do

      if (k_deposite == -1) k_deposite = khi

      if (k_deposite < klo .or. k_deposite > khi) then
         write(*,'(a)') 'find_k_deposite: deposit region out of bounds'
         write(*,11) 'klo, khi, k_deposite', klo, khi, k_deposite
         ierr = 1
         return
      end if

      if (debug) write(*,11) 'find_k_deposite: klo, k_rain, k_deposite, khi', klo, k_rain, k_deposite, khi

   end subroutine find_k_deposite


   subroutine find_k_rain(s, Xi_miscible_max_monotonic, is_miscible, klo, khi, k_rain, ierr)
      ! last cell in klo:khi that is immiscible outside the core; this is
      ! where the deposit region starts
      type(star_info), pointer :: s
      real(dp), allocatable, intent(in) :: Xi_miscible_max_monotonic(:)  ! (working range)
      logical, allocatable, intent(in) :: is_miscible(:)  ! (1:nz)
      integer, intent(in) :: klo, khi
      integer, intent(out) :: k_rain
      integer, intent(out) :: ierr

      integer :: k
      logical :: inside_core
      include 'formats'

      ierr = 0

      if (all(is_miscible(klo:khi))) then
         write(*,'(a)') 'find_k_rain: model is miscible everywhere in the region klo:khi; no k_rain can be found'
         write(*,11) 'klo, khi', klo, khi
         ierr = 1
         return
      end if

      ! starting from the bottom, skip the core (if there is one)
      inside_core = is_inside_core(khi)

      do k = khi, klo, -1
         if (is_miscible(k) .and. inside_core .and. .not. is_saturated(k)) then
            inside_core = .false.  ! now outside the core
         else if (.not. is_miscible(k) .and. .not. inside_core) then
            k_rain = k
            if (debug) write(*,11) 'find_k_rain: klo, khi, k_rain', klo, khi, k_rain
            return
         end if
      end do

      ! miscible everywhere outside the core
      k_rain = klo - 1
      if (debug) write(*,11) 'find_k_rain: miscible everywhere outside the core; klo, khi, k_rain', klo, khi, k_rain

   contains

      logical function is_inside_core(k)
         ! a cell is inside the core if it is immiscible (a core is about to
         ! form) or if it is miscible but above the miscibility maximum (a
         ! core already exists)
         integer, intent(in) :: k
         is_inside_core = (.not. is_miscible(k)) .or. is_saturated(k)
      end function is_inside_core

      logical function is_saturated(k)
         integer, intent(in) :: k
         is_saturated = (is_miscible(k) .and. Xi_miscible_max_monotonic(k) >= 0d0 &
            .and. Xi_miscible_max_monotonic(k) < s% xa(i_Xi, k))
      end function is_saturated

   end subroutine find_k_rain


   subroutine adjust_miscibility_gap(s, klo, khi, k_start, k_end, is_miscible, Xi_miscible_min, Xi_miscible_max, &
         Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
      ! monotonic miscibility limits over klo:khi (monotonicity enforced in k_start:k_end)
      type(star_info), pointer :: s
      integer, intent(in) :: klo, khi        ! region in which we transfer the element
      integer, intent(in) :: k_start, k_end  ! subregion in which monotonicity is enforced
      logical, allocatable, intent(in) :: is_miscible(:)  ! (1:nz)
      real(dp), allocatable, intent(in) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! (1:nz)
      real(dp), allocatable, intent(out) :: Xi_miscible_min_monotonic(:)  ! (klo:khi) monotonically increasing inwards
      real(dp), allocatable, intent(out) :: Xi_miscible_max_monotonic(:)  ! (klo:khi) monotonically increasing inwards
      integer, intent(out) :: ierr

      integer :: k

      ierr = 0

      allocate(Xi_miscible_min_monotonic(klo:khi), source=Xi_miscible_min(klo:khi))
      allocate(Xi_miscible_max_monotonic(klo:khi), source=Xi_miscible_max(klo:khi))

      call enforce_monotonicity(Xi_miscible_min_monotonic, k_start, k_end, enforce_monotonicity_bottom_up, ierr)
      if (ierr /= 0) return
      call debugging_abundance_array_check( &
         'adjust_miscibility_gap: error in check_abundance_array after enforce_monotonicity (min)', &
         Xi_miscible_min_monotonic(k_start:k_end), k_start, k_end, .true., monotonic_increase, ierr, &
         Xi_miscible_min_monotonic(k_start:k_end) /= -1d0)

      call enforce_monotonicity(Xi_miscible_max_monotonic, k_start, k_end, enforce_monotonicity_top_down, ierr)
      if (ierr /= 0) return
      call debugging_abundance_array_check( &
         'adjust_miscibility_gap: error in check_abundance_array after enforce_monotonicity (max)', &
         Xi_miscible_max_monotonic(k_start:k_end), k_start, k_end, .true., monotonic_increase, ierr, &
         Xi_miscible_max_monotonic(k_start:k_end) /= -1d0)

      ! never increase the abundance in the region above the precipitating region
      ! (we start, typically, from the RCB and stop at the first immiscible cell)
      do k = k_start, k_end
         if (.not. is_miscible(k)) then
            exit
         else if (Xi_miscible_min_monotonic(k) > 0d0 .and. Xi_miscible_min_monotonic(k) > s% xa(i_Xi, k)) then
            if (debug) write(*,'(a,i0,a,1pd26.16,a,1pd26.16)') &
               'adjust_miscibility_gap: reducing Xi_miscible_min_monotonic at k = ', &
               k, ' from ', Xi_miscible_min_monotonic(k), ' to ', s% xa(i_Xi, k)
            Xi_miscible_min_monotonic(k) = s% xa(i_Xi, k)
         end if
      end do

      if (damp_sedimentation) then
         call reduce_Xi_difference(s, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, klo, khi, ierr)
         if (ierr /= 0) return
         call debugging_abundance_array_check( &
            'adjust_miscibility_gap: error in check_abundance_array after damping Xi (min)', &
            Xi_miscible_min_monotonic(k_start:k_end), k_start, k_end, .true., monotonic_increase, ierr, &
            Xi_miscible_min_monotonic(k_start:k_end) /= -1d0)
         call debugging_abundance_array_check( &
            'adjust_miscibility_gap: error in check_abundance_array after damping Xi (max)', &
            Xi_miscible_max_monotonic(k_start:k_end), k_start, k_end, .true., monotonic_increase, ierr, &
            Xi_miscible_max_monotonic(k_start:k_end) /= -1d0)
      end if

   end subroutine adjust_miscibility_gap


   subroutine enforce_monotonicity(Xi, klo, khi, method, ierr)
      ! makes the miscibility limit Xi monotonically increasing inwards in klo:khi.
      ! Entries of -1d0 (no gap) are filled from the neighbour in the direction of the sweep.
      real(dp), allocatable, intent(inout) :: Xi(:)
      integer, intent(in) :: klo, khi  ! cells to be treated
      integer, intent(in) :: method    ! enforce_monotonicity_bottom_up (for the minimum) or _top_down (for the maximum)
      integer, intent(out) :: ierr

      integer :: k

      ierr = 0

      if (debug) write(*,'(a,4i8)') 'enforce_monotonicity: size(Xi), klo, khi, method: ', size(Xi), klo, khi, method

      if (method == enforce_monotonicity_bottom_up) then
         do k = khi, klo+1, -1
            if (Xi(k) >= 0d0 .and. Xi(k-1) < 0d0) then
               Xi(k-1) = Xi(k)
            else if (Xi(k) < 0d0) then
               ! out of phase diagram range, nothing to do
            else if (Xi(k) < Xi(k-1)) then
               if (debug) write(*,'(a,i0,2(1pd26.16))') 'enforce_monotonicity: at k, Xi(k-1), Xi(k) = ', k, Xi(k-1), Xi(k)
               Xi(k-1) = Xi(k)
            end if
         end do
      else if (method == enforce_monotonicity_top_down) then
         do k = klo+1, khi
            if (Xi(k) < 0d0 .and. Xi(k-1) >= 0d0) then
               Xi(k) = Xi(k-1)
            else if (Xi(k) < 0d0) then
               ! out of phase diagram range, nothing to do
            else if (Xi(k) < Xi(k-1)) then
               if (debug) write(*,'(a,i0,2(1pd26.16))') 'enforce_monotonicity: at k, Xi(k), Xi(k-1) = ', k, Xi(k), Xi(k-1)
               Xi(k) = Xi(k-1)
            end if
         end do
      else
         write(*,'(a,i0)') 'enforce_monotonicity: unknown method: ', method
         write(*,'(a,2i8)') 'klo, khi: ', klo, khi
         ierr = 1
         return
      end if

      call check_abundance_array(Xi(klo:khi), klo, khi, .false., monotonic_increase, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'enforce_monotonicity: error in check_abundance_array'
         return
      end if

   end subroutine enforce_monotonicity


   subroutine reduce_Xi_difference(s, Xi_miscible_min, Xi_miscible_max, klo, khi, ierr)
      ! limits the change of the mass fraction per cell to
      ! s% element_sedimentation_maximum_dXi_per_cell relative to the start of the step
      type(star_info), pointer :: s
      real(dp), allocatable, intent(inout) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! (working range)
      integer, intent(in) :: klo, khi
      integer, intent(out) :: ierr

      real(dp), allocatable :: Xi_old(:)
      real(dp) :: dXi_max

      ierr = 0

      allocate(Xi_old(klo:khi))
      Xi_old(klo:khi) = s% xa_old(i_Xi, klo:khi)

      dXi_max = s% element_sedimentation_maximum_dXi_per_cell  ! > 0

      ! is_miscible already guarantees Xi_miscible_min < Xi < Xi_miscible_max, hence no absolute values
      Xi_miscible_min(klo:khi) = merge(Xi_old(klo:khi) - dXi_max, Xi_miscible_min(klo:khi), &
         (Xi_miscible_min(klo:khi) >= 0d0) .and. (Xi_miscible_min(klo:khi) < Xi_old(klo:khi) - dXi_max))
      Xi_miscible_max(klo:khi) = merge(Xi_old(klo:khi) + dXi_max, Xi_miscible_max(klo:khi), &
         (Xi_miscible_max(klo:khi) >= 0d0) .and. (Xi_miscible_max(klo:khi) > Xi_old(klo:khi) + dXi_max))

   end subroutine reduce_Xi_difference


   subroutine transfer_Xi(s, Xi, Xi_ref, k, klo, khi, ierr)
      ! sets Xi(k) = Xi_ref(k) and spreads the mass difference uniformly over the cells klo:khi
      type(star_info), pointer :: s
      real(dp), allocatable, intent(inout) :: Xi(:)   ! (working range) mass fractions
      real(dp), allocatable, intent(in) :: Xi_ref(:)  ! (working range) target mass fractions
      integer, intent(in) :: k                        ! cell that is set to Xi_ref(k)
      integer, intent(in) :: klo, khi                 ! cells that receive the difference
      integer, intent(out) :: ierr

      real(dp) :: dXi    ! change of the mass fraction in cell k
      real(dp) :: dM_Xi  ! transferred mass of the species
      real(dp) :: M      ! total mass of the cells klo:khi
      include 'formats'

      ierr = 0

      if (klo > khi) then
         write(*,'(a)') 'transfer_Xi: invalid bounds, klo > khi'
         write(*,11) 'k, klo, khi', k, klo, khi
         ierr = 1
         return
      end if

      dXi = Xi_ref(k) - Xi(k)
      Xi(k) = Xi_ref(k)

      dM_Xi = s% dm(k)*dXi
      M = sum(s% dm(klo:khi))

      if (M <= 0d0) then
         write(*,'(a)') 'transfer_Xi: total mass of the cells klo:khi is zero or negative'
         write(*,11) 'klo, khi', klo, khi
         write(*,1) 'M', M
         ierr = 1
         return
      end if

      Xi(klo:khi) = Xi(klo:khi) - dM_Xi/M

      if (any(Xi(klo:khi) < 0d0) .or. any(Xi(klo:khi) > 1d0)) then
         write(*,'(a)') 'transfer_Xi: mass fraction out of [0,1] after transfer'
         write(*,11) 'k, klo, khi', k, klo, khi
         write(*,1) 'dXi', dXi
         write(*,1) 'min Xi(klo:khi)', minval(Xi(klo:khi))
         write(*,1) 'max Xi(klo:khi)', maxval(Xi(klo:khi))
         ierr = 1
         return
      end if

      if (debug) write(*,'(a)') 'transfer_Xi: finished.'

   end subroutine transfer_Xi


   subroutine homogenize_region(s, Xi, klo, khi, ierr)
      ! sets Xi in the cells klo:khi to the mass-weighted mean
      type(star_info), pointer :: s
      real(dp), allocatable, intent(inout) :: Xi(:)  ! (working range) mass fractions
      integer, intent(in) :: klo, khi
      integer, intent(out) :: ierr

      real(dp) :: M, M_Xi  ! total mass and mass of the species in klo:khi
      include 'formats'

      ierr = 0

      if (klo > khi) then
         write(*,'(a)') 'homogenize_region: invalid bounds, klo > khi'
         write(*,11) 'klo, khi', klo, khi
         ierr = 1
         return
      end if

      M = sum(s% dm(klo:khi))
      M_Xi = dot_product(s% dm(klo:khi), Xi(klo:khi))

      if (M <= 0d0) then
         write(*,'(a)') 'homogenize_region: total mass of the cells klo:khi is zero or negative'
         write(*,11) 'klo, khi', klo, khi
         write(*,1) 'M', M
         ierr = 1
         return
      end if

      if (M_Xi < 0d0) then
         write(*,'(a)') 'homogenize_region: mass of the species in the cells klo:khi is negative'
         write(*,11) 'klo, khi', klo, khi
         write(*,1) 'M_Xi', M_Xi
         write(*,1) 'M', M
         write(*,1) 'mean Xi', M_Xi/M
         ierr = 1
         return
      end if

      Xi(klo:khi) = M_Xi/M

   end subroutine homogenize_region


   ! ----------------------------------------------------------------------
   ! application to the model
   ! ----------------------------------------------------------------------

   subroutine adjust_mass_fractions(s, dXi, ierr)
      ! applies the change dXi of species i_Xi to s% xa. The other species
      ! are displaced by the counter-flow: the net mass flux of species i_Xi
      ! across every face is balanced by an upwind flux of the remaining
      ! mixture, applied in explicit sub-steps that respect the CFL condition.
      type(star_info), pointer :: s
      real(dp), intent(in) :: dXi(:)  ! (1:nz)
      integer, intent(out) :: ierr

      integer :: k, n_species, step, max_steps, steps_needed
      real(dp), allocatable :: dX_array(:,:)   ! (species, nz)
      real(dp), allocatable :: F_He_total(:)   ! (0:nz) net mass flux of species i_Xi across the face between
                                               !        cells k and k+1, > 0 outwards
      real(dp), allocatable :: F_He_sub(:)     ! (0:nz) per sub-step
      real(dp), allocatable :: F_j(:,:)        ! (species, 0:nz) counter-flow mass fluxes
      real(dp), allocatable :: dXi_sub(:)      ! (nz) per sub-step
      real(dp) :: upwind_factor, donor_mass
      include 'formats'

      ierr = 0

      n_species = size(s% xa, dim=1)

      if (size(dXi) /= s% nz) then
         write(*,'(a)') 'adjust_mass_fractions: dXi does not have nz entries'
         write(*,11) 'size(dXi), nz', size(dXi), s% nz
         ierr = 1
         return
      end if

      if (.not. is_mass_conserved(dXi(1:s% nz), s% dm(1:s% nz))) then
         write(*,'(a)') 'adjust_mass_fractions: mass is not conserved in dXi'
         ierr = 1
         return
      end if

      allocate(dX_array(1:n_species, 1:s% nz), source=0d0)
      allocate(F_He_total(0:s% nz), source=0d0)
      allocate(F_He_sub(0:s% nz), source=0d0)
      allocate(F_j(1:n_species, 0:s% nz), source=0d0)
      allocate(dXi_sub(1:s% nz), source=0d0)

      ! 1. net mass flux of species i_Xi across the faces
      F_He_total(0) = 0d0
      do k = 1, s% nz
         F_He_total(k) = F_He_total(k-1) - dXi(k)*s% dm(k)
      end do
      F_He_total(s% nz) = 0d0

      ! 2. number of sub-steps needed for the counter-flow to satisfy the CFL condition (< 0.5)
      max_steps = 1
      do k = 1, s% nz-1
         if (F_He_total(k) > 0d0) then
            donor_mass = s% dm(k+1)*(1d0 - s% xa(i_Xi, k+1))
         else if (F_He_total(k) < 0d0) then
            donor_mass = s% dm(k)*(1d0 - s% xa(i_Xi, k))
         else
            donor_mass = 1d0  ! arbitrary, the flux is zero
         end if
         if (abs(F_He_total(k)) > 0d0) then
            steps_needed = ceiling(abs(F_He_total(k))/(0.5d0*donor_mass))
            if (steps_needed > max_steps) max_steps = steps_needed
         end if
      end do

      ! 3. explicit sub-stepped advection
      dXi_sub = dXi/dble(max_steps)
      F_He_sub = F_He_total/dble(max_steps)

      do step = 1, max_steps
         ! counter-flow fluxes from the instantaneous abundances of the donor cell
         do k = 1, s% nz-1
            if (F_He_sub(k) > 0d0) then
               upwind_factor = 1d0/(1d0 - s% xa(i_Xi, k+1))
               F_j(:, k) = -F_He_sub(k)*s% xa(:, k+1)*upwind_factor
            else
               upwind_factor = 1d0/(1d0 - s% xa(i_Xi, k))
               F_j(:, k) = -F_He_sub(k)*s% xa(:, k)*upwind_factor
            end if
         end do
         F_j(:, 0) = 0d0
         F_j(:, s% nz) = 0d0

         ! flux divergence
         do k = 1, s% nz
            dX_array(1:n_species, k) = (F_j(1:n_species, k-1) - F_j(1:n_species, k))/s% dm(k)
         end do

         ! the sedimenting species follows the prescribed change
         dX_array(i_Xi, 1:s% nz) = dXi_sub(1:s% nz)

         s% xa(1:n_species, 1:s% nz) = s% xa(1:n_species, 1:s% nz) + dX_array(1:n_species, 1:s% nz)
      end do

      if (debug) write(*,11) 'adjust_mass_fractions: sub-steps', max_steps

   end subroutine adjust_mass_fractions


   subroutine adjust_mixing_in_sedimenting_region(s, is_miscible, ierr)
      ! switches off mixing in the immiscible cells
      type(star_info), pointer :: s
      logical, intent(in) :: is_miscible(:)  ! (1:nz)
      integer, intent(out) :: ierr

      ierr = 0

      if (all(is_miscible(1:s% nz))) then
         if (debug) write(*,'(a)') 'adjust_mixing_in_sedimenting_region: all cells are miscible, no adjustment needed'
         return
      end if

      if (s% do_conv_premix) then
         write(*,'(a)') 'You are using convective premixing and instant helium rain simultaneously. This will likely cause crashes.'
         write(*,'(a)') 'In CPM, regions that rain out will flip-flop between homogeneous and rainout mixing.'
         write(*,'(a)') 'This check sets do_conv_premix = .false.'
         write(*,'(a)') 'If you know what you are doing, feel free to remove this check.'
         s% do_conv_premix = .false.
      end if

      ! simplest approach: cells that are immiscible are not mixed
      s% D_mix(1:s% nz) = merge(s% D_mix(1:s% nz), 0d0, is_miscible(1:s% nz))

   end subroutine adjust_mixing_in_sedimenting_region


   ! ----------------------------------------------------------------------
   ! checks and diagnostics
   ! ----------------------------------------------------------------------

   subroutine check_after_rain(s, klo, khi, xa_pre, ierr)
      ! checks that every species conserved its mass in klo:khi and that the
      ! changed mass fractions are within [0,1]. On failure the messages say
      ! which species and cells are affected and ierr /= 0 is returned, so that
      ! MESA retries the step (the original implementation stopped the run).
      type(star_info), pointer :: s
      integer, intent(in) :: klo, khi
      real(dp), allocatable, intent(in) :: xa_pre(:,:)  ! composition before the rain
      integer, intent(out) :: ierr

      integer :: i, op_err
      real(dp) :: M_chem, M_chem_pre
      include 'formats'

      ierr = 0

      do i = 1, s% species

         M_chem = dot_product(s% xa(i, klo:khi), s% dm(klo:khi))
         M_chem_pre = dot_product(xa_pre(i, klo:khi), s% dm(klo:khi))
         if (abs(M_chem - M_chem_pre) > small*max(abs(M_chem_pre), tiny(1d0))) then
            ierr = 1
            write(*,'(a,i0,a)') 'check_after_rain: mass of species ', i, ' not conserved'
            write(*,3) 'klo, khi, M_chem/M_chem_pre, M_chem_pre, M_chem', klo, khi, M_chem/M_chem_pre, M_chem_pre, M_chem
         end if

         ! no monotonicity check: composition gradients are allowed
         call check_abundance_array(s% xa(i, klo:khi), klo, khi, .true., no_monotonicity, op_err, &
            mask=(xa_pre(i, klo:khi) /= s% xa(i, klo:khi)), Xi_ref=xa_pre(i, klo:khi))
         if (op_err /= 0) then
            write(*,'(a,i0)') 'check_after_rain: abundance check failed for species ', i
            ierr = op_err
         end if

      end do

      call print_maximum_abundance_change(s% xa(i_Xi, klo:khi), xa_pre(i_Xi, klo:khi), klo, do_check=verbose)

      if (ierr /= 0) write(*,'(a)') 'check_after_rain: legacy rain step failed; MESA will retry the step'

   end subroutine check_after_rain


   subroutine print_maximum_abundance_change(Xi, Xi_ref, klo, do_check)
      ! reports the largest change of the mass fraction; Xi and Xi_ref hold the cells klo:khi
      real(dp), intent(in) :: Xi(:)
      real(dp), intent(in) :: Xi_ref(:)
      integer, intent(in) :: klo
      logical, intent(in) :: do_check

      integer :: k_max
      real(dp) :: diff_max
      include 'formats'

      if (.not. do_check) return

      k_max = maxloc(abs(Xi - Xi_ref), dim=1)
      diff_max = abs(Xi(k_max) - Xi_ref(k_max))

      write(*,2) 'print_maximum_abundance_change: k, maximum difference', klo + k_max - 1, diff_max

   end subroutine print_maximum_abundance_change


   subroutine check_abundance_array(Xi, klo, khi, check_range, monotonicity_method, ierr, mask, Xi_ref)
      ! Xi, mask and Xi_ref hold the cells klo:khi (element 1 is cell klo).
      ! check_range: all masked values must be in [0,1].
      ! monotonicity_method: monotonic_increase/monotonic_decrease require the masked
      ! values to be monotonic inwards within the relative tolerance 'small', unless
      ! Xi_ref shows the same violation; no_monotonicity skips the check.
      real(dp), intent(in) :: Xi(:)
      integer, intent(in) :: klo, khi
      logical, intent(in) :: check_range
      integer, intent(in) :: monotonicity_method
      integer, intent(out) :: ierr
      logical, intent(in), optional :: mask(:)
      real(dp), intent(in), optional :: Xi_ref(:)

      integer :: j, n
      logical, allocatable :: local_mask(:), bad(:)
      include 'formats'

      ierr = 0

      n = khi - klo + 1
      if (n <= 0) return

      if (size(Xi) < n) then
         write(*,'(a)') 'check_abundance_array: Xi has fewer entries than cells klo:khi'
         write(*,11) 'size(Xi), klo, khi', size(Xi), klo, khi
         ierr = 1
         return
      end if

      allocate(local_mask(n), bad(n))
      if (present(mask)) then
         if (size(mask) < n) then
            write(*,'(a)') 'check_abundance_array: mask has fewer entries than cells klo:khi'
            write(*,11) 'size(mask), klo, khi', size(mask), klo, khi
            ierr = 1
            return
         end if
         local_mask(1:n) = mask(1:n)
      else
         local_mask(1:n) = .true.
      end if
      if (present(Xi_ref)) then
         if (size(Xi_ref) < n) then
            write(*,'(a)') 'check_abundance_array: Xi_ref has fewer entries than cells klo:khi'
            write(*,11) 'size(Xi_ref), klo, khi', size(Xi_ref), klo, khi
            ierr = 1
            return
         end if
      end if

      if (check_range) then
         bad(1:n) = (Xi(1:n) < 0d0 .or. Xi(1:n) > 1d0) .and. local_mask(1:n)
         if (any(bad(1:n))) then
            ierr = 1
            write(*,'(a)') 'check_abundance_array: Xi contains values outside of [0,1]'
            write(*,'(a)') '   k, Xi(k)'
            do j = 1, n
               if (bad(j)) write(*,'(i12,1pd26.16)') klo + j - 1, Xi(j)
            end do
            return
         end if
      end if

      select case (monotonicity_method)
      case (monotonic_increase, monotonic_decrease)
         bad(1:n) = .false.
         do j = 1, n-1
            if (.not. (local_mask(j) .and. local_mask(j+1))) cycle
            if (violates(Xi(j), Xi(j+1))) then
               if (present(Xi_ref)) then
                  ! skip if the reference array shows the same violation
                  if (violates(Xi_ref(j), Xi_ref(j+1))) cycle
               end if
               bad(j) = .true.
            end if
         end do
         if (any(bad(1:n))) then
            ierr = 1
            if (monotonicity_method == monotonic_increase) then
               write(*,'(a)') 'check_abundance_array: Xi is not monotonically increasing'
            else
               write(*,'(a)') 'check_abundance_array: Xi is not monotonically decreasing'
            end if
            write(*,'(a)') '   k, Xi(k), Xi(k+1)'
            do j = 1, n-1
               if (bad(j)) write(*,'(i12,2(1pd26.16))') klo + j - 1, Xi(j), Xi(j+1)
            end do
            return
         end if
      case (no_monotonicity)
         ! nothing to check
      case default
         write(*,'(a,i0)') 'check_abundance_array: unknown monotonicity method: ', monotonicity_method
         ierr = 1
         return
      end select

   contains

      logical function violates(a, b)
         ! relative decrease (increase) from a to b beyond the tolerance 'small'
         real(dp), intent(in) :: a, b
         if (monotonicity_method == monotonic_increase) then
            violates = (a - b)/b > small
         else
            violates = (a - b)/b < -small
         end if
      end function violates

   end subroutine check_abundance_array


   subroutine debugging_abundance_array_check(message, Xi, klo, khi, check_range, monotonicity_method, ierr, mask)
      ! check_abundance_array that is only active in debug mode and aborts on failure
      character(len=*), intent(in) :: message
      real(dp), intent(in) :: Xi(:)  ! cells klo:khi
      integer, intent(in) :: klo, khi
      logical, intent(in) :: check_range
      integer, intent(in) :: monotonicity_method
      integer, intent(out) :: ierr
      logical, intent(in), optional :: mask(:)  ! cells klo:khi

      ierr = 0
      if (.not. debug) return

      call check_abundance_array(Xi, klo, khi, check_range, monotonicity_method, ierr, mask)
      if (ierr /= 0) then
         write(*,'(a)') trim(message)
         call mesa_error(__FILE__,__LINE__)
      end if

   end subroutine debugging_abundance_array_check


   logical function is_mass_conserved(dXi, dm)
      ! relative change of the total mass of the species below 'small'
      real(dp), intent(in) :: dXi(:), dm(:)
      real(dp) :: dM_species
      dM_species = dot_product(dXi, dm)
      is_mass_conserved = abs(dM_species/sum(dm)) < small
      if (.not. is_mass_conserved) then
         write(*,'(a,2(1pd26.16))') 'is_mass_conserved: mass is not conserved, dM, dM/M = ', dM_species, dM_species/sum(dm)
      else if (debug) then
         write(*,'(a,2(1pd26.16))') 'is_mass_conserved: mass is conserved, dM, dM/M = ', dM_species, dM_species/sum(dm)
      end if
   end function is_mass_conserved


   subroutine print_miscibility_info(s, is_miscible, Xi_miscible_min, Xi_miscible_max, &
         print_only_immiscible_zones, file_path)
      ! writes the miscibility arrays to stdout or to file_path.
      ! The arrays may be allocated over any range of cells; the common range is printed.
      use iso_fortran_env, only: output_unit
      type(star_info), pointer :: s
      logical, allocatable, intent(in) :: is_miscible(:)
      real(dp), allocatable, intent(in) :: Xi_miscible_min(:), Xi_miscible_max(:)
      logical, intent(in) :: print_only_immiscible_zones
      character(len=*), intent(in), optional :: file_path

      integer :: k, klo, khi, out_unit, ios
      logical :: write_to_file

      write_to_file = present(file_path)

      if (write_to_file) then
         open(newunit=out_unit, file=file_path, status='replace', action='write', iostat=ios)
         if (ios /= 0) then
            write(*,'(a)') 'print_miscibility_info: cannot open ' // trim(file_path)
            return
         end if
      else
         out_unit = output_unit
      end if

      klo = max(lbound(is_miscible, dim=1), lbound(Xi_miscible_min, dim=1), lbound(Xi_miscible_max, dim=1), 1)
      khi = min(ubound(is_miscible, dim=1), ubound(Xi_miscible_min, dim=1), ubound(Xi_miscible_max, dim=1), s% nz)

      write(out_unit,'(a)') 'Miscibility information:'
      write(out_unit,'(a)') 'zone index, mass, X_i, miscibility, Xi_miscible_min, Xi_miscible_max, Xi_miscible_min - Xi'
      do k = klo, khi
         if (.not. print_only_immiscible_zones .or. .not. is_miscible(k)) then
            write(out_unit,'(i12,2(1pd26.16),a12,3(1pd26.16))') &
               k, s% m(k), s% xa(i_Xi, k), &
               merge('  miscible', 'immiscible', is_miscible(k)), &
               Xi_miscible_min(k), Xi_miscible_max(k), &
               merge(Xi_miscible_min(k) - s% xa(i_Xi, k), 0d0, .not. is_miscible(k))
         end if
      end do

      if (write_to_file) close(out_unit)

   end subroutine print_miscibility_info

end module element_sedimentation_legacy
