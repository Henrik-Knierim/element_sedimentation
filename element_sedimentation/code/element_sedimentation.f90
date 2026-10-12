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

module element_sedimentation

   ! MESA-facing driver for element sedimentation ("helium rain").
   !
   ! A sedimenting species (element_sedimentation_chemical_species, default he4)
   ! is transported through the model by a conservative advection-diffusion
   ! equation in the Lagrangian mass coordinate: convective mixing acts as
   ! diffusion (MESA's sig built from D_mix) and rain as upwind advection with
   ! a carrier mass flow rate w = 4 pi r^2 rho v. Only the part of the species
   ! above the He-poor miscibility limit of a cell rains out of it, and inflow
   ! into a cell is throttled as the cell approaches its He-rich limit, so a
   ! species-rich layer grows upward from the centre. The other species are
   ! moved by the displaced counter-flow so that every cell keeps sum(xa) = 1
   ! and every species conserves its total mass.
   !
   ! The numerics live in advection_diffusion (no MESA dependency), the phase
   ! diagram in sedimentation_phase_diagram, and the composition arithmetic in
   ! sedimentation_miscibility. The old instant_rain and bottom_up_rain
   ! algorithms are kept in element_sedimentation_legacy.
   !
   ! Called from evolve.f90 after the mixing information of the step has been
   ! set (okay_to_set_mixing_info = .false.); the evolve.f90 snippet installed
   ! by setup.py records energy_start before the call and computes
   ! eps_element_sedimentation from the energy change afterwards, exactly like
   ! MESA's element_diffusion.

   use star_private_def
   use const_def, only: dp, pi4
   use utils_lib, only: is_bad
   use chem_def, only: chem_isos
   use chem_lib, only: chem_get_iso_id
   use advection_diffusion
   use sedimentation_phase_diagram, only: phase_diagram, load_phase_diagram, free_phase_diagram, &
      phase_diagram_error_message
   use sedimentation_miscibility, only: evaluate_miscibility, miscibility_error_message
   use element_sedimentation_legacy, only: init_element_sedimentation_legacy, &
      shutdown_element_sedimentation_legacy, do_instant_rain, do_bottom_up_rain

   implicit none

   private
   public :: do_element_sedimentation, shutdown_element_sedimentation

   character(len=*), parameter :: method_advection_diffusion = 'advection_diffusion_rain'
   character(len=*), parameter :: method_instant = 'instant_rain'
   character(len=*), parameter :: method_bottom_up = 'bottom_up_rain'
   character(len=*), parameter :: velocity_constant = 'constant'

   real(dp), parameter :: sum_renormalize_tol = 1d-10  ! renormalize a cell when |sum xa - 1| exceeds this
   real(dp), parameter :: sum_error_tol = 1d-8         ! error when |sum xa - 1| exceeds this
   real(dp), parameter :: species_mass_rtol = 1d-8     ! error when the total mass of a species changes by more
   real(dp), parameter :: min_one_minus_x = 1d-30      ! floor of (1 - x_i) in the passive counter-flow scale
   real(dp), parameter :: passive_sum_error_tol = 1d-3 ! error when the passive species miss 1 - x_i by more after a step

   ! module state: the phase diagram is loaded once per run (freed by shutdown);
   ! species indices and atomic weights are refreshed on every call (cheap, and
   ! robust against a change of the nuclear net)
   type(phase_diagram), target :: pd
   integer :: i_sed = 0       ! index of the sedimenting species in s% xa
   integer :: i_partner = 0   ! index of the partner species of the two-component phase diagram
   real(dp), allocatable :: A(:)  ! (species) atomic weights
   integer :: legacy_i_sed = 0, legacy_i_partner = 0, legacy_species = 0  ! state the legacy module was initialized with
   integer :: pd_load_count = 0          ! incremented on every (re)load of the phase diagram
   integer :: legacy_pd_load_count = -1  ! pd_load_count when the legacy module copied the phase diagram
   logical :: legacy_notice_printed = .false.
   logical :: conv_premix_notice_printed = .false.
   ! retry damping (element_sedimentation_retry_change_factor): the largest change of the
   ! mass fraction of the sedimenting species in the last attempt and the model it belongs to
   integer :: last_attempt_model_number = -1
   real(dp) :: last_attempt_max_change = 0d0

contains

   ! ---------------------------------------------------------------------------
   ! public entry points
   ! ---------------------------------------------------------------------------

   subroutine do_element_sedimentation(s, dt, ierr)
      ! transport the sedimenting species for one MESA step of length dt [s];
      ! ierr /= 0 (with messages on stdout) makes MESA retry the step
      type (star_info), pointer :: s
      real(dp), intent(in) :: dt
      integer, intent(out) :: ierr

      ierr = 0

      call initialize(s, ierr)
      if (ierr /= 0) then
         write(*,'(a)') ' element_sedimentation: initialization failed; see the messages above'
         return
      end if

      s% eps_element_sedimentation(1:s% nz) = 0d0

      select case (trim(s% element_sedimentation_method))
      case (method_advection_diffusion)
         call do_advection_diffusion_rain(s, dt, ierr)
      case (method_instant)
         call prepare_legacy(s, ierr)
         if (ierr == 0) call do_instant_rain(s, ierr)
      case (method_bottom_up)
         call prepare_legacy(s, ierr)
         if (ierr == 0) call do_bottom_up_rain(s, ierr)
      case default
         ierr = 1
         write(*,'(a)') ' element_sedimentation: unknown element_sedimentation_method ''' // &
            trim(s% element_sedimentation_method) // ''''
         write(*,'(a)') '    valid options are ''' // method_advection_diffusion // ''' (default), ''' // &
            method_instant // ''' and ''' // method_bottom_up // ''' (legacy)'
      end select

      s% need_to_setvars = .true.

   end subroutine do_element_sedimentation


   subroutine shutdown_element_sedimentation()
      ! free the phase diagram and forget the species bookkeeping (called from init.f90)
      call shutdown_element_sedimentation_legacy()
      call free_phase_diagram(pd)
      if (allocated(A)) deallocate(A)
      i_sed = 0
      i_partner = 0
      legacy_i_sed = 0
      legacy_i_partner = 0
      legacy_species = 0
      legacy_pd_load_count = -1
      last_attempt_model_number = -1
      last_attempt_max_change = 0d0
   end subroutine shutdown_element_sedimentation


   ! ---------------------------------------------------------------------------
   ! initialization
   ! ---------------------------------------------------------------------------

   subroutine initialize(s, ierr)
      ! load the phase diagram (once, or again if the phase diagram controls changed)
      ! and resolve the species indices and atomic weights for the current net
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      integer :: j, species
      logical :: reload

      ierr = 0

      if (len_trim(s% element_sedimentation_phase_diagram) == 0) then
         ierr = 1
         write(*,'(a)') ' element_sedimentation: element_sedimentation_phase_diagram is not set'
         return
      end if

      reload = .not. pd% loaded
      if (.not. reload) reload = &
         trim(pd% filename) /= trim(s% element_sedimentation_phase_diagram) .or. &
         pd% T_shift /= s% element_sedimentation_phase_diagram_temperature_shift .or. &
         (pd% extrapolate_in_T .neqv. s% element_sedimentation_extrapolate_phase_diagram)

      if (reload) then
         call load_phase_diagram(pd, trim(s% element_sedimentation_phase_diagram), &
            s% element_sedimentation_phase_diagram_temperature_shift, &
            s% element_sedimentation_extrapolate_phase_diagram, ierr)
         if (ierr /= 0) then
            write(*,'(a)') ' element_sedimentation: failed to load the phase diagram ''' // &
               trim(s% element_sedimentation_phase_diagram) // ''': ' // trim(phase_diagram_error_message(ierr))
            return
         end if
         pd_load_count = pd_load_count + 1
         write(*,'(a,i0,a,i0,a)') ' element_sedimentation: loaded phase diagram ''' // &
            trim(s% element_sedimentation_phase_diagram) // ''' (', pd% n_points, ' points, ', &
            pd% num_unique_P, ' pressure blocks)'
         write(*,'(a,1pd12.4,a,l1)') '    temperature shift [K] = ', pd% T_shift, &
            ', extrapolate in T = ', pd% extrapolate_in_T
      end if

      call resolve_species(s, s% element_sedimentation_chemical_species, &
         'element_sedimentation_chemical_species', i_sed, ierr)
      if (ierr /= 0) return
      call resolve_species(s, s% element_sedimentation_partner_species, &
         'element_sedimentation_partner_species', i_partner, ierr)
      if (ierr /= 0) return
      if (i_sed == i_partner) then
         ierr = 1
         write(*,'(a)') ' element_sedimentation: element_sedimentation_chemical_species and ' // &
            'element_sedimentation_partner_species must be different isotopes'
         return
      end if

      species = s% species
      if (allocated(A)) then
         if (size(A) /= species) deallocate(A)
      end if
      if (.not. allocated(A)) allocate(A(species))
      do j = 1, species
         A(j) = chem_isos% W(s% chem_id(j))
      end do

   end subroutine initialize


   subroutine resolve_species(s, name, control, i, ierr)
      ! index in s% xa of the isotope called name; ierr /= 0 with a message naming the control if it is unknown
      type (star_info), pointer :: s
      character(len=*), intent(in) :: name, control
      integer, intent(out) :: i, ierr
      integer :: cid
      ierr = 0
      i = 0
      cid = chem_get_iso_id(trim(name))
      if (cid <= 0) then
         ierr = 1
         write(*,'(a)') ' element_sedimentation: ' // trim(control) // ' = ''' // trim(name) // &
            ''' is not a known isotope'
         return
      end if
      i = s% net_iso(cid)
      if (i <= 0) then
         ierr = 1
         i = 0
         write(*,'(a)') ' element_sedimentation: ' // trim(control) // ' = ''' // trim(name) // &
            ''' is not in the current nuclear net'
      end if
   end subroutine resolve_species


   subroutine prepare_legacy(s, ierr)
      ! one-time notice and (re)initialization of the legacy module
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      if (.not. legacy_notice_printed) then
         write(*,'(a)') ' element_sedimentation: NOTE: element_sedimentation_method = ''' // &
            trim(s% element_sedimentation_method) // ''' is a legacy method kept for reference;'
         write(*,'(a)') '    ''' // method_advection_diffusion // ''' is the recommended method.'
         legacy_notice_printed = .true.
      end if
      ! the legacy module keeps its own copy of the phase diagram: re-initialize it when the
      ! species or the phase diagram (file, temperature shift, extrapolation) changed
      if (legacy_i_sed /= i_sed .or. legacy_i_partner /= i_partner .or. legacy_species /= s% species &
            .or. legacy_pd_load_count /= pd_load_count) then
         call init_element_sedimentation_legacy(pd, i_sed, i_partner, A, ierr)
         if (ierr /= 0) then
            write(*,'(a)') ' element_sedimentation: init_element_sedimentation_legacy failed'
            return
         end if
         legacy_i_sed = i_sed
         legacy_i_partner = i_partner
         legacy_species = s% species
         legacy_pd_load_count = pd_load_count
      end if
   end subroutine prepare_legacy


   ! ---------------------------------------------------------------------------
   ! advection-diffusion rain
   ! ---------------------------------------------------------------------------

   subroutine do_advection_diffusion_rain(s, dt, ierr)
      ! Transport of the sedimenting species by the advection_diffusion solver
      ! (to steady state, or for the physical time dt), with the passive
      ! counter-flow of all other species applied after every accepted step.
      ! On any failure s% xa is restored to its state at entry and ierr /= 0.
      type (star_info), pointer :: s
      real(dp), intent(in) :: dt
      integer, intent(out) :: ierr

      integer :: nz, species, nvar, k, j, jj, donor, n_gap, k_first, k_last, num_renormalized
      real(dp) :: total_time, v0, rho_face, area, f, xsum, dev, max_sum_dev, mass_end, rel, &
         max_species_dev, mass_moved, cons_err, max_passive_defect, change_limit, &
         max_change, damping
      logical, allocatable :: has_gap(:), in_gap(:)
      real(dp), allocatable :: x_lo(:), x_hi(:), x(:), xa_start(:,:), mass_start(:), xp_old(:,:), xp_new(:,:)
      type(adv_diff_problem) :: p, pp
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats
      include 'formats'

      ierr = 0
      nz = s% nz
      species = s% species
      nvar = species - 1

      ! --- miscibility of every cell at its current (P, T, composition)
      allocate(has_gap(nz), in_gap(nz), x_lo(nz), x_hi(nz))
      call evaluate_miscibility(pd, nz, species, s% xa, A, i_sed, i_partner, s% Peos, s% T, &
         has_gap, in_gap, x_lo, x_hi, ierr)
      if (ierr /= 0) then
         write(*,'(a,i0,a)') ' element_sedimentation: evaluate_miscibility failed with ierr = ', ierr, &
            ' (' // trim(miscibility_error_message(ierr)) // ')'
         return
      end if

      n_gap = 0
      k_first = 0
      k_last = 0
      do k = 1, nz
         if (in_gap(k)) then
            n_gap = n_gap + 1
            if (k_first == 0) k_first = k
            k_last = k
         end if
      end do

      ! nothing can rain unless at least one cell is inside the gap (the threshold
      ! control can only raise this minimum)
      if (n_gap < max(1, s% element_sedimentation_cell_threshold_for_precipitation)) then
         if (s% element_sedimentation_report) write(*,'(a,i0,a,i0,a,i0,a)') &
            ' element_sedimentation: model ', s% model_number, ': ', n_gap, &
            ' cells inside the miscibility gap (threshold ', &
            max(1, s% element_sedimentation_cell_threshold_for_precipitation), '), nothing to do'
         return
      end if

      ! --- sedimentation velocity (cm/s, positive = sinking)
      select case (trim(s% element_sedimentation_velocity_method))
      case (velocity_constant)
         v0 = s% element_sedimentation_velocity
      case default
         ierr = 1
         write(*,'(a)') ' element_sedimentation: unknown element_sedimentation_velocity_method ''' // &
            trim(s% element_sedimentation_velocity_method) // ''''
         write(*,'(a)') '    valid options are ''' // velocity_constant // ''''
         return
      end select
      if (is_bad(v0)) then
         ierr = 1
         write(*,'(a,1pd26.16)') ' element_sedimentation: element_sedimentation_velocity is NaN or Inf: ', v0
         return
      end if

      ! --- the transport problem on MESA's grid (cells 1:nz, faces 1:nz+1, zero boundary flux)
      call alloc_adv_diff_problem(p, nz, ierr)
      if (ierr /= 0) then
         write(*,'(a)') ' element_sedimentation: alloc_adv_diff_problem failed'
         return
      end if
      p% dm(1:nz) = s% dm(1:nz)
      p% cap_width = s% element_sedimentation_cap_width
      do k = 2, nz
         rho_face = (s% dq(k-1)*s% rho(k) + s% dq(k)*s% rho(k-1))/(s% dq(k-1) + s% dq(k))
         area = pi4*s% r(k)*s% r(k)
         f = area*rho_face
         ! MESA's sig (mix_info: get_convection_sigmas) scaled by the user factor
         p% sig(k) = f*f*s% D_mix(k)*s% element_sedimentation_D_mix_factor/(0.5d0*(s% dm(k-1) + s% dm(k)))
         ! carrier mass flow rate; cells without a gap never rain (their floor is
         ! adv_diff_huge, so the flux would be zero anyway), hence such faces are
         ! left inactive and do not enter the timescale bookkeeping of the solver
         if (v0 > 0d0) then
            donor = k - 1
         else
            donor = k
         end if
         if (has_gap(donor)) p% w(k) = f*v0
      end do
      do k = 1, nz
         if (has_gap(k)) then
            p% x_floor(k) = x_lo(k)   ! only the species above the He-poor limit rains out
            p% x_cap(k) = x_hi(k)     ! inflow stops at the He-rich limit
         else
            p% x_floor(k) = adv_diff_huge
            p% x_cap(k) = 1d0
         end if
      end do

      ! largest allowed change of the sedimenting species per cell in this step (damping,
      ! applied after the solve; <= 0: no limit)
      change_limit = change_limit_for_this_attempt(s)

      opts% tol_steady = s% element_sedimentation_steady_state_tolerance
      opts% max_steps = s% element_sedimentation_max_steps
      opts% dt_initial_factor = s% element_sedimentation_advection_diffusion_time_step_factor
      opts% x_upper_bound = 1d0
      opts% report = s% element_sedimentation_report
      if (s% element_sedimentation_solve_to_steady_state) then
         total_time = 0d0
      else
         total_time = dt
      end if

      ! --- state at entry (restored on failure; reference for the conservation checks)
      allocate(x(nz), xa_start(species,nz), mass_start(species), xp_old(nz,nvar), xp_new(nz,nvar))
      xa_start(1:species,1:nz) = s% xa(1:species,1:nz)
      do j = 1, species
         mass_start(j) = dot_product(s% dm(1:nz), s% xa(j,1:nz))
      end do
      x(1:nz) = s% xa(i_sed,1:nz)

      ! --- passive problem of the other species: same grid and diffusion, carrier
      ! flow and donor scale set after every accepted step by after_step; no floor, no cap
      call alloc_adv_diff_problem(pp, nz, ierr)
      if (ierr /= 0) then
         write(*,'(a)') ' element_sedimentation: alloc_adv_diff_problem failed (passive problem)'
         call cleanup()
         return
      end if
      pp% dm(1:nz) = p% dm(1:nz)
      pp% sig(1:nz+1) = p% sig(1:nz+1)
      max_passive_defect = 0d0

      call adv_diff_integrate(p, x, total_time, opts, stats, ierr, after_step)
      max_change = maxval(abs(x(1:nz) - xa_start(i_sed,1:nz)))
      last_attempt_model_number = s% model_number
      last_attempt_max_change = max_change
      if (ierr /= 0) then
         call report_solver_failure(s, p, x, stats, ierr)
         s% xa(1:species,1:nz) = xa_start(1:species,1:nz)
         call cleanup()
         return
      end if

      ! --- damping (element_sedimentation_maximum_dXi_per_cell and the retry damping):
      ! if the steady state changes the sedimenting species in some cell by more than
      ! change_limit, move only part of the way towards it,
      !    xa = xa_start + damping*(xa_steady - xa_start),  damping = change_limit/max_change,
      ! for all species. A linear blend of two states that both conserve every species,
      ! sum to one and are non-negative has the same properties, convective regions that
      ! are uniform in both stay uniform, and the largest change equals change_limit
      ! exactly. The rest of the rain happens in the following steps. (Stopping the
      ! transport integration early, or raising floors, does not have these properties.)
      damping = 1d0
      if (change_limit > 0d0 .and. max_change > change_limit) then
         damping = change_limit/max_change
         s% xa(1:species,1:nz) = xa_start(1:species,1:nz) + &
            damping*(s% xa(1:species,1:nz) - xa_start(1:species,1:nz))
         last_attempt_max_change = change_limit
      end if

      ! --- post checks: bounds, sum of mass fractions, conservation of every species
      do k = 1, nz
         do j = 1, species
            if (is_bad(s% xa(j,k)) .or. s% xa(j,k) < 0d0) then
               ierr = 1
               write(*,'(a)') ' element_sedimentation: negative or NaN mass fraction after sedimentation'
               write(*,3) 'k, j, xa(j,k) for ' // trim(chem_isos% name(s% chem_id(j))), k, j, s% xa(j,k)
               s% xa(1:species,1:nz) = xa_start(1:species,1:nz)
               call cleanup()
               return
            end if
         end do
      end do

      max_species_dev = 0d0
      do j = 1, species
         mass_end = dot_product(s% dm(1:nz), s% xa(j,1:nz))
         if (mass_start(j) > 0d0) then
            rel = abs(mass_end - mass_start(j))/mass_start(j)
         else if (mass_end /= 0d0) then
            rel = 1d0  ! a species with zero mass must stay exactly zero
         else
            rel = 0d0
         end if
         max_species_dev = max(max_species_dev, rel)
         if (rel > species_mass_rtol) then
            ierr = 1
            write(*,'(a)') ' element_sedimentation: total mass of a species changed beyond tolerance'
            write(*,2) 'species index, mass before, mass after [g] for ' // trim(chem_isos% name(s% chem_id(j))), &
               j, mass_start(j), mass_end
            write(*,1) 'relative change, tolerance', rel, species_mass_rtol
            s% xa(1:species,1:nz) = xa_start(1:species,1:nz)
            call cleanup()
            return
         end if
      end do

      max_sum_dev = 0d0
      num_renormalized = 0
      do k = 1, nz
         xsum = sum(s% xa(1:species,k))
         dev = abs(xsum - 1d0)
         max_sum_dev = max(max_sum_dev, dev)
         if (dev > sum_error_tol) then
            ierr = 1
            write(*,'(a)') ' element_sedimentation: mass fractions do not sum to 1 after sedimentation'
            write(*,2) 'k, sum(xa(:,k)) - 1, tolerance', k, xsum - 1d0, sum_error_tol
            s% xa(1:species,1:nz) = xa_start(1:species,1:nz)
            call cleanup()
            return
         else if (dev > sum_renormalize_tol) then
            s% xa(1:species,k) = s% xa(1:species,k)/xsum
            num_renormalized = num_renormalized + 1
         end if
      end do

      ! --- optionally switch off MESA's mixing where the mixture is immiscible
      if (s% element_sedimentation_suppress_mixing_in_immiscible_regions) then
         if (s% do_conv_premix) then
            if (.not. conv_premix_notice_printed) then
               write(*,'(a)') ' element_sedimentation: WARNING: convective premixing and element sedimentation ' // &
                  'with suppressed mixing are both on.'
               write(*,'(a)') '    Regions that rain out would flip-flop between homogeneous and rained-out mixing,'
               write(*,'(a)') '    so do_conv_premix is set to .false. for the rest of the run.'
               conv_premix_notice_printed = .true.
            end if
            s% do_conv_premix = .false.
         end if
         do k = 1, nz
            if (.not. in_gap(k)) cycle
            s% D_mix(k) = 0d0
            s% sig(k) = 0d0
            if (k < nz) then
               s% D_mix(k+1) = 0d0
               s% sig(k+1) = 0d0
            end if
         end do
      end if

      ! --- summary
      if (s% element_sedimentation_report) then
         mass_moved = 0d0
         do k = 1, nz
            mass_moved = mass_moved + s% dm(k)*max(s% xa(i_sed,k) - xa_start(i_sed,k), 0d0)
         end do
         if (mass_start(i_sed) > 0d0) then
            cons_err = abs(dot_product(s% dm(1:nz), s% xa(i_sed,1:nz)) - mass_start(i_sed))/mass_start(i_sed)
         else
            cons_err = 0d0
         end if
         write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,a,a,1pd10.3,a,1pd10.3,a,1pd9.2,a,1pd9.2,a,1pd9.2,a,1pd9.2,a,i0)') &
            ' element_sedimentation: model ', s% model_number, &
            ' n_gap ', n_gap, ' (k ', k_first, '-', k_last, ') steps ', stats% steps, &
            ' rejected ', stats% rejected_steps, &
            ' newton ', stats% newton_iters, ' stop ', trim(stop_reason_name(stats% stop_reason)), &
            ' max_dX ', maxval(abs(s% xa(i_sed,1:nz) - xa_start(i_sed,1:nz))), ' mass_moved[g] ', mass_moved, &
            ' cons_err ', cons_err, ' damping ', damping, ' passive_defect ', max_passive_defect, &
            ' max|sum-1| ', max_sum_dev, ' renormalized ', num_renormalized
         if (stats% num_clipped > 0) write(*,11) 'element_sedimentation: tiny negatives clipped', stats% num_clipped
         if (max_species_dev > 0d0) write(*,'(a,1pd12.4)') &
            ' element_sedimentation: max relative change of a species total mass ', max_species_dev
      end if

      call cleanup()

   contains

      subroutine after_step(dt_step, x_new, flux_adv, flux_diff, change_other, ierr_step)
         ! Passive counter-flow of all other species after an accepted step of the
         ! active one: each passive species is diffused with the same sig and
         ! advected with the carrier flow w' = -flux_adv (the displaced mixture)
         ! and donor scale 1/(1 - x_i) evaluated at the new x_i. The implicit
         ! solve with the same dt keeps sum(xa) = 1 to round-off.
         ! Host association provides s, pp, xp_old, xp_new, nz, nvar, species.
         real(dp), intent(in) :: dt_step
         real(dp), intent(in) :: x_new(:)
         real(dp), intent(in) :: flux_adv(:)
         real(dp), intent(in) :: flux_diff(:)  ! not needed: the passive diffusion uses pp% sig
         real(dp), intent(out) :: change_other ! largest change of a passive species in this step
         integer, intent(out) :: ierr_step
         integer :: kk, j1, jj1
         real(dp) :: psum, x_rest, defect

         ierr_step = 0
         change_other = 0d0

         pp% w(1) = 0d0
         pp% w(nz+1) = 0d0
         do kk = 2, nz
            pp% w(kk) = -flux_adv(kk)
         end do
         do kk = 1, nz
            pp% scale(kk) = 1d0/max(1d0 - x_new(kk), min_one_minus_x)
         end do

         jj1 = 0
         do j1 = 1, species
            if (j1 == i_sed) cycle
            jj1 = jj1 + 1
            xp_old(1:nz,jj1) = s% xa(j1,1:nz)
         end do

         call adv_diff_passive_step(pp, dt_step, nvar, xp_old, xp_new, ierr_step)
         if (ierr_step /= 0) then
            write(*,'(a,i0,a)') ' element_sedimentation: passive counter-flow step failed with ierr = ', &
               ierr_step, ' (' // trim(adv_diff_error_message(ierr_step)) // ')'
            return
         end if

         ! In exact arithmetic the passive species sum to 1 - x_i after the solve. In floating
         ! point they miss it by the residual of the accepted active state, which is limited by
         ! the representability of x: where the active Jacobian is large compared with the
         ! passive coupling (a rain donor just above its floor, a receiver in the cap ramp) the
         ! defect is ~ epsilon * x * dt/tau_cell, i.e. 1d-11 for dt/tau ~ 1d4 and up to ~1d-6
         ! for dt/tau ~ 1d13 (seen in a forced-rain test). Rescale the passive species of every
         ! cell to exactly 1 - x_i, which distributes the defect in proportion to their
         ! abundances (their total masses change only at second order, checked below) and
         ! keeps it from accumulating over the steps; a defect beyond passive_sum_error_tol
         ! means something else went wrong.
         do kk = 1, nz
            psum = sum(xp_new(kk,1:nvar))
            x_rest = max(1d0 - x_new(kk), 0d0)
            defect = abs(psum - x_rest)
            max_passive_defect = max(max_passive_defect, defect)
            if (defect > passive_sum_error_tol .or. is_bad(psum)) then
               ierr_step = 1
               write(*,'(a,i0,a,1pd16.8,a,1pd16.8)') &
                  ' element_sedimentation: passive species do not sum to 1 - x_i in cell ', kk, &
                  ': sum = ', psum, ', 1 - x_i = ', x_rest
               return
            end if
            if (psum > 0d0) xp_new(kk,1:nvar) = xp_new(kk,1:nvar)*(x_rest/psum)
         end do
         ! the passive species enter the steady-state test of the solver, so that it does not
         ! stop while they are still being mixed
         change_other = maxval(abs(xp_new(1:nz,1:nvar) - xp_old(1:nz,1:nvar)))

         jj1 = 0
         do j1 = 1, species
            if (j1 == i_sed) cycle
            jj1 = jj1 + 1
            s% xa(j1,1:nz) = xp_new(1:nz,jj1)
         end do
         s% xa(i_sed,1:nz) = x_new(1:nz)

      end subroutine after_step


      subroutine cleanup()
         call free_adv_diff_problem(p)
         call free_adv_diff_problem(pp)
      end subroutine cleanup

   end subroutine do_advection_diffusion_rain


   real(dp) function change_limit_for_this_attempt(s) result(limit)
      ! Largest allowed change of the mass fraction of the sedimenting species in a cell
      ! for this attempt (applied by damping the step): the user limit
      ! element_sedimentation_maximum_dXi_per_cell, and, when MESA is retrying the same
      ! model after an attempt that transported material, element_sedimentation_retry_change_factor
      ! times the largest change applied in that attempt (a steady-state solve does not
      ! depend on the time step, so without this every retry would reproduce the
      ! composition change that made the step fail). Returns a value <= 0 when there is
      ! no limit.
      type (star_info), pointer :: s
      real(dp) :: retry_limit
      logical :: is_retry

      limit = s% element_sedimentation_maximum_dXi_per_cell

      is_retry = s% model_number_for_last_retry == s% model_number .and. &
         last_attempt_model_number == s% model_number .and. &
         last_attempt_max_change > 0d0 .and. &
         s% element_sedimentation_retry_change_factor > 0d0
      if (.not. is_retry) return

      retry_limit = s% element_sedimentation_retry_change_factor*last_attempt_max_change
      if (limit <= 0d0 .or. retry_limit < limit) then
         limit = retry_limit
         write(*,'(a,i0,a,1pd10.3,a,1pd10.3,a)') &
            ' element_sedimentation: model ', s% model_number, &
            ' is a retry; the change per cell is limited to ', limit, &
            ' (previous attempt: ', last_attempt_max_change, &
            '; see element_sedimentation_retry_change_factor)'
      end if
   end function change_limit_for_this_attempt


   subroutine report_solver_failure(s, p, x, stats, ierr)
      ! diagnostic block printed when adv_diff_integrate fails
      type (star_info), pointer :: s
      type(adv_diff_problem), intent(in) :: p
      real(dp), intent(in) :: x(:)
      type(adv_diff_stats), intent(in) :: stats
      integer, intent(in) :: ierr
      integer :: k
      include 'formats'
      write(*,'(a)') ' element_sedimentation: advection_diffusion_rain FAILED'
      write(*,11) 'model_number', s% model_number
      write(*,'(a55,2x,a)') 'method', trim(s% element_sedimentation_method)
      write(*,'(a55,i12,2x,a)') 'ierr', ierr, '(' // trim(adv_diff_error_message(ierr)) // ')'
      write(*,'(a55,2x,a)') 'message', trim(stats% message)
      write(*,'(a55,2x,a)') 'stop_reason', trim(stop_reason_name(stats% stop_reason))
      write(*,11) 'steps', stats% steps
      write(*,11) 'rejected_steps', stats% rejected_steps
      write(*,11) 'newton_iters', stats% newton_iters
      write(*,1) 'time integrated [s]', stats% time
      write(*,1) 'dt_last [s]', stats% dt_last
      write(*,1) 'tau_min [s]', stats% tau_min
      write(*,1) 'tau_max [s]', stats% tau_max
      write(*,1) 'residual', stats% residual
      write(*,1) 'max_change_last', stats% max_change_last
      write(*,1) 'max_change_total', stats% max_change_total
      write(*,11) 'k_fail', stats% k_fail
      k = stats% k_fail
      if (k >= 1 .and. k <= p% nz) then
         write(*,2) 'x(k)', k, x(k)
         write(*,2) 'x_floor(k)', k, p% x_floor(k)
         write(*,2) 'x_cap(k)', k, p% x_cap(k)
         write(*,2) 'dm(k) [g]', k, p% dm(k)
         write(*,2) 'sig(k), sig(k+1) [g/s]', k, p% sig(k), p% sig(k+1)
         write(*,2) 'w(k), w(k+1) [g/s]', k, p% w(k), p% w(k+1)
         write(*,2) 'Peos(k) [dyn/cm^2], T(k) [K]', k, s% Peos(k), s% T(k)
      end if
      write(*,'(a)') ' element_sedimentation: s% xa restored to its state before the call'
   end subroutine report_solver_failure


   function stop_reason_name(stop_reason) result(name)
      integer, intent(in) :: stop_reason
      character(len=16) :: name
      select case (stop_reason)
      case (adv_diff_stop_steady_state)
         name = 'steady_state'
      case (adv_diff_stop_total_time)
         name = 'total_time'
      case (adv_diff_stop_change_limit)
         name = 'change_limit'
      case (adv_diff_stop_nothing_to_do)
         name = 'nothing_to_do'
      case default
         name = 'none'
      end select
   end function stop_reason_name

end module element_sedimentation
