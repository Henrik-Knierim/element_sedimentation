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
   ! sedimentation_miscibility. Layers that the transport leaves convectively
   ! unstable (a heavy deposit on top of cells that were not mixed at the start
   ! of the step) are mixed at once by sedimentation_overturn with the Ledoux
   ! criterion evaluated here.
   !
   ! Called from evolve.f90 after the mixing information of the step has been
   ! set (okay_to_set_mixing_info = .false.). The evolve.f90 snippet installed
   ! by setup.py records energy_start before the call, re-evaluates the EOS with
   ! the new composition afterwards (set_vars) and then calls
   ! finish_element_sedimentation, which sets eps_element_sedimentation for the
   ! energy equation: the change of the internal energy at fixed (rho, T),
   ! exactly like MESA's element_diffusion, plus the divergence of the enthalpy
   ! flux carried by the moved species (sedimentation_energy).

   use star_private_def
   use const_def, only: dp, pi4, ln10, crad, boltzm, mp
   use utils_lib, only: is_bad
   use chem_def, only: chem_isos
   use chem_lib, only: chem_get_iso_id
   use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_lnE, i_lnPgas
   use eos_support, only: get_eos
   use advection_diffusion
   use sedimentation_phase_diagram, only: phase_diagram, load_phase_diagram, free_phase_diagram, &
      phase_diagram_error_message
   use sedimentation_miscibility, only: evaluate_miscibility, miscibility_error_message
   use sedimentation_energy, only: species_face_flow_rates, face_energy_flux, &
      flux_divergence_heating, partial_enthalpy_increment
   use sedimentation_overturn, only: overturn_unstable_layers

   implicit none

   private
   public :: do_element_sedimentation, finish_element_sedimentation, shutdown_element_sedimentation

   character(len=*), parameter :: method_advection_diffusion = 'advection_diffusion_rain'
   character(len=*), parameter :: method_instant = 'instant_rain'
   character(len=*), parameter :: method_bottom_up = 'bottom_up_rain'
   character(len=*), parameter :: velocity_constant = 'constant'

   real(dp), parameter :: sum_renormalize_tol = 1d-10  ! renormalize a cell when |sum xa - 1| exceeds this
   real(dp), parameter :: sum_error_tol = 1d-8         ! error when |sum xa - 1| exceeds this
   real(dp), parameter :: species_mass_rtol = 1d-8     ! error when the total mass of a species changes by more
   real(dp), parameter :: min_one_minus_x = 1d-30      ! floor of (1 - x_i) in the passive counter-flow scale
   real(dp), parameter :: passive_sum_error_tol = 1d-3 ! error when the passive species miss 1 - x_i by more after a step
   real(dp), parameter :: saturation_tol = 1d-10       ! a cell counts as saturated if x_i >= x_lo - saturation_tol
   real(dp), parameter :: transfer_round_tol = 1d-12   ! bottom_up_rain stops its rounds when a transfer moves less than this
   ! finite differences of the EOS for the partial specific enthalpies (enthalpy_flux_heating)
   real(dp), parameter :: dxa_eos = 1d-3      ! composition step
   real(dp), parameter :: dxa_eos_min = 1d-6  ! below this abundance of both species the derivative is not needed
   ! overturn of unstable layers: two blocks whose compositions agree to this are never tested
   real(dp), parameter :: overturn_composition_tol = 1d-12

   ! module state: the phase diagram is loaded once per run (freed by shutdown);
   ! species indices and atomic weights are refreshed on every call (cheap, and
   ! robust against a change of the nuclear net)
   type(phase_diagram), target :: pd
   integer :: i_sed = 0       ! index of the sedimenting species in s% xa
   integer :: i_partner = 0   ! index of the partner species of the two-component phase diagram
   real(dp), allocatable :: A(:)  ! (species) atomic weights
   logical :: conv_premix_notice_printed = .false.
   ! retry damping (element_sedimentation_retry_change_factor): the largest change of the
   ! mass fraction of the sedimenting species in the last attempt and the model it belongs to
   integer :: last_attempt_model_number = -1
   real(dp) :: last_attempt_max_change = 0d0
   real(dp) :: last_attempt_dt = 0d0
   ! energy carried by the moved species (element_sedimentation_include_enthalpy_flux):
   ! heating rate of every cell from the divergence of the enthalpy flux of the last
   ! call, added to eps_element_sedimentation by finish_element_sedimentation
   real(dp), allocatable :: eps_enthalpy_flux(:)
   logical :: enthalpy_flux_valid = .false.

contains

   ! ---------------------------------------------------------------------------
   ! public entry points
   ! ---------------------------------------------------------------------------

   subroutine do_element_sedimentation(s, dt, ierr)
      ! transport the sedimenting species for one MESA step of length dt [s];
      ! ierr > 0 (with messages on stdout) makes MESA retry the step, ierr < 0 is a
      ! configuration error (unknown method, phase diagram, species) that ends the run
      ! (evolve.f90 snippet)
      type (star_info), pointer :: s
      real(dp), intent(in) :: dt
      integer, intent(out) :: ierr
      real(dp), allocatable :: xa_entry(:,:)
      integer(8) :: clock_start, clock_end, clock_rate
      real(dp) :: seconds, max_dx

      ierr = 0
      enthalpy_flux_valid = .false.
      call system_clock(clock_start, clock_rate)

      call initialize(s, ierr)
      if (ierr /= 0) then
         write(*,'(a)') ' element_sedimentation: initialization failed (phase diagram or species); see the messages above'
         ierr = -1  ! configuration error: the run ends
         return
      end if

      s% eps_element_sedimentation(1:s% nz) = 0d0
      allocate(xa_entry(s% species, s% nz))
      xa_entry(1:s% species,1:s% nz) = s% xa(1:s% species,1:s% nz)

      select case (trim(s% element_sedimentation_method))
      case (method_advection_diffusion, method_instant, method_bottom_up)
         call do_advection_diffusion_rain(s, dt, ierr)
      case default
         ierr = -1  ! configuration error: the run ends
         write(*,'(a)') ' element_sedimentation: unknown element_sedimentation_method ''' // &
            trim(s% element_sedimentation_method) // ''''
         write(*,'(a)') '    valid options are ''' // method_advection_diffusion // ''' (default), ''' // &
            method_instant // ''' and ''' // method_bottom_up // ''''
      end select
      if (ierr /= 0) return

      ! energy carried by the moved species, from the composition change of the step
      ! (any method)
      if (s% element_sedimentation_include_enthalpy_flux) then
         call enthalpy_flux_heating(s, dt, xa_entry, ierr)
         if (ierr /= 0) then
            write(*,'(a)') ' element_sedimentation: enthalpy flux of the moved species failed; see the messages above'
            s% xa(1:s% species,1:s% nz) = xa_entry(1:s% species,1:s% nz)
            return
         end if
      end if

      s% need_to_setvars = .true.

      if (s% element_sedimentation_report) then
         max_dx = maxval(abs(s% xa(i_sed,1:s% nz) - xa_entry(i_sed,1:s% nz)))
         if (max_dx > 0d0) then
            call system_clock(clock_end)
            seconds = dble(clock_end - clock_start)/dble(max(clock_rate, 1_8))
            write(*,'(a,i0,a,f8.3,a)') ' element_sedimentation: model ', s% model_number, &
               ' wall time of the module ', seconds, ' s'
         end if
      end if

   end subroutine do_element_sedimentation


   subroutine finish_element_sedimentation(s, dt)
      ! eps_element_sedimentation for the energy equation; called from evolve.f90 after
      ! set_vars has re-evaluated the EOS with the new composition. Like MESA's
      ! finish_element_diffusion: the change of the specific internal energy at fixed
      ! (rho, T) (energy_start was recorded before the composition change), so that
      ! every cell keeps its internal energy while its composition changes and the
      ! total energy is conserved. On top of that, when
      ! element_sedimentation_include_enthalpy_flux is true, the divergence of the
      ! enthalpy flux of the moved species (see sedimentation_energy), which redistributes
      ! energy between the cells and sums to zero over the star.
      type (star_info), pointer :: s
      real(dp), intent(in) :: dt
      integer :: k
      do k = 1, s% nz
         s% eps_element_sedimentation(k) = (s% energy_start(k) - s% energy(k))/dt
      end do
      if (enthalpy_flux_valid) then
         do k = 1, s% nz
            s% eps_element_sedimentation(k) = s% eps_element_sedimentation(k) + eps_enthalpy_flux(k)
         end do
      end if
   end subroutine finish_element_sedimentation


   subroutine shutdown_element_sedimentation()
      ! free the phase diagram and forget the species bookkeeping (called from init.f90)
      call free_phase_diagram(pd)
      if (allocated(A)) deallocate(A)
      if (allocated(eps_enthalpy_flux)) deallocate(eps_enthalpy_flux)
      enthalpy_flux_valid = .false.
      i_sed = 0
      i_partner = 0
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


   ! ---------------------------------------------------------------------------
   ! advection-diffusion rain
   ! ---------------------------------------------------------------------------

   subroutine do_advection_diffusion_rain(s, dt, ierr)
      ! Transport of the sedimenting species on MESA's grid with the floors (He-poor
      ! limits) and caps (He-rich limits) of the phase diagram, for the three methods
      ! that share this path:
      !   'advection_diffusion_rain'  the advection_diffusion solver (to steady state,
      !                               or for the physical time dt) with the passive
      !                               counter-flow of all other species after every
      !                               accepted step;
      !   'instant_rain'              the same steady state in the limit of infinitely
      !                               fast rain (adv_diff_instant_transfer), one
      !                               transfer from the composition homogenized over
      !                               the mixed regions, then the counter-flow;
      !   'bottom_up_rain'            'instant_rain' repeated with the miscibility
      !                               limits re-evaluated from the new composition
      !                               after every transfer, until nothing moves
      !                               (element_sedimentation_transfer_rounds).
      ! Damping, overturn of unstable layers, checks, mixing takeover and report are
      ! common. On any failure s% xa is restored to its state at entry and ierr /= 0.
      type (star_info), pointer :: s
      real(dp), intent(in) :: dt
      integer, intent(out) :: ierr

      integer :: nz, species, nvar, k, j, jj, donor, n_gap, n_sat, k_first, k_last, num_renormalized, &
         n_merged_faces, n_mixed_blocks, k_overturn_outer, k_overturn_inner, n_contrast, n_rounds, &
         n_dummy1, n_dummy2, k_dummy1, k_dummy2
      logical :: rain_active, instant_method, bottom_up_method
      real(dp), allocatable :: w_face(:), flux_inst(:), flux_zero(:)
      real(dp) :: change_other_inst
      real(dp) :: total_time, v0, rho_face, area, f, xsum, dev, max_sum_dev, mass_end, rel, &
         max_species_dev, mass_moved, cons_err, max_passive_defect, change_limit, retry_limit, &
         max_change, max_change_total, damping
      logical, allocatable :: has_gap(:), in_gap(:)
      real(dp), allocatable :: x_lo(:), x_hi(:), x(:), xa_start(:,:), xa_base(:,:), mass_start(:), &
         xp_old(:,:), xp_new(:,:)
      type(adv_diff_problem) :: p, pp, p_mix
      type(adv_diff_options) :: opts
      type(adv_diff_stats) :: stats, stats_mix
      include 'formats'

      ierr = 0
      nz = s% nz
      species = s% species
      nvar = species - 1
      instant_method = (trim(s% element_sedimentation_method) == method_instant .or. &
         trim(s% element_sedimentation_method) == method_bottom_up)
      bottom_up_method = (trim(s% element_sedimentation_method) == method_bottom_up)
      n_rounds = 0

      ! --- miscibility of every cell at its current (P, T, composition)
      allocate(has_gap(nz), in_gap(nz), x_lo(nz), x_hi(nz))
      call evaluate_miscibility(pd, nz, species, s% xa, A, i_sed, i_partner, s% Peos, s% T, &
         has_gap, in_gap, x_lo, x_hi, ierr)
      if (ierr /= 0) then
         write(*,'(a,i0,a)') ' element_sedimentation: evaluate_miscibility failed with ierr = ', ierr, &
            ' (' // trim(miscibility_error_message(ierr)) // ')'
         return
      end if

      ! cells inside the gap (strictly between the limits) and saturated cells (at or above
      ! the He-poor limit, which includes the cells a previous step left exactly at their
      ! limit). The module acts as long as some cell is saturated: otherwise a rain region
      ! that sits at its limit would be handed back to MESA's mixing in the next step, be
      ! homogenized, and rain again in the step after.
      n_gap = 0
      n_sat = 0
      k_first = 0
      k_last = 0
      do k = 1, nz
         if (in_gap(k)) n_gap = n_gap + 1
         if (has_gap(k)) then
            if (s% xa(i_sed,k) >= x_lo(k) - saturation_tol) then
               n_sat = n_sat + 1
               if (k_first == 0) k_first = k
               k_last = k
            end if
         end if
      end do

      ! nothing can rain unless at least one cell is saturated (the threshold control can
      ! only raise this minimum)
      rain_active = (n_sat >= max(1, s% element_sedimentation_cell_threshold_for_precipitation))

      ! Faces with mixing across which the species differs: MESA's mixing would change
      ! the composition there, and MESA's mixing carries no enthalpy flux (every cell
      ! keeps its internal energy while its composition changes). After rain the
      ! species-poor envelope sits on the species-rich interior, and when MESA's
      ! convection reaches through that interface the exchange is of the same kind as
      ! the rain itself; so the module takes such steps over as well (mixing only, with
      ! its energy terms), see element_sedimentation_act_on_mixing_contrast.
      n_contrast = 0
      if (s% element_sedimentation_act_on_mixing_contrast) then
         do k = 2, nz
            if (s% D_mix(k) <= 0d0) cycle
            if (abs(s% xa(i_sed,k-1) - s% xa(i_sed,k)) > s% element_sedimentation_mixing_contrast_tolerance) &
               n_contrast = n_contrast + 1
         end do
      end if

      if (.not. rain_active .and. n_contrast == 0) then
         if (s% element_sedimentation_report) write(*,'(a,i0,a,i0,a,i0,a)') &
            ' element_sedimentation: model ', s% model_number, ': ', n_sat, &
            ' cells at or above the He-poor limit (threshold ', &
            max(1, s% element_sedimentation_cell_threshold_for_precipitation), &
            ') and no mixed face with a composition contrast, nothing to do'
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
      allocate(w_face(nz+1))
      w_face(1) = 0d0
      w_face(nz+1) = 0d0
      do k = 2, nz
         rho_face = (s% dq(k-1)*s% rho(k) + s% dq(k)*s% rho(k-1))/(s% dq(k-1) + s% dq(k))
         area = pi4*s% r(k)*s% r(k)
         f = area*rho_face
         ! MESA's sig (mix_info: get_convection_sigmas) scaled by the user factor
         p% sig(k) = f*f*s% D_mix(k)*s% element_sedimentation_D_mix_factor/(0.5d0*(s% dm(k-1) + s% dm(k)))
         w_face(k) = f*v0   ! carrier mass flow rate [g/s] if the face rains
      end do
      call set_limits_and_flow()

      ! largest allowed change of the sedimenting species per cell in this step (damping,
      ! applied after the solve; <= 0: no limit): the user limit, relative to the
      ! mixing-only state, and on a MESA retry a limit on the whole change of the step
      call change_limits_for_this_attempt(s, change_limit, retry_limit)

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
      allocate(x(nz), xa_start(species,nz), xa_base(species,nz), mass_start(species), xp_old(nz,nvar), xp_new(nz,nvar))
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
      ! the instant transfer hands after_step the whole change of the active species as
      ! fluxes with dt = 1 (masses); the other species were already homogenized over the
      ! mixed regions, so their passive step must carry the counter-flow only (a sig with
      ! dt = 1 would add one second of diffusion, a unit-dependent term)
      if (instant_method) pp% sig(1:nz+1) = 0d0
      max_passive_defect = 0d0

      ! --- reference state of the damping: the same transport without rain (mixing only).
      ! Only needed when a change limit applies. In mixed regions this is the homogenized
      ! composition, in unmixed (e.g. radiative) cells the composition at the start.
      xa_base(1:species,1:nz) = xa_start(1:species,1:nz)
      if (instant_method) then
         ! the instant methods start from the composition homogenized over the mixed
         ! regions (every species, like the solver's diffusion in the steady state);
         ! that state is also the reference of the damping
         call overturn_unstable_layers(nz, species, s% dm, xa_base, mixed_face, &
            n_dummy1, n_dummy2, k_dummy1, k_dummy2)
      else if (change_limit > 0d0) then
         call alloc_adv_diff_problem(p_mix, nz, ierr)
         if (ierr /= 0) then
            write(*,'(a)') ' element_sedimentation: alloc_adv_diff_problem failed (mixing-only problem)'
            call cleanup()
            return
         end if
         p_mix% dm(1:nz) = p% dm(1:nz)
         p_mix% sig(1:nz+1) = p% sig(1:nz+1)
         call adv_diff_integrate(p_mix, x, total_time, opts, stats_mix, ierr, after_step)
         if (ierr /= 0) then
            write(*,'(a)') ' element_sedimentation: mixing-only reference state for the damping failed'
            call report_solver_failure(s, p_mix, x, stats_mix, ierr)
            s% xa(1:species,1:nz) = xa_start(1:species,1:nz)
            call cleanup()
            return
         end if
         xa_base(1:species,1:nz) = s% xa(1:species,1:nz)
         s% xa(1:species,1:nz) = xa_start(1:species,1:nz)
         x(1:nz) = xa_start(i_sed,1:nz)
      end if

      if (instant_method) then
         s% xa(1:species,1:nz) = xa_base(1:species,1:nz)
         allocate(flux_inst(nz+1), flux_zero(nz+1))
         flux_zero(1:nz+1) = 0d0
         do
            n_rounds = n_rounds + 1
            x(1:nz) = s% xa(i_sed,1:nz)
            call adv_diff_instant_transfer(p, x, flux_inst, s% element_sedimentation_max_steps, stats, ierr)
            if (ierr /= 0) exit
            ! the counter-flow of the other species for the whole transfer (dt = 1: the
            ! fluxes are masses)
            call after_step(1d0, x, flux_inst, flux_zero, change_other_inst, ierr)
            if (ierr /= 0) then
               ierr = adv_diff_err_callback
               stats% message = 'passive counter-flow of the instant transfer failed'
               exit
            end if
            if (.not. bottom_up_method) exit
            if (n_rounds >= max(1, s% element_sedimentation_transfer_rounds)) exit
            if (stats% max_change_total <= transfer_round_tol) exit
            ! the composition changed: re-evaluate the limits and rain faces ("status
            ! update") and transfer again
            call evaluate_miscibility(pd, nz, species, s% xa, A, i_sed, i_partner, s% Peos, s% T, &
               has_gap, in_gap, x_lo, x_hi, ierr)
            if (ierr /= 0) then
               write(*,'(a,i0,a)') ' element_sedimentation: evaluate_miscibility failed in transfer round ', &
                  n_rounds + 1, ' (' // trim(miscibility_error_message(ierr)) // ')'
               s% xa(1:species,1:nz) = xa_start(1:species,1:nz)
               call cleanup()
               return
            end if
            call set_limits_and_flow()
         end do
      else
         call adv_diff_integrate(p, x, total_time, opts, stats, ierr, after_step)
      end if
      max_change = maxval(abs(x(1:nz) - xa_base(i_sed,1:nz)))          ! what the rain did
      max_change_total = maxval(abs(x(1:nz) - xa_start(i_sed,1:nz)))   ! rain and mixing together
      last_attempt_model_number = s% model_number
      last_attempt_max_change = max_change_total
      last_attempt_dt = dt
      if (ierr /= 0) then
         call report_solver_failure(s, p, x, stats, ierr)
         s% xa(1:species,1:nz) = xa_start(1:species,1:nz)
         call cleanup()
         return
      end if

      ! --- damping (element_sedimentation_maximum_dXi_per_cell and the retry damping):
      ! if rain changes the sedimenting species in some cell by more than change_limit
      ! relative to the mixing-only state xa_base, move only part of the way,
      !    xa = xa_base + damping*(xa_rain - xa_base),  damping = change_limit/max_change,
      ! for all species. This is a ceiling on the change caused by rain in every cell,
      ! mixed or not (in unmixed cells xa_base is the composition at the start). A linear
      ! blend of two states that both conserve every species, sum to one and are
      ! non-negative has the same properties. Blending from the mixing-only state rather
      ! than from the start keeps mixed regions mixed: a layer that was stable and has
      ! become convective is homogenized as MESA's mixing would do, instead of keeping
      ! a fraction of its old (possibly inverted) composition jump. The rest of the rain
      ! happens in the following steps.
      ! On a MESA retry of the same model the whole change of the step (mixing included) is
      ! limited instead, relative to the composition at the start: MESA could not converge
      ! with the previous attempt, and the mixing-only state can itself be a large change
      ! (e.g. when the mixing coefficients of the start of the step would homogenize a
      ! helium ocean that the previous attempt had just formed), so a limit relative to it
      ! would not reduce what MESA sees.
      damping = 1d0
      if (change_limit > 0d0 .and. max_change > change_limit) then
         damping = change_limit/max_change
         s% xa(1:species,1:nz) = xa_base(1:species,1:nz) + &
            damping*(s% xa(1:species,1:nz) - xa_base(1:species,1:nz))
      end if
      if (retry_limit > 0d0) then
         max_change_total = maxval(abs(s% xa(i_sed,1:nz) - xa_start(i_sed,1:nz)))
         if (max_change_total > retry_limit) then
            damping = damping*retry_limit/max_change_total
            s% xa(1:species,1:nz) = xa_start(1:species,1:nz) + &
               (retry_limit/max_change_total)*(s% xa(1:species,1:nz) - xa_start(1:species,1:nz))
         end if
      end if
      last_attempt_max_change = maxval(abs(s% xa(i_sed,1:nz) - xa_start(i_sed,1:nz)))

      ! --- overturn of layers the transport left convectively unstable. The rain is
      ! deposited below the immiscible region with the mixing of the start of the step;
      ! where that mixing was zero (a layer that had become radiative) the deposit sits
      ! as a heavier layer on top of lighter material. The mixing coefficients of this
      ! step are fixed, so MESA could not mix it away before its convection treatment
      ! reacts to the inverted composition gradient (gradT -> gradL < grada, an inverted
      ! temperature gradient that pumps luminosity into the layer). Physically such a
      ! layer overturns on a dynamical timescale, so every run of cells whose faces are
      ! unstable by the Ledoux criterion (with B from the equation of state at the new
      ! composition) is homogenized here, every species mass-weighted. Blocks grow until
      ! the face above them is stable, e.g. down to the centre of a uniform interior.
      n_merged_faces = 0
      n_mixed_blocks = 0
      k_overturn_outer = 0
      k_overturn_inner = 0
      if (s% element_sedimentation_overturn_unstable_layers) then
         call overturn_unstable_layers(nz, species, s% dm, s% xa, ledoux_unstable_face, &
            n_merged_faces, n_mixed_blocks, k_overturn_outer, k_overturn_inner)
         if (ierr /= 0) then  ! set by ledoux_unstable_face if the EOS failed
            write(*,'(a)') ' element_sedimentation: overturn of unstable layers failed (equation of state)'
            s% xa(1:species,1:nz) = xa_start(1:species,1:nz)
            call cleanup()
            return
         end if
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

      ! --- MESA's composition mixing for the rest of the step: the transport above already
      ! mixed every species with D_mix together with the rain, so switch MESA's mixing off.
      ! D_mix is not recomputed within the step after this hook (okay_to_set_mixing_info is
      ! false), but do_struct_burn_mix rebuilds sig from D_mix (get_convection_sigmas), so
      ! D_mix itself is set to zero (sig follows); MESA recomputes it at the next step.
      if (s% element_sedimentation_replace_mesa_mixing) then
         s% D_mix(1:nz) = 0d0
         s% sig(1:nz) = 0d0
      else if (s% element_sedimentation_suppress_mixing_in_immiscible_regions) then
         ! --- optionally switch off MESA's mixing where the mixture is immiscible
         if (s% do_conv_premix) then
            if (.not. conv_premix_notice_printed) then
               write(*,'(a)') ' element_sedimentation: WARNING: convective premixing and element sedimentation ' // &
                  'with suppressed mixing are both on.'
               write(*,'(a)') '    Regions that rain out would flip-flop between homogeneous and rained-out mixing,'
               write(*,'(a)') '    so do_conv_premix is set to .false. for the rest of the run (this also switches off'
               write(*,'(a)') '    the convective premixing of the gentle_mixing extension, if installed).'
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
         write(*,'(a,i0,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,a,a,1pd10.3,a,1pd10.3,a,1pd9.2,a,1pd9.2,a,1pd9.2,a,1pd9.2,a,i0)') &
            ' element_sedimentation: model ', s% model_number, trim(merge(' rain', ' mix ', rain_active)), &
            ' n_gap ', n_gap, ' n_sat ', n_sat, ' (k ', k_first, '-', k_last, ') n_contrast ', n_contrast, &
            ' rounds ', n_rounds, ' steps ', stats% steps, &
            ' rejected ', stats% rejected_steps, &
            ' newton ', stats% newton_iters, ' stop ', trim(stop_reason_name(stats% stop_reason)), &
            ' max_dX ', maxval(abs(s% xa(i_sed,1:nz) - xa_start(i_sed,1:nz))), ' mass_moved[g] ', mass_moved, &
            ' cons_err ', cons_err, ' damping ', damping, ' passive_defect ', max_passive_defect, &
            ' max|sum-1| ', max_sum_dev, ' renormalized ', num_renormalized
         if (n_merged_faces > 0) write(*,'(a,i0,a,i0,a,i0,a,i0,a)') &
            ' element_sedimentation: overturn of unstable layers: ', n_merged_faces, ' faces mixed away in ', &
            n_mixed_blocks, ' block(s) (cells ', k_overturn_outer, '-', k_overturn_inner, ')'
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


      logical function ledoux_unstable_face(k, xa_above, xa_below)
         ! Face k (between cell k-1 above and cell k below) is convectively unstable by
         ! the Ledoux criterion, gradr > grada + B, for the given compositions of the
         ! two sides. B is MESA's composition term (brunt.f90, get_brunt_B): the change
         ! of ln P at the face (rho, T) from the composition of the upper to that of
         ! the lower side, over the pressure difference of the cells, divided by chiT.
         ! B < 0 for a heavier layer on top. gradr and grada are the face values of the
         ! start of the step. Host association provides s, nz, species and ierr; an EOS
         ! failure sets ierr and returns .false.
         integer, intent(in) :: k
         real(dp), intent(in) :: xa_above(:), xa_below(:)
         real(dp) :: alfa, rho_f, T_f, chiT_f, Prad_f, Ppoint, delta_lnP, lnP_above, lnP_below, &
            e_f, Pgas_f, B
         integer :: ierr_eos

         ledoux_unstable_face = .false.
         if (ierr /= 0) return
         if (k < 2 .or. k > nz) return
         if (maxval(abs(xa_above(1:species) - xa_below(1:species))) < overturn_composition_tol) return

         alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         rho_f = alfa*s% rho(k) + (1d0 - alfa)*s% rho(k-1)
         T_f = alfa*s% T(k) + (1d0 - alfa)*s% T(k-1)
         chiT_f = alfa*s% chiT(k) + (1d0 - alfa)*s% chiT(k-1)
         if (chiT_f <= 0d0 .or. rho_f <= 0d0 .or. T_f <= 0d0) return
         Prad_f = crad*T_f*T_f*T_f*T_f/3d0

         delta_lnP = s% lnPeos(k-1) - s% lnPeos(k)
         if (delta_lnP > -1d-50) then  ! hydrostatic estimate, as in MESA
            Ppoint = alfa*s% Peos(k) + (1d0 - alfa)*s% Peos(k-1)
            delta_lnP = -s% cgrav(k)*s% m(k)/(pi4*s% r(k)*s% r(k)*s% r(k)*s% r(k)*Ppoint)*s% dm_bar(k)
            if (delta_lnP >= 0d0) return
         end if

         call eos_at_rho_T(s, 0, xa_above, rho_f, T_f, e_f, Pgas_f, ierr_eos)
         if (ierr_eos /= 0) then
            ierr = ierr_eos
            write(*,'(a,i0)') ' element_sedimentation: EOS failed at the face of the overturn test, k = ', k
            return
         end if
         lnP_above = log(Prad_f + Pgas_f)
         call eos_at_rho_T(s, 0, xa_below, rho_f, T_f, e_f, Pgas_f, ierr_eos)
         if (ierr_eos /= 0) then
            ierr = ierr_eos
            write(*,'(a,i0)') ' element_sedimentation: EOS failed at the face of the overturn test, k = ', k
            return
         end if
         lnP_below = log(Prad_f + Pgas_f)

         B = (lnP_below - lnP_above)/delta_lnP/chiT_f
         if (is_bad(B)) return
         ledoux_unstable_face = (B < 0d0 .and. s% gradr(k) > s% grada_face(k) + B)
      end function ledoux_unstable_face


      subroutine set_limits_and_flow()
         ! floors and caps from the miscibility limits, carrier flow on the faces whose
         ! donor cell has a gap (cells without a gap never rain: their floor is
         ! adv_diff_huge, so the flux would be zero anyway, and such faces are left
         ! inactive and do not enter the timescale bookkeeping of the solver). In a
         ! mixing-only step (no saturated cell) nothing rains.
         ! The instant methods have no velocity: a face rains if the cell above it has a
         ! gap (the species always sinks), and p% w only marks those faces. In a
         ! mixing-only step no limits are installed at all (the instant transfer would
         ! otherwise clamp gap cells to their floors without any rain face).
         integer :: kk, dn
         p% w(1:nz+1) = 0d0
         do kk = 2, nz
            if (v0 > 0d0 .or. instant_method) then
               dn = kk - 1
            else
               dn = kk
            end if
            if (rain_active .and. has_gap(dn)) then
               if (instant_method) then
                  p% w(kk) = 1d0
               else
                  p% w(kk) = w_face(kk)
               end if
            end if
         end do
         do kk = 1, nz
            if (has_gap(kk) .and. rain_active) then
               p% x_floor(kk) = x_lo(kk)   ! only the species above the He-poor limit rains out
               p% x_cap(kk) = x_hi(kk)     ! inflow stops at the He-rich limit
            else
               p% x_floor(kk) = adv_diff_huge
               p% x_cap(kk) = 1d0
            end if
         end do
      end subroutine set_limits_and_flow


      logical function mixed_face(kk, xa_above, xa_below)
         ! for overturn_unstable_layers: a face with mixing joins its two cells
         integer, intent(in) :: kk
         real(dp), intent(in) :: xa_above(:), xa_below(:)
         mixed_face = (p% sig(kk) > 0d0)
         if (.false.) mixed_face = mixed_face .or. (size(xa_above) + size(xa_below) < 0)
      end function mixed_face


      subroutine cleanup()
         call free_adv_diff_problem(p)
         call free_adv_diff_problem(pp)
         call free_adv_diff_problem(p_mix)
      end subroutine cleanup

   end subroutine do_advection_diffusion_rain


   subroutine change_limits_for_this_attempt(s, user_limit, retry_limit)
      ! Largest allowed changes of the mass fraction of the sedimenting species in a cell
      ! for this attempt (applied by damping the step), <= 0 when there is no limit:
      ! user_limit = element_sedimentation_maximum_dXi_per_cell (relative to the mixing-only
      ! state of the step), and, when MESA is retrying the same model after an attempt
      ! that transported material, retry_limit = element_sedimentation_retry_change_factor
      ! times the largest change (relative to the start of the step) applied in that
      ! attempt, times the ratio of the new to the previous time step. A steady-state
      ! solve does not depend on the time step, so without the retry limit every retry
      ! would reproduce the composition change that made the step fail; and since the
      ! heating of the step is the change divided by dt, a limit that only halves the
      ! change while MESA halves dt would leave the heating rate unchanged, retry after
      ! retry. Scaling with dt makes the change, and so the heating rate, shrink with
      ! the time step.
      type (star_info), pointer :: s
      real(dp), intent(out) :: user_limit, retry_limit
      logical :: is_retry
      real(dp) :: dt_ratio

      user_limit = s% element_sedimentation_maximum_dXi_per_cell
      retry_limit = -1d0

      is_retry = s% model_number_for_last_retry == s% model_number .and. &
         last_attempt_model_number == s% model_number .and. &
         last_attempt_max_change > 0d0 .and. &
         s% element_sedimentation_retry_change_factor > 0d0
      if (.not. is_retry) return

      dt_ratio = 1d0
      if (last_attempt_dt > 0d0 .and. s% dt > 0d0) dt_ratio = min(1d0, s% dt/last_attempt_dt)
      retry_limit = s% element_sedimentation_retry_change_factor*dt_ratio*last_attempt_max_change
      write(*,'(a,i0,a,1pd10.3,a,1pd10.3,a,1pd10.3,a)') &
         ' element_sedimentation: model ', s% model_number, &
         ' is a retry; the change per cell of the whole step is limited to ', retry_limit, &
         ' (previous attempt: ', last_attempt_max_change, ', dt ratio ', dt_ratio, &
         '; see element_sedimentation_retry_change_factor)'
   end subroutine change_limits_for_this_attempt


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


   ! ---------------------------------------------------------------------------
   ! energy carried by the moved species
   ! ---------------------------------------------------------------------------

   subroutine enthalpy_flux_heating(s, dt, xa_old, ierr)
      ! Heating rate of every cell from the divergence of the enthalpy flux carried
      ! by the species the step moved (sedimentation_energy): the mass flow rates
      ! through the faces follow from the composition change, the partial specific
      ! enthalpies from finite differences of the EOS at the (rho, T) of the cells
      ! (values at the start of the step, where the composition change is applied),
      ! averaged to the faces with the cell masses as weights. The result is kept in
      ! eps_enthalpy_flux for finish_element_sedimentation.
      type (star_info), pointer :: s
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: xa_old(:,:)  ! (species, nz) composition at entry of the step
      integer, intent(out) :: ierr
      integer :: nz, species, j, k, k_top, k_bottom, k_max
      real(dp), allocatable :: flow(:,:), h_cell(:,:), h_face(:,:), energy_flux(:)
      logical, allocatable :: needed(:,:)
      real(dp) :: dh, total
      include 'formats'

      ierr = 0
      nz = s% nz
      species = s% species
      enthalpy_flux_valid = .false.
      if (allocated(eps_enthalpy_flux)) then
         if (size(eps_enthalpy_flux) < nz) deallocate(eps_enthalpy_flux)
      end if
      if (.not. allocated(eps_enthalpy_flux)) allocate(eps_enthalpy_flux(nz))
      eps_enthalpy_flux(1:nz) = 0d0

      allocate(flow(species,nz+1), h_cell(species,nz), h_face(species,nz+1), energy_flux(nz+1), &
         needed(species,nz))
      call species_face_flow_rates(nz, species, s% dm, xa_old, s% xa, dt, flow)
      if (maxval(abs(flow(1:species,1:nz+1))) == 0d0) then
         enthalpy_flux_valid = .true.
         return
      end if

      ! partial enthalpies relative to the partner species: the flow rates of a face add
      ! up to zero, so the reference drops out; needed only in the cells next to a face
      ! with a flow of that species
      needed = .false.
      do k = 2, nz
         do j = 1, species
            if (j == i_partner .or. flow(j,k) == 0d0) cycle
            needed(j,k-1) = .true.
            needed(j,k) = .true.
         end do
      end do
      h_cell = 0d0
      do k = 1, nz
         do j = 1, species
            if (.not. needed(j,k)) cycle
            call partial_enthalpy(s, k, j, xa_old(1:species,k), dh, ierr)
            if (ierr /= 0) return
            h_cell(j,k) = dh
         end do
      end do
      h_face = 0d0
      do k = 2, nz
         do j = 1, species
            if (j == i_partner) cycle
            h_face(j,k) = (s% dm(k)*h_cell(j,k-1) + s% dm(k-1)*h_cell(j,k))/(s% dm(k-1) + s% dm(k))
         end do
      end do

      call face_energy_flux(nz, species, flow, h_face, energy_flux)
      call flux_divergence_heating(nz, s% dm, energy_flux, eps_enthalpy_flux)
      do k = 1, nz
         if (is_bad(eps_enthalpy_flux(k))) then
            ierr = 1
            write(*,2) ' element_sedimentation: NaN or Inf in the enthalpy flux heating of cell', k, eps_enthalpy_flux(k)
            return
         end if
      end do
      enthalpy_flux_valid = .true.

      if (s% element_sedimentation_report) then
         total = dot_product(s% dm(1:nz), eps_enthalpy_flux(1:nz))
         write(*,'(a,1pd10.3,a,1pd10.3,a,1pd10.3,a,1pd10.3,a,1pd10.3)') &
            ' element_sedimentation: enthalpy flux of the moved species: max |E| [erg/s] ', &
            maxval(abs(energy_flux(1:nz+1))), ' eps min ', minval(eps_enthalpy_flux(1:nz)), &
            ' max ', maxval(eps_enthalpy_flux(1:nz)), ' [erg/g/s] sum(dm eps) ', total, &
            ' /(dm max|eps|) ', total/max(maxval(s% dm(1:nz)*abs(eps_enthalpy_flux(1:nz))), tiny(total))
         ! partial enthalpy difference of the sedimenting species and its partner in units
         ! of kT/m_p at the outermost and innermost face the species crosses and at the
         ! face with the largest energy flux (an ideal H2/He gas gives -9/8 for He - H;
         ! the dense H/He fluid of a giant planet several times that)
         k_top = 0
         k_bottom = 0
         do k = 2, nz
            if (flow(i_sed,k) == 0d0) cycle
            if (k_top == 0) k_top = k
            k_bottom = k
         end do
         k_max = maxloc(abs(energy_flux(2:nz)), dim=1) + 1
         if (k_top > 0) write(*,'(a,i0,a,1pd10.3,a,i0,a,1pd10.3,a,i0,a,1pd10.3)') &
            ' element_sedimentation: (h_' // trim(chem_isos% name(s% chem_id(i_sed))) // ' - h_' // &
            trim(chem_isos% name(s% chem_id(i_partner))) // ')/(kT/m_p): outermost crossed face ', k_top, ' ', &
            h_face(i_sed,k_top)/(boltzm*face_T(k_top)/mp), ', innermost ', k_bottom, ' ', &
            h_face(i_sed,k_bottom)/(boltzm*face_T(k_bottom)/mp), ', face of max |E| ', k_max, ' ', &
            h_face(i_sed,k_max)/(boltzm*face_T(k_max)/mp)
      end if

   contains

      real(dp) function face_T(k)
         integer, intent(in) :: k
         face_T = (s% dm(k)*s% T(k-1) + s% dm(k-1)*s% T(k))/(s% dm(k-1) + s% dm(k))
      end function face_T

   end subroutine enthalpy_flux_heating


   subroutine partial_enthalpy(s, k, j, xa, dh, ierr)
      ! Partial specific enthalpy of species j relative to the partner species in cell
      ! k, (dh/dX_j)_{T,P} with X_partner compensating, from a finite difference of the
      ! EOS at fixed (rho, T) around the composition xa (sedimentation_energy for the
      ! conversion to fixed pressure). Central differences of size dxa_eos where both
      ! species allow it, one-sided ones at the edges of the composition range.
      type (star_info), pointer :: s
      integer, intent(in) :: k, j
      real(dp), intent(in) :: xa(:)  ! (species)
      real(dp), intent(out) :: dh
      integer, intent(out) :: ierr
      real(dp) :: xa_plus(s% species), xa_minus(s% species), delta_plus, delta_minus, &
         e_plus, e_minus, Pgas_plus, Pgas_minus
      include 'formats'

      ierr = 0
      dh = 0d0

      ! steps towards more (plus) and less (minus) of species j, each taken from the
      ! partner; limited by the abundances so that no mass fraction becomes negative
      delta_plus = min(dxa_eos, xa(i_partner))
      delta_minus = min(dxa_eos, xa(j))
      if (delta_plus + delta_minus < dxa_eos_min) return  ! both species absent: nothing moves here

      xa_plus(1:s% species) = xa(1:s% species)
      xa_plus(j) = xa(j) + delta_plus
      xa_plus(i_partner) = xa(i_partner) - delta_plus
      xa_minus(1:s% species) = xa(1:s% species)
      xa_minus(j) = xa(j) - delta_minus
      xa_minus(i_partner) = xa(i_partner) + delta_minus

      call eos_e_and_Pgas(s, k, xa_plus, e_plus, Pgas_plus, ierr)
      if (ierr /= 0) return
      call eos_e_and_Pgas(s, k, xa_minus, e_minus, Pgas_minus, ierr)
      if (ierr /= 0) return

      dh = partial_enthalpy_increment((e_plus - e_minus)/(delta_plus + delta_minus), &
         (Pgas_plus - Pgas_minus)/(delta_plus + delta_minus), s% rho(k), s% chiT(k), s% chiRho(k))
      if (is_bad(dh)) then
         ierr = 1
         write(*,3) ' element_sedimentation: NaN or Inf partial enthalpy for species, cell', j, k, dh
      end if

   end subroutine partial_enthalpy


   subroutine eos_e_and_Pgas(s, k, xa, e, Pgas, ierr)
      ! specific internal energy and gas pressure of cell k at its (rho, T) for the composition xa
      type (star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(in) :: xa(:)
      real(dp), intent(out) :: e, Pgas
      integer, intent(out) :: ierr
      include 'formats'
      call eos_at_rho_T(s, k, xa, s% rho(k), s% T(k), e, Pgas, ierr)
      if (ierr /= 0) write(*,2) ' element_sedimentation: EOS failed for a perturbed composition in cell', k
   end subroutine eos_e_and_Pgas


   subroutine eos_at_rho_T(s, k, xa, rho, T, e, Pgas, ierr)
      ! specific internal energy and gas pressure at (rho, T) for the composition xa;
      ! k is the cell the call belongs to (0: a face or no particular cell)
      type (star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(in) :: xa(:), rho, T
      real(dp), intent(out) :: e, Pgas
      integer, intent(out) :: ierr
      real(dp), dimension(num_eos_basic_results) :: res, dres_dlnd, dres_dlnT
      real(dp) :: dres_dxa(num_eos_d_dxa_results, s% species)
      e = 0d0
      Pgas = 0d0
      call get_eos(s, k, xa, rho, log10(rho), T, log10(T), res, dres_dlnd, dres_dlnT, dres_dxa, ierr)
      if (ierr /= 0) return
      e = exp(res(i_lnE))
      Pgas = exp(res(i_lnPgas))
   end subroutine eos_at_rho_T

end module element_sedimentation
