! ***********************************************************************
!   Copyright (C) 2025  Henrik Knierim
! ***********************************************************************
module element_sedimentation

   use star_lib
   use const_def
   use chem_lib

   implicit none

   private

   public :: init_element_sedimentation
   public :: shutdown_phase_diagram
   public :: do_element_sedimentation

   integer, parameter :: instant_rain = 0
   integer, parameter :: bottom_up_rain = 1
   integer, parameter :: advection_diffusion_rain = 2

   real(dp), parameter :: very_tiny = 1d-12
   real(dp), parameter :: tiny = 1d-10
   real(dp), parameter :: small = 1d-4
   real(dp), parameter :: sedimentation_threshold = 1d-4   ! ignore sedimentation if the change in the mass fraction is smaller than this value
   real(dp), parameter :: temperature_shift = 0d0

   integer :: i_Xi = -1
   logical :: debug = .false.
   logical :: debug_advection_diffusion_solver = .false.
   logical :: debug_phase_diagram = .false.
   logical :: verbose = .false.
   logical :: extrapolate = .false.
   logical :: damp_sedimentation ! originally, damping was determined via do_gentle_mixing. Since gentle_mixing is not officially supported, I decided to make this module self-contained.

   ! TODO: make this a stand-alone module that can call other star modules
   ! TODO: once done, implement the _update method.

   type phase_diagram
      character(len=1024) :: phase_diagram_filename
      integer :: n_points
      real(kind=dp), allocatable :: He_rich(:), He_poor(:), T(:), P(:)
      integer :: num_unique_Ps
      integer, allocatable :: unique_P_index(:), num_Ts_at_unique_Ps(:)
      real(kind=dp), allocatable :: unique_P(:)
   end type phase_diagram

   type immiscibility_region
      integer :: klo, khi  ! the precipitating region
      integer :: k_reservoir_lo, k_reservoir_hi  ! the convective element reservoir
      integer :: k_deposite_lo, k_deposite_hi  ! the region where the element is deposited
      real(kind=dp), allocatable :: dXi(:) ! change in the mass fractions
   end type immiscibility_region

   type advection_diffusion_state
      integer :: nz  ! number of zones
      real(kind=dp) :: dt  ! time step
      real(kind=dp), allocatable :: D(:), v(:)
      real(kind=dp), allocatable :: dr(:), Xi(:), dXi(:)  ! dXi is the change in the element's mass fraction
   end type advection_diffusion_state

   type(phase_diagram) :: immiscibility_phase_diagram
   logical :: phase_diagram_loaded = .false.


contains

   subroutine init_element_sedimentation(s, ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr

      ierr = 0 ! 0 means AOK

      ! check if the phase diagram is already loaded
      ! if so, throw an error. This rountine should only be called once
      if (phase_diagram_loaded) then
         ierr = 1
         write(*,*) 'init_element_sedimentation: phase diagram already loaded'
         return
      end if

      ! load the phase diagram
      immiscibility_phase_diagram%phase_diagram_filename = trim(s% element_sedimentation_phase_diagram)
      call load_phase_diagram(immiscibility_phase_diagram, ierr)
      phase_diagram_loaded = .true.
      if (ierr /= 0) return

      extrapolate = s% element_sedimentation_extrapolate_phase_diagram

      if (debug_phase_diagram) call print_phase_diagram(immiscibility_phase_diagram)
      
      ! set index of sedimenting species
      i_Xi = s% net_iso(chem_get_iso_id(trim(s% element_sedimentation_chemical_species)))

   end subroutine init_element_sedimentation

   subroutine do_element_sedimentation(s, ierr)

      type (star_info), pointer :: s
      integer :: method ! 0: instant rain, 1: constant velocity rain
      integer, intent(out) :: ierr

      ierr = 0 ! 0 means AOK

      ! get the method to use for the element sedimentation
      ! we first assign the method and then select the correct subroutine
      ! to be able to reference the method in the subroutine later if needed (currently not the case)
      if (trim(s% element_sedimentation_method) == 'instant_rain') then
         method = instant_rain
      else if (trim(s% element_sedimentation_method) == 'bottom_up_rain') then
         method = bottom_up_rain
      else if (trim(s% element_sedimentation_method) == 'advection_diffusion_rain') then
         method = advection_diffusion_rain
      else
         ierr = 1
         write(*,*) 'do_element_sedimentation: unknown method for element sedimentation: ', trim(s% element_sedimentation_method)
         return
      end if

      damp_sedimentation = (s% element_sedimentation_maximum_dXi_per_cell > 0d0)

      select case (method)
       case (instant_rain)
         call do_instant_rain(s, ierr)
       case (bottom_up_rain)
         call do_bottom_up_rain(s, ierr)
       case (advection_diffusion_rain)
         call do_advection_diffusion_rain(s, ierr)
       case default
         ierr = 1
         write(*,*) 'do_element_sedimentation: unknown method for element sedimentation: ', trim(s% element_sedimentation_method)
         return
      end select

   end subroutine do_element_sedimentation

   subroutine do_advection_diffusion_rain(s, ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr

      ! local
      integer, parameter :: max_iterations = 100
      integer :: i, j
      type(advection_diffusion_state) :: adv_diff_state

      ierr = 0 ! 0 means AOK


      do i = 1, max_iterations
         if (verbose) write(*,*) 'do_advection_diffusion_rain: iteration: ', i

         ! initialize the advection-diffusion state
         call init_advection_diffusion_state(s, adv_diff_state, ierr)
         if (ierr /= 0) return

         ! if adv_diff_state % v is not allocated, no precipitating regions were found; exit the loop
         if (.not. allocated(adv_diff_state % v)) then
            if (debug_advection_diffusion_solver) write(*,*) 'do_advection_diffusion_rain: no precipitating regions found, skipping advection-diffusion step'
            call free_advection_diffusion_state(adv_diff_state)
            exit
         end if

         ! perform the advection-diffusion step
         call backward_Euler_step(adv_diff_state, ierr)
         if (ierr /= 0) return
         
         if (debug_advection_diffusion_solver) then
            do j = 1, adv_diff_state%nz
               write(*,*) 'do_advection_diffusion_rain: (k, Xi(k), dXi(k)) = ', j, adv_diff_state%Xi(j), adv_diff_state%dXi(j)
            end do
         end if

         call adjust_mass_fractions(s, adv_diff_state% dXi, ierr)

         ! check if the change in the element's mass fraction is below a threshold
         if (all(abs(adv_diff_state%dXi) < sedimentation_threshold)) then
            if (verbose) write(*,*) 'do_advection_diffusion_rain: no significant change in element mass fraction'
            if (verbose) write(*,*) 'do_advection_diffusion_rain: max(|dXi|) = ', maxval(abs(adv_diff_state%dXi))
            exit
         end if

         call free_advection_diffusion_state(adv_diff_state)
      end do

      ! don't mix afterwards
      s% D_mix(1:s% nz) = 0d0


   end subroutine do_advection_diffusion_rain

   subroutine init_advection_diffusion_state(s, adv_diff_state, ierr)
      type (star_info), pointer :: s
      type(advection_diffusion_state), intent(out) :: adv_diff_state
      integer, intent(out) :: ierr

      ! local
      real(kind=dp), allocatable :: Xi_miscible(:)  ! miscible element mass fraction
      logical, allocatable :: is_miscible(:)  ! miscible regions
      real(kind=dp), allocatable :: v_settle(:)  ! settling velocity
      integer :: k

      ierr = 0 ! 0 means AOK

      ! function to populate a advection_diffusion_state type used for the advection-diffusion solver
      adv_diff_state % nz = s% nz
      allocate(adv_diff_state % D(1:s% nz), source=s% D_mix(1: s% nz))  ! diffusion coefficient
      allocate(adv_diff_state % Xi(1:s% nz), source=s% xa(i_Xi, 1: s% nz))  ! diffusion coefficient
      allocate(adv_diff_state % dXi(1:s% nz), source=0d0)  ! diffusion coefficient

      ! radial step size
      allocate(adv_diff_state % dr(1:s% nz))  
      !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(guided)
      do k = 1, s % nz
         adv_diff_state % dr(k) = get_dr(s, k)
      end do
      !$OMP END PARALLEL DO

      call get_raw_miscibility(s, is_miscible, Xi_miscible,  ierr)
      if (ierr /= 0) then
         write(*,*) 'init_advection_diffusion_state: error in get_raw_miscibility'
         return
      end if

      ! if everything is miscible, we can skip the advection-diffusion step
      if (all(is_miscible)) then
         if (debug_advection_diffusion_solver) write(*,*) 'init_advection_diffusion_state: all regions are miscible, skipping advection-diffusion step'
         return
      end if
      
      call get_settling_velocity(s, adv_diff_state % v, is_miscible, Xi_miscible, ierr)
      if (ierr /= 0) then
         write(*,*) 'init_advection_diffusion_state: error in get_settling_velocity'
         return
      end if

      call get_time_step(s, adv_diff_state, ierr)
      if (ierr /= 0) then
         write(*,*) 'init_advection_diffusion_state: error in get_time_step'
         return
      end if

      if (debug_advection_diffusion_solver) then
         write(*,*) 'init_advection_diffusion_state: advection-diffusion state initialized'
         write(*,*) 'init_advection_diffusion_state: nz = ', adv_diff_state % nz
         write(*,*) 'init_advection_diffusion_state: dt = ', adv_diff_state % dt
         do k = 1, adv_diff_state % nz
            write(*,'(a,i12,99(1pd26.16))') 'init_advection_diffusion_state: (k, D(k), v(k), dr(k), Xi(k)) = ', k, adv_diff_state % D(k), adv_diff_state % v(k), adv_diff_state % dr(k), adv_diff_state % Xi(k)
         end do
      end if

      deallocate(Xi_miscible, is_miscible)
      
   end subroutine init_advection_diffusion_state

   subroutine get_time_step(s, adv_diff_state, ierr)
      type (star_info), pointer :: s
      type(advection_diffusion_state), intent(inout) :: adv_diff_state
      integer, intent(out) :: ierr

      ! local
      real(kind=dp) :: tau_adv ! advection timescale
      real(kind=dp) :: tau_diff ! diffusion timescale
      real(kind=dp) :: tau_adv_div_tau_diff  ! ratio of advection to diffusion timescale
      real(kind=dp) :: L, v, D
      logical, allocatable :: mask(:)
      real(kind=dp) :: fraction_of_timescale
      
      ierr = 0 ! 0 means AOK
      
      fraction_of_timescale = s% element_sedimentation_advection_diffusion_time_step_factor
      ! HHK: Check beta. To see some change at all, I had to set this value to 1d7 last time.

      allocate(mask(adv_diff_state%nz))
      mask = (adv_diff_state%v /= 0d0)

      v = dot_product(pack(adv_diff_state% v, mask), pack(adv_diff_state% dr, mask))
      D = dot_product(pack(adv_diff_state% D, mask), pack(adv_diff_state% dr, mask))
      L = sum(pack(adv_diff_state% dr, mask))

      if (D <= 0d0) then
         ierr = 1
         write(*,*) 'get_time_step: diffusion coefficient is zero or negative'
         return
      end if
      if (v <= 0d0) then
         ierr = 1
         write(*,*) 'get_time_step: settling velocity is zero or negative'
         return
      end if
      if (L <= 0d0) then
         ierr = 1
         write(*,*) 'get_time_step: radial step size is zero or negative'
         return
      end if

      tau_adv = L/v
      tau_diff = L**2 / D
      tau_adv_div_tau_diff = tau_adv/tau_diff
      adv_diff_state% dt = min(tau_adv, tau_diff) * fraction_of_timescale

      if (debug_advection_diffusion_solver) then
         write(*,*) 'get_time_step: tau_adv = ', tau_adv
         write(*,*) 'get_time_step: tau_diff = ', tau_diff
         write(*,*) 'get_time_step: tau_adv_div_tau_diff = ', tau_adv_div_tau_diff
         write(*,*) 'get_time_step: dt = ', adv_diff_state% dt
      end if

   
      deallocate(mask)

   end subroutine get_time_step

   subroutine free_advection_diffusion_state(adv_diff_state)
      type(advection_diffusion_state), intent(inout) :: adv_diff_state

      ! free the allocated arrays
      if (allocated(adv_diff_state % D)) deallocate(adv_diff_state % D)
      if (allocated(adv_diff_state % v)) deallocate(adv_diff_state % v)
      if (allocated(adv_diff_state % dr)) deallocate(adv_diff_state % dr)
      if (allocated(adv_diff_state % Xi)) deallocate(adv_diff_state % Xi)
      if (allocated(adv_diff_state % dXi)) deallocate(adv_diff_state % dXi)

      ! reset the state
      adv_diff_state % nz = 0
      adv_diff_state % dt = 0d0

      if (debug_advection_diffusion_solver) then
         write(*,*) 'free_advection_diffusion_state: advection-diffusion state freed'
      end if
   
      
   end subroutine free_advection_diffusion_state

   subroutine backward_Euler_step(adv_diff_state, ierr)
      type(advection_diffusion_state), intent(inout) :: adv_diff_state
      integer, intent(out) :: ierr

      ! local variables
      integer :: i, nz, k
      real(kind=dp), allocatable :: a(:), b(:), c(:), f_n(:), f_np1(:)

      ierr = 0 ! 0 means AOK
      nz = adv_diff_state % nz

      ! use tridiagonal matrix solver in the range [2, N-1]
      ! the other values are set by boundary conditions
      allocate(a(nz-2), b(nz-2), c(nz-2), f_n(nz-2), f_np1(nz-2))

      ! Set tridiagonal matrix coefficients
      ! a(1) and c(-1) are not used
      call setup_tridiagonal(adv_diff_state, a, b, c)  ! hhk
      
      f_n(1:nz-2) = adv_diff_state% Xi(2:nz-1)
      call tridag(a, b, c, f_n, f_np1)
      
      ! set change in element abundance
      adv_diff_state% dXi(1) = f_np1(1) - f_n(1) ! Neumann BC at r=0
      adv_diff_state% dXi(2:nz-1) = f_np1(1:nz-2) - f_n(1:nz-2) ! interior points
      adv_diff_state% dXi(nz) = f_np1(nz-2) - f_n(nz-2) ! Neumann BC at r=L

      if ( debug_advection_diffusion_solver ) then
         write(*,'(a16,a16,a16,a16)') 'k','f_n(k)', 'f_np1(k)', 'dXi(k)'
         do k = 1, nz-2
            write(*,'(i12,99(1pd26.16))') k+1, f_n(k), f_np1(k), adv_diff_state% dXi(k+1)
         end do  
      end if

      deallocate(a, b, c, f_n, f_np1)
   end subroutine backward_Euler_step

   subroutine do_bottom_up_rain(s,  ierr)
      type (star_info), pointer :: s
      integer :: i
      integer, intent(out) :: ierr

      integer, parameter :: max_iterations = 10000
      integer :: j
      integer :: klo, khi, k_mix, k

      ! debugging
      logical :: debug = .false.
      integer :: k_debug
      real(kind=dp), allocatable :: xa_pre(:,:)
      real(kind=dp) :: X, Y, Z

      type(immiscibility_region), allocatable :: precipitating_regions(:)

      include 'formats'

      ierr = 0 ! 0 means AOK

      allocate(xa_pre(size(s%xa, 1), size(s%xa, 2)))
      
      do j=1, max_iterations
         
         if (verbose) write(*,*) 'do_bottom_up_rain: iteration: ', j
         xa_pre = s% xa

         ! determine the precipitating region
         ! if allocated, deallocate precipitating regions
         if (allocated(precipitating_regions)) deallocate(precipitating_regions)
         call determine_precipitating_regions(s, precipitating_regions, ierr)
         if (ierr /= 0) return
         if (.not. allocated(precipitating_regions)) then
            if (verbose) write(*,1) 'do_bottom_up_rain: no precipitating region found.'
            return
         end if



         ! determine the upper convective element reservoir
         call determine_top_reservoirs(s, precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_bottom_up_rain: error in determine_top_reservoirs'
            return
         endif

         ! determine deposit region
         call determine_deposites(s, precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_bottom_up_rain: error in determine_deposites'
            return
         endif

         call check_precipitating_regions(precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_bottom_up_rain: error in check_precipitating_regions'
            return
         endif

         do i = 1, size(precipitating_regions)

            ! the differnce to instant rain: we change the abundances cell-by-cell
            klo = precipitating_regions(i)% klo
            khi = precipitating_regions(i)% khi
            
            ! as long as we are raining down, we don't need to determine the top reservoir
            ! case 1: We are inside the raining region (klo /= khi)
            ! Set all changes in dXi except the lowest cell equal to zero
            ! no need to adjust the top reservoir (yet)
            if (klo /= khi) then
               write(*,*) 'precipitating_regions(i)%dXi(klo:khi): ', precipitating_regions(i)%dXi(klo:khi)
               k_mix = khi
               do k = khi, klo, -1
                  if (abs(precipitating_regions(i)%dXi(k)) < sedimentation_threshold) then
                     precipitating_regions(i)%dXi(k) = 0d0
                     k_mix = k - 1
                  else
                     exit
                  end if
               end do
               precipitating_regions(i)%dXi(1:k_mix-1) = 0d0
               write(*,*) 'precipitating_regions(i)%dXi(klo:khi): ', precipitating_regions(i)%dXi(klo:khi)
            end if

            call set_top_reservoir_abundance_change(s,  precipitating_regions(i), ierr)
            if (ierr /= 0) then
               write(*,*) 'do_bottom_up_rain: error in set_top_reservoir_abundance_change'
               return
            endif

            call set_deposit_region_abundance_change(s,  precipitating_regions(i), ierr)
            if (ierr /= 0) then
               write(*,*) 'do_bottom_up_rain: error in set_deposit_region_abundance_change'
               return
            endif

            if ( debug ) then
               do k_debug = 1, s%nz
                  call get_mass_fractions(s, s% xa(:,k_debug), X, Y, Z)
                  write(*,2) 'do_bottom_up_rain: (k, X, Y, Z, dXi): ', &
                     k_debug, X, Y, Z, &
                     precipitating_regions(i)%dXi(k_debug)
               end do
            end if

            call adjust_mass_fractions(s, precipitating_regions(i)% dXi, ierr)
            if (ierr /= 0) then
               write(*,*) 'do_bottom_up_rain: error in adjust_mass_fractions'
               return
            end if

            if ( debug ) then
               write(*,*) 'do_bottom_up_rain: after setting the abundances:'
               do k_debug = 1, s%nz
                  call get_mass_fractions(s, s% xa(:,k_debug), X, Y, Z)
                  write(*,2) 'do_bottom_up_rain: (k, X, Y, Z): ', &
                     k_debug, X, Y, Z
               end do
            end if

            ! check mass fraction conservation
            if ( debug ) then
               call monitor_mass_fractions(s,  precipitating_regions(i), ierr)
            end if

         end do
         
         if ( verbose ) then
            k_debug = maxloc(abs(s%xa(i_Xi, 1:s%nz)-xa_pre(i_Xi, 1:s%nz)), dim=1)
            if (verbose) write(*,*) 'maximum difference: (k, val)', k_debug,  abs(s%xa(i_Xi, k_debug)-xa_pre(i_Xi, k_debug))
         end if
         
         if (all(abs(s%xa(i_Xi, 1:s%nz)-xa_pre(i_Xi,1:s%nz)) < small)) exit

      end do

      call check_after_rain(s, precipitating_regions, xa_pre, ierr)
      deallocate(xa_pre, precipitating_regions)
   
   end subroutine do_bottom_up_rain

   subroutine do_instant_rain(s, ierr)
      type (star_info), pointer :: s
      integer :: i
      integer, intent(out) :: ierr

      integer, parameter :: max_iterations = 1000
      integer :: j

      ! debugging
      logical :: debug = .false.
      integer :: k_debug
      real(kind=dp), allocatable :: xa_pre(:,:)
      real(kind=dp) :: X, Y, Z

      type(immiscibility_region), allocatable :: precipitating_regions(:)

      include 'formats'

      ierr = 0 ! 0 means AOK

      allocate(xa_pre(size(s%xa, 1), size(s%xa, 2)))

      do j=1, max_iterations
         
         if (verbose) write(*,*) 'do_instant_rain: iteration: ', j
         xa_pre = s%xa

         ! determine the precipitating region
         call determine_precipitating_regions(s, precipitating_regions, ierr)
         if (ierr /= 0) return
         if (.not. allocated(precipitating_regions)) then
            if (verbose) write(*,1) 'do_instant_rain: no precipitating region found.'
            return
         end if

         ! determine the upper convective element reservoir
         call determine_top_reservoirs(s, precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_instant_rain: error in determine_top_reservoirs'
            return
         endif

         ! determine deposit region
         call determine_deposites(s, precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_instant_rain: error in determine_deposites'
            return
         endif

         call check_precipitating_regions(precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_instant_rain: error in check_precipitating_regions'
            return
         endif

         do i = 1, size(precipitating_regions)

            call set_top_reservoir_abundance_change(s,  precipitating_regions(i), ierr)
            if (ierr /= 0) then
               write(*,*) 'do_instant_rain: error in set_top_reservoir_abundance_change'
               return
            endif

            call set_deposit_region_abundance_change(s,  precipitating_regions(i), ierr)
            if (ierr /= 0) then
               write(*,*) 'do_instant_rain: error in set_deposit_region_abundance_change'
               return
            endif

            if ( debug ) then
               do k_debug = 1, s%nz
                  call get_mass_fractions(s, s% xa(:,k_debug), X, Y, Z)
                  write(*,2) 'do_instant_rain: (k, X, Y, Z, dXi): ', &
                     k_debug, X, Y, Z, &
                     precipitating_regions(i)%dXi(k_debug)
               end do
            end if

            call adjust_mass_fractions(s, precipitating_regions(i)% dXi, ierr)
            if (ierr /= 0) then
               write(*,*) 'do_bottom_up_rain: error in adjust_mass_fractions'
               return
            end if

            if ( debug ) then
               write(*,*) 'do_instant_rain: after setting the abundances:'
               do k_debug = 1, s%nz
                  call get_mass_fractions(s, s% xa(:,k_debug), X, Y, Z)
                  write(*,2) 'do_instant_rain: (k, X, Y, Z): ', &
                     k_debug, X, Y, Z
               end do
            end if

            ! check mass fraction conservation
            if ( debug ) then
               call monitor_mass_fractions(s,  precipitating_regions(i), ierr)
            end if

         end do

         call check_after_rain(s, precipitating_regions, xa_pre, ierr)
         deallocate(precipitating_regions)

         if ( verbose ) then
            k_debug = maxloc(abs(s%xa(i_Xi, 1:s%nz)-xa_pre(i_Xi, 1:s%nz)), dim=1)
            if (verbose) write(*,*) 'maximum difference: (k, val)', k_debug,  abs(s%xa(i_Xi, k_debug)-xa_pre(i_Xi, k_debug))
         end if
         if (all(abs(s%xa(i_Xi, 1:s%nz)-xa_pre(i_Xi,1:s%nz)) < small)) exit

      end do

      deallocate(xa_pre)

   end subroutine do_instant_rain

   subroutine check_after_rain(s, precipitating_regions, xa_pre, ierr)
      type (star_info), pointer :: s
      type (immiscibility_region), allocatable, intent(in) :: precipitating_regions(:)
      real(kind=dp), allocatable, intent(in) :: xa_pre(:,:)
      integer, intent(out) :: ierr

      integer :: i, k, j
      integer :: klo, khi
      integer :: k_debug
      real(kind=dp) :: dY_max
      include 'formats'

      ! test that the sedimenting element's mass fraction is monotonically increasing and the other mass fractions are monotonically decreasing
      ! in the precipitating region
      do i = 1, size(precipitating_regions)
         klo = precipitating_regions(i)%klo
         khi = precipitating_regions(i)%khi
         !$OMP PARALLEL DO PRIVATE(k, j) SCHEDULE(guided)
         do k = klo, khi-1

         if (s%xa(i_Xi, k) - small > s%xa(i_Xi, k+1)) then
            ierr = 1
            write(*,*) 'check_after_rain: Xi is not monotonically increasing'
            write(*,*) 'klo, khi: ', klo, khi
            write(*,*) 'k, Xi(k), Xi(k+1): ', k, s%xa(i_Xi, k), s%xa(i_Xi, k+1)
            write(*,*) 'k, Y_pre(k), Y_pre(k+1): ', k, xa_pre(i_Xi, k), xa_pre(i_Xi, k+1)
         end if

         do j = 1, i_Xi - 1
            if (s%xa(j, k) + small < s%xa(j, k+1) .and. xa_pre(j, k) > xa_pre(j, k+1)) then
            ierr = 1
            write(*,*) 'check_after_rain: X_i is not monotonically decreasing for i = ', j
            write(*,*) 'klo, khi: ', klo, khi
            write(*,2) 'k, X_i(k): ', k, s% xa(j, k)
            write(*,2) 'k, X_i(k+1): ', k, s% xa(j, k+1)
            write(*,2) 'k, X_i_pre(k): ', k, xa_pre(j, k)
            write(*,2) 'k, X_i_pre(k+1): ', k, xa_pre(j, k+1)
            end if
         end do

         do j = i_Xi + 1, s%species
            if (s%xa(j, k) + small < s%xa(j, k+1) .and. xa_pre(j, k) > xa_pre(j, k+1)) then
            ierr = 1
            write(*,*) 'check_after_rain: X_i is not monotonically decreasing for i = ', j
            write(*,*) 'klo, khi: ', klo, khi
            write(*,2) 'k, X_i(k): ', k, s% xa(j, k)
            write(*,2) 'k, X_i(k+1): ', k, s% xa(j, k+1)
            write(*,2) 'k, X_i_pre(k): ', k, xa_pre(j, k)
            write(*,2) 'k, X_i_pre(k+1): ', k, xa_pre(j, k+1)
            end if
         end do
         end do
         !$OMP END PARALLEL DO

         ! check that Xi in the top reservoir either decreases or stays the same
         if (any(s%xa(i_Xi, precipitating_regions(i)%k_reservoir_lo:precipitating_regions(i)%k_reservoir_hi) - small > xa_pre(i_Xi, precipitating_regions(i)%k_reservoir_lo:precipitating_regions(i)%k_reservoir_hi))) then
            ierr = 1
            write(*,*) 'check_after_rain: Xi in the top reservoir is not decreasing'
            write(*,*) 'klo, khi: ', klo, khi
            write(*,*) 'k_reservoir_lo, k_reservoir_hi: ', precipitating_regions(i)%k_reservoir_lo, precipitating_regions(i)%k_reservoir_hi
            write(*,*) 'Xi(k_reservoir_lo), Xi(k_reservoir_hi): ', s%xa(i_Xi, precipitating_regions(i)%k_reservoir_lo), s%xa(i_Xi, precipitating_regions(i)%k_reservoir_hi)
            write(*,*) 'Y_pre(k_reservoir_lo), Y_pre(k_reservoir_hi): ', xa_pre(i_Xi, precipitating_regions(i)%k_reservoir_lo), xa_pre(i_Xi, precipitating_regions(i)%k_reservoir_hi)
         end if

         ! check that Xi in the top reservoir and precipitating region did not change by more than s% element_sedimentation_maximum_dXi_per_cell
         k_debug = maxloc(abs(s%xa(i_Xi, precipitating_regions(i)%k_reservoir_lo:precipitating_regions(i)%khi) - xa_pre(i_Xi, precipitating_regions(i)%k_reservoir_lo:precipitating_regions(i)%khi)), dim=1)
         dY_max = abs(s%xa(i_Xi, k_debug) - xa_pre(i_Xi, k_debug))
         if (dY_max > abs(s% element_sedimentation_maximum_dXi_per_cell) + small) then
            ierr = 1
            write(*,*) 'check_after_rain: Xi in the top reservoir and precipitating region changed by more than s% element_sedimentation_maximum_dXi_per_cell'
            write(*,*) 'k, abs(dY_max), abs(s% element_sedimentation_maximum_dXi_per_cell): ', k_debug, dY_max, abs(s% element_sedimentation_maximum_dXi_per_cell)
         end if

      end do
      if ( ierr /= 0 ) then
         call mesa_error(__FILE__,__LINE__)
      end if
   end subroutine check_after_rain

   subroutine monitor_mass_fractions(s, precipitating_region, ierr)
      type (star_info), pointer :: s
      type (immiscibility_region) :: precipitating_region
      integer, intent(out) :: ierr

      ! local
      real(kind=dp) :: M_deposite, M_rain
      integer :: k_reservoir_lo, k_reservoir_hi, k_deposite_lo, k_deposite_hi
      real(kind=dp), allocatable :: dXi(:)
      integer :: klo, khi


      ierr = 0 ! 0 means AOK
      k_reservoir_lo = precipitating_region%k_reservoir_lo
      k_reservoir_hi = precipitating_region%k_reservoir_hi
      k_deposite_lo = precipitating_region%k_deposite_lo
      k_deposite_hi = precipitating_region%k_deposite_hi
      klo = precipitating_region%klo
      khi = precipitating_region%khi

      ! allocate the mass fraction changes
      allocate(dXi(1:s%nz), source=0d0)
      dXi(1:s%nz) = precipitating_region%dXi(1:s%nz)

      M_rain = sum(s%dm(k_reservoir_lo:khi))
      M_deposite = sum(s%dm(k_deposite_lo:k_deposite_hi))
      write(*,*) '----------------------------------------'
      write(*,*) 'Helium rain convergence test:'
      write(*,*) '----------------------------------------'
      write(*,*) 'rain region: ', k_reservoir_lo, k_reservoir_hi, klo, khi
      write(*,*) 'integrated dY_rain/mstar: ', dot_product(s%dm(k_reservoir_lo:khi), dXi(k_reservoir_lo:khi))/M_rain
      write(*,*) ''
      write(*,*) 'depositing region: ', k_deposite_lo, k_deposite_hi
      write(*,*) 'do_instant_rain: integrated dXi/M_deposite: ', dot_product(s%dm(k_deposite_lo:k_deposite_hi), dXi(k_deposite_lo:k_deposite_hi))/M_deposite
      write(*,*) ''
      write(*,*) 'full region: ', 1, s%nz
      write(*,*) 'do_instant_rain: integrated dXi/mstar: ', dot_product(s%dm(1:s%nz), dXi(1:s%nz))/s%mstar
      write(*,*) '----------------------------------------'
      deallocate(dXi)
   end subroutine monitor_mass_fractions

   subroutine dummy_get_precipitating_bounds(s,  klo, khi, n_precipitating_regions, dXi)
      ! quick foo to isolate effect of phase diagram and element sedimentation
      type (star_info), pointer :: s
      integer, allocatable, intent(out) :: klo(:), khi(:)  ! bounds of the precipitating regions
      real(kind = dp), allocatable, intent(out) :: dXi(:)  ! difference between the miscible and the actual element mass fraction
      integer, intent(out) :: n_precipitating_regions  ! number of precipitating regions
      real(kind=dp)  :: Xi_crit, Xi_min, Xi_max, logt_min, logt_max, logt

      logt = log10(s% time_years)
      if (logt < 9d0) return

      allocate(dXi(1:s%nz), source=0d0)

      allocate(klo(1), khi(1))
      klo(1) = 400
      khi(1) = 420
      n_precipitating_regions = 1      

      Xi_min = 0.15d0
      Xi_max = 2.6d-1
      logt_min = 9d0
      logt_max = 9.6580113967d0 !9.6580113967d0
      
      Xi_crit = (logt_max*Xi_max - logt_min*Xi_min + logt*(-Xi_max + Xi_min))/(logt_max - logt_min)
      Xi_crit = max(Xi_crit, 0d0)
      
      dXi(klo(1):khi(1)) = Xi_crit - s%xa(s%net_iso(chem_get_iso_id('he4')), klo(1):khi(1))

   end subroutine dummy_get_precipitating_bounds

   subroutine enforce_Xi_gradient(Xi_miscible, klo, khi, ierr)
      real(kind=dp), allocatable, intent(inout) :: Xi_miscible(:)
      integer, allocatable, intent(in) :: klo(:), khi(:)
      integer, intent(out) :: ierr

      integer :: i, k

      logical :: debug = .false.

      if ( debug ) then
         write(*,*) 'enforce_Xi_gradient: before enforcing gradient: '
         do i = 1, size(klo)
            write(*,*) 'i, klo, khi: ', i, klo(i), khi(i)
            do k = klo(i), khi(i)
               write(*,*) 'k, Xi_miscible(k): ', k, Xi_miscible(k)
            end do
         end do
      end if

      ! make sure that the element mass fraction is monotonically increasing with depth
      do i = 1, size(klo)
         do k = khi(i), klo(i), -1
            if (Xi_miscible(k) < Xi_miscible(k-1)) then
               if ( debug ) write(*,*) 'enforce_Xi_gradient: enforcing gradient at k, Xi_miscible(k-1), Xi_miscible(k) :',&
                  k, Xi_miscible(k-1), Xi_miscible(k)
               Xi_miscible(k-1) = Xi_miscible(k)
            end if
         end do
      end do

      ! test if the element mass fraction is monotonically increasing
      do i = 1, size(klo)
         do k = klo(i), khi(i)-1
            if (Xi_miscible(k) > Xi_miscible(k+1)) then
               write(*,*) 'enforce_Xi_gradient: element mass fraction is not monotonically increasing'
               write(*,*) 'klo, khi: ', klo(i), khi(i)
               write(*,*) 'k, Xi_miscible(k), Xi_miscible(k+1): ', k, Xi_miscible(k), Xi_miscible(k+1)
               ierr = 1
               return
            end if
         end do
      end do

      if ( debug ) then
         write(*,*) 'enforce_Xi_gradient: after enforcing gradient: '
         do i = 1, size(klo)
            write(*,*) 'i, klo, khi: ', i, klo(i), khi(i)
            do k = klo(i), khi(i)
               write(*,*) 'k, Xi_miscible(k): ', k, Xi_miscible(k)
            end do
         end do
      end if


   end subroutine enforce_Xi_gradient

   subroutine get_raw_miscibility(s, is_miscible, Xi_miscible,  ierr)
      ! miscibility without any post-processing
      type(star_info), pointer :: s
      real(kind=dp), allocatable, intent(out) :: Xi_miscible(:)
      logical, allocatable, intent(out) :: is_miscible(:)  ! miscible regions
      
      ! local vars
      real(kind=dp), allocatable :: Xi(:)
      integer :: k
      integer :: op_err, ierr

      ierr = 0 ! 0 means AOK
      ! set up arrays
      allocate(Xi_miscible(s%nz), Xi(s%nz), is_miscible(s%nz))
      
      ! shorthand for the element mass fraction
      Xi = s% xa(i_Xi, 1:s%nz)

      !$OMP PARALLEL DO PRIVATE(k, op_err) SCHEDULE(guided)
      do k = 1, s%nz
         op_err = 0  ! Each thread gets its own local copy of ierr
         call Xi_is_miscible(s%lnPgas(k), s%lnT(k), Xi(k), is_miscible(k), Xi_miscible(k), immiscibility_phase_diagram, op_err)
         if (op_err /= 0) ierr = op_err
      end do
      !$OMP END PARALLEL DO
      if (ierr /= 0) return
      
   end subroutine get_raw_miscibility

   subroutine get_precipitating_bounds(s, klo, khi, n_precipitating_regions, dXi, ierr)
      type (star_info), pointer :: s
      integer, allocatable, intent(out) :: klo(:), khi(:)  ! bounds of the precipitating regions
      real(kind = dp), allocatable, intent(out) :: dXi(:)  ! difference between the miscible and the actual element mass fraction
      integer, intent(out) :: n_precipitating_regions  ! number of precipitating regions


      ! local
      integer, allocatable :: klo_dummy(:), khi_dummy(:)
      integer :: k, i, op_err
      real(kind=dp), allocatable :: Xi_miscible(:), Xi(:)
      logical, allocatable :: is_miscible(:)  ! miscible regions

      ! finding the precipitating regions
      logical :: in_region

      ! utils
      integer, intent(out) :: ierr
      integer :: k_debug

      include 'formats'

      call get_raw_miscibility(s, is_miscible, Xi_miscible,  ierr)

      ! set up arrays
      allocate(Xi(s%nz))
      
      Xi = s%xa(i_Xi, 1:s%nz)


      ! if is_miscible is .true. everywhere, we don't need to do anything
      if (all(is_miscible)) then
         deallocate(Xi_miscible, Xi, is_miscible)
         return
      end if


      ! find the precipitating regions
      allocate(klo_dummy(s%nz), khi_dummy(s%nz), source=-1)
      in_region = .false.
      do k = 1, s%nz
         if (.not. is_miscible(k) .and. .not. in_region) then
            klo_dummy(k) = k
            in_region = .true.
         else if (is_miscible(k) .and. in_region) then
            khi_dummy(k) = k - 1
            in_region = .false.
         end if
      end do

      ! if the last zone is not miscible, set the end of the region
      if (in_region) then
         khi_dummy(k) = s%nz
         in_region = .false.
      end if

      ! allocate the precipitating regions
      n_precipitating_regions = size(pack(klo_dummy, klo_dummy > 0))
      allocate(klo(n_precipitating_regions), khi(n_precipitating_regions), source=-1)
      klo = pack(klo_dummy(1:s%nz), klo_dummy(1:s%nz) > 0)
      khi = pack(khi_dummy(1:s%nz), khi_dummy(1:s%nz) > 0)

      if ( debug ) then
         write(*,*) 'get_precipitating_bounds: before joining regions: '
         do i = 1, n_precipitating_regions
            write(*,*) 'i, klo, khi: ', i, klo(i), khi(i)
            do k = klo(i), khi(i)
               write(*,*) 'k, is_miscible(k), Xi(k), Xi_miscible(k)', k, is_miscible(k), Xi(k), Xi_miscible(k)
            end do
         end do
      end if

      call join_small_precipitating_regions(s, klo, khi, n_precipitating_regions, Xi_miscible, is_miscible, ierr)

      ! check that size(klo) == size(khi)
      if (size(klo) /= size(khi)) then
         ierr = 1
         write(*,*) 'get_precipitating_bounds: size(klo) /= size(khi)'
         write(*,*) 'size(klo), size(khi): ', size(klo), size(khi)
         write(*,*) 'klo: ', klo
         write(*,*) 'khi: ', khi
         deallocate(Xi_miscible, Xi, is_miscible, klo_dummy, khi_dummy, klo, khi)
         return
      end if

      ! check that klo and khi are strictly increasing
      do i = 1, size(klo) - 1
         if (klo(i) >= klo(i+1)) then
            ierr = 1
            write(*,*) 'get_precipitating_bounds: klo is not strictly increasing'
            write(*,*) 'klo_clean(i), klo_clean(i+1): ', klo(i), klo(i+1)
            deallocate(Xi_miscible, Xi, is_miscible, klo_dummy, khi_dummy, klo, khi)
            return
         else if (khi(i) >= khi(i+1)) then
            ierr = 1
            write(*,*) 'get_precipitating_bounds: khi is not strictly increasing'
            write(*,*) 'khi_clean(i), khi_clean(i+1): ', khi(i), khi(i+1)
            deallocate(Xi_miscible, Xi, is_miscible, klo_dummy, khi_dummy, klo, khi)
            return
         end if
      end do

      if (damp_sedimentation) then
         do i = 1, size(klo)
            call reduce_Xi_difference(s, Xi_miscible, is_miscible, klo(i), khi(i), ierr)
            if (ierr /= 0) return
            ! call trim_region(s, Xi_miscible, is_miscible, klo(i), khi(i), ierr)
            ! if (ierr /= 0) return
         end do
      end if

      call enforce_Xi_gradient(Xi_miscible, klo, khi, ierr)
      if (ierr /= 0) return

      ! get the difference between the miscible and the actual element mass fraction in the precipitating region
      allocate(dXi(1:s%nz), source=0d0)
      dXi(1:s%nz) = merge(0d0, Xi_miscible(1:s%nz) - Xi(1:s%nz), is_miscible(1:s%nz))

      if (debug) then
         do k = 1, s%nz
            write(*,'(a,i12,l5,99(1pd26.16))') 'get_precipitating_bounds: (k, is_miscible, Xi, Xi_miscible, dXi): ', &
               k, is_miscible(k), Xi(k), Xi_miscible(k), dXi(k)
         end do
      end if

   end subroutine get_precipitating_bounds

   subroutine reduce_Xi_difference(s, Xi_miscible, is_miscible, klo, khi, ierr)
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(inout) :: Xi_miscible(:)
      real(kind=dp), allocatable :: Xi_old(:)
      logical, allocatable, intent(inout) :: is_miscible(:)
      integer, intent(inout) :: klo, khi  ! region in which we reduce the element difference
      integer, intent(out) :: ierr
      
      integer :: k, khi_new

      real(kind=dp) :: dXi_max

      allocate(Xi_old(klo:khi), source=0d0)
      Xi_old = s% xa_old(i_Xi, klo:khi)

      ! s% element_sedimentation_maximum_dXi_per_cell is positive for convience and to have a turn-off switch
      ! However, in reality it is negative, because material is raining out from a region
      ! To make this clear, we introduce the negative dXi_max here.
      dXi_max =  - s% element_sedimentation_maximum_dXi_per_cell
      
      ! the is_miscible check already makes sure that Xi_miscible is less than Xi; hence, we don't need absolute values
      Xi_miscible(klo:khi) = merge(Xi_old(klo:khi) + dXi_max, Xi_miscible(klo:khi), (.not. is_miscible(klo:khi)) .and. (Xi_miscible(klo:khi) < Xi_old(klo:khi) + dXi_max))
      
      deallocate(Xi_old)
   
   end subroutine reduce_Xi_difference

   subroutine trim_region(s, Xi_miscible, is_miscible, klo, khi, ierr)
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(inout) :: Xi_miscible(:)
      logical, allocatable, intent(inout) :: is_miscible(:)
      integer, intent(inout) :: klo, khi  ! region in which we reduce the element difference
      integer, intent(out) :: ierr
      
      ! local
      integer :: k, khi_new
      real(kind=dp), allocatable :: Xi(:)

      if (klo == khi) return  ! no trimming if the zone has size 1
      ! make sure that dXi is not zero for the bottom zone

      allocate(Xi(klo:khi), source=s%xa(i_Xi, klo:khi))
      
      
      khi_new = khi
      write(*,*) 'reduce_Xi_difference: klo+1, klo, khi, khi_new: ', klo+1, klo, khi_new  ! HHK
         do k = khi, klo+1, -1
            write(*,*) 'k, Xi_miscible(k), Xi(k): ', k, Xi_miscible(k), Xi(k)
            if (abs(Xi_miscible(k) - Xi(k)) > very_tiny) then
               exit
            else
               write(*,*) 'reduce_Xi_difference: Xi_miscible(k) is too small' ! HHK
               Xi_miscible(k) = -1d0
               khi_new = k - 1
               is_miscible(k) = .true.
            end if
         end do

      khi = khi_new
      write(*,*) 'reduce_Xi_difference: klo, khi, khi_new: ', klo, khi,  khi_new
      deallocate(Xi)
   end subroutine trim_region

   subroutine join_small_precipitating_regions(s, klo, khi, n_precipitating_regions, Xi_miscible, is_miscible, ierr)
      ! if regions are very small any only differ by a few zones, artifacts form.
      ! These regions are likely also not physical, because we interpolate the (uncertain) phase diagram linearly.

      type (star_info), pointer :: s
      integer, allocatable, intent(inout) :: klo(:), khi(:)
      integer, intent(inout) :: n_precipitating_regions
      logical, allocatable, intent(inout) :: is_miscible(:)
      real(kind=dp), allocatable, intent(inout) :: Xi_miscible(:)
      integer, intent(out) :: ierr

      ! local
      integer :: i, j, k

      logical, allocatable :: was_joined(:)
      integer, allocatable :: klo_temp(:), khi_temp(:)
      real(kind=dp) :: minimal_pressure_scale_height = 1d0
      real(kind=dp) :: dHp

      ! utils
      logical :: debug = .false.

      ierr = 0 ! 0 means AOK

      if (n_precipitating_regions < 2) return   ! nothing to do

      if (debug) then
         write(*,*) 'join_small_precipitating_regions: old regions: '
         do i = 1, n_precipitating_regions
            write(*,*) 'i, klo, khi: ', i, klo(i), khi(i)
         end do
      end if


      allocate(klo_temp(n_precipitating_regions), khi_temp(n_precipitating_regions), source=-1)
      allocate(was_joined(n_precipitating_regions), source=.false.)

      do i = 1, n_precipitating_regions-1

         if (was_joined(i)) cycle

         do j = n_precipitating_regions, i+1, -1

            dHp = pressure_scale_height_distance(s, klo(i), khi(i))

            if (dHp < minimal_pressure_scale_height) then
               if (debug) write(*,*) 'join_small_precipitating_regions: joining regions: ', i, j
               klo_temp(i) = klo(i)
               khi_temp(i) = khi(j)
               was_joined(i:j) = .true.
               exit
            else
               klo_temp(i) = klo(i)
               khi_temp(i) = khi(i)
            end if
         end do
      end do


      deallocate(klo, khi)
      n_precipitating_regions = size(pack(klo_temp, klo_temp > 0))
      allocate(klo(n_precipitating_regions), khi(n_precipitating_regions))
      klo = pack(klo_temp, klo_temp > 0)
      khi = pack(khi_temp, khi_temp > 0)

      if (debug) then
         write(*,*) 'join_small_precipitating_regions: new regions: '
         do i = 1, n_precipitating_regions
            write(*,*) 'i, klo, khi: ', i, klo(i), khi(i)
         end do
      end if

      ! adjust the miscible parameters accordingly
      do i = 1, size(khi)
         is_miscible(klo(i):khi(i)) = .false.
         Xi_miscible(klo(i):khi(i)) = merge(s%xa(s%net_iso(chem_get_iso_id('he4')),klo(i):khi(i)), Xi_miscible(klo(i):khi(i)), Xi_miscible(klo(i):khi(i)) < 0d0)
      end do

      deallocate(klo_temp, khi_temp, was_joined)

   end subroutine join_small_precipitating_regions

   subroutine determine_precipitating_regions(s, precipitating_regions, ierr)
      type (star_info), pointer :: s
      type (immiscibility_region), allocatable, intent(out) :: precipitating_regions(:)
      integer :: i, k, n_precipitating_regions

      ! precipitating region vars
      integer, allocatable :: klo(:), khi(:)
      real(kind = dp), allocatable :: dXi(:)

      ! utils
      integer, intent(out) :: ierr

      real(kind=dp) :: X, Y, Z
      include 'formats'

      ierr = 0 ! 0 means AOK


      call get_precipitating_bounds(s, klo, khi, n_precipitating_regions, dXi, ierr)
      ! call dummy_get_precipitating_bounds(s, klo, khi, n_precipitating_regions, dXi)
      if (ierr /= 0) return
      if (.not. allocated(klo)) then
         write(*,*) 'determine_precipitating_region: no precipitating region found.'
         return
      end if

      allocate(precipitating_regions(n_precipitating_regions))
      !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(guided)
      do i = 1, n_precipitating_regions

         precipitating_regions(i)%klo = klo(i)
         precipitating_regions(i)%khi = khi(i)

         ! at this point, we haven't determined the top reservoir and the depositing region yet. Hence, we consider the largest possible region, setting it to zero everywhere.
         allocate(precipitating_regions(i)%dXi(1:s%nz), source=0d0)

         precipitating_regions(i)%dXi(klo(i):khi(i)) = dXi(klo(i):khi(i))

      end do
      !$OMP END PARALLEL DO

      if ( debug ) then
         do i = 1, n_precipitating_regions
            write(*,*) 'determine_precipitating_region: (i, klo, khi): ', i, klo(i), khi(i)
            do k = klo(i), khi(i)
               call get_mass_fractions(s, s% xa(:, k), X, Y, Z)
               write(*,*) 'determine_precipitating_region: (k, X, Y, Z, dXi): ', &
                  i, k, X, Y, Z, precipitating_regions(i)%dXi(k)
            end do
         end do

      end if

      if ( verbose ) then
         write(*,11) 'determine_precipitating_region: number of precipitating regions = ', n_precipitating_regions
         write(*,11) 'determine_precipitating_region: precipitating regions: '
         do i = 1, n_precipitating_regions
            write(*,11) 'determine_precipitating_region (region, klo, khi): ', i, precipitating_regions(i)%klo, precipitating_regions(i)%khi
         end do
      end if

      deallocate(dXi, klo, khi)
   end subroutine determine_precipitating_regions


   subroutine determine_deposites(s, precipitating_regions, ierr)
      type (star_info), pointer :: s
      type (immiscibility_region), intent(inout) :: precipitating_regions(:)
      integer, intent(out) :: ierr

      integer :: k, klo_local, khi_max, k_deposite_hi, k_deposite_lo, khi
      integer :: i
      real(kind=dp) :: dH
      real(kind=dp), parameter :: minimal_pressure_scale_height = 1d0

      include 'formats'

      ierr = 0
      ! TODO: implement an algorithm in case the layer below the precipitating region is not convective
      ! I thought a bit about this. I could of course add the element to the radiative region, which would then trigger a convective instability.
      ! However, I think one can also make the argument that the rain droplets would fall until they reach a convective region, where they would be destroyed by convection.
      ! Both are too simplistic and I should solve the advection equation for the element. However, I think as a starting point its fine.
      ! I'll implement the first option for now; including the radiative region though.


      do i = 1, size(precipitating_regions)

         precipitating_regions(i)% k_deposite_hi = -1

         precipitating_regions(i)% k_deposite_lo = precipitating_regions(i)% khi + 1

         klo_local = precipitating_regions(i)% k_deposite_lo

         if (i == size(precipitating_regions)) then
            khi_max = s%nz
         else
            khi_max = precipitating_regions(i+1)%klo - 1
         end if

         if ( debug ) write(*,*) 'determine_deposites: (i, k_deposite_lo, k_deposite_hi, klo_local, khi_max): ', i, &
            precipitating_regions(i)%k_deposite_lo, precipitating_regions(i)%k_deposite_hi, klo_local, khi_max

         if (s%mixing_type(precipitating_regions(i)%k_deposite_lo) /= convective_mixing) then
            do k = (precipitating_regions(i)%k_deposite_lo+1), khi_max
               if (s%mixing_type(k) == convective_mixing) then
                  klo_local = k
                  exit
               end if
            end do

            if (klo_local == precipitating_regions(i)%k_deposite_lo) then
               if (verbose) write(*,*) 'determine_deposites: depositing region not found. Using the entire region.'
               precipitating_regions(i)%k_deposite_lo = precipitating_regions(i)%khi + 1
               precipitating_regions(i)%k_deposite_hi = khi_max
               cycle
            end if
         end if

         if ( debug ) write(*,*) 'determine_deposites: (i, klo_local, khi_max): ', i, klo_local, khi_max

         ! pressure scale height condition

         do k = (klo_local+1), khi_max

            ! to avoid very small deposites, we require that the element spreads at least a certain distance
            ! this is a bit of a hack, but it works for now

            dH = pressure_scale_height_distance(s, klo_local, k)

            if (s%mixing_type(k) /= convective_mixing .and. dH > minimal_pressure_scale_height) then
               precipitating_regions(i)%k_deposite_hi = k - 1
               exit
            end if
         end do

         ! if the region is not found, the entire region until the next precipitating region is used
         if (precipitating_regions(i)%k_deposite_hi == -1) then
            precipitating_regions(i)%k_deposite_hi = khi_max
         end if
      end do

      if ( verbose ) then
         do i = 1, size(precipitating_regions)
            write(*,*) 'determine_deposites (region, k_reservoir_lo, k_reservoir_hi, klo, khi, k_deposite_lo, k_deposite_hi): ', &
               i, &
               precipitating_regions(i)%k_reservoir_lo, precipitating_regions(i)%k_reservoir_hi,&
               precipitating_regions(i)%klo, precipitating_regions(i)%khi, &
               precipitating_regions(i)%k_deposite_lo, precipitating_regions(i)%k_deposite_hi
         end do
      endif

   end subroutine determine_deposites

   function pressure_scale_height_distance(s, k_top, k_bot) result(delta_Hp)
      type (star_info), pointer :: s
      integer, intent(in) :: k_top, k_bot
      real(kind=dp) :: delta_Hp

      real(kind=dp) :: top_r, bot_r, dr, Hp
      real(kind=dp) :: top_Hp, bot_Hp


      top_r = s%r(k_top)
      bot_r = s%r(k_bot)
      dr = top_r - bot_r

      top_Hp = s%scale_height(k_top)
      bot_Hp = s%scale_height(k_bot)
      Hp = (top_Hp + bot_Hp)/2

      delta_Hp = dr/Hp

   end function pressure_scale_height_distance

   subroutine determine_top_reservoirs(s,  precipitating_regions, ierr)
      type (star_info), pointer :: s
      type (immiscibility_region), intent(inout) :: precipitating_regions(:)
      integer, intent(out) :: ierr
      integer :: k, i, k_max
      include 'formats'

      ierr = 0

      do i = 1, size(precipitating_regions)

         precipitating_regions(i)% k_reservoir_lo = -1
         precipitating_regions(i)% k_reservoir_hi = precipitating_regions(i)%klo - 1

         if ( i>1 ) then
            k_max = precipitating_regions(i-1)%khi + 1
         else
            k_max = 1
         end if

         if ( debug ) write(*,11) 'determine_top_reservoirs: (i, k_reservoir_hi, k_max): ', i, precipitating_regions(i)%k_reservoir_hi, k_max

         ! find the top convective element reservoir
         do k = precipitating_regions(i)% k_reservoir_hi, k_max, -1
            if (s%mixing_type(k) /= convective_mixing) then

               if (debug) write(*,11) 'determine_top_reservoirs: (i, k, mixing_type): ', i, k, s%mixing_type(k)

               ! in the edge case where the first zone of the interation is not convective, we need to set the k_reservoir_hi to the first zone
               if ( k == precipitating_regions(i)% k_reservoir_hi) then
                  precipitating_regions(i)% k_reservoir_lo = k
               else
                  precipitating_regions(i)% k_reservoir_lo = k + 1
               end if

               exit
            else if (k == k_max) then
               precipitating_regions(i)% k_reservoir_lo = k
            end if
         end do

         ! if the region is not found, return an error
         if (precipitating_regions(i)%k_reservoir_lo == -1) then
            ierr = 1
            write(*,*) 'determine_top_reservoir: reservoir not found'
            return
         end if
      end do

      if ( verbose ) then
         do i = 1, size(precipitating_regions)
            write(*,*) 'determine_top_reservoir: precipitating region (i, klo, khi, k_reservoir_lo, k_reservoir_hi): ',&
               i, precipitating_regions(i)%klo, precipitating_regions(i)%khi,&
               precipitating_regions(i)%k_reservoir_lo, precipitating_regions(i)%k_reservoir_hi
         end do
      endif

   end subroutine determine_top_reservoirs

   subroutine check_precipitating_regions(precipitating_regions, ierr)
      type (immiscibility_region), allocatable, intent(in) :: precipitating_regions(:)
      integer, intent(out) :: ierr
      integer :: i, k
      logical :: strictly_increasing
      include 'formats'

      ierr = 0 ! 0 means AOK

      ! test that k_reservoir_lo <= k_reservoir_hi < klo <= khi < k_deposite_lo <= k_deposite_hi

      do i = 1, size(precipitating_regions)
         strictly_increasing = precipitating_regions(i)%k_reservoir_lo <= (precipitating_regions(i)%k_reservoir_hi) .and. &
            precipitating_regions(i)%k_reservoir_hi < precipitating_regions(i)%klo .and. &
            precipitating_regions(i)%klo <= precipitating_regions(i)%khi .and. &
            precipitating_regions(i)%khi < precipitating_regions(i)%k_deposite_lo .and. &
            precipitating_regions(i)%k_deposite_lo <= precipitating_regions(i)%k_deposite_hi
         if (.not. strictly_increasing) then
            ierr = 1
            write(*,*) 'check_precipitating_regions: region bounds are not strictly increasing'
            write(*,*) '(region, k_reservoir_lo, k_reservoir_hi, klo, khi, k_deposite_lo, k_deposite_hi): ', &
               i, precipitating_regions(i)%k_reservoir_lo, precipitating_regions(i)%k_reservoir_hi, &
               precipitating_regions(i)%klo, precipitating_regions(i)%khi, &
               precipitating_regions(i)%k_deposite_lo, precipitating_regions(i)%k_deposite_hi
            return
         end if
      end do

   end subroutine check_precipitating_regions

   subroutine set_top_reservoir_abundance_change(s, precipitating_region, ierr)
      type (star_info), pointer :: s
      type (immiscibility_region) :: precipitating_region
      integer, intent(out) :: ierr
      integer :: i
      real(kind=dp), allocatable :: Xi_crit(:)
      logical, allocatable :: is_miscible(:)
      integer :: klo, k_reservoir_lo, k_reservoir_hi
      ierr = 0 ! 0 means AOK

      k_reservoir_lo = precipitating_region%k_reservoir_lo
      k_reservoir_hi = precipitating_region%k_reservoir_hi

      ! if k_reservoir_lo == k_reservoir_hi, we don't need to do anything because the zone above the precipitating region is not convective
      if (k_reservoir_lo == k_reservoir_hi) return

      allocate(Xi_crit(k_reservoir_lo:k_reservoir_hi), source=s%xa(i_Xi, precipitating_region%klo) + precipitating_region%dXi(precipitating_region%klo))

      ! make sure cells change only by a maximum of s% element_sedimentation_maximum_dXi_per_cell if gentle mixing is on
      if ( damp_sedimentation ) then
         allocate(is_miscible(k_reservoir_lo:k_reservoir_hi), source=.false.)
         call reduce_Xi_difference(s, Xi_crit, is_miscible, k_reservoir_lo, k_reservoir_hi, ierr)
         deallocate(is_miscible)
         if (ierr /= 0) return
      end if

      precipitating_region%dXi(k_reservoir_lo:k_reservoir_hi) = merge(0d0, Xi_crit - s%xa(i_Xi, k_reservoir_lo:k_reservoir_hi), s%xa(i_Xi, k_reservoir_lo:k_reservoir_hi) < Xi_crit(k_reservoir_lo:k_reservoir_hi))
      deallocate(Xi_crit)

      if ( debug ) then
         write(*,*) 'set_top_reservoir_abundance_change: Xi(klo) = ', s%xa(i_Xi, precipitating_region%klo)
         write(*,*) 'set_top_reservoir_abundance_change: dXi(klo) = ', precipitating_region%dXi(precipitating_region%klo)
         write(*,*) 'set_top_reservoir_abundance_change: Xi_crit(klo) = ', s%xa(i_Xi, precipitating_region%klo) + precipitating_region%dXi(precipitating_region%klo)
         write(*,*) 'set_top_reservoir_abundance_change: Xi(k_reservoir_lo:k_reservoir_hi) = ', s%xa(i_Xi, k_reservoir_lo:k_reservoir_hi)
         write(*,*) 'set_top_reservoir_abundance_change: dXi(k_reservoir_lo:k_reservoir_hi) = ', precipitating_region%dXi(k_reservoir_lo:k_reservoir_hi)
      end if

   end subroutine set_top_reservoir_abundance_change

   subroutine set_deposit_region_abundance_change(s, precipitating_region, ierr)
      type (star_info), pointer :: s
      type (immiscibility_region) :: precipitating_region
      integer, intent(out) :: ierr

      ! local
      real(kind=dp) :: delta_M_Xi_rain, M_deposite
      integer :: k_reservoir_lo, k_deposite_lo, k_deposite_hi, khi
      logical :: debug = .false.
      ierr = 0 ! 0 means AOK

      k_reservoir_lo = precipitating_region%k_reservoir_lo
      khi = precipitating_region%khi
      k_deposite_lo = precipitating_region%k_deposite_lo
      k_deposite_hi = precipitating_region%k_deposite_hi

      if ( debug ) then
         write(*,*) 'set_deposit_region_abundance_change: k_reservoir_lo = ', k_reservoir_lo
         write(*,*) 'set_deposit_region_abundance_change: khi = ', khi
         write(*,*) 'set_deposit_region_abundance_change: k_deposite_lo = ', k_deposite_lo
         write(*,*) 'set_deposit_region_abundance_change: k_deposite_hi = ', k_deposite_hi
      endif

      ! get the rain mass of Xi
      delta_M_Xi_rain = - dot_product(s%dm(k_reservoir_lo:khi), precipitating_region%dXi(k_reservoir_lo:khi))

      ! mass abundance change in deeper layers
      M_deposite = sum(s%dm(k_deposite_lo:k_deposite_hi))
      precipitating_region%dXi(k_deposite_lo:k_deposite_hi) = delta_M_Xi_rain/M_deposite

   end subroutine set_deposit_region_abundance_change

   subroutine print_phase_diagram(pd)
      type(phase_diagram), intent(in) :: pd
      integer :: i
      include 'formats'

      write(*,*) '----------------------------------------'
      write(*,*) 'Phase diagram:', trim(pd%phase_diagram_filename)
      write(*,*) 'P, T, He-poor, He-rich'
      do i = 1, pd%n_points
         write(*,'(99(1pd26.16))') pd%P(i), pd%T(i), pd%He_poor(i), pd%He_rich(i)
      end do
      write(*,*) '----------------------------------------'
   end subroutine print_phase_diagram

   subroutine load_phase_diagram(s, pd, ierr)
      type(star_info), pointer :: s
      type(phase_diagram), intent(inout) :: pd
      integer, intent(out) :: ierr
      integer :: i, ios, unit
      character(len=256) :: line
      real(kind=dp), allocatable :: temp_data(:,:)

      ierr = 0

      ! Open the file
      open(newunit=unit, file=trim(pd%phase_diagram_filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 1
         write(*,*) 'load_phase_diagram: Unable to open file ', pd%phase_diagram_filename
         return
      end if

      ! Count the number of data points
      pd%n_points = -1
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         pd%n_points = pd%n_points + 1
      end do

      ! Rewind the file to read data
      rewind(unit)

      ! Allocate temporary storage
      allocate(temp_data(4, pd%n_points))

      ! Read the data
      read(unit, '(a)', iostat=ios) ! skip the header
      do i = 1, pd%n_points
         read(unit, *, iostat=ios) temp_data(1, i), temp_data(2, i), temp_data(3, i), temp_data(4, i)
         if (ios /= 0) then
            ierr = 1
            write(*,*) 'load_phase_diagram: Unable to read data from file ', pd%phase_diagram_filename
            close(unit)
            return
         end if
      end do

      ! Close the file
      close(unit)

      ! Allocate and assign data to the phase diagram
      allocate(pd%He_rich(pd%n_points), pd%He_poor(pd%n_points), pd%T(pd%n_points), pd%P(pd%n_points))
      pd%He_poor = temp_data(1, :)
      pd%He_rich = temp_data(2, :)
      pd%T = temp_data(3, :)
      pd%P = temp_data(4, :) * 1d12 ! convert from Mbars to cgs

      ! apply temperature shift
      pd%T = pd%T + s% element_sedimentation_phase_diagram_temperature_shift

      ! Deallocate temporary storage
      deallocate(temp_data)

      call unique_Ps_and_Ts(pd, ierr)
      if (ierr /= 0) return

      ! Check the phase diagram
      call check_phase_diagram(pd, ierr)
      if (ierr /= 0) return
   end subroutine load_phase_diagram

   subroutine unique_Ps_and_Ts(pd,  ierr)
      type(phase_diagram), intent(inout) :: pd
      integer, intent(out) :: ierr
      integer :: i, j
      integer, allocatable :: P_indices(:), T_indices(:)
      include 'formats'

      ierr = 0 ! 0 means AOK

      ! integer, allocatable :: num_Ps(:), num_Ts_at_P(:,:)
      ! real(kind=dp), allocatable :: P_unique(:), T_unique_at_P(:)

      ! at maximum, the number of unique P and T points is the same as the number of points
      ! -1 is used to indicate that the point is not unique
      allocate(P_indices(pd%n_points), source=-1)

      P_indices(1) = 1
      do i = 2, pd%n_points
         if ( abs(pd%P(i-1)-pd%P(i)) > tiny ) then
            P_indices(i) = i
         end if
      end do

      pd%num_unique_Ps = size(pack(P_indices, P_indices > 0))
      allocate(pd%unique_P_index(pd%num_unique_Ps), pd%unique_P(pd%num_unique_Ps), pd%num_Ts_at_unique_Ps(pd%num_unique_Ps))
      pd%unique_P_index = pack(P_indices, P_indices > 0)
      pd%unique_P = pd%P(pd%unique_P_index)

      pd%num_Ts_at_unique_Ps(:pd%num_unique_Ps-1) = pd%unique_P_index(2:) - pd%unique_P_index(:pd%num_unique_Ps-1)
      pd%num_Ts_at_unique_Ps(pd%num_unique_Ps) = pd%n_points - pd%unique_P_index(pd%num_unique_Ps) + 1

      if ( debug_phase_diagram ) then
         write(*,11) 'unique_Ps_and_Ts: num_Ps = ', pd%num_unique_Ps
         write(*,11) 'unique_Ps_and_Ts: unique_P_index = ', pd%unique_P_index
         write(*,1)  'unique_Ps_and_Ts: unique_P = ', pd%unique_P
         write(*,11) 'unique_Ps_and_Ts: num_Ts_at_unique_Ps = ', pd%num_Ts_at_unique_Ps
      end if

   end subroutine unique_Ps_and_Ts

   subroutine check_phase_diagram(pd, ierr)
      type(phase_diagram), intent(inout) :: pd
      integer, intent(out) :: ierr
      integer :: i

      do i = 1, pd%n_points-1

         ! check that no (P, T) point is repeated
         if (pd%P(i) == pd%P(i+1) .and. pd%T(i) == pd%T(i+1)) then
            ierr = 1
            write(*,*) 'check_phase_diagram: repeated (P, T) point in the phase diagram'
            return
         end if

         ! check that He-rich is greater than He-poor
         if (pd%He_rich(i) <= pd%He_poor(i)) then
            ierr = 1
            write(*,*) 'check_phase_diagram: He-rich is not greater than He-poor'
            return
         end if
      end do

      do i = 1, pd%num_unique_Ps-1

         ! check that unique_P_index is strictly monotonically increasing
         if (pd%unique_P_index(i) >= pd%unique_P_index(i+1)) then
            ierr = 1
            write(*,*) 'check_phase_diagram: unique_P_index is not strictly monotonically increasing'
            return
         end if

         ! check that unique_P is strictly monotonically increasing
         if (pd%unique_P(i) >= pd%unique_P(i+1)) then
            ierr = 1
            write(*,*) 'check_phase_diagram: unique_P is not strictly monotonically increasing'
            return
         end if

      end do


   end subroutine check_phase_diagram

   subroutine shutdown_phase_diagram(pd)
      type(phase_diagram), intent(inout) :: pd
      

      if (allocated(pd%He_rich)) deallocate(pd%He_rich)
      if (allocated(pd%He_poor)) deallocate(pd%He_poor)
      if (allocated(pd%T)) deallocate(pd%T)
      if (allocated(pd%P)) deallocate(pd%P)
      if (allocated(pd%unique_P)) deallocate(pd%unique_P)
      if (allocated(pd%unique_P_index)) deallocate(pd%unique_P_index)
      if (allocated(pd%num_Ts_at_unique_Ps)) deallocate(pd%num_Ts_at_unique_Ps)

   end subroutine shutdown_phase_diagram

   function linear_interpolation(x, x1, x2, y1, y2) result(y)
      real(kind=dp), intent(in) :: x, x1, x2, y1, y2
      real(kind=dp) :: y, fraction
      if (x2 == x1) then
         y = y1  ! oder Fehlerbehandlung
      else
         fraction = (x - x1) / (x2 - x1)
         y = y1 * (1d0 - fraction) + y2 * fraction
      end if
   end function linear_interpolation

   function linear_list_interpolation(x, x_list, y_list, outside_bounds, ierr) result(y)
      real(kind=dp), intent(in) :: x
      real(kind=dp), dimension(:), intent(in) :: x_list, y_list
      real(kind=dp) :: y
      logical, intent(out) :: outside_bounds
      integer, intent(out) :: ierr
      integer :: i
      ierr = 0 ! 0 means AOK

      ! check that x_list and y_list are the same size
      if (size(x_list) /= size(y_list)) then
         ierr = 1
         write(*,*) 'linear_list_interpolation: x_list and y_list are not the same size'
         return
      end if

      call bisection_locate(x, x_list, i, outside_bounds, ierr)
      if (ierr/=0) then
         write(*,*) 'linear_list_interpolation: bisection_locate failed'
         return
      end if

      ! linear interpolation
      y = linear_interpolation(x, x_list(i), x_list(i+1), y_list(i), y_list(i+1))
   end function linear_list_interpolation

   subroutine immiscible_Xi_fraction_interpolation(P0, T0, pd, Xi_min, Xi_max, ierr)
      real(kind=dp), intent(in) :: P0, T0
      type(phase_diagram), intent(in) :: pd
      real(kind=dp), intent(out) :: Xi_min, Xi_max
      real(kind=dp) :: P_lo, P_hi
      real(kind=dp) :: Xi_min_lo_interp, Xi_min_hi_interp, Xi_max_lo_interp, Xi_max_hi_interp
      real(kind=dp), allocatable :: T_lo(:), T_hi(:), Y_min_lo(:), Y_min_hi(:), Y_max_lo(:), Y_max_hi(:)
      integer, intent(out) :: ierr
      integer :: i, j
      logical :: outside_bounds
      include 'formats'
      ierr = 0 ! 0 means AOK

      ! -1d0 means out of bounds and thus by default miscible
      Xi_min = -1d0
      Xi_max = -1d0
      ! in principle, you should check this. However, when you paralalize the code, this
      ! will create a lot of overhead. So, I just assume that the phase diagram is correct for now
      ! call check_phase_diagram(pd, ierr)
      ! if (ierr /= 0) return
      outside_bounds = .false.
      ! find the pressure values that enclose the P0
      call bisection_locate(P0, pd%unique_P, j, outside_bounds, ierr)
      if (ierr /= 0) then
         write(*,*) 'immiscible_Xi_fraction_interpolation: bisection_locate failed'
         return
      end if

      ! if we are outside the phase diagram, everything is miscible by default
      if (outside_bounds) then
         if ( debug_phase_diagram ) write(*,*) 'immiscible_Xi_fraction_interpolation: P0 is outside the phase diagram. P0 = ', P0
         return
      end if

      P_lo = pd%unique_P(j)
      P_hi = pd%unique_P(j+1)

      ! now that we have the pressure values, find the respective temperature and Xi_min/Xi_max values
      allocate(T_lo(pd%num_Ts_at_unique_Ps(j)), T_hi(pd%num_Ts_at_unique_Ps(j+1)), &
         Y_min_lo(pd%num_Ts_at_unique_Ps(j)), Y_min_hi(pd%num_Ts_at_unique_Ps(j+1)), &
         Y_max_lo(pd%num_Ts_at_unique_Ps(j)), Y_max_hi(pd%num_Ts_at_unique_Ps(j+1)), source=-1d0)

      T_lo = pd%T(pd%unique_P_index(j):pd%unique_P_index(j+1)-1)
      Y_min_lo = pd%He_poor(pd%unique_P_index(j):pd%unique_P_index(j+1)-1)
      Y_max_lo = pd%He_rich(pd%unique_P_index(j):pd%unique_P_index(j+1)-1)

      T_hi = pd%T(pd%unique_P_index(j+1):pd%num_Ts_at_unique_Ps(j+1)+pd%unique_P_index(j+1)-1)  ! slightly cumbersome way because of the edge case where P(n-1) < P0 < P(n)
      Y_min_hi = pd%He_poor(pd%unique_P_index(j+1):pd%num_Ts_at_unique_Ps(j+1)+pd%unique_P_index(j+1)-1)
      Y_max_hi = pd%He_rich(pd%unique_P_index(j+1):pd%num_Ts_at_unique_Ps(j+1)+pd%unique_P_index(j+1)-1)

      ! check that the temperature values are strictly increasing
      do i = 1, size(T_lo)-1
         if (T_lo(i) >= T_lo(i+1)) then
            ierr = 1
            write(*,*) 'immiscible_Xi_fraction_interpolation: Tj is not strictly increasing'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
      end do

      do i = 1, size(T_hi)-1
         if (T_hi(i) >= T_hi(i+1)) then
            ierr = 1
            write(*,*) 'immiscible_Xi_fraction_interpolation: Tjp1 is not strictly increasing'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
      end do

      ! if T0 below the temperature range, return -1d0
      if (T0 < T_lo(1) .or. T0 < T_hi(1)) then
         if ( debug_phase_diagram ) write(*,*) 'immiscible_Xi_fraction_interpolation: T0 below the data range T0 = ', T0
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      ! if T0 is above the temperature range, extrapolate linearly (done in lineater_list_interpolation by default)

      ! Xi_min
      Xi_min_lo_interp = linear_list_interpolation(T0, T_lo, Y_min_lo, outside_bounds, ierr)
      if (ierr /= 0) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      ! if we are outside the bounds and don't extrapolate, return -1d0
      if (outside_bounds .and. .not. extrapolate) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      Xi_min_hi_interp = linear_list_interpolation(T0, T_hi, Y_min_hi, outside_bounds, ierr)
      if (ierr /= 0) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if
      ! if we are outside the bounds and don't extrapolate, return -1d0
      if (outside_bounds .and. .not. extrapolate) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      Xi_min = linear_interpolation(P0, P_lo, P_hi, Xi_min_lo_interp, Xi_min_hi_interp)

      if ( extrapolate ) then
         Xi_min = max(min(1d0, Xi_min), 0d0)  ! make sure we are in the range of [0, 1]
      else
         if (Xi_min < 0d0) then
            ierr = 1
            write(*,*) 'immiscible_Xi_fraction_interpolation: Y_min_0 < 0'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
   
         if (Xi_min > 1d0) then
            ierr = 1
            write(*,*) 'immiscible_Xi_fraction_interpolation: Y_min_0 > 1'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
      end if
      
      ! Xi_max
      Xi_max_lo_interp = linear_list_interpolation(T0, T_lo, Y_max_lo, outside_bounds, ierr)
      if (ierr /= 0) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if
      ! if we are outside the bounds and don't extrapolate, return -1d0
      if (outside_bounds .and. .not. extrapolate) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      Xi_max_hi_interp = linear_list_interpolation(T0, T_hi, Y_max_hi, outside_bounds, ierr)
      if (ierr /= 0) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if
      ! if we are outside the bounds and don't extrapolate, return -1d0
      if (outside_bounds .and. .not. extrapolate) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if
      
      ! linear interpolation for Xi_max
      Xi_max = linear_interpolation(P0, P_lo, P_hi, Xi_max_lo_interp, Xi_max_hi_interp)
      
      ! if extrapolation is allowed, the linear functions might go outside the range of [0, 1]
      ! in this case, we set the values to 0 and 1
      ! if extrapolation is not allowed, we need to check that the values are in the range of [0, 1]
      ! and return an error if they are not
      if (extrapolate) then
         Xi_max = min(max(0d0, Xi_max),1d0)
      else
         if (Xi_max < 0d0) then
            ierr = 1
            write(*,*) 'immiscible_Xi_fraction_interpolation: Y_max_0 < 0'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
         
         if (Xi_max > 1d0) then
            ierr = 1
            write(*,*) 'immiscible_Xi_fraction_interpolation: Y_max_0 > 1'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if

         if (Xi_min > Xi_max) then
            ierr = 1
            write(*,*) 'immiscible_Xi_fraction_interpolation: Xi_min > Xi_max'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
      end if

      if ( debug_phase_diagram ) then
         write(*,*) '----'
         write(*,*) 'immiscible_Xi_fraction_interpolation: P0(Mbar) = ', P0/1d12
         write(*,*) 'immiscible_Xi_fraction_interpolation: T0 = ', T0
         write(*,*) 'immiscible_Xi_fraction_interpolation: P_lo(Mbar) = ', P_lo/1d12
         write(*,*) 'immiscible_Xi_fraction_interpolation: P_hi(Mbar) = ', P_hi/1d12
         write(*,*) '----'
         write(*,*) 'immiscible_Xi_fraction_interpolation: T_lo = ', T_lo
         write(*,*) 'immiscible_Xi_fraction_interpolation: Y_min_lo = ', Y_min_lo
         write(*,*) 'immiscible_Xi_fraction_interpolation: Y_lo_interp = ', Xi_min_lo_interp
         write(*,*) '----'
         write(*,*) 'immiscible_Xi_fraction_interpolation: T_hi = ', T_hi
         write(*,*) 'immiscible_Xi_fraction_interpolation: Y_min_hi = ', Y_min_hi
         write(*,*) 'immiscible_Xi_fraction_interpolation: Y_hi_interp = ', Xi_min_hi_interp
         write(*,*) '----'
      end if

      deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)

   end subroutine immiscible_Xi_fraction_interpolation

   subroutine Xi_is_miscible(lnP, lnT, Xi, is_miscible, Xi_miscible, pd, ierr)
      ! gives the highest element mass fraction that is miscible before the immiscibility gap
      real(kind=dp), intent(in) :: lnP, lnT, Xi
      type(phase_diagram), intent(in) :: pd
      real(kind=dp), intent(out) :: Xi_miscible
      logical, intent(out) :: is_miscible
      real(kind=dp) :: P0, T0
      integer, intent(out) :: ierr
      real(kind=dp) :: Xi_max

      ierr = 0 ! 0 means AOK

      ! convert to lin (already in cgs)
      P0 = exp(lnP)
      T0 = exp(lnT)

      call immiscible_Xi_fraction_interpolation(P0, T0, pd, Xi_miscible, Xi_max, ierr)
      if (ierr /= 0) return

      ! check that Xi is in the miscible range
      if (Xi_max < Xi_miscible) then  ! outside of the miscible range
         is_miscible = .true.
      else if (Xi <= Xi_miscible .or. Xi >= Xi_max) then
         is_miscible = .true.
      else
         is_miscible = .false.
      end if

      if ( debug .and. .not. is_miscible) then
         write(*,*) 'Xi_is_miscible: P0 = ', P0/1d12
         write(*,*) 'Xi_is_miscible: T0 = ', T0
         write(*,*) 'Xi_is_miscible: Xi = ', Xi
         write(*,*) 'Xi_is_miscible: Xi_min = ', Xi_miscible
         write(*,*) 'Xi_is_miscible: Xi_max = ', Xi_max
      end if

   end subroutine Xi_is_miscible

   subroutine adjust_mass_fractions(s,  dXi, ierr)
      type (star_info), pointer :: s
      real(kind=dp), intent(in) :: dXi(:)
      integer, intent(out) :: ierr
      integer :: i, k, n_species
      real(kind=dp) :: f
      ierr = 0 ! 0 means AOK
      
      n_species = size(s%xa, 1)
      ! test that dXi is the same size as s%xa
      if (size(dXi) /= s%nz) then
         write(*,*) 'adjust_mass_fractions: dXi is not the same size as s%xa'
         ierr = 1
         return
      end if

      !$OMP PARALLEL DO PRIVATE(i, k, f) SCHEDULE(guided)
      do i = 1, n_species
         if (i == i_Xi) cycle
         do k = 1, s% nz
         f = 1d0-dXi(k)/(1d0 - s% xa(i_Xi, k))
         s% xa(i, k) = f * s% xa(i, k)
         end do
      end do
      !$OMP END PARALLEL DO

      s% xa(i_Xi, 1:s% nz) = s% xa(i_Xi, 1:s% nz) + dXi(1:s% nz)

      ! HHK: test that again we conserved mass

   end subroutine adjust_mass_fractions

   subroutine get_mass_fractions(s, xa, X, Y, Z)
      type (star_info), pointer :: s
      real(kind=dp), intent(in) :: xa(:)
      real(kind=dp), intent(out) :: X, Y, Z

      ! fluff that we don't need
      real(dp) :: abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
      
      ! get the composition info
      call basic_composition_info( &
         s% species, s% chem_id, xa, X, Y, Z, &
         abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
   end subroutine get_mass_fractions

   subroutine get_settling_velocity(s, v_settle, is_miscible, Xi_miscible, ierr)
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(out) :: v_settle(:)
      logical, allocatable, intent(in) :: is_miscible(:)
      real(kind=dp), allocatable, intent(in) :: Xi_miscible(:)
      integer, intent(out) :: ierr

      ! local variables
      real(kind=dp) :: v_settle_0
      ierr = 0 ! 0 means AOK

      allocate(v_settle(s%nz), source=0d0)

      v_settle_0 = s% element_sedimentation_velocity
      v_settle(1:s%nz) = merge(0d0, v_settle_0, is_miscible(1:s%nz))
      
   end subroutine get_settling_velocity

   subroutine tridag(a, b, c, f_n, f_np1)
      real(kind=dp), intent(in) :: a(:), b(:), c(:)
      real(kind=dp), intent(in) :: f_n(:)
      real(kind=dp), intent(inout) :: f_np1(:)
      
      ! local variables
      integer :: j, n
      real(kind=dp), allocatable :: gam(:)
      real(kind=dp) :: bet

      ! check if the input arrays are of the same size
      if (size(a) /= size(b) .or. size(b) /= size(c) .or. size(a) /= size(f_n) .or. size(f_n) /= size(f_np1)) then
         write(*,*) "tridag: Error: Input arrays must be of the same size."
         return
      end if
      if (size(a) == 0) then
         write(*,*) "tridag: Error: Input arrays are empty."
         return
      end if
      
      ! if (debug) write(*,*) "tridag: f_n(1) = ", f_n(1)
      ! if (debug) write(*,*) "tridag: f_np1(1) = ", f_np1(1)

      n = size(b)
      ! if (debug) write(*,*) 'tridag: n = ', n 
      
      allocate(gam(n))
      ! if (debug) write(*,*) 'tridag: gam allocated'

      if(b(1) == 0.0d0) then
         write(*,*) "tridag: Error: b[0] is zero."
         return
      end if

      ! forward elimination
      ! if (debug) write(*,*) "tridag: start forward elimination"
      bet = b(1)
      f_np1(1) = f_n(1) / bet
      

      ! if (debug) write(*,*) "tridag: j = 1 finished"

      do j = 2, n
         ! if (debug) write(*,*) "tridag: j = ", j
         gam(j) = c(j-1) / bet
         bet = b(j) - a(j) * gam(j)
         if (bet == 0.0d0) then
            write(*,*) "Error: Division by zero in forward elimination."
            return
         end if
         f_np1(j) = (f_n(j) - a(j) * f_np1(j-1)) / bet
      end do
      
      ! back substitution
      ! if (debug) write(*,*) "tridag: start backward elimination"
      do j = n-1, 1, -1
         ! if (debug) write(*,*) "tridag: j = ", j
         f_np1(j) = f_np1(j) - gam(j+1) * f_np1(j+1)
      end do

      ! deallocate temporary arrays
      deallocate(gam)

      
   end subroutine tridag

   subroutine setup_tridiagonal(adv_diff_state, a, b, c)
      type(advection_diffusion_state), intent(in) :: adv_diff_state
      real(kind=dp), allocatable, intent(out) :: a(:), b(:), c(:)
      integer :: i, j
      integer :: nz
   
      nz = adv_diff_state% nz
      
      ! Allocate the arrays
      allocate(a(nz-2), b(nz-2), c(nz-2))
      
      ! c(1) is the second index of 1, 2, …, N-1, N, i.e., 2

      do i = 1, nz-2
         ! we need to add one to i because a, b, and c, range from 2 to N-1 instead of 1 to N
         a(i) = -alpha_half_step_down(adv_diff_state, i+1)
         c(i) = -alpha_half_step_up(adv_diff_state, i+1) + beta_half_step_up(adv_diff_state, i+1)
         b(i) = 1d0 + alpha_half_step_down(adv_diff_state, i+1) + alpha_half_step_up(adv_diff_state, i+1) - beta_half_step_down(adv_diff_state, i+1)
      end do
      
      ! Set the boundary conditions
      b(1) = 1d0 + alpha_half_step_up(adv_diff_state, 2) - beta_half_step_down(adv_diff_state, 2)
      b(nz-2) = 1d0 + alpha_half_step_down(adv_diff_state, nz-1) - beta_half_step_down(adv_diff_state, nz-1) + beta_half_step_up(adv_diff_state, nz-1)

      if (debug_advection_diffusion_solver) then
         write(*,*) "--------------------------------"
         write(*,'(a,i12,a,1pd26.16)') "Tridiagonal matrix setup with N = ", nz, " dt = ", adv_diff_state% dt
         write(*,*) "Tridiagonal matrix coefficients:"
         do i = 1, nz-2
            write(*,'(i12,99(1pd26.16))') i, a(i), b(i), c(i)
         end do
         write(*,*) "Tridiagonal matrix setup successfully."
         write(*,*) "--------------------------------"
      end if
   end subroutine setup_tridiagonal

   real(kind=dp) function alpha_half_step_down(adv_diff_state, j)
      type(advection_diffusion_state), intent(in) :: adv_diff_state
      integer, intent(in) :: j
      
      real(kind=dp) :: dxCell, dxSurf, dt, D

      dxCell = (adv_diff_state% dr(j+1) + adv_diff_state% dr(j)) / 2.0d0
      dxSurf = adv_diff_state% dr(j)
      dt = adv_diff_state% dt
      D = (adv_diff_state% D(j-1) + adv_diff_state% D(j)) / 2.0d0
      
      ! Calculate the stability parameter
      alpha_half_step_down = alpha(dxCell, dxSurf, dt, D)
  
   end function alpha_half_step_down

   real(kind=dp) function alpha_half_step_up(adv_diff_state, j)
      type(advection_diffusion_state), intent(in) :: adv_diff_state
      integer, intent(in) :: j
      
      real(kind=dp) :: dxCell, dxSurf, dt, D

      dxCell = (adv_diff_state% dr(j+1) + adv_diff_state% dr(j)) / 2.0d0
      dxSurf = adv_diff_state% dr(j+1)
      dt = adv_diff_state% dt
      D = (adv_diff_state% D(j) + adv_diff_state% D(j+1)) / 2.0d0
      
      ! Calculate the stability parameter
      alpha_half_step_up = alpha(dxCell, dxSurf, dt, D)
  
   end function alpha_half_step_up

   function alpha(dxCell, dxSurf, dt, D)
      real(kind=dp), intent(in) :: dxCell, dxSurf, dt, D
      real(kind=dp) :: alpha
      
      ! Calculate the stability parameter
      alpha = D * dt / (dxCell * dxSurf)
      
   end function alpha

   real(kind=dp) function beta_half_step_down(adv_diff_state, j)
      type(advection_diffusion_state), intent(in) :: adv_diff_state
      integer, intent(in) :: j

      real(kind=dp) :: dxCell, dt, v
      dxCell = (adv_diff_state% dr(j+1) + adv_diff_state% dr(j)) / 2.0d0
      dt = adv_diff_state% dt
      v = (adv_diff_state% v(j-1) + adv_diff_state% v(j)) / 2.0d0
      ! Calculate the stability parameter
      beta_half_step_down = beta(dxCell, dt, v)

   end function beta_half_step_down

   real(kind=dp) function beta_half_step_up(adv_diff_state, j)
      type(advection_diffusion_state), intent(in) :: adv_diff_state
      integer, intent(in) :: j

      real(kind=dp) :: dxCell, dt, v
      dxCell = (adv_diff_state% dr(j+1) + adv_diff_state% dr(j)) / 2.0d0
      dt = adv_diff_state% dt
      v = (adv_diff_state% v(j) + adv_diff_state% v(j+1)) / 2.0d0
      ! Calculate the stability parameter
      beta_half_step_up = beta(dxCell, dt, v)
   
   end function beta_half_step_up

   real(kind=dp) function beta(dxCell, dt, v)
      real(kind=dp), intent(in) :: dxCell, dt, v
      
      ! Calculate the beta parameter
      beta = dt * v / dxCell

      ! throw an error if beta is greater than 1
      if (abs(beta) > 1.0d0) then
         write(*,*) "Error: beta is greater than 1.0, beta = ", beta
         stop
      end if
   
      
   end function beta

   function get_dr(s, k) result(dr)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      real(kind=dp) :: dr

      if (k == s% nz) then
         dr = s% r(k) - s% R_center
      else
         dr = s% r(k) - s% r(k+1)
      end if
      
   end function get_dr

   !* phase diagram utilities
   ! (1) locating the value in the phase diagram; use bisection because of robustness; alternatively, think about the hunt method

   subroutine check_locate_input(x, stop_interpolation)
      real(kind=dp), intent(in) :: x(:)   ! Array of values to search
      logical, intent(out) :: stop_interpolation  ! Flag to stop interpolation
      
      integer ::  n, i
      n = size(x)
      
      stop_interpolation = .false.
      if (n < 2) then
         write(*,*) "Error: Input array must have at least two elements."
         stop_interpolation = .true.
         return
      end if
      
      ! Check if the array is strictly increasing
      do i = 1, n-1
         if (x(i) >= x(i+1)) then
            write(*,*) "Error: Input array is not strictly increasing."
            stop_interpolation = .true.
            return
         end if
      end do
      
   end subroutine check_locate_input

   subroutine bisection_locate(x0, x, i, outside_bounds, ierr)
      real(kind=dp), intent(in) :: x0  ! Initial value to locate
      real(kind=dp), intent(in) :: x(:)   ! Array of values to search
      logical, intent(out) :: outside_bounds  ! Flag to indicate if x0 is outside the bounds of x
      integer, intent(out) :: i  ! Index of the located value
      integer, intent(out) :: ierr  ! Error code: 0 means AOK, 1 means error
      integer :: low, high, mid
      logical :: stop_interpolation

      ierr = 0  ! 0 means AOK
      low = 1
      high = size(x)
      
      call check_locate_input(x, stop_interpolation)
      if (stop_interpolation) then
         ierr = 1 
         return
      end if

      outside_bounds = .false.  ! Initialize the flag

      if (x0 < minval(x)) then
         i = 1  ! If x0 is less than the minimum value, return the first index
         outside_bounds = .true.
         return
      end if
      if (x0 > maxval(x)) then
         i = size(x) - 1  ! If x0 is greater than the maximum value, return one minus the last index
         outside_bounds = .true.
         return
      end if

      do while (high-low > 1)
         mid = (low + high) / 2
         if (x(mid) < x0) then
            low = mid
         else if (x(mid) > x0) then
            high = mid
         else
            i = mid
            return
         end if
      end do

      i = low
      
   end subroutine bisection_locate

   ! copied from phase_separation.F90
   ! subroutine update_model_ (s, kc_t, kc_b, do_brunt)

   !      use turb_info, only: set_mlt_vars
   !      use brunt, only: do_brunt_B
   !      use micro
        
   !      type(star_info), pointer :: s
   !      integer, intent(in)      :: kc_t
   !      integer, intent(in)      :: kc_b
   !      logical, intent(in)      :: do_brunt
        
   !      integer  :: ierr
   !      integer  :: kf_t
   !      integer  :: kf_b

   !      logical :: mask(s%nz)

   !      mask(:) = .true.

   !      ! Update the model to reflect changes in the abundances across
   !      ! cells kc_t:kc_b (the mask part of this call is unused, mask=true for all zones).
   !      ! Do updates at constant (P,T) rather than constant (rho,T).
   !      s%fix_Pgas = .true.
   !      call set_eos_with_mask(s, kc_t, kc_b, mask, ierr)
   !      if (ierr /= 0) then
   !         write(*,*) 'phase_separation: error from call to set_eos_with_mask'
   !         stop
   !      end if
   !      s%fix_Pgas = .false.
        
   !      ! Update opacities across cells kc_t:kc_b (this also sets rho_face
   !      ! and related quantities on faces kc_t:kc_b)        
   !      call set_micro_vars(s, kc_t, kc_b, &
   !           skip_eos=.TRUE., skip_net=.TRUE., skip_neu=.TRUE., skip_kap=.FALSE., ierr=ierr)
   !      if (ierr /= 0) then
   !         write(*,*) 'phase_separation: error from call to set_micro_vars'
   !         stop
   !      end if

   !      ! This is expensive, so only do it if we really need to.
   !      if(do_brunt) then
   !         ! Need to make sure we can set brunt for mix_outward calculation.
   !         if(.not. s% calculate_Brunt_B) then
   !            stop "phase separation requires s% calculate_Brunt_B = .true."
   !         end if
   !         call do_brunt_B(s, kc_t, kc_b, ierr) ! for unsmoothed_brunt_B
   !         if (ierr /= 0) then
   !            write(*,*) 'phase_separation: error from call to do_brunt_B'
   !            stop
   !         end if
   !      end if

   !      ! Finally update MLT for interior faces
        
   !      kf_t = kc_t
   !      kf_b = kc_b + 1
        
   !      call set_mlt_vars(s, kf_t+1, kf_b-1, ierr)
   !      if (ierr /= 0) then
   !         write(*,*) 'phase_separation: failed in call to set_mlt_vars during update_model_'
   !         stop
   !      endif
        
   !      ! Finish
        
   !      return
        
   !    end subroutine update_model_


end module element_sedimentation
