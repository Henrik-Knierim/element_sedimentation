! ***********************************************************************
!   Copyright (C) 2025  Henrik Knierim
! ***********************************************************************
module helium_rain

   use star_lib
   use const_def
   use chem_lib

   implicit none

   private

   public :: init_helium_rain
   public :: free_phase_diagram
   public :: do_helium_rain
   public :: determine_precipitating_regions
   public :: instant_rain
   public :: bottom_up_rain
   public :: advection_diffusion_rain
   public :: phase_diagram
   public :: immiscibility_region
   public :: helium_immiscibility_phase_diagram

   integer, parameter :: instant_rain = 0
   integer, parameter :: bottom_up_rain = 1
   integer, parameter :: advection_diffusion_rain = 2

   real(dp), parameter :: very_tiny = 1d-12
   real(dp), parameter :: tiny = 1d-10
   real(dp), parameter :: small = 1d-4
   real(dp), parameter :: helium_rain_threshold = 1d-4   ! ignore helium rain if the change in the helium mass fraction is smaller than this value

   real(kind=dp) :: gentle_mixing_maximum_helium_rain_dY ! maximum change allowed in the helium mass fraction
   integer :: i_Y = -1
   logical :: debug = .true.
   logical :: debug_advection_diffusion_solver = .false.
   logical :: debug_phase_diagram = .true.
   logical :: verbose = .true.

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
      integer :: k_reservoir_lo, k_reservoir_hi  ! the convective helium reservoir
      integer :: k_deposite_lo, k_deposite_hi  ! the region where the helium is deposited
      real(kind=dp), allocatable :: dY(:) ! change in the mass fractions
   end type immiscibility_region

   type advection_diffusion_state
      integer :: nz  ! number of zones
      real(kind=dp) :: dt  ! time step
      real(kind=dp), allocatable :: D(:), v(:)
      real(kind=dp), allocatable :: dr(:), Y(:), dY(:)  ! dY is the change in the helium mass fraction
   end type advection_diffusion_state

   type(phase_diagram) :: helium_immiscibility_phase_diagram
   logical :: phase_diagram_loaded = .false.


contains

   subroutine init_helium_rain(s, phase_diagram_filename, maximum_dY_change_per_step, ierr)
      type (star_info), pointer :: s
      character(len=*), intent(in) :: phase_diagram_filename
      real(kind=dp), intent(in) :: maximum_dY_change_per_step
      integer, intent(out) :: ierr

      include 'formats'

      ierr = 0 ! 0 means AOK

      gentle_mixing_maximum_helium_rain_dY = maximum_dY_change_per_step

      ! check if the phase diagram is already loaded
      ! if so, throw an error. This rountine should only be called once
      if (phase_diagram_loaded) then
         ierr = 1
         write(*,*) 'init_helium_rain: phase diagram already loaded'
         return
      end if

      ! load the phase diagram
      helium_immiscibility_phase_diagram%phase_diagram_filename = trim(phase_diagram_filename)
      call load_phase_diagram(helium_immiscibility_phase_diagram, ierr)
      phase_diagram_loaded = .true.
      if (ierr /= 0) return

      if (debug_phase_diagram) call print_phase_diagram(helium_immiscibility_phase_diagram)
      
   end subroutine init_helium_rain

   subroutine do_helium_rain(s,  method, ierr)

      type (star_info), pointer :: s
      integer :: method ! 0: instant rain, 1: constant velocity rain
      integer, intent(out) :: ierr

      ierr = 0 ! 0 means AOK

      if (i_Y == -1) i_Y = s% net_iso(chem_get_iso_id('he4'))

      select case (method)
       case (instant_rain)
         call do_instant_rain(s, ierr)
       case (bottom_up_rain)
         call do_bottom_up_rain(s, ierr)
       case (advection_diffusion_rain)
         call do_advection_diffusion_rain(s, ierr)
       case default
         ierr = 1
         write(*,*) 'do_helium_rain: unknown method for helium rain'
         return
      end select

   end subroutine do_helium_rain

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
               write(*,*) 'do_advection_diffusion_rain: (k, Y(k), dY(k)) = ', j, adv_diff_state%Y(j), adv_diff_state%dY(j)
            end do
         end if

         call adjust_mass_fractions(s, adv_diff_state% dY, ierr)

         ! check if the change in the helium mass fraction is below a threshold
         if (all(abs(adv_diff_state%dY) < helium_rain_threshold)) then
            if (verbose) write(*,*) 'do_advection_diffusion_rain: no significant change in helium mass fraction'
            if (verbose) write(*,*) 'do_advection_diffusion_rain: max(|dY|) = ', maxval(abs(adv_diff_state%dY))
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
      real(kind=dp), allocatable :: Y_miscible(:)  ! miscible helium mass fraction
      logical, allocatable :: is_miscible(:)  ! miscible regions
      real(kind=dp), allocatable :: v_settle(:)  ! settling velocity
      integer :: k

      ierr = 0 ! 0 means AOK

      ! function to populate a advection_diffusion_state type used for the advection-diffusion solver
      adv_diff_state % nz = s% nz
      allocate(adv_diff_state % D(1:s% nz), source=s% D_mix(1: s% nz))  ! diffusion coefficient
      allocate(adv_diff_state % Y(1:s% nz), source=s% xa(i_Y, 1: s% nz))  ! diffusion coefficient
      allocate(adv_diff_state % dY(1:s% nz), source=0d0)  ! diffusion coefficient

      ! radial step size
      allocate(adv_diff_state % dr(1:s% nz))  
      !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(guided)
      do k = 1, s % nz
         adv_diff_state % dr(k) = get_dr(s, k)
      end do
      !$OMP END PARALLEL DO

      call get_raw_miscibility(s, is_miscible, Y_miscible,  ierr)
      if (ierr /= 0) then
         write(*,*) 'init_advection_diffusion_state: error in get_raw_miscibility'
         return
      end if

      ! if everything is miscible, we can skip the advection-diffusion step
      if (all(is_miscible)) then
         if (debug_advection_diffusion_solver) write(*,*) 'init_advection_diffusion_state: all regions are miscible, skipping advection-diffusion step'
         return
      end if
      
      call get_settling_velocity(s, adv_diff_state % v, is_miscible, Y_miscible, ierr)
      if (ierr /= 0) then
         write(*,*) 'init_advection_diffusion_state: error in get_settling_velocity'
         return
      end if

      call get_time_step(adv_diff_state, ierr)
      if (ierr /= 0) then
         write(*,*) 'init_advection_diffusion_state: error in get_time_step'
         return
      end if

      if (debug_advection_diffusion_solver) then
         write(*,*) 'init_advection_diffusion_state: advection-diffusion state initialized'
         write(*,*) 'init_advection_diffusion_state: nz = ', adv_diff_state % nz
         write(*,*) 'init_advection_diffusion_state: dt = ', adv_diff_state % dt
         do k = 1, adv_diff_state % nz
            write(*,'(a,i12,99(1pd26.16))') 'init_advection_diffusion_state: (k, D(k), v(k), dr(k), Y(k)) = ', k, adv_diff_state % D(k), adv_diff_state % v(k), adv_diff_state % dr(k), adv_diff_state % Y(k)
         end do
      end if

      deallocate(Y_miscible, is_miscible)
      
   end subroutine init_advection_diffusion_state

   subroutine get_time_step(adv_diff_state, ierr)
      type(advection_diffusion_state), intent(inout) :: adv_diff_state
      integer, intent(out) :: ierr

      ! local
      real(kind=dp) :: tau_adv ! advection timescale
      real(kind=dp) :: tau_diff ! diffusion timescale
      real(kind=dp) :: tau_adv_div_tau_diff  ! ratio of advection to diffusion timescale
      real(kind=dp) :: L, v, D
      logical, allocatable :: mask(:)
      real(kind=dp) :: fraction_of_timescale = 1d7
      
      ierr = 0 ! 0 means AOK

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
      if (allocated(adv_diff_state % Y)) deallocate(adv_diff_state % Y)
      if (allocated(adv_diff_state % dY)) deallocate(adv_diff_state % dY)

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
      
      f_n(1:nz-2) = adv_diff_state% Y(2:nz-1)
      call tridag(a, b, c, f_n, f_np1)
      
      ! set change in helium abundance
      adv_diff_state% dY(1) = f_np1(1) - f_n(1) ! Neumann BC at r=0
      adv_diff_state% dY(2:nz-1) = f_np1(1:nz-2) - f_n(1:nz-2) ! interior points
      adv_diff_state% dY(nz) = f_np1(nz-2) - f_n(nz-2) ! Neumann BC at r=L

      if ( debug_advection_diffusion_solver ) then
         write(*,'(a16,a16,a16,a16)') 'k','f_n(k)', 'f_np1(k)', 'dY(k)'
         do k = 1, nz-2
            write(*,'(i12,99(1pd26.16))') k+1, f_n(k), f_np1(k), adv_diff_state% dY(k+1)
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



         ! determine the upper convective helium reservoir
         call determine_top_reservoirs(s, precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_bottom_up_rain: error in determine_top_reservoirs'
            return
         endif

         ! determine deposit region
         call determine_helium_deposites(s, precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_bottom_up_rain: error in determine_helium_deposites'
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
            ! Set all changes in dY except the lowest cell equal to zero
            ! no need to adjust the top reservoir (yet)
            if (klo /= khi) then
               write(*,*) 'precipitating_regions(i)%dY(klo:khi): ', precipitating_regions(i)%dY(klo:khi)
               k_mix = khi
               do k = khi, klo, -1
                  if (abs(precipitating_regions(i)%dY(k)) < helium_rain_threshold) then
                     precipitating_regions(i)%dY(k) = 0d0
                     k_mix = k - 1
                  else
                     exit
                  end if
               end do
               precipitating_regions(i)%dY(1:k_mix-1) = 0d0
               write(*,*) 'precipitating_regions(i)%dY(klo:khi): ', precipitating_regions(i)%dY(klo:khi)
            end if

            call set_top_reservoir_abundance_change(s,  precipitating_regions(i), ierr)
            if (ierr /= 0) then
               write(*,*) 'do_bottom_up_rain: error in set_top_reservoir_abundance_change'
               return
            endif

            call set_Helium_deposit_region_abundance_change(s,  precipitating_regions(i), ierr)
            if (ierr /= 0) then
               write(*,*) 'do_bottom_up_rain: error in set_Helium_deposit_region_abundance_change'
               return
            endif

            if ( debug ) then
               do k_debug = 1, s%nz
                  call get_mass_fractions(s, s% xa(:,k_debug), X, Y, Z)
                  write(*,2) 'do_bottom_up_rain: (k, X, Y, Z, dY): ', &
                     k_debug, X, Y, Z, &
                     precipitating_regions(i)%dY(k_debug)
               end do
            end if

            call adjust_mass_fractions(s, precipitating_regions(i)% dY, ierr)
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
            k_debug = maxloc(abs(s%xa(i_Y, 1:s%nz)-xa_pre(i_Y, 1:s%nz)), dim=1)
            if (verbose) write(*,*) 'maximum difference: (k, val)', k_debug,  abs(s%xa(i_Y, k_debug)-xa_pre(i_Y, k_debug))
         end if
         
         if (all(abs(s%xa(i_Y, 1:s%nz)-xa_pre(i_Y,1:s%nz)) < small)) exit

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

         ! determine the upper convective helium reservoir
         call determine_top_reservoirs(s, precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_instant_rain: error in determine_top_reservoirs'
            return
         endif

         ! determine deposit region
         call determine_helium_deposites(s, precipitating_regions, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_instant_rain: error in determine_helium_deposites'
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

            call set_Helium_deposit_region_abundance_change(s,  precipitating_regions(i), ierr)
            if (ierr /= 0) then
               write(*,*) 'do_instant_rain: error in set_Helium_deposit_region_abundance_change'
               return
            endif

            if ( debug ) then
               do k_debug = 1, s%nz
                  call get_mass_fractions(s, s% xa(:,k_debug), X, Y, Z)
                  write(*,2) 'do_instant_rain: (k, X, Y, Z, dY): ', &
                     k_debug, X, Y, Z, &
                     precipitating_regions(i)%dY(k_debug)
               end do
            end if

            call adjust_mass_fractions(s, precipitating_regions(i)% dY, ierr)
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
            k_debug = maxloc(abs(s%xa(i_Y, 1:s%nz)-xa_pre(i_Y, 1:s%nz)), dim=1)
            if (verbose) write(*,*) 'maximum difference: (k, val)', k_debug,  abs(s%xa(i_Y, k_debug)-xa_pre(i_Y, k_debug))
         end if
         if (all(abs(s%xa(i_Y, 1:s%nz)-xa_pre(i_Y,1:s%nz)) < small)) exit

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

      ! test that the helium mass fraction is monotonically increasing and the X and Z mass fractions are monotonically decreasing
      ! in the precipitating region
      do i = 1, size(precipitating_regions)
         klo = precipitating_regions(i)%klo
         khi = precipitating_regions(i)%khi
         !$OMP PARALLEL DO PRIVATE(k, j) SCHEDULE(guided)
         do k = klo, khi-1

         if (s%xa(i_Y, k) - small > s%xa(i_Y, k+1)) then
            ierr = 1
            write(*,*) 'check_after_rain: Y is not monotonically increasing'
            write(*,*) 'klo, khi: ', klo, khi
            write(*,*) 'k, Y(k), Y(k+1): ', k, s%xa(i_Y, k), s%xa(i_Y, k+1)
            write(*,*) 'k, Y_pre(k), Y_pre(k+1): ', k, xa_pre(i_Y, k), xa_pre(i_Y, k+1)
         end if

         do j = 1, i_Y - 1
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

         do j = i_Y + 1, s%species
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

         ! check that Y in the top reservoir either decreases or stays the same
         if (any(s%xa(i_Y, precipitating_regions(i)%k_reservoir_lo:precipitating_regions(i)%k_reservoir_hi) - small > xa_pre(i_Y, precipitating_regions(i)%k_reservoir_lo:precipitating_regions(i)%k_reservoir_hi))) then
            ierr = 1
            write(*,*) 'check_after_rain: Y in the top reservoir is not decreasing'
            write(*,*) 'klo, khi: ', klo, khi
            write(*,*) 'k_reservoir_lo, k_reservoir_hi: ', precipitating_regions(i)%k_reservoir_lo, precipitating_regions(i)%k_reservoir_hi
            write(*,*) 'Y(k_reservoir_lo), Y(k_reservoir_hi): ', s%xa(i_Y, precipitating_regions(i)%k_reservoir_lo), s%xa(i_Y, precipitating_regions(i)%k_reservoir_hi)
            write(*,*) 'Y_pre(k_reservoir_lo), Y_pre(k_reservoir_hi): ', xa_pre(i_Y, precipitating_regions(i)%k_reservoir_lo), xa_pre(i_Y, precipitating_regions(i)%k_reservoir_hi)
         end if

         ! check that Y in the top reservoir and precipitating region did not change by more than gentle_mixing_maximum_helium_rain_dY
         k_debug = maxloc(abs(s%xa(i_Y, precipitating_regions(i)%k_reservoir_lo:precipitating_regions(i)%khi) - xa_pre(i_Y, precipitating_regions(i)%k_reservoir_lo:precipitating_regions(i)%khi)), dim=1)
         dY_max = abs(s%xa(i_Y, k_debug) - xa_pre(i_Y, k_debug))
         if (dY_max > abs(gentle_mixing_maximum_helium_rain_dY) + small) then
            ierr = 1
            write(*,*) 'check_after_rain: Y in the top reservoir and precipitating region changed by more than gentle_mixing_maximum_helium_rain_dY'
            write(*,*) 'k, abs(dY_max), abs(gentle_mixing_maximum_helium_rain_dY): ', k_debug, dY_max, abs(gentle_mixing_maximum_helium_rain_dY)
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
      real(kind=dp), allocatable :: dY(:)
      integer :: klo, khi


      ierr = 0 ! 0 means AOK
      k_reservoir_lo = precipitating_region%k_reservoir_lo
      k_reservoir_hi = precipitating_region%k_reservoir_hi
      k_deposite_lo = precipitating_region%k_deposite_lo
      k_deposite_hi = precipitating_region%k_deposite_hi
      klo = precipitating_region%klo
      khi = precipitating_region%khi

      ! allocate the mass fraction changes
      allocate(dY(1:s%nz), source=0d0)
      dY(1:s%nz) = precipitating_region%dY(1:s%nz)

      M_rain = sum(s%dm(k_reservoir_lo:khi))
      M_deposite = sum(s%dm(k_deposite_lo:k_deposite_hi))
      write(*,*) '----------------------------------------'
      write(*,*) 'Helium rain convergence test:'
      write(*,*) '----------------------------------------'
      write(*,*) 'rain region: ', k_reservoir_lo, k_reservoir_hi, klo, khi
      write(*,*) 'integrated dY_rain/mstar: ', dot_product(s%dm(k_reservoir_lo:khi), dY(k_reservoir_lo:khi))/M_rain
      write(*,*) ''
      write(*,*) 'depositing region: ', k_deposite_lo, k_deposite_hi
      write(*,*) 'do_instant_rain: integrated dY/M_deposite: ', dot_product(s%dm(k_deposite_lo:k_deposite_hi), dY(k_deposite_lo:k_deposite_hi))/M_deposite
      write(*,*) ''
      write(*,*) 'full region: ', 1, s%nz
      write(*,*) 'do_instant_rain: integrated dY/mstar: ', dot_product(s%dm(1:s%nz), dY(1:s%nz))/s%mstar
      write(*,*) '----------------------------------------'
      deallocate(dY)
   end subroutine monitor_mass_fractions

   subroutine dummy_get_precipitating_bounds(s,  klo, khi, n_precipitating_regions, dY)
      ! quick foo to isolate effect of phase diagram and helium rain
      type (star_info), pointer :: s
      integer, allocatable, intent(out) :: klo(:), khi(:)  ! bounds of the precipitating regions
      real(kind = dp), allocatable, intent(out) :: dY(:)  ! difference between the miscible and the actual helium mass fraction
      integer, intent(out) :: n_precipitating_regions  ! number of precipitating regions
      real(kind=dp)  :: Y_crit, Y_min, Y_max, logt_min, logt_max, logt

      logt = log10(s% time_years)
      if (logt < 9d0) return

      allocate(dY(1:s%nz), source=0d0)

      allocate(klo(1), khi(1))
      klo(1) = 400
      khi(1) = 420
      n_precipitating_regions = 1      

      Y_min = 0.15d0
      Y_max = 2.6d-1
      logt_min = 9d0
      logt_max = 9.6580113967d0 !9.6580113967d0
      
      Y_crit = (logt_max*Y_max - logt_min*Y_min + logt*(-Y_max + Y_min))/(logt_max - logt_min)
      Y_crit = max(Y_crit, 0d0)
      
      dY(klo(1):khi(1)) = Y_crit - s%xa(s%net_iso(chem_get_iso_id('he4')), klo(1):khi(1))

   end subroutine dummy_get_precipitating_bounds

   subroutine enforce_helium_gradient(Y_miscible, klo, khi, ierr)
      real(kind=dp), allocatable, intent(inout) :: Y_miscible(:)
      integer, allocatable, intent(in) :: klo(:), khi(:)
      integer, intent(out) :: ierr

      integer :: i, k

      logical :: debug = .false.

      if ( debug ) then
         write(*,*) 'enforce_helium_gradient: before enforcing gradient: '
         do i = 1, size(klo)
            write(*,*) 'i, klo, khi: ', i, klo(i), khi(i)
            do k = klo(i), khi(i)
               write(*,*) 'k, Y_miscible(k): ', k, Y_miscible(k)
            end do
         end do
      end if

      ! make sure that the helium mass fraction is monotonically increasing with depth
      do i = 1, size(klo)
         do k = khi(i), klo(i), -1
            if (Y_miscible(k) < Y_miscible(k-1)) then
               if ( debug ) write(*,*) 'enforce_helium_gradient: enforcing gradient at k, Y_miscible(k-1), Y_miscible(k) :',&
                  k, Y_miscible(k-1), Y_miscible(k)
               Y_miscible(k-1) = Y_miscible(k)
            end if
         end do
      end do

      ! test if the helium mass fraction is monotonically increasing
      do i = 1, size(klo)
         do k = klo(i), khi(i)-1
            if (Y_miscible(k) > Y_miscible(k+1)) then
               write(*,*) 'enforce_helium_gradient: helium mass fraction is not monotonically increasing'
               write(*,*) 'klo, khi: ', klo(i), khi(i)
               write(*,*) 'k, Y_miscible(k), Y_miscible(k+1): ', k, Y_miscible(k), Y_miscible(k+1)
               ierr = 1
               return
            end if
         end do
      end do

      if ( debug ) then
         write(*,*) 'enforce_helium_gradient: after enforcing gradient: '
         do i = 1, size(klo)
            write(*,*) 'i, klo, khi: ', i, klo(i), khi(i)
            do k = klo(i), khi(i)
               write(*,*) 'k, Y_miscible(k): ', k, Y_miscible(k)
            end do
         end do
      end if


   end subroutine enforce_helium_gradient

   subroutine get_raw_miscibility(s, is_miscible, Y_miscible,  ierr)
      ! miscibility without any post-processing
      type(star_info), pointer :: s
      real(kind=dp), allocatable, intent(out) :: Y_miscible(:)
      logical, allocatable, intent(out) :: is_miscible(:)  ! miscible regions
      
      ! local vars
      real(kind=dp), allocatable :: Y(:)
      integer :: k
      integer :: op_err, ierr

      ierr = 0 ! 0 means AOK
      ! set up arrays
      allocate(Y_miscible(s%nz), Y(s%nz), is_miscible(s%nz))

      ! if i_Y is not set (=-1), call chem_get_iso_id('he4') to get the index
      if (i_Y == -1) i_Y = s% net_iso(chem_get_iso_id('he4'))
      
      ! shorthand for the helium mass fraction
      Y = s% xa(i_Y, 1:s%nz)

      !$OMP PARALLEL DO PRIVATE(k, op_err) SCHEDULE(guided)
      do k = 1, s%nz
         op_err = 0  ! Each thread gets its own local copy of ierr
         call Y_is_miscible(s%lnPgas(k), s%lnT(k), Y(k), is_miscible(k), Y_miscible(k), helium_immiscibility_phase_diagram, op_err)
         if (op_err /= 0) ierr = op_err
      end do
      !$OMP END PARALLEL DO
      if (ierr /= 0) return
      
   end subroutine get_raw_miscibility

   subroutine get_precipitating_bounds(s, klo, khi, n_precipitating_regions, dY, ierr)
      type (star_info), pointer :: s
      integer, allocatable, intent(out) :: klo(:), khi(:)  ! bounds of the precipitating regions
      real(kind = dp), allocatable, intent(out) :: dY(:)  ! difference between the miscible and the actual helium mass fraction
      integer, intent(out) :: n_precipitating_regions  ! number of precipitating regions


      ! local
      integer, allocatable :: klo_dummy(:), khi_dummy(:)
      integer :: k, i, op_err
      real(kind=dp), allocatable :: Y_miscible(:), Y(:)
      logical, allocatable :: is_miscible(:)  ! miscible regions

      ! finding the precipitating regions
      logical :: in_region

      ! utils
      integer, intent(out) :: ierr
      integer :: k_debug

      include 'formats'

      call get_raw_miscibility(s, is_miscible, Y_miscible,  ierr)

      ! set up arrays
      allocate(Y(s%nz))
      
      Y = s%xa(i_Y, 1:s%nz)


      ! if is_miscible is .true. everywhere, we don't need to do anything
      if (all(is_miscible)) then
         deallocate(Y_miscible, Y, is_miscible)
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
               write(*,*) 'k, is_miscible(k), Y(k), Y_miscible(k)', k, is_miscible(k), Y(k), Y_miscible(k)
            end do
         end do
      end if

      call join_small_precipitating_regions(s, klo, khi, n_precipitating_regions, Y_miscible, is_miscible, ierr)

      ! check that size(klo) == size(khi)
      if (size(klo) /= size(khi)) then
         ierr = 1
         write(*,*) 'get_precipitating_bounds: size(klo) /= size(khi)'
         write(*,*) 'size(klo), size(khi): ', size(klo), size(khi)
         write(*,*) 'klo: ', klo
         write(*,*) 'khi: ', khi
         deallocate(Y_miscible, Y, is_miscible, klo_dummy, khi_dummy, klo, khi)
         return
      end if

      ! check that klo and khi are strictly increasing
      do i = 1, size(klo) - 1
         if (klo(i) >= klo(i+1)) then
            ierr = 1
            write(*,*) 'get_precipitating_bounds: klo is not strictly increasing'
            write(*,*) 'klo_clean(i), klo_clean(i+1): ', klo(i), klo(i+1)
            deallocate(Y_miscible, Y, is_miscible, klo_dummy, khi_dummy, klo, khi)
            return
         else if (khi(i) >= khi(i+1)) then
            ierr = 1
            write(*,*) 'get_precipitating_bounds: khi is not strictly increasing'
            write(*,*) 'khi_clean(i), khi_clean(i+1): ', khi(i), khi(i+1)
            deallocate(Y_miscible, Y, is_miscible, klo_dummy, khi_dummy, klo, khi)
            return
         end if
      end do

      if (s% do_gentle_mixing) then
         do i = 1, size(klo)
            call reduce_helium_difference(s, Y_miscible, is_miscible, klo(i), khi(i), ierr)
            if (ierr /= 0) return
            ! call trim_region(s, Y_miscible, is_miscible, klo(i), khi(i), ierr)
            ! if (ierr /= 0) return
         end do
      end if

      call enforce_helium_gradient(Y_miscible, klo, khi, ierr)
      if (ierr /= 0) return

      ! get the difference between the miscible and the actual helium mass fraction in the precipitating region
      allocate(dY(1:s%nz), source=0d0)
      dY(1:s%nz) = merge(0d0, Y_miscible(1:s%nz) - Y(1:s%nz), is_miscible(1:s%nz))

      if (debug) then
         do k = 1, s%nz
            write(*,'(a,i12,l5,99(1pd26.16))') 'get_precipitating_bounds: (k, is_miscible, Y, Y_miscible, dY): ', &
               k, is_miscible(k), Y(k), Y_miscible(k), dY(k)
         end do
      end if

   end subroutine get_precipitating_bounds

   subroutine reduce_helium_difference(s, Y_miscible, is_miscible, klo, khi, ierr)
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(inout) :: Y_miscible(:)
      real(kind=dp), allocatable :: Y_old(:)
      logical, allocatable, intent(inout) :: is_miscible(:)
      integer, intent(inout) :: klo, khi  ! region in which we reduce the helium difference
      integer, intent(out) :: ierr
      
      integer :: k, khi_new

      allocate(Y_old(klo:khi), source=0d0)
      Y_old = s% xa_old(i_Y, klo:khi)
      
      ! the is_miscible check already makes sure that Y_miscible is less than Y; hence, we don't need absolute values
      Y_miscible(klo:khi) = merge(Y_old(klo:khi) + gentle_mixing_maximum_helium_rain_dY, Y_miscible(klo:khi), (.not. is_miscible(klo:khi)) .and. (Y_miscible(klo:khi) < Y_old(klo:khi) + gentle_mixing_maximum_helium_rain_dY))
      
      deallocate(Y_old)
   
   end subroutine reduce_helium_difference

   subroutine trim_region(s, Y_miscible, is_miscible, klo, khi, ierr)
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(inout) :: Y_miscible(:)
      logical, allocatable, intent(inout) :: is_miscible(:)
      integer, intent(inout) :: klo, khi  ! region in which we reduce the helium difference
      integer, intent(out) :: ierr
      
      ! local
      integer :: k, khi_new
      real(kind=dp), allocatable :: Y(:)

      if (klo == khi) return  ! no trimming if the zone has size 1
      ! make sure that dY is not zero for the bottom zone

      allocate(Y(klo:khi), source=s%xa(i_Y, klo:khi))
      
      
      khi_new = khi
      write(*,*) 'reduce_helium_difference: klo+1, klo, khi, khi_new: ', klo+1, klo, khi_new  ! HHK
         do k = khi, klo+1, -1
            write(*,*) 'k, Y_miscible(k), Y(k): ', k, Y_miscible(k), Y(k)
            if (abs(Y_miscible(k) - Y(k)) > very_tiny) then
               exit
            else
               write(*,*) 'reduce_helium_difference: Y_miscible(k) is too small' ! HHK
               Y_miscible(k) = -1d0
               khi_new = k - 1
               is_miscible(k) = .true.
            end if
         end do

      khi = khi_new
      write(*,*) 'reduce_helium_difference: klo, khi, khi_new: ', klo, khi,  khi_new
      deallocate(Y)
   end subroutine trim_region

   subroutine join_small_precipitating_regions(s, klo, khi, n_precipitating_regions, Y_miscible, is_miscible, ierr)
      ! if regions are very small any only differ by a few zones, artifacts form.
      ! These regions are likely also not physical, because we interpolate the (uncertain) phase diagram linearly.

      type (star_info), pointer :: s
      integer, allocatable, intent(inout) :: klo(:), khi(:)
      integer, intent(inout) :: n_precipitating_regions
      logical, allocatable, intent(inout) :: is_miscible(:)
      real(kind=dp), allocatable, intent(inout) :: Y_miscible(:)
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
         Y_miscible(klo(i):khi(i)) = merge(s%xa(s%net_iso(chem_get_iso_id('he4')),klo(i):khi(i)), Y_miscible(klo(i):khi(i)), Y_miscible(klo(i):khi(i)) < 0d0)
      end do

      deallocate(klo_temp, khi_temp, was_joined)

   end subroutine join_small_precipitating_regions

   subroutine determine_precipitating_regions(s, precipitating_regions, ierr)
      type (star_info), pointer :: s
      type (immiscibility_region), allocatable, intent(out) :: precipitating_regions(:)
      integer :: i, k, n_precipitating_regions

      ! precipitating region vars
      integer, allocatable :: klo(:), khi(:)
      real(kind = dp), allocatable :: dY(:)

      ! utils
      integer, intent(out) :: ierr

      real(kind=dp) :: X, Y, Z
      include 'formats'

      ierr = 0 ! 0 means AOK


      call get_precipitating_bounds(s, klo, khi, n_precipitating_regions, dY, ierr)
      ! call dummy_get_precipitating_bounds(s, klo, khi, n_precipitating_regions, dY)
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
         allocate(precipitating_regions(i)%dY(1:s%nz), source=0d0)

         precipitating_regions(i)%dY(klo(i):khi(i)) = dY(klo(i):khi(i))

      end do
      !$OMP END PARALLEL DO

      if ( debug ) then
         do i = 1, n_precipitating_regions
            write(*,*) 'determine_precipitating_region: (i, klo, khi): ', i, klo(i), khi(i)
            do k = klo(i), khi(i)
               call get_mass_fractions(s, s% xa(:, k), X, Y, Z)
               write(*,*) 'determine_precipitating_region: (k, X, Y, Z, dY): ', &
                  i, k, X, Y, Z, precipitating_regions(i)%dY(k)
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

      deallocate(dY, klo, khi)
   end subroutine determine_precipitating_regions


   subroutine determine_helium_deposites(s, precipitating_regions, ierr)
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
      ! I thought a bit about this. I could of course add the helium to the radiative region, which would then trigger a convective instability.
      ! However, I think one can also make the argument that the rain droplets would fall until they reach a convective region, where they would be destroyed by convection.
      ! Both are too simplistic and I should solve the advection equation for the helium. However, I think as a starting point its fine.
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

         if ( debug ) write(*,*) 'determine_helium_deposites: (i, k_deposite_lo, k_deposite_hi, klo_local, khi_max): ', i, &
            precipitating_regions(i)%k_deposite_lo, precipitating_regions(i)%k_deposite_hi, klo_local, khi_max

         if (s%mixing_type(precipitating_regions(i)%k_deposite_lo) /= convective_mixing) then
            do k = (precipitating_regions(i)%k_deposite_lo+1), khi_max
               if (s%mixing_type(k) == convective_mixing) then
                  klo_local = k
                  exit
               end if
            end do

            if (klo_local == precipitating_regions(i)%k_deposite_lo) then
               if (verbose) write(*,*) 'determine_helium_deposites: depositing region not found. Using the entire region.'
               precipitating_regions(i)%k_deposite_lo = precipitating_regions(i)%khi + 1
               precipitating_regions(i)%k_deposite_hi = khi_max
               cycle
            end if
         end if

         if ( debug ) write(*,*) 'determine_helium_deposites: (i, klo_local, khi_max): ', i, klo_local, khi_max

         ! pressure scale height condition

         do k = (klo_local+1), khi_max

            ! to avoid very small deposites, we require that the helium spreads at least a certain distance
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
            write(*,*) 'determine_helium_deposites (region, k_reservoir_lo, k_reservoir_hi, klo, khi, k_deposite_lo, k_deposite_hi): ', &
               i, &
               precipitating_regions(i)%k_reservoir_lo, precipitating_regions(i)%k_reservoir_hi,&
               precipitating_regions(i)%klo, precipitating_regions(i)%khi, &
               precipitating_regions(i)%k_deposite_lo, precipitating_regions(i)%k_deposite_hi
         end do
      endif

   end subroutine determine_helium_deposites

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

         ! find the top convective helium reservoir
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
      real(kind=dp), allocatable :: Y_crit(:)
      logical, allocatable :: is_miscible(:)
      integer :: klo, k_reservoir_lo, k_reservoir_hi
      ierr = 0 ! 0 means AOK

      k_reservoir_lo = precipitating_region%k_reservoir_lo
      k_reservoir_hi = precipitating_region%k_reservoir_hi

      ! if k_reservoir_lo == k_reservoir_hi, we don't need to do anything because the zone above the precipitating region is not convective
      if (k_reservoir_lo == k_reservoir_hi) return

      allocate(Y_crit(k_reservoir_lo:k_reservoir_hi), source=s%xa(i_Y, precipitating_region%klo) + precipitating_region%dY(precipitating_region%klo))

      ! make sure cells change only by a maximum of gentle_mixing_maximum_helium_rain_dY if gentle mixing is on
      if ( s% do_gentle_mixing ) then
         allocate(is_miscible(k_reservoir_lo:k_reservoir_hi), source=.false.)
         call reduce_helium_difference(s, Y_crit, is_miscible, k_reservoir_lo, k_reservoir_hi, ierr)
         deallocate(is_miscible)
         if (ierr /= 0) return
      end if

      precipitating_region%dY(k_reservoir_lo:k_reservoir_hi) = merge(0d0, Y_crit - s%xa(i_Y, k_reservoir_lo:k_reservoir_hi), s%xa(i_Y, k_reservoir_lo:k_reservoir_hi) < Y_crit(k_reservoir_lo:k_reservoir_hi))
      deallocate(Y_crit)

      if ( debug ) then
         write(*,*) 'set_top_reservoir_abundance_change: Y(klo) = ', s%xa(i_Y, precipitating_region%klo)
         write(*,*) 'set_top_reservoir_abundance_change: dY(klo) = ', precipitating_region%dY(precipitating_region%klo)
         write(*,*) 'set_top_reservoir_abundance_change: Y_crit(klo) = ', s%xa(i_Y, precipitating_region%klo) + precipitating_region%dY(precipitating_region%klo)
         write(*,*) 'set_top_reservoir_abundance_change: Y(k_reservoir_lo:k_reservoir_hi) = ', s%xa(i_Y, k_reservoir_lo:k_reservoir_hi)
         write(*,*) 'set_top_reservoir_abundance_change: dY(k_reservoir_lo:k_reservoir_hi) = ', precipitating_region%dY(k_reservoir_lo:k_reservoir_hi)
      end if

   end subroutine set_top_reservoir_abundance_change

   subroutine set_Helium_deposit_region_abundance_change(s, precipitating_region, ierr)
      type (star_info), pointer :: s
      type (immiscibility_region) :: precipitating_region
      integer, intent(out) :: ierr

      ! local
      real(kind=dp) :: delta_M_Y_rain, M_deposite
      integer :: k_reservoir_lo, k_deposite_lo, k_deposite_hi, khi
      logical :: debug = .false.
      ierr = 0 ! 0 means AOK

      k_reservoir_lo = precipitating_region%k_reservoir_lo
      khi = precipitating_region%khi
      k_deposite_lo = precipitating_region%k_deposite_lo
      k_deposite_hi = precipitating_region%k_deposite_hi

      if ( debug ) then
         write(*,*) 'set_Helium_deposit_region_abundance_change: k_reservoir_lo = ', k_reservoir_lo
         write(*,*) 'set_Helium_deposit_region_abundance_change: khi = ', khi
         write(*,*) 'set_Helium_deposit_region_abundance_change: k_deposite_lo = ', k_deposite_lo
         write(*,*) 'set_Helium_deposit_region_abundance_change: k_deposite_hi = ', k_deposite_hi
      endif

      ! get the rain mass of He (Y)
      delta_M_Y_rain = - dot_product(s%dm(k_reservoir_lo:khi), precipitating_region%dY(k_reservoir_lo:khi))

      ! mass abundance change in deeper layers
      M_deposite = sum(s%dm(k_deposite_lo:k_deposite_hi))
      precipitating_region%dY(k_deposite_lo:k_deposite_hi) = delta_M_Y_rain/M_deposite

   end subroutine set_Helium_deposit_region_abundance_change

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

   subroutine load_phase_diagram(pd, ierr)
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

   subroutine free_phase_diagram(pd)
      type(phase_diagram), intent(inout) :: pd
      

      if (allocated(pd%He_rich)) deallocate(pd%He_rich)
      if (allocated(pd%He_poor)) deallocate(pd%He_poor)
      if (allocated(pd%T)) deallocate(pd%T)
      if (allocated(pd%P)) deallocate(pd%P)
      if (allocated(pd%unique_P)) deallocate(pd%unique_P)
      if (allocated(pd%unique_P_index)) deallocate(pd%unique_P_index)
      if (allocated(pd%num_Ts_at_unique_Ps)) deallocate(pd%num_Ts_at_unique_Ps)

   end subroutine free_phase_diagram

   function linear_interpolation(x, x1, x2, y1, y2) result(y)
      real(kind=dp), intent(in) :: x, x1, x2, y1, y2
      real(kind=dp) :: y, fraction
      fraction = (x - x1) / (x2 - x1)
      y = y1 * (1d0 - fraction) + y2 * fraction
   end function linear_interpolation

   function linear_list_interpolation(x, x_list, y_list, extrapolate, ierr) result(y)
      real(kind=dp), intent(in) :: x
      real(kind=dp), dimension(:), intent(in) :: x_list, y_list
      logical, intent(in) :: extrapolate
      real(kind=dp) :: y
      integer, intent(out) :: ierr
      integer :: i
      ierr = 0 ! 0 means AOK

      ! default value
      y = 0d0

      ! check that x_list and y_list are the same size
      if (size(x_list) /= size(y_list)) then
         ierr = 1
         write(*,*) 'linear_list_interpolation: x_list and y_list are not the same size'
         return
      end if

      ! check that x_list containts at least 2 points
      if (size(x_list) < 2) then
         ierr = 1
         write(*,*) 'linear_list_interpolation: x_list is too small'
         return
      end if

      ! test that x is strictly increasing
      do i = 1, size(x_list)-1
         if (x_list(i) >= x_list(i+1)) then
            ierr = 1
            write(*,*) 'linear_list_interpolation: x_list is not strictly increasing'
            return
         end if
      end do

      if (x < x_list(1)) then
         
         if ( extrapolate ) then
            y = linear_interpolation(x, x_list(1), x_list(2), y_list(1), y_list(2))
         else
            y = y_list(1)
         end if
         return
      
      else if (x > x_list(size(x_list))) then
         
         if ( extrapolate ) then
            y = linear_interpolation(x, x_list(size(x_list)-1), x_list(size(x_list)), y_list(size(y_list)-1), y_list(size(y_list)))
         else
            y = y_list(size(y_list))
         end if
         return
      end if
      
      ! find the interval
      i = 1
      do while (x_list(i) < x)
         i = i + 1
      end do

      ! linear interpolation
      y = linear_interpolation(x, x_list(i-1), x_list(i), y_list(i-1), y_list(i))
   end function linear_list_interpolation

   subroutine immiscible_helium_fraction_interpolation(P0, T0, pd, Y_min, Y_max, ierr)
      real(kind=dp), intent(in) :: P0, T0
      type(phase_diagram), intent(in) :: pd
      real(kind=dp), intent(out) :: Y_min, Y_max
      real(kind=dp) :: P_lo, P_hi
      real(kind=dp) :: Y_min_lo_interp, Y_min_hi_interp, Y_max_lo_interp, Y_max_hi_interp
      real(kind=dp), allocatable :: T_lo(:), T_hi(:), Y_min_lo(:), Y_min_hi(:), Y_max_lo(:), Y_max_hi(:)
      integer, intent(out) :: ierr
      integer :: i, j
      logical :: extrapolate = .false.
      include 'formats'
      ierr = 0 ! 0 means AOK

      ! -1d0 means out of bounds and thus by default miscible
      Y_min = -1d0
      Y_max = -1d0
      ! in principle, you should check this. However, when you paralalize the code, this
      ! will create a lot of overhead. So, I just assume that the phase diagram is correct for now
      ! call check_phase_diagram(pd, ierr)
      ! if (ierr /= 0) return

      ! find the pressure values that enclose the P0

      ! edge cases
      if (P0 < pd%unique_P(1) .or. P0 > pd%unique_P(pd%num_unique_Ps)) then
         if ( debug_phase_diagram ) write(*,*) 'immiscible_helium_fraction_interpolation: P0 is outside the phase diagram. P0 = ', P0
         return
      end if

      j = -1
      do i = 2, pd%num_unique_Ps
         if ( P0 < pd%unique_P(i) ) then
            j = i - 1
            exit
         end if
      end do

      if (j == -1) then
         ierr = 1
         write(*,*) 'immiscible_helium_fraction_interpolation: error in finding the pressure index'
         return
      end if

      P_lo = pd%unique_P(j)
      P_hi = pd%unique_P(j+1)

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
            write(*,*) 'immiscible_helium_fraction_interpolation: Tj is not strictly increasing'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
      end do

      do i = 1, size(T_hi)-1
         if (T_hi(i) >= T_hi(i+1)) then
            ierr = 1
            write(*,*) 'immiscible_helium_fraction_interpolation: Tjp1 is not strictly increasing'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
      end do

      ! if T0 below the temperature range, return -1d0
      if (T0 < T_lo(1) .or. T0 < T_hi(1)) then
         if ( debug_phase_diagram ) write(*,*) 'immiscible_helium_fraction_interpolation: T0 below the data range T0 = ', T0
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      ! if T0 is above the temperature range, extrapolate linearly (done in lineater_list_interpolation by default)

      ! Y_min
      Y_min_lo_interp = linear_list_interpolation(T0, T_lo, Y_min_lo, extrapolate, ierr)
      if (ierr /= 0) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      Y_min_hi_interp = linear_list_interpolation(T0, T_hi, Y_min_hi, extrapolate, ierr)
      if (ierr /= 0) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      Y_min = linear_interpolation(P0, P_lo, P_hi, Y_min_lo_interp, Y_min_hi_interp)

      if ( extrapolate ) then
         Y_min = min(1d0, Y_min)
      else
         if (Y_min < 0d0) then
            ierr = 1
            write(*,*) 'immiscible_helium_fraction_interpolation: Y_min_0 < 0'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
   
         if (Y_min > 1d0) then
            ierr = 1
            write(*,*) 'immiscible_helium_fraction_interpolation: Y_min_0 > 1'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
      end if
      
      ! Y_max
      Y_max_lo_interp = linear_list_interpolation(T0, T_lo, Y_max_lo, extrapolate, ierr)
      if (ierr /= 0) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      Y_max_hi_interp = linear_list_interpolation(T0, T_hi, Y_max_hi, extrapolate, ierr)
      if (ierr /= 0) then
         deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
         return
      end if

      Y_max = linear_interpolation(P0, P_lo, P_hi, Y_max_lo_interp, Y_max_hi_interp)
      
      ! if extrapolation is allowed, the linear functions might go outside the range of [0, 1]
      ! in this case, we set the values to 0 and 1
      ! if extrapolation is not allowed, we need to check that the values are in the range of [0, 1]
      ! and return an error if they are not
      if (extrapolate) then
         Y_max = max(0d0, Y_max)

      else
         if (Y_max < 0d0) then
            ierr = 1
            write(*,*) 'immiscible_helium_fraction_interpolation: Y_max_0 < 0'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
         
         if (Y_max > 1d0) then
            ierr = 1
            write(*,*) 'immiscible_helium_fraction_interpolation: Y_max_0 > 1'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if

         if (Y_min > Y_max) then
            ierr = 1
            write(*,*) 'immiscible_helium_fraction_interpolation: Y_min > Y_max'
            deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)
            return
         end if
      end if

      if ( debug_phase_diagram ) then
         write(*,*) '----'
         write(*,*) 'immiscible_helium_fraction_interpolation: P0(Mbar) = ', P0/1d12
         write(*,*) 'immiscible_helium_fraction_interpolation: T0 = ', T0
         write(*,*) 'immiscible_helium_fraction_interpolation: P_lo(Mbar) = ', P_lo/1d12
         write(*,*) 'immiscible_helium_fraction_interpolation: P_hi(Mbar) = ', P_hi/1d12
         write(*,*) '----'
         write(*,*) 'immiscible_helium_fraction_interpolation: T_lo = ', T_lo
         write(*,*) 'immiscible_helium_fraction_interpolation: Y_min_lo = ', Y_min_lo
         write(*,*) 'immiscible_helium_fraction_interpolation: Y_lo_interp = ', Y_min_lo_interp
         write(*,*) '----'
         write(*,*) 'immiscible_helium_fraction_interpolation: T_hi = ', T_hi
         write(*,*) 'immiscible_helium_fraction_interpolation: Y_min_hi = ', Y_min_hi
         write(*,*) 'immiscible_helium_fraction_interpolation: Y_hi_interp = ', Y_min_hi_interp
         write(*,*) '----'
      end if

      deallocate(T_lo, T_hi, Y_min_lo, Y_min_hi, Y_max_lo, Y_max_hi)

   end subroutine immiscible_helium_fraction_interpolation

   subroutine Y_is_miscible(lnP, lnT, Y, is_miscible, Y_miscible, pd, ierr)
      ! gives the highest helium mass fraction that is miscible before the immiscibility gap
      real(kind=dp), intent(in) :: lnP, lnT, Y
      type(phase_diagram), intent(in) :: pd
      real(kind=dp), intent(out) :: Y_miscible
      logical, intent(out) :: is_miscible
      real(kind=dp) :: P0, T0
      integer, intent(out) :: ierr
      real(kind=dp) :: Y_max

      ierr = 0 ! 0 means AOK

      ! convert to lin (already in cgs)
      P0 = exp(lnP)
      T0 = exp(lnT)

      call immiscible_helium_fraction_interpolation(P0, T0, pd, Y_miscible, Y_max, ierr)
      if (ierr /= 0) return

      ! check that Y is in the miscible range
      if (Y_max < Y_miscible) then  ! outside of the miscible range
         is_miscible = .true.
      else if (Y <= Y_miscible .or. Y >= Y_max) then
         is_miscible = .true.
      else
         is_miscible = .false.
      end if

      if ( debug .and. .not. is_miscible) then
         write(*,*) 'Y_is_miscible: P0 = ', P0/1d12
         write(*,*) 'Y_is_miscible: T0 = ', T0
         write(*,*) 'Y_is_miscible: Y = ', Y
         write(*,*) 'Y_is_miscible: Y_min = ', Y_miscible
         write(*,*) 'Y_is_miscible: Y_max = ', Y_max
      end if

   end subroutine Y_is_miscible

   subroutine adjust_mass_fractions(s,  dY, ierr)
      type (star_info), pointer :: s
      real(kind=dp), intent(in) :: dY(:)
      integer, intent(out) :: ierr
      integer :: i, k, n_species
      real(kind=dp) :: f
      ierr = 0 ! 0 means AOK
      
      n_species = size(s%xa, 1)
      ! test that dY is the same size as s%xa
      if (size(dY) /= s%nz) then
         write(*,*) 'adjust_mass_fractions: dY is not the same size as s%xa'
         ierr = 1
         return
      end if

      !$OMP PARALLEL DO PRIVATE(i, k, f) SCHEDULE(guided)
      do i = 1, n_species
         if (i == i_Y) cycle
         do k = 1, s% nz
         f = 1d0-dY(k)/(1d0 - s% xa(i_Y, k))
         s% xa(i, k) = f * s% xa(i, k)
         end do
      end do
      !$OMP END PARALLEL DO

      s% xa(i_Y, 1:s% nz) = s% xa(i_Y, 1:s% nz) + dY(1:s% nz)

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

   subroutine get_settling_velocity(s, v_settle, is_miscible, Y_miscible, ierr)
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(out) :: v_settle(:)
      logical, allocatable, intent(in) :: is_miscible(:)
      real(kind=dp), allocatable, intent(in) :: Y_miscible(:)
      integer, intent(out) :: ierr

      ! local variables
      real(kind=dp) :: v_settle_0 = 1d3 ! fixed setteling veliocity in cgs
      ierr = 0 ! 0 means AOK

      allocate(v_settle(s%nz), source=0d0)

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


end module helium_rain
