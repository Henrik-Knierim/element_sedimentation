! ***********************************************************************
!   Copyright (C) 2025  Henrik Knierim
! ***********************************************************************
module element_sedimentation

   use star_private_def
   use const_def
   use chem_lib

   implicit none

   private

   public :: shutdown_phase_diagram
   public :: do_element_sedimentation

   integer :: sedimentation_method = -1
   integer, parameter :: instant_rain = 0
   integer, parameter :: bottom_up_rain = 1
   integer, parameter :: advection_diffusion_rain = 2

   real(dp), parameter :: very_tiny = 1d-12
   real(dp), parameter :: tiny = 1d-10
   real(dp), parameter :: very_small = 1d-6
   real(dp), parameter :: small = 1d-4

   real(dp), parameter :: temperature_shift = 0d0
   character(len=32), parameter :: phase_diagram_elements(2) = [character(len=32) :: 'h1', 'he4']

   real(kind=dp), parameter :: minimal_pressure_scale_height = 1d0   ! the minimal size below which we join precipitating regions.

   integer :: i_Xi = -1
   logical :: debug = .false.
   logical :: debug_advection_diffusion_solver = .false.
   logical :: debug_number_fraction_scaling = .false.
   logical :: debug_phase_diagram = .false.
   logical :: debug_set_dXi = .false.
   logical :: debug_bottom_up_rain = .false.
   logical :: verbose = .false.
   logical :: custom_debug = .false. ! for dirty debugging
   logical :: extrapolate = .false.
   logical :: damp_sedimentation ! originally, damping was determined via do_gentle_mixing. Since gentle_mixing is not officially supported, I decided to make this module self-contained.

   ! monotonicity methods
   integer, parameter :: enforce_monotonicity_bottom_up = -1
   integer, parameter :: enforce_no_monotonicity = 0
   integer, parameter :: enforce_monotonicity_top_down = 1

   integer, parameter :: monotonic_increase = 1
   integer, parameter :: no_monotonicity = 0
   integer, parameter :: monotonic_decrease = -1

   type phase_diagram
      character(len=1024) :: phase_diagram_filename
      integer :: n_points
      real(kind=dp), allocatable :: He_rich(:), He_poor(:), T(:), P(:)
      integer :: num_unique_Ps
      integer, allocatable :: unique_P_index(:), num_Ts_at_unique_Ps(:)
      real(kind=dp), allocatable :: unique_P(:)
   end type phase_diagram

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
      immiscibility_phase_diagram% phase_diagram_filename = trim(s% element_sedimentation_phase_diagram)
      call load_phase_diagram(s, immiscibility_phase_diagram, ierr)
      if (ierr /= 0) then
         write(*,*) 'init_element_sedimentation: error in load_phase_diagram'
         return
      end if
      phase_diagram_loaded = .true.

      extrapolate = s% element_sedimentation_extrapolate_phase_diagram

      if (debug_phase_diagram) call print_phase_diagram(immiscibility_phase_diagram)

      ! set index of sedimenting species
      i_Xi = s% net_iso(chem_get_iso_id(trim(s% element_sedimentation_chemical_species)))

   end subroutine init_element_sedimentation

   subroutine do_element_sedimentation(s, ierr)

      type (star_info), pointer :: s
      integer, intent(out) :: ierr

      ierr = 0 ! 0 means AOK

      ! if not already done, initialize the element sedimentation
      if (.not. phase_diagram_loaded) call init_element_sedimentation(s, ierr)
      if (ierr /= 0) then
         write(*,*) 'do_element_sedimentation: error in init_element_sedimentation'
         return
      end if
      ! get the method to use for the element sedimentation
      ! we first assign the method and then select the correct subroutine
      ! to be able to reference the method in the subroutine later if needed (currently not the case)
      if (trim(s% element_sedimentation_method) == 'instant_rain') then
         sedimentation_method = instant_rain
      else if (trim(s% element_sedimentation_method) == 'bottom_up_rain') then
         sedimentation_method = bottom_up_rain
      else if (trim(s% element_sedimentation_method) == 'advection_diffusion_rain') then
         sedimentation_method = advection_diffusion_rain
      else
         ierr = 1
         write(*,*) 'do_element_sedimentation: unknown method for element sedimentation: ', trim(s% element_sedimentation_method)
         return
      end if

      damp_sedimentation = (s% element_sedimentation_maximum_dXi_per_cell > 0d0)

      select case (sedimentation_method)
       case (instant_rain)
         call do_instant_rain(s, ierr)
       case (bottom_up_rain)
         call do_bottom_up_rain(s, ierr)
       case (advection_diffusion_rain)
         ! call do_advection_diffusion_rain(s, ierr)
       case default
         ierr = 1
         write(*,*) 'do_element_sedimentation: unknown method for element sedimentation: ', trim(s% element_sedimentation_method)
         return
      end select

      call update_model_(s, 1, s% nz, .true.)
      s% need_to_setvars = .true.

   end subroutine do_element_sedimentation

   ! subroutine do_advection_diffusion_rain(s, ierr)
   !    type (star_info), pointer :: s
   !    integer, intent(out) :: ierr

   !    ! local
   !    integer, parameter :: max_iterations = 100
   !    integer :: i, j
   !    type(advection_diffusion_state) :: adv_diff_state

   !    ierr = 0 ! 0 means AOK


   !    do i = 1, max_iterations
   !       if (verbose) write(*,*) 'do_advection_diffusion_rain: iteration: ', i

   !       ! initialize the advection-diffusion state
   !       call init_advection_diffusion_state(s, adv_diff_state, ierr)
   !       if (ierr /= 0) return

   !       ! if adv_diff_state % v is not allocated, no precipitating regions were found; exit the loop
   !       if (.not. allocated(adv_diff_state % v)) then
   !          if (debug_advection_diffusion_solver) write(*,*) 'do_advection_diffusion_rain: no precipitating regions found, skipping advection-diffusion step'
   !          call free_advection_diffusion_state(adv_diff_state)
   !          exit
   !       end if

   !       ! perform the advection-diffusion step
   !       call backward_Euler_step(adv_diff_state, ierr)
   !       if (ierr /= 0) return

   !       if (debug_advection_diffusion_solver) then
   !          do j = 1, adv_diff_state%nz
   !             write(*,*) 'do_advection_diffusion_rain: (k, Xi(k), dXi(k)) = ', j, adv_diff_state%Xi(j), adv_diff_state%dXi(j)
   !          end do
   !       end if

   !       call adjust_mass_fractions(s, adv_diff_state% dXi, ierr)

   !       ! check if the change in the element's mass fraction is below a threshold
   !       if (all(abs(adv_diff_state%dXi) < sedimentation_threshold)) then
   !          if (verbose) write(*,*) 'do_advection_diffusion_rain: no significant change in element mass fraction'
   !          if (verbose) write(*,*) 'do_advection_diffusion_rain: max(|dXi|) = ', maxval(abs(adv_diff_state%dXi))
   !          exit
   !       end if

   !       call free_advection_diffusion_state(adv_diff_state)
   !    end do

   !    ! TODO: Adjust the mixing in the sedimenting region
   !    stop
   !    ! call adjust_mixing_in_sedimenting_region(s, precipitating_regions, ierr)


   ! end subroutine do_advection_diffusion_rain

   ! subroutine init_advection_diffusion_state(s, adv_diff_state, ierr)
   !    type (star_info), pointer :: s
   !    type(advection_diffusion_state), intent(out) :: adv_diff_state
   !    integer, intent(out) :: ierr

   !    ! local
   !    real(kind=dp), allocatable :: Xi_miscible(:)  ! miscible element mass fraction
   !    logical, allocatable :: is_miscible(:)  ! miscible regions
   !    real(kind=dp), allocatable :: v_settle(:)  ! settling velocity
   !    integer :: k

   !    ierr = 0 ! 0 means AOK

   !    ! function to populate a advection_diffusion_state type used for the advection-diffusion solver
   !    adv_diff_state % nz = s% nz
   !    allocate(adv_diff_state % D(1:s% nz), source=s% D_mix(1: s% nz))  ! diffusion coefficient
   !    allocate(adv_diff_state % Xi(1:s% nz), source=s% xa(i_Xi, 1: s% nz))  ! diffusion coefficient
   !    allocate(adv_diff_state % dXi(1:s% nz), source=0d0)  ! diffusion coefficient

   !    ! radial step size
   !    allocate(adv_diff_state % dr(1:s% nz))
   !    !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(guided)
   !    do k = 1, s % nz
   !       adv_diff_state % dr(k) = get_dr(s, k)
   !    end do
   !    !$OMP END PARALLEL DO

   !    call get_raw_miscibility(s, is_miscible, Xi_miscible,  ierr)
   !    if (ierr /= 0) then
   !       write(*,*) 'init_advection_diffusion_state: error in get_raw_miscibility'
   !       return
   !    end if

   !    ! if everything is miscible, we can skip the advection-diffusion step
   !    if (all(is_miscible)) then
   !       if (debug_advection_diffusion_solver) write(*,*) 'init_advection_diffusion_state: all regions are miscible, skipping advection-diffusion step'
   !       return
   !    end if

   !    call get_settling_velocity(s, adv_diff_state % v, is_miscible, Xi_miscible, ierr)
   !    if (ierr /= 0) then
   !       write(*,*) 'init_advection_diffusion_state: error in get_settling_velocity'
   !       return
   !    end if

   !    call get_time_step(s, adv_diff_state, ierr)
   !    if (ierr /= 0) then
   !       write(*,*) 'init_advection_diffusion_state: error in get_time_step'
   !       return
   !    end if

   !    if (debug_advection_diffusion_solver) then
   !       write(*,*) 'init_advection_diffusion_state: advection-diffusion state initialized'
   !       write(*,*) 'init_advection_diffusion_state: nz = ', adv_diff_state % nz
   !       write(*,*) 'init_advection_diffusion_state: dt = ', adv_diff_state % dt
   !       do k = 1, adv_diff_state % nz
   !          write(*,'(a,i12,99(1pd26.16))') 'init_advection_diffusion_state: (k, D(k), v(k), dr(k), Xi(k)) = ', k, adv_diff_state % D(k), adv_diff_state % v(k), adv_diff_state % dr(k), adv_diff_state % Xi(k)
   !       end do
   !    end if

   !    deallocate(Xi_miscible, is_miscible)

   ! end subroutine init_advection_diffusion_state

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
      integer :: k

      ! debugging
      integer :: k_debug
      real(kind=dp), allocatable :: xa_pre(:,:)
      real(kind=dp), allocatable :: dXi(:)  ! change in the element's mass fraction
      real(kind=dp) :: dXi_max
      logical, allocatable :: is_miscible(:)  ! miscible regions
      integer :: k_update
      real(kind=dp) :: dXi_debug

      include 'formats'

      ierr = 0 ! 0 means AOK

      allocate(xa_pre(size(s%xa, 1), size(s%xa, 2)))

      ! in the instant rain method, we transfer all the abundaces at once
      ! this can lead to a region raining out that previously was not raining out
      ! thus, we need to loop until we find no significant change in the element's mass fraction

      xa_pre = s%xa

      do k_update = s% nz, 1, -1

         if ( debug_bottom_up_rain ) write(*,*) 'do_bottom_up_rain: k_update: ', k_update

         call get_dXi(s, k_update, s% nz, dXi, is_miscible, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_bottom_up_rain: error in get_dXi'
            if (allocated(dXi)) deallocate(dXi)
            if (allocated(is_miscible)) deallocate(is_miscible)
            return
         end if

         if (.not. allocated(dXi)) then
            cycle
         end if

         call adjust_mass_fractions(s, dXi(1:s% nz), ierr)
         if (ierr /= 0) then
            write(*,*) 'do_bottom_up_rain: error in adjust_mass_fractions'
            deallocate(dXi, is_miscible)
            return
         end if

         ! check mass fraction conservation
         call check_after_rain(s, 1, s% nz, xa_pre, ierr)

         ! adjust model after mixing
         call update_model_(s, k_update, s% nz, .true.)

         deallocate(dXi, is_miscible)

      end do


      ! call adjust_mixing_in_sedimenting_region(s, is_miscible, ierr)




      deallocate(xa_pre)

   end subroutine do_bottom_up_rain

   integer function first_nonzero_element(x)
      ! finds the first non-zero element in the array x
      real(kind=dp), intent(in) :: x(:)
      integer :: i

      first_nonzero_element = -1  ! default value if no non-zero element is found
      do i = 1, size(x)
         if (abs(x(i)) > 0d0) then
            first_nonzero_element = i
            return
         end if
      end do
   end function first_nonzero_element

   integer function last_nonzero_element(x)
      ! finds the last non-zero element in the array x
      real(kind=dp), intent(in) :: x(:)
      integer :: i

      last_nonzero_element = -1  ! default value if no non-zero element is found
      do i = size(x), 1, -1
         if (abs(x(i)) > 0d0) then
            last_nonzero_element = i
            return
         end if
      end do
   end function last_nonzero_element




   subroutine do_instant_rain(s, ierr)
      type (star_info), pointer :: s
      integer :: i
      integer, intent(out) :: ierr

      integer, parameter :: max_iterations = 10
      integer :: j

      ! debugging
      integer :: k_debug
      real(kind=dp), allocatable :: xa_pre(:,:)
      real(kind=dp), allocatable :: dXi(:)  ! change in the element's mass fraction
      logical, allocatable :: is_miscible(:)  ! miscible regions
      include 'formats'

      ierr = 0 ! 0 means AOK

      allocate(xa_pre(size(s%xa, 1), size(s%xa, 2)))

      ! in the instant rain method, we transfer all the abundaces at once
      ! this can lead to a region raining out that previously was not raining out
      ! thus, we need to loop until we find no significant change in the element's mass fraction
      do j=1, max_iterations

         if (verbose) write(*,*) 'do_instant_rain: iteration: ', j
         xa_pre = s%xa

         call get_dXi(s, 1, s%nz, dXi, is_miscible, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_instant_rain: error in get_dXi'
            if (allocated(dXi)) deallocate(dXi)
            if (allocated(is_miscible)) deallocate(is_miscible)
            return
         end if

         ! if dXi is not allocated, no precipitating regions were found; exit the loop
         if (.not. allocated(dXi)) then
            if (verbose) write(*,1) 'do_instant_rain: no precipitating region found.'
            return
         end if

         call adjust_mass_fractions(s, dXi(1:s% nz), ierr)

         ! check mass fraction conservation
         call check_after_rain(s, 1, s% nz, xa_pre, ierr)

         ! call adjust_mixing_in_sedimenting_region(s, is_miscible, ierr)

         deallocate(dXi, is_miscible)

         if (all(abs(s%xa(i_Xi, 1:s% nz)-xa_pre(i_Xi,1:s% nz)) < small)) exit

      end do

      deallocate(xa_pre)

   end subroutine do_instant_rain

   subroutine check_after_rain(s, klo, khi, xa_pre, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: klo, khi
      real(kind=dp), allocatable, intent(in) :: xa_pre(:,:)
      integer, intent(out) :: ierr

      integer :: i, k, j
      integer :: k_debug
      integer :: k_RCB
      real(kind=dp) :: M_chem, M_chem_pre
      integer :: monotonicity_method
      logical :: check_range = .true.


      ierr = 0 ! 0 means AOK


      ! test if the total mass fraction was conserved
      do i = 1, s% species

         if ( i == i_Xi ) then
            monotonicity_method = monotonic_increase
         else
            monotonicity_method = monotonic_decrease
         end if

         M_chem = dot_product(s% xa(i, klo:khi), s% dm(klo:khi))
         M_chem_pre = dot_product(xa_pre(i, klo:khi), s% dm(klo:khi))
         if (abs((M_chem-M_chem_pre)/M_chem_pre) > small) then
            ierr = 1
            write(*,*) 'check_after_rain: mass fraction of species ', i, ' not conserved'
            write(*,*) 'klo, khi, M_chem/M_chem_pre, M_chem_pre, M_chem: ', &
               klo, khi, M_chem/M_chem_pre, M_chem_pre, M_chem
         end if

         ! TODO: Only check monotonicity in radiative regions.
         ! ? The error probably arises from ommiting mixing
         ! call find_RCB(s, klo, khi, k_RCB, ierr)
         call check_abundance_array(s% xa(i, 1:s% nz), klo, khi, check_range, monotonicity_method, ierr, mask = (s% mixing_type(klo:khi)/=no_mixing))
      end do

      if ( verbose ) then
         k_debug = maxloc(abs(s%xa(i_Xi, klo:khi)-xa_pre(i_Xi, klo:khi)), dim=1)
         if (verbose) write(*,*) 'maximum difference: (k, dXi)', k_debug,  abs(s%xa(i_Xi, k_debug)-xa_pre(i_Xi, k_debug))
      end if

      if ( ierr /= 0 ) then
         call mesa_error(__FILE__,__LINE__)
      end if
   end subroutine check_after_rain

   subroutine get_raw_miscibility(s, klo, khi, is_miscible, Xi_miscible_min, Xi_miscible_max, ierr)
      ! miscibility without any post-processing
      type(star_info), pointer :: s
      integer, intent(in) :: klo, khi  ! range of zones to consider
      real(kind=dp), allocatable, intent(out) :: Xi_miscible_min(:), Xi_miscible_max(:)   ! miscible element mass fraction
      logical, allocatable, intent(out) :: is_miscible(:)  ! miscible regions

      ! local vars
      real(kind=dp), allocatable :: yi_miscible_min(:), yi_miscible_max(:) ! miscible element _number_ fraction
      real(kind=dp), allocatable :: yi(:)  ! element number fractions
      real(kind=dp) :: X, Y, Z
      integer :: k
      integer :: op_err, ierr

      ierr = 0 ! 0 means AOK

      ! set up arrays
      ! we need to evaluate the entire range later
      ! however, only values ranging from klo to khi are checked for miscibility
      allocate(Xi_miscible_min(1:s% nz), Xi_miscible_max(1:s% nz), yi_miscible_min(1:s% nz), yi_miscible_max(1:s% nz), source = -1d0)
      allocate(is_miscible(1:s% nz), source =.true.)  ! miscible regions

      !$OMP PARALLEL DO PRIVATE(k, yi,  op_err) SCHEDULE(guided)
      do k = klo, khi

         op_err = 0  ! Each thread gets its own local copy of ierr
         call get_scaled_number_fractions(s, s% xa(:, k), yi, op_err)
         if (op_err /= 0) ierr = op_err

         call yi_is_miscible(s%lnPgas(k), s%lnT(k), yi(i_Xi), is_miscible(k), yi_miscible_min(k), yi_miscible_max(k), immiscibility_phase_diagram, op_err)
         if (op_err /= 0) ierr = op_err

         ! for y_min
         call from_scaled_number_fraction_to_mass_fraction(s, yi, yi_miscible_min(k), is_miscible(k), Xi_miscible_min(k), op_err)
         if (op_err /= 0) ierr = op_err

         ! for y_max
         call from_scaled_number_fraction_to_mass_fraction(s, yi, yi_miscible_max(k), is_miscible(k), Xi_miscible_max(k), op_err)
         if (op_err /= 0) ierr = op_err

         deallocate(yi)
      end do
      !$OMP END PARALLEL DO
      if (ierr /= 0) return

   end subroutine get_raw_miscibility

   subroutine get_dXi(s, klo, khi, dXi, is_miscible, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: klo, khi  ! range of zones to consider
      real(kind = dp), allocatable, intent(out) :: dXi(:)  ! difference between the miscible and the actual element mass fraction
      logical, allocatable, intent(out) :: is_miscible(:)  ! miscible regions
      integer, intent(out) :: ierr

      ! local
      integer :: k, i, op_err
      real(kind=dp), allocatable :: Xi_miscible_min(:), Xi_miscible_max(:)
      integer :: k_debug ! utils

      call get_raw_miscibility(s, klo, khi, is_miscible, Xi_miscible_min, Xi_miscible_max, ierr)
      if (ierr /= 0) then
         write(*,*) 'get_precipitating_bounds: error in get_raw_miscibility'
         return
      end if

      ! if is_miscible is .true. everywhere, we don't need to do anything
      if (all(is_miscible)) then
         deallocate(Xi_miscible_min, Xi_miscible_max, is_miscible)
         return
      end if

      ! get the difference between the miscible and the actual element mass fraction in the precipitating region
      call set_dXi(s, is_miscible, Xi_miscible_min, Xi_miscible_max, 1, s% nz, dXi, ierr)

   end subroutine get_dXi

   logical function is_at_miscibility_edge(Xi, is_miscible, Xi_miscible_min, Xi_miscible_max, dXi_tolerance)
      ! checks whether a zone is either immsicible or very close to the miscibility edge
      ! the idea is that any change in this cell would likely lead to a rainout
      real(kind=dp), intent(in) :: Xi  ! element mass fraction
      logical, intent(in) :: is_miscible  ! miscible regions
      real(kind=dp), intent(in) :: Xi_miscible_min  ! miscible element mass fraction
      real(kind=dp), intent(in) :: Xi_miscible_max  ! miscible element mass fraction
      real(kind=dp), intent(in), optional :: dXi_tolerance  ! tolerance for the difference between the miscible and the actual element mass fraction

      ! local
      real(kind=dp) :: dXi_tolerance_local

      is_at_miscibility_edge = .false.

      if ( .not. is_miscible ) then
         is_at_miscibility_edge = .true.  ! immiscible regions are always at the miscibility edge
         return
      end if

      ! Set up local mask - if mask is not present, create a mask that's all .true.
      if (present(dXi_tolerance)) then
         dXi_tolerance_local = dXi_tolerance
      else
         dXi_tolerance_local = 0.01d0
      end if

      if (Xi + dXi_tolerance_local > Xi_miscible_min .and. Xi_miscible_min >= 0d0) then
         is_at_miscibility_edge = .true.
         return
      end if

      if (Xi - dXi_tolerance_local < Xi_miscible_max .and. Xi_miscible_max >= 0d0 ) then
         is_at_miscibility_edge = .true.
         return
      end if

   end function is_at_miscibility_edge

   subroutine set_dXi(s, is_miscible, Xi_miscible_min, Xi_miscible_max, klo, khi, dXi, ierr)
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(in) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! miscible element mass fraction
      logical, allocatable, intent(in) :: is_miscible(:)  ! miscible regions
      integer, intent(in) :: klo, khi
      real(kind=dp), allocatable, intent(out) :: dXi(:)  ! difference between the miscible and the actual element mass fraction
      integer, intent(out) :: ierr

      ! local vars
      real(kind=dp), allocatable :: Xi_temp(:)  ! temporary array to hold the mass fractions
      real(kind=dp), allocatable :: Xi_miscible_min_monotonic(:)  ! miscible element mass fraction, monotonically increasing
      real(kind=dp), allocatable :: Xi_miscible_max_monotonic(:)  ! miscible element mass fraction, monotonically decreasing
      integer :: k, j, k_core, k_rain, k_RCB
      logical :: it_rains
      logical :: has_core  ! whether we have a core

      ierr = 0 ! 0 means AOK

      allocate(Xi_temp(klo:khi), source=s% xa(i_Xi, klo:khi))

      call init_set_dXi(s, klo, khi, is_miscible, Xi_miscible_min, Xi_miscible_max, k_RCB, k_rain, k_core, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
      if (ierr /= 0) then
         write(*,*) 'set_dXi: error in init_set_dXi'
         return
      end if

      has_core = .false.
      if (custom_debug) write(*,*) 'set_dXi: (k, is_miscible(k), Xi_temp(k), Xi_miscible_min_monotonic(k), Xi_miscible_max_monotonic(k)): '
      do k = k_RCB, khi

         if (custom_debug) write(*,'(i12,l5,99(1pd26.16))') &
            k, is_miscible(k), Xi_temp(k), Xi_miscible_min_monotonic(k), Xi_miscible_max_monotonic(k)

         if (should_exit_loop(k, k_core, k_rain, has_core)) exit

         it_rains = (Xi_temp(k) > Xi_miscible_min_monotonic(k)) & ! we are above the minimum mixture
            .and. (Xi_miscible_min_monotonic(k) >= 0d0)  ! the miscible limit is not -1
         !.and. (Xi_temp(k) < Xi_miscible_max_monotonic(k) .or. Xi_miscible_max_monotonic(k) <= 0d0) & ! we are below the maximum mixture

         if (it_rains) then
            if (debug_set_dXi) write(*,*) 'set_dXi: raining at k = ', k

            if (k + 1 == k_core) then
               ! we are one cell above the core and still immiscible.
               ! the only way forward is to saturate the core
               call transfer_Xi(s, Xi_temp, Xi_miscible_max_monotonic, k, k_rain+1, k-1, ierr)
               if (ierr /= 0) then
                  write(*,*) 'set_dXi: error in transfer_Xi for edge case k + 1 == k_core'
                  return
               end if
               exit
            end if

            ! if the mass fraction is in the immiscible region, transfer it to the deposite below
            call transfer_Xi(s, Xi_temp, Xi_miscible_min_monotonic, k, k+1, k_core-1, ierr)
            if (ierr /= 0) then
               write(*,*) 'set_dXi: error in transfer_Xi for rain condition'
               return
            end if

            ! make sure that after the transfer, none of the mass fractions are greater than the miscible maximum
            if (any((Xi_temp(k+1:k_core-1) > Xi_miscible_max_monotonic(k+1:k_core-1)) .and. (Xi_miscible_max_monotonic(k+1:k_core-1) >= 0d0))) then
               has_core = .true.
               call saturate_core(s, k, klo, khi, k_RCB, k_rain, k_core, is_miscible, Xi_miscible_min, Xi_miscible_max, Xi_temp, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
               if (ierr /= 0) return
            end if

         end if
      end do

      ! if there is a core, we homogenize the mass fraction one last time
      if (has_core) then
         if (debug_set_dXi) write(*,*) 'set_dXi: homogenizing region (k_rain+1, k_core-1): ', k_rain+1, k_core-1
         call homogenize_region(s, Xi_temp, k_rain+1, k_core-1, ierr)
         if ( ierr /= 0 ) return

         if (debug_set_dXi) then
            call check_abundance_array(Xi_temp, k_rain+1, khi, .true., monotonic_increase, ierr)
            if (ierr /= 0) then
               write(*,*) 'set_dXi: error in check_abundance_array after homogenization of the core'
               call mesa_error(__FILE__,__LINE__)
            end if
         end if
      end if

      ! homogenize convective envelope
      if (sedimentation_method == instant_rain .or. klo == 1 ) then
         call homogenize_convective_envelope(s, is_miscible, Xi_miscible_min, Xi_miscible_max, klo, khi, Xi_temp, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_dXi: error in homogenize_convective_envelope'
            return
         end if
      end if

      ! now, we compute the difference between the miscible and the actual element mass fraction in the region klo:khi
      if (allocated(dXi)) deallocate(dXi)
      allocate(dXi(klo:khi), source = 0d0)
      dXi(klo:khi) = Xi_temp(klo:khi) - s% xa(i_Xi, klo:khi)

      if ( custom_debug ) then
         write(*,*) ''
         write(*,*) 'After set_dXi:'
         write(*,*) 'klo, khi: ', klo, khi
         write(*,*) 'k, is_miscible, mixing_type, Xi, Xi_temp, dXi, Xi_miscible_min, Xi_miscible_min_monotonic, Xi_miscible_max, Xi_miscible_max_monotonic:'
         do k = klo, khi
            write(*,'(i12,l5,i12,99(1pd26.16))') k, is_miscible(k), s% mixing_type(k), s% xa(i_Xi, k), Xi_temp(k), dXi(k), &
               Xi_miscible_min(k), Xi_miscible_min_monotonic(k), Xi_miscible_max(k), Xi_miscible_max_monotonic(k)
         end do
         write(*,*) ''
      end if

      deallocate(Xi_temp)
      deallocate(Xi_miscible_min_monotonic)
      deallocate(Xi_miscible_max_monotonic)

   end subroutine set_dXi

   subroutine homogenize_convective_envelope(s, is_miscible, Xi_miscible_min, Xi_miscible_max, klo, khi, Xi, ierr)
      type (star_info), pointer :: s
      logical, allocatable, intent(in) :: is_miscible(:)  ! miscible regions
      real(kind=dp), allocatable, intent(in) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! miscible element mass fraction
      integer, intent(in) :: klo, khi  ! range of zones to consider
      real(kind=dp), allocatable, intent(inout) :: Xi(:)  ! element mass fraction
      integer, intent(out) :: ierr  ! error code

      integer :: k_RCB  ! index of the radiative-convective boundary
      integer :: k_core  ! index of the core
      integer :: k  ! loop index

      call find_RCB(s, klo, khi, k_RCB, ierr)
      if (ierr /= 0) then
         write(*,*) 'homogenize_convective_envelope: error in find_RCB'
         return
      end if

      do k = k_RCB+1, khi
         if (is_at_miscibility_edge(Xi(k), is_miscible(k), Xi_miscible_min(k), Xi_miscible_max(k))) then
            if (debug_set_dXi) write(*,*) 'set_dXi: homogenizing region (k_RCB, k): ', k_RCB, k-1
            call homogenize_region(s, Xi, k_RCB, k-1, ierr)
            if (ierr /= 0) return

            if (debug_set_dXi) then
               call check_abundance_array(Xi, k_RCB, khi, .true., monotonic_increase, ierr)
               if (ierr /= 0) then
                  write(*,*) 'set_dXi: error in check_abundance_array after homogenization'
                  call mesa_error(__FILE__,__LINE__)
               end if
            end if
            exit
         end if
      end do

   end subroutine homogenize_convective_envelope

   subroutine saturate_core(s, k, klo, khi, k_RCB, k_rain, k_core, is_miscible, Xi_miscible_min, Xi_miscible_max, Xi, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
      ! saturates the core of the star, i.e. the region where the element mass fraction is above the miscible maximum
      type (star_info), pointer :: s
      integer, intent(in) :: k  ! index of the last zone in the core
      integer, intent(in) :: klo, khi  ! entire modelling range
      integer, intent(in) :: k_RCB  ! index of the radiative-convective boundary
      integer, intent(inout) :: k_rain  ! index of the first zone in the rain
      integer, intent(inout) :: k_core  ! index of the core
      logical, allocatable, intent(in) :: is_miscible(:)  ! miscible regions
      real(kind=dp), allocatable, intent(in) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! miscible element mass fraction
      real(kind=dp), allocatable, intent(inout) :: Xi(:)  ! element mass fraction
      real(kind=dp), allocatable, intent(inout) :: Xi_miscible_min_monotonic(:), Xi_miscible_max_monotonic(:)  ! miscible element mass fraction, monotonically increasing/decreasing
      integer, intent(out) :: ierr  ! error code

      ! local variables
      integer :: j  ! loop index
      integer :: k_core_old  ! old value of k_core

      ierr = 0 ! 0 means AOK
      k_core_old = k_core

      if (debug_set_dXi) write(*,*) 'saturate_core: k_rain, k+1, k_core_old', k_rain, k+1, k_core_old
      do j = k_core_old-1, k+1, -1
         if (Xi(j) > Xi_miscible_max_monotonic(j) .and. Xi_miscible_max_monotonic(j) >= 0d0) then
            k_core = j
            if (debug) write(*,*) 'core oversaturated: j, Xi_temp, Xi_miscible_max_monotonic, k_core', j, Xi(j), Xi_miscible_max_monotonic(j), k_core

            if (k_rain + 1 == j) then
               call find_k_rain(s, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, is_miscible, k_RCB, k_rain - 1, k_rain, ierr)
               if ( ierr /= 0 ) return
               if (debug) write(*,*) 'set_dXi: k_rain adjusted to ', k_rain, ' because k_rain == j = ', j
            end if

            ! transfer the excess mass fraction to the miscible maximum
            call transfer_Xi(s, Xi, Xi_miscible_max_monotonic, j, k_rain+1, j-1, ierr)
            if (ierr /= 0) then
               write(*,*) 'set_dXi: error in transfer_Xi for oversaturated core'
               return
            end if
         else

            if (debug_set_dXi) write(*,*) 'exit core: k_core', k_core

            if (k_rain + 1 == j) then
               call find_k_rain(s, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, is_miscible, k_RCB, k_rain - 1, k_rain, ierr)
               if ( ierr /= 0 ) return
               if (debug_set_dXi) write(*,*) 'set_dXi: k_rain adjusted to ', k_rain, ' because k_rain == j = ', j
            end if

            call adjust_miscibility_gap(s, klo, khi, k_RCB, k_core-1, is_miscible, Xi_miscible_min, Xi_miscible_max, &
               Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)

            exit
         end if
      end do

   end subroutine saturate_core

   subroutine homogenize_region(s, Xi, klo, khi, ierr)
      ! homogenizes the element mass fraction Xi in the region klo:khi
      type (star_info), pointer :: s
      real(kind=dp), intent(inout) :: Xi(:)  ! element mass fraction
      integer, intent(in) :: klo, khi  ! region to be homogenized
      integer, intent(out) :: ierr  ! error code

      ! local variables
      real(kind=dp) :: M, M_Xi  ! total mass and total mass of the element in the region klo:khi
      integer :: k

      ierr = 0 ! 0 means AOK

      M = sum(s% dm(klo:khi))  ! total mass in the region klo:khi
      M_Xi = dot_product(s% dm(klo:khi), Xi(klo:khi))  ! total mass of the element in the region klo:khi

      if (M <= 0d0) then
         ierr = 1
         write(*,*) 'homogenize_region: total mass in the region klo:khi is zero or negative'
         write(*,*) 'klo, khi: ', klo, khi
         write(*,*) 'M: ', M
         return
      end if

      if ( M_Xi < 0d0 ) then
         ierr = 1
         write(*,*) 'homogenize_region: total mass of the element in the region klo:khi is negative'
         write(*,*) 'klo, khi: ', klo, khi
         write(*,*) 'M_Xi: ', M_Xi
         return
      end if

      ! now, we homogenize the element mass fraction in the region klo:khi
      Xi(klo:khi) = M_Xi / M

   end subroutine homogenize_region

   subroutine enforce_monotonicity(Xi, klo, khi, method, ierr)
      ! this subroutine enforces monotonicity of the miscible element mass fraction Xi
      real(kind=dp), allocatable, intent(inout) :: Xi(:)  ! miscible element mass fraction
      integer, intent(in) :: klo, khi  ! indices of Xi to be checked
      integer, intent(in) :: method  ! if true, enforce monotonicity from bottom to top, otherwise from top to bottom
      ! typically, you want to start from the bottom for Xi_min and from the top for Xi_max
      integer, intent(out) :: ierr

      ! local variables
      integer :: k
      logical :: check_range = .false.
      ierr = 0 ! 0 means AOK

      if (debug_set_dXi) write(*,*) 'enforce_monotonicity: size(Xi), klo, khi, method: ', size(Xi), klo, khi, method

      ! TODO: Optimize by finding the minimum/maximum value
      if (method == enforce_monotonicity_bottom_up) then
         do k = khi, klo+1, -1
            if ( Xi(k) >= 0d0 .and. Xi(k-1) < 0d0 ) then
               Xi(k-1) = Xi(k)
            else if ( Xi(k) < 0d0 ) then
               ! continue
            else if (Xi(k) < Xi(k-1)) then
               if ( debug ) write(*,*) 'enforce_monotonicity: enforcing monotonicity at k, Xi(k-1), Xi(k) :', &
                  k, Xi(k-1), Xi(k)
               Xi(k-1) = Xi(k)
            end if
         end do
      else if (method == enforce_monotonicity_top_down) then
         do k = klo+1, khi
            if (Xi(k) < 0d0 .and. Xi(k-1) >= 0d0 ) then
               Xi(k) = Xi(k-1)
            else if (Xi(k) < 0d0) then
               ! continue
            else if (Xi(k) < Xi(k-1)) then
               if ( debug ) write(*,*) 'enforce_monotonicity: enforcing monotonicity at k, Xi(k), Xi(k-1):', &
                  k, Xi(k), Xi(k-1)
               Xi(k) = Xi(k-1)
            end if
         end do
      else
         ierr = 1
         write(*,*) 'enforce_monotonicity: unknown method: ', method
         write(*,*) 'klo, khi: ', klo, khi
         return
      end if

      call check_abundance_array(Xi, klo, khi, check_range, monotonic_increase, ierr)
      if (ierr /= 0) then
         write(*,*) 'enforce_monotonicity: error in check_abundance_array'
         return
      end if

   end subroutine enforce_monotonicity

   subroutine adjust_miscibility_gap(s, klo, khi, k_start, k_end, is_miscible, Xi_miscible_min, Xi_miscible_max, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
      ! sets the abundance profiles for the miscible element
      type (star_info), pointer :: s
      integer, intent(in) :: klo, khi  ! region in which we transfer the element
      integer, intent(in) :: k_start, k_end  ! subregion for which we ensure monotonicity
      logical, allocatable,  intent(in) :: is_miscible(:)  ! miscible regions
      real(kind=dp), allocatable, intent(in) :: Xi_miscible_min(:)  ! miscible element mass fraction
      real(kind=dp), allocatable, intent(in) :: Xi_miscible_max(:)  ! miscible element mass fraction
      real(kind=dp), allocatable, intent(out) :: Xi_miscible_min_monotonic(:)  ! miscible element mass fraction, monotonically increasing
      real(kind=dp), allocatable, intent(out) :: Xi_miscible_max_monotonic(:)  ! miscible element mass fraction, monotonically decreasing
      integer, intent(out) :: ierr  ! error code

      ierr = 0 ! 0 means AOK

      if (allocated(Xi_miscible_min_monotonic)) deallocate(Xi_miscible_min_monotonic)
      if (allocated(Xi_miscible_max_monotonic)) deallocate(Xi_miscible_max_monotonic)
      allocate(Xi_miscible_min_monotonic(klo:khi), source=Xi_miscible_min(klo:khi))
      allocate(Xi_miscible_max_monotonic(klo:khi), source=Xi_miscible_max(klo:khi))

      call enforce_monotonicity(Xi_miscible_min_monotonic, k_start, k_end, enforce_monotonicity_bottom_up, ierr)
      if (ierr /= 0) return

      call debugging_abundance_array_check(debug_set_dXi,&
         'adjust_miscibility_gap: error in check_abundance_array after enforce_monotonicity (min)', &
         Xi_miscible_min_monotonic, k_start, k_end, .true., monotonic_increase, ierr, Xi_miscible_min_monotonic /= -1d0)

      call enforce_monotonicity(Xi_miscible_max_monotonic, k_start, k_end, enforce_monotonicity_top_down, ierr)
      if (ierr /= 0) return

      call debugging_abundance_array_check(debug_set_dXi,&
         'adjust_miscibility_gap: error in check_abundance_array after damping Xi (max)', &
         Xi_miscible_max_monotonic, k_start, k_end, .true., monotonic_increase, ierr, Xi_miscible_max_monotonic /= -1d0)

      if (damp_sedimentation) then
         call reduce_Xi_difference(s, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, is_miscible, k_start, k_end, ierr)
         if (ierr /= 0) return

         call debugging_abundance_array_check(debug_set_dXi,&
            'adjust_miscibility_gap: error in check_abundance_array after damping Xi (min)', &
            Xi_miscible_min_monotonic, k_start, k_end, .true., monotonic_increase, ierr, Xi_miscible_min_monotonic /= -1d0)

         call debugging_abundance_array_check(debug_set_dXi,&
            'adjust_miscibility_gap: error in check_abundance_array after damping Xi (max)', &
            Xi_miscible_max_monotonic, k_start, k_end, .true., monotonic_increase, ierr, Xi_miscible_max_monotonic /= -1d0)

      end if

   end subroutine adjust_miscibility_gap

   subroutine debugging_abundance_array_check(debug,  message, Xi, klo, khi, check_range, monotonicity_method, ierr, mask)
      logical, intent(in) :: debug  ! whether to print debugging information
      character(len=*), intent(in) :: message  ! message to print
      real(kind=dp), intent(in) :: Xi(:)  ! element mass fraction
      integer, intent(in) :: klo, khi  ! region in which we check the abundance array
      logical, intent(in) :: check_range  ! whether to check the range of the abundance array
      integer, intent(in) :: monotonicity_method  ! method to check monotonicity
      integer, intent(out) :: ierr  ! error code
      logical, intent(in), optional :: mask(:)  ! mask to apply to the abundance array

      if (.not. debug) return

      call check_abundance_array(Xi, klo, khi, .true., monotonicity_method, ierr, mask)
      if (ierr /= 0) then
         write(*,*) trim(message)
         call mesa_error(__FILE__,__LINE__)
      end if


   end subroutine debugging_abundance_array_check

   logical function should_exit_loop(k, k_core, k_rain, has_core)
      ! checks whether we should exit the loop in set_dXi
      integer, intent(in) :: k  ! current index
      integer, intent(in) :: k_core  ! index of the core
      integer, intent(in) :: k_rain  ! index of the rain threshold
      logical, intent(in) :: has_core  ! whether we have a core

      should_exit_loop = .false.

      ! we reached the core. Nothing more to do.
      if (k == k_core) then
         if (debug_set_dXi) write(*,*) 'set_dXi: k_core reached, exiting loop'
         should_exit_loop = .true.
         ! If there is no core and we are below the rain threshold,
         ! we already transferred all mass fractions to the miscible value
         ! TODO: fix this check to robustly handle k_core

         ! We rained out all the material above the rain threshold. Nothing more to do.
      ! elseif (.not. has_core .and. k > k_rain) then
      !    if (debug_set_dXi) write(*,*) 'set_dXi: no core and k > k_rain, exiting loop'
      !    should_exit_loop = .true.
      end if
   end function

   subroutine init_set_dXi(s, klo, khi, is_miscible, Xi_miscible_min, Xi_miscible_max, k_RCB, k_rain, k_core, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
      ! initializes the parameters for set_dXi
      type (star_info), pointer :: s
      integer, intent(in) :: klo, khi  ! region in which we transfer the element
      logical, allocatable, intent(in) :: is_miscible(:)  ! miscible regions
      real(kind=dp), allocatable, intent(in) :: Xi_miscible_min(:)
      real(kind=dp), allocatable, intent(in) :: Xi_miscible_max(:)  ! miscible element mass fraction
      integer, intent(out) :: k_RCB  ! index of the radiative-convective boundary
      integer, intent(out) :: k_rain  ! index of the rain threshold
      integer, intent(out) :: k_core  ! index of the core
      real(kind=dp), allocatable, intent(inout) :: Xi_miscible_min_monotonic(:)  ! miscible element mass fraction, monotonically increasing
      real(kind=dp), allocatable, intent(inout) :: Xi_miscible_max_monotonic(:)  ! miscible element mass fraction, monotonically decreasing
      integer, intent(out) :: ierr  ! error code

      call find_RCB(s, klo, khi, k_RCB, ierr)
      if (ierr/=0) return

      call adjust_miscibility_gap(s, klo, khi, k_RCB, khi, is_miscible, Xi_miscible_min, Xi_miscible_max, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, ierr)
      if (ierr /= 0) then
         write(*,*) 'init_set_dXi: error in adjust_miscibility_gap'
         return
      end if

      call find_k_rain(s, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, is_miscible, k_RCB, khi, k_rain, ierr)
      if (ierr /= 0) then
         write(*,*) 'set_dXi: error in find_k_rain'
         return
      end if

      k_core = khi+1

   end subroutine init_set_dXi

   subroutine transfer_Xi(s, Xi, Xi_ref, k, klo, khi, ierr)
      ! transfers the mass fraction Xi_ref(k)-Xi(k) to the region klo:khi
      ! Although many of these parameters are present in the parent scope, we pass them to avoid confusion
      type (star_info), pointer :: s
      real(kind=dp), intent(inout) :: Xi(:)  ! element mass fraction
      real(kind=dp), intent(in) :: Xi_ref(:)  ! reference element mass fraction
      integer, intent(in) :: k  ! index of the element mass fraction to transfer
      integer, intent(in) :: klo, khi  ! region in which we transfer the element mass fraction
      integer, intent(out) :: ierr  ! error code

      ! local variables
      real(kind=dp) :: dXi    ! difference between the reference and the current element mass fraction
      real(kind=dp) :: dM_Xi  ! mass of the element that is transferred
      real(kind=dp) :: M      ! total mass in the region klo:khi

      ierr = 0 ! 0 means AOK

      ! first, compute the difference between the reference and the current element mass fraction
      dXi = Xi_ref(k) - Xi(k)

      ! next, adjust the element mass fraction in cell k
      Xi(k) = Xi_ref(k)

      ! Now, we transfer the difference to the region klo:khi
      dM_Xi = s% dm(k) * dXi     ! mass of the element that is transferred
      M = sum(s% dm(klo:khi))    ! total mass in the region klo:khi

      ! test that M is not zero
      if (M <= 0d0) then
         ierr = 1
         write(*,*) 'transfer_Xi: total mass in the region klo:khi is zero or negative'
         write(*,*) 'klo, khi: ', klo, khi
         write(*,*) 'M: ', M
         return
      end if

      ! now, we adjust the element mass fraction in the region klo:khi
      Xi(klo:khi) = Xi(klo:khi) - dM_Xi / M

   end subroutine transfer_Xi

   subroutine find_RCB(s, klo, khi, k_RCB, ierr)
      ! finds the radiative-convective boundary (RCB) in the region (klo, khi)
      type (star_info), pointer :: s
      integer, intent(in) :: klo, khi  ! region in which we search for the RCB
      integer, intent(out) :: k_RCB  ! index of the RCB
      integer, intent(out) :: ierr  ! error code

      ! local variables
      integer :: k

      ierr = 0 ! 0 means AOK
      k_RCB = -1  ! initialize to -1, meaning that we haven't found the RCB yet
      do k = klo, khi
         if ( s% mixing_type(k) /= no_mixing) then
            k_RCB = k  ! we found the RCB
            exit
         end if
      end do

      if ( k_RCB == -1 ) then
         ! no RCB found, use the entire region as the RCB
         k_RCB = klo
      end if

      if ( k_RCB < klo .or.k_RCB > khi) then
         ierr = 1
         write(*,*) 'find_RCB: RCB out of bounds'
         write(*,*) 'klo, khi, k_RCB: ', klo, khi, k_RCB
         return
      end if

      if ( debug_set_dXi ) then
         if (k_RCB > 0) then
            write(*,*) 'find_RCB: found RCB at k: ', k_RCB
         end if
      end if

   end subroutine find_RCB

   subroutine find_k_rain(s, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, is_miscible, klo, khi, k_rain, ierr)
      ! finds the last index k_rain inside (klo,khi) where the the model is immiscible outside the core
      ! this is used to determine where to deposite helium
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(in) :: Xi_miscible_min_monotonic(:)  ! miscible element mass fraction, monotonically increasing
      real(kind=dp), allocatable, intent(in) :: Xi_miscible_max_monotonic(:)  ! miscible element mass fraction, monotonically decreasing
      logical, allocatable, intent(in) :: is_miscible(:)  ! miscible regions
      integer, intent(in) :: klo, khi  ! region in which we search for the last index k_rain
      integer, intent(out) :: k_rain  ! index of the last cell where the model is immiscible outside the core
      integer, intent(out) :: ierr  ! error code

      ! local variables
      integer :: k, k_end
      logical :: inside_core, is_cell_saturated
      ierr = 0 ! 0 means AOK

      ! if there is no k_core, we have to start from the bottom
      ! first, we need to find out if we are inside the core

      inside_core = is_inside_core(s, is_miscible, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, khi)

      do k = khi, klo, -1
         is_cell_saturated = is_saturated(s, is_miscible, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, k)
         if (is_miscible(k) .and. inside_core .and. .not. is_cell_saturated) then
            inside_core = .false.  ! we are now outside the core
         else if (.not. is_miscible(k) .and. .not. inside_core) then
            k_rain = k  ! we found the last index where the model is immiscible outside the core
            if (debug) then
               write(*,*) 'find_k_rain: found k_rain: ', k_rain
               write(*,*) 'klo, khi: ', klo, khi
            end if
            return
         end if
      end do

      ! if you reach this point, it means that the model is miscible everywhere outside the core

      k_rain = klo-1  ! set k_rain to klo-1, meaning that there is no immiscible region outside the core

      if (debug_set_dXi) then
         write(*,*) 'find_k_rain: model is miscible everywhere outside the core'
         write(*,*) 'klo, khi: ', klo, khi
         write(*,*) 'k_rain set to: ', k_rain
      end if

   contains

      logical function is_inside_core(s, is_miscible, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, k)
         ! Returns .true. if cell k is inside the core
         ! there are two ways of being inside the core:
         ! 1. the last cell is immiscible. In this case, there is no core yet, but one will form
         ! 2. the cell is miscible, but greater than Xi_miscible_max. In that case, a core already exists.
         type (star_info), intent(in) :: s
         logical, intent(in) :: is_miscible(:)
         real(kind=dp), intent(in) :: Xi_miscible_min_monotonic(:)
         real(kind=dp), intent(in) :: Xi_miscible_max_monotonic(:)
         integer, intent(in) :: k

         is_inside_core = (.not. is_miscible(k)) .or. is_saturated(s, is_miscible, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, k)
      end function is_inside_core

      logical function is_saturated(s, is_miscible, Xi_miscible_min_monotonic, Xi_miscible_max_monotonic, k)
         type (star_info), intent(in) :: s
         logical, intent(in) :: is_miscible(:)
         real(kind=dp), intent(in) :: Xi_miscible_min_monotonic(:)
         real(kind=dp), intent(in) :: Xi_miscible_max_monotonic(:)
         integer, intent(in) :: k

         is_saturated = (is_miscible(k) .and. Xi_miscible_max_monotonic(k) >= 0d0 .and. Xi_miscible_max_monotonic(k) < s% xa(i_Xi, k))
      end function is_saturated


   end subroutine find_k_rain

   logical function both_radiative(s, k, j)
      type(star_info), pointer :: s
      integer, intent(in) :: k, j

      both_radiative = (s%mixing_type(k) == no_mixing) .and. (s% mixing_type(j) == no_mixing)
   end function both_radiative

   subroutine check_abundance_array(Xi, klo, khi, check_range, monotonicity_method, ierr, mask)
      real(kind=dp), intent(in) :: Xi(:)
      integer, intent(in) :: klo, khi  ! indices of Xi to be checked
      logical, intent(in) :: check_range
      integer, intent(in) :: monotonicity_method
      logical, intent(in), optional :: mask(:)
      integer, intent(out) :: ierr
      integer, allocatable :: op_err(:)

      ! local
      integer :: k
      logical, allocatable :: local_mask(:)

      ierr = 0 ! 0 means AOK

      ! Set up local mask - if mask is not present, create a mask that's all .true.
      if (present(mask)) then
         local_mask = mask
      else
         allocate(local_mask(size(Xi(klo:khi))))
         local_mask = .true.
      end if

      if (check_range .and. any((Xi(klo:khi) < 0d0 .or. Xi(klo:khi) > 1d0) .and. local_mask(klo:khi))) then
         ierr = 1
         write(*,*) 'check_abundance_array: Xi contains values outside of [0, 1]'

         if (allocated(op_err)) deallocate(op_err)
         allocate(op_err(size(Xi)), source=0)

         !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(guided)
         do k = klo, khi-1
            if ((Xi(k) < 0d0 .or. (Xi(k) > 1d0)) .and. local_mask(k)) then
               op_err(k) = 1
            end if
         end do
         !$OMP END PARALLEL DO

         write(*,*) achar(9),'k, Xi(k), /= -1d0:'
         do k = klo, khi
            if (op_err(k) /= 0) write(*,*) achar(9), k, Xi(k), Xi(k) /= -1d0, local_mask(k)
         end do

         return
      end if

      if (monotonicity_method==monotonic_increase) then

         if (allocated(op_err)) deallocate(op_err)
         allocate(op_err(size(Xi)), source=0)

         !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(guided)
         do k = klo, khi-1
            if ((Xi(k)-Xi(k+1))/Xi(k+1) > very_small .and. (local_mask(k) .and. local_mask(k+1))) then
               op_err(k) = 1
            end if
         end do
         !$OMP END PARALLEL DO

         if (any(op_err > 0)) then

            write(*,*) 'check_abundance_array: Xi is not monotonically increasing'
            write(*,*) 'k, Xi(k), Xi(k+1): '
            do k = klo, khi-1
               if (op_err(k) > 0) then
                  write(*,*) k, Xi(k), Xi(k+1)
               end if
            end do
            ierr = 1
            deallocate(op_err)
            return
         end if

      else if (monotonicity_method==monotonic_decrease) then

         if (allocated(op_err)) deallocate(op_err)
         allocate(op_err(size(Xi)), source=0)

         !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(guided)
         do k = klo, khi-1
            if ((Xi(k)-Xi(k+1))/Xi(k+1) < -very_small .and. (local_mask(k) .and. local_mask(k+1))) then
               op_err(k) = 1
            end if
         end do
         !$OMP END PARALLEL DO

         if (any(op_err > 0)) then

            write(*,*) 'check_abundance_array: Xi is not monotonically decreasing'
            write(*,*) 'k, Xi(k), Xi(k+1): '
            do k = klo, khi-1
               if (op_err(k) > 0) then
                  write(*,*) k, Xi(k), Xi(k+1)
               end if
            end do
            ierr = 1
            deallocate(op_err)
            return
         end if
      else if (monotonicity_method == no_monotonicity) then
         ! do nothing, we don't check for monotonicity
      else
         ierr = 1
         write(*,*) 'check_abundance_array: unknown monotonicity method: ', monotonicity_method
         return
      end if

      ! Clean up local mask if it was allocated
      if (allocated(local_mask)) deallocate(local_mask)
      if (allocated(op_err)) deallocate(op_err)

   end subroutine check_abundance_array

   subroutine reduce_Xi_difference(s, Xi_miscible_min, Xi_miscible_max, is_miscible, klo, khi, ierr)
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(inout) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! miscible element mass fraction
      logical, allocatable, intent(in) :: is_miscible(:)
      integer, intent(in) :: klo, khi  ! region in which we reduce the element difference
      integer, intent(out) :: ierr

      ! local variables
      integer :: k, khi_new
      real(kind=dp), allocatable :: Xi_old(:)
      real(kind=dp) :: dXi_max
      ierr = 0 ! 0 means AOK

      allocate(Xi_old(klo:khi), source=0d0)
      Xi_old = s% xa_old(i_Xi, klo:khi)

      ! define short-hand for the maximum element mass fraction change per cell (is positive)
      dXi_max = s% element_sedimentation_maximum_dXi_per_cell

      ! the is_miscible check already makes sure that Xi_miscible is less than Xi; hence, we don't need absolute values
      Xi_miscible_min(klo:khi) = merge(Xi_old(klo:khi) - dXi_max, Xi_miscible_min(klo:khi), (Xi_miscible_min(klo:khi) >= 0d0) .and. (Xi_miscible_min(klo:khi) < Xi_old(klo:khi) - dXi_max))
      Xi_miscible_max(klo:khi) = merge(Xi_old(klo:khi) + dXi_max, Xi_miscible_max(klo:khi), (Xi_miscible_max(klo:khi) >= 0d0) .and. (Xi_miscible_max(klo:khi) > Xi_old(klo:khi) + dXi_max))
      deallocate(Xi_old)

   end subroutine reduce_Xi_difference

   ! subroutine join_small_precipitating_regions(s, klo, khi, Xi_miscible_min, Xi_miscible_max, is_miscible, ierr)
   !    ! if regions are very small any only differ by a few zones, artifacts form.
   !    ! These regions are likely also not physical, because we interpolate the (uncertain) phase diagram linearly.

   !    type (star_info), pointer :: s
   !    integer, intent(in) :: klo, khi
   !    logical, allocatable, intent(inout) :: is_miscible(:)
   !    real(kind=dp), allocatable, intent(inout) :: Xi_miscible_min(:), Xi_miscible_max(:)  ! miscible element mass fraction
   !    integer, intent(out) :: ierr

   !    ! local
   !    integer :: k,j, k_start, k_end

   !    logical :: inside
   !    real(kind=dp) :: minimal_pressure_scale_height = 1d0
   !    real(kind=dp) :: dHp

   !    ierr = 0 ! 0 means AOK


   !    inside = .false.
   !    k_end = khi
   !    k_star = -1

   !    do k = khi, klo, -1
   !       if ( .not. is_miscible(k)) then
   !          k_end = k
   !       end if

   !       dHp = pressure_scale_height_distance(s, k_end, khi)
   !       if (dHp < minimal_pressure_scale_height) then

   !       endif

   !       if (debug) write(*,*) 'join_small_precipitating_regions: joining zones: ', k_start, k_end

   !    end do


   !    if (dHp < minimal_pressure_scale_height) then

   !       klo_temp(i) = klo(i)
   !       khi_temp(i) = khi(j)
   !       was_joined(i:j) = .true.
   !       exit
   !    endif


   !    ! adjust the miscible parameters accordingly
   !    do i = 1, size(khi)
   !       is_miscible(klo(i):khi(i)) = .false.
   !       Xi_miscible_min(klo(i):khi(i)) = merge(s%xa(i_Xi,klo(i):khi(i)), Xi_miscible_min(klo(i):khi(i)), Xi_miscible_min(klo(i):khi(i)) < 0d0)
   !       Xi_miscible_max(klo(i):khi(i)) = merge(s%xa(i_Xi,klo(i):khi(i)), Xi_miscible_max(klo(i):khi(i)), Xi_miscible_max(klo(i):khi(i)) < 0d0)
   !    end do

   !    deallocate(klo_temp, khi_temp, was_joined)

   ! end subroutine join_small_precipitating_regions

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

   subroutine shutdown_phase_diagram()
      type(phase_diagram) :: pd
      pd = immiscibility_phase_diagram

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

   subroutine immiscible_yi_fraction_interpolation(P0, T0, pd, yi_min, yi_max, ierr)
      real(kind=dp), intent(in) :: P0, T0
      type(phase_diagram), intent(in) :: pd
      real(kind=dp), intent(out) :: yi_min, yi_max
      real(kind=dp) :: P_lo, P_hi
      real(kind=dp) :: yi_min_lo_interp, yi_min_hi_interp, yi_max_lo_interp, yi_max_hi_interp
      real(kind=dp), allocatable :: T_lo(:), T_hi(:), yi_min_lo(:), yi_min_hi(:), yi_max_lo(:), yi_max_hi(:)
      integer, intent(out) :: ierr
      integer :: i, j
      logical :: outside_bounds
      include 'formats'
      ierr = 0 ! 0 means AOK

      ! -1d0 means out of bounds and thus by default miscible
      yi_min = -1d0
      yi_max = -1d0
      ! in principle, you should check this. However, when you paralalize the code, this
      ! will create a lot of overhead. So, I just assume that the phase diagram is correct for now
      ! call check_phase_diagram(pd, ierr)
      ! if (ierr /= 0) return
      outside_bounds = .false.
      ! find the pressure values that enclose the P0
      call bisection_locate(P0, pd%unique_P, j, outside_bounds, ierr)
      if (ierr /= 0) then
         write(*,*) 'immiscible_yi_fraction_interpolation: bisection_locate failed'
         return
      end if

      ! if we are outside the phase diagram, everything is miscible by default
      if (outside_bounds) then
         if ( debug_phase_diagram ) write(*,*) 'immiscible_yi_fraction_interpolation: P0 is outside the phase diagram. P0 = ', P0
         return
      end if

      P_lo = pd%unique_P(j)
      P_hi = pd%unique_P(j+1)

      ! now that we have the pressure values, find the respective temperature and yi_min/yi_max values
      allocate(T_lo(pd%num_Ts_at_unique_Ps(j)), T_hi(pd%num_Ts_at_unique_Ps(j+1)), &
         yi_min_lo(pd%num_Ts_at_unique_Ps(j)), yi_min_hi(pd%num_Ts_at_unique_Ps(j+1)), &
         yi_max_lo(pd%num_Ts_at_unique_Ps(j)), yi_max_hi(pd%num_Ts_at_unique_Ps(j+1)), source=-1d0)

      T_lo = pd%T(pd%unique_P_index(j):pd%unique_P_index(j+1)-1)
      yi_min_lo = pd%He_poor(pd%unique_P_index(j):pd%unique_P_index(j+1)-1)
      yi_max_lo = pd%He_rich(pd%unique_P_index(j):pd%unique_P_index(j+1)-1)

      T_hi = pd%T(pd%unique_P_index(j+1):pd%num_Ts_at_unique_Ps(j+1)+pd%unique_P_index(j+1)-1)  ! slightly cumbersome way because of the edge case where P(n-1) < P0 < P(n)
      yi_min_hi = pd%He_poor(pd%unique_P_index(j+1):pd%num_Ts_at_unique_Ps(j+1)+pd%unique_P_index(j+1)-1)
      yi_max_hi = pd%He_rich(pd%unique_P_index(j+1):pd%num_Ts_at_unique_Ps(j+1)+pd%unique_P_index(j+1)-1)

      ! check that the temperature values are strictly increasing
      do i = 1, size(T_lo)-1
         if (T_lo(i) >= T_lo(i+1)) then
            ierr = 1
            write(*,*) 'immiscible_yi_fraction_interpolation: Tj is not strictly increasing'
            deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
            return
         end if
      end do

      do i = 1, size(T_hi)-1
         if (T_hi(i) >= T_hi(i+1)) then
            ierr = 1
            write(*,*) 'immiscible_yi_fraction_interpolation: Tjp1 is not strictly increasing'
            deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
            return
         end if
      end do

      ! if T0 below the temperature range, return -1d0
      if (T0 < T_lo(1) .or. T0 < T_hi(1)) then
         if ( debug_phase_diagram ) write(*,*) 'immiscible_yi_fraction_interpolation: T0 below the data range T0 = ', T0
         deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
         return
      end if

      ! if T0 is above the temperature range, extrapolate linearly (done in lineater_list_interpolation by default)

      ! yi_min
      yi_min_lo_interp = linear_list_interpolation(T0, T_lo, yi_min_lo, outside_bounds, ierr)
      if (ierr /= 0) then
         write(*,*) 'immiscible_yi_fraction_interpolation: linear_list_interpolation for yi_min_lo failed'
         deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
         return
      end if

      ! if we are outside the bounds and don't extrapolate, return -1d0
      if (outside_bounds .and. .not. extrapolate) then
         deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
         return
      end if

      yi_min_hi_interp = linear_list_interpolation(T0, T_hi, yi_min_hi, outside_bounds, ierr)
      if (ierr /= 0) then
         write(*,*) 'immiscible_yi_fraction_interpolation: linear_list_interpolation for yi_min_hi failed'
         deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
         return
      end if
      ! if we are outside the bounds and don't extrapolate, return -1d0
      if (outside_bounds .and. .not. extrapolate) then
         deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
         return
      end if

      yi_min = linear_interpolation(P0, P_lo, P_hi, yi_min_lo_interp, yi_min_hi_interp)

      if ( extrapolate ) then
         yi_min = max(min(1d0, yi_min), 0d0)  ! make sure we are in the range of [0, 1]
      else
         if (yi_min < 0d0) then
            ierr = 1
            write(*,*) 'immiscible_yi_fraction_interpolation: Y_min_0 < 0'
            deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
            return
         end if

         if (yi_min > 1d0) then
            ierr = 1
            write(*,*) 'immiscible_yi_fraction_interpolation: Y_min_0 > 1'
            deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
            return
         end if
      end if

      ! yi_max
      yi_max_lo_interp = linear_list_interpolation(T0, T_lo, yi_max_lo, outside_bounds, ierr)
      if (ierr /= 0) then
         write(*,*) 'immiscible_yi_fraction_interpolation: linear_list_interpolation for yi_max_lo failed'
         deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
         return
      end if
      ! if we are outside the bounds and don't extrapolate, return -1d0
      if (outside_bounds .and. .not. extrapolate) then
         deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
         return
      end if

      yi_max_hi_interp = linear_list_interpolation(T0, T_hi, yi_max_hi, outside_bounds, ierr)
      if (ierr /= 0) then
         write(*,*) 'immiscible_yi_fraction_interpolation: linear_list_interpolation for yi_max_hi failed'
         deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
         return
      end if
      ! if we are outside the bounds and don't extrapolate, return -1d0
      if (outside_bounds .and. .not. extrapolate) then
         deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
         return
      end if

      ! linear interpolation for yi_max
      yi_max = linear_interpolation(P0, P_lo, P_hi, yi_max_lo_interp, yi_max_hi_interp)

      ! if extrapolation is allowed, the linear functions might go outside the range of [0, 1]
      ! in this case, we set the values to 0 and 1
      ! if extrapolation is not allowed, we need to check that the values are in the range of [0, 1]
      ! and return an error if they are not
      if (extrapolate) then
         yi_max = min(max(0d0, yi_max),1d0)
      else
         if (yi_max < 0d0) then
            ierr = 1
            write(*,*) 'immiscible_yi_fraction_interpolation: Y_max_0 < 0'
            deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
            return
         end if

         if (yi_max > 1d0) then
            ierr = 1
            write(*,*) 'immiscible_yi_fraction_interpolation: Y_max_0 > 1'
            deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
            return
         end if

         if (yi_min > yi_max) then
            ierr = 1
            write(*,*) 'immiscible_yi_fraction_interpolation: yi_min > yi_max'
            deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)
            return
         end if
      end if

      if ( debug_phase_diagram ) then
         write(*,*) '----'
         write(*,*) 'immiscible_yi_fraction_interpolation: P0(Mbar) = ', P0/1d12
         write(*,*) 'immiscible_yi_fraction_interpolation: T0 = ', T0
         write(*,*) 'immiscible_yi_fraction_interpolation: P_lo(Mbar) = ', P_lo/1d12
         write(*,*) 'immiscible_yi_fraction_interpolation: P_hi(Mbar) = ', P_hi/1d12
         write(*,*) '----'
         write(*,*) 'immiscible_yi_fraction_interpolation: T_lo = ', T_lo
         write(*,*) 'immiscible_yi_fraction_interpolation: yi_min_lo = ', yi_min_lo
         write(*,*) 'immiscible_yi_fraction_interpolation: Y_lo_interp = ', yi_min_lo_interp
         write(*,*) '----'
         write(*,*) 'immiscible_yi_fraction_interpolation: T_hi = ', T_hi
         write(*,*) 'immiscible_yi_fraction_interpolation: yi_min_hi = ', yi_min_hi
         write(*,*) 'immiscible_yi_fraction_interpolation: Y_hi_interp = ', yi_min_hi_interp
         write(*,*) '----'
      end if

      deallocate(T_lo, T_hi, yi_min_lo, yi_min_hi, yi_max_lo, yi_max_hi)

   end subroutine immiscible_yi_fraction_interpolation

   subroutine yi_is_miscible(lnP, lnT, yi, is_miscible, yi_miscible_min, yi_miscible_max, pd, ierr)
      ! gives the highest element number fraction that is miscible before the immiscibility gap
      real(kind=dp), intent(in) :: lnP, lnT, yi
      type(phase_diagram), intent(in) :: pd
      real(kind=dp), intent(out) :: yi_miscible_min, yi_miscible_max
      logical, intent(out) :: is_miscible
      real(kind=dp) :: P0, T0
      integer, intent(out) :: ierr

      ierr = 0 ! 0 means AOK

      ! convert to lin (already in cgs)
      P0 = exp(lnP)
      T0 = exp(lnT)

      call immiscible_yi_fraction_interpolation(P0, T0, pd, yi_miscible_min, yi_miscible_max, ierr)
      if (ierr /= 0) return

      ! check that yi is in the miscible range
      if (yi_miscible_min < 0d0 .or. yi_miscible_max < 0d0) then
         is_miscible = .true.
      else if (yi_miscible_max < yi_miscible_min) then  ! outside of the miscible range;

         is_miscible = .true.
         ! unphysical results; thus set to -1
         yi_miscible_min = -1d0
         yi_miscible_max = -1d0
      else if (yi <= yi_miscible_min .or. yi >= yi_miscible_max) then
         is_miscible = .true.
      else
         is_miscible = .false.
      end if

      if ( custom_debug .and. .not. is_miscible) then
         write(*,*) 'yi_is_miscible: P0 = ', P0/1d12
         write(*,*) 'yi_is_miscible: T0 = ', T0
         write(*,*) 'yi_is_miscible: yi = ', yi
         write(*,*) 'yi_is_miscible: yi_min = ', yi_miscible_min
         write(*,*) 'yi_is_miscible: yi_max = ', yi_miscible_max
      end if

   end subroutine yi_is_miscible

   subroutine adjust_mass_fractions(s,  dXi, ierr)
      type (star_info), pointer :: s
      real(kind=dp), intent(in) :: dXi(:)
      integer, intent(out) :: ierr
      integer :: i, k, n_species
      real(kind=dp) :: f
      ierr = 0 ! 0 means AOK

      n_species = size(s%xa, 1)
      ! test that dXi is the same size as s%xa
      if (size(dXi) /= s% nz) then
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

   subroutine get_scaled_number_fractions(s, xa, yi, ierr)
      use chem_def, only: chem_isos
      ! converts the mass fractions to number fractions
      ! and allowing for renormalization, in case the phase diagram is not complete
      type (star_info), pointer :: s
      real(kind=dp), intent(in) :: xa(:)
      real(kind=dp), intent(out), allocatable :: yi(:)

      integer, intent(out) :: ierr
      integer :: i, j, cid

      real(kind=dp) :: yi_scaling


      ierr = 0 ! 0 means AOK
      allocate(yi(s% species), source=0d0)

      do i = 1, s% species
         cid = s% chem_id(i)
         yi(i) = xa(i) / dble(chem_isos% Z_plus_N(cid))
      end do

      yi_scaling = 0d0
      do i = 1, size(phase_diagram_elements)
         j = s% net_iso(chem_get_iso_id(trim(phase_diagram_elements(i))))
         yi_scaling = yi_scaling + yi(j)
      end do

      yi = yi / yi_scaling ! note that this is not normalized at the moment and does thus not sum to 1

      if (any(yi < 0d0)) then
         ierr = 1
         write(*,*) 'get_scaled_number_fractions: yi contains a negative value: ', yi
         return
      end if

      if ( debug_number_fraction_scaling ) then
         write(*,*) 'get_scaled_number_fractions: xa = ', xa
         write(*,*) 'get_scaled_number_fractions: yi_scaling = ', yi_scaling
         write(*,*) 'get_scaled_number_fractions: yi'
         do i = 1, s% species
            write(*, '(i3,1pd26.16)') i, yi(i)
         end do
      end if

   end subroutine get_scaled_number_fractions

   subroutine from_scaled_number_fraction_to_mass_fraction(s, yi, yi_miscible, is_miscible,  Xi_miscible, ierr)
      ! converts the number fractions to mass fractions
      use chem_def, only: chem_isos
      type (star_info), pointer :: s
      real(kind=dp), intent(inout), allocatable :: yi(:)
      real(kind=dp), intent(in) :: yi_miscible
      logical, intent(in) :: is_miscible
      real(kind=dp), intent(out) :: Xi_miscible
      integer, intent(out) :: ierr
      integer :: i, j, cid
      real(kind=dp) :: dyi, f, mean_molecular_weight

      ierr = 0 ! 0 means AOK

      if (is_miscible .and. yi_miscible < 0d0) then
         if (debug_number_fraction_scaling) write(*,*) 'from_scaled_number_fraction_to_mass_fraction: is_miscible is true, returning without changes'
         Xi_miscible = -1d0  ! indicate that we are miscible
         return
      end if


      dyi = yi_miscible - yi(i_Xi)  ! the difference to the miscible fraction

      f = 1d0 - dyi / (1d0 - yi(i_Xi))  ! the factor to scale the number fractions

      yi = f * yi ! adjust number fractions to refelct the change in yi_miscible
      yi(i_Xi) = yi_miscible  ! set the miscible fraction to the desired value
      yi = yi / sum(yi)  ! normalize the mass fractions
      mean_molecular_weight = 0d0
      do i = 1, s% species
         cid = s% chem_id(i)
         mean_molecular_weight = mean_molecular_weight + yi(i) * dble(chem_isos% Z_plus_N(cid))
      end do

      cid = s% chem_id(i_Xi)
      Xi_miscible = yi(i_Xi) * dble(chem_isos% Z_plus_N(cid)) / mean_molecular_weight

      if ( debug_number_fraction_scaling ) then
         !$OMP CRITICAL
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: yi_miscible = ', yi_miscible
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: is_miscible = ', is_miscible
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: dyi = ', dyi
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: f = ', f
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: yi after normalization = '
         do i = 1, s% species
            write(*, '(i3, 1pd26.16)') i, yi(i)
         end do
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: mean_molecular_weight = ', mean_molecular_weight
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: Xi_miscible = ', Xi_miscible
         !$OMP END CRITICAL
      end if

      ! some tests
      if (any(yi < 0d0) .or. any(yi > 1d0)) then
         ierr = 1
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: yi is out of bounds [0, 1]: ', yi
         return
      end if
      if (abs(sum(yi) - 1d0) > tiny) then
         ierr = 1
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: yi does not sum to 1: ', sum(yi)
         return
      end if
      if (Xi_miscible < 0d0 .or. Xi_miscible > 1d0) then
         ierr = 1
         write(*,*) 'from_scaled_number_fraction_to_mass_fraction: Xi_miscible is out of bounds [0, 1]: ', Xi_miscible
         return
      end if

   end subroutine from_scaled_number_fraction_to_mass_fraction

   subroutine get_settling_velocity(s, v_settle, is_miscible, Xi_miscible, ierr)
      type (star_info), pointer :: s
      real(kind=dp), allocatable, intent(out) :: v_settle(:)
      logical, allocatable, intent(in) :: is_miscible(:)
      real(kind=dp), allocatable, intent(in) :: Xi_miscible(:)
      integer, intent(out) :: ierr

      ! local variables
      real(kind=dp) :: v_settle_0
      ierr = 0 ! 0 means AOK

      allocate(v_settle(s% nz), source=0d0)

      v_settle_0 = s% element_sedimentation_velocity
      v_settle(1:s% nz) = merge(0d0, v_settle_0, is_miscible(1:s% nz))

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

   subroutine adjust_mixing_in_sedimenting_region(s, is_miscible, ierr)
      type(star_info), pointer :: s
      logical, intent(in) :: is_miscible(:)
      integer, intent(out) :: ierr

      ! local variables
      integer :: k
      integer :: klo, khi  ! lower and upper bounds of the sedimenting region that we exclude from mixing
      real(kind=dp) :: m_min, m_max  ! lower and upper bounds in mass of the sedimenting region that we exclude from mixing
      logical :: boundaries_were_set = .false.
      save :: m_min, m_max, boundaries_were_set
      ierr = 0 ! 0 means AOK

      if ( s% do_conv_premix ) then
         write(*,*) 'You are using convective premixing and instant helium rain simultaneously. This will likely cause crashes.'
         write(*,*) 'In CPM, regions that rain out will flip-flop between homogeneous and rainout mixing.'
         write(*,*) 'This check sets do_conv_premix = .false.'
         write(*,*) 'If you know what you are doing, feel free to remove this check.'
         s% do_conv_premix = .false.
      end if

      ! let's try the easiest solution first: all cells that are not miscible will not be mixed
      s% D_mix(1:s% nz) = merge(s% D_mix(1:s% nz), 0d0, is_miscible(1:s% nz))
      ! idea: I should check if a cell is immiscible or--alternatively--if it is in a precipiating region.

      ! if the region is not set yet, initialize it with the first region
      ! if ( .not. boundaries_were_set ) then
      !    if (debug) write(*,*) 'adjust_mixing_in_sedimenting_region: set initial boundaries'

      !     klo = -1
      !     khi = -1
      !     do k = 1, size(is_miscible)
      !       if (.not. is_miscible(k)) then
      !          if (klo == -1) klo = k
      !          khi = i
      !       end if
      !     end do
      !     if (klo == -1 .or. khi == -1) then
      !       ierr = 1
      !       write(*,*) 'adjust_mixing_in_sedimenting_region: no immiscible region found in is_miscible'
      !       return
      !     end if

      !    if (verbose) write(*,*) 'initializing klo, khi = ', klo, khi

      !    ! ? Is this the best way to determine the non-mixing region? Just take the largest possible one? I should probably introduce multiple regions.
      !    ! ? For now, we typically only have one region anyways, so this is fine. But for more complex cases, we might want to generalize this.
      !    if (size(precipitating_regions) > 1) then
      !       if (verbose) write(*,*) 'zone_boundaries_by_mass: multiple regions found, adjusting boundaries'
      !       do i = 2, size(precipitating_regions)
      !          klo = min(precipitating_regions(i)% klo, klo)
      !          khi = max(precipitating_regions(i)% khi, khi)
      !       end do
      !    end if

      !    m_max = s%m(klo)
      !    m_min = s%m(khi)
      !    boundaries_were_set = .true.
      ! else

      !    ! find klo and khi, the indices closest to m_max and m_min, respectively
      !    khi = minloc(abs(s%m(1:s%nz)-m_min), dim=1)
      !    klo = minloc(abs(s%m(1:s%nz)-m_max), dim=1)

      !    ! check if the regions expanded
      !    do i = 1, size(precipitating_regions)
      !       klo = min(klo, precipitating_regions(i)% klo)
      !       khi = max(khi, precipitating_regions(i)% khi)
      !    end do
      !    m_max = s%m(klo)
      !    m_min = s%m(khi)
      ! end if

      ! s% D_mix(klo:khi) = 0d0  ! turn everything off for now; HHK

      if (verbose) write(*,'(a,2i12,3f26.16)') 'zone_boundaries_by_mass: klo, khi, m_max, m_min, dm = ', klo, khi, m_min/m_jupiter, m_max/m_jupiter, (m_max-m_min)/m_jupiter

   end subroutine adjust_mixing_in_sedimenting_region

   ! copied from phase_separation.F90
   subroutine update_model_ (s, kc_t, kc_b, do_brunt)

      use turb_info, only: set_mlt_vars
      use brunt, only: do_brunt_B
      use micro

      type(star_info), pointer :: s
      integer, intent(in)      :: kc_t
      integer, intent(in)      :: kc_b
      logical, intent(in)      :: do_brunt

      integer  :: ierr
      integer  :: kf_t
      integer  :: kf_b

      logical :: mask(s% nz)

      mask(:) = .true.

      ! Update the model to reflect changes in the abundances across
      ! cells kc_t:kc_b (the mask part of this call is unused, mask=true for all zones).
      ! Do updates at constant (P,T) rather than constant (rho,T).
      ! TODO: Add another mode that does not hold T constant
      s%fix_Pgas = .true.
      call set_eos_with_mask(s, kc_t, kc_b, mask, ierr)
      if (ierr /= 0) then
         write(*,*) 'element_sedimentation: error from call to set_eos_with_mask'
         stop
      end if
      s%fix_Pgas = .false.

      ! Update opacities across cells kc_t:kc_b (this also sets rho_face
      ! and related quantities on faces kc_t:kc_b)
      call set_micro_vars(s, kc_t, kc_b, &
         skip_eos=.TRUE., skip_net=.TRUE., skip_neu=.TRUE., skip_kap=.FALSE., ierr=ierr)
      if (ierr /= 0) then
         write(*,*) 'element_sedimentation: error from call to set_micro_vars'
         stop
      end if

      ! This is expensive, so only do it if we really need to.
      if(do_brunt) then
         ! Need to make sure we can set brunt for mix_outward calculation.
         if(.not. s% calculate_Brunt_B) then
            stop "phase separation requires s% calculate_Brunt_B = .true."
         end if
         call do_brunt_B(s, kc_t, kc_b, ierr) ! for unsmoothed_brunt_B
         if (ierr /= 0) then
            write(*,*) 'element_sedimentation: error from call to do_brunt_B'
            stop
         end if
      end if

      ! Finally update MLT for interior faces

      kf_t = kc_t
      kf_b = kc_b + 1

      call set_mlt_vars(s, kf_t+1, kf_b-1, ierr)
      if (ierr /= 0) then
         write(*,*) 'element_sedimentation: failed in call to set_mlt_vars during update_model_'
         stop
      endif

      ! Finish

      return

   end subroutine update_model_


end module element_sedimentation
