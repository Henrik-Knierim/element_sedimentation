! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_star_extras

   implicit none

   ! Globale Variable für das Regenmodell
   integer, parameter :: helium_rain_method = instant_rain
   real(kind=dp), parameter :: maximum_dY_change_per_step = -1d-2

   character(len=*), parameter :: phase_diagram_filename = 'src/Schötter+2018.csv'

contains

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

   
      s% extras_after_evolve => extras_after_evolve
 
      ! Helium Rain
      call init_helium_rain(phase_diagram_filename, maximum_dY_change_per_step, ierr)
      s% other_diffusion => my_other_diffusion
      s% other_D_mix => my_other_D_mix

   end subroutine extras_controls


   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      call free_phase_diagram(helium_immiscibility_phase_diagram)
   end subroutine extras_after_evolve

   subroutine my_other_diffusion(id, dt, ierr)
      use star_def
      integer, intent(in) :: id
      real(dp), intent(in) :: dt
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      real(dp), allocatable :: y(:)
      real(dp) :: logP_lo, logP_hi, dy
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      call do_helium_rain(s, helium_rain_method, ierr)
      if (ierr /= 0) then
         write(*,*) 'my_other_diffusion: error in do_helium_rain'
         return
      end if

   end subroutine my_other_diffusion

   subroutine determine_no_mixing_zone(s,  precipitating_regions, klo, khi, ierr)
      ! for now, assume only one zone
      type(star_info), pointer :: s
      type(immiscibility_region), allocatable, intent(in) :: precipitating_regions(:)
      integer, intent(out) :: ierr

      integer :: klo, khi
      integer :: klo_local, khi_local
      integer :: i

      logical :: first_call = .true.

      real(kind=dp) :: m_min, m_max
      save :: m_min, m_max, first_call
      logical :: verbose = .false.

      ierr = 0 ! 0 means AOK

      ! Check if precipitating_regions is empty
      if (size(precipitating_regions) == 0) then
         ierr = 1
         write(*,*) 'zone_boundaries_by_mass: no precipitating regions provided'
         return
      end if

      if ( size(precipitating_regions) > 1) then
         ierr = 1
         write(*,*) 'zone_boundaries_by_mass: more than one precipitating region'
         return
      end if

      if ( first_call ) then
         if (verbose) write(*,*) 'zone_boundaries_by_mass: first call'
         klo = precipitating_regions(1)% klo
         khi = precipitating_regions(1)% khi
         ! do i = 1, size(precipitating_regions)
         ! you need to initalize klo and khi first before taking the min and max.
         !    klo = min(precipitating_regions(i)% klo, klo)
         !    khi = max(precipitating_regions(i)% khi, khi)
         ! end do
         m_max = s%m(klo)
         m_min = s%m(khi)
         first_call = .false.
         return
      end if

      ! find klo and khi, the indices closest to m_max and m_min, respectively
      khi = minloc(abs(s%m(1:s%nz)-m_min), dim=1)
      klo = minloc(abs(s%m(1:s%nz)-m_max), dim=1)

      do i = 1, size(precipitating_regions)
         klo_local = precipitating_regions(i)% klo
         khi_local = precipitating_regions(i)% khi

         klo = min(klo, klo_local)
         khi = max(khi, khi_local)

      end do

      m_max = s%m(klo)
      m_min = s%m(khi)

      if (verbose) write(*,'(a,2i12,3f26.16)') 'zone_boundaries_by_mass: klo, khi, m_max, m_min, dm = ', klo, khi, m_min/m_jupiter, m_max/m_jupiter, (m_max-m_min)/m_jupiter

   end subroutine determine_no_mixing_zone

   subroutine my_other_D_mix(id, ierr)
      use star_def
      integer, intent(out) :: ierr
      integer, intent(in) :: id
      type (star_info), pointer :: s
      type(immiscibility_region), allocatable :: precipitating_regions(:)


      integer :: i, klo, khi
      logical :: verbose = .false.
      logical :: boundaries_were_set = .false.
      save :: klo, khi, boundaries_were_set
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! if the helium rain mode is instant_rain, we set the D_mix to zero in the precipitation zone
      select case (helium_rain_method)
       case (instant_rain)
         call determine_precipitating_regions(s, precipitating_regions, ierr)
         if (ierr /= 0) return

         if (.not. allocated(precipitating_regions)) then
            if ( verbose ) write(*,*) 'my_other_D_mix: no precipitating regions found'
            if (.not. boundaries_were_set) return
            if ( verbose ) write(*,*) 'my_other_D_mix: take old boundaries.'
         else
            if ( verbose ) write(*,*) 'my_other_D_mix: call determine_no_mixing_zone'
            call determine_no_mixing_zone(s, precipitating_regions, klo, khi, ierr)
            if (ierr /= 0) return
            boundaries_were_set = .true.
            deallocate(precipitating_regions)
         end if
         if ( verbose ) write(*,'(a,2i12)') 'my_other_D_mix: klo, khi', klo, khi

         if ( s%do_conv_premix ) then
            write(*,*) 'You are using convective premixing and instant helium rain simultaneously. This will likely cause crashes.'
            write(*,*) 'In CPM, regions that rain out will flip-flop between homogeneous and rainout mixing.'
            write(*,*) 'If you know what you are doing, feel free to remove this check.'
            ierr = 1
            return
         end if

         s% D_mix(klo:khi) = 0d0
       case (constant_velocity_rain)
         ! call do_constant_velocity_rain(s, ierr)
       case (advenction_equation_rain)
         ! call do_advection_equation_rain(s, ierr)
       case default
         ierr = 1
         write(*,*) 'my_other_D_mix: unknown method for helium rain'
         return
      end select

   end subroutine my_other_D_mix

end module run_star_extras


