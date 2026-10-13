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

module test_sedimentation_energy

   ! unit tests for sedimentation_energy: flow rates from a composition change,
   ! energy flux and its divergence, and the partial enthalpy at fixed (T, P)
   ! from (rho, T) derivatives checked against ideal gases of H2 and He.

   use const_def, only: dp, kerg, mp
   use sedimentation_energy
   use test_support, only: assert_true, assert_close

   implicit none

   private
   public :: run_sedimentation_energy_tests

   real(dp), parameter :: rtol = 1d-12, atol = 1d-300

contains

   subroutine run_sedimentation_energy_tests(nfail)
      integer, intent(inout) :: nfail
      write(*,'(a)') 'test_sedimentation_energy'
      call test_flow_rates(nfail)
      call test_flux_and_divergence(nfail)
      call test_partial_enthalpy_ideal_gas(nfail)
      call test_partial_enthalpy_nonideal(nfail)
   end subroutine run_sedimentation_energy_tests


   subroutine test_flow_rates(nfail)
      ! three cells, three species: He sinks from cell 1 to cell 3, replaced by H and Z
      ! in proportion to their abundances; check the face flow rates against the
      ! masses moved and the zero sums
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 3, species = 3
      real(dp) :: dm(nz), xa_old(species,nz), xa_new(species,nz), flow(species,nz+1), dt
      real(dp) :: lost_he
      integer :: j, k

      dm = [2d0, 1d0, 4d0]
      dt = 10d0
      ! species 1 = H, 2 = He, 3 = Z
      xa_old(:,1) = [0.70d0, 0.28d0, 0.02d0]
      xa_old(:,2) = [0.70d0, 0.28d0, 0.02d0]
      xa_old(:,3) = [0.60d0, 0.38d0, 0.02d0]
      xa_new = xa_old
      ! cell 1 loses 0.04 of He (mass 0.08), cell 2 unchanged, cell 3 gains it (0.08/4 = 0.02)
      xa_new(2,1) = 0.24d0
      xa_new(2,3) = 0.40d0
      ! the counter-flow fills the donor with H and Z in the proportions 0.7:0.02 and takes
      ! the same masses from the receiver (every species conserves its total mass)
      xa_new(1,1) = xa_old(1,1) + 0.04d0*0.70d0/0.72d0
      xa_new(3,1) = xa_old(3,1) + 0.04d0*0.02d0/0.72d0
      xa_new(1,3) = xa_old(1,3) - 0.02d0*0.70d0/0.72d0
      xa_new(3,3) = xa_old(3,3) - 0.02d0*0.02d0/0.72d0

      call species_face_flow_rates(nz, species, dm, xa_old, xa_new, dt, flow)

      lost_he = dm(1)*0.04d0  ! 0.08 g of He crossed faces 2 and 3 inward
      call assert_close('flow He face 2 (inward, negative)', flow(2,2), -lost_he/dt, rtol, atol, nfail)
      call assert_close('flow He face 3 (inward, negative)', flow(2,3), -lost_he/dt, rtol, atol, nfail)
      call assert_close('flow H face 2 (outward)', flow(1,2), lost_he/dt*0.70d0/0.72d0, rtol, atol, nfail)
      call assert_close('flow Z face 2 (outward)', flow(3,2), lost_he/dt*0.02d0/0.72d0, rtol, atol, nfail)
      ! what crosses face 3 is what cell 3 gained
      call assert_close('flow H face 3 from the receiver', flow(1,3), dm(3)*0.02d0*0.70d0/0.72d0/dt, rtol, atol, nfail)
      call assert_true('no flow through the outer boundary', all(flow(:,1) == 0d0), nfail)
      call assert_true('no flow through the inner boundary', all(flow(:,nz+1) == 0d0), nfail)
      do k = 2, nz
         call assert_close('flow rates of a face add up to zero', sum(flow(:,k)), 0d0, 0d0, 1d-15, nfail)
      end do
      ! cell balance: dm*(xa_new - xa_old)/dt = flow(k+1) - flow(k)
      do k = 1, nz
         do j = 1, species
            call assert_close('cell balance', dm(k)*(xa_new(j,k) - xa_old(j,k))/dt, &
               flow(j,k+1) - flow(j,k), 1d-12, 1d-15, nfail)
         end do
      end do
   end subroutine test_flow_rates


   subroutine test_flux_and_divergence(nfail)
      ! the energy flux is independent of a common offset of the partial enthalpies and
      ! its divergence conserves energy
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 4, species = 2
      real(dp) :: dm(nz), flow(species,nz+1), h(species,nz+1), h_shifted(species,nz+1), &
         energy_flux(nz+1), energy_flux_shifted(nz+1), eps(nz)
      integer :: k

      dm = [1d0, 2d0, 3d0, 4d0]
      flow = 0d0
      do k = 2, nz
         flow(2,k) = -1d3*dble(k)       ! He inward
         flow(1,k) = -flow(2,k)         ! H outward
      end do
      do k = 1, nz + 1
         h(1,k) = 1d12 + 1d10*dble(k)   ! H, hotter inward
         h(2,k) = 3d11 + 3d9*dble(k)    ! He
      end do
      call face_energy_flux(nz, species, flow, h, energy_flux)
      call assert_close('energy flux face 2', energy_flux(2), 2d3*(h(1,2) - h(2,2)), rtol, atol, nfail)
      call assert_close('energy flux face 4', energy_flux(4), 4d3*(h(1,4) - h(2,4)), rtol, atol, nfail)
      call assert_true('outward energy flux when hot H replaces cold He', all(energy_flux(2:nz) > 0d0), nfail)

      h_shifted = h + 7d11
      call face_energy_flux(nz, species, flow, h_shifted, energy_flux_shifted)
      do k = 1, nz + 1
         call assert_close('flux independent of the enthalpy reference', energy_flux_shifted(k), &
            energy_flux(k), 1d-12, 1d-6, nfail)
      end do

      call flux_divergence_heating(nz, dm, energy_flux, eps)
      call assert_close('energy conserved: sum(dm eps) = 0', dot_product(dm, eps), 0d0, 0d0, &
         1d-12*maxval(abs(energy_flux)), nfail)
      call assert_close('cell 1 gains the flux through face 2', eps(1), energy_flux(2)/dm(1), rtol, atol, nfail)
      call assert_close('cell 4 loses the flux through face 4', eps(4), -energy_flux(4)/dm(4), rtol, atol, nfail)
      call assert_close('cell 2 net', eps(2), (energy_flux(3) - energy_flux(2))/dm(2), rtol, atol, nfail)
   end subroutine test_flux_and_divergence


   subroutine test_partial_enthalpy_ideal_gas(nfail)
      ! mixture of ideal gases H2 (5/2 kT per molecule, mass 2 m_p) and He (3/2 kT, 4 m_p):
      ! replacing He by H2 at fixed (T, P) changes h by 7/4 - 5/8 = 9/8 kT/m_p per gram,
      ! while the (rho, T) differences alone give 5/4 - 3/8 = 7/8 for e; the pressure
      ! correction supplies the remaining 1/4
      integer, intent(inout) :: nfail
      real(dp), parameter :: T = 1d4, rho = 1d-1, delta = 1d-3
      real(dp) :: x, de, dPgas, dh, e_plus, e_minus, P_plus, P_minus, kT_mp

      kT_mp = kerg*T/mp
      x = 0.7d0  ! H2 mass fraction, the rest He
      e_plus = e_ideal(x + delta, T)
      e_minus = e_ideal(x - delta, T)
      P_plus = P_ideal(x + delta, rho, T)
      P_minus = P_ideal(x - delta, rho, T)
      de = (e_plus - e_minus)/(2d0*delta)
      dPgas = (P_plus - P_minus)/(2d0*delta)
      call assert_close('ideal gas: de/dX at fixed rho,T = 7/8 kT/m_p', de, 0.875d0*kT_mp, 1d-10, atol, nfail)
      call assert_close('ideal gas: dP/dX/rho at fixed rho,T = 1/4 kT/m_p', dPgas/rho, 0.25d0*kT_mp, 1d-10, atol, nfail)
      dh = partial_enthalpy_increment(de, dPgas, rho, 1d0, 1d0)
      call assert_close('ideal gas: dh/dX at fixed T,P = 9/8 kT/m_p', dh, 1.125d0*kT_mp, 1d-10, atol, nfail)
      ! the sign: hydrogen carries more enthalpy per gram than helium
      call assert_true('hydrogen carries more enthalpy per gram', dh > 0d0, nfail)
   contains
      real(dp) function e_ideal(xh2, T) result(e)
         real(dp), intent(in) :: xh2, T
         e = kerg*T/mp*(xh2*1.25d0 + (1d0 - xh2)*0.375d0)
      end function e_ideal
      real(dp) function P_ideal(xh2, rho, T) result(P)
         real(dp), intent(in) :: xh2, rho, T
         P = rho*kerg*T/mp*(xh2/2d0 + (1d0 - xh2)/4d0)
      end function P_ideal
   end subroutine test_partial_enthalpy_ideal_gas


   subroutine test_partial_enthalpy_nonideal(nfail)
      ! a gas with a composition-dependent van der Waals-like correction: the formula
      ! must agree with a direct finite difference of h at fixed (T, P), where the
      ! density is solved for at the perturbed composition
      integer, intent(inout) :: nfail
      real(dp), parameter :: T = 2d4, rho0 = 0.3d0, x0 = 0.6d0, delta = 1d-4
      real(dp) :: P0, de, dPgas, chiT, chiRho, dh_formula, dh_direct, rho_plus, rho_minus, h_plus, h_minus

      P0 = P_of(x0, rho0, T)
      de = (e_of(x0 + delta, rho0, T) - e_of(x0 - delta, rho0, T))/(2d0*delta)
      dPgas = (P_of(x0 + delta, rho0, T) - P_of(x0 - delta, rho0, T))/(2d0*delta)
      chiT = T/P0*(P_of(x0, rho0, T*(1d0 + 1d-6)) - P_of(x0, rho0, T*(1d0 - 1d-6)))/(2d-6*T)
      chiRho = rho0/P0*(P_of(x0, rho0*(1d0 + 1d-6), T) - P_of(x0, rho0*(1d0 - 1d-6), T))/(2d-6*rho0)
      dh_formula = partial_enthalpy_increment(de, dPgas, rho0, chiT, chiRho)

      rho_plus = rho_at_P(x0 + delta, P0, T)
      rho_minus = rho_at_P(x0 - delta, P0, T)
      h_plus = e_of(x0 + delta, rho_plus, T) + P0/rho_plus
      h_minus = e_of(x0 - delta, rho_minus, T) + P0/rho_minus
      dh_direct = (h_plus - h_minus)/(2d0*delta)
      call assert_close('non-ideal gas: formula = direct (T,P) difference', dh_formula, dh_direct, 1d-8, atol, nfail)
      call assert_true('non-ideal gas: pressure correction matters', abs(dh_formula - de) > 1d-3*abs(dh_formula), nfail)
   contains
      real(dp) function e_of(x, rho, T) result(e)
         real(dp), intent(in) :: x, rho, T
         ! ideal part plus an attraction term -a(x) rho
         e = kerg*T/mp*(x*1.25d0 + (1d0 - x)*0.375d0) - a_of(x)*rho
      end function e_of
      real(dp) function P_of(x, rho, T) result(P)
         real(dp), intent(in) :: x, rho, T
         real(dp) :: n
         n = rho/mp*(x/2d0 + (1d0 - x)/4d0)
         P = n*kerg*T/(1d0 - b_of(x)*rho) - a_of(x)*rho*rho
      end function P_of
      real(dp) function a_of(x)
         real(dp), intent(in) :: x
         a_of = 1d10*(1d0 + 0.5d0*x)
      end function a_of
      real(dp) function b_of(x)
         real(dp), intent(in) :: x
         b_of = 0.2d0*(1d0 + x)
      end function b_of
      real(dp) function rho_at_P(x, P, T) result(rho)
         ! Newton iteration for the density at given (x, P, T)
         real(dp), intent(in) :: x, P, T
         real(dp) :: f, dfdrho, r
         integer :: it
         rho = rho0
         do it = 1, 50
            f = P_of(x, rho, T) - P
            r = rho*1d-7
            dfdrho = (P_of(x, rho + r, T) - P_of(x, rho - r, T))/(2d0*r)
            rho = rho - f/dfdrho
            if (abs(f) < 1d-13*P) exit
         end do
      end function rho_at_P
   end subroutine test_partial_enthalpy_nonideal

end module test_sedimentation_energy
