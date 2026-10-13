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

module sedimentation_energy

   ! Energy carried by the species that the sedimentation moves between cells.
   !
   ! In the Lagrangian (barycentric) frame every cell keeps its mass, so the
   ! species fluxes through a face add up to zero: the sinking species is
   ! replaced by an equal mass of the others. Gravity accelerates all species
   ! alike, hence the relative motion does no net gravitational work in this
   ! frame and the gravitational energy released by the rain enters the energy
   ! equation through the contraction of the mass shells (P dV work in MESA's
   ! dedt form, -T ds/dt in the eps_grav form). What the relative motion does
   ! carry is enthalpy: the energy flux through a face is
   !
   !    E = sum_i hbar_i J_i,
   !
   ! with J_i the mass flow rate of species i and hbar_i its partial specific
   ! enthalpy at the temperature and pressure of the face (de Groot & Mazur,
   ! Non-Equilibrium Thermodynamics, ch. III; Bird, Stewart & Lightfoot,
   ! Transport Phenomena, eq. 19.3-3). Its divergence heats a cell by the
   ! enthalpy difference of what arrives and what leaves. It contains the
   ! buoyancy work g sum_i vbar_i J_i of the sinking species (the gravitational
   ! energy dissipated by drag where the species falls, vbar_i the partial
   ! specific volumes) and the heat carried by the exchange of hot material
   ! from below against cold material from above. Summed over the star the
   ! term vanishes: it redistributes energy and conserves the total.
   !
   ! The partial enthalpy differences follow from the equation of state by a
   ! composition derivative at fixed (T, P), which the (rho, T) derivatives
   ! MESA's EOS returns give as
   !
   !    (dh/dX)_{T,P} = (de/dX)_{rho,T} + (chiT/chiRho) (dP/dX)_{rho,T} / rho
   !
   ! (from h = e + P/rho and the Maxwell relation (de/dv)_T = T (dP/dT)_v - P).
   ! For an ideal gas of H2 and He this gives hbar_H2 - hbar_He = 9/8 kT/m_p.
   !
   ! Everything here is plain arithmetic on arrays (no MESA star_info), so it
   ! can be unit tested; the driver supplies the EOS evaluations.

   use const_def, only: dp

   implicit none

   private
   public :: species_face_flow_rates, face_energy_flux, flux_divergence_heating, &
      partial_enthalpy_increment

contains

   subroutine species_face_flow_rates(nz, species, dm, xa_old, xa_new, dt, flow)
      ! Net mass flow rate of every species through every face [g/s], positive
      ! outward, from the composition change of a step: the mass of species j that
      ! left the cells above face k (1..k-1) must have crossed face k. Face k is the
      ! outer face of cell k; faces 1 and nz+1 carry no flow. The rates of all
      ! species add up to zero at every face when the cells keep their mass
      ! (sum of xa = 1 before and after).
      integer, intent(in) :: nz, species
      real(dp), intent(in) :: dm(:)          ! (nz) cell masses [g]
      real(dp), intent(in) :: xa_old(:,:)    ! (species, nz) mass fractions before
      real(dp), intent(in) :: xa_new(:,:)    ! (species, nz) mass fractions after
      real(dp), intent(in) :: dt             ! length of the step [s]
      real(dp), intent(out) :: flow(:,:)     ! (species, nz+1) [g/s]
      integer :: j, k
      real(dp) :: lost
      do j = 1, species
         flow(j,1) = 0d0
         lost = 0d0
         do k = 2, nz + 1
            ! mass of species j that the cells 1..k-1 lost during the step
            lost = lost + dm(k-1)*(xa_old(j,k-1) - xa_new(j,k-1))
            flow(j,k) = -lost/dt
         end do
      end do
      flow(1:species,nz+1) = 0d0  ! exact when the species conserves its mass; drop the round-off
   end subroutine species_face_flow_rates


   subroutine face_energy_flux(nz, species, flow, h_partial, energy_flux)
      ! Energy flow rate through every face [erg/s], positive outward:
      ! E(k) = sum_j hbar_j(k) J_j(k) with the partial specific enthalpies of the
      ! species at the face. Because the flow rates add up to zero, the hbar_j may
      ! be given relative to any common reference (e.g. the partner species).
      integer, intent(in) :: nz, species
      real(dp), intent(in) :: flow(:,:)        ! (species, nz+1) [g/s]
      real(dp), intent(in) :: h_partial(:,:)   ! (species, nz+1) [erg/g]
      real(dp), intent(out) :: energy_flux(:)  ! (nz+1) [erg/s]
      integer :: k
      do k = 1, nz + 1
         energy_flux(k) = dot_product(h_partial(1:species,k), flow(1:species,k))
      end do
   end subroutine face_energy_flux


   subroutine flux_divergence_heating(nz, dm, energy_flux, eps)
      ! Heating rate per unit mass [erg/g/s] of every cell from the energy flux
      ! through its faces: cell k gains what enters through its inner face k+1 and
      ! loses what leaves through its outer face k. Sum(dm*eps) telescopes to
      ! E(nz+1) - E(1), i.e. zero for closed boundaries.
      integer, intent(in) :: nz
      real(dp), intent(in) :: dm(:)           ! (nz) [g]
      real(dp), intent(in) :: energy_flux(:)  ! (nz+1) [erg/s]
      real(dp), intent(out) :: eps(:)         ! (nz) [erg/g/s]
      integer :: k
      do k = 1, nz
         eps(k) = (energy_flux(k+1) - energy_flux(k))/dm(k)
      end do
   end subroutine flux_divergence_heating


   real(dp) function partial_enthalpy_increment(de, dPgas, rho, chiT, chiRho) result(dh)
      ! Change of the specific enthalpy at fixed (T, P) along a composition
      ! direction, from the changes of the specific internal energy de [erg/g] and
      ! of the gas pressure dPgas [dyn/cm^2] along the same direction at fixed
      ! (rho, T), the density rho [g/cm^3] and the logarithmic derivatives
      ! chiT = (dlnP/dlnT)_rho, chiRho = (dlnP/dlnrho)_T of the total pressure:
      !    (dh/dX)_{T,P} = (de/dX)_{rho,T} + (chiT/chiRho) (dP/dX)_{rho,T} / rho.
      ! Radiation pressure does not depend on the composition, so dP = dPgas.
      real(dp), intent(in) :: de, dPgas, rho, chiT, chiRho
      dh = de + (chiT/chiRho)*dPgas/rho
   end function partial_enthalpy_increment

end module sedimentation_energy
