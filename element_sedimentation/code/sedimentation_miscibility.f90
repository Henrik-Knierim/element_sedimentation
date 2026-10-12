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

module sedimentation_miscibility

   ! Composition arithmetic for the sedimenting species and per-cell
   ! evaluation of the immiscibility limits. No star_info, no chem_def:
   ! atomic weights are passed in.
   !
   ! The phase diagram is tabulated in the number fraction y of the
   ! sedimenting species i within the two-component (i, partner) subsystem,
   !    y = n_i/(n_i + n_p),  n_j = xa(j)/A(j).
   ! Other species are carried along without entering the phase diagram.

   use const_def, only: dp
   use utils_lib, only: is_bad
   use sedimentation_phase_diagram, only: phase_diagram, phase_diagram_limits, &
      phase_diagram_err_not_loaded, phase_diagram_error_message

   implicit none

   private
   public :: pair_number_fraction, mass_fraction_for_pair_number_fraction, evaluate_miscibility
   public :: miscibility_error_message

   ! error codes (ierr = 0 means OK). Codes of sedimentation_phase_diagram
   ! (1..5) are passed through unchanged by evaluate_miscibility.
   integer, parameter, public :: miscibility_ok = 0
   integer, parameter, public :: miscibility_err_bad_input = 11    ! sizes, indices, atomic weights
   integer, parameter, public :: miscibility_err_composition = 12  ! NaN or negative mass fraction in a cell

contains

   pure real(dp) function pair_number_fraction(xa, A, i_sed, i_partner) result(y)
      ! y = (xa(i)/A(i)) / (xa(i)/A(i) + xa(p)/A(p)); 0 if both are 0
      real(dp), intent(in) :: xa(:)  ! mass fractions
      real(dp), intent(in) :: A(:)   ! atomic weights
      integer, intent(in) :: i_sed, i_partner
      real(dp) :: n_sed, n_partner
      n_sed = xa(i_sed)/A(i_sed)
      n_partner = xa(i_partner)/A(i_partner)
      if (n_sed + n_partner > 0d0) then
         y = n_sed/(n_sed + n_partner)
      else
         y = 0d0
      end if
   end function pair_number_fraction


   pure real(dp) function mass_fraction_for_pair_number_fraction(xa, A, i_sed, i_partner, y_target) result(x)
      ! The mass fraction of species i_sed that the cell would have if its pair
      ! number fraction were y_target, keeping the number ratios of all other
      ! species to the partner fixed.
      !
      ! Derivation: with n_j = xa(j)/A(j), scale all n_j (j /= i) by a common
      ! factor f and choose n_i such that n_i/(n_i + f n_p) = y_target, i.e.
      ! n_i = y_target f n_p/(1 - y_target). Then
      !    X_i = n_i A_i / (n_i A_i + f sum_{j/=i} n_j A_j)
      !        = y_target A_i n_p / (y_target A_i n_p + (1 - y_target) S),
      ! with S = sum_{j/=i} xa(j); f cancels. The old implementation used
      ! f = (1 - y_target)/(1 - y_current) and renormalized, which gives the
      ! same X_i but divides by (1 - y_current); this form has no such division.
      !
      ! Degenerate case (denominator 0, i.e. y_target = 0 and S = 0, or
      ! xa(partner) = 0 and (y_target = 1 or S = 0)): for y_target = 0 the
      ! answer is 0 (X_i = 0 always gives y = 0); otherwise no composition
      ! change can reach y_target and the current normalized xa(i_sed) is
      ! returned (so that in_gap is false for such a cell).
      !
      ! pair_number_fraction of the composition (xa(j) scaled by (1 - X_i)/S
      ! for j /= i, xa(i) = X_i) equals y_target to round-off.
      real(dp), intent(in) :: xa(:)  ! mass fractions
      real(dp), intent(in) :: A(:)   ! atomic weights
      integer, intent(in) :: i_sed, i_partner
      real(dp), intent(in) :: y_target
      real(dp) :: n_partner, others, numerator, denominator, total
      integer :: j

      n_partner = xa(i_partner)/A(i_partner)
      others = 0d0
      do j = 1, size(xa)
         if (j /= i_sed) others = others + xa(j)
      end do

      numerator = y_target*A(i_sed)*n_partner
      denominator = numerator + (1d0 - y_target)*others

      if (denominator > 0d0) then
         x = numerator/denominator
      else if (y_target <= 0d0) then
         x = 0d0
      else
         total = others + xa(i_sed)
         if (total > 0d0) then
            x = xa(i_sed)/total
         else
            x = 0d0
         end if
      end if

   end function mass_fraction_for_pair_number_fraction


   subroutine evaluate_miscibility(pd, nz, species, xa, A, i_sed, i_partner, P, T, &
         has_gap, in_gap, x_lo, x_hi, ierr)
      ! Per-cell miscibility of species i_sed.
      !   has_gap(k)      the phase diagram has an immiscibility gap at (P(k), T(k))
      !   x_lo(k), x_hi(k) mass-fraction limits of species i_sed in cell k
      !                   (valid only where has_gap(k); -1d0 elsewhere)
      !   in_gap(k)       has_gap(k) .and. y_lo < y(k) < y_hi, with y(k) the pair number
      !                   fraction of the cell (equivalent to x_lo(k) < xa(i_sed,k) < x_hi(k), but
      !                   decided in the variable of the phase diagram, so that degenerate cells
      !                   (no partner: y = 1; no pair at all: y = 0) are never in the gap because
      !                   of round-off in the mass-fraction conversion)
      ! OpenMP parallel over k; per-cell error codes are collected and the first
      ! failing cell is reported in ierr and on stdout.
      type(phase_diagram), intent(in) :: pd
      integer, intent(in) :: nz, species
      real(dp), intent(in) :: xa(:,:)  ! (species, nz) mass fractions; may be larger
      real(dp), intent(in) :: A(:)     ! (species) atomic weights
      integer, intent(in) :: i_sed, i_partner
      real(dp), intent(in) :: P(:), T(:)  ! (nz) [dyn/cm^2], [K]
      logical, intent(out) :: has_gap(:), in_gap(:)  ! (nz)
      real(dp), intent(out) :: x_lo(:), x_hi(:)      ! (nz)
      integer, intent(out) :: ierr
      integer, allocatable :: cell_ierr(:)
      integer :: k, k_fail, op_err
      real(dp) :: y_lo, y_hi, y
      include 'formats'

      ierr = 0

      if (nz < 0 .or. species < 2) then
         write(*,'(a,i0,a,i0)') 'evaluate_miscibility: need nz >= 0 and species >= 2; nz, species = ', nz, ', ', species
         ierr = miscibility_err_bad_input
         return
      end if
      if (i_sed < 1 .or. i_sed > species .or. i_partner < 1 .or. i_partner > species .or. i_sed == i_partner) then
         write(*,'(a,i0,a,i0,a,i0)') 'evaluate_miscibility: invalid species indices i_sed, i_partner = ', &
            i_sed, ', ', i_partner, ' for species = ', species
         ierr = miscibility_err_bad_input
         return
      end if
      if (size(xa,dim=1) < species .or. size(xa,dim=2) < nz .or. size(A) < species &
            .or. size(P) < nz .or. size(T) < nz .or. size(has_gap) < nz .or. size(in_gap) < nz &
            .or. size(x_lo) < nz .or. size(x_hi) < nz) then
         write(*,'(a)') 'evaluate_miscibility: array sizes are inconsistent with nz and species'
         write(*,11) 'nz, species', nz, species
         write(*,11) 'shape(xa)', shape(xa)
         write(*,11) 'size(A), size(P), size(T)', size(A), size(P), size(T)
         write(*,11) 'size(has_gap), size(in_gap), size(x_lo), size(x_hi)', &
            size(has_gap), size(in_gap), size(x_lo), size(x_hi)
         ierr = miscibility_err_bad_input
         return
      end if
      if (any(is_bad(A(1:species))) .or. any(A(1:species) <= 0d0)) then
         write(*,'(a)') 'evaluate_miscibility: atomic weights must be finite and > 0'
         write(*,1) 'A', A(1:species)
         ierr = miscibility_err_bad_input
         return
      end if
      if (.not. pd% loaded) then
         write(*,'(a)') 'evaluate_miscibility: phase diagram is not loaded'
         ierr = phase_diagram_err_not_loaded
         return
      end if

      has_gap(1:nz) = .false.
      in_gap(1:nz) = .false.
      x_lo(1:nz) = -1d0
      x_hi(1:nz) = -1d0
      if (nz == 0) return

      allocate(cell_ierr(nz))
      cell_ierr(1:nz) = 0

      !$OMP PARALLEL DO PRIVATE(k, y_lo, y_hi, y, op_err) SCHEDULE(static)
      do k = 1, nz
         if (any(is_bad(xa(1:species,k))) .or. any(xa(1:species,k) < 0d0)) then
            cell_ierr(k) = miscibility_err_composition
            cycle
         end if
         call phase_diagram_limits(pd, P(k), T(k), y_lo, y_hi, has_gap(k), op_err)
         if (op_err /= 0) then
            cell_ierr(k) = op_err
            has_gap(k) = .false.
            cycle
         end if
         if (has_gap(k)) then
            x_lo(k) = mass_fraction_for_pair_number_fraction(xa(1:species,k), A(1:species), i_sed, i_partner, y_lo)
            x_hi(k) = mass_fraction_for_pair_number_fraction(xa(1:species,k), A(1:species), i_sed, i_partner, y_hi)
            y = pair_number_fraction(xa(1:species,k), A(1:species), i_sed, i_partner)
            in_gap(k) = (y_lo < y) .and. (y < y_hi)
         end if
      end do
      !$OMP END PARALLEL DO

      k_fail = 0
      do k = 1, nz
         if (cell_ierr(k) /= 0) then
            k_fail = k
            exit
         end if
      end do

      if (k_fail > 0) then
         ierr = cell_ierr(k_fail)
         write(*,'(a,i0,a,i0,a)') 'evaluate_miscibility: failed for cell k = ', k_fail, ' with ierr = ', ierr, &
            ' (' // trim(miscibility_error_message(ierr)) // ')'
         write(*,11) 'number of failing cells', count(cell_ierr(1:nz) /= 0)
         write(*,1) 'P [dyn/cm^2]', P(k_fail)
         write(*,1) 'T [K]', T(k_fail)
         write(*,1) 'xa(i_sed)', xa(i_sed,k_fail)
         write(*,1) 'xa(i_partner)', xa(i_partner,k_fail)
         write(*,1) 'min xa in cell', minval(xa(1:species,k_fail))
      end if

   end subroutine evaluate_miscibility


   function miscibility_error_message(ierr) result(msg)
      integer, intent(in) :: ierr
      character(len=64) :: msg
      select case (ierr)
      case (miscibility_ok)
         msg = 'ok'
      case (miscibility_err_bad_input)
         msg = 'invalid sizes, species indices or atomic weights'
      case (miscibility_err_composition)
         msg = 'NaN or negative mass fraction in a cell'
      case default
         msg = phase_diagram_error_message(ierr)
      end select
   end function miscibility_error_message

end module sedimentation_miscibility
