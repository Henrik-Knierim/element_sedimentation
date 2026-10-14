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

! Unit tests for the element_sedimentation modules (run from test/ via ./rn).

program tester

   use test_support, only: print_test_summary
   use test_advection_diffusion, only: run_advection_diffusion_tests
   use test_instant_transfer, only: run_instant_transfer_tests
   use test_phase_diagram, only: run_phase_diagram_tests
   use test_miscibility, only: run_miscibility_tests
   use test_sedimentation_energy, only: run_sedimentation_energy_tests
   use test_sedimentation_overturn, only: run_sedimentation_overturn_tests

   implicit none

   integer :: nfail

   nfail = 0
   call run_advection_diffusion_tests(nfail)
   call run_instant_transfer_tests(nfail)
   call run_phase_diagram_tests(nfail)
   call run_miscibility_tests(nfail)
   call run_sedimentation_energy_tests(nfail)
   call run_sedimentation_overturn_tests(nfail)
   call print_test_summary(nfail)
   if (nfail > 0) error stop 'element_sedimentation tests FAILED'

end program tester
