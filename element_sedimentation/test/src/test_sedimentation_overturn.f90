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

module test_sedimentation_overturn

   ! unit tests for sedimentation_overturn: pool-adjacent-violators mixing of
   ! composition inversions with a caller-supplied face test

   use const_def, only: dp
   use sedimentation_overturn
   use test_support, only: assert_true, assert_close, assert_equal_int

   implicit none

   private
   public :: run_sedimentation_overturn_tests

   real(dp), parameter :: rtol = 1d-12, atol = 1d-300
   integer, parameter :: i_heavy = 2   ! the species whose inversion counts as unstable
   integer :: barrier_face = 0         ! a face that is never unstable (0: none)

contains

   subroutine run_sedimentation_overturn_tests(nfail)
      integer, intent(inout) :: nfail
      write(*,'(a)') 'test_sedimentation_overturn'
      call test_stable_profile_unchanged(nfail)
      call test_single_inversion(nfail)
      call test_cascade_to_the_centre(nfail)
      call test_upward_propagation(nfail)
      call test_barrier_face(nfail)
   end subroutine run_sedimentation_overturn_tests


   logical function heavy_above(k, xa_above, xa_below)
      ! unstable if the block above has more of the heavy species, unless k is the barrier
      integer, intent(in) :: k
      real(dp), intent(in) :: xa_above(:), xa_below(:)
      heavy_above = xa_above(i_heavy) > xa_below(i_heavy) + 1d-14
      if (k == barrier_face) heavy_above = .false.
   end function heavy_above


   subroutine set_two_species(xa, x_heavy)
      ! species 1 = light, 2 = heavy, 3 = trace that follows the light one
      real(dp), intent(out) :: xa(:,:)
      real(dp), intent(in) :: x_heavy(:)
      integer :: k
      do k = 1, size(x_heavy)
         xa(2,k) = x_heavy(k)
         xa(3,k) = 0.02d0*(1d0 - x_heavy(k))
         xa(1,k) = 1d0 - xa(2,k) - xa(3,k)
      end do
   end subroutine set_two_species


   subroutine check_invariants(name, dm, xa_old, xa_new, nfail)
      ! every species keeps its mass, every cell sums to one
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: dm(:), xa_old(:,:), xa_new(:,:)
      integer, intent(inout) :: nfail
      integer :: j, k, nz, species
      nz = size(dm)
      species = size(xa_old, 1)
      do j = 1, species
         call assert_close(name // ': mass of species conserved', dot_product(dm, xa_new(j,1:nz)), &
            dot_product(dm, xa_old(j,1:nz)), rtol, atol, nfail)
      end do
      do k = 1, nz
         call assert_close(name // ': sum of mass fractions', sum(xa_new(1:species,k)), 1d0, rtol, atol, nfail)
      end do
   end subroutine check_invariants


   subroutine test_stable_profile_unchanged(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 5, species = 3
      real(dp) :: dm(nz), xa(species,nz), xa0(species,nz)
      integer :: n_faces, n_blocks, k_outer, k_inner
      dm = [1d0, 2d0, 1d0, 3d0, 1d0]
      call set_two_species(xa, [0.25d0, 0.25d0, 0.27d0, 0.30d0, 0.30d0])
      xa0 = xa
      barrier_face = 0
      call overturn_unstable_layers(nz, species, dm, xa, heavy_above, n_faces, n_blocks, k_outer, k_inner)
      call assert_equal_int('stable: no face merged', n_faces, 0, nfail)
      call assert_equal_int('stable: no block', n_blocks, 0, nfail)
      call assert_equal_int('stable: k_outer 0', k_outer, 0, nfail)
      call assert_equal_int('stable: k_inner 0', k_inner, 0, nfail)
      call assert_true('stable: composition unchanged', all(xa == xa0), nfail)
   end subroutine test_stable_profile_unchanged


   subroutine test_single_inversion(nfail)
      ! a heavy cell on top of a lighter one, the rest stable: only the pair is mixed
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 5, species = 3
      real(dp) :: dm(nz), xa(species,nz), xa0(species,nz), x_mix
      integer :: n_faces, n_blocks, k_outer, k_inner
      dm = [1d0, 2d0, 1d0, 3d0, 1d0]
      call set_two_species(xa, [0.20d0, 0.30d0, 0.28d0, 0.40d0, 0.45d0])
      xa0 = xa
      barrier_face = 0
      call overturn_unstable_layers(nz, species, dm, xa, heavy_above, n_faces, n_blocks, k_outer, k_inner)
      x_mix = (2d0*0.30d0 + 1d0*0.28d0)/3d0
      call assert_equal_int('single: one face merged', n_faces, 1, nfail)
      call assert_equal_int('single: one block', n_blocks, 1, nfail)
      call assert_equal_int('single: k_outer', k_outer, 2, nfail)
      call assert_equal_int('single: k_inner', k_inner, 3, nfail)
      call assert_close('single: cell 2 mixed', xa(2,2), x_mix, rtol, atol, nfail)
      call assert_close('single: cell 3 mixed', xa(2,3), x_mix, rtol, atol, nfail)
      call assert_true('single: other cells untouched', &
         all(xa(:,1) == xa0(:,1)) .and. all(xa(:,4) == xa0(:,4)) .and. all(xa(:,5) == xa0(:,5)), nfail)
      call assert_true('single: trace species uniform in the block', xa(3,2) == xa(3,3), nfail)
      call check_invariants('single', dm, xa0, xa, nfail)
   end subroutine test_single_inversion


   subroutine test_cascade_to_the_centre(nfail)
      ! a deposit on top of a uniform interior: the block grows to the centre
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 8, species = 3
      real(dp) :: dm(nz), xa(species,nz), xa0(species,nz), x_mix, m_int
      integer :: n_faces, n_blocks, k_outer, k_inner, k
      dm = [1d0, 1d0, 1d0, 2d0, 2d0, 2d0, 2d0, 2d0]
      call set_two_species(xa, [0.20d0, 0.20d0, 0.21d0, 0.31d0, 0.28d0, 0.28d0, 0.28d0, 0.28d0])
      xa0 = xa
      barrier_face = 0
      call overturn_unstable_layers(nz, species, dm, xa, heavy_above, n_faces, n_blocks, k_outer, k_inner)
      m_int = sum(dm(4:nz))
      x_mix = dot_product(dm(4:nz), xa0(2,4:nz))/m_int
      call assert_equal_int('cascade: four faces merged', n_faces, 4, nfail)
      call assert_equal_int('cascade: one block', n_blocks, 1, nfail)
      call assert_equal_int('cascade: block starts at the deposit', k_outer, 4, nfail)
      call assert_equal_int('cascade: block ends at the centre', k_inner, nz, nfail)
      do k = 4, nz
         call assert_close('cascade: uniform interior', xa(2,k), x_mix, rtol, atol, nfail)
      end do
      call assert_true('cascade: envelope untouched', all(xa(:,1:3) == xa0(:,1:3)), nfail)
      call assert_true('cascade: interior heavier than the envelope (stable)', xa(2,4) > xa(2,3), nfail)
      call check_invariants('cascade', dm, xa0, xa, nfail)
   end subroutine test_cascade_to_the_centre


   subroutine test_upward_propagation(nfail)
      ! merging lowers the block mean below the block above, which must then be merged too
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 4, species = 3
      real(dp) :: dm(nz), xa(species,nz), xa0(species,nz), x_mix
      integer :: n_faces, n_blocks, k_outer, k_inner, k
      dm = 1d0
      call set_two_species(xa, [0.28d0, 0.30d0, 0.20d0, 0.20d0])
      xa0 = xa
      barrier_face = 0
      call overturn_unstable_layers(nz, species, dm, xa, heavy_above, n_faces, n_blocks, k_outer, k_inner)
      x_mix = sum(xa0(2,1:nz))/nz
      call assert_equal_int('upward: three faces merged', n_faces, 3, nfail)
      call assert_equal_int('upward: one block', n_blocks, 1, nfail)
      call assert_equal_int('upward: whole column', k_outer, 1, nfail)
      call assert_equal_int('upward: whole column (inner)', k_inner, nz, nfail)
      do k = 1, nz
         call assert_close('upward: uniform', xa(2,k), x_mix, rtol, atol, nfail)
      end do
      call check_invariants('upward', dm, xa0, xa, nfail)
   end subroutine test_upward_propagation


   subroutine test_barrier_face(nfail)
      ! the face test can refuse a face (e.g. a Ledoux-stable one): blocks stop there
      integer, intent(inout) :: nfail
      integer, parameter :: nz = 6, species = 3
      real(dp) :: dm(nz), xa(species,nz), xa0(species,nz), x_mix
      integer :: n_faces, n_blocks, k_outer, k_inner
      dm = [1d0, 1d0, 1d0, 1d0, 1d0, 1d0]
      call set_two_species(xa, [0.20d0, 0.20d0, 0.30d0, 0.28d0, 0.26d0, 0.26d0])
      xa0 = xa
      barrier_face = 5   ! between cells 4 and 5
      call overturn_unstable_layers(nz, species, dm, xa, heavy_above, n_faces, n_blocks, k_outer, k_inner)
      x_mix = 0.5d0*(0.30d0 + 0.28d0)
      call assert_equal_int('barrier: one face merged', n_faces, 1, nfail)
      call assert_equal_int('barrier: one block', n_blocks, 1, nfail)
      call assert_equal_int('barrier: block 3-4 (outer)', k_outer, 3, nfail)
      call assert_equal_int('barrier: block 3-4 (inner)', k_inner, 4, nfail)
      call assert_close('barrier: mixed above the barrier', xa(2,3), x_mix, rtol, atol, nfail)
      call assert_close('barrier: mixed above the barrier (4)', xa(2,4), x_mix, rtol, atol, nfail)
      call assert_true('barrier: below the barrier untouched', all(xa(:,5:6) == xa0(:,5:6)), nfail)
      call assert_true('barrier: inversion across the barrier remains', xa(2,4) > xa(2,5), nfail)
      call check_invariants('barrier', dm, xa0, xa, nfail)
      barrier_face = 0
   end subroutine test_barrier_face

end module test_sedimentation_overturn
