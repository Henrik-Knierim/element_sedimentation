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

module sedimentation_overturn

   ! Instantaneous mixing of composition inversions left behind by a transport
   ! step (no MESA dependency).
   !
   ! A transport step that deposits a heavy species into cells that were not
   ! mixed at the start of the step leaves a layer that is denser than the one
   ! below it. Such a layer overturns on a dynamical timescale, long before the
   ! next step, so it is mixed here at once instead of being handed to the
   ! stellar code, whose mixing coefficients of the step were fixed before the
   ! transport and whose convection treatment would otherwise react to the
   ! inversion with an inverted temperature gradient.
   !
   ! The cells are gathered from the surface inward into blocks of uniform
   ! composition by the pool-adjacent-violators algorithm: a cell starts its
   ! own block; while the face between the newest block and the block above it
   ! is unstable (decided by a caller-supplied test, e.g. the Ledoux criterion
   ! with the equation of state), the two blocks are merged, every species
   ! mass-weighted, and the face above the merged block is tested again. Every
   ! species conserves its mass exactly and the mass fractions of a mixed cell
   ! still sum to one. Cells, faces and masses follow the MESA convention: cell
   ! k has outer face k and inner face k+1, k = 1 is the surface.

   use const_def, only: dp

   implicit none

   private
   public :: overturn_unstable_layers, face_is_unstable_interface

   abstract interface
      logical function face_is_unstable_interface(k, xa_above, xa_below)
         ! .true. if the face k between the block ending in cell k-1 (composition
         ! xa_above) and the block starting in cell k (composition xa_below) is
         ! unstable and the two blocks should be mixed
         import dp
         integer, intent(in) :: k
         real(dp), intent(in) :: xa_above(:)  ! (species)
         real(dp), intent(in) :: xa_below(:)  ! (species)
      end function face_is_unstable_interface
   end interface

contains

   subroutine overturn_unstable_layers(nz, species, dm, xa, face_is_unstable, &
         n_merged_faces, n_mixed_blocks, k_outer, k_inner)
      ! Mix every maximal run of cells connected by unstable faces to a uniform
      ! composition. n_merged_faces is the number of faces mixed away, n_mixed_blocks
      ! the number of resulting blocks with more than one cell, k_outer/k_inner the
      ! outermost/innermost cell of any mixed block (0 if nothing was mixed).
      integer, intent(in) :: nz, species
      real(dp), intent(in) :: dm(:)         ! (nz) cell masses
      real(dp), intent(inout) :: xa(:,:)    ! (species, nz) mass fractions
      procedure(face_is_unstable_interface) :: face_is_unstable
      integer, intent(out) :: n_merged_faces, n_mixed_blocks, k_outer, k_inner

      integer :: k, nb, b
      integer, allocatable :: bstart(:), bend(:)
      real(dp), allocatable :: bmass(:), bxa(:,:)

      n_merged_faces = 0
      n_mixed_blocks = 0
      k_outer = 0
      k_inner = 0
      if (nz < 2) return

      allocate(bstart(nz), bend(nz), bmass(nz), bxa(species,nz))
      nb = 0
      do k = 1, nz
         nb = nb + 1
         bstart(nb) = k
         bend(nb) = k
         bmass(nb) = dm(k)
         bxa(1:species,nb) = xa(1:species,k)
         do while (nb >= 2)
            if (.not. face_is_unstable(bstart(nb), bxa(1:species,nb-1), bxa(1:species,nb))) exit
            bxa(1:species,nb-1) = (bmass(nb-1)*bxa(1:species,nb-1) + bmass(nb)*bxa(1:species,nb)) &
               /(bmass(nb-1) + bmass(nb))
            bmass(nb-1) = bmass(nb-1) + bmass(nb)
            bend(nb-1) = bend(nb)
            nb = nb - 1
            n_merged_faces = n_merged_faces + 1
         end do
      end do

      do b = 1, nb
         if (bend(b) == bstart(b)) cycle
         n_mixed_blocks = n_mixed_blocks + 1
         if (k_outer == 0) k_outer = bstart(b)
         k_inner = bend(b)
         do k = bstart(b), bend(b)
            xa(1:species,k) = bxa(1:species,b)
         end do
      end do

   end subroutine overturn_unstable_layers

end module sedimentation_overturn
