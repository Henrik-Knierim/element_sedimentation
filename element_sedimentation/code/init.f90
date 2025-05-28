! anchor = 'use paquette_coeffs, only: free_collision_integrals'
         use element_sedimentation, only: shutdown_phase_diagram
! anchor = 'call free_collision_integrals()'
         call shutdown_phase_diagram()
! anchor = 'use star_utils, only: init_random'
         use element_sedimentation, only: init_element_sedimentation
! anchor = 'call init_random(s)'

         ! HHK: load phase diagram
         call init_element_sedimentation(s, ierr)