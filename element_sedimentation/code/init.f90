! anchor = 'use paquette_coeffs, only: free_collision_integrals'
        use element_sedimentation, only: shutdown_element_sedimentation
! anchor = 'call free_collision_integrals()'
        call shutdown_element_sedimentation()
