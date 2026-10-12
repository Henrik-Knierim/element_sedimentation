! anchor = 'real(dp), pointer :: eps_phase_separation(:)'

      ! element sedimentation (extension, not part of the official MESA release)
      real(dp), pointer :: eps_element_sedimentation(:) ! (nz)
! anchor = 'total_energy_from_phase_separation, &'
            total_energy_from_element_sedimentation, &  ! element sedimentation (extension)
