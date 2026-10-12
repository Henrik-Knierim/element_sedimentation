! anchor = 'call set_to_NaN(s% total_energy_from_phase_separation)'
            call set_to_NaN(s% total_energy_from_element_sedimentation)
! anchor = 'use phase_separation, only: do_phase_separation'
         use element_sedimentation, only: do_element_sedimentation
! anchor = 's% okay_to_set_mixing_info = .false.'

               ! element sedimentation (extension, not part of the official MESA release)
               if (s% do_element_sedimentation) then
                  do k = 1, s% nz
                     ! energy_start is reset later in do_struct_burn_mix (same use as in do_element_diffusion)
                     s% energy_start(k) = s% energy(k)
                  end do
                  call do_element_sedimentation(s, dt, ierr)
                  if (failed('do_element_sedimentation')) return
                  call set_vars_if_needed(s, dt, 'after element sedimentation', ierr)
                  if (failed('set_vars_if_needed after element sedimentation')) return
                  do k = 1, s% nz  ! for use by the energy equation
                     s% eps_element_sedimentation(k) = (s% energy_start(k) - s% energy(k))/dt
                  end do
               end if
! anchor = 's% total_energy_from_phase_separation = 0d0'; prepend
            s% total_energy_from_element_sedimentation = 0d0
            if (s% do_element_sedimentation) then
               s% total_energy_from_element_sedimentation = &
                  dt*dot_product(s% dm(1:nz), s% eps_element_sedimentation(1:nz))
            end if

! anchor = '+ s% total_energy_from_phase_separation &'
               + s% total_energy_from_element_sedimentation &
! anchor = '- s% total_energy_from_phase_separation &'; all
               - s% total_energy_from_element_sedimentation &
