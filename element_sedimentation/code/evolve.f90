! anchor = 'call set_to_NaN(s% total_energy_from_phase_separation)'
            call set_to_NaN(s% total_energy_from_element_sedimentation)
! anchor = 'use phase_separation, only: do_phase_separation'
         use element_sedimentation, only: do_element_sedimentation, finish_element_sedimentation
! anchor = 's% okay_to_set_mixing_info = .false.'

               ! element sedimentation (extension, not part of the official MESA release)
               if (s% do_element_sedimentation) then
                  do k = 1, s% nz
                     ! energy_start is reset later in do_struct_burn_mix (same use as in do_element_diffusion)
                     s% energy_start(k) = s% energy(k)
                  end do
                  call do_element_sedimentation(s, dt, ierr)
                  if (ierr < 0) then
                     ! configuration error (unknown method, phase diagram, species): no retry can help
                     write(*,'(a)') ' do_step_part2: do_element_sedimentation: configuration error, terminating'
                     termination_code_str(t_xtra1) = 'element_sedimentation configuration error'
                     s% termination_code = t_xtra1
                     s% result_reason = nonzero_ierr
                     do_step_part2 = terminate
                     return
                  end if
                  if (failed('do_element_sedimentation')) return
                  call set_vars_if_needed(s, dt, 'after element sedimentation', ierr)
                  if (failed('set_vars_if_needed after element sedimentation')) return
                  ! eps_element_sedimentation for the energy equation (like finish_element_diffusion)
                  call finish_element_sedimentation(s, dt)
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
