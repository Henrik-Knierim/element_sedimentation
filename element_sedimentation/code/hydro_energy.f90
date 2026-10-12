! anchor = 'others_ad%val = others_ad%val + s% eps_phase_separation(k)'
            if (s% do_element_sedimentation .and. s% element_sedimentation_heating) &
               others_ad%val = others_ad%val + s% eps_element_sedimentation(k)
