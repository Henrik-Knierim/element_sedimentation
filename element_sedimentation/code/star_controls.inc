! anchor = 'character (len=32) :: phase_separation_option'

      ! element sedimentation (extension, not part of the official MESA release)

         logical :: do_element_sedimentation
         character (len=strlen) :: element_sedimentation_method
         character (len=strlen) :: element_sedimentation_chemical_species
         character (len=strlen) :: element_sedimentation_partner_species
         character (len=strlen) :: element_sedimentation_phase_diagram
         real(dp) :: element_sedimentation_phase_diagram_temperature_shift
         logical :: element_sedimentation_extrapolate_phase_diagram
         character (len=strlen) :: element_sedimentation_velocity_method
         real(dp) :: element_sedimentation_velocity
         real(dp) :: element_sedimentation_D_mix_factor
         real(dp) :: element_sedimentation_cap_width
         logical :: element_sedimentation_solve_to_steady_state
         real(dp) :: element_sedimentation_steady_state_tolerance
         integer :: element_sedimentation_max_steps
         real(dp) :: element_sedimentation_advection_diffusion_time_step_factor
         real(dp) :: element_sedimentation_maximum_dXi_per_cell
         real(dp) :: element_sedimentation_retry_change_factor
         integer :: element_sedimentation_cell_threshold_for_precipitation
         logical :: element_sedimentation_suppress_mixing_in_immiscible_regions
         logical :: element_sedimentation_disable_conv_premix
         logical :: element_sedimentation_replace_mesa_mixing
         logical :: element_sedimentation_overturn_unstable_layers
         logical :: element_sedimentation_act_on_mixing_contrast
         real(dp) :: element_sedimentation_mixing_contrast_tolerance
         integer :: element_sedimentation_transfer_rounds
         logical :: element_sedimentation_heating
         logical :: element_sedimentation_include_enthalpy_flux
         logical :: element_sedimentation_report
