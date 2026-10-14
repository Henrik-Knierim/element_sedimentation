! anchor = 'phase_separation_no_diffusion, &'

    ! element sedimentation (extension, not part of the official MESA release)
    do_element_sedimentation, &
    element_sedimentation_method, &
    element_sedimentation_chemical_species, &
    element_sedimentation_partner_species, &
    element_sedimentation_phase_diagram, &
    element_sedimentation_phase_diagram_temperature_shift, &
    element_sedimentation_extrapolate_phase_diagram, &
    element_sedimentation_velocity_method, &
    element_sedimentation_velocity, &
    element_sedimentation_D_mix_factor, &
    element_sedimentation_cap_width, &
    element_sedimentation_solve_to_steady_state, &
    element_sedimentation_steady_state_tolerance, &
    element_sedimentation_max_steps, &
    element_sedimentation_advection_diffusion_time_step_factor, &
    element_sedimentation_maximum_dXi_per_cell, &
    element_sedimentation_retry_change_factor, &
    element_sedimentation_cell_threshold_for_precipitation, &
    element_sedimentation_suppress_mixing_in_immiscible_regions, &
    element_sedimentation_replace_mesa_mixing, &
    element_sedimentation_overturn_unstable_layers, &
    element_sedimentation_act_on_mixing_contrast, &
    element_sedimentation_mixing_contrast_tolerance, &
    element_sedimentation_transfer_rounds, &
    element_sedimentation_heating, &
    element_sedimentation_include_enthalpy_flux, &
    element_sedimentation_report, &
! anchor = 's% phase_separation_no_diffusion = phase_separation_no_diffusion'

 ! element sedimentation (extension, not part of the official MESA release)
 s% do_element_sedimentation = do_element_sedimentation
 s% element_sedimentation_method = element_sedimentation_method
 s% element_sedimentation_chemical_species = element_sedimentation_chemical_species
 s% element_sedimentation_partner_species = element_sedimentation_partner_species
 s% element_sedimentation_phase_diagram = element_sedimentation_phase_diagram
 s% element_sedimentation_phase_diagram_temperature_shift = element_sedimentation_phase_diagram_temperature_shift
 s% element_sedimentation_extrapolate_phase_diagram = element_sedimentation_extrapolate_phase_diagram
 s% element_sedimentation_velocity_method = element_sedimentation_velocity_method
 s% element_sedimentation_velocity = element_sedimentation_velocity
 s% element_sedimentation_D_mix_factor = element_sedimentation_D_mix_factor
 s% element_sedimentation_cap_width = element_sedimentation_cap_width
 s% element_sedimentation_solve_to_steady_state = element_sedimentation_solve_to_steady_state
 s% element_sedimentation_steady_state_tolerance = element_sedimentation_steady_state_tolerance
 s% element_sedimentation_max_steps = element_sedimentation_max_steps
 s% element_sedimentation_advection_diffusion_time_step_factor = element_sedimentation_advection_diffusion_time_step_factor
 s% element_sedimentation_maximum_dXi_per_cell = element_sedimentation_maximum_dXi_per_cell
 s% element_sedimentation_retry_change_factor = element_sedimentation_retry_change_factor
 s% element_sedimentation_cell_threshold_for_precipitation = element_sedimentation_cell_threshold_for_precipitation
 s% element_sedimentation_suppress_mixing_in_immiscible_regions = element_sedimentation_suppress_mixing_in_immiscible_regions
 s% element_sedimentation_replace_mesa_mixing = element_sedimentation_replace_mesa_mixing
 s% element_sedimentation_overturn_unstable_layers = element_sedimentation_overturn_unstable_layers
 s% element_sedimentation_act_on_mixing_contrast = element_sedimentation_act_on_mixing_contrast
 s% element_sedimentation_mixing_contrast_tolerance = element_sedimentation_mixing_contrast_tolerance
 s% element_sedimentation_transfer_rounds = element_sedimentation_transfer_rounds
 s% element_sedimentation_heating = element_sedimentation_heating
 s% element_sedimentation_include_enthalpy_flux = element_sedimentation_include_enthalpy_flux
 s% element_sedimentation_report = element_sedimentation_report
! anchor = 'phase_separation_no_diffusion = s% phase_separation_no_diffusion'

 ! element sedimentation (extension, not part of the official MESA release)
 do_element_sedimentation = s% do_element_sedimentation
 element_sedimentation_method = s% element_sedimentation_method
 element_sedimentation_chemical_species = s% element_sedimentation_chemical_species
 element_sedimentation_partner_species = s% element_sedimentation_partner_species
 element_sedimentation_phase_diagram = s% element_sedimentation_phase_diagram
 element_sedimentation_phase_diagram_temperature_shift = s% element_sedimentation_phase_diagram_temperature_shift
 element_sedimentation_extrapolate_phase_diagram = s% element_sedimentation_extrapolate_phase_diagram
 element_sedimentation_velocity_method = s% element_sedimentation_velocity_method
 element_sedimentation_velocity = s% element_sedimentation_velocity
 element_sedimentation_D_mix_factor = s% element_sedimentation_D_mix_factor
 element_sedimentation_cap_width = s% element_sedimentation_cap_width
 element_sedimentation_solve_to_steady_state = s% element_sedimentation_solve_to_steady_state
 element_sedimentation_steady_state_tolerance = s% element_sedimentation_steady_state_tolerance
 element_sedimentation_max_steps = s% element_sedimentation_max_steps
 element_sedimentation_advection_diffusion_time_step_factor = s% element_sedimentation_advection_diffusion_time_step_factor
 element_sedimentation_maximum_dXi_per_cell = s% element_sedimentation_maximum_dXi_per_cell
 element_sedimentation_retry_change_factor = s% element_sedimentation_retry_change_factor
 element_sedimentation_cell_threshold_for_precipitation = s% element_sedimentation_cell_threshold_for_precipitation
 element_sedimentation_suppress_mixing_in_immiscible_regions = s% element_sedimentation_suppress_mixing_in_immiscible_regions
 element_sedimentation_replace_mesa_mixing = s% element_sedimentation_replace_mesa_mixing
 element_sedimentation_overturn_unstable_layers = s% element_sedimentation_overturn_unstable_layers
 element_sedimentation_act_on_mixing_contrast = s% element_sedimentation_act_on_mixing_contrast
 element_sedimentation_mixing_contrast_tolerance = s% element_sedimentation_mixing_contrast_tolerance
 element_sedimentation_transfer_rounds = s% element_sedimentation_transfer_rounds
 element_sedimentation_heating = s% element_sedimentation_heating
 element_sedimentation_include_enthalpy_flux = s% element_sedimentation_include_enthalpy_flux
 element_sedimentation_report = s% element_sedimentation_report
