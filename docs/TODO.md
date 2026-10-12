# TODO

## Physics
- [ ] Velocity models beyond `'constant'`: droplet settling velocity from droplet size,
      viscosity and density contrast (Stokes/turbulent regimes), possibly depending on the
      local supersaturation; `element_sedimentation_velocity_method` is the hook.
- [ ] Include the latent heat of demixing (currently only the internal-energy change of
      the composition change enters the energy equation).
- [ ] Mass transfer on top of abundance transfer (droplets carry mass, not only composition).
- [ ] Decide whether cells that are saturated at the He-poor limit should stay unmixed
      after the step (`element_sedimentation_suppress_mixing_in_immiscible_regions` uses the
      pre-step gap cells and does not update `mixing_type`).

## Numerics
- [ ] The steady-state residual is normalized with the longest face timescale; purely
      diffusive modes relax on ~nz^2 face timescales, so the slowest mode may be ~1e-4 from
      its limit at the default tolerance (tighten `element_sedimentation_steady_state_tolerance`
      or raise the solver's `max_dt_factor`). Consider a residual that uses the diffusive
      relaxation time of the whole convective region.
- [ ] Time-accurate mode (`solve_to_steady_state = .false.`) has only been run on the
      `planet_helium_rain` example, where the MESA steps (~1e7 yr) are far longer than the
      transport timescales, so it reproduces the steady state there. A case with
      comparable timescales (slow rain, short steps) is still untested in MESA.

## Example run (`planet_helium_rain`)
- [ ] With the default settings the 1 M_J model rains from ~4.75 Gyr on (model ~497). Once
      the convective envelope above the rain region has been drained to Y ~ 0.15-0.2,
      MESA's structure solver fails for good (default settings: model ~512;
      `element_sedimentation_maximum_dXi_per_cell = 0.01`: model ~524, 4.87 Gyr). The
      retries fail even when the composition change of the retried step is 1e-15, so the
      failure comes from the depleted state, not from the transport: the CMS equation of
      state returns `chiT < 0` and `grada < 0` for the He-poor cells near logT ~ 3.5,
      logRho ~ -0.2 (P ~ 0.6 Mbar) and `store_eos_for_cell` (micro.f90) rejects them.
      Options: blend to another EOS in that region (`inlist_eos` CMS limits), or revisit the
      phase diagram shift (-350 K) that puts the He-poor limit this low.
      (An earlier version of the damping stopped the transport integration early; that run
      reached 5 Gyr only because the envelope was never properly drained.)
- [ ] The He-poor envelope value is set by the He-poor limit at the top of the rain region
      (the convective envelope drains to it); the jump between steps is large when the rain
      region deepens quickly. A velocity law or time-accurate runs with short steps would
      show the transient.
- [ ] In steady-state mode the whole composition change of a step, and its heat, is applied
      in every step regardless of `dt`; when MESA cuts the time step on retries,
      `eps_element_sedimentation` = dE/dt grows accordingly (the retry damping reduces dE in
      step with dt). Check whether this hinders convergence near the failure above.

## Code and project
- [ ] Modules other than `star` that use `star_info` (binary, astero) are only rebuilt by
      the full `install`; `soft_install` rebuilds `star_data` and `star`.
- [ ] Legacy `bottom_up_rain`: `dXi` can exceed `element_sedimentation_maximum_dXi_per_cell`
      during core saturation (old known issue, kept with the legacy code).
