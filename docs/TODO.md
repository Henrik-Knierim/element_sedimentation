# TODO

## Physics
- [ ] Velocity models beyond `'constant'`: droplet settling velocity from droplet size,
      viscosity and density contrast (Stokes/turbulent regimes), possibly depending on the
      local supersaturation; `element_sedimentation_velocity_method` is the hook.
- [ ] Latent heat of demixing. The energy of the moved species is the partial enthalpy from
      the equation of state (`sedimentation_energy.f90`), so a heat of mixing enters only if
      the EOS has one; the CMS tables are combined with the additive-volume rule, so the
      enthalpy of forming and dissolving the droplets is missing. Its size from the phase
      diagram itself (van 't Hoff slope of the He-poor limit in `Schoettler+2018.csv`, or
      `T dsbar = (kT/m_He) ln(x_rich/x_poor)` with the ideal mixing entropy): 0.5-0.8 kT
      per proton mass per gram of He where the 1 M_J model rains (1-2 Mbar, 4000-6000 K),
      i.e. 20-50 % of the enthalpy-exchange term, not a small residual (Markham & Guillot
      2024 adopt 2e12 erg/g). It is released where droplets form and absorbed where they
      dissolve, so it is a second energy flux of opposite sign, nearly zero net; implement
      it as the face enthalpy of He in the He-rich phase, `hbar_He(poor) - L`, inside the
      same flux (a divergence), never as a one-sided heating.
- [ ] Sedimentation velocity relative to what. The carrier flow `w = 4 pi r^2 rho v` moves
      the mass fraction with `v` relative to the mass shells (the barycentre), and in that
      frame the mass transport is complete (the He flux is balanced by an equal upward mass
      flux of the rest; "droplets carry mass" needs nothing more). A droplet's terminal
      velocity is relative to the surrounding fluid, and the droplets are the He-rich phase,
      so the droplet mass flux is `J_He/(Y_rich - Y_mix)` (`J_He/(1 - Y)` for pure He).
      Irrelevant for the steady state per step, which the floor and cap set, but it matters
      for time-accurate runs and a velocity law.
- [ ] The helium gradient at the base of the drained envelope is a Ledoux-stable layer
      (one or two cells with `gradr < gradL`, `mixing_type = 0`) that throttles the
      luminosity from the interior to 0.1-0.2 of its value in the run without rain. The
      interior below then stores heat and the envelope radiates from its own content.
      Decide whether such a layer should get a thermohaline/double-diffusive coefficient
      (MESA's `thermohaline_coeff`, or an ODDC prescription as in Mankovich et al. 2016)
      instead of being either fully convective or not mixed at all; this, not the energy
      bookkeeping, now decides the sign and size of the T_eff difference.
- [ ] The rain comes in bursts (a limit cycle): the envelope drains until its helium
      reaches the He-poor limit at the top of the rain region, the layer below warms (heat
      stored under the stable layer plus the compression work of the deposit), the limit at
      the gap rises above the envelope value, no cell is saturated any more and the module
      stops raining; once MESA's Ledoux criterion finds the interface convective
      (`gradr > grada + B`, the warmed interior under the cooled envelope), the interior is
      mixed back into the envelope over a few steps (the module does that mixing with its
      energy terms, `element_sedimentation_act_on_mixing_contrast`), T_eff rises by a few K,
      and the envelope rains again ~100 Myr later. Decide whether this is physical (it is
      what the phase diagram, the heating and the Ledoux criterion give) or whether the
      interface should stay stratified (thermohaline instead of full convective remixing
      of a stabilizing helium gradient). The remixing hands the whole deposit back to the
      envelope, so by 5 Gyr the planet has rained out only ~0.7 % of its helium net.
- [ ] Long run (20 Gyr requested, `dedt`, 2026-10-14, without velocity time centering, see
      `docs/index.md` Physics): the planet rains from 4.65 Gyr on, Y_env 0.2765 -> 0.145 at
      9.7 Gyr, where the CMS EOS failure of the first item ends the run
      (`store_eos_for_cell` for the He-poor envelope cells). The released gravitational
      energy (W - W_base = -7.8e40 erg) is stored below the Ledoux-stable helium gradient:
      the interior is up to ~1300 K warmer at the centre than the no-rain planet (an earlier
      run with velocity time centering showed the centre 230 K cooler, an artefact of the
      half-charged work term), the envelope cooler. L/L_base falls to 0.68 at 7 Gyr
      (T_eff -6.4 K at 6 Gyr) and climbs back to 0.88 at 9.7 Gyr: the time-integrated
      luminosity deficit, 4.5e40 erg, equals the stored (U + W) - base, i.e. the excess
      luminosity of the Stevenson & Salpeter / Fortney & Hubbard picture has not started yet
      5 Gyr after onset in this model, because the barrier is a radiative layer with no
      thermohaline/ODDC transport. In the Saturn example, where the gap reaches the centre
      and a helium ocean forms, the planet is 36 % more luminous and 4.4 K warmer than
      without rain at 5 Gyr. Needed to continue the Jupiter run: the EOS blend for He-poor
      cells (first item) and the barrier transport (item above). The `eps_grav` form is not
      usable for this experiment: its cumulative energy error (2.4e-3 of the total energy,
      ~2e40 erg) is as large as the stored energy.
- [ ] The radiative atmosphere above the convective envelope keeps its initial helium
      while the envelope below drains (a heavier layer on top of a lighter one, with a
      stable thermal stratification: thermohaline regime, so the overturn leaves it to
      MESA). Check whether MESA's thermohaline mixing should be on in the example, and
      whether the He-poor limit should be applied there.
- [ ] The `eps_grav` form of the energy equation carries a cumulative energy error of 2e-3
      in this planet even without rain (`dedt`: 7e-8); the rain adds nothing to it
      (2.366e-3 vs 2.365e-3), but comparisons must be made within one energy form.

## Numerics
- [ ] `'instant_rain'` deviates from the infinite-velocity limit of the solver in three
      documented ways (docs/index.md, section 8): cells at their cap are fixed and separate
      mixed regions (the solver's diffusion crosses them), a cell above its cap never rains,
      and the limit is taken from the homogenized state. A region-wise water filling with
      soft caps (caps throttle advective inflow only) would remove the first; whether a
      cell above the He-rich limit should rain at all is a physics question for the solver
      too (its advective term ignores the donor's cap).
- [ ] The composition change of a steady-state step does not depend on `dt`, so the
      heating `eps_element_sedimentation` grows as 1/dt when MESA shortens the step (the
      Saturn ocean formation triggers this). The retry damping now scales with the
      time-step ratio, and a per-cell change limit (`element_sedimentation_maximum_dXi_per_cell`,
      0.02 in the Saturn example) spreads the ocean formation over many steps; a transport
      integrated over the physical `dt` (the solver's time-accurate mode) would bound the
      rate physically, but only for steps shorter than the rain time scale (~1 yr at
      1 m/s), far below MESA's steps.
- [ ] The steady-state residual is normalized with the longest face timescale; purely
      diffusive modes relax on ~nz^2 face timescales, so the slowest mode may be ~1e-4 from
      its limit at the default tolerance (tighten `element_sedimentation_steady_state_tolerance`
      or raise the solver's `max_dt_factor`). Consider a residual that uses the diffusive
      relaxation time of the whole convective region.
- [ ] Energy terms: the partial enthalpies are cell-centred finite differences of the EOS
      (`dxa_eos = 1e-3` on the coarse composition grid of the CMS tables) averaged to the
      faces with the cell masses; donor-cell values would be an alternative (both first
      order, both conserve energy). Neither the step nor the averaging has been tested for
      resolution sensitivity. `eps_element_sedimentation` is explicit (no Jacobian partials,
      like MESA's `eps_diffusion`), and in steady-state mode it scales as 1/dt when MESA
      cuts the step; the planet runs needed no retries, larger per-step changes are untested.
- [ ] Time-accurate mode (`solve_to_steady_state = .false.`) has only been run on the
      `planet_helium_rain` example, where the MESA steps (~1e7 yr) are far longer than the
      transport timescales, so it reproduces the steady state there. A case with
      comparable timescales (slow rain, short steps) is still untested in MESA.

## Example run (`planet_helium_rain`)
- [ ] With the default settings the 1 M_J model rains from ~4.75 Gyr on (model ~497) and
      now reaches 5 Gyr in both energy forms without retries. Before the enthalpy flux was
      added, the spurious cooling of the drained envelope (see `docs/index.md`, Physics)
      lowered the He-poor limit further (positive feedback), the envelope drained to
      Y ~ 0.2, and MESA's structure solver failed for good at ~4.82 Gyr: the CMS equation
      of state returns `chiT < 0` and `grada < 0` for He-poor cells near logT ~ 3.5,
      logRho ~ -0.2 (P ~ 0.6 Mbar) and `store_eos_for_cell` (micro.f90) rejects them.
      That state is no longer reached, but a model that drains further (other phase
      diagram shift, other mass) will meet it again: blend to another EOS in that region
      (`inlist_eos` CMS limits) or revisit the phase diagram shift (-350 K).
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
