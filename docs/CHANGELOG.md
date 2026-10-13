# Changelog

Details of what changed and why; `docs/index.md` describes the current state, `docs/TODO.md`
the open points.

## 2026-10: refactor into modules, advection-diffusion transport, mixing and energy

Branch `develop`, five commits (split, unit tests, installer, docs, mixing/energy).

### Architecture
- The 3485-line `element_sedimentation.f90` is split into library modules without MESA
  dependency, `advection_diffusion` (solver), `sedimentation_phase_diagram` (CSV reader,
  interpolation), `sedimentation_miscibility` (composition arithmetic, per-cell limits),
  `sedimentation_energy` (energy carried by the moved species) and `sedimentation_overturn`
  (mixing of unstable layers), plus the MESA driver `element_sedimentation` and
  `element_sedimentation_legacy` with the old `instant_rain` and `bottom_up_rain`
  algorithms, unchanged.
- The two components of the phase diagram are the controls
  `element_sedimentation_chemical_species` and `element_sedimentation_partner_species`;
  atomic weights come from `chem_isos% W` for every species in the net. The old code had
  `h1`/`he4` hard-coded and used integer mass numbers.

### Transport (`'advection_diffusion_rain'`, the new default method)
- Conservative advection-diffusion equation in the Lagrangian mass coordinate: convective
  mixing as diffusion (MESA's `sig` from `D_mix`), rain as upwind advection with the carrier
  mass flow `4 pi r^2 rho v`. Backward Euler, semismooth Newton for the per-cell floor (only
  the species above the He-poor limit rains) and the receiver cap (inflow stops at the
  He-rich limit), pseudo-transient continuation to the steady state of each MESA step
  (`element_sedimentation_solve_to_steady_state`) or integration over `dt`.
- All other species are moved by the displaced counter-flow, so every cell keeps
  `sum(xa) = 1` and every species conserves its mass; both are checked after every step.
- Damping (`element_sedimentation_maximum_dXi_per_cell`, retry damping) blends from the
  mixing-only state. Every failure returns a named error with a diagnostic block and makes
  MESA retry the step; nothing calls `stop`.

### Mixing
- `element_sedimentation_replace_mesa_mixing` (default true): in steps where the module acts
  it sets `D_mix` and `sig` to zero after its transport, so MESA's convection no longer
  undoes the rain; `mix_factor` can stay at 1 (the example used `mix_factor = 0`).
- `element_sedimentation_act_on_mixing_contrast` (default true): the module also takes over
  steps in which MESA would mix a composition contrast of the species, because MESA's own
  mixing carries no enthalpy flux (it cooled the core by ~160 K when its convection remixed
  the helium-rich interior into the envelope).
- `element_sedimentation_overturn_unstable_layers` (default true): layers the transport
  leaves Ledoux-unstable (the deposit on cells whose `D_mix` was zero) are mixed at once
  by the pool-adjacent-violators algorithm with MESA's Ledoux composition term. Without it
  MESA's MLT put an inverted temperature gradient at the composition inversion.

### Energy
- `eps_element_sedimentation` = change of the internal energy at fixed (rho, T) by the
  composition change (needed in both energy forms, since MESA resets `xa_start` and
  `energy_start` after the hook) plus, with `element_sedimentation_include_enthalpy_flux`
  (default true), the divergence of the enthalpy flux `sum_i hbar_i J_i` of the moved
  species, with the partial specific enthalpies at fixed (T, P) from (rho, T) derivatives
  of the EOS. The term sums to zero over the star. No separate gravitational source: the
  release of gravitational energy enters through the structure equations in both forms.
- Without the flux every cell kept its internal energy while its composition changed,
  which cooled the drained envelope by `hbar_H - hbar_He` per gram exchanged (9/8 kT/m_p
  for an ideal gas, 6-7 kT/m_p in the dense H/He fluid of the planet), halved L, lowered
  T_eff by 15 K and crashed the example at 4.82 Gyr. With it the example reaches 5 Gyr in
  the `dedt` and the `eps_grav` form without retries, T_eff within a few K of the run
  without rain. The latent heat of demixing is still missing (TODO).

### Installer, tests, documentation
- `setup.py` rewritten: anchor-based patching of a pristine MESA with backups, `check`
  (anchors, snippets, control consistency), `status`, `install`, `soft_install`
  (`star_data` and `star` only), `uninstall`, `reinstall`; sets `OMP_STACKSIZE=512m` for
  the MESA build (bus error in MESA's `colors` test on macOS otherwise).
- Unit tests of the five library modules (`element_sedimentation/test`, 505 checks),
  linked against the MESA SDK libraries without installing anything.
- `docs/index.md` (physics, numerics, controls, installation, tests), `README.md`,
  `docs/TODO.md`; the example `planet_helium_rain` updated (`mix_factor = 1`).
- Sign convention: `element_sedimentation_velocity` is now positive for sinking.
