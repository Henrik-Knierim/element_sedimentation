# Changelog

Details of what changed and why; `docs/index.md` describes the current state, `docs/TODO.md`
the open points.

## 2026-10: refactor into modules, advection-diffusion transport, mixing and energy

Branch `develop`; commits: split, unit tests, installer, docs, mixing/energy, the other
methods with the Saturn example, removal of the old implementations.

### Architecture
- The 3485-line `element_sedimentation.f90` is split into library modules without MESA
  dependency, `advection_diffusion` (solver), `sedimentation_phase_diagram` (CSV reader,
  interpolation), `sedimentation_miscibility` (composition arithmetic, per-cell limits),
  `sedimentation_energy` (energy carried by the moved species) and `sedimentation_overturn`
  (mixing of unstable layers), plus the MESA driver `element_sedimentation`. The old
  `instant_rain` and `bottom_up_rain` algorithms were kept in `element_sedimentation_legacy`
  for the comparison below and removed afterwards (2026-10-14).
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
  MESA retry the step; nothing calls `stop`. Configuration errors (an unknown method, a
  phase diagram that cannot be loaded, a species missing from the net) end the run at once
  with the message instead of a retry cascade (2026-10-14).

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

### The other methods (`'instant_rain'`, `'bottom_up_rain'`)
- Both re-implemented on the same floors, caps, counter-flow, damping, overturn, energy
  terms, checks and mixing takeover as the solver; only the transport core differs.
  `'instant_rain'` is the steady state for infinitely fast rain, computed directly by
  `adv_diff_instant_transfer` in the solver module (drain to the floors, fill the cells
  below up to their limits, back-fill to the caps when the way down is blocked, which
  makes a helium ocean grow upward from the centre; mixed regions settle to a water
  filling). `'bottom_up_rain'` repeats it with the limits re-evaluated after every
  transfer (`element_sedimentation_transfer_rounds`), which is what the old bottom-up
  sweep was for.
- The old implementations were kept as `'instant_rain_legacy'` and `'bottom_up_rain_legacy'`
  for the comparison (`docs/index.md`, section 8) and then removed, together with
  `element_sedimentation_legacy.f90`. Their known problems: the deposit region
  and the helium ocean are found by heuristics, the He-rich limit is not enforced when
  depositing, the bottom-up sweep costs n_z evaluations per step, and the bottom-up
  variant leaves MESA's convection on, which refills the raining cells every step; with
  the energy terms this futile cycle pumps energy (T_eff +33 K in the Jupiter example).
- New: a Saturn example (`saturn_helium_rain`), the integration test
  `test/integration/compare_methods.py` that runs an example with several methods and
  compares them, and unit tests of the instant transfer.
- Found with the Saturn example, where the helium gap reaches the centre and a helium
  ocean forms: MESA's `dedt` work term with velocity time centering charges a cell that
  contracts within one step only half its P dV work, which cooled the forming ocean
  from 6700 K to 2900 K and crashed every method; the examples now run without
  velocity time centering (`docs/index.md`, Physics). The retry damping now scales the
  allowed change with the time-step ratio of the retry (otherwise the heating rate of
  the step stays constant through a retry cascade), and it limits the whole change of
  the step relative to its start rather than to the mixing-only state.
- An independent review of the instant transfer (three reviewers, adversarial
  verification) led to: no limits in mixing-only steps, no stray diffusion in the
  counter-flow of the instant transfer, floors acting only on cells that advect, exact
  filling to the cap, and the instant methods ignoring the velocity; the remaining
  documented deviations from the solver are listed in `docs/index.md`, section 8.

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
