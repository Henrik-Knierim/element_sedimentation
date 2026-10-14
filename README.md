# element_sedimentation

Element sedimentation ("rain") for [MESA](https://docs.mesastar.org), developed for
helium rain in giant planets but written for any species that becomes immiscible in the
surrounding mixture. The species is transported by a conservative advection-diffusion
equation in MESA's Lagrangian grid: convective mixing acts as diffusion, rain as advection
limited by the miscibility gap of a phase diagram. All other species are moved by the
displaced counter-flow, so mass fractions keep summing to one and every species conserves
its mass.

The extension is not part of the official MESA release. `setup.py` copies the sources
into a MESA tree, patches a handful of MESA files and recompiles. Tested with MESA r26.4.1.

## Quickstart

```
export MESA_DIR=/path/to/mesa          # a built MESA
cd element_sedimentation
python3 setup.py check                  # verify all patch anchors in this MESA and the consistency of the controls
python3 setup.py install -y             # patch + copy + full MESA rebuild (./clean; ./install)
python3 setup.py status                 # what is installed
```

Then, in your work directory, run `./clean` and `./mk` and add to `&controls`:

```fortran
   do_element_sedimentation = .true.
   element_sedimentation_phase_diagram = 'src/Schoettler+2018.csv'
   element_sedimentation_velocity = 1d2            ! cm/s, positive = sinking
   element_sedimentation_report = .true.
```

`planet_helium_rain/` is a complete example (Jupiter-like planet; see its `rn*` scripts).
`python3 setup.py uninstall` restores MESA. `soft_install` rebuilds only `star_data` and
`star` and is the quick way to recompile after editing the sources.

Unit tests of the solver, phase diagram and composition modules (no installation needed):

```
cd element_sedimentation/test && ./mk && ./rn
```

## Options

All controls are documented in `element_sedimentation/code/controls.defaults` (installed
into `$MESA_DIR/star/defaults/controls.defaults`):

| control | default |
|---|---|
| `do_element_sedimentation` | `.false.` |
| `element_sedimentation_method` | `'advection_diffusion_rain'` (or `'instant_rain'`, `'bottom_up_rain'`) |
| `element_sedimentation_chemical_species` | `'he4'` |
| `element_sedimentation_partner_species` | `'h1'` |
| `element_sedimentation_phase_diagram` | `''` (CSV: `y_lo, y_hi, T [K], P [Mbar]`) |
| `element_sedimentation_phase_diagram_temperature_shift` | `0d0` |
| `element_sedimentation_extrapolate_phase_diagram` | `.false.` |
| `element_sedimentation_velocity_method` | `'constant'` |
| `element_sedimentation_velocity` | `1d2` (cm/s, positive = sinking) |
| `element_sedimentation_D_mix_factor` | `1d0` |
| `element_sedimentation_cap_width` | `1d-3` |
| `element_sedimentation_solve_to_steady_state` | `.true.` |
| `element_sedimentation_steady_state_tolerance` | `1d-8` |
| `element_sedimentation_max_steps` | `1000` |
| `element_sedimentation_advection_diffusion_time_step_factor` | `0.5d0` |
| `element_sedimentation_maximum_dXi_per_cell` | `-1d0` (off) |
| `element_sedimentation_retry_change_factor` | `0.5d0` |
| `element_sedimentation_cell_threshold_for_precipitation` | `1` |
| `element_sedimentation_suppress_mixing_in_immiscible_regions` | `.false.` |
| `element_sedimentation_replace_mesa_mixing` | `.true.` |
| `element_sedimentation_overturn_unstable_layers` | `.true.` |
| `element_sedimentation_act_on_mixing_contrast` | `.true.` |
| `element_sedimentation_mixing_contrast_tolerance` | `1d-6` |
| `element_sedimentation_transfer_rounds` | `5` |
| `element_sedimentation_heating` | `.true.` |
| `element_sedimentation_include_enthalpy_flux` | `.true.` |
| `element_sedimentation_report` | `.false.` |

## Structure

```
element_sedimentation/
   setup.py                 installer: install | soft_install | uninstall | reinstall | single_* | status | check | clean
   code/                    Fortran sources and MESA patch snippets
      advection_diffusion.f90          generic implicit advection-diffusion solver (no MESA dependency)
      sedimentation_phase_diagram.f90  phase diagram reader and interpolation
      sedimentation_miscibility.f90    composition arithmetic and per-cell miscibility limits
      sedimentation_energy.f90         energy carried by the moved species (enthalpy flux and its divergence)
      sedimentation_overturn.f90       instantaneous mixing of layers left convectively unstable
      element_sedimentation.f90        MESA driver (do_element_sedimentation, finish_element_sedimentation, shutdown_element_sedimentation)
      makefile_base, controls.defaults, star_controls.inc, ctrls_io.f90,
      star_data_step_work.inc, alloc.f90, hydro_energy.f90, evolve.f90, init.f90   patch snippets
   test/                    unit tests (./mk, ./rn, ./clean), data/Schoettler+2018.csv,
                            integration/compare_methods.py (runs an example with several methods)
   backup/                  pristine MESA files (created by setup.py, not in git)
docs/index.md               physics, numerics, controls, installation, tests
docs/CHANGELOG.md           what changed and why (2026-10 refactor)
docs/TODO.md                open points
planet_helium_rain/         example work directory (1 M_J)
saturn_helium_rain/         example work directory (1 M_Saturn; the helium gap reaches deeper)
```

See `docs/index.md` for the equations, the numerical method and the details of the
installer.
