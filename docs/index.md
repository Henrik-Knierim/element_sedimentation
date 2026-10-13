# element_sedimentation

Element sedimentation ("rain") for MESA. The extension models the transport of a
chemical species that becomes immiscible in the surrounding mixture, e.g. helium
rain in giant planets, as a conservative advection-diffusion problem in MESA's
Lagrangian grid. It is not part of the official MESA release: `setup.py` copies
the sources into a MESA tree, patches a handful of MESA files by anchor insertion
and rebuilds `star_data` and `star`.

Tested with MESA r26.4.1 (mesasdk gfortran 15.2). `docs/CHANGELOG.md` summarizes what the
2026-10 refactor changed and why; `docs/TODO.md` lists the open points.

## 1. Architecture

```
element_sedimentation/code/
   advection_diffusion.f90          generic conservative 1-D advection-diffusion solver (no MESA dependency)
   sedimentation_phase_diagram.f90  immiscibility phase diagram: CSV reader, bilinear interpolation (no MESA dependency)
   sedimentation_miscibility.f90    composition arithmetic, per-cell miscibility limits (no MESA dependency)
   sedimentation_energy.f90         energy carried by the moved species: face flow rates, enthalpy flux, its divergence (no MESA dependency)
   sedimentation_overturn.f90       instantaneous mixing of layers left convectively unstable, pool-adjacent-violators (no MESA dependency)
   element_sedimentation.f90        MESA-facing driver: builds the problem from star_info, applies the result, Ledoux test, energy terms
   element_sedimentation_legacy.f90 the old 'instant_rain' and 'bottom_up_rain' algorithms
   <patch snippets>                 makefile_base, controls.defaults, star_controls.inc, ctrls_io.f90,
                                    star_data_step_work.inc, alloc.f90, hydro_energy.f90, evolve.f90, init.f90
element_sedimentation/setup.py      installer (see section 6)
element_sedimentation/test/         self-checking unit tests of the five library modules (see section 7)
planet_helium_rain/                 example MESA work directory (Jupiter-like planet with helium rain)
```

Call flow per MESA step (`evolve.f90`, right after `okay_to_set_mixing_info = .false.`,
i.e. after the mixing information of the step is fixed):

```
do_element_sedimentation(s, dt, ierr)                         element_sedimentation.f90
   initialize                                                 load phase diagram once, resolve species, atomic weights
   eps_element_sedimentation = 0
   select element_sedimentation_method
      'advection_diffusion_rain' -> do_advection_diffusion_rain
            evaluate_miscibility(pd, ..., s% xa, A, s% Peos, s% T)   sedimentation_miscibility.f90
            build adv_diff_problem (dm, sig, w, x_floor, x_cap)
            adv_diff_integrate(p, x, total_time, opts, stats, ierr, after_step)   advection_diffusion.f90
                after every accepted step: after_step -> adv_diff_passive_step (other species)
            damping (change limit), overturn_unstable_layers (Ledoux-unstable faces, EOS)   sedimentation_overturn.f90
            post checks (bounds, sum = 1, species masses), MESA's mixing switched off for the step, report
      'instant_rain' / 'bottom_up_rain' -> element_sedimentation_legacy.f90
   enthalpy_flux_heating (flow rates from the composition change, partial enthalpies from the EOS)
   need_to_setvars = .true.
set_vars_if_needed                                            evolve.f90 snippet
finish_element_sedimentation(s, dt)                           eps_element_sedimentation = (energy_start - energy)/dt
                                                                 + divergence of the enthalpy flux
```

The library modules depend only on `const_def`, `utils_lib` (and `math_lib`), so
they are unit-tested without installing anything into MESA.

## 2. Physics

Species i (`element_sedimentation_chemical_species`, default `he4`) with mass fraction
$X_k$ in cell $k$ ($k = 1$ at the surface, $k = n_z$ at the centre). Faces are numbered
like cells: face $k$ is the outer face of cell $k$; faces 1 and $n_z + 1$ carry zero flux.

Conservative transport in the Lagrangian mass coordinate:

$$
\mathrm{d}m_k \frac{\mathrm{d}X_k}{\mathrm{d}t} = F_k - F_{k+1}, \qquad
F_k = \underbrace{\sigma_k (X_{k-1} - X_k)}_{\text{convective mixing}}
    + \underbrace{w_k\, \max(X_d - X_{\mathrm{floor},d},\,0)\, \theta_r}_{\text{rain}}
$$

* $\sigma_k = (4\pi r_k^2 \rho_{\mathrm{face},k})^2 D_k / \overline{\mathrm{d}m}_k$ [g/s] is
  MESA's `sig` built from `D_mix(k)` (`mix_info.f90: get_convection_sigmas`), multiplied
  by `element_sedimentation_D_mix_factor`; $\overline{\mathrm{d}m}_k = (\mathrm{d}m_{k-1} + \mathrm{d}m_k)/2$,
  $\rho_{\mathrm{face},k} = (\mathrm{d}q_{k-1}\rho_k + \mathrm{d}q_k \rho_{k-1})/(\mathrm{d}q_{k-1} + \mathrm{d}q_k)$.
* $w_k = 4\pi r_k^2 \rho_{\mathrm{face},k} v_k$ [g/s] is the mass flow rate of the carrier
  (the droplets), $v_k > 0$ meaning sinking. The donor cell is $d = k-1$ for $w_k > 0$ and
  $d = k$ otherwise; the receiver $r$ is the other cell. Only `'constant'`
  ($v_k = $ `element_sedimentation_velocity`) is implemented.
* Floor: only the part of the species above the He-poor miscibility limit
  $X_{\mathrm{lo}}(P_d, T_d)$ of the donor rains. Cells whose $(P, T)$ lie outside the
  immiscibility region never rain (floor $= \infty$).
* Throttle: $\theta_r = \min(1, \max(0, (X_{\mathrm{cap},r} - X_r)/\Delta))$ shuts the inflow
  off as the receiver approaches its cap, the He-rich limit $X_{\mathrm{hi}}(P_r, T_r)$
  (cap $= 1$ where there is no gap); $\Delta$ = `element_sedimentation_cap_width`. This is
  how a helium-rich layer grows upward from the centre.

Summing over cells telescopes to $\sum_k \mathrm{d}m_k\, \mathrm{d}X_k/\mathrm{d}t = 0$: the
total mass of the species is conserved exactly by construction.

The miscibility limits come from the phase diagram, tabulated in the number fraction
$y = n_i/(n_i + n_p)$ of the two-component (i, partner) subsystem, and are converted to
mass fractions of the cell keeping the number ratios of all other species to the
partner fixed (`mass_fraction_for_pair_number_fraction`).

**Other species.** Falling droplets displace an equal mass of the surrounding mixture
upward. Each other species $j$ is diffused with the same $\sigma$ and advected with the
face rate $w'_k = -F_{\mathrm{adv},k}$ (the actual rain flux of the active species) and
donor scale $1/(1 - X_{i,d'})$ at the new time level:

$$
G_{j,k} = \sigma_k (X_{j,k-1} - X_{j,k}) + w'_k\, \frac{X_{j,d'}}{1 - X_{i,d'}}
$$

Solved implicitly with the same time step, $S = \sum_{j \ne i} X_j$ satisfies exactly
the equation solved by $1 - X_i$, so $\sum_j X_j = 1$ holds to round-off after every
step, every species conserves its total mass, and positivity is preserved (M-matrix
with upwinding).

**Energy.** The composition change is applied at fixed $(\rho, T)$ before MESA's
structure solve. Two things have to be right: the total energy and where the heat goes.

* *Total energy.* As for `element_diffusion` and `phase_separation`, the internal energy
  before and after the composition change gives the term
  `(energy_start - energy)/dt`: every cell keeps its internal energy while its
  composition changes. MESA resets `energy_start` and `xa_start` after this hook
  (`struct_burn_mix`), so neither the `dedt` form nor the composition term of the
  `eps_grav` form sees the sedimentation, and the term is needed, and counted once,
  in both forms of the energy equation. The gravitational energy released by the rain
  needs no term of its own: in the Lagrangian frame the sinking helium is replaced by
  an equal mass of the other species, gravity does no net work on that exchange
  ($\sum_i J_i = 0$, so $g \cdot \sum_i J_i = 0$), and the release
  $-dW/dt = -\int P\, d(1/\rho)/dt\, dm$ is delivered by the structure equations in
  both forms, as $P\,dV$ work in `dedt` and inside the full $-T\,ds/dt$ in `eps_grav`.
  Locally that work is booked where the volume changes: the helium-receiving layers
  are compressed, the helium-losing ones expand. A separate gravitational source
  (a $g v J$ term or a constant) on top of this would count the release twice.
* *Where the heat goes.* The moved species carry energy. For a multicomponent fluid the
  energy flux relative to the barycentre is $\sum_i \bar h_i J_i$ with the partial
  specific enthalpies $\bar h_i$ at fixed $(T, P)$ (de Groot & Mazur ch. III; Bird,
  Stewart & Lightfoot eq. 19.3-3). `sedimentation_energy.f90` builds the mass flow rates
  of every species through every face from the composition change of the step, the
  driver evaluates $\partial h/\partial X_i|_{T,P}$ by finite differences of the EOS at
  fixed $(\rho, T)$,
  $(\partial h/\partial X)_{T,P} = (\partial e/\partial X)_{\rho,T} +
  (\chi_T/\chi_\rho)\,(\partial P/\partial X)_{\rho,T}/\rho$,
  and the divergence of the flux is added to `eps_element_sedimentation`
  (`element_sedimentation_include_enthalpy_flux`). Written out, the term contains the
  buoyancy work $g \sum_i \bar v_i J_i$ of the sinking species (the gravitational
  energy dissipated by drag where it falls), the flow work, and the heat carried by the
  exchange of hot material from below against cold material from above. It sums to
  zero over the star: it moves the released energy from where the hydrostatic work
  books it (the compressed deposit, the expanding source) to where the droplets fall.
  Without it every cell keeps its enthalpy at its pressure while its composition
  changes: the helium-depleted cells cool by $(\bar h_{\rm H} - \bar h_{\rm He})\,
  \delta X / c_p$ and the helium-enriched ones heat by the same rule. For an ideal
  H$_2$/He gas $\bar h_{\rm H} - \bar h_{\rm He} = 9/8\,kT/m_p$ (internal energy
  $7/8$ plus flow work $1/4$); in the dense fluid at the rain region of a Jupiter-mass
  planet (1-2 Mbar, 4000 K) the CMS equation of state gives 6-7 $kT/m_p$, split about
  equally between the internal-energy and the $P\bar v$ part at fixed $(T, P)$, or, in
  the $(\rho, T)$ form the code evaluates, ~88 % from $(\partial e/\partial X)_{\rho,T}$
  (pressure-dissociated hydrogen against neutral helium). That is why the old
  bookkeeping cooled the rain-top cells by 25-30 % of $T$ for $\delta Y = -0.066$,
  halved $L$ and lowered $T_{\rm eff}$ by 15 K; `element_sedimentation_report` prints
  the ratio to $kT/m_p$ at the top and bottom of the rain region. The partial
  enthalpies of the additive-volume CMS mixture contain no enthalpy of mixing, so the
  latent heat of the demixing itself is not included (`docs/TODO.md`).

`eps_element_sedimentation` enters the energy equation when `element_sedimentation_heating`
is true and the total energy accounting through `total_energy_from_element_sedimentation`.
As for MESA's phase separation, the accounting term is always included, so
`element_sedimentation_heating = .false.` makes the reported energy-conservation error
reflect the omitted heating and triggers `hard_limit_for_rel_error_in_energy_conservation`
retries unless that limit is relaxed.

## 3. Numerics (`advection_diffusion.f90`)

* Backward Euler (fully implicit) in time, first-order upwind for advection, central
  differences for diffusion. The implicit matrix is tridiagonal with positive diagonal and
  non-positive off-diagonals (an M-matrix after scaling by $\mathrm{d}m_k$), so the Thomas
  algorithm without pivoting is stable, the scheme is unconditionally stable and
  positivity preserving, and no minimum-diffusion hack or post-hoc mass rescaling is needed.
* Floor and cap make a step piecewise linear; it is solved by a semismooth Newton
  iteration (generalized derivatives at the kinks) with a backtracking line search. A step
  counts as converged when the corrections are negligible and every cell-equation defect is
  below 1e-8 (in units of the mass fraction) or at the round-off floor of its Jacobian row;
  otherwise the step is rejected and retried with a smaller pseudo time step. The
  tridiagonal solves use pivots computed from the exact column sums ($\mathrm{d}m_k$), so
  they stay accurate for any stiffness ($\Delta t\,\sigma/\mathrm{d}m \gg 10^{16}$).
* The steady state is only accepted once the pseudo time step exceeds the slowest face
  timescale, and the change of the other species (reported by the driver after every step)
  enters the steady-state test, so the solver does not stop while they are still mixing.
* Steady state (default, `element_sedimentation_solve_to_steady_state = .true.`) by
  pseudo-transient continuation: the pseudo time step starts at
  `element_sedimentation_advection_diffusion_time_step_factor` times the shortest face
  timescale $\tau_k = \min(\mathrm{d}m_{k-1}, \mathrm{d}m_k)/(\sigma_k + |w_k|)$ and grows
  geometrically (switched evolution relaxation), so the number of steps is logarithmic
  in $\tau_{\max}/\tau_{\min}$. Converged when
  $\max_k |\Delta X_k| / \Delta t \cdot \tau_{\max} \le$ `element_sedimentation_steady_state_tolerance`.
  Note that diffusion-dominated modes relax on $\sim n^2 \tau_{\mathrm{face}}$; tighten the
  tolerance if the slowest mode matters.
* Time-accurate mode (`solve_to_steady_state = .false.`) integrates exactly to the MESA
  time step, with the same step cap ($\Delta t \le 10^4\,\tau_{\max}$). It ends early,
  reporting `steady_state`, once the profile is steady by the same criterion, because the
  coefficients are frozen during the step; for MESA steps much longer than the rain and
  mixing timescales it therefore reproduces the steady-state result.
* `element_sedimentation_maximum_dXi_per_cell > 0` is a ceiling on the change caused by
  rain: the transport is also solved without rain (mixing only, same `D_mix`), and if rain
  changes $X_i$ in some cell by more than the limit relative to that mixing-only state, the
  step is scaled down for all species, $X \leftarrow X_{\rm mix} + f\,(X_{\rm rain} -
  X_{\rm mix})$ with $f$ = limit / largest change. In unmixed (e.g. radiative) cells
  $X_{\rm mix}$ is the composition at the start of the step, so the ceiling applies there
  directly; in mixed regions it applies to the deviation from the homogenized state. Both
  states conserve every species, sum to one and are non-negative, so the blend does too,
  and mixed regions stay mixed (a layer that was stable and has become convective is
  homogenized instead of keeping part of a possibly inverted composition jump). The
  remaining rain happens in the following steps. (Stopping the integration early, or
  raising the He-poor floors, was tried and rejected: the first leaves gradients inside
  convective zones, the second does not bound the change of convective envelopes.)
* The module acts in every step in which at least one cell is saturated (mass fraction at
  or above the He-poor limit; cells left exactly at their limit by the previous step count).
  With `element_sedimentation_replace_mesa_mixing = .true.` (default) it then replaces
  MESA's composition mixing for the rest of the step: its transport has already mixed every
  species with `D_mix` together with the rain, so `D_mix` and `sig` are set to zero. MESA
  does not recompute `D_mix` within a step after the hook (`okay_to_set_mixing_info` is
  false) but rebuilds `sig` from `D_mix` in `do_struct_burn_mix`, which is why `D_mix`
  itself is zeroed. In steps without saturated cells MESA mixes as usual, so `mix_factor`
  can stay at 1 (`mix_factor = 0` was a workaround needed before this existed).
* The transport spreads the deposited rain with the mixing of the start of the step.
  Where that mixing is zero (a layer below the rain region that has become radiative
  because the helium gradient throttles the luminosity) the deposit sits as a heavier
  layer on top of lighter material, and MESA, whose mixing coefficients are fixed for the
  step, would react with an inverted temperature gradient (`gradT = gradL < grada`, even
  negative) that pumps luminosity into the deposit and out of the core. With
  `element_sedimentation_overturn_unstable_layers = .true.` (default) such layers are
  mixed at once after the damping (`sedimentation_overturn.f90`): the cells are gathered
  from the surface inward into blocks of uniform composition (pool-adjacent-violators)
  as long as the face between a block and the block above is convectively unstable by
  the Ledoux criterion, `gradr > grada + B`, with `B` evaluated like MESA's `brunt.f90`
  from the equation of state at the new compositions of the two sides and `gradr`,
  `grada` from the start of the step. A deposit on top of an interior of uniform
  composition is thus mixed through the whole interior, as a dynamical overturn would
  do; faces that are only thermohaline-unstable are left to MESA. Every species is
  mass-weighted, so masses and sums are preserved exactly; the energy of the overturn
  is part of the composition change the energy terms see.
* MESA's own composition mixing carries no enthalpy flux either (every cell keeps its
  internal energy while its composition changes). After rain the helium-poor envelope
  lies on the helium-rich interior; once MESA's convection reaches through that interface
  (in the planet run this happens when the interior has warmed), the exchange of
  hydrogen-rich and helium-rich material is of the same kind as the rain, and with
  $\bar h_{\rm H} - \bar h_{\rm He}$ of several $kT/m_p$ MESA's bookkeeping cooled the core
  by ~160 K for a 0.3 % change of $Y$. With
  `element_sedimentation_act_on_mixing_contrast = .true.` (default) the module therefore
  also acts in steps without saturated cells whenever a face with `D_mix > 0` separates
  cells whose mass fractions of the species differ by more than
  `element_sedimentation_mixing_contrast_tolerance`: mixing only (no rain), with the
  passive species, the overturn and the energy terms, and MESA's mixing switched off for
  the step. The report line marks such steps with `mix` instead of `rain`.
* A steady-state solve does not depend on MESA's time step, so when MESA retries a step
  the same composition change would be reproduced on every retry. With
  `element_sedimentation_retry_change_factor > 0` (default 0.5) each retry of a model
  limits the largest change per cell to that factor times the largest change applied in
  the previous attempt (in the same way as `element_sedimentation_maximum_dXi_per_cell`),
  compounding with further retries; the deferred rain happens in later steps.
  A notice is printed whenever this limit is applied.
* Every accepted step is checked for mass conservation and bounds ($0 \le X \le 1$).
  All failures return a named error code, a message (`stats% message`) and the offending
  cell (`stats% k_fail`); nothing in the solver calls `stop`.

On a solver failure the driver prints a diagnostic block (model number, error, message,
steps, Newton iterations, last dt, timescales, residual, and for the failing cell $X$,
floor, cap, $\sigma$, $w$, $P$, $T$), restores `s% xa`, and returns `ierr /= 0` so that
MESA retries the step. After every implicit step the passive species of each cell are
rescaled to exactly $1 - X_i$ (in exact arithmetic they already are; in floating point they
miss it by the round-off floor of the accepted step, ~ epsilon times the ratio of the step
to the fastest cell timescale: ~1e-11 in the example, up to ~1e-6 in a forced-rain test;
a defect above 1e-3 is an error), and their change enters the steady-state test. After a successful solve the driver checks that every mass
fraction is non-negative, that every cell sums to one (renormalized if the deviation
exceeds `1e-10`, error above `1e-8`) and that the total mass of every species changed by
less than `1e-8` relative.

## 4. Controls (`&controls`)

| control | default | meaning |
|---|---|---|
| `do_element_sedimentation` | `.false.` | switch |
| `element_sedimentation_method` | `'advection_diffusion_rain'` | or legacy `'instant_rain'`, `'bottom_up_rain'` |
| `element_sedimentation_chemical_species` | `'he4'` | sedimenting isotope |
| `element_sedimentation_partner_species` | `'h1'` | other component of the phase diagram |
| `element_sedimentation_phase_diagram` | `''` | CSV file `y_lo, y_hi, T [K], P [Mbar]` |
| `element_sedimentation_phase_diagram_temperature_shift` | `0d0` | K, added to the table temperatures |
| `element_sedimentation_extrapolate_phase_diagram` | `.false.` | linear extrapolation in T above the table |
| `element_sedimentation_velocity_method` | `'constant'` | only option so far |
| `element_sedimentation_velocity` | `1d2` | cm/s, **positive = sinking** |
| `element_sedimentation_D_mix_factor` | `1d0` | multiplies `D_mix` in the diffusive term (independent of MESA's `mix_factor`) |
| `element_sedimentation_cap_width` | `1d-3` | width of the inflow throttle below the He-rich limit |
| `element_sedimentation_solve_to_steady_state` | `.true.` | steady state vs. integration over `dt` |
| `element_sedimentation_steady_state_tolerance` | `1d-8` | dimensionless residual |
| `element_sedimentation_max_steps` | `1000` | implicit steps per MESA step |
| `element_sedimentation_advection_diffusion_time_step_factor` | `0.5d0` | first step / shortest face timescale |
| `element_sedimentation_maximum_dXi_per_cell` | `-1d0` | ceiling on the change of $X_i$ caused by rain per cell and step, relative to the mixing-only state (off if negative) |
| `element_sedimentation_retry_change_factor` | `0.5d0` | on a MESA retry of the same model, limit the change per cell to this factor times the previous attempt's largest change (off if <= 0) |
| `element_sedimentation_cell_threshold_for_precipitation` | `1` | minimum number of saturated cells (at or above the He-poor limit) for the module to act |
| `element_sedimentation_suppress_mixing_in_immiscible_regions` | `.false.` | zero `D_mix`/`sig` around cells in the gap; disables `do_conv_premix` (advection-diffusion method only; `instant_rain` always does this, `bottom_up_rain` never) |
| `element_sedimentation_replace_mesa_mixing` | `.true.` | in steps where the module acts, set `D_mix` and `sig` to zero after its transport (MESA's mixing would undo the rain) |
| `element_sedimentation_overturn_unstable_layers` | `.true.` | mix layers the transport leaves Ledoux-unstable (heavy deposit on lighter material) at once, every species mass-weighted |
| `element_sedimentation_act_on_mixing_contrast` | `.true.` | also act (mixing only, with the energy terms) in steps where a mixed face separates cells of different species mass fraction; MESA's mixing carries no enthalpy flux |
| `element_sedimentation_mixing_contrast_tolerance` | `1d-6` | contrast across a mixed face above which the previous control applies |
| `element_sedimentation_heating` | `.true.` | include `eps_element_sedimentation` in the energy equation |
| `element_sedimentation_include_enthalpy_flux` | `.true.` | the moved species carry their partial specific enthalpy (EOS at fixed T, P); off: every cell keeps its internal energy, as for MESA's element diffusion |
| `element_sedimentation_report` | `.false.` | summary line per MESA step and one line per solver step |

The full documentation of each control is in `code/controls.defaults` and, after
installation, in `$MESA_DIR/star/defaults/controls.defaults`.

Note: the sign convention of `element_sedimentation_velocity` changed with the
refactoring; positive values now sink (the old code used negative values for sinking).

## 5. Phase diagram file

CSV with one header line followed by rows `y_lo, y_hi, T [K], P [Mbar]`
(`planet_helium_rain/src/Schoettler+2018.csv` is an example). Rows are grouped in blocks
of constant P with T strictly increasing inside a block and P strictly increasing between
blocks; at least two blocks and two temperatures per block are needed. Interpolation is
linear in T inside the two bracketing P blocks and then linear in P. A cell has a gap
when P is inside the table range, T is inside the T range of both bracketing blocks (or
extrapolation is on) and the interpolated `y_lo < y_hi`.

## 6. Installation (`setup.py`)

```
export MESA_DIR=/path/to/mesa
cd element_sedimentation
python3 setup.py check          # dry run: are all anchors found in this MESA version?
python3 setup.py install -y     # patch, copy, full MESA rebuild (./clean; ./install); sets OMP_STACKSIZE=512m
                                # unless defined (MESA's colors test crashes on macOS with the default thread stack,
                                # and a failed ./install leaves MESA without libraries)
python3 setup.py status
```

| command | effect |
|---|---|
| `install` | patch MESA, copy the five sources, full MESA rebuild |
| `soft_install` | same, but rebuild only `star_data` and `star` (`./clean; ./mk; ./export`) |
| `uninstall` / `soft_uninstall` | restore the backups, remove the sources, rebuild `star_data` and `star` |
| `reinstall` / `soft_reinstall` | uninstall + install / soft_install |
| `single_install` / `single_uninstall` / `single_reinstall -f NAME` | one file and its module |
| `status` | patched/pristine and up-to-date/STALE state of every file, MESA version, backup version |
| `check` | parse all snippets, verify every anchor without modifying MESA, and check that every control is declared, initialized and read consistently |
| `clean` | delete the backup directory (refused while MESA is patched) |

Options: `-y/--yes` (no prompts), `--no-build` (patch/restore only), `--force` (ignore
the MESA version recorded with the backups), `-v`, `--mesa-dir DIR`.

Patch snippets consist of blocks `! anchor = '<text of a MESA line>'[; replace|; prepend|; all]`
followed by the lines to insert after (default), instead of, before, or after every
occurrence of the anchor line. Backups of the pristine MESA files are kept in
`element_sedimentation/backup/` (created on install, never committed), each with a record
of the MESA version it was taken from (`<file>.bak.version`; a `backup/MESA_VERSION` from an
older version of the script is only read as a fallback). A MESA file counts as pristine when
it does not contain the substring `element_sedimentation`; backups are retaken whenever the
MESA file is pristine, so stale backups from an older MESA are harmless. Before anything is
modified, `install`, `uninstall`, `reinstall` and the `single_*` commands check that every
patched file has a pristine backup recorded for the current MESA version, and the install
commands do a dry run of every anchor, so a failing command leaves MESA untouched.
Backups without a version record (made by an older version of the script) are adopted by
running `python3 setup.py soft_install --force` once. `install` re-patches from the
pristine version, so it is idempotent. `status` reports every file as up to date or STALE
relative to `code/` (or as unverifiable when its backup is missing), so it shows when MESA
needs a `soft_install` after editing the sources. Restored and copied files get the current
modification time, so builds without `./clean` recompile them, and the build scripts run
with `MESA_DIR` set to the installation being patched (also with `--mesa-dir`).

Modules other than `star` that use `star_info` (e.g. `binary`, `astero`) are only rebuilt
by the full `install`; use that if you need them.

## 7. Tests

```
cd element_sedimentation/test
export MESA_DIR=/path/to/mesa     # MESA must be built (libconst, libmath, libutils)
./mk && ./rn                      # builds ../tester, runs it, exits nonzero on failure
./clean
```

The tests link the five library modules directly against the MESA SDK libraries; nothing
is installed into MESA. They cover input validation, pure diffusion and pure advection
steady states (analytic), floor/cap behaviour, positivity, the passive counter-flow
(sum to one, species masses), time-accurate integration against the analytic two-cell
solution, the change limit, failure handling, the phase diagram interpolation on
`test/data/Schoettler+2018.csv` (nodes, midpoints, ranges, extrapolation, temperature
shift, malformed files, OpenMP), the composition arithmetic, and the energy bookkeeping
(face flow rates from a composition change, the enthalpy flux and its divergence, and
the fixed-(T, P) partial enthalpy from (rho, T) derivatives against ideal H2/He gases and
a non-ideal mock EOS), and the overturn of unstable layers (stable profiles untouched, a
single inversion, a deposit cascading to the centre, merges that propagate upward, a
face the stability test refuses; species masses and sums preserved).

## 8. Legacy methods

`'instant_rain'` and `'bottom_up_rain'` (`element_sedimentation_legacy.f90`) are the
previous algorithms, kept for comparison with older runs. They use the same phase diagram
and miscibility evaluation as the new method but transfer mass by their own rules; a
one-time notice is printed when they are selected.
