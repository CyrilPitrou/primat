# Changelog

All notable changes to `primat` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project has
not yet reached a stable public API (`0.x`), so backwards-incompatible
changes may still land in minor bumps.

Entries here are coarse-grained by design — one line per user-visible change
or theme, not a mirror of `git log`. For full detail on any entry, `git log`
in this repository is the authoritative source.

## [Unreleased]

### Fixed
- **QED plasma-pressure tables are now read with the same interpolant they are
  written with** (both backends). `Plasma._load_tables` documents three
  interchangeable modes (load from file / compute analytically / recompute),
  but the file path built a *linear* `interp1d` while the analytic path built a
  cubic spline. On the shipped 500-point log grid over [1e-3, 1e2] MeV
  (`d(lnT) = 0.023`, `δP ~ T⁴`) that cost ~8e-4 relative on δP, so — since
  `δP/ρ_pl ~ 4e-4` during BBN — Neff shifted in its 6th decimal according only
  to whether the cache files happened to exist. Both paths now go through
  `plasma._qed_spline` (and `cpr_cubic_spline_fit_notaknot` in C). Residual
  file-vs-analytic disagreement is 3.9e-06, set by the tables' own `%.6E` write
  precision rather than by interpolation, and pinned by
  `tests/test_qed_pressure.py::test_file_and_analytic_paths_agree`.
- **The NEVO heating function is clamped to N ≥ 0** (both backends). Heating is
  entropy flowing from the EM plasma into the neutrinos and cannot reverse, but
  74 of the 600 rows of the shipped `NEVOPRIMAT_col_1_7.csv` carry a negative
  residual from the NEVO solve itself (all within T_γ ∈ [0.0315, 0.0835] MeV,
  reaching −4.2e-06 against a peak N of 4.2e-03), which the a(T_γ) ODE
  integrated as a spurious reverse transfer.
- Together these two shift the default run by ~4e-6 relative in D/H and Li7/H
  (`D/H` 2.4358955e-05 → 2.4359049e-05 on the Python backend) — two orders of
  magnitude inside the ±3e-9 D/H regression tolerance, and cross-backend
  agreement stays at ~5e-6 against the documented 5e-5 bound.
- **`wnEDE ≤ 1/3` is now rejected** when `fEDE > 0`, on both backends. The EDE
  peak scale factor solves `u^(3(1+wnEDE)) = 4/(3·wnEDE − 1)`, which has no root
  for `wnEDE ≤ 1/3`: such a component dilutes no faster than radiation, so its
  energy *fraction* never peaks during radiation domination and `fEDE` (defined
  at that peak) is meaningless. Previously `wnEDE = 1/3` raised a bare
  `ZeroDivisionError` and `wnEDE = 0` silently produced a *complex* scale factor
  that surfaced hundreds of lines later as solve_ivp's "`y0` is complex", while
  the C backend produced a NaN background without complaint. Both are standard
  axion-like values (`wn = (n−1)/(n+1)` for n = 2 and n = 1).
- `AnalyticDistortion` no longer omits `x_of_Tg`, which made a documented
  `NeutrinoHistory` protocol attribute raise `AttributeError` instead of
  returning `None`. Latent only — `PRIMATConfig` currently forbids the flag
  combination that would reach it.
- Custom NEVO tables: overriding `nevo_spectral_file` *without* also overriding
  `nevo_grid_file` is now validated against the shipped `NEVOGrid.csv`, instead
  of being computed and then never compared — a width mismatch used to surface
  as a shape error deep inside `RegularGridInterpolator`.
- `Constants.erg` was missing a square on `second` (it read `gram·cm²/second`
  against its own docstring). Numerically inert under the natural-units
  convention, where all three base units are 1.
- `primat/qed_pressure.py`'s local `_ME_MEV` was the CODATA 2014 electron mass
  while `CONST.me` is CODATA 2018 — and the C backend already used the latter,
  so the two backends generated QED tables at different electron masses. Now
  identical (an 8e-9 relative change to freshly computed tables).
- **n↔p weak-rate cache keys** (`WEAK_RATE_FORMAT_VERSION` 1 → 4, both
  backends). Three configuration fields changed the rates but were absent from
  the fingerprint, so runs that differed only in one of them silently shared a
  cache file:
  - `munuOverTnu`/`munuOverTnu_e` was missing from the **thermal** (CCRTh)
    fingerprint, although the thermal integrands carry an explicit
    `exp(−sgnq·ξ_ν)` neutrino occupation. Degenerate-BBN runs were reusing the
    ξ=0 table — worth ~4e-3 of the base rate at ξ_e = 0.3, T = 1e10 K, i.e. far
    above anything YP tolerates — and, on a cold cache, writing their own
    ξ-specific numbers under the filename standard runs load.
  - `nevo_grid_file` was missing from the weak-rate fingerprint while its
    partner `nevo_spectral_file` was present; the two jointly define the
    tabulated distortion the SD term integrates.
  - `sampling_temperature_per_decade` was missing from the weak-rate
    fingerprint. It sets the node spacing of the linear T_ν(T_γ) interpolant
    every rate integrand reads: coarsening it moves the rates by ~1e-3 (40
    points/decade) down to ~1e-5 at the default 600.

  The version constant is also bumped past the v2/v3 generations that were
  documented in the changelog comment but never actually applied, so pre-v3
  cache files (whose `nTOp_*.txt` still included CCRTh, and whose thermal table
  was unclamped below 10^8.2 K) can no longer be loaded. The shipped tables were
  re-keyed in place — same numbers, new hash-named filenames — so default runs
  still hit the cache. **Existing `cache_dir` trees and editable installs:**
  stale files are simply never loaded again (delete them, or `primat
  --cache-clear`), but a compiled C extension built before this change computes
  the old hashes and will miss the re-keyed tables — rebuild it
  (`python setup.py build_ext --inplace`).
- No observable changes: with a rebuilt extension both backends reproduce the
  previous D/H, YP and Neff bit-for-bit.

### Documented
- **`Omeganurel` and `OneOverOmeganunr` are per neutrino flavour** (ν + ν̄), not
  summed over the three — now stated in `Background.Omeganuh2_relnu`/`_nrnu`,
  `docs/howto/output.md` and `primat-c/include/background.h`. The values are
  unchanged (multiply `Omeganurel` by 3 for the usual quoted total ≈ 17); the
  per-flavour convention is the natural one for `OneOverOmeganunr`, whose ≈ 93
  reproduces the standard Σm_ν / 93.1 eV normalisation.
- **The ΔNeff extra species deliberately uses a different "instantaneous
  decoupling" normalisation from the SM neutrinos** when `QED_corrections=True`:
  `T_nu_decoupling`'s free-gas σ_∞ rather than the QED-corrected `_sbar_ref`,
  leaving it ~0.31 % low in energy density during BBN. That is what makes the
  reported `Neff` come out as `Neff_SM + ΔNeff` to machine precision, i.e. makes
  the knob mean what its name says; `rho_nu_extra` now spells out the trade-off
  and the measured cost of the alternative, and `T_nu_decoupling` no longer
  claims to be the SM neutrino temperature in that mode.
- `plasma.rho_SM`/`p_SM` are labelled as the ξ = 0, no-spectral-distortion
  reference quantities they are — **not** the Friedmann source, which is
  `StandardBackground.Hubble` (it adds each flavour's own ξ and `rho_nu_SD`).
- `Background.t_of_T` documents its valid range: outside
  `[T_end, T_start_cosmo]` it extrapolates linearly and can return a negative
  time. This is distinct from the radiation-domination extrapolation below the
  NEVO table's edge, which is inside the solved span and correct in both
  `external_scale_factor` modes.
- The `external_scale_factor` True/False agreement is quoted at its measured
  ~1e-5 (per-observable figures given) instead of "~1e-6".
- Every constant in `primat/constants.py` now records its edition (SI 2019
  exact / CODATA 2018 / CODATA 2010 / PDG 2020 / PDG 2018 / AME2020 / AME2016 /
  Fixsen 2009); the set is deliberately not single-vintage, and `gA`/`Vud` are
  flagged as the group that feeds the n↔p rates.
- Assorted `plasma.py` docstring drift corrected: a reference to a
  `_setup_qed_pressure` method that does not exist, a recompute-mode paragraph
  naming the wrong output file, stale `rates/plasma/` paths, and a "no
  module-level mutable state" claim that overlooked the four numba integrand
  handles (harmless, but now described accurately).
- `primat/weak_rates/corrections.py` now cites the Phys. Rep. **equation**
  numbers for every correction term, instead of section numbers alone: the
  relativistic Fermi function (Eq. 100), the resummed radiative factor
  (Eq. B35, with g = B32, Spence L = B33, constants B31/B36), the finite-mass
  terms (Eqs. 114, 115a/115b, χ_FM from App. B.3) and the four CCRTh
  sub-integrands (Eqs. 107, 109, 112a, 113, with the F_± kernels of Eq. B41 and
  the B51a/B51b kernel). Three section references were also off against
  `biblio/Pitrou_etal_PhysReptArxivVersion.pdf` and are corrected: the T=0
  radiative corrections are §III.E (was §III.D) and the finite-temperature ones
  §III.F (was §III.H, which is "Weak magnetism").
- The `F_+`/`F_−` asymmetry in the bremsstrahlung soft subtraction is now
  documented as deliberate on both backends, citing Phys. Rep. Eq. B43 where it
  is printed, plus the measured consequence of "correcting" it (the CCRTh sum
  would grow to ~0.8% of the base rate at 3e10 K).
- Two accepted-but-not-self-consistent flag combinations are called out in
  `config.py`: `thermal_corrections=True` with `radiative_corrections=False`,
  and the absence of the SD-FM term outside analytic-distortion mode.
- `background=`: documented that the weak-rate cache is keyed on the config
  alone and cannot see a custom background's temperature grid — use
  `weak_rate_cache=False`/`save_nTOp=False` for a non-standard history.

### Changed
- Python backend: the ten weak-rate Fermi-Dirac integrand kernels
  (`weak_rates/integrands.py`) and the four e± electron-thermo integrands
  (`plasma.py`) are now numba-compiled with `cache=True`, so a fresh process
  (joblib MC worker, Streamlit server, re-run CLI) loads the compiled machine
  code from numba's on-disk cache instead of recompiling — ~2.3 s of cold-start
  JIT saved per process (measured). The plasma integrands were moved to module
  level and now take the electron mass `me` as an explicit argument (rather than
  closing over `cfg.me`), which is what makes the on-disk cache safe. No effect
  on any observable.
- C backend (`primat-c`) now routes its ~250 unrecoverable heap allocations
  (ODE work vectors, spline tables, loaded network, result arrays, …) through
  new checked helpers `cpr_xmalloc`/`cpr_xcalloc`/`cpr_xrealloc`
  (`include/xalloc.h`): a failed allocation now prints
  `primat: out of memory (<bytes>) at <file>:<line>` and exits, instead of
  dereferencing NULL and crashing anonymously. Sites that intentionally
  degrade gracefully (cache writers with their own NULL checks, `errmsg`-return
  paths) are unchanged. No effect on any observable — purely OOM-diagnostic
  robustness.

### Removed
- `rate_interp_order` `DEFAULT_PARAMS`/C config parameter: never consumed by
  any solver, resampler, or rate-lookup path on either backend (rate-table
  resampling always hardcoded log-log cubic, and the per-step master-grid
  lookup always used linear `searchsorted` regardless of its value), so
  setting it to `quadratic`/`cubic` changed neither observables nor runtime.
  Removed together with its C field/default/validation, both param
  templates' entries, both backends' enum tests, and the now-empty
  `_PARAM_CHOICES` machinery it was the last user of.

## [0.3.2] - 2026-07-11

### Added
- `output_rates_time_evolution` now writes per-reaction forward-rate columns
  (`<reaction>_frwrd`, small/small_parthenope networks), on both backends —
  previously a no-op. Populated in `EvolutionResult.rates` and round-tripped
  by `primat.evolution.load_evolution`.
- `primat.sensitivity.sensitivity_table` — a one-call API returning the
  logarithmic-sensitivity matrix ∂ln(observable)/∂ln(parameter) as a
  `SensitivityTable` dataclass (`.to_markdown()`/`.to_dataframe()` views),
  with `notebooks/Sensitivity.ipynb` reduced to a thin demo and a new
  *How-to → Sensitivity tables* docs page.
- `primat --list-params` and `--version` (with backend build status) for CLI
  discoverability.
- Type hints across the public API, `py.typed` marker, and a lenient mypy CI
  job.
- `primat.__citation__` (ready-made BibTeX string) and top-level
  `run_bbn`/`run_mc`/`HAS_C_BACKEND` aliases, so common entry points no
  longer require reaching into `primat.backend`.
- `CITATION.cff` and Zenodo archival instructions in `PyPiGuide.md`.
- This changelog.
- Config validation with clearer error messages, shared between both
  backends.
- MC covariance/correlation output (`MCResult.cov()`/`.corr()`) and a
  `runfiles/primat_mc.py` demo script.
- Sphinx documentation site (`docs/`, published to Read the Docs), migrating
  and superseding the old `README.md`/`EXTENDING.md` prose content.
- `primat-gui`'s Final abundances tab gained a single reproduction-bundle
  download (`.py`/`.ini` + README, `primat/gui/export_params.py`): prints the
  full standard-ratio `run_bbn` centrals plus a `run_mc(seed=0)` std-only
  block, pins `force_backend` to whichever backend actually ran, and for a
  custom network embeds the exact `custom_network` dict (Python) or a
  `nuclear/` overlay directory (`.ini`, via `user_nuclear_dir`) — including
  any uploaded/edited rate tables — so a downloaded bundle reproduces the
  GUI run bit-exactly.
- `notebooks/ReactionRates.ipynb` (⟨σv⟩(T9) of any reaction vs. the Hubble
  rate, with the master-grid reinterpolation overlaid) and
  `notebooks/AnimatedAbundances.ipynb` (animated GIFs of the small-network
  abundance evolution vs. ΔNeff/Ω_b h²), both wired into the docs tutorial
  gallery with a guard test against future gallery drift.

### Changed
- The C backend now supports `extra_rho` and `decay_era`, closing two of the
  three former Python-only feature gaps. `extra_rho` callables are
  sampled onto a dense temperature grid and splined into the C Friedmann
  equation; `decay_era`'s long-lived-isotope Decay-Time propagation is ported
  via a scaling-and-squaring Padé matrix exponential, writing an identical
  `output_decay_evolution` TSV. Only `background=` (a custom `Background`
  object) remains Python-only.
- Default `Omegabh2` changed to the Planck 2018 + BAO value (`0.02242`),
  from the previous default.
- `plotly` and `joblib` moved from hard dependencies to the `mc`/`plots`/
  `gui` extras — a plain `pip install primat` install is now lighter.
- Nuclear rate tables (`large` network) regenerated at 1000 points instead
  of being reinterpolated onto the master T9 grid from a coarser source
  grid, for both backends.
- Cache trees consolidated under `cache_plasma_weak/{weak,plasma}/` with an
  additive `cache_dir` overlay redirect for read-only installs, and
  non-fatal (warn, not crash) cache-write failures on both backends
  (80 → 81 `DEFAULT_PARAMS` keys).
- Per-flavour neutrino degeneracies `munuOverTnu_e/mu/tau`.

### Fixed
- Assorted release-blocker fixes ahead of the first PyPI publish: packaging
  metadata, wheel build matrix, and related polish.
- The Python backend's `run_bbn` now also exposes a `Y_final` sub-dict,
  matching the C backend and restoring result-dict parity (CLAUDE.md).
- Windows editable-install C-extension shadowing, and assorted
  Windows-portability failures in the CI Tests matrix.
- MSVC POSIX-header/pthreads build failures on the Windows leg of
  `wheels.yml`, so Windows users now get binary wheels (with the fast C
  backend) from PyPI instead of falling back to a source build (64-bit only).
- `GN` (Newton's constant) default corrected to the exact CODATA literal
  `6.6743e-11`; the previous default (`6.674299257609439e-11`) was off at
  the ~1.1e-7 relative level.

### Documented
- Electron-thermo / QED-pressure table extrapolation behaviour beyond the
  tabulated range, with a runtime warning when a run goes past it.

## [0.3.1] - 2026-07-02

### Added
- `show_progress` config flag to control `[primat]`/`[MC]` stderr progress
  messages, wired through both the Python and C backends' CLIs
  (including `--flag`/`--no-flag` boolean parsing).
- `mc_rate_rescale_cap` parameter to cap Monte Carlo rate-rescaling factors
  (default lowered from `1e3` to `30` after further validation).
- Flat `sigma_<name>` fields in MC results on both backends.
- Ctrl-C abort support for a running `primat-c` Monte Carlo (`run_mc`)
  sample.
- Memory-leak checking (`make leak-test`) and an ASan/UBSan CI job for the
  C backend.
- A pre-computed `PRIMAT_Yp_DH_ErrorMC_1000_2026.dat` table for CLASS/CAMB
  consumption.

### Changed
- Renamed `rates_dir`/`user_rates_dir` config fields to `data_dir`/
  `user_nuclear_dir` (clearer overlay semantics — see `CLAUDE.md`'s "Rates
  directory resolution" section).
- Split the three separate QED plasma-pressure correction table files into
  one consolidated `QED_tables.txt`.
- Integrated three background-ODE performance branches: dense-output RK45,
  a combined 2D background ODE, and monotone spline lookups for rate
  interpolation.
- Smoothed MC progress reporting (was jumping straight from 0% to 100%).
- Moved C backend headers from `include/cprimat/` to `include/`.

### Fixed
- `--no-show_progress` being silently ignored by the `primat-c` CLI's
  `--mc` path.
- `delta_<rxn>` rate perturbations not applying when
  `rescale_nuclear_rates` was left at its default.
- A NaN issue in the C backend's `electron_thermo` plasma cache.
- Several gcc-14 warnings (`-Wformat-truncation`, `-Wmaybe-uninitialized`,
  unused-parameter) in the C backend.

## [0.3.0] - 2026-06-25

Initial PyPI-track release. By this point the project already had its
current two-backend architecture in place:

### Added
- Dual backend: a pure-Python implementation (`primat/`) and a fast C99
  port (`primat-c/`), exposed to Python via a compiled extension and
  dispatched through `primat.backend.run_bbn(force_backend={"auto","c","python"})`.
- `primat-gui`, a Streamlit application (four usage modes: Python API, CLI,
  GUI, notebooks), including a "Customise Reactions" flow to build and
  import/export custom nuclear networks.
- Monte Carlo uncertainty propagation (`run_mc`) with rate-key resolution
  and incremental sample reuse (`prev`) on the Python backend.
- Non-instantaneous neutrino decoupling via the NEVO tables, with
  overridable `nevo_file`/`nevo_spectral_file`/`nevo_grid_file`/
  `nevo_file_prefix` parameters.
- Analytic QED plasma-pressure corrections and per-reaction weak-rate
  correction terms (Born, CCR, finite-mass, thermal, spectral-distortion).
- Grid-agnostic nuclear rate loading (tables resampled onto a configurable
  master T9 grid at load time) and the `large`/`amax` network-filtering
  mechanism.
- Unified log-log cubic not-a-knot interpolation for n↔p weak rates on both
  backends, collapsing the former C-vs-Python D/H gap to ≲1e-5.
- Unified time-evolution TSV schema (`primat/evolution.py`,
  `EvolutionResult`/`load_evolution`), implemented identically by both
  backends.

[Unreleased]: https://github.com/CyrilPitrou/primat/compare/v0.3.2...HEAD
[0.3.2]: https://github.com/CyrilPitrou/primat/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/CyrilPitrou/primat/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/CyrilPitrou/primat/releases/tag/v0.3.0
