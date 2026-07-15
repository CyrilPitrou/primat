# Changelog

All notable changes to `primat` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project has
not yet reached a stable public API (`0.x`), so backwards-incompatible
changes may still land in minor bumps.

Entries here are coarse-grained by design — one line per user-visible change
or theme, not a mirror of `git log`. For full detail on any entry, `git log`
in this repository is the authoritative source.

## [Unreleased]

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
