# Changelog

All notable changes to `primat` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project has
not yet reached a stable public API (`0.x`), so backwards-incompatible
changes may still land in minor bumps.

Entries here are coarse-grained by design — one line per user-visible change
or theme, not a mirror of `git log`. For full detail on any entry, `git log`
in this repository is the authoritative source.

## [Unreleased]

### Added
- `primat.sensitivity.sensitivity_table` — a one-call API returning the
  logarithmic-sensitivity matrix ∂ln(observable)/∂ln(parameter) as a
  `SensitivityTable` dataclass (`.to_markdown()`/`.to_dataframe()` views),
  with `notebooks/Sensitivity.ipynb` reduced to a thin demo and a new
  *How-to → Sensitivity tables* docs page (O-10).
- `primat --list-params` and `--version` (with backend build status) for CLI
  discoverability (S-11).
- Type hints across the public API, `py.typed` marker, and a lenient mypy CI
  job (S-10).
- `primat.__citation__` (ready-made BibTeX string) and top-level
  `run_bbn`/`run_mc`/`HAS_C_BACKEND` aliases, so common entry points no
  longer require reaching into `primat.backend` (S-9).
- `CITATION.cff`, a `paper/` JOSS submission stub, and Zenodo archival
  instructions in `PyPiGuide.md` (S-12).
- This changelog (S-13).
- Config validation with clearer error messages, shared between both
  backends (FABLEADVICE O-1).
- MC covariance/correlation output (`MCResult.cov()`/`.corr()`) and a
  `runfiles/primat_mc.py` demo script (F-1, F-2).
- Sphinx documentation site (`docs/`, published to Read the Docs), migrating
  and superseding the old `README.md`/`EXTENDING.md` prose content
  (FABLEADVICE O-3).

### Changed
- Default `Omegabh2` changed to the Planck 2018 + BAO value (`0.02242`),
  from the previous default (S-8).
- `plotly` and `joblib` moved from hard dependencies to the `mc`/`plots`/
  `gui` extras — a plain `pip install primat` install is now lighter
  (FABLEADVICE O-2).

### Fixed
- Assorted release-blocker fixes ahead of the first PyPI publish: packaging
  metadata, wheel build matrix, and related polish (S-1 through S-6).

## [0.3.2] - 2026-07-02

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
