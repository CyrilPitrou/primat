# primat Test Suite

This test suite is meant to be read as documentation as much as run as a check:
every test group below states the physics or software property it pins down.

## Why these tests exist

primat is a numerical physics code whose output (primordial abundances) is a
single set of numbers that is easy to get *subtly* wrong — a mistyped Jacobian
entry, a rate read from the wrong column, a refactor that shifts a result by a
few parts in 1e5. The suite is built in layers so that:

1. **fast unit tests** catch gross breakage in seconds (config, plasma
   thermodynamics, the public API, the Monte-Carlo machinery, weak-rate
   helpers);
2. **structural cross-checks** prove the compiled, stoichiometry-driven network
   kernels agree with the declarative reaction table to machine precision (this
   is what caught four latent Jacobian bugs — see `test_network_builder.py`);
3. **regression tests** pin the final abundances, both loosely at default
   precision and tightly (the `reference` marker) at the exact settings that
   produced the published CLAUDE.md values.

## Running the tests

From the repository root:

```bash
pytest tests/                          # everything (~4 min)
pytest tests/ -m "not slow"            # fast lane: config/plasma/structural unit tests, <10 s
pytest tests/ -m "not slow or solve"   # fast lane + default-precision solves (CI: every push/PR, .github/workflows/tests.yml)
pytest tests/ -m "not reference"       # skip only the ~1 min high-precision reference runs
pytest tests/ -m reference             # only the tight CLAUDE.md regression (CI: nightly, .github/workflows/tests.yml)
pytest tests/test_plasma.py -v         # a single file, verbose
```

## Markers

Three speed tiers:

| Marker | Meaning |
|--------|---------|
| `slow` | any test excluded from the fast lane: a full primat solve (or a Monte-Carlo loop of solves), a weak-rate recompute (~1.8 s, bypassing the fingerprinted cache), or packaging checks. Deselect with `-m "not slow"`. |
| `solve` | the "solve" tier: tests that run >=1 full primat solve at *default* (non-reference) precision; always also marked `slow`. `-m "not slow or solve"` selects the fast lane plus this tier. |
| `reference` | high-precision runs (numerical_precision=1e-10, sampling_temperature_per_decade=2000, sampling_nTOp_per_decade=125, T_start_cosmo=100 MeV) that reproduce the documented reference values to YP ±1e-5, D/H ±3e-9; ~60 s total; always also marked `slow`. |
| `wheel` | builds a wheel and `pip install`s it into a clean venv before running a smoke solve; always also marked `slow`. |
| `gui` | drives the optional Streamlit GUI (`primat.gui`) via `AppTest`; skipped if the `gui` extra is not installed; always also marked `slow` and `solve`. |

`tests/test_gui.py` (`gui` marker) is skipped automatically unless the
optional `gui` extra is installed (`pip install -e ".[gui]"`); install it to
also exercise the Streamlit GUI end-to-end.

The fast lane (`-m "not slow"`) does include *one* cheap solve: the
`solved_small` session fixture (`conftest.py`), used by most of
`test_api.py`. It uses the default config (`weak_rate_cache=True`), so it
loads the n<->p rates from the fingerprinted cache instead of recomputing
them (~1 s total for `__init__` + `solve()`). Anything that needs *more* than
this single default-precision solve -- a second solve with different flags, a
Monte-Carlo loop, etc. -- is tagged `solve` (and `solved_large`, used only by
`test_regression.py`, is entirely in the `slow`/`solve`/`reference` tiers).

**Deferred**: replacing `solve`-tier tests with era-level tests (e.g.
seeding the LT era directly via Saha, instead of integrating the HT+MT eras
first) would shrink the `solve` tier further, but needs the era integrations
to be exposed as callable units first (a future architecture phase) -- not
yet done, so the `solve` tier still runs full three-era solves.

## Structure

| File | What it checks |
|------|----------------|
| `conftest.py` | Session-scoped fixtures: pre-solved small- and large-network `PRIMAT` instances reused across tests (built once, not per test). |
| `test_config.py` | `PyPRConfig`: defaults, user overrides, unknown-key warnings, p_*/delta_* reaction-name typo warnings, the `Nuclides` table, that `eta0b` tracks `Omegabh2`, and that there is exactly one MCMC weight per network reaction. |
| `test_constants.py` | `primat.constants.CONST`'s derived electroweak values: `sW2` (sin²θ_W) against an independent hand-computation of the on-shell relation, the `geL`/`geR`/`gmuL` effective couplings derived from it, and `T_weak`/`T_nucl` against `MeV_to_Kelvin`. |
| `test_plasma.py` | Plasma/neutrino thermodynamics: `rho_g`, `rho_e`/`p_e` positivity and the e± cutoff, `spl`/`dspl_dT` self-consistency (combined vs separate evaluation, vs finite differences), `T_nu_decoupling` high- and low-T limits. |
| `test_decoupling_qed.py` | The `incomplete_decoupling` × `QED_corrections` 2×2 flag matrix: that `PofT`/`dPdT`/`d2PdT2` vanish when `QED_corrections=False`; that `spl/T³` equals `11π²/45` (free-gas) or differs from it (QED) at high T; that the instantaneous-decoupling $(T_\gamma/T_\nu)^3$ ratio equals `11/4` without QED and the Dodelson–Turner–Heckler perturbative formula with QED; that the correct NEVO file is loaded for each combination; and Neff reference values pinned for all four combinations. |
| `test_api.py` | Public API: `A/N/Z` dicts, `__getitem__` abundance interpolators (scalar and array input, non-negativity), `get_quantity`, lazy `solve()`, `T_of_t`/`t_of_T`. |
| `test_mc.py` | Monte-Carlo machinery: `MCResult`/`MCQuantityResult` shapes and attributes, mean/std consistency, reproducibility for a fixed seed, that varying rates gives `std > 0`. |
| `test_weak_rates.py` | n↔p weak rates: Fermi-Dirac helpers, the Fermi-Coulomb correction, the neutron-decay phase-space integral `ComputeFn`, that the two loaded rate interpolants (forward/backward) are positive and obey detailed balance (ratio → 1 at high T); (`slow`-tier) that `RecomputeWeakRates`'s recompute path (`weak_rate_cache=False`) agrees with the fingerprinted cache it normally loads; and that `_setup_fd_impls` re-wraps the module-level `FD_*` implementations (jitted vs pure-Python) every time `numba_installed` changes, rather than latching on the first call. |
| `test_cache_utils.py` | The fingerprinted-cache helpers (`primat.cache_utils`): `fingerprint_hash` is order-independent and value-sensitive; write/read round-trips; a missing file or a header-less/corrupt file is reported as unknown fingerprint (`None`) rather than raising; and that `write_cache_with_fingerprint` writes atomically (temp file + `os.replace`, no leftover `.tmp.<pid>` file, safe to overwrite an existing cache). |
| `test_refactor_invariants.py` | Properties the performance refactor relies on: MC results independent of `n_jobs`, `eta0b` recomputed on reassignment, GN/tau_n overridable, electron-thermo tabulation ≈ exact integrals, `_LinearRate` ≈ `interp1d(kind='linear')`. |
| `test_custom_loader.py` | The `small_parthenope` custom network file: verifies the reaction set and species match the standard small network; that reactions routed to non-default files (e.g. `ddTOHe3n_parthenope.txt`) actually use different rate values; that the loaded network passes N/Z/Q conservation; and that a full BBN solve gives physically reasonable YP and D/H. |
| `test_qed_pressure.py` | The analytical QED plasma-pressure module (`primat.qed_pressure`): Fermi-Dirac integral analytic limits (UR limit → π²/12, non-relativistic cutoff), sign conventions (δP_a < 0, δP_{e3} > 0), agreement with the PRIMAT-generated tables to 0.5% at T ≥ 2 MeV, numerical derivative consistency, and a save/load round-trip check. |
| `test_network_generation.py` | The offline generation layer (`generate_rates/`): token resolution, the formal baryon/charge conservation check, that `nuclides.csv` agrees with the hard-coded table, that the deduced reaction list is a superset of the 12-/62-reaction networks, and that the computed detailed-balance coefficients reproduce the published values. |
| `test_network_builder.py` | The generic stoichiometry-driven kernels (`network_builder`), the single network path: compiled RHS/Jacobian equal the `reactions` reference to machine precision; the formal N/Z conservation check (passes for real nets, fires on a broken one); numerical baryon-number conservation; the full `UpdateNuclearRates` driver methods (rhs/rhsMT/rhsLT + Jacobians); era-independent table invariants (buffer-order lengths, per-reaction A/Z conservation); and the `amax` mass-cutoff filter (correct count, nuclide bound, conservation, and invalid-value rejection). |
| `test_large_network.py` | The large network: it loads (~59 nuclides, ~433 reactions) and passes the formal conservation check; the vectorised rate buffer stays finite/bounded across the LT range; and a full solve conserves baryon number exactly while matching the medium network on the light elements. |
| `test_nuclear_qed.py` | QED corrections to radiative-capture rates (Pitrou & Pospelov 2020): correction factors are > 1 and sub-percent; the npTOdg polynomial matches its T9→0 cap; the four Kroll-formula reactions increase monotonically with T9; reference magnitudes at T9=0.1 GK are pinned to ±2e-6; non-QED reactions are unchanged; p_* variations stack correctly on the corrected median; and a full solve with the flag on shifts D/H by a detectable but sub-percent amount. |
| `test_regression.py` | Final abundances: loose default-precision sanity checks, tight `reference`-marked checks against the published CLAUDE.md values, no-numba full solve checks (pure-Python kernels must match JIT to 1e-4), and the `amax` cutoff verification (large network filtered to A ≤ 20 matches medium light elements to ~1e-3). |
| `test_wheel_smoke.py` | The `wheel`-marked "pip install" smoke test: builds a wheel, installs it into a clean venv, and runs a small-network solve there to catch package-data/path regressions (e.g. `rates/` not shipped, or a path computed relative to the source tree instead of the installed package) that an editable install would not reveal. |
| `test_docs_consistency.py` | Guards README claims that aren't checked anywhere else: `PRIMATConfig`'s `save_nTOp`/`save_nTOp_thermal` defaults match what README states, and the parameter names/values documented for `runfiles/primat_reference_run.py` (`sampling_temperature_per_decade`, `numerical_precision`, `sampling_nTOp_per_decade`, `T_start_cosmo_MeV`) still exist verbatim in that script and are recognised by `PRIMATConfig` (no "unknown parameter" warning). |
| `test_gui.py` | The optional Streamlit GUI (`primat.gui`): `import primat` does not pull in `primat.gui`/streamlit; the parameter-form metadata covers `amax` (the one `None`-default key) and the network choices; an end-to-end `AppTest` run of `primat/gui/app.py` reproduces `test_cli.py`'s pinned default-run values (Neff/YPBBN/D-H and the per-nuclide table) -- i.e. the GUI drives `PRIMAT` identically to the CLI; the abundance-evolution panel renders with its default "light elements" nuclide selection; the `amax` widget appears only for `network='large'`; and an invalid flag combination (`spectral_distortions=True` with `incomplete_decoupling=False`) is shown as a clean `st.error` rather than a traceback. Skipped entirely if the `gui` extra is not installed. |

| `test_backend_parity.py` | Backend parity: `primat._primat_c` (C) vs `primat.main.PRIMAT` (Python) — the result-dict shape and numerical agreement across the two backends. |
| `test_background.py` | Direct unit tests for `primat.background.StandardBackground`. |
| `test_custom_background.py` | The `custom_background` mode: a user-supplied (T, t, a) table drives the cosmological background while the nuclear network is solved with instantaneous-decoupling n↔p weak rates. |
| `test_decay_rates.py` | Radioactive-decay treatment in the `large` network. |
| `test_detailed_balance.py` | `compute_detailed_balance_coefficients` reproduces the values in `detailed_balance.csv`. |
| `test_deuterium_network.py` | The `network="large", amax=2` configuration reproduces the old standalone `deuterium` network (single-reaction `n_p__d_g`), the clean-slate starting point for custom networks. |
| `test_evolution.py` | Unit tests for `primat.evolution` (the unified time-evolution schema), independent of any BBN solve: round-trips a synthetic `EvolutionResult` through `dump_evolution`/`load_evolution`, both to a string and to a file on disk. |
| `test_gui_custom_network.py` | The GUI's "Manage networks" dialog and the "Create custom network" dialog it gates (`primat.gui.params_form`). |
| `test_gui_run_view.py` | Direct parity test for `primat.gui.run_view.GuiRun`. |
| `test_neutrino_history.py` | The pluggable neutrino-sector background. |
| `test_notebooks.py` | Notebook smoke test: papermill-executes the example notebooks. |
| `test_nuclear.py` | The auto-derivation fallback in `reaction_stoichiometry` and the duplicate-entry check in `load_network`. |
| `test_rate_variations.py` | Nuclear rate variation and MC uncertainty propagation. |
| `test_runfiles.py` | Smoke test for the example scripts in `runfiles/`. |
| `test_sensitivity.py` | `primat.sensitivity` — the logarithmic-sensitivity API. |
| `test_spectral_distortions.py` | Two full solves with `spectral_distortions` on/off: a small but non-zero effect on D/H, and zero distortion energy density in NEVO by construction (`solve` tier). |
| `reference_values.py` | (helper, not a test) Centralised default-run reference observables shared by test_cli/test_gui/test_regression, and the validation-reference constants (single source, see test_docs_consistency). |
| `_oracles.py` | (helper, not a test) Test-only reference RHS/Jacobian oracle implementations the nuclear-network tests compare against. |

## Validation reference (authoritative copy)

After any modification, run `python runfiles/primat_run.py` from the repo
root and check the output against these tables. A result outside these
bounds indicates a regression. (This section moved here from the untracked
CLAUDE.md so that CI and public clones carry it; tests parse THESE tables
— see tests/test_docs_consistency.py.)

The following values must hold (`Omegabh2=0.02242`, `spectral_distortions=True`
and `nuclear_qed_corrections=True`, both defaults (the latter from the
radiative-capture QED corrections of Pitrou & Pospelov 2020).
References were produced by `runfiles/primat_reference_run.py` (high-precision
run, `numerical_precision=1e-10`, `sampling_temperature_per_decade=2000`,
`sampling_nTOp_per_decade=125`, `T_start_cosmo_MeV=100`, `rate_grid_npts=4000`
explicit so this reference stays decoupled from the routine-run default below).

**Small network** (`network="small"`):

| Observable | Expected | Tolerance |
|------------|----------|-----------|
| YP (BBN) | 0.24699819 | ±1e-5 |
| D/H | 2.43588e-5 | ±3e-9 |

**`large, amax=8`** (the old "medium" network's exact 68-reaction equivalent):

| Observable | Expected | Tolerance |
|------------|----------|-----------|
| YP (BBN) | 0.24700154 | ±1e-5 |
| D/H | 2.43659e-5 | ±3e-9 |

A result outside these bounds indicates a regression.

### Per-nuclide final abundances (small / large+amax=8 / large)

Final mass-fraction abundances `Y` of the small-network nuclides at the end of
BBN, from the default run (`numerical_precision=1e-7`, generic kernels). The
full large network must match `large, amax=8` (the old "medium" network's
exact equivalent) on every light element to ≲1e-3, and now also on the free
neutron `n` (≲1e-3): the reverse-rate clamp in `primat/network_data.py`
(see its module docstring, "exothermic blow-up") removed the spurious low-T
flooding of heavy nuclides such as B10, which previously fed β-delayed
neutron emission and inflated `n` to ~7e-13. With the clamp the heavy tail is
negligible and `n` tracks the `amax=8` value.

| Nuclide | small | large, amax=8 | large |
|---------|-------|----------------|-------|
| n   | 3.997435e-16 | 3.996517e-16 | 3.996503e-16 |
| p   | 7.529426e-01 | 7.529393e-01 | 7.529393e-01 |
| H2  | 1.834099e-05 | 1.834603e-05 | 1.834601e-05 |
| H3  | 5.852032e-08 | 5.839099e-08 | 5.839082e-08 |
| He4 | 6.174931e-02 | 6.175013e-02 | 6.175013e-02 |
| Li7 | 2.181386e-11 | 9.178083e-11 | 9.178037e-11 |
| Be7 | 3.966369e-10 | 3.223616e-10 | 3.223599e-10 |

Baryon number is conserved exactly by every network (`sum_s A_s Y_s = 1` to ~1e-10).

**`small`-only exception (H3/Li7/Be7).** Commit `6221e43` added the
`tTOHe3Bm` and `Be7TOLi7Bp` analytic beta-decay/electron-capture reactions to
`primat/data/nuclear/networks/large.txt` (present at every `amax` that keeps A ≤ 8,
including `amax=8` itself) but not to `small`'s hard-coded `ORDER_SMALL`.
Their rates are the *laboratory* decay constants (`ln2/τ_{1/2}`,
τ_{1/2}(H3) = 12.32 yr, τ_{1/2}(Be7) = 53.29 d) and act over the full
integration window down to `T_end = 0.001 MeV` (`t_end ≈ 1.3e6 s ≈ 15 days`,
see `nuclear_network.py`'s `_setup_evolution`), so by `t_end` they convert a
non-negligible fraction of H3 -> He3 (~0.23%) and Be7 -> Li7 (~18%) in both
the `large, amax=8` and full `large` columns above (hence both differing
sharply from `small` on these three rows, while agreeing closely with each
other). `tests/test_large_network.py` nonetheless excludes H3/Li7/Be7 from
its large-vs-amax=8 ≲1e-3 assertions, out of caution against this shared
decay-channel sensitivity.
