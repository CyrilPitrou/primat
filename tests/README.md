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
   produced the published values in "Validation reference" below.

## Running the tests

From the repository root:

```bash
pytest tests/                          # everything (~12 min)
pytest tests/ -m "not slow"            # fast lane: config/plasma/structural unit tests, ~30 s
pytest tests/ -m "not slow or solve"   # fast lane + default-precision solves (CI: every push/PR, .github/workflows/tests.yml)
pytest tests/ -m "not reference"       # skip only the ~1 min high-precision reference runs
pytest tests/ -m reference             # only the tight high-precision regression
pytest tests/test_plasma.py -v         # a single file, verbose
```

CI (`.github/workflows/tests.yml`) runs `-m "not slow or solve"` on every
push/PR across an OS × Python matrix, and the **full** suite (`pytest tests/`,
i.e. including the `reference`, `wheel` and `notebook` tiers) nightly.

## Markers

Every marker below is registered in `pytest.ini`; the list here and that one
must stay in step (`pytest --strict-markers` is the check).

| Marker | Meaning |
|--------|---------|
| `slow` | any test excluded from the fast lane: a full primat solve (or a Monte-Carlo loop of solves), a weak-rate recompute (~1.8 s, bypassing the fingerprinted cache), or packaging checks. Deselect with `-m "not slow"`. |
| `solve` | the "solve" tier: tests that run >=1 full primat solve at *default* (non-reference) precision; always also marked `slow`. `-m "not slow or solve"` selects the fast lane plus this tier. |
| `reference` | high-precision runs (numerical_precision=1e-10, sampling_temperature_per_decade=2000, sampling_nTOp_per_decade=125, T_start_cosmo=100 MeV) that reproduce the documented reference values to YP ±1e-5, D/H ±3e-9; ~60 s total; always also marked `slow`. |
| `wheel` | builds a wheel and `pip install`s it into a clean venv before running a smoke solve; always also marked `slow`. |
| `gui` | drives the optional Streamlit GUI (`primat.gui`) via `AppTest`; skipped if the `gui` extra is not installed; always also marked `slow` and `solve`. |
| `notebook` | papermill-executes a demonstration notebook end-to-end; skipped if the optional `notebooks` extra is not installed; always also marked `slow`. |
| `backend` | compares the `primat._primat_c` C-extension backend against the pure-Python one; skipped if the C extension is not built; always also marked `slow` and `solve`. |

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
| `test_config.py` | `PRIMATConfig`: defaults, user overrides, unknown-key warnings, p_*/delta_* reaction-name typo warnings, the `Nuclides` table, that `eta0b` tracks `Omegabh2`, and that there is exactly one MCMC weight per network reaction. |
| `test_constants.py` | `primat.constants.CONST`'s derived electroweak values: `sW2` (sin²θ_W) against an independent hand-computation of the on-shell relation, the `geL`/`geR`/`gmuL` effective couplings derived from it, and `T_weak`/`T_nucl` against `MeV_to_Kelvin`. |
| `test_plasma.py` | Plasma/neutrino thermodynamics: `rho_g`, `rho_e`/`p_e` positivity and the e± cutoff, `spl`/`dspl_dT` self-consistency (combined vs separate evaluation, vs finite differences), `T_nu_decoupling` high- and low-T limits. |
| `test_decoupling_qed.py` | The `incomplete_decoupling` × `QED_corrections` 2×2 flag matrix: that `PofT`/`dPdT`/`d2PdT2` vanish when `QED_corrections=False`; that `spl/T³` equals `11π²/45` (free-gas) or differs from it (QED) at high T; that the instantaneous-decoupling $(T_\gamma/T_\nu)^3$ ratio equals `11/4` without QED and the Dodelson–Turner–Heckler perturbative formula with QED; that the correct NEVO file is loaded for each combination; and Neff reference values pinned for all four combinations. |
| `test_api.py` | Public API: `A/N/Z` dicts, `__getitem__` abundance interpolators (scalar and array input, non-negativity), `get_quantity`, lazy `solve()`, `T_of_t`/`t_of_T`; and solver-failure reporting — `nuclear_network._check_solver` raises `RuntimeError` (quoting the era and scipy's own `sol.message`) rather than letting an unconverged `solve_ivp` return silently-wrong abundances, verified by injecting a failure into each of the HT/MT/LT eras in turn. |
| `test_mc.py` | Monte-Carlo machinery: `MCResult`/`MCQuantityResult` shapes and attributes, mean/std consistency, reproducibility for a fixed seed, that varying rates gives `std > 0`. |
| `test_weak_rates.py` | n↔p weak rates: Fermi-Dirac helpers, the Fermi-Coulomb correction, the neutron-decay phase-space integral `ComputeFn`, that the two loaded rate interpolants (forward/backward) are positive and obey detailed balance (ratio → 1 at high T); (`slow`-tier) that `RecomputeWeakRates`'s recompute path (`weak_rate_cache=False`) reproduces the shipped, git-tracked cache to 1e-6 relative — the only pin on the *contents* of the n↔p rate tables, and so on the individual correction terms; and that `_setup_fd_impls` re-wraps the module-level `FD_*` implementations (jitted vs pure-Python) every time `numba_installed` changes, rather than latching on the first call. |
| `test_cache_utils.py` | The fingerprinted-cache helpers (`primat.cache_utils`): `fingerprint_hash` is order-independent and value-sensitive; write/read round-trips; a missing file or a header-less/corrupt file is reported as unknown fingerprint (`None`) rather than raising; and that `write_cache_with_fingerprint` writes atomically (temp file + `os.replace`, no leftover `.tmp.<pid>` file, safe to overwrite an existing cache). |
| `test_refactor_invariants.py` | Properties the performance refactor relies on: MC results independent of `n_jobs`, `eta0b` recomputed on reassignment, GN/tau_n overridable, electron-thermo tabulation ≈ exact integrals, `_LinearRate` ≈ `interp1d(kind='linear')`. |
| `test_custom_loader.py` | The `small_parthenope` custom network file: verifies the reaction set and species match the standard small network; that reactions routed to non-default files (e.g. `ddTOHe3n_parthenope.txt`) actually use different rate values; that the loaded network passes N/Z/Q conservation; and that a full BBN solve gives physically reasonable YP and D/H. |
| `test_qed_pressure.py` | The analytical QED plasma-pressure module (`primat.qed_pressure`): Fermi-Dirac integral analytic limits (UR limit → π²/12, non-relativistic cutoff), sign conventions (δP_a < 0, δP_{e3} > 0), agreement with the PRIMAT-generated tables to 0.5% at T ≥ 2 MeV, numerical derivative consistency, and a save/load round-trip check. |
| `test_network_generation.py` | The offline generation layer (`generate_rates/`): token resolution, the formal baryon/charge conservation check, that `nuclides.csv` agrees with the hard-coded table, that the deduced reaction list is a superset of the 12- and 68-reaction networks, that the computed detailed-balance coefficients reproduce the published values, and that the NUBASE half-life field is read at its documented (1-based) column offsets — limits included, since a bound is not a measurement. Plus (`slow` tier) the two generation *commands* end to end: `convert_ac2024_rates.py`, run from a throwaway cwd, must regenerate all 395 shipped artifacts (390 rate tables + `decays.txt` + the three CSVs + `networks/large.txt`) **byte for byte** and emit no warning, and `generate_qed_tables.py` must still run and write fingerprinted tables — into the directory the solver actually loads, a path that had silently gone stale before anything executed the script. |
| `test_network_builder.py` | The generic stoichiometry-driven kernels (`network_builder`), the single network path: compiled RHS/Jacobian equal the `reactions` reference to machine precision; the formal N/Z conservation check (passes for real nets, fires on a broken one); numerical baryon-number conservation; the full `UpdateNuclearRates` driver methods (rhs/rhsMT/rhsLT + Jacobians); era-independent table invariants (buffer-order lengths, per-reaction A/Z conservation); and the `amax` mass-cutoff filter (correct count, nuclide bound, conservation, and invalid-value rejection). |
| `test_large_network.py` | The large network: it loads (59 nuclides, 429 LT reactions) and passes the formal conservation check; the vectorised rate buffer stays finite/bounded across the LT range; and a full solve conserves baryon number (to 1e-10) while matching `large, amax=8` on the light elements *and on the free neutron* `n` (the reverse-rate-clamp guard, see "Per-nuclide final abundances" below). |
| `test_nuclear_qed.py` | QED corrections to radiative-capture rates (Pitrou & Pospelov 2020): correction factors are > 1 and sub-percent; the npTOdg polynomial matches its T9→0 cap; the four Kroll-formula reactions increase monotonically with T9; reference magnitudes at T9=0.1 GK are pinned to ±2e-6; non-QED reactions are unchanged; p_* variations stack correctly on the corrected median; and a full solve with the flag on shifts D/H by a detectable but sub-percent amount. |
| `test_regression.py` | Final abundances: loose default-precision sanity checks, tight `reference`-marked checks against the "Validation reference" values below, the per-nuclide table check (every cell of it, at 1e-4 relative), the physical-effect guards for the n↔p correction flags (Born mode and CCRTh both lowering YP), no-numba full solve checks (pure-Python kernels must match JIT to 1e-4), and the `amax` cutoff verification (the large network filtered to A ≤ 20 matches the full large network's light elements to ~1e-3). |
| `test_wheel_smoke.py` | The `wheel`-marked "pip install" smoke test: builds a wheel, installs it into a clean venv, and runs a small-network solve there to catch package-data/path regressions (e.g. `rates/` not shipped, or a path computed relative to the source tree instead of the installed package) that an editable install would not reveal. |
| `test_docs_consistency.py` | Every documentation claim that no other test would catch, all as static file reads (no solve, so it stays in the fast lane): the `CPRIMAT_VERSION` ↔ `pyproject.toml` version sync; the Streamlit-Cloud chain (`requirements.txt` → an existing `wheels/*.whl` at the right version); that both generated parameter templates are byte-for-byte what `primat.tools.gen_param_templates` produces *and* list exactly the `DEFAULT_PARAMS` keys; that `runfiles/primat_reference_run.py` still declares the four reference-run parameters verbatim and `PRIMATConfig` recognises every key it sets; `PRIMATConfig`'s `save_nTOp`/`save_nTOp_thermal` defaults; a batch of README guards (the `--set KEY=VALUE` form, the MC key names, the Python-only-features list vs `backend.py`'s actual gate, the unified evolution-schema columns, the result-dict key table, the per-reaction rate columns); that `notebooks/README.md`, `docs/tutorials/index.md` and this file each list every file they index; and that this README's "Validation reference" tables — observables *and* the per-nuclide table — quote exactly the constants in `tests/reference_values.py`. |
| `test_gui.py` | The optional Streamlit GUI (`primat.gui`): `import primat` does not pull in `primat.gui`/streamlit; the parameter-form metadata covers `amax` (the one `None`-default key) and the network choices; an end-to-end `AppTest` run of `primat/gui/app.py` reproduces `test_cli.py`'s pinned default-run values (Neff/YPBBN/D-H and the per-nuclide table) -- i.e. the GUI drives `PRIMAT` identically to the CLI; the abundance-evolution panel renders with its default "light elements" nuclide selection; the `amax` widget appears only for `network='large'`; and an invalid flag combination (`spectral_distortions=True` with `incomplete_decoupling=False`) is shown as a clean `st.error` rather than a traceback. Skipped entirely if the `gui` extra is not installed. |
| `test_backend_parity.py` | Backend parity: `primat._primat_c` (C) vs `primat.main.PRIMAT` (Python) — the result-dict shape, the written evolution-TSV header, and numerical agreement across the two backends. |
| `test_cache_parity.py` | Cross-backend **cache** parity, the companion to `test_backend_parity.py`: the two backends *share* every on-disk cache, and this is what makes sharing safe rather than merely convenient. Drives both backends into separate `cache_dir`s on a coarse grid and asserts (a) **hash identity** — both emit the same `nTOp_<hash>.txt` / `electron_thermo_<hash>.txt` filename and the same QED fingerprint header, which pins `cpr_weak_rate_fingerprint`/`cpr_constants_hash`/`cpr_qed_fingerprint` against their Python counterparts field-for-field — and (b) **column agreement** at documented, measured tolerances (n↔p rates 1e-8; electron-thermo 1e-9; QED e2/e3 1e-15 on a column-scale metric, since pointwise-relative is dominated by the Boltzmann-suppressed tail). The CCRTh thermal table's *contents* are deliberately excluded (independent Monte-Carlo streams, minutes to recompute); its fingerprint is still pinned via a both-backends-hit-the-shipped-file assertion. Skipped without the C extension. ~13 s. |
| `test_background.py` | Direct unit tests for `primat.background.StandardBackground`. |
| `test_custom_background.py` | The `custom_background` mode: a user-supplied (T, t, a) table drives the cosmological background while the nuclear network is solved with instantaneous-decoupling n↔p weak rates. |
| `test_decay_rates.py` | Radioactive-decay treatment in the `large` network. |
| `test_detailed_balance.py` | `compute_detailed_balance_coefficients` reproduces every (α, β, γ) row of `detailed_balance.csv` from the nuclide data alone — i.e. the shipped reverse-rate coefficients are a derivable consequence of spins/masses/Q-values, not hand-maintained magic numbers. |
| `test_deuterium_network.py` | The `network="large", amax=2` configuration reproduces the old standalone `deuterium` network (single-reaction `n_p__d_g`), the clean-slate starting point for custom networks. |
| `test_evolution.py` | The unified time-evolution schema (`primat.evolution`) on its own, with no BBN solve: round-trips a synthetic `EvolutionResult` through `dump_evolution`/`load_evolution` (to a string and to a file) and pins the exact header. This is the Python-writer half of the schema's guard; the C writer's header is compared against it in `test_backend_parity.py`. |
| `test_gui_custom_network.py` | The GUI's "Manage networks" dialog and the "Create custom network" dialog it gates (`primat.gui.params_form`). |
| `test_gui_run_view.py` | Direct parity test for `primat.gui.run_view.GuiRun`. |
| `test_neutrino_history.py` | The pluggable neutrino-sector background. |
| `test_notebooks.py` | Notebook smoke test: papermill-executes the example notebooks. |
| `test_nuclear.py` | The auto-derivation fallback in `reaction_stoichiometry` and the duplicate-entry check in `load_network`. |
| `test_rate_variations.py` | Nuclear rate variation and MC uncertainty propagation. |
| `test_runfiles.py` | The example scripts in `runfiles/`: each runs as a real subprocess from a throwaway cwd and must exit 0 with no traceback; and `primat_run.py`'s *printed* YP/(D/H) are checked against the "Validation reference" below at the routine tolerance, so the documented "run this after any modification" workflow is automated rather than honour-system. |
| `test_sensitivity.py` | `primat.sensitivity` — the logarithmic-sensitivity API. |
| `test_spectral_distortions.py` | Non-thermal neutrino spectra, each pinned as a *difference* between two full solves (`solve` tier): `spectral_distortions` on/off (small but non-zero on D/H, zero distortion energy in NEVO by construction); the analytic y-type (`y_SZ`) and gray (`y_gray`) distortions shifting Neff, and `finite_mass_corrections` genuinely gating the SD-FM term; and neutrino chemical potentials (`munuOverTnu`, per-flavour `xi_*`) — Neff even in the sign, per-flavour knobs reducing to the common one, and the discriminating case that `xi_mu` alone gravitates but must *not* shift the n↔p rates. |
| `reference_values.py` | (helper, not a test) Centralised default-run reference observables shared by test_cli/test_gui/test_regression, and the validation-reference constants (single source, see test_docs_consistency). |
| `_oracles.py` | (helper, not a test) Test-only reference RHS/Jacobian oracle implementations the nuclear-network tests compare against. |

## Validation reference (authoritative copy)

(This section moved here from the untracked CLAUDE.md so that CI and public
clones carry it; tests parse THESE tables — see
`tests/test_docs_consistency.py`, which fails if they drift from
`tests/reference_values.py`.)

The values below hold at the defaults `Omegabh2=0.02242`,
`spectral_distortions=True` and `nuclear_qed_corrections=True` (the last from
the radiative-capture QED corrections of Pitrou & Pospelov 2020). They were
produced by `runfiles/primat_reference_run.py` — a **high-precision** run with
`numerical_precision=1e-10`, `sampling_temperature_per_decade=2000`,
`sampling_nTOp_per_decade=125`, `T_start_cosmo_MeV=100` and
`rate_grid_npts=4000` explicit, so this reference stays decoupled from the
routine-run defaults.

### Which tolerance applies to which command

The two columns below are **not** interchangeable — mixing them up makes the
documented check fail on a perfectly healthy tree:

| You ran | Precision | Use column | Automated by |
|---------|-----------|-----------|--------------|
| `runfiles/primat_reference_run.py` | `numerical_precision=1e-10` (+ the settings above) | **Reference** | `pytest tests/ -m reference` |
| `runfiles/primat_run.py` | `numerical_precision=1e-7` (the default) | **Routine** | `tests/test_runfiles.py::test_primat_run_matches_the_validation_reference` |

The routine column is looser because a default-precision solve carries ~1e-8
of adaptive-step jitter in D/H, and because the two backends differ by a few
parts in 1e6 (see `tests/test_backend_parity.py`). Measured 2026-08-05,
`primat_run.py` lands 3.8e-9 (C backend) / 5.7e-10 (Python backend) below the
reference D/H — i.e. *outside* the ±3e-9 reference bound, and correctly inside
the ±2e-8 routine one.

Note also that `runfiles/primat_run.py` runs **only** `large, amax=8`; the
small-network table is checked by the `reference` tier and by
`tests/test_regression.py`, not by that script.

**Small network** (`network="small"`):

| Observable | Expected | Reference tol. | Routine tol. |
|------------|----------|----------------|--------------|
| YP (BBN) | 0.24699819 | ±1e-5 | ±1e-5 |
| D/H | 2.4358800e-05 | ±3e-9 | ±2e-8 |

**`large, amax=8`** (the 68-reaction subset of `large`):

| Observable | Expected | Reference tol. | Routine tol. |
|------------|----------|----------------|--------------|
| YP (BBN) | 0.24700154 | ±1e-5 | ±1e-5 |
| D/H | 2.4365900e-05 | ±3e-9 | ±2e-8 |

A result outside the applicable bound indicates a regression.

### Per-nuclide final abundances (small / large+amax=8 / large)

Final mass-fraction abundances `Y` of the small-network nuclides at the end of
BBN, from the **default** run (`numerical_precision=1e-7`) on the auto backend
(C when built). Snapshotted 2026-08-05.

**Read these to 5 significant figures, not 7.** The two backends agree on them
to ≤2.2e-5 relative, and an ordinary numerics improvement moves them at the
1e-5 level (review passes 6–7 moved `n` by 6.2e-5). The pinned bound is
therefore **1e-4 relative** — see `NUCLIDE_REL_TOL` in
`tests/reference_values.py`, which holds these same numbers as the single
source and is checked live by
`tests/test_regression.py::test_per_nuclide_abundances_match_the_reference_table`.

| Nuclide | small | large, amax=8 | large |
|---------|-------|----------------|-------|
| n   | 3.997234e-16 | 3.996332e-16 | 3.996365e-16 |
| p   | 7.529405e-01 | 7.529372e-01 | 7.529372e-01 |
| H2  | 1.834071e-05 | 1.834570e-05 | 1.834580e-05 |
| H3  | 5.851937e-08 | 5.838985e-08 | 5.839019e-08 |
| He4 | 6.174982e-02 | 6.175066e-02 | 6.175066e-02 |
| Li7 | 2.181375e-11 | 9.178225e-11 | 9.178156e-11 |
| Be7 | 3.966446e-10 | 3.223685e-10 | 3.223652e-10 |

The full large network must match `large, amax=8` on every light element to
≲1e-3, **and also on the free neutron `n`**: the reverse-rate clamp in
`primat/network_data.py` (see its module docstring, "exothermic blow-up")
removed the spurious low-T flooding of heavy nuclides such as B10, which
previously fed β-delayed neutron emission and inflated `n` to ~7e-13. With the
clamp the heavy tail is negligible and `n` tracks the `amax=8` value (measured
8.3e-6 relative). `n` is included in
`test_large_network.py::test_large_solve_conserves_baryon_and_matches_amax8`'s
comparison loop precisely so that a regression of that clamp — a factor ~1750
on `n` — fails a test rather than passing unnoticed.

Baryon number is conserved by every network to `sum_s A_s Y_s - 1 ≈ 1.6e-12`
(measured 2026-08-05); the same test pins it at 1e-10.

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

## Known cross-backend divergences

README.md's "Backend parity contract" says the two backends mirror each
other's physics and numerics. Three places knowingly do not, and this is the
tracked list of them, so that a future reviewer neither re-discovers them nor
"fixes" one that was measured and kept on purpose. The *magnitudes* are pinned
and kept current in `tests/test_backend_parity.py`'s module docstring — this
section states the causes and the decisions.

| Divergence | Status |
|---|---|
| **HT-era integrator.** Python integrates the n↔p-only HT era with `LSODA`, C with Dormand-Prince RK45. | **Intentional.** This one mismatch is the whole YP gap: patching Python's HT method to RK45 reproduces C's `YPBBN` exactly. It therefore does not shrink with tighter `numerical_precision`. Aligning both on BDF was tried and *degraded* YP parity, so alignment would buy reviewability, not accuracy. Documented in place in both backends' `nuclear_network` sources. |
| **`external_scale_factor` interpolant.** Python reads T(a) linearly inside its time-integration RHS; C fits a not-a-knot cubic over the same nodes. | **Intentional.** The C cubic is a performance workaround — its RK45 stepper rejected ~65 % of steps on the kinks — and leaves the solution unchanged; LSODA has no such problem. Making Python match was tried and measured *worse* in both self-convergence and cross-backend YP, because in this mode a(T) is itself a table read (NEVO's `x` column), so a cubic through those nodes manufactures curvature the data does not contain. |
| **Residual D/H gap.** At converged tolerance the two settle on different D/H, with Li7/H and He3/H following and YP/Neff agreeing far better. | **Open.** It does not shrink with `numerical_precision` (swept 1e-6 … 1e-10), so it is structural, not round-off. The background is ruled out: `t(T)` agrees cross-backend and each backend self-converges, and the weak-rate tables are tabulated on identical grids (see `test_cache_parity.py`). The cause is downstream — the nuclear network itself, or how the n↔p rates couple into it. It sits far below observational significance (observed D/H is known to ~1 %). |

Everything else that round-1 review found divergent between the two backends
was closed rather than documented; `git log` is the record of those.
