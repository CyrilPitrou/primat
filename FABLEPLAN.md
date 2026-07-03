# FABLEPLAN.md — Improvement plan for primat / primat-c

> **Status (2026-07-03): all actionable items complete.** 1.1, 1.3, 1.5, 2.1,
> 3.1, 3.2, 3.3 are done (commit refs on each heading); 1.2 was dropped by the
> author and 1.4 was already moot. Nothing outstanding.

Scope agreed with the author (2026-07-03): the plan retains the **structural
cleanups** and the **investigation of the un-root-caused C-vs-Python D/H
difference**, plus a short list of cheap hygiene/safety items. Explicitly
dropped after review: relocating runtime caches out of the package tree
(caches stay in `primat/data/`), the declarative parameter registry (the
three-way manual sync of the 73 parameters stays, as documented in
CLAUDE.md), and compiled-network/numba caching for the Python MC path
(lengthy computations use the C backend, so Python-side JIT warm-up is not
worth optimizing).

Each task carries an effort estimate (S = hours, M = a day or two) and a
**model assignment**: *Sonnet* for well-scoped mechanical refactors with a
clear correctness check (the test suite + the CLAUDE.md reference values),
*Opus* for tasks needing sustained numerical/physics judgment or open-ended
diagnosis.

Ground rules for every task below (from CLAUDE.md):

- After any change, `python runfiles/primat_run.py` must reproduce the
  reference observables within tolerance (YP 0.24700028 ±1e-5,
  D/H 2.43500e-5 ±3e-9 for the small network).
- Pure refactors (renames, splits, no numerical effect) do **not** need to
  be mirrored in `primat-c/`; anything touching numbers does.
- Preserve the heavy-commenting convention: every extracted helper keeps or
  gains a physics-explaining docstring with paper citations.

---

## 1. Structural cleanup

### 1.1 Split the longest functions — **Sonnet** (M total, independent sub-tasks) — **DONE (2026-07-03)**

Measured by AST scan, the outliers:

| Lines | Function | File | Sub-task model | Status |
|-------|----------|------|----------------|--------|
| 359 | `_build_analytic_distortion` | `primat/neutrino_history.py` | Sonnet | done, commit `7602d5a` |
| 357 | `_L_CCRTh_interpolants` | `primat/weak_rates/corrections.py` | Sonnet | done, commit `e568b77` |
| 342 | `NuclearNetwork.solve` | `primat/nuclear_network.py` | Sonnet | done, commit `2ee88b7` |
| 273 | `PRIMATConfig.__init__` | `primat/config.py` | Sonnet | done, commit `dcfd95b` |
| 239 | `mc_uncertainty` | `primat/main.py` | Sonnet | done, commit `7daf836` |
| 232 | `StandardBackground._setup_background_and_cosmo` | `primat/background.py` | Sonnet | done, commit `e706096` |
| 213 | `render_sidebar_form` | `primat/gui/params_form.py` | Sonnet | done, commit `9b2636f` |

Each split verified against `runfiles/primat_run.py`'s reference values and
the relevant test file(s) before committing; the physics-dense pair
(`_build_analytic_distortion`, `_L_CCRTh_interpolants`) additionally got a
targeted numerical-equivalence check (git-stash before/after comparison)
since their normal test coverage doesn't exercise every code path.

These are all *pure* refactors — extract helpers, move code, change no
numbers — with the full test suite as the safety net, hence Sonnet. Do them
one function per commit so a regression bisects trivially. Priorities:

- **`NuclearNetwork.solve`** → `_solve_HT()` / `_solve_MT()` / `_solve_LT()`
  plus a short orchestrator. The three eras are already conceptually
  separate in the module docstring; the code should mirror that. Each era
  method returns `(t, Y)`; the orchestrator stitches them. Bonus: single
  eras become unit-testable in isolation.
- **`PRIMATConfig.__init__`** → per-group validators
  (`_validate_nevo_files()`, `_validate_data_dirs()`,
  `_validate_physics_flag_combos()`, ...). Behaviour identical (same
  errors, same messages, same order where tests depend on it); the
  constructor becomes a readable checklist.
- **`mc_uncertainty`** → separate the prev-reuse guard, the sample-loop
  dispatch, and result assembly (see also 1.3).
- **`_L_CCRTh_interpolants` / `_build_analytic_distortion`** — physics-dense;
  split only at natural physics boundaries (per correction term / per
  distortion type), each helper carrying its own equation citation. Still
  Sonnet (no numbers change), but the reviewer should check the citations
  landed on the right fragments.
- **`render_sidebar_form`** — split per sidebar group; the GUI tests
  (`tests/test_gui*.py`) pin the behaviour.

### 1.2 ~~Split `network_data.py` (2127 lines) into a subpackage~~ — **dropped**

Author decision (2026-07-03, after 1.1 completed): not worth doing. The
subpackage split below is left for reference only.

### 1.2 (reference only, not being done) Split `network_data.py` (2127 lines) into a subpackage — **Sonnet** (M)

`primat/network_data.py` mixes six concerns. Suggested layout:

```
primat/network/
  __init__.py      # re-export everything network_data exports today (no API break)
  loader.py        # load_network, _parse_network_entries, amax filter, era restriction
  tables.py        # rate-table loading, _resample_rate_table, _LinearRate, overlay resolution
  balance.py       # compute_detailed_balance_coefficients, reverse-rate cap/clamp
  qed.py           # _apply_nuclear_qed / _qed_nuclear_rescale
  display.py       # nuclide_latex, reaction_display_name, reaction_category, grouping
  custom.py        # custom-network injection (kept_to_custom_network glue)
  definition.py    # NetworkDefinition dataclass + UpdateNuclearRates
```

Keep `primat/network_data.py` as a thin re-export shim for one release so
external imports (`from primat.network_data import load_network`) keep
working. Pure code motion → Sonnet; the module's `__all__` (42 names) is
the checklist that nothing was dropped.

### 1.3 Deduplicate the MC `prev`-reuse guard — **Sonnet** (S) — **DONE** (commit `3e9d217`)

The seed/quantities/params/custom_network compatibility check exists both
in `primat.main.mc_uncertainty` and in `primat.backend.run_mc` (which adds
the backend-match condition). Extract one `mc_prev_is_reusable(prev, ...)`
helper used by both so the two guards can never drift. `tests/test_mc.py`
pins the reuse semantics.

### 1.4 ~~Break the `weak_rates` import cycle~~ — already resolved, moot

Checked 2026-07-03: no cycle exists. `corrections.py` and `api.py` both
import `_setup_fd_impls` and the FD kernels straight from the leaf module
`integrands.py`, never from `weak_rates/__init__.py`; a static AST scan of
`__init__.py`/`api.py`/`corrections.py`/`integrands.py`/`cache.py` shows no
back-edges, and importing each submodule independently (and via
`importlib`) works cleanly. This was fixed by commit `7abbfa7` ("weak_rates:
kill scalar/_v Fermi-Dirac kernel duplication, split into subpackage"),
before this plan was written — the description above was stale.

### 1.5 Repo hygiene — **Sonnet** (S) — **DONE** (commit `27c7f4d`)

- Delete the committed editor backup `.github/workflows/build_linux.yml~`.
- The runtime weak-rate caches stay in `primat/data/weak/` (author's
  decision), so add the runtime-generated `nTOp_*.txt` pattern to
  `.gitignore` while keeping the shipped cache files tracked (ignore by
  pattern, `git add -f` the shipped set once), so `git status` stays clean
  after routine runs.
- Decide tracked-vs-ignored for `.codex/` and `.graphifyignore` and commit
  that decision.

---

## 2. C-vs-Python D/H gap investigation

### 2.1 Root-cause the ~1.7e-8 small-network D/H discrepancy — **Opus** (M) — **DONE (2026-07-03)**

**Root cause: a weak-rate interpolation-*scheme* mismatch.** Python interpolated
the cached n↔p rate table with a linear-space quadratic spline
(scipy `interp1d(kind='quadratic')`); the C backend used a local 3-point
Lagrange quadratic (`cpr_interp_quadratic_local`). On the *same* cached nodes
the two curves differed by up to ~1e-4 relative through the n/p freeze-out
window (T ~ 0.2–2 MeV), C systematically above Python — propagating to a
systematic ~2.5e-5 YP / ~1.8e-5 D/H offset (Neff always identical). Ruled out
BDF controller noise: the gap *plateaued* (didn't scale) as
`numerical_precision` tightened 1e-7 → 1e-10. Measured (subsampling the shipped
80-pt/decade table) that log10-log10 cubic is the most accurate scheme at every
density and that more points/decade is the wrong lever — the fix is a better
*scheme*, not a denser grid.

**Fix (mirrored, both backends):** both now interpolate the weak rates with the
same log10-log10 not-a-knot cubic spline — Python
`weak_rates.api._weak_rate_loglog_interp`, C `cpr_weak_rate_nTOp`/`cpr_weak_rate_pTOn`
via a new `CPRWeakInterp` (`cpr_cubic_spline_fit_notaknot`) — the scheme the
nuclear rate tables already share, and ~1–2 orders more accurate than either old
scheme at the shipped density. The backward p→n rate's exp(−Q/T)-suppressed
zero prefix is handled by splining only the positive suffix and returning 0
below it (forward rate is positive throughout and extrapolates). This collapsed
the YP gap to ~1e-6 and the D/H gap to ~1e-6 (`small`) / ~6e-6 (`large`+amax=8);
the residual is the separate, smaller nuclear-rate-interpolation + BDF
step-sequence difference (doesn't touch YP). `test_backend_parity.py`'s
cross-backend D/H tolerance tightened 1e-3 → 5e-5; CLAUDE.md and the test
docstring updated. Both reference configs stay within the CLAUDE.md tolerances.

*Historical framing (superseded by the finding above):* `tests/test_backend_parity.py`
budgeted a ~1.7e-8 absolute (~7e-4 relative) C-vs-Python D/H gap for
`network="small"` as an unexplained cross-backend tolerance — larger in spirit
than the ±3e-9 same-backend regression tolerance CLAUDE.md enforces. Genuine
numerical detective work (two solver stacks, several plausible culprits, no
single failing assert), hence **Opus**.

Suggested bisection strategy:

1. **Localize in time.** Run both backends with
   `output_time_evolution=True` and diff `Y_H2` (and `Y_n`, `Y_H3`,
   `Y_He3`) at matched times/temperatures — the unified evolution schema
   makes the two TSVs directly comparable. Determine whether the gap is
   already present at the MT→LT boundary or grows during LT.
2. **Localize in the rates.** If the gap is born in MT/LT, compare the two
   backends' per-reaction forward/backward rate buffers at a few
   temperatures (add a temporary debug dump on the C side). Likely
   suspects, in order: interpolation-kind mismatch (Python's log-log cubic
   `_resample_rate_table` / `_LinearRate` vs the C spline/interp path),
   the reverse-rate cap and `_EXP_CAP`/`_FLOOR` clamps
   (`primat/network_data.py` vs their C mirrors), and the weak-rate
   interpolants' grids.
3. **Distinguish physics from solver noise.** If rates agree, tighten
   `numerical_precision` on both sides (1e-8, 1e-9, 1e-10) and check
   whether the gap shrinks proportionally. If it does, it is BDF
   step-sequence noise between the two controllers; document that
   conclusion (with the scaling evidence) next to the tolerance in
   `test_backend_parity.py`. If it does not, there is a real formula or
   clamp mismatch to fix — and per CLAUDE.md the fix must be mirrored so
   both backends move together.

Deliverable either way: the parity test's cross-backend tolerance is
replaced by (a) a tightened tolerance after a fix, or (b) a comment citing
the measured tolerance-scaling that proves it is controller noise.

---

## 3. Cheap safety items (optional, kept for their cost/benefit)

### 3.1 Automated `CPRIMAT_VERSION` vs `pyproject.toml` sync test — **Sonnet** (S) — **DONE** (commit `020dbd8`)

CLAUDE.md notes there is no automated check for the version-macro sync. A
five-line test in `tests/test_docs_consistency.py` (regex the macro out of
`primat-c/include/config.h`, compare to the pyproject `version`) closes
that gap permanently.

### 3.2 Run the C unit tests under ASan/UBSan in CI — **Sonnet** (S) — **DONE** (commit `eb890b8`)

`make debug` already builds with `-fsanitize=address,undefined`; nothing
runs it automatically. Add a CI job building the unit-test binaries with
the debug flags and running them. Cheapest possible memory-bug detector
for a C codebase with ~250 manual allocation sites and few NULL checks.

### 3.3 Pytest markers for the slow tests — **Sonnet** (S) — **DONE** (commit `2ff287c`)

Mark the wheel-install smoke test, notebook executions, and full
large-network solves `@pytest.mark.slow`; document `pytest -m "not slow"`
as the inner-loop command in `tests/README.md`. Faster feedback during the
1.x refactors above.

---

## 4. Suggested ordering

1. **3.3** pytest markers (makes everything after it cheaper to iterate on)
2. **1.5** repo hygiene, **3.1** version-sync test, **3.2** sanitizer CI — one small commit each
3. **1.3** MC guard dedup, **1.4** import cycle
4. **1.1** long-function splits, one function per commit, `solve()` and
   `PRIMATConfig.__init__` first
5. **1.2** `network_data.py` subpackage split (after 1.1 so the moved code
   is already in its final shape)
6. **2.1** D/H gap investigation (Opus) — independent of the refactors;
   can run in parallel, but land it on a clean tree, not mid-split
