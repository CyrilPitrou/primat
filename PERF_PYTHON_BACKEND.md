# Handover: speed up the pure-Python backend (~3× target, no solver change)

> **STATUS (implemented 2026-07-16).** Tasks 1–3 done. Warm small-network run
> **1.61 s → ~1.08 s** (≈1.5×). Tasks 2 (njit `fill_buffer`) and 3 (per-`t`
> memo) are bit-identical (DoH `2.435900710117168e-05` unchanged). Task 1 (njit
> weak-rate evaluator, `primat/weak_rates/fast_eval.py`) shifts DoH to
> `2.43589158770391e-05` — a **5e-11 absolute / 2e-6 relative** move that is a
> pure BDF-`rtol` tolerance artifact (the LT solver runs at `rtol≈1e-6`, so a
> ~1e-15 change in how the *same* spline is evaluated nudges step selection);
> it sits 1.7e-9 from the pinned reference `2.435909e-5`, inside the `abs=3e-9`
> regression tolerance. **83 tests pass** (regression + parity + weak-rates).
> The remaining runtime is scipy's pure-Python BDF stepper (the out-of-scope
> item below); the background *linear* interps were deliberately **left on
> scipy** because a njit linear reimplementation was not bit-identical (~3e-14)
> and they are already memoised by Task 3. The `~0.5 s` target in the original
> plan assumed also removing the BDF Python overhead, which was ruled out.

## Context and goal

The Python backend (`force_backend="python"`) is ~35× slower than the C
backend on a warm default run (small network): **C ≈ 0.05 s, Python ≈ 1.6–1.9 s**
on the reference machine. Profiling shows the compiled numba RHS kernel is
only ~2 % of the runtime — the rest is Python-level glue per ODE step.

**Goal of this task: remove that glue (tasks 1–3 below), taking the warm
Python run from ~1.6 s to ~0.5–0.7 s, while keeping results numerically
unchanged.** Explicitly **out of scope**: replacing `scipy.integrate.solve_ivp`
/ BDF with numbalsoda/SUNDIALS (that changes the stepper and all regression
pins; a separate decision).

Reproduce the baseline first:

```bash
python -c "
import time
from primat.backend import run_bbn
run_bbn(force_backend='python')          # warm caches/JIT
t0=time.time(); r=run_bbn(force_backend='python')
print('PY warm', time.time()-t0, 's  DoH =', repr(r['DoH']))
"
```

Record the warm time and the full-precision `DoH` before touching anything.

## Profile summary (cProfile, warm run, ≈3 s under profiler, 6.3 M calls)

| Hot spot | Cumulative | Nature |
|---|---|---|
| scipy interpolator `__call__`s (`interp1d`/`PPoly`/`BSpline`) | ~1.4 s | 105 k **scalar** evaluations of `T_of_t`, `rhoB_BBN`, weak-rate `_eval`, thermal-correction lambda; each pays ~10 µs of `asarray`/validation overhead for a trivial spline lookup |
| scipy BDF stepper internals (`bdf.py`) | ~0.8 s | pure-Python stepping logic — the out-of-scope floor |
| `NetworkDefinition.fill_buffer` (`primat/network_data.py:960`) | ~0.75 s | Python+numpy rate-buffer fill, 24 k calls, slice copies + temporaries |
| numba RHS kernel (`primat/network_builder.py:431`) | 0.06 s | the actual physics — already fast |

Call sites that drive all of this: `Y_prime_MT`/`Jacobian_MT`
(`primat/nuclear_network.py:247-253`) and `Y_prime_LT`/`Jacobian_LT`
(`primat/nuclear_network.py:303-309`). Each call does
`rhoB_BBN(t)`, `T_of_t(t)`, then inside `rhsLT`/`JacobianLT` →
`fill_buffer` → `nTOp_frwrd(T)`/`nTOp_bkwrd(T)` (each of which is a scipy
spline chain, see `primat/weak_rates/api.py:102` `_eval` and the thermal
recombination lambdas around `primat/weak_rates/api.py:355` and
`primat/weak_rates/corrections.py:1324`).

## Task 1 — njit scalar evaluators replacing scipy interpolators on the hot path

For each hot interpolant, extract the knot/coefficient arrays **once at
setup** and evaluate with a small `@njit(cache=True)` scalar function
(binary search via `np.searchsorted` + polynomial/linear eval). Targets:

- **Weak rates** `_eval` (`primat/weak_rates/api.py:102`): it is
  `10**cubic_spline(log10 T)` with a clamp and a zero-mask below
  `T_zero_below`. The spline is a not-a-knot cubic `interp1d` in
  log10–log10 space — **this exact interpolant is a cross-backend parity
  contract** (see CLAUDE.md "share the same log10-log10 not-a-knot cubic
  weak-rate interpolant"), so reproduce its values bit-for-bit (or to
  float ulp): pull the `PPoly`-equivalent coefficients out of the fitted
  `interp1d` (`spline._spline` / convert via `scipy.interpolate.splrep`
  breakpoints) rather than re-fitting with a different algorithm.
  Extrapolation below the first node must keep working (the forward rate
  relies on it).
- **Thermal-correction recombination** lambdas
  (`primat/weak_rates/api.py:355-356`, `primat/weak_rates/corrections.py:1324`):
  fold them into the same njit evaluator so `nTOp_frwrd(T_K)` becomes ONE
  compiled scalar call.
- **Background** `T_of_t` and `rhoB_BBN` (attributes of `Background`,
  `primat/background.py`; the frwrd/bkwrd wrappers are at
  `primat/background.py:1096` and `:1107`): same treatment — they are fixed
  1-D interpolants built once per run.

Keep the existing scipy-based construction as the *fitting* step; only the
*evaluation* path changes. Preserve vector-input behaviour where callers use
it (grep for call sites; the hot path is scalar, but e.g. table dumps may
pass arrays — an `np.ndim` dispatch wrapper is fine, only the scalar branch
must be njit-fast).

## Task 2 — fuse `fill_buffer` into a compiled kernel

`fill_buffer` (`primat/network_data.py:960`) is searchsorted + linear blend
of `self._fwd` columns + detailed-balance prefactor
(`alpha * T9**beta * exp(min(gamma/T9, _EXP_CAP)) * fwd`) + floors/clamps,
all on plain contiguous arrays. Move the body into an
`@njit(cache=True)` function taking `(T9, grid, fwd_table, abg, buf, clamp,
nTOp_f, nTOp_b)` — with the n↔p rates passed as already-evaluated floats
(computed via Task 1's njit evaluators). Keep:

- the exact-`T_t` one-entry cache semantics (docstring at
  `network_data.py:961-972` — callers rely on the returned buffer being the
  same object, read-only, consumed immediately);
- the clamp behaviour and `_EXP_CAP`/`_FLOOR` constants exactly;
- the "slice copies are intentional" invariant (`network_data.py:990`): the
  njit kernel must not mutate the cached table columns.

If ambitious, pass the buffer straight into the existing njit RHS/Jacobian
kernels (`network_builder.py`) so one compiled call covers
fill_buffer→RHS; but the simple split already removes most of the cost.

## Task 3 — one-entry memo on `t` for background quantities

scipy's BDF evaluates `fun(t, Y)` and `jac(t, Y)` at the same `t`;
`Y_prime_LT` and `Jacobian_LT` (and the MT pair) currently redo identical
`T_of_t(t)`/`rhoB_BBN(t)`/weak-rate work. Add a one-entry cache keyed on
bit-identical `t` (same pattern as `fill_buffer`'s `_cache_T_t`) computing
`(rho, T_K, nTOp_f, nTOp_b)` once, shared by the RHS and Jacobian closures
in `_solve_MT`/`_solve_LT` (`primat/nuclear_network.py:228` and `:283`).
Exact-equality keying only — no tolerance-based caching (that would change
numerics).

## Hard constraints (from CLAUDE.md — read it first)

- **Results must not change.** These are overhead removals; the target is
  bit-identical output, and at minimum within the same-backend regression
  tolerance (±3e-9 on D/H). If a task can't be done bit-identically
  (e.g. spline coefficient extraction reorders a sum), justify the ulp-level
  diff explicitly in the commit message.
- **No physics/numerics change ⇒ nothing to mirror in `primat-c/`.** If any
  step turns out to alter a formula, clamp, or tolerance, STOP — that's out
  of scope for this task.
- Comment heavily (physicist-readable), per repo conventions. Docstrings
  explain what/why; cite `biblio/Pitrou_etal_PhysReptArxivVersion.pdf` if a
  formula is touched (it shouldn't be).
- numba is *recommended*, not mandatory: the package must keep working
  without it. Follow the existing pattern in `network_builder.py` /
  `weak_rates/integrands.py` for the numba-optional fallback (pure-numpy
  path when numba is absent), and keep `@njit(cache=True)` so the
  compilation cost is paid once (see the cached-integrands commit
  `4e2cb02`).
- Run `graphify update .` after code changes; the repo hooks require
  `graphify query` before exploratory greps/reads.

## Validation checklist (in order)

1. `python runfiles/primat_run.py` from repo root — compare against the
   pinned tables in `tests/README.md` ("Validation reference"). Report
   observables at full precision: Neff to 8 decimals, YP to 8, D/H to 7,
   Li7/H to 6.
2. `DoH` from the baseline snippet above: bit-identical (print `repr`).
3. `pytest tests/test_backend_parity.py tests/test_docs_consistency.py` and
   then the full suite.
4. Re-run the timing snippet warm: expect ~0.5–0.7 s (from ~1.6 s). Re-run
   the cProfile below and confirm scipy interpolator `__call__`s and
   `fill_buffer` have left the top of the profile:

```bash
python -c "
import cProfile, pstats, io
from primat.backend import run_bbn
run_bbn(force_backend='python')
pr=cProfile.Profile(); pr.enable(); run_bbn(force_backend='python'); pr.disable()
s=io.StringIO(); pstats.Stats(pr,stream=s).sort_stats('tottime').print_stats(20)
print(s.getvalue())
"
```

5. Also validate a `network="large", amax=8` run and a
   `spectral_distortions=False` run (different weak-rate correction chain)
   to exercise both thermal-recombination code paths.

## Explicitly out of scope

- Replacing the scipy BDF stepper (numbalsoda / SUNDIALS). That's the
  remaining ~0.8 s floor and the only route to matching C (~0.05 s), but it
  changes the integrator and every regression pin. Do not attempt here.
- Any change to `primat-c/`, the evolution schema, or DEFAULT_PARAMS.
