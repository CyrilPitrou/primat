# Performance & benchmarks

Typical wall-clock times for a single BBN solve and for an MC
uncertainty-propagation run, both backends, plus guidance on the knobs that
trade accuracy for speed.

## Reference machine

Apple Silicon Mac (macOS, `arm64`), Python 3.12, `primat` built with its
default compiler flags (`-O2`/`-O3`). Absolute numbers will differ on other
hardware; the *relative* C-vs-Python speedup and the shape of the
network/backend cost curve should transfer.

## Typical timings (warm cache)

All timings below are *warm-cache*: the shipped
`primat/data/cache_plasma_weak/weak/nTOp_*.txt` (n<->p weak rates) and
`primat/data/cache_plasma_weak/plasma/QED_*.txt` (QED pressure corrections) tables already
match the run's configuration fingerprint, so initialisation skips
recomputing them. This is the timing every run after the first sees; see
[Cold-cache cost](#cold-cache-cost) below for the one-time cost when they
don't match.

| Run | Wall time |
|-----|-----------|
| small (c) | 0.040 s |
| small (python) | 0.460 s |
| large, amax=8 (c) | 0.134 s |
| large, amax=8 (python) | 0.684 s |
| MC-100 (small, c) | 0.906 s |

- **small** — the 12-reaction default network.
- **large, amax=8** — the large network's reactions filtered to A <= 8 (68
  reactions, the old "medium" network's exact equivalent); the LT-era
  solve dominates the extra cost over `small`.
- **MC-100** — `run_mc(100, params={"network": "small"}, force_backend="c")`,
  100 samples, default `n_jobs=-1` (parallel across CPU cores); reported as
  a single wall-clock measurement rather than a per-solve minimum, since MC
  parallelism/scheduling overhead is itself part of what the number
  documents.
- The C backend's advantage grows with network size (~12x for `small`,
  ~5x for `large, amax=8` here) because the Python solver's per-step
  overhead (mostly `scipy.integrate.solve_ivp` callback/array-marshalling
  cost) is roughly independent of the RHS's cost, so it amortises worse as
  the RHS (nuclide count) grows; `numba`-JIT-compiled kernels reduce but
  don't eliminate this gap on the Python side.
- The full `large` network (~429 reactions, ~59 tracked nuclides) is not
  included above — it is a further ~5-10x slower than `large, amax=8` on
  either backend (LT-era stiffness scales with reaction count) and is
  primarily useful for including the heavy-nuclide tail, not for
  fast iteration; use `amax=8` for MC/scan workloads that only need the
  light elements at ≲1e-4 accuracy (see `tests/README.md`'s per-nuclide
  abundance table).

## Cold-cache cost

The n<->p weak-rate cache (`primat/data/cache_plasma_weak/weak/nTOp_<hash>.txt` and
`nTOp_thermal_<hash>.txt`) is the most expensive part of initialisation when
it must be (re)computed rather than loaded — see
{doc}`howto/weak-rate-cache` for the fingerprinting mechanism that makes
this a one-time cost per distinct configuration:

- Non-thermal rate (Born + finite-mass + Coulomb/radiative + spectral-distortion
  corrections): ~1.8 s, numerical integration, no external solver.
- Finite-temperature radiative correction (`thermal_corrections=True`,
  the default): an additional multi-minute `vegas` Monte-Carlo integration,
  cached separately in `nTOp_thermal_<hash>.txt`.

Both are skipped (loaded from cache instead) whenever the current
configuration's fingerprint matches an existing cache file — the common case
for repeated runs at default settings, or across an MC batch (every sample
shares the same background/weak-rate configuration, so only the *first*
sample in a fresh process pays this cost). Set `thermal_corrections=False`
during exploratory work if the multi-minute `vegas` cost is a bottleneck and
the ~1e-3 Neff-level correction it provides is not needed for that
particular check.

## Tuning knobs

| Knob | Default | Effect on speed |
|------|---------|------------------|
| `numerical_precision` | `1e-7` | `solve_ivp` `rtol` for every era; raising it (looser tolerance, e.g. `1e-5`) speeds up the LT-era solve roughly linearly in step count but degrades the `D/H`/`YPBBN` precision pinned in `tests/reference_values.py` — do not raise it past the point where results are no longer distinguishable from the reference-run tolerances. |
| `rate_grid_npts` | 1000 | Master T9 grid size every rate table is resampled onto at load time; this is a fixed one-time cost per solve (not per RHS call), so it matters more for `small`/fast networks than for `large`. It also sets an accuracy floor — see "What the default grids cost" below. |
| `sampling_nTOp_per_decade` / `sampling_nTOp_thermal_per_decade` | see `primat/config.py` | Only affects *cold-cache* weak-rate computation time (see above); irrelevant once cached. |
| `amax` | `None` | Filters any named network to reactions with `A <= amax` — the single biggest lever for LT-era solve time on the `large` network (see the `large` vs `large, amax=8` gap above). |
| `n_jobs` (`run_mc`) | `-1` (all cores) | MC samples are embarrassingly parallel (independent solves); `n_jobs=1` is useful for reproducible profiling of a single sample's cost, not for a production MC run. |
| `force_backend` | `None`/`"auto"` (prefers C) | The only `PRIMAT.__init__` feature the C backend cannot express is `background=` (a custom `Background` object) — seeing an unexpected Python-backend fallback in `log_backend=True` output usually means that is set. (`extra_rho` and `decay_era` *are* supported on the C backend.) |

The high-precision reference-run settings (`numerical_precision=1e-10`,
`sampling_temperature_per_decade=2000`, `sampling_nTOp_per_decade=125`,
`rate_grid_npts=4000` — all set by `runfiles/primat_reference_run.py`) trade
several minutes of wall time for the extra digits behind the tolerance bands
in `tests/README.md`'s "Validation reference" — not something to use for
routine runs or MC sampling. Most of that time is the thermal (CCRTh) `vegas`
recompute: those settings miss both shipped weak-rate caches.

## What the default grids cost

`numerical_precision` is not the only thing limiting a run's accuracy, and past
about `1e-8` it stops being the binding one. Two *fixed sampling grids* take
over: the master T9 grid every rate table is resampled onto (`rate_grid_npts`)
and the background's photon-temperature grid
(`sampling_temperature_per_decade`). Both converge at second order, so their
limits are well defined; the default's distance from them, measured on the
`small` network at `numerical_precision=1e-10` and identical on both backends:

| Grid | Default | Limit measured at | `YPBBN` | `D/H` | `He3/H` | `Li7/H` |
|------|---------|-------------------|---------|-------|---------|---------|
| `rate_grid_npts` | 1000 | 16000 | 1.4e-08 | −6.3e-06 | −1.1e-05 | **+9.7e-05** |
| `sampling_temperature_per_decade` | 600 | 4800 | +2.9e-07 | −7.2e-06 | −2.7e-06 | +9.7e-06 |

For comparison, one more decade of `numerical_precision` (`1e-10` → `1e-11`)
moves `D/H` by 8.0e-09 and `Li7/H` by 4.6e-09 — three orders of magnitude less.

What this means in practice:

* `YPBBN` and `D/H` are unaffected at the level anything is pinned to: the two
  terms together are 1.4e-05 relative in `D/H`, i.e. 3.3e-10 absolute, about
  a tenth of the ±3e-9 regression tolerance in `tests/README.md`.
* `Li7/H` is the exception. Its default value carries ~1e-04 of grid error, so
  the last two of the six decimals it is quoted to are grid artefacts rather
  than physics. Comparisons *between* runs at the same grid settings are
  unaffected — the error is systematic, not noise.
* Both backends carry the same error, so no cross-backend comparison can
  reveal it. `tests/test_backend_parity.py` is not the check for this.
* `runfiles/primat_reference_run.py` already sits close to both limits
  (`rate_grid_npts=4000`, `sampling_temperature_per_decade=2000`), which is
  part of why its numbers are the ones `tests/reference_values.py` pins.

Raising the defaults to the limit column would cost roughly 16× the rate
resampling and 8× the background setup for a shift that stays inside the
existing tolerance bands, which is why the defaults are where they are.

## Where the Python backend's time goes

The pure-Python backend is ~12× slower than C on a warm `small` run, and that
gap is *not* the physics: the compiled numba RHS kernel is a few per cent of
the runtime. The rest is Python-level glue paid once per ODE step, plus the
one-off background construction.

Share of a warm `small` run under `cProfile` (the profiler roughly doubles the
absolute time, so read the column as proportions):

| Where | Share | Nature |
|---|---|---|
| scipy's BDF stepping logic | ~35 % | pure-Python predictor/corrector, error norms and order selection, 4.4 k steps |
| RHS + Jacobian evaluation | ~35 % | 15.5 k calls: background lookups, `fill_buffer`, then the numba kernel (~6 % on its own) |
| background construction | ~22 % | the `a(T)` and `t(T)` LSODA solves, once per run |
| weak-rate setup | ~7 % | reading and fitting the cached n↔p tables, once per run |
| dense LU (LAPACK `getrf`/`getrs`) | ~2 % | the linear algebra itself |

Four things keep the per-step glue down. All four are **bit-exact** — they
change how a value is computed, never which value:

1. njit scalar evaluators (`primat/weak_rates/fast_eval.py`) instead of scipy
   interpolator calls for the n↔p rates. The weak-rate interpolant is a
   cross-backend parity contract, so its coefficients are pulled out of the
   fitted spline rather than re-fitted.
2. `fill_buffer` fused into a compiled kernel.
3. A one-entry memo keyed on bit-identical `t`, so the RHS and Jacobian
   closures do not redo the same `T_of_t(t)` / `rhoB_BBN(t)` / weak-rate work
   at the same step.
4. Scalar fast paths for the background's own linear interpolants
   (`background._scalar_linear_eval`, used by `T_of_t`/`t_of_T` and
   `rhoB_BBN`), and a BDF subclass calling LAPACK's `getrf`/`getrs` directly
   instead of through `scipy.linalg.lu_factor`/`lu_solve`
   (`nuclear_network._bdf_method`). Both check at build time that they
   reproduce what they replace, and fall back to scipy if they cannot;
   `tests/test_refactor_invariants.py` pins both.

**What remains is scipy's pure-Python BDF stepper** — the top row above, and
the only route left to matching C. Replacing it (a numba port, or
numbalsoda/SUNDIALS) changes the integrator itself and therefore every
regression pin, so it is a deliberate decision rather than an optimisation.
`primat-c/src/ode_bdf.c` is a term-for-term transcription of scipy's
`_ivp/bdf.py`, which gives such a port a checked reference implementation to
follow.

To re-profile:

```bash
python -c "
import cProfile, pstats, io
from primat.backend import run_bbn
run_bbn(force_backend='python')          # warm caches/JIT
pr=cProfile.Profile(); pr.enable(); run_bbn(force_backend='python'); pr.disable()
s=io.StringIO(); pstats.Stats(pr,stream=s).sort_stats('tottime').print_stats(20)
print(s.getvalue())
"
```

Note that `numba` is recommended but not mandatory: every njit kernel has a
pure-numpy fallback, so the package keeps working (more slowly) without it.

## Regenerating this table

```bash
python runfiles/benchmark.py
```

Runs each backend/network combination `--repeats` times (default 3) and
reports the minimum wall time (steady-state, avoiding first-call
JIT/import warm-up noise), plus a 100-sample MC run on the C backend; prints
a ready-to-paste Markdown table. Pass `--quick` for a fast (few-second)
smoke test whose numbers are not meaningful as a benchmark — only used by
`tests/test_runfiles.py` to check the script still runs. Re-run and update
the table above whenever the physics/numerics change enough to shift these
numbers by more than noise, or when benchmarking on a new reference machine
(update the "Reference machine" description too).
