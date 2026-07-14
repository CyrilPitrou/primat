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
| small (c) | 0.041 s |
| small (python) | 1.468 s |
| large, amax=8 (c) | 0.140 s |
| large, amax=8 (python) | 2.116 s |
| MC-100 (small, c) | 0.885 s |

- **small** — the 12-reaction default network.
- **large, amax=8** — the large network's reactions filtered to A <= 8 (68
  reactions, the old "medium" network's exact equivalent); the LT-era
  solve dominates the extra cost over `small`.
- **MC-100** — `run_mc(100, params={"network": "small"}, force_backend="c")`,
  100 samples, default `n_jobs=-1` (parallel across CPU cores); reported as
  a single wall-clock measurement rather than a per-solve minimum, since MC
  parallelism/scheduling overhead is itself part of what the number
  documents.
- The C backend's advantage grows with network size (~36x for `small`,
  ~15x for `large, amax=8` here) because the Python solver's per-step
  overhead (mostly `scipy.integrate.solve_ivp` callback/array-marshalling
  cost) is roughly independent of the RHS's cost, so it amortises worse as
  the RHS (nuclide count) grows; `numba`-JIT-compiled kernels reduce but
  don't eliminate this gap on the Python side.
- The full `large` network (~429 reactions, ~59 tracked nuclides) is not
  included above — it is a further ~5-10x slower than `large, amax=8` on
  either backend (LT-era stiffness scales with reaction count) and is
  primarily useful for including the heavy-nuclide tail, not for
  fast iteration; use `amax=8` for MC/scan workloads that only need the
  light elements at ≲1e-4 accuracy (see `CLAUDE.md`'s per-nuclide table).

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
| `numerical_precision` | `1e-7` | `solve_ivp` `rtol` for every era; raising it (looser tolerance, e.g. `1e-5`) speeds up the LT-era solve roughly linearly in step count but degrades the `D/H`/`YPBBN` precision documented in `CLAUDE.md` — do not raise it past the point where results are no longer distinguishable from the reference-run tolerances. |
| `rate_grid_npts` | 1000 | Master T9 grid size every rate table is resampled onto at load time; this is a fixed one-time cost per solve (not per RHS call), so it matters more for `small`/fast networks than for `large`. |
| `sampling_nTOp_per_decade` / `sampling_nTOp_thermal_per_decade` | see `primat/config.py` | Only affects *cold-cache* weak-rate computation time (see above); irrelevant once cached. |
| `amax` | `None` | Filters any named network to reactions with `A <= amax` — the single biggest lever for LT-era solve time on the `large` network (see the `large` vs `large, amax=8` gap above). |
| `n_jobs` (`run_mc`) | `-1` (all cores) | MC samples are embarrassingly parallel (independent solves); `n_jobs=1` is useful for reproducible profiling of a single sample's cost, not for a production MC run. |
| `force_backend` | `None`/`"auto"` (prefers C) | The only `PRIMAT.__init__` feature the C backend cannot express is `background=` (a custom `Background` object) — seeing an unexpected Python-backend fallback in `log_backend=True` output usually means that is set. (`extra_rho` and `decay_era` *are* supported on the C backend.) |

The high-precision reference-run settings documented in `CLAUDE.md`
(`numerical_precision=1e-10`, `sampling_temperature_per_decade=2000`,
`sampling_nTOp_per_decade=125`, `rate_grid_npts=4000`,
`runfiles/primat_reference_run.py`) trade several minutes of wall time for
the extra digits used to derive the tolerance bands in that file's
validation table — not something to use for routine runs or MC sampling.

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
