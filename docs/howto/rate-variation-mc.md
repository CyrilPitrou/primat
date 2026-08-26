# Rate variation and Monte-Carlo uncertainty

primat provides two distinct mechanisms for varying nuclear reaction rates,
plus a Monte-Carlo driver that uses the first one to propagate uncertainty
onto any observable.

## 1. Log-normal rate variations: `p_<reaction>` parameters

Each nuclear reaction has a corresponding parameter `p_<name>` (e.g.
`p_n_p__d_g` for the n + p → d + γ reaction). This varies the rate as:

**Rate = median × exp(p × σ)**

where σ is the rate's log-normal uncertainty width (from the rate table's
error column).

- **Primary use case**: Monte-Carlo uncertainty propagation (below) — use
  `run_mc()`/`mc_uncertainty()` to automatically sample `p_*` from N(0,1) for
  each reaction.
- **Manual use**: set `p_<name>` directly to explore a fixed variation, e.g.
  `p_n_p__d_g = 1` increases the rate by roughly +1σ, `p_n_p__d_g = -2`
  decreases it by roughly −2σ. Also settable from the CLI:
  ```bash
  primat --set p_n_p__d_g=1
  ```

## 2. Additive rate rescaling: `delta_<reaction>`

For deterministic sensitivity studies, set `delta_<name>`. When `p_<name>=0`
(the default), the rate becomes:

**Rate = median × (1 + delta_&lt;name&gt;)**

This allows uniform or per-reaction rescaling. When `p_<name>≠0` as well, the
combined formula is:

**Rate = median × (exp(p × σ) + delta_&lt;name&gt;)**

```python
from primat.backend import run_bbn

# Sensitivity study: vary n+p->d+gamma rate by +10%
result = run_bbn({
    "delta_n_p__d_g": 0.1,
})
```

:::{important}
`p_<reaction>` is designed for MC uncertainty propagation (log-normal
variations); `delta_<reaction>` is designed for deterministic sensitivity
studies (additive variations). They can be used together, but interpret the
combined effect carefully.
:::

## 3. Computing the uncertainty: `run_mc()` and `--mc N`

`run_mc()` (or its pure-Python counterpart `primat.main.mc_uncertainty()`)
computes the propagated nuclear-rate/τ_n uncertainty on any observable: it
runs many independent BBN solves, each with randomly-sampled reaction rates
(and neutron lifetime), and reports the spread of results as the
uncertainty.

```python
from primat.backend import run_mc

mc = run_mc(100, ["YPBBN", "DoH"], params={"Omegabh2": 0.02242})

mc["DoH"].central   # nominal (best-estimate) value
mc["DoH"].mean      # mean over the 100 MC samples
mc["DoH"].std       # MC uncertainty (1-sigma, sample std/ddof=1) -- "the error"
mc["DoH"].values    # full array of per-sample values, length 100
```

Regardless of which `quantities` you ask for, the result also includes every
standard observable the run produced — `Neff`, `YPBBN`, `YPCMB`, `DoH`,
`He4oH`, `He3oH`, `He3oHe4`, `Li7oH`, plus `Li6oLi7`/`YCNO` on the networks
that track those nuclides — and every tracked nuclide's final abundance, at no
extra cost. `mc.quantity_names()` lists what a given run actually has; on the
default `small` network that is the eight observables above and eight
nuclides. See {doc}`output` for which keys are conditional and why.

### Joint uncertainty (covariance and correlation)

When you constrain cosmology with several abundances at once (typically
`YPBBN` *and* `DoH`), you need how they *co-vary* across the same MC
samples, not just their individual σ's. `MCResult` exposes both matrices
(sample estimators, `ddof=1`, so `diag(cov) == std**2`):

```python
mc.cov()                  # full (n_q, n_q) covariance matrix, quantity_names() order
mc.corr()                 # full correlation matrix (unit diagonal)
mc.cov("YPBBN", "DoH")    # scalar covariance between two named quantities
mc.corr("YPBBN", "DoH")   # scalar correlation, e.g. -0.13 (mildly anti-correlated)
```

The `YPBBN`–`DoH` correlation is *weak* because the two are set by largely
different physics: `YPBBN` is fixed by the n/p ratio at freeze-out (τ_n, the
weak rates) while `D/H` is fixed by the deuterium-burning reactions. The
strongly correlated pairs are the ones sharing a burning chain — e.g.
`DoH`–`He3oHe4` at about −0.4.

A quantity identical in every sample (zero variance) gives `NaN`
off-diagonal correlations (never a warning storm).

### From the command line

Add `--mc N`; the summary prints the 4×4 correlation and covariance matrices
of the four main products (`YPBBN`, `DoH`, `He3oHe4`, `Li7oH`):

```bash
primat --Omegabh2 0.02242 --mc 300 --mc-seed 0
# YP (BBN)   = 0.24699907 +/- 0.00010631
# D/H        = 2.4358767e-05 +/- 2.6358624e-07
# ...
# Correlation matrix (YPBBN, DoH, He3oHe4, Li7oH):
#             YPBBN      DoH  He3oHe4    Li7oH
#    YPBBN    1.000   -0.129   -0.064    0.071
#      DoH   -0.129    1.000   -0.421   -0.388
#  He3oHe4   -0.064   -0.421    1.000    0.460
#    Li7oH    0.071   -0.388    0.460    1.000
```

Those numbers are a real run (300 samples, `seed=0`, C backend, default
`small` network), not an illustration — but they are still a *finite* sample:
with N = 300 each σ carries ~4 % statistical error of its own and the
off-diagonal correlations rather more, so expect the last digits to move
between seeds. The value before each `+/-` is the unperturbed (central)
solve, so it depends on neither N nor the seed. Raise `--mc` for anything you
intend to quote.

`--mc-seed` sets the random seed (use the same seed to reproduce a run) and
`--mc-jobs` the number of parallel workers. The three MC output files share
one filename stem `--output_mc_file_prefix PREFIX` (default
`results/output_mc`), each gated by its own flag:

- `--output_mc_samples` → `PREFIX_samples.tsv` (every raw per-sample value,
  one column per quantity)
- `--output_mc_covariance` → `PREFIX_covariance.tsv` (the covariance matrix)
- `--output_mc_correlation` → `PREFIX_correlation.tsv` (the correlation
  matrix)

All three, and both printed matrices, are identical (same shape, same
header wording) whether the C or the pure-Python backend runs.
Programmatically the writers are `primat.backend.dump_mc_samples` /
`dump_mc_covariance` / `dump_mc_correlation`.

See `runfiles/primat_mc.py` for a runnable end-to-end demo, and the
{doc}`../tutorials/MonteCarloRates` notebook for the full uncertainty budget
(nuclear rates + τ_n + Ω_b h² simultaneously) with corner plots.

## Backend parallelisation strategies

Both backends support MC with parallel acceleration, using different
strategies:

**Python backend** (`primat.main.mc_uncertainty`):
- Parallelises MC samples across worker processes using `joblib.Parallel`.
- `n_jobs`: `-1` uses all available CPU cores, or an explicit integer count.
- Each worker draws independent rate/τ_n samples and runs a full PRIMAT
  solve; results are aggregated in the main process.
- Uses NumPy's `default_rng` for the per-sample RNG.

**C backend** (`cpr_mc_uncertainty`):
- Parallelises internally using POSIX threads (pthread) within a single
  process.
- `n_jobs` is passed to the C extension and controls the thread count.
- Uses xoshiro256** (a fast parallel RNG) for per-thread sample generation.
- Typically faster than joblib for small-to-moderate sample counts
  (N < a few thousand) due to lower inter-process overhead.

The two backends produce **statistically equivalent but not bit-for-bit
identical** MC samples — different RNG streams mean individual samples
differ, but mean/std/correlations converge to the same values as N_MC
increases (within numerical precision ~1e-7 to 1e-8). When comparing MC
results across backends, verify *statistics* agree rather than expecting
exact sample values.
