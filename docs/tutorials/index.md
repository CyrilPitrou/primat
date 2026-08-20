# Worked examples

Ten notebooks, each a complete piece of analysis: the standard figures,
parameter scans, rate plots and uncertainty budgets. They assume you know what
a BBN run is and what it produces — {doc}`first-run` is the guided
introduction if you do not.

They are rendered directly from the repository's
[`notebooks/`](https://github.com/CyrilPitrou/primat/tree/master/notebooks)
directory and ship with their stored outputs, so the plots on this site are
the ones their author saw.

To run any notebook yourself, clone the repo and install the `notebooks`
extra:

```bash
git clone https://github.com/CyrilPitrou/primat
cd primat
pip install -e ".[notebooks,recommended]"
jupyter lab notebooks/
```

## Standard results

- **{doc}`StandardPlots`** — the Schramm diagram: primordial abundances vs.
  η_b with 1σ nuclear-rate uncertainty bands and observational constraints
  (YP, D/H, ³He/⁴He, ⁷Li/H).
- **{doc}`AbundanceEvolution`** — time evolution of `A_i Y_i(t)` for every
  nuclide from 1 s to 10⁵ s, for both the small (12-reaction) and large
  (~429-reaction) networks.
- **{doc}`CompareSmallNetworks`** — `small` vs. `small_parthenope` head to
  head.
- **{doc}`AnimatedAbundances`** — animated GIFs of the small-network
  abundance evolution `A_i Y_i(t)`, sweeping ΔNeff and Ω_b h² in turn.

## Nuclear rates

- **{doc}`ReactionRates`** — plots the tabulated rate ⟨σv⟩(T9) of any
  reaction in the network, with the master-grid reinterpolation overlaid,
  alongside the n↔p weak rates, all compared against the Hubble rate H(T)
  to show freeze-out.

## Parameter scans

- **{doc}`PosteriorBaryons`** — a posterior on Ω_b h² from YP and D/H:
  scans Ω_b h² ∈ [0.020, 0.024] and computes Gaussian likelihoods from each
  observable.
- **{doc}`AbundancesNrelat`** — abundances vs. ΔNeff: scans
  ΔNeff ∈ [−2, +2] to show how extra relativistic species shift YP and D/H.
- **{doc}`AbundancesXi`** — abundances vs. neutrino degeneracy
  ξ = μ_ν/T_ν: scans ξ ∈ [−0.05, +0.05].

## Uncertainty analysis

- **{doc}`MonteCarloRates`** — the full MC uncertainty budget: draws nuclear
  rates, τ_n, and Ω_b h² simultaneously, with histograms and a corner plot of
  the joint distribution of every observable. See also
  {doc}`../howto/rate-variation-mc`.
- **{doc}`Sensitivity`** — sensitivity tables: the logarithmic derivative
  ∂ ln(observable) / ∂ ln(parameter) for each of the 12 nuclear rates, τ_n,
  G_N, Ω_b h², and ΔNeff, as formatted tables and a heat-map.

:::{note}
Common conventions across the parameter-scan notebooks: a fixed MC seed
(`MC_SEED = 0`) at every grid point so finite-sample MC bias cancels across
the grid; observational constraints shown as grey bands; the baryon
density Ω_b h² = 0.02242 ± 0.00014 (Planck 2018 + BAO) shown as a red vertical band; and
`num_mc = 500`+ for publication-quality uncertainty bands (the notebooks
default to 50 for speed).
:::

```{toctree}
:maxdepth: 1
:hidden:

StandardPlots
AbundanceEvolution
CompareSmallNetworks
AnimatedAbundances
ReactionRates
PosteriorBaryons
AbundancesNrelat
AbundancesXi
MonteCarloRates
Sensitivity
```
