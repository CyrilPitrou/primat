# Tutorials

:::{note}
*(stub — FABLEADVICE O-3)* The `notebooks/` gallery will be rendered here via
`myst-nb`. Notebook execution is disabled in the ordinary `-W` docs build (the
committed notebooks ship with stored outputs); CI's nightly lane re-executes
them with `NB_EXECUTION_MODE=cache`.

To wire a notebook into the site, add its path to the toctree below, e.g.
`../../notebooks/StandardPlots.ipynb` (Sonnet's migration step decides whether
to include notebooks in place or copy curated copies into `docs/`).
:::

The tutorial notebooks live in the repository's
[`notebooks/`](https://github.com/CyrilPitrou/primat/tree/master/notebooks)
directory. Highlights:

- **StandardPlots** — the canonical abundance plots.
- **AbundanceEvolution** — time evolution of every nuclide.
- **MonteCarloRates** — nuclear-rate uncertainty propagation.
- **Sensitivity** — response of observables to input parameters.
- **CompareSmallNetworks** — `small` vs `small_parthenope`.
- **PosteriorBaryons** — a baryon-density posterior from BBN + CMB.

```{toctree}
:maxdepth: 1
:hidden:
```
