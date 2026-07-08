# Installation

:::{note}
*(stub — FABLEADVICE O-3)* This page will absorb the README's install section:
wheels, the optional extras (`recommended` / `mc` / `plots` / `gui` /
`notebooks` / `all`), from-source builds of the C backend, and conda-forge once
it exists.
:::

## From PyPI

```bash
pip install primat
```

This installs the lean core (numpy + scipy) with the fast C backend compiled
where a wheel is available, and a graceful pure-Python fallback otherwise.

## Optional extras

| Extra | Pulls in | For |
|-------|----------|-----|
| `recommended` | numba, vegas, joblib | JIT kernels, thermal weak-rate MC integration, parallel Monte-Carlo |
| `mc` | joblib | parallel Monte-Carlo only |
| `plots` | plotly | interactive GUI figures |
| `gui` | streamlit, pandas, plotly | the `primat-gui` app |
| `notebooks` | matplotlib, pandas, papermill | running the tutorial notebooks |
| `all` | everything above | one-shot full install |
| `docs` | sphinx + theme + extensions | building this documentation |

```bash
pip install "primat[recommended]"
```

## From source

```bash
git clone https://github.com/CyrilPitrou/primat
cd primat
pip install -e ".[recommended]"
```
