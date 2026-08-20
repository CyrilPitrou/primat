# Installation

## For most users

```bash
pip install primat
```

That's it. The package includes a fast C backend compiled for your platform,
with a pure-Python fallback if no compiled binary is available — both give
identical results, just different speed. To get started, just type

```bash
primat --help
```

## Core dependencies

| Package | Role |
|---------|------|
| `numpy`, `scipy` | **Mandatory** (installed by `pip install primat`) — everything needed for a single `run_bbn`/`solve` |
| `numba` (`recommended`) | JIT compilation gives ~5× speedup on rate kernels |
| `vegas` (`recommended`) | Monte Carlo integration for thermal weak-rate corrections |
| `joblib` (`recommended`, `mc`) | Parallel Monte-Carlo (`run_mc(..., n_jobs != 1)`); a core install can still run MC serially with `n_jobs=1` or on the C backend |

`joblib` and `plotly` are **not** mandatory — they live in extras so a lean
`pip install primat` core (numpy + scipy) is enough to install and solve.

## Optional extras

Pick an extra by use case:

| Extra | `pip install "primat[…]"` | Adds |
|-------|---------------------------|------|
| `recommended` | `numba`, `vegas`, `joblib` | fast kernels + thermal-rate integrator + parallel MC |
| `mc` | `joblib` | parallel Monte-Carlo only |
| `plots` | `plotly` | interactive plotting (used by the GUI figures) |
| `gui` | `streamlit`, `pandas`, `plotly` | the `primat-gui` web app |
| `notebooks` | `matplotlib`, `pandas`, `papermill` | the example notebooks |
| `docs` | `sphinx`, `furo`, `myst-nb`, ... | building this documentation site |
| `all` | everything above except `docs` | one-shot full install |

```bash
pip install "primat[recommended]"
```

### GUI

```bash
pip install "primat[gui]"
primat-gui
```

| Package | Role |
|---------|------|
| `streamlit` | **Required for `primat-gui`** — the web app framework |
| `pandas` | **Required for `primat-gui`** — final-abundance table |
| `plotly` | **Required for `primat-gui`** — the abundance-evolution figures |

See {doc}`howto/custom-networks` for what the GUI adds beyond the CLI/API.

### Notebooks

```bash
pip install -e ".[notebooks]"
```

| Package | Role |
|---------|------|
| `matplotlib`, `pandas` | Plotting and tabular display in the notebooks |
| `papermill` | Headless notebook re-execution (used by the nightly docs build) |

See {doc}`tutorials/index` for the notebook gallery.

## From source

```bash
git clone https://github.com/CyrilPitrou/primat
cd primat
pip install -e ".[recommended]"
```

For users who prefer to work directly with the C code, `primat-c/` is a
standalone C99 implementation that builds independently via `make` — see
`primat-c/README.md` in the repository for compilation instructions.

If you intend to change the code rather than only run it, read
[`docs/development.md`](https://github.com/CyrilPitrou/primat/blob/master/docs/development.md)
in the repository: the conventions this code is written to, what the test
suite enforces, and the list of numerical decisions already measured and
settled. It is addressed to whoever edits the source, so it is not part of
this site.
