---
sd_hide_title: false
---

# primat

**primat** is a precise Big Bang Nucleosynthesis (BBN) solver, distributed as a
single pip-installable package with two interchangeable backends: a fast C
engine (the default) and a pure-Python implementation used as a fallback and for
development. It integrates the coupled ODEs of the cosmological background
(photon/neutrino temperatures, scale factor) together with a nuclear reaction
network to predict the primordial abundances of H, D, ³He, ⁴He, ⁷Li and heavier
nuclides.

:::{note}
This documentation site is being migrated from the project's README, LaTeX
manual, and notebook galleries. Pages marked *(stub)* are scaffolding for the
in-progress content migration (FABLEADVICE O-3).
:::

## Quick start

```bash
pip install primat
```

```python
from primat import PRIMAT

# A standard-model run at the Planck baryon density.
bbn = PRIMAT({"Omegabh2": 0.022425, "network": "small"})
results = bbn.solve()
print(f"YP  = {results['YPBBN']:.8f}")
print(f"D/H = {results['DoH']:.7e}")
```

Or straight from the command line:

```bash
primat --Omegabh2 0.022425 --network large --amax 8
```

## Citing primat

If you use primat in published work, please cite:

> Pitrou, Coc, Uzan, Vangioni, *Physics Reports* **04** (2018) 005
> ([arXiv:1801.08023](https://arxiv.org/abs/1801.08023)).

See {doc}`citing` for the BibTeX entry.

```{toctree}
:maxdepth: 2
:caption: Getting started
:hidden:

installation
tutorials/index
```

```{toctree}
:maxdepth: 2
:caption: How-to guides
:hidden:

howto/index
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

physics
api/index
cli
```

```{toctree}
:maxdepth: 1
:caption: About
:hidden:

changelog
citing
```
