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

## Quick start

```bash
pip install primat
```

```python
from primat.backend import run_bbn

result = run_bbn({"Omegabh2": 0.02242})

print(f"YP  (BBN) = {result['YPBBN']:.8f}")  # 0.24699907
print(f"D/H       = {result['DoH']:.7e}")    # 2.4358767e-05
```

```text
[primat]  HT.  MT.  LT.  done.
YP  (BBN) = 0.24699907
D/H       = 2.4358767e-05
```

New to BBN codes? {doc}`tutorials/first-run` walks through what those two
numbers mean and how to change them.

(Those are the C backend's values for the default `small` network; the
pure-Python backend gives `0.24699896` / `2.4358605e-05`, within the
cross-backend tolerance. The authoritative reference values, with the
tolerance that applies to each, live in
[`tests/README.md`](https://github.com/CyrilPitrou/primat/blob/master/tests/README.md)'s
"Validation reference".)

`run_bbn()` is the main entry point and automatically selects the best
available backend (fast C engine by default, pure-Python fallback if
needed). Pass an optional parameter dict to override defaults; {doc}`parameters`
lists all 95 keys with their defaults, and every one of them can equally be set
on the command line or in an INI file.

## Four ways to use primat

All four produce identical results (to the cross-backend tolerance described
in {doc}`api/backend`).

1. **Python API** (recommended) — `primat.backend.run_bbn(params)`, shown
   above. Force a specific backend with `force_backend="c"` /
   `force_backend="python"`.
2. **Command-line interface**:
   ```bash
   primat --Omegabh2 0.02242 --network large --amax 8
   ```
   See {doc}`cli` for the full flag reference (auto-generated from
   `primat --help`); anything not exposed as a flag is reachable via the
   repeatable `--set KEY=VALUE` escape hatch.
3. **Graphical interface** — after `pip install "primat[gui]"`, run
   `primat-gui` for a browser-based parameter form, interactive
   abundance-evolution plot, and final-abundances panel. See {doc}`howto/gui`;
   it is also the only interface that can build a network interactively
   ({doc}`howto/custom-networks`).
4. **Example scripts** (development/source-only) — `python
   runfiles/primat_run.py`, `primat_compare.py`, `primat_reference_run.py`,
   `primat_mc.py`; run from a source checkout's repo root.

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
tutorials/first-run
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
extending
performance
parameters
api/index
cli
glossary
```

```{toctree}
:maxdepth: 1
:caption: About
:hidden:

changelog
citing
```
