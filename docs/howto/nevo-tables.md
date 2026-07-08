# Custom NEVO (neutrino-decoupling) tables

The neutrino-decoupling history is read from `data/NEVO/`. Three optional
parameters point at alternative tables instead (filenames resolved relative
to `data/NEVO/`, or absolute paths); each defaults to `None`, which keeps
the shipped file selected by `QED_corrections`/`incomplete_decoupling`:

| Parameter | Replaces | Expected shape |
|-----------|----------|-----------------|
| `nevo_file` | `NEVOPRIMAT[_NoQED]_col_1_7.csv` | 6 or 7 columns |
| `nevo_spectral_file` | `NEVOPRIMAT[_NoQED].csv` (only read when `spectral_distortions=True` and `analytic_distortions=False`) | 6 thermo columns + N spectral columns (86 in the shipped tables) |
| `nevo_grid_file` | `NEVOGrid.csv` | length = N (the spectral table's column count minus 6) |

```python
from primat.backend import run_bbn

result = run_bbn({"nevo_file": "my_NEVO_col_1_7.csv"})
```

A custom file is checked for existence and shape at construction time
(`PRIMATConfig.__init__` raises `ValueError` early on a mismatch), and is
included in the n↔p weak-rate cache fingerprint so pointing at a different
table correctly triggers a recompute (see {doc}`weak-rate-cache`).
`nevo_grid_file` is currently assumed to be a Gauss-Laguerre quadrature; the
format for handing NEVO results to primat will evolve in future releases.

## `nevo_file_prefix`

A single knob that rebuilds all four default filenames at once:
`<prefix>[_NoQED]_col_1_7.csv`, `<prefix>[_NoQED].csv`, and the shared
`NEVOGrid.csv` (not prefixed). It is ignored for any file selected
explicitly via `nevo_file`/`nevo_spectral_file` (those still win), and has
no effect when `incomplete_decoupling=False`. A custom prefix is your
responsibility to keep self-consistent with the shipped conventions;
`PRIMATConfig.__init__` still checks the derived files exist with the
expected shape.

```python
result = run_bbn({"nevo_file_prefix": "MyNEVO"})
```

## `external_scale_factor`

When `True`, `a(T_γ)` is read directly from the NEVO table's `x` column
(`x ∝ a` by the NEVO convention) instead of solved from the
entropy-conservation ODE driven by the NEVO heating function `N_NEVO(T_γ)`.
`t(a)` is unchanged (Hubble integration in both modes). Requires
`incomplete_decoupling=True`. The two modes agree to ~1e-5 in `a(T)`/`t(T)`
and reproduce the same `Neff`/`YPBBN`/`D/H` to the quoted precision.

```python
result = run_bbn({"incomplete_decoupling": True, "external_scale_factor": True})
```

This mode is mutually exclusive with {doc}`backgrounds`'s
`custom_background` mechanism.
