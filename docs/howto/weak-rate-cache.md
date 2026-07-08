# The n↔p weak-rate cache workflow

The n↔p weak rates are the most expensive part of initialisation (~1.8 s).
The non-thermal rate (Born+FM+CCR+SD) is cached in
`data/weak/nTOp_<hash>.txt` (forward and backward columns together); the
finite-temperature radiative correction (CCRTh) is cached separately in
`data/weak/nTOp_thermal_<hash>.txt`.

Each file is tagged with a *fingerprint* header: a hash of every config field
that affects its numeric content (background thermodynamics,
`sampling_nTOp_per_decade`/`sampling_nTOp_thermal_per_decade`,
`radiative_corrections`, `finite_mass_corrections`, `thermal_corrections`,
etc. — see `primat.weak_rates`). At every run:

- If `weak_rate_cache=True` (default) and a cache file's fingerprint matches
  the current configuration, the corresponding rates are loaded directly —
  initialisation is effectively instantaneous.
- Otherwise (fingerprint mismatch, missing file, or `weak_rate_cache=False`),
  the rates are recomputed from scratch by numerical integration (~1.8 s).
- `save_nTOp` and `save_nTOp_thermal` (both default **`True`**) write the
  (re)computed rates back to `data/weak/` with a fresh fingerprint header, so
  future runs with the same configuration load the cache. The hash is part
  of the filename, so different configurations coexist without overwriting
  each other — set either flag to `False` only to avoid littering
  `data/weak/` during throwaway experiments.

Recomputing the thermal correction (`thermal_corrections=True`) requires a
`vegas` Monte-Carlo integration that can take a few minutes; the fingerprint
mechanism above is what makes this avoidable across runs that share the same
configuration.

## Typical workflow for a high-precision study

```python
from primat.backend import run_bbn

# Step 1 -- compute and save high-precision rates once (non-default
# sampling_nTOp_per_decade gives a fingerprint that the shipped cache won't
# match, so this recomputes; save_nTOp=True is the default)
result1 = run_bbn({"save_nTOp": True, "sampling_nTOp_per_decade": 160})

# Step 2 -- all subsequent runs with the same sampling_nTOp_per_decade reuse
# the saved tables
result2 = run_bbn({"sampling_nTOp_per_decade": 160})
```

See `primat.weak_rates.api.ComputeWeakRates`/`RecomputeWeakRates` for the
full cache-loading algorithm.
