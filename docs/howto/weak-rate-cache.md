# The n↔p weak-rate cache workflow

The n↔p weak rates are the most expensive part of initialisation. The
non-thermal rate — Born plus the finite-mass, radiative and
spectral-distortion corrections, `FM`/`CCR`/`SD` in the {doc}`../glossary` —
is cached in
`data/cache_plasma_weak/weak/nTOp_<hash>.txt` (forward and backward columns
together); the finite-temperature radiative correction (CCRTh) is cached
separately in `data/cache_plasma_weak/weak/nTOp_thermal_<hash>.txt`.

Each file is named after a *fingerprint*: a hash of every config field that
affects its numeric content (background thermodynamics,
`sampling_nTOp_per_decade`/`sampling_nTOp_thermal_per_decade`,
`radiative_corrections`, `finite_mass_corrections`, the neutrino degeneracy
`munuOverTnu`, the NEVO table selection, etc. — the authoritative lists are
`_WEAK_RATE_BG_FIELDS` and `_THERMAL_BG_FIELDS` in
`primat/weak_rates/cache.py`, each field carrying a comment on why it is in
or out). The same hash is also written into the file as a
`# fingerprint_hash:`/`# fingerprint:` header, together with the full field
dict — that header is for humans (and migration scripts) to read, not for the
loader, which matches on the *filename* alone. At every run:

- If `weak_rate_cache=True` (default) and a cache file's fingerprint matches
  the current configuration, the corresponding rates are loaded directly —
  initialisation is effectively instantaneous.
- Otherwise (fingerprint mismatch or missing file) the rates are recomputed
  from scratch by numerical integration — 0.4 s on the C backend, 8 s on the
  pure-Python one.
- `weak_rate_cache=False` forces that recompute for the **non-thermal** rate
  only. The thermal (CCRTh) table is loaded whenever a file matching its own
  fingerprint exists, on both backends, so this flag will not make a run
  repeat the expensive integration below. To redo that one, change a field its
  fingerprint covers, or delete the `nTOp_thermal_<hash>.txt` file.
- `save_nTOp` and `save_nTOp_thermal` (both default **`True`**) write the
  (re)computed rates back to `data/cache_plasma_weak/weak/` with a fresh
  fingerprint header, so future runs with the same configuration load the
  cache. The hash is part of the filename, so different configurations
  coexist without overwriting each other — set either flag to `False` only
  to avoid littering `cache_plasma_weak/weak/` during throwaway experiments.

Recomputing the thermal correction (`thermal_corrections=True`) requires a
`vegas` Monte-Carlo integration that can take a few minutes; the fingerprint
mechanism above is what makes this avoidable across runs that share the same
configuration.

Because that recompute is expensive, the thermal fingerprint deliberately
depends on *fewer* fields than the non-thermal one: the CCRTh table is reused
across runs that differ only in something the thermal integral cannot see
(`T_end_MeV`, the spectral-distortion settings, the T-grid density). Anything
the integral *does* see is keyed — including the electron-neutrino degeneracy
`munuOverTnu`/`munuOverTnu_e`, which enters the integrands' neutrino
occupation directly. Each exclusion is justified in a comment next to
`_THERMAL_BG_FIELDS`; if you add a term to the thermal integrand, check
whether it reads a config field that list does not yet cover.

The physical constants are keyed the same way, through a `constants_hash`
field: each cache hashes only the constants it reads
(`cache_utils.CACHE_CONSTANTS`), so `--gA 1.276` or `--me 0.511` re-keys the
rate table while `--T0CMB 2.7250` — which the integrands cannot see — leaves
it valid and skips the multi-minute thermal recompute.
`tests/test_cache_constant_deps.py` proves those lists exact in both
directions, so a term that starts reading a new constant fails the suite
rather than serving a stale table.

## Format version

`WEAK_RATE_FORMAT_VERSION` (`primat/weak_rates/cache.py`, mirrored in
`primat-c/src/cache.c`) is part of both fingerprints, so bumping it
invalidates every cache file at once. Bump it whenever a code change alters
the *numeric content* of these files for an unchanged configuration — a new
term, a changed formula, a new clamp. The shipped tables under
`primat/data/cache_plasma_weak/weak/` must then be re-keyed in the same
commit (their filenames and headers carry the old hash), or every default run
silently misses the cache and pays a fresh integration;
`tests/test_weak_rates.py::test_shipped_weak_caches_carry_current_format_version`
fails if that step is forgotten.

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

Cache files are written under the shipped `primat/data/` tree by default; on
a read-only install, point `cache_dir` at a writable directory instead — see
{doc}`data-overlays`.

See `primat.weak_rates.api.ComputeWeakRates`/`RecomputeWeakRates` for the
full cache-loading algorithm.
