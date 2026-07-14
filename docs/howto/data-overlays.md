# Data overlays: `data_dir` and `user_nuclear_dir`

The shipped default `data/` tree lives inside the `primat` package
(`primat/data/`). Two `DEFAULT_PARAMS` fields let you override where data is
read from — one full-replacement, one additive:

## `data_dir` — full replacement

When set, `data_dir` completely replaces the entire `primat/data/` tree. The
supplied directory must contain `NEVO/`, `nuclear/`, `csv/`, and
`cache_plasma_weak/` subdirectories (the n↔p weak-rate and plasma caches
live together under `cache_plasma_weak/{weak,plasma}/`) — every data file is
then read from there instead of the package's shipped copy.

```python
from primat.backend import run_bbn

result = run_bbn({"data_dir": "/path/to/my_data_tree"})
```

On the C CLI, the equivalent is the `--data_dir` flag or the
`CPRIMAT_DATA_DIR` environment variable.

## `user_nuclear_dir` — additive overlay

`user_nuclear_dir` points at a directory with the same `networks/` and/or
`tables/<name>/` layout as the shipped `data/nuclear/` folder. Any network
file or per-reaction table found there is used *instead of* the shipped
one, while everything not overridden still falls back to the shipped
default — an additive overlay, not a takeover, so `small`/`large` remain
accessible even when `user_nuclear_dir` doesn't contain them.

```python
result = run_bbn({"user_nuclear_dir": "/path/to/my_nuclear_overlay"})
```

Currently wired through this resolver: the `network` validation/loading
path (`nuclear/networks/<name>.txt`), `available_rate_tables` (the GUI's
per-reaction table dropdown), and each individual reaction's rate-table
file — resolved per-file, so one overlaid table doesn't shadow the rest of
the shipped `tables/` tree.

Not yet routed through either override: NEVO tables (which have their own
mechanism — see {doc}`nevo-tables`), or the reaction catalog
(`nuclides.csv`/`reactions_large.csv`/`detailed_balance.csv`/`decays.txt`,
always read from the data root on both backends). The n↔p weak-rate cache
and the QED/electron-thermo plasma caches (`cache_plasma_weak/{weak,plasma}/`)
have their own separate additive overlay, `cache_dir` — see
{doc}`weak-rate-cache`.

## Common ground

- Both fields default to `None` and are eagerly validated as existing
  directories at construction time.
- `user_nuclear_dir` is not part of the n↔p weak-rate fingerprint machinery
  (`_WEAK_RATE_BG_FIELDS`/`_THERMAL_BG_FIELDS`) since it only affects
  network/rate-table resolution, not anything those fingerprints cover.
- The C side (`primat-c/`) mirrors this overlay exactly — same lookup
  order, same two wired call sites (`cpr_load_network` and each reaction's
  rate-table open) via `cpr_config_resolve_rates_path`.
