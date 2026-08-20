# Data overlays: `data_dir`, `user_nuclear_dir` and `cache_dir`

The shipped default `data/` tree lives inside the `primat` package
(`primat/data/`). Three `DEFAULT_PARAMS` fields let you override where data is
read from and written to — one full replacement, two additive overlays:

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

The example above works on either backend: `data_dir` is an ordinary `params`
key, which `primat.backend` additionally hands to the C extension as its data
root (the C side loads `csv/nuclides.csv` before any parameter is applied, so it
needs the directory up front).

On the C CLI, the equivalent is the `--data_dir` flag, the
`CPRIMAT_DATA_DIR` environment variable, or a `data_dir = …` line in an
`--ini` file.

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
always read from the data root on both backends). The regenerable caches have
their own overlay, below.

## `cache_dir` — where the regenerable caches are written

The n↔p weak-rate cache and the QED/electron-thermo plasma caches
(`cache_plasma_weak/{weak,plasma}/`, see {doc}`weak-rate-cache`) are the only
data primat *writes*. By default they go into the shipped tree, which fails on
a read-only install — a system-wide `pip install`, a container image, a
read-only home. `cache_dir` moves both:

```python
result = run_bbn({"cache_dir": "~/.cache/primat"})
```

Caches are then written to `<cache_dir>/{weak,plasma}/` and read from there
first, falling back to the shipped copies on a miss — so the tables that ship
with the package are never shadowed, only added to. Unlike the other two
fields, `cache_dir` is *not* validated at construction: it is a write target,
created on demand, so a directory that does not exist yet is normal. A failed
write warns (naming `cache_dir`) rather than raising, since a cache is an
optimisation and the run can proceed without it.

`cache_dir` is deliberately absent from every cache fingerprint: where a table
is stored cannot change what is in it.

## Common ground

- All three default to `None`. `data_dir` and `user_nuclear_dir` are eagerly
  validated as existing directories at construction time (`cache_dir` is not,
  for the reason above). Because it is a *takeover*, `data_dir` is
  additionally checked for the `csv/` and `nuclear/` subdirectories — before
  `nuclides.csv` is read from it, so a typo'd path is reported as a bad
  `data_dir` rather than as a missing CSV. A leading `~` is expanded in all
  three, on both backends.
- `user_nuclear_dir` is not part of the n↔p weak-rate fingerprint machinery
  (`_WEAK_RATE_BG_FIELDS`/`_THERMAL_BG_FIELDS`) since it only affects
  network/rate-table resolution, not anything those fingerprints cover.
- The C side (`primat-c/`) mirrors this overlay exactly — same lookup
  order, same two wired call sites (`cpr_load_network` and each reaction's
  rate-table open) via `cpr_config_resolve_rates_path`.
