# Custom networks (GUI zip and API)

There are two ways to customise the reaction set: interactively through the
GUI, or permanently by adding a rate table + network entry to the source tree.

## Interactively, via the GUI

After `pip install "primat[gui]"`, `primat-gui`'s sidebar "Nuclear reactions"
group offers two buttons:

- **"Create custom network"** — a popup to start from any named network,
  toggle reactions in/out grouped by mass-number category
  (`reaction_category`/`group_reactions_by_category`), and substitute or
  upload an alternate rate table per reaction, or add brand-new reactions.
- **"Import custom network"** — re-load a previously exported `.zip`.

Either save the result as a `.zip` (re-importable later) or apply it and run
BBN directly. The exported zip always contains the *reinterpolated*
(on-grid) version of every user-supplied table, so it stays consistent with
whatever `rate_grid_npts`/`rate_grid_T9_min`/`rate_grid_T9_max` the run used.

## Programmatically, via `custom_network=`

The GUI's export format is also the API's input format: pass a dict with the
`{"removed": [...], "replaced": {...}, "added": {...}}` schema to
`run_bbn`/`PRIMAT` directly — see `primat.network_data.UpdateNuclearRates`
for the full schema, and `primat.main.PRIMAT.__init__` for worked examples.
It is supported on **both** backends (not one of the Python-only features
listed in {doc}`../api/backend`).

```python
from primat import PRIMAT

# Drop one reaction, override another's rate table, and add a brand-new
# reaction (its stoichiometry is read from the name):
PRIMAT({"network": "small"}, custom_network={
    "removed": ["d_d__t_p"],
    "replaced": {"n_p__d_g": "0.001 1.2e3\n10.0 4.5e1\n"},
    "added": {"t_t__He4_n_n": "0.001 1.0e2\n10.0 1.0e2\n"},
})
```

`custom_network` is not a `PRIMATConfig` field (it carries bulk table data
rather than a fingerprintable scalar), so it does not participate in any
rate-cache fingerprint.

## Permanently, via the source tree

For a reaction you want to keep long-term rather than toggle at runtime, add
it directly:

1. Drop a rate table under `primat/data/nuclear/tables/<name>/<name>.txt`.
2. Add `<name>` to the relevant network file under
   `primat/data/nuclear/networks/`.
3. `reaction_stoichiometry` (`primat.network_data`) auto-derives the
   stoichiometry from the reaction name's `TO`-separated tokens, falling
   back to a manual `reactions_large.csv`/`detailed_balance.csv` row only if
   the name can't be tokenised. `load_network` validates A/Z conservation
   and rejects duplicate entries, so a malformed addition fails fast and
   loudly rather than silently mis-integrating.

For a one-off sensitivity study rather than a permanent addition, use the
existing `p_<reaction>`/`delta_<reaction>` config knobs instead — see
{doc}`rate-variation-mc` — no file changes needed.
