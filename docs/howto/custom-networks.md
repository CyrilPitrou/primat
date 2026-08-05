# Custom networks (GUI zip and API)

There are two ways to customise the reaction set: interactively through the
GUI, or permanently by adding a rate table + network entry to the source tree.

## Interactively, via the GUI

After `pip install "primat[gui]"`, `primat-gui`'s sidebar "Nuclear reactions"
group offers a single **"Manage networks"** button — the one gateway to every
network action. The dialog it opens lists the networks built or imported this
session (select / rename / remove), loads one from a previously exported
`.zip`, and hands off to:

- **"Create new network"** — start from any named network, toggle reactions
  in/out grouped by mass-number category
  (`reaction_category`/`group_reactions_by_category`), substitute or upload
  an alternate rate table per reaction, override a decay rate, or add
  brand-new reactions.

That dialog's footer only *saves* the network (and offers it as a
re-importable `.zip`); running it is left to the main "Run BBN" button, the
same as for an imported one.

### What the exported zip contains

Every user-supplied table is written **verbatim — on its own original grid,
at full `%.17e` precision** — deliberately *not* pre-resampled onto the
master T9 grid. That is what makes a round trip bit-for-bit: `load_network`'s
resampling then runs exactly once, on the same data the GUI's own live run
resampled. Pre-resampling at export would round the values and extrapolate a
coarse upload onto the wider master grid, which a re-import would then
resample a *second* time, drifting ~1e-6 from the run being reproduced.

A reaction you did **not** customise is copied verbatim from the table its
network actually pins — `small_parthenope`'s `*_parthenope3.0.txt`, say —
never from an assumed `<name>_primat.txt`.

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
