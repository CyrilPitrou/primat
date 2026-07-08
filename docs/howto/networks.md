# Choose a network and `amax`

:::{note}
*(stub — FABLEADVICE O-3)* Migrate from the README "Key configuration flags"
table and CLAUDE.md's network discussion: `small` / `small_parthenope` /
`large`, the `amax` A≤N filter (e.g. `network="large", amax=8` = the old
"medium" 68-reaction network), the shared HT/MT eras, and the accuracy notes
for the heavy-nuclide tail.
:::

Set the network at construction time:

```python
from primat import PRIMAT

bbn = PRIMAT({"network": "large", "amax": 8})
results = bbn.solve()
```

| `network` | Reactions | Notes |
|-----------|-----------|-------|
| `small` (default) | 12 | fast, light elements only |
| `small_parthenope` | 12 | Parthenope 3.0 rate tables |
| `large` | ~429 | full network, ~59 nuclides |
| any name | — | loads `data/nuclear/networks/<name>.txt` |

The `amax` flag filters *any* network to reactions with A ≤ `amax`.
