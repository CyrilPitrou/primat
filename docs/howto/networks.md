# Choose a network and `amax`

Set the network at construction time — via the `network` parameter (Python
API/CLI) or the `--network`/`--amax` flags:

```python
from primat.backend import run_bbn

result = run_bbn({"network": "large", "amax": 8})
```

```bash
primat --network large --amax 8
```

## Available networks

Two named networks (plus a Parthenope-rates variant of the small one) are
available via `network`; `amax` (any positive integer) further restricts
*any* of them to reactions whose nuclides all have mass number A ≤ `amax`:

| `network` | Reactions (nuclear + n↔p) | Nuclides | Notes |
|-----------|---------------------------|----------|-------|
| `"small"` (default) | 12 + 1 = 13 | 8 | the key reactions; fastest |
| `"small_parthenope"` | 12 + 1 = 13 | 8 | same reactions, Parthenope 3.0 rate tables (comparison runs) |
| `"large"` | 428 + 1 = 429 | 59 | from the AC2024 compilation; LT era only |
| `"large"`, `amax=8` | 67 + 1 = 68 | 12 | the old "medium" network's exact equivalent |
| `"large"`, `amax=2` | 2 + 1 = 3 | 3 | the old "deuterium" network's equivalent (n↔p + n_p__d_g + p_p_n__d_p) |
| any other name | — | — | loads `data/nuclear/networks/<name>.txt` |

The first number counts the *nuclear* reactions (what a network file lists,
or what survives the `amax` filter); the total is what
`load_network(...).n_reac` reports, since every network additionally carries
the n↔p weak reaction, which no file lists. Elsewhere in these docs `small`
is called "the 12-reaction network" and `large` "~429 reactions" — the same
two networks, counted the two different ways. Nuclide counts include `n` and
`p`.

All networks share the HT (n↔p) and MT eras — the MT era always uses a fixed
18-reaction subset, too stiff to run the full network; only the LT reaction
set is filtered by `network`/`amax`. The light-element abundances of the
full `large` network match the `amax=8` restriction to ≲1e-4; its
heavy-nuclide tail (B, C, N, O, …) is approximate (limited by the AC2024
rate floors). See the {doc}`../tutorials/AbundanceEvolution` notebook for
evolution plots across networks.

## Related physics-configuration flags

| Flag | Default | Effect |
|------|---------|--------|
| `radiative_corrections` | `True` | Coulomb + T=0 resummed radiative corrections to n↔p (CCR) |
| `finite_mass_corrections` | `True` | Fokker-Planck finite-nucleon-mass correction (FM) |
| `thermal_corrections` | `True` | Finite-temperature radiative corrections to n↔p (CCRTh) |
| `spectral_distortions` | `True` | Correct n↔p rates for non-FD neutrino distributions (SD) |
| `tau_n_normalization` | `True` | Normalise weak rates using τ_n (neutron lifetime) |
| `numerical_precision` | `1e-7` | `solve_ivp` relative tolerance (rtol) for all ODE integration |

See {doc}`weak-rate-cache` for how changing these flags interacts with the
n↔p rate cache.
