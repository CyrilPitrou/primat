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
| `"large"` | 427 + 1 = 428 | 59 | from the AC2024 compilation; LT era only |
| `"large"`, `amax=8` | 66 + 1 = 67 | 12 | the old "medium" network, less its double-counted ⁷Be+³He→2p+2α |
| `"large"`, `amax=2` | 2 + 1 = 3 | 3 | the old "deuterium" network's equivalent (n↔p + n_p__d_g + p_p_n__d_p) |
| any other name | — | — | loads `data/nuclear/networks/<name>.txt` |

The first number counts the *nuclear* reactions (what a network file lists,
or what survives the `amax` filter); the total is what
`load_network(...).n_reac` reports, since every network additionally carries
the n↔p weak reaction, which no file lists. Elsewhere in these docs `small`
is called "the 12-reaction network" and `large` "~428 reactions" — the same
two networks, counted the two different ways. Nuclide counts include `n` and
`p`.

Every network is integrated across the same three temperature eras — high
(HT), middle (MT) and low (LT), defined in the {doc}`../glossary`. HT evolves
n and p only, and MT intersects the chosen network with a fixed list of 18
reactions, the full network being too stiff to run there, so it holds 18
reactions for `large` and 13 for `small`. Only the LT reaction set is filtered
by `network`/`amax`. The light-element abundances of the
full `large` network match the `amax=8` restriction to ≲1e-4; its
heavy-nuclide tail (B, C, N, O, …) is approximate (limited by the AC2024
rate floors). See the {doc}`../tutorials/AbundanceEvolution` notebook for
evolution plots across networks.

## Two channels with the same reactants and products

Two reactions with the same reactants and products are one channel to the
solver: their rates add up. Products listed in a different order
(`Be7_He3__a_a_p_p` vs `Be7_He3__p_p_a_a`) or the reverse direction
(`d_p__He3_g` vs `He3_g__d_p`) still count as the same channel. Photons do
count, so `Li7_p__a_a` and `Li7_p__a_a_g` are different channels.

Loading a network that lists a channel twice prints one warning per run on
both backends, naming the reactions:

```text
warning: reactions d_p__He3_g, He3_g__d_p of network 'small' are the same channel (H2 + p -> He3 + g), so their rates add up. ...
```

To list a second channel on purpose (for example a resonant and a direct
contribution, each with its own table), add a `___<label>` suffix (three
underscores) to its name in the network file:

```text
d_p__He3_g,       d_p__He3_g_primat.txt
d_p__He3_g___dir, d_p__He3_g_direct.txt
```

The label is free text, and a third channel just takes another label
(`___dir2`). A labelled reaction is never reported as a duplicate. It
borrows from the reaction named before `___`:

- its reactants and products;
- its detailed-balance (reverse-rate) coefficients;
- its nuclear-QED correction;
- its table folder: `d_p__He3_g_direct.txt` above is read from
  `nuclear/tables/d_p__He3_g/`, the base reaction's folder.

So the second column must name a table file in that folder. The labelled
reaction gets its own `p_<name>`/`delta_<name>` rate-variation parameters,
e.g. `p_d_p__He3_g___dir`. The base reaction must be in the reaction catalog
(`csv/reactions_large.csv`) or be added in the same network.

The rate-table generator (`generate_rates/convert_ac2024_rates.py`) is
stricter: it refuses to write the shipped catalog if any channel appears
twice, since shipped data must have one source per channel.

## Related physics-configuration flags

| Flag | Default | Effect |
|------|---------|--------|
| `radiative_corrections` | `True` | Coulomb + T=0 resummed radiative corrections to n↔p (CCR) |
| `finite_mass_corrections` | `True` | Fokker-Planck finite-nucleon-mass correction (FM) |
| `thermal_corrections` | `True` | Finite-temperature radiative corrections to n↔p (CCRTh) |
| `spectral_distortions` | `True` | Correct n↔p rates for non-FD neutrino distributions (SD) |
| `tau_n_normalization` | `True` | Normalise weak rates using τ_n (neutron lifetime) |
| `numerical_precision` | `1e-7` | `solve_ivp` relative tolerance (rtol) for all ODE integration |

The four correction abbreviations (CCR, FM, CCRTh, SD) each have a line in the
{doc}`../glossary`. See {doc}`weak-rate-cache` for how changing these flags
interacts with the n↔p rate cache.
