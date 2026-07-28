# Reading the result dict and time-evolution output

## The result dict

`run_bbn()` (and `PRIMAT.solve()`) returns a dict with these keys:

| Key | Description |
|-----|-------------|
| `YPBBN` | Helium-4 mass fraction (BBN convention) |
| `YPCMB` | Helium-4 mass fraction (CMB convention) |
| `He4oH` | He4/H (by number) |
| `DoH` | D/H |
| `He3oH` | (He3+H3)/H |
| `Li7oH` | (Li7+Be7)/H |
| `Neff` | Effective number of neutrino species |
| `Omeganurel` | Ω_ν h² × 10⁶ (relativistic) |
| `OneOverOmeganunr` | 1 / (Ω_ν h² × 10⁻⁶) (non-relativistic) |

The neutrino-sector keys (`Neff`, `Omeganurel`, `OneOverOmeganunr`) are only
present if the background actually provides that information — see
`primat.main.PRIMAT.solve`.

## Monte-Carlo uncertainties

When a Monte-Carlo run is requested (`--mc N` on the CLI, or
`run_mc()`/`mc_uncertainty()` via `to_flat_dict()`), every observable above
also gets a matching `sigma_<key>` entry with its 1-sigma MC uncertainty,
e.g. `sigma_DoH` alongside `DoH`. See {doc}`rate-variation-mc` for how the
underlying samples are drawn.

## Time-evolution output

Set `output_time_evolution=True` to make the full cosmic-time evolution
available via the `"evolution"` key of the result dict (an
`EvolutionResult`, see `primat.evolution`). If `output_file` is also set to a
path, a TSV file is written with the unified schema:

```text
t_s, a, T_gamma_MeV, T_nue_MeV, T_numu_MeV, T_nutau_MeV, [Nheating],
Y_<nuclide> (one column per tracked species),
n_to_p_weak_rate, p_to_n_weak_rate, [nuclear rates]
```

- `output_file` defaults to `results/output_tables.tsv` (relative to the
  current directory); set it to `None` to skip the disk write entirely — the
  time-evolution data is still accessible via `result["evolution"]` either
  way, on both backends (see `primat.backend`).
- `Nheating` is included only for `incomplete_decoupling=True` (a real NEVO
  heating table).
- `[abundances]` is one `Y<species>` column per nuclide of the chosen
  network — 8 for `small`/`small_parthenope`, ~59 for `large`, fewer with an
  `amax` cutoff.
- `[nuclear rates]` (`output_rates_time_evolution=True`) is only available
  for `small`/`small_parthenope` — omitted (with a printed note) for
  `network="large"`.
- `output_n_points` (default 500) controls how many interpolated rows the
  file has.

`primat.evolution` (`load_evolution`, `dump_evolution`) and `primat.plotting`
provide tools for loading and plotting this data — see the
{doc}`../tutorials/AbundanceEvolution` notebook for a worked example.
