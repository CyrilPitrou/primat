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
| `Omeganurel` | Ω_ν h² × 10⁶ **per flavour** (relativistic) |
| `OneOverOmeganunr` | 1 / (Ω_ν h² × 10⁻⁶) **per flavour** (non-relativistic) |

Both Ω_ν keys are **per neutrino flavour** (ν + ν̄), not summed over the three.
Multiply `Omeganurel` by 3 to compare with the usual quoted total (≈ 17 for
Neff ≈ 3.044). The per-flavour convention is the natural one for
`OneOverOmeganunr`, whose value ≈ 93 reproduces the standard
Σm_ν / 93.1 eV normalisation, and `Omeganurel` follows it for consistency.

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
t_s  a  T_gamma_MeV  T_nue_MeV  T_numu_MeV  T_nutau_MeV
Y_<nuclide> (one column per tracked species)
<reaction>_frwrd  (optional trailing block)
```

The authoritative contract for this schema — column names, order and
semantics — is `primat.evolution`'s module docstring; both backends' writers
conform to it, and `load_evolution` reads the header dynamically rather than
assuming a column count.

- The six leading columns are always present, in that order.
- `output_file` defaults to `results/output_tables.tsv` (relative to the
  current directory); set it to `None` to skip the disk write entirely — the
  time-evolution data is still accessible via `result["evolution"]` either
  way, on both backends (see `primat.backend`).
- Each `Y_<nuclide>` is an abundance per baryon, `Y_i = n_i/n_b`, not a mass
  fraction: the mass fraction is `A_i Y_i`, so `YPBBN` is `4 * Y_He4`.
- The `Y_<nuclide>` block is one column per nuclide of the chosen network —
  8 for `small`/`small_parthenope`, ~59 for `large`, fewer with an `amax`
  cutoff.
- The `<reaction>_frwrd` block (`output_rates_time_evolution=True`) is
  appended after the abundances: one column per reaction in the active LT
  network, lexicographically sorted, carrying the forward rate evaluated at
  each row's photon temperature. It is only available for
  `small`/`small_parthenope` — omitted (with a printed note) for
  `network="large"`.
- `a` and the three `T_nu` columns are `NaN` when the active background has
  no scale-factor / neutrino-sector tracking (e.g. a minimal custom
  background).
- `output_n_points` (default 500) controls how many interpolated rows the
  file has.

## Background time-evolution output

`output_background_evolution=True` writes a *separate* file — the background
thermodynamics rather than the nuclear network — to
`results/output_background.tsv`:

```text
T [MeV]  t [s]  a [1]  H [s^-1]  Tnue [MeV]  Tnumu [MeV]  Tnutau [MeV]
Nheating [1]  rho_plasma [MeV^4]  rho_nu_tot [MeV^4]  rho_extra [MeV^4]
rho_tot [MeV^4]
```

`Nheating` is the NEVO heating function (meaningful only with
`incomplete_decoupling=True`) and `rho_extra` the summed `extra_rho`
contribution (0 without one). This is the file to read for "what did the
expansion history actually do", as opposed to "what did the abundances do";
`runfiles/primat_run.py` turns both flags on.

`primat.evolution` (`load_evolution`, `dump_evolution`) and `primat.plotting`
provide tools for loading and plotting this data — see the
{doc}`../tutorials/AbundanceEvolution` notebook for a worked example.
