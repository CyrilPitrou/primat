# Design: `generate_rates/thermal_average.ipynb`

Date: 2026-07-28

## Goal

Give a user a Jupyter notebook that turns an astrophysical S-factor — or a
bare cross-section — into a primat-format thermonuclear rate table, complete
with a Monte-Carlo-propagated 1σ uncertainty column and the detailed-balance
header line. It is the Python replacement for the Mathematica notebook
`PRIMAT-mma-dev/AlternateRates/Thermal-Average.nb` (whose repo-local stub is
`generate_rates/Thermal-Average.m`, containing only the `<<PRIMAT-Main.m`
loader).

The user should be able to drop in their own `S(E)` with an uncertainty and
get a usable rate file, without editing anything outside a single input cell.

## What the Mathematica original does

Extracted from the `.nb`'s input cells:

- `MassRed[n1,n2]` — reduced mass from PRIMAT's `Mass[]` nuclide data.
- `PhiMB[n1,n2,T,v] = √(2/π) (μc²/kT)^{3/2} exp(−μv²/2kT) v²/c³`
- `σ[n1,n2,Sfun,v] = S(E)/E · exp(−2πη)`, with `η = α_FS Z₁Z₂ c/v` and
  `E = ½μv²`.
- `Rate[T] = N_A ∫₀^{20 v_th} σ(v) Φ_MB(T,v) v dv`, evaluated with
  `NIntegrate` at `AccuracyGoal → PrecisionGoal → 40`.
- Tabulated over `ListTWagoner` (60 points, `PRIMAT-Main.m:2136`).
- Two worked S-factors, in MeV·barn with `E` in MeV:
  - `S_ddp(E) = 0.05520 + 0.2151 E − 0.02555 E²`
  - `S_ddn(E) = 0.05225 + 0.3655 E − 0.1799 E² + 0.05832 E³ − 0.007393 E⁴`

It has **no** uncertainty propagation and **no** file output. Both are added
here.

## Design decisions (settled with the user)

| Question | Decision |
|---|---|
| S(E) input form | Uniform interface accepting *either* a Python callable *or* a tabulated file |
| Uncertainty | Monte Carlo over the S(E) parameter vector |
| Temperature grid | Wagoner 60-point grid (as in the Mathematica original) |
| Code layout | Self-contained notebook — no new importable module |
| Output location | User-chosen `OUTDIR`, defaulting to an untracked overlay under `generate_rates/` |
| Cross-section entry points | `S(E)` for charged particles **and** direct `σ(E)` |
| Validation | Reproduce d+d as a worked, pre-filled example |

## Physics

The notebook integrates in energy rather than velocity. The two forms are
identical under `E = ½μv²`, but the energy form is better conditioned near the
Gamow peak:

```
σ(E)      = S(E)/E · exp(−2πη),      η = α_FS Z₁Z₂ √(μc²/2E)
N_A⟨σv⟩   = N_A (8/πμ)^{1/2} (kT)^{−3/2} ∫₀^∞ σ(E) E e^{−E/kT} dE
```

Units: `S` in MeV·barn, `σ` in barn, `E` in MeV, output `N_A⟨σv⟩` in
cm³ mol⁻¹ s⁻¹ — the unit of every shipped rate table.

Reference: Pitrou, Coc, Uzan & Vangioni, *Phys. Rept.* **04** (2018) 005
(`biblio/Pitrou_etal_PhysReptArxivVersion.pdf`), nuclear-rates section.

## Nuclide data provenance

Nothing is hard-coded. The chain is:

```
generate_rates/nubase_4.mas20.txt          NUBASE2020 evaluation
  ↓ generate_rates/nuclide_table.py        (offline, run once)
primat/data/csv/nuclides.csv               name,N,Z,A,Q,mass_excess_keV,spin
  ↓ primat/config.py  _load_nuclides()
PRIMATConfig.Nuclides / .NuclExcessMass / .NuclSpin
```

The notebook instantiates a live `PRIMATConfig()`, so it follows any
`data_dir` override and cannot drift from what the solver uses. Masses are
reconstructed exactly as `compute_detailed_balance_coefficients` does —
**nuclear**, not atomic:

```
M(s) = A·m_u + Δ(s)·keV − Z·m_e
```

with `m_u = cfg.ma`, `m_e = cfg.me`, and `Z` from `cfg.Nuclides[s] = [N, Z]`.
`α_FS` comes from `primat.constants.CONST.alphaem`; Avogadro's number is
derived as `1/m_u[g]` rather than typed in.

`nuclides.csv` only holds nuclides present in primat's reaction catalog. A
species outside it fails at `reaction_species(REACTION)`; the notebook catches
this and points the user at `nuclide_table.py`.

## Notebook structure

Single file, `generate_rates/thermal_average.ipynb`. Sections:

### §1 Physics preamble (markdown)
The formulas above, the unit conventions, the reference, and a one-paragraph
map of which cell the user edits.

### §2 Nuclide data and constants
Puts the repo root on `sys.path`, builds `cfg = PRIMATConfig()`, and defines
`nuclear_mass(name)`, `charge(name)`, `reduced_mass(n1, n2)`, `N_A`.

### §3 USER INPUT — the only cell to edit

```python
REACTION  = "d_d__t_p"     # primat reaction name
REF       = "MyFit2026"    # → d_d__t_p_MyFit2026.txt
MODE      = "S"            # "S" (MeV·barn) or "sigma" (barn)

def cross_section(E_MeV, theta):        # vectorized in E
    return theta[0] + theta[1]*E_MeV + theta[2]*E_MeV**2

THETA0    = np.array([0.05520, 0.2151, -0.02555])
COV       = None           # covariance, or None → error column ≡ 1.0
N_MC      = 300
OUTDIR    = REPO / "generate_rates/rate_tables_out/tables"
OVERWRITE = False
```

A user may instead assign `cross_section = from_table(path, kind=...)` (§4)
and/or supply `sample_theta(rng)` in place of `COV`.

### §4 Uniform interface adapter
Collapses every input form into one signature
`sigma_of_E(E_MeV, theta) → barn`:

- `MODE="S"` → divide by `E`, multiply by `exp(−2πη)`.
- `MODE="sigma"` → pass through, needed for neutron-induced reactions where
  `Z₁Z₂ = 0` and `S(E)` is not the natural variable.
- `from_table(path, kind)` → a `cross_section`-shaped callable that log-log
  interpolates a two- or three-column file, with `theta` acting as a
  log-normal overall normalisation so the Monte Carlo still applies.

No code downstream of this cell branches on the input form.

### §5 Thermal-average kernel
`Nsv(T9, theta)` implements the energy integral. Substituting `E = kT·x`, the
domain `x ∈ [0, 200]` reproduces the Mathematica cutoff `v < 20 v_th`. The
domain is split into log-spaced panels, each integrated with 20-point
Gauss–Legendre, so the Gamow peak is resolved at every temperature. The whole
evaluation is vectorized over `(60 temperatures × N_MC samples)`, making the
Monte Carlo one array operation.

### §6 Monte-Carlo uncertainty
Draws `theta` from `multivariate_normal(THETA0, COV)` (or the user's
`sample_theta`). Central value is `Nsv(T9, THETA0)`. The error column is
`sqrt(p84/p16)` of the sampled rate distribution at each temperature — the
multiplicative 1σ envelope. This matches both the convention documented in
the shipped `d_d__*_parthenope3.0.txt` headers and primat's own rate-variation
model, in which `p_<reaction>` samples the rate at `median · exp(p·expsigma)`.
`COV=None` yields a column of 1.0. Diagnostic plots show the S(E) band and
the resulting rate band.

### §7 Header and file write
`reaction_species(REACTION)` → `compute_detailed_balance_coefficients(...)`
produce the two header lines used by every shipped table:

```
# d + d > t + p   [d_d__t_p]   ref=MyFit2026
# detailed balance: alpha=1.73492 beta=0 gamma=-46.7971  Q=4.03266
# T9                 rate                error
```

Columns are written on the 60-point Wagoner T9 grid transcribed from
`PRIMAT-Main.m:2136`, `%.6e` formatted, to
`OUTDIR/<reaction>/<reaction>_<REF>.txt`. Existing files are not clobbered
unless `OVERWRITE=True`.

`OUTDIR` defaults to `generate_rates/rate_tables_out/tables`, which is
gitignored: the notebook never writes into the shipped `primat/data/` tree
unless the user deliberately redirects it there. The default layout is chosen
so that `rate_tables_out` works directly as a `user_nuclear_dir` overlay
(which is resolved per-file, so one overlaid table does not shadow the rest of
the shipped tree).

The Wagoner grid is coarse below T9 = 0.01 relative to primat's 1000-point
master grid, onto which `network_data.py` resamples at load time. This is the
same grid the Mathematica notebook used; the notebook states the trade-off in
a markdown note so a user who cares can widen it.

### §8 Worked example and validation
Pre-filled with `S_ddp` and `S_ddn` from the Mathematica notebook. A cell
overlays the generated rate on the shipped `d_d__t_p_primat.txt` and
`d_d__He3_n_primat.txt`, with a ratio panel — this is simultaneously the
user's template and the correctness check on §5. A closing cell shows how to
select the new table variant in a primat run.

### §9 Documentation pointer

`generate_rates/README.md`'s "Pipeline map" gains an entry for
`thermal_average.ipynb`, describing it as the user-facing entry point for
adding a rate from one's own S(E)/σ(E) — distinct from the bulk
`convert_ac2024_rates.py` regeneration — and noting that it writes a
per-reaction table variant into a `user_nuclear_dir`-shaped overlay rather
than rebuilding or modifying the shipped tree. The entry also records that it
supersedes the Mathematica `Thermal-Average.nb`.

## Guardrails

- Unknown `REACTION` → error listing near-matches from the reaction catalog,
  plus the `nuclide_table.py` hint when the species itself is missing.
- `MODE="S"` with `Z₁Z₂ = 0` → warning that `S(E)` is then merely `σE`.
- Non-positive-semidefinite `COV` → error before any integration runs.
- Negative sampled cross-sections clipped to zero, with the count reported so
  a badly-extrapolated polynomial does not fail silently.

## Out of scope

- No changes to `primat/` or `primat-c/`: the notebook only writes data files,
  so the backend-parity requirement in `CLAUDE.md` does not apply.
- No pytest regression, per the self-contained-notebook decision. §8's
  comparison against the shipped d+d tables is the validation.
- No electron screening, no resonance-by-resonance narrow-resonance
  formalism. A user needing either can encode it inside their own
  `cross_section` callable.
