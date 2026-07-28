# Thermal-Average Rate-Generation Notebook — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task, inline in the current session. Do **not** dispatch subagents. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `generate_rates/thermal_average.ipynb`, a self-contained Jupyter notebook that turns a user-supplied astrophysical S-factor or cross-section — with a parameter covariance — into a primat-format rate table with a Monte-Carlo 1σ error column.

**Architecture:** One notebook, nine sections. All nuclide data (masses, charges, spins, Q-values, detailed-balance coefficients) is read from a live `PRIMATConfig()` and `primat.network_data`; nothing is hard-coded. The user edits exactly one cell. The thermal average is an energy-space integral evaluated with log-spaced Gauss–Legendre panels, vectorized over temperature and Monte-Carlo sample.

**Tech Stack:** Python 3, numpy, scipy (`roots_legendre`, `quad`), matplotlib, nbformat/nbclient (for building and executing the notebook), primat itself.

**Spec:** `docs/superpowers/specs/2026-07-28-thermal-average-notebook-design.md`

## Global Constraints

- The deliverable is a **single self-contained notebook**. Do not create an importable helper module under `generate_rates/` or `primat/`. No pytest files.
- **No changes to `primat/` or `primat-c/`** — source *or* data. The notebook reads primat APIs and writes only into its own untracked output directory `generate_rates/rate_tables_out/`; the shipped `primat/data/` tree is never modified, so `CLAUDE.md`'s backend-parity rule does not apply.
- **No hard-coded nuclide data.** Masses, charges, spins, Q-values and detailed-balance coefficients come from `PRIMATConfig()` / `primat.network_data`. `α_FS` comes from `primat.constants.CONST.alphaem`. Avogadro's number is derived as `1/m_u[g]`.
- **Comment heavily** (`CLAUDE.md`): every function gets a docstring stating what it computes and why, the meaning and units of each argument and of the return value, and a usage example. Every non-obvious step gets an inline comment. Every magic number is named and explained.
- Units throughout: energies in **MeV**, cross-sections in **barn** at the user interface and **cm²** inside the kernel, temperatures as **T9** (10⁹ K), output rate `N_A⟨σv⟩` in **cm³ mol⁻¹ s⁻¹**.
- The output temperature grid is the **60-point Wagoner grid** transcribed from `generate_rates/PRIMAT-Main.m:2136`.
- Output file format must byte-match the conventions of `generate_rates/convert_ac2024_rates.py:331` `write_reaction_file`: two `#` header lines then a `# T9 rate error` column line, `%.6e` values, `"   "` delimiter.
- Verify every task by **executing the notebook end to end** with the harness from Task 1. A task is not done until the notebook runs clean.

## File Structure

| File | Responsibility |
|---|---|
| `generate_rates/thermal_average.ipynb` | **Create.** The entire deliverable: nine sections, §3 being the only user-edited cell. |
| `generate_rates/README.md` | **Modify.** Add a "Pipeline map" entry pointing at the notebook. |
| `.gitignore` | **Modify.** Ignore `generate_rates/rate_tables_out/`. |
| `generate_rates/rate_tables_out/` | Generated, untracked. The notebook's default output: a `user_nuclear_dir`-shaped overlay holding `tables/<reaction>/<reaction>_<ref>.txt` and `networks/custom_rate.txt`. |
| `<scratchpad>/build_nb.py` | Throwaway. One-off nbformat script that creates the empty notebook in Task 1. Not committed. |
| `<scratchpad>/run_nb.py` | Throwaway. Executes the notebook and prints cell outputs; the verification harness for every task. Not committed. |

`<scratchpad>` means the session scratchpad directory. After Task 1 the notebook is edited with the `NotebookEdit` tool, never regenerated from a script.

---

### Task 1: Notebook scaffold, verification harness, and §1–§2 (nuclide data)

**Files:**
- Create: `generate_rates/thermal_average.ipynb`
- Create (throwaway): `<scratchpad>/build_nb.py`, `<scratchpad>/run_nb.py`

**Interfaces:**
- Consumes: nothing.
- Produces: notebook globals `cfg` (`PRIMATConfig`), `REPO` (`pathlib.Path`), `N_A` (float, mol⁻¹), `MEV_PER_GK` (float, MeV per unit T9), `BARN_CM2` (float), `ALPHA_FS` (float), and functions `nuclear_mass_MeV(name: str) -> float`, `charge(name: str) -> int`, `reduced_mass_MeV(n1: str, n2: str) -> float`.

- [ ] **Step 1: Write the verification harness**

Create `<scratchpad>/run_nb.py`:

```python
"""Execute generate_rates/thermal_average.ipynb and print every cell's output.

Usage:  python <scratchpad>/run_nb.py
Exit code 0 means the notebook ran to completion with no exception.
"""
import sys
import nbformat
from nbclient import NotebookClient

NB = "generate_rates/thermal_average.ipynb"
nb = nbformat.read(NB, as_version=4)
client = NotebookClient(nb, timeout=900, kernel_name="python3",
                        resources={"metadata": {"path": "generate_rates/"}})
try:
    client.execute()
except Exception as exc:                      # noqa: BLE001 - we want the traceback
    print("NOTEBOOK FAILED:", type(exc).__name__, exc)
    sys.exit(1)

for i, cell in enumerate(nb.cells):
    if cell.cell_type != "code":
        continue
    for out in cell.get("outputs", []):
        text = out.get("text") or out.get("data", {}).get("text/plain")
        if text:
            print(f"--- cell {i} ---")
            print(text if isinstance(text, str) else "".join(text))
print("NOTEBOOK OK")
```

- [ ] **Step 2: Create the empty notebook**

Create `<scratchpad>/build_nb.py` and run it from the repo root:

```python
import nbformat as nbf

nb = nbf.v4.new_notebook()
nb.cells = [nbf.v4.new_markdown_cell("# placeholder")]
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
}
nbf.write(nb, "generate_rates/thermal_average.ipynb")
print("created")
```

Run: `python <scratchpad>/build_nb.py`
Expected: `created`

- [ ] **Step 3: Write §1, the physics preamble (markdown cell)**

Replace the placeholder cell with a markdown cell containing exactly this content:

````markdown
# Thermal-averaged reaction rates from S(E) or σ(E)

This notebook converts a nuclear cross-section — supplied either as an
astrophysical S-factor `S(E)` or directly as `σ(E)` — into a **primat-format
thermonuclear rate table**, with a Monte-Carlo-propagated 1σ uncertainty
column and the correct detailed-balance header.

It is the Python replacement for the Mathematica notebook
`Thermal-Average.nb` (see `Thermal-Average.m`, the `<<PRIMAT-Main.m` stub).

## Physics

For charged particles the cross-section is factorised into the slowly-varying
astrophysical S-factor and the Coulomb (Gamow) penetration factor:

$$\sigma(E) = \frac{S(E)}{E}\,e^{-2\pi\eta},
\qquad \eta = \alpha_{\rm FS} Z_1 Z_2 \sqrt{\frac{\mu c^2}{2E}} $$

The thermal average over a Maxwell–Boltzmann distribution at temperature $T$ is

$$N_A\langle\sigma v\rangle
  = N_A \left(\frac{8}{\pi\mu}\right)^{1/2} (k_BT)^{-3/2}
    \int_0^\infty \sigma(E)\, E\, e^{-E/k_BT}\, {\rm d}E $$

This is the energy-space form of the velocity integral used in the Mathematica
notebook; the two are identical under $E = \tfrac12\mu v^2$, but the energy
form is better conditioned near the Gamow peak.

Reference: Pitrou, Coc, Uzan & Vangioni, *Physics Reports* **04** (2018) 005
(`biblio/Pitrou_etal_PhysReptArxivVersion.pdf`), nuclear-rates section.

## Units

| Quantity | Unit |
|---|---|
| Energy `E` | MeV |
| S-factor `S(E)` | MeV·barn |
| Cross-section `σ(E)` | barn |
| Temperature `T9` | 10⁹ K |
| Output `N_A⟨σv⟩` | cm³ mol⁻¹ s⁻¹ |

## How to use this notebook

**Edit §3 only.** Everything else is machinery. In §3 you declare the
reaction, a reference label for the output filename, whether you are giving
`S(E)` or `σ(E)`, the function itself, its parameter vector and covariance,
and where to write the result. Then run all cells.

All nuclide data — masses, charges, spins, Q-values, detailed-balance
coefficients — is read from primat itself, so this notebook cannot drift from
the solver's own nuclear data.
````

- [ ] **Step 4: Write §2, nuclide data and constants (markdown + code cell)**

Append a markdown cell:

```markdown
## §2 — Nuclide data and constants (from primat)

Nothing here is hard-coded. `PRIMATConfig` reads
`primat/data/csv/nuclides.csv` (generated offline by
`generate_rates/nuclide_table.py` from the NUBASE2020 evaluation
`nubase_4.mas20.txt`), giving `(N, Z)`, mass excess in keV, and spin for
every nuclide in primat's reaction catalog.
```

Then a code cell:

```python
import sys
from pathlib import Path

import numpy as np

# The notebook lives in generate_rates/; primat is importable from the repo root.
REPO = Path.cwd().parent if Path.cwd().name == "generate_rates" else Path.cwd()
sys.path.insert(0, str(REPO))

from primat.config import PRIMATConfig
from primat.constants import CONST

cfg = PRIMATConfig()

# --- Unit conversions and constants, all taken from primat ------------------
ALPHA_FS = CONST.alphaem              # fine-structure constant (dimensionless)
M_U_MEV  = cfg.ma                     # atomic mass unit          [MeV]
M_E_MEV  = cfg.me                     # electron mass             [MeV]
BARN_CM2 = 1.0e-24                    # 1 barn in cm^2 (definition of the barn)

# k_B * 1e9 K expressed in MeV, i.e. kT[MeV] = MEV_PER_GK * T9.  cfg.kB is in
# erg/K and cfg.MeV is 1 MeV in erg, so the ratio converts erg -> MeV.
MEV_PER_GK = cfg.kB * 1.0e9 / cfg.MeV

# Avogadro's number is *derived*, not typed in: it is the reciprocal of the
# atomic mass unit expressed in grams (m_u[MeV] * erg/MeV / c^2).
N_A = 1.0 / (M_U_MEV * cfg.MeV / cfg.clight**2)   # [mol^-1]


def nuclear_mass_MeV(name):
    """Rest-mass energy of a nuclide, in MeV.

    This is the *nuclear* mass (bare nucleus), not the atomic mass: the Z bound
    electrons are subtracted.  It is built exactly the way
    ``primat.network_data.compute_detailed_balance_coefficients`` builds it, so
    the reduced mass used here is consistent with the reverse-rate coefficients
    written into the output file's header:

        M = A * m_u + Delta - Z * m_e

    with ``Delta`` the mass excess from ``nuclides.csv`` (stored in keV, hence
    the 1e-3 conversion to MeV).

    Args:
        name: primat nuclide key, e.g. ``"n"``, ``"p"``, ``"H2"``, ``"He4"``.

    Returns:
        Rest-mass energy in MeV.

    Example:
        >>> round(nuclear_mass_MeV("H2"), 4)     # deuteron
        1875.6128
    """
    N, Z = cfg.Nuclides[name]
    A = N + Z
    return A * M_U_MEV + cfg.NuclExcessMass[name] * 1.0e-3 - Z * M_E_MEV


def charge(name):
    """Atomic number Z of a nuclide, from primat's ``Nuclides`` table.

    Z enters the Gamow penetration factor as the product Z1*Z2; it is zero for
    the neutron, which is why neutron-induced reactions have no Coulomb barrier.

    Args:
        name: primat nuclide key, e.g. ``"He4"``.

    Returns:
        Atomic number (int).

    Example:
        >>> charge("He4"), charge("n")
        (2, 0)
    """
    return cfg.Nuclides[name][1]


def reduced_mass_MeV(n1, n2):
    """Reduced rest-mass energy mu*c^2 of the entrance channel, in MeV.

    The thermal average is an integral over the *relative* kinetic energy of
    the two reactants, whose inertia is the reduced mass
    ``mu = m1 m2 / (m1 + m2)``.  Working with mu*c^2 in MeV keeps every energy
    in the notebook in the same unit.

    Args:
        n1, n2: primat nuclide keys of the two reactants.

    Returns:
        mu*c^2 in MeV.

    Example:
        >>> round(reduced_mass_MeV("H2", "H2"), 4)   # half the deuteron mass
        937.8064
    """
    m1, m2 = nuclear_mass_MeV(n1), nuclear_mass_MeV(n2)
    return m1 * m2 / (m1 + m2)


print(f"N_A        = {N_A:.6e} mol^-1        (expect 6.022141e+23)")
print(f"MEV_PER_GK = {MEV_PER_GK:.7f} MeV/T9  (expect 0.0861733)")
print(f"M(d)       = {nuclear_mass_MeV('H2'):.4f} MeV   (expect 1875.6128)")
print(f"mu(d,d)    = {reduced_mass_MeV('H2', 'H2'):.4f} MeV   (expect 937.8064)")
```

- [ ] **Step 5: Run the notebook and check the four printed values**

Run: `python <scratchpad>/run_nb.py`

Expected output includes:
```
N_A        = 6.022141e+23 mol^-1        (expect 6.022141e+23)
MEV_PER_GK = 0.0861733 MeV/T9  (expect 0.0861733)
M(d)       = 1875.6128 MeV   (expect 1875.6128)
mu(d,d)    = 937.8064 MeV   (expect 937.8064)
NOTEBOOK OK
```

If any value disagrees, the constant chain is wrong — fix before proceeding.

- [ ] **Step 6: Commit**

```bash
git add generate_rates/thermal_average.ipynb
git commit -m "generate_rates: scaffold thermal-average notebook with primat nuclide data"
```

---

### Task 2: §3 user-input cell and §4 uniform cross-section adapter

**Files:**
- Modify: `generate_rates/thermal_average.ipynb`

**Interfaces:**
- Consumes: `cfg`, `REPO`, `ALPHA_FS`, `BARN_CM2`, `charge`, `reduced_mass_MeV` (Task 1).
- Produces: user globals `REACTION`, `REF`, `MODE`, `cross_section`, `THETA0`, `COV`, `SAMPLE_THETA`, `N_MC`, `SEED`, `OUTDIR`, `OVERWRITE`; derived globals `REACTANTS: list[str]`, `PRODUCTS: list[str]`, `MU_MEV: float`, `Z1Z2: int`; functions `from_table(path, kind) -> callable`, `gamow_exponent(E_MeV, mu_MeV, Z1Z2) -> ndarray`, `sigma_of_E_cm2(E_MeV, theta) -> ndarray`, and the counter `NEGATIVE_SIGMA_COUNT: list[int]`.

- [ ] **Step 1: Append §3, the user-input markdown cell**

```markdown
## §3 — USER INPUT (edit this section, and nothing else)

`cross_section(E_MeV, theta)` must be vectorized in `E_MeV` and take a
1-D parameter vector `theta`, so that the Monte Carlo of §6 can resample it.
Return **MeV·barn** when `MODE="S"`, **barn** when `MODE="sigma"`.

Use `MODE="sigma"` for neutron-induced reactions (Z₁Z₂ = 0), where the
S-factor factorisation buys you nothing.

Instead of writing a formula you can read a tabulated cross-section from a
file with `cross_section = from_table(path, kind="S")` (see §4); `theta[0]`
then acts as a log-normal overall normalisation, so an overall systematic
uncertainty is set with `THETA0 = np.array([0.0])` and `COV = [[0.05**2]]`
for 5%.

Uncertainties: give `COV`, the covariance matrix of `theta`, for a Gaussian
Monte Carlo. Set `COV = None` for no uncertainty (the error column becomes
1.0). For a non-Gaussian prior, set `COV = None` and provide
`SAMPLE_THETA = lambda rng, n: <(n, len(THETA0)) array>` instead.
```

- [ ] **Step 2: Append the §3 code cell**

```python
# ===========================================================================
#  USER INPUT
# ===========================================================================
# The reaction, by primat's reaction name (the directory name under
# primat/data/nuclear/tables/).  Species short names are joined by "_", and
# reactants are separated from products by "__":  "d_d__t_p" is d + d -> t + p.
REACTION = "d_d__t_p"

# Reference label; becomes the filename suffix and appears in the file header.
REF = "Mathematica-Sddp"

# "S"     -> cross_section returns the astrophysical S-factor in MeV*barn
# "sigma" -> cross_section returns the cross-section directly, in barn
MODE = "S"


def cross_section(E_MeV, theta):
    """S-factor of d + d -> t + p, in MeV*barn.

    Worked example: the polynomial fit used in the Mathematica notebook
    ``Thermal-Average.nb``,
    ``S(E) = 0.05520 + 0.2151 E - 0.02555 E^2`` with E in MeV.  Replace this
    body with your own parametrisation.

    Args:
        E_MeV: centre-of-mass kinetic energy [MeV]; may be an array.
        theta: 1-D parameter vector; here the three polynomial coefficients.

    Returns:
        S(E) in MeV*barn, same shape as ``E_MeV``.

    Example:
        >>> float(cross_section(np.array([0.1]), THETA0))
        0.0765...
    """
    return theta[0] + theta[1] * E_MeV + theta[2] * E_MeV**2


# Central parameter values.
THETA0 = np.array([0.05520, 0.2151, -0.02555])

# Covariance matrix of theta, or None for "no uncertainty".  The Mathematica
# notebook quoted no uncertainty on this fit, so the worked example uses a
# 2% fully-correlated normalisation error as an illustration: a 2% error on
# the constant term alone would be inconsistent, so it is applied by scaling
# the whole vector (see SAMPLE_THETA below for the general case).
COV = np.diag((0.02 * np.abs(THETA0))**2)

# Optional non-Gaussian sampler: SAMPLE_THETA(rng, n) -> (n, len(THETA0)).
# Leave as None to use the Gaussian COV above.
SAMPLE_THETA = None

N_MC = 300          # Monte-Carlo samples; 300 is ample for a 16/84 percentile
SEED = 20260728     # fixed seed so the written table is reproducible

# Where to write.  The default is an untracked overlay directory, so this
# notebook never touches the shipped primat/data/ tree.  Its layout is exactly
# what primat's `user_nuclear_dir` parameter expects, so a run can pick the new
# table up with  user_nuclear_dir="<REPO>/generate_rates/rate_tables_out"
# (see §8).  Point OUTDIR at primat/data/nuclear/tables only if you really do
# intend to modify the shipped data.
OUTDIR = REPO / "generate_rates" / "rate_tables_out" / "tables"

OVERWRITE = False   # refuse to clobber an existing file unless True
# ===========================================================================
```

- [ ] **Step 3: Append §4, the adapter markdown cell**

```markdown
## §4 — Uniform cross-section interface

Everything the user can supply — an S-factor formula, a σ(E) formula, or a
tabulated file — is collapsed here into a single function
`sigma_of_E_cm2(E_MeV, theta)` returning cm². No code below this point
branches on `MODE`.
```

- [ ] **Step 4: Append the §4 code cell**

```python
from primat.network_data import reaction_species

# --- Resolve the reaction against primat's catalog --------------------------
try:
    REACTANTS, PRODUCTS = reaction_species(REACTION)
except Exception as exc:
    # A bad reaction name is by far the most common user error, so say exactly
    # what went wrong and offer the near-misses from the shipped tables tree.
    known = sorted(p.name for p in (REPO / "primat/data/nuclear/tables").iterdir()
                   if p.is_dir())
    stem = REACTION.split("__")[0]
    close = [k for k in known if k.startswith(stem[:3])]
    raise ValueError(
        f"REACTION={REACTION!r} is not a reaction primat knows: {exc}\n"
        f"Did you mean one of: {close[:10]}\n"
        "If the *nuclide* itself is missing from primat/data/csv/nuclides.csv, "
        "extend the catalog with generate_rates/nuclide_table.py first."
    ) from exc

if len(REACTANTS) != 2:
    raise ValueError(
        f"{REACTION} has {len(REACTANTS)} reactants ({REACTANTS}); the "
        "two-body thermal average implemented here needs exactly 2."
    )

MU_MEV = reduced_mass_MeV(*REACTANTS)          # entrance-channel reduced mass
Z1Z2 = charge(REACTANTS[0]) * charge(REACTANTS[1])

if MODE == "S" and Z1Z2 == 0:
    print(f"WARNING: {REACTION} has Z1*Z2 = 0, so there is no Coulomb barrier "
          "and the Gamow factor is 1.  S(E) then means nothing more than "
          "sigma(E)*E; MODE='sigma' is the natural choice.")

# Counter for clipped unphysical (negative) cross-sections; a badly
# extrapolated polynomial must not fail silently.
NEGATIVE_SIGMA_COUNT = [0]


def gamow_exponent(E_MeV, mu_MeV, z1z2):
    """The exponent 2*pi*eta of the Coulomb penetration factor.

    The Sommerfeld parameter is eta = Z1 Z2 alpha c / v; with the relative
    velocity written in terms of the centre-of-mass energy,
    v = sqrt(2E/mu), this becomes eta = Z1 Z2 alpha sqrt(mu c^2 / 2E).  The
    tunnelling probability through the Coulomb barrier is exp(-2 pi eta),
    which is what suppresses charged-particle rates at low temperature.

    Args:
        E_MeV: centre-of-mass energy [MeV], array-like.
        mu_MeV: reduced rest-mass energy mu*c^2 [MeV].
        z1z2: product of the reactants' atomic numbers (0 for neutrons).

    Returns:
        2*pi*eta, dimensionless, same shape as ``E_MeV``.

    Example:
        >>> float(gamow_exponent(np.array([0.1]), 937.8064, 1))
        9.93...
    """
    return 2.0 * np.pi * ALPHA_FS * z1z2 * np.sqrt(mu_MeV / (2.0 * E_MeV))


def from_table(path, kind="S"):
    """Build a ``cross_section``-shaped callable from a tabulated file.

    Reads a whitespace- or comma-separated file whose first column is the
    energy in MeV and whose second column is S(E) in MeV*barn (``kind="S"``)
    or sigma(E) in barn (``kind="sigma"``); further columns are ignored.
    Interpolation is linear in log-log, which is the right choice for
    quantities spanning decades, and is clamped (not extrapolated) outside
    the tabulated range — extrapolating a measured cross-section is a
    physics decision the user should make deliberately, not a side effect.

    ``theta[0]`` acts as a log-normal overall normalisation
    (sigma -> sigma * exp(theta[0])), so the Monte Carlo of §6 propagates an
    overall systematic uncertainty with ``THETA0 = np.array([0.0])`` and
    ``COV = [[rel_err**2]]``.

    Args:
        path: path to the two-or-more-column table.
        kind: ``"S"`` or ``"sigma"``, matching ``MODE``.

    Returns:
        A callable ``f(E_MeV, theta)`` with the same contract as the
        hand-written ``cross_section``.

    Example:
        >>> cross_section = from_table("my_sfactor.dat", kind="S")
        >>> THETA0 = np.array([0.0]); COV = np.array([[0.05**2]])
    """
    data = np.loadtxt(path, delimiter=None, comments="#", ndmin=2)
    E_tab, y_tab = data[:, 0], data[:, 1]
    order = np.argsort(E_tab)
    logE, logy = np.log(E_tab[order]), np.log(y_tab[order])

    def f(E_MeV, theta):
        # np.interp clamps to the end values outside [E_tab[0], E_tab[-1]].
        y = np.exp(np.interp(np.log(E_MeV), logE, logy))
        return y * np.exp(theta[0])

    f.kind = kind
    return f


def sigma_of_E_cm2(E_MeV, theta):
    """Cross-section in cm^2, whatever form the user supplied it in.

    This is the single interface the thermal-average kernel of §5 sees.  For
    ``MODE="S"`` it applies the Gamow factorisation
    ``sigma = S(E)/E * exp(-2 pi eta)``; for ``MODE="sigma"`` it merely
    converts barn to cm^2.  Negative values (an over-extrapolated polynomial,
    say) are unphysical and are clipped to zero, with a running count kept in
    ``NEGATIVE_SIGMA_COUNT`` so the clipping is reported rather than hidden.

    Args:
        E_MeV: centre-of-mass energy [MeV], array-like.
        theta: parameter vector passed straight through to ``cross_section``.

    Returns:
        sigma(E) in cm^2, same shape as ``E_MeV``.

    Example:
        >>> float(sigma_of_E_cm2(np.array([0.1]), THETA0))    # d+d at 100 keV
        3.7e-26...
    """
    y = np.asarray(cross_section(E_MeV, theta), dtype=float)
    if MODE == "S":
        sigma_barn = y / E_MeV * np.exp(-gamow_exponent(E_MeV, MU_MEV, Z1Z2))
    elif MODE == "sigma":
        sigma_barn = y
    else:
        raise ValueError(f"MODE must be 'S' or 'sigma', got {MODE!r}")
    n_neg = int(np.count_nonzero(sigma_barn < 0.0))
    if n_neg:
        NEGATIVE_SIGMA_COUNT[0] += n_neg
        sigma_barn = np.clip(sigma_barn, 0.0, None)
    return sigma_barn * BARN_CM2


print(f"{REACTION}: {' + '.join(REACTANTS)} -> {' + '.join(PRODUCTS)}")
print(f"  mu = {MU_MEV:.4f} MeV, Z1*Z2 = {Z1Z2}")
print(f"  2*pi*eta at E = 100 keV: "
      f"{float(gamow_exponent(np.array([0.1]), MU_MEV, Z1Z2)):.4f}")
print(f"  sigma(100 keV) = {float(sigma_of_E_cm2(np.array([0.1]), THETA0)):.4e} cm^2")
```

- [ ] **Step 5: Run and check the Gamow exponent by hand**

Run: `python <scratchpad>/run_nb.py`

Expected: `NOTEBOOK OK`, and the printed values must satisfy the hand check

```
2*pi*eta = 2*pi*(1/137.036)*1*sqrt(937.8064/(2*0.1)) = 9.9314
```

so the line must read `2*pi*eta at E = 100 keV: 9.9314` (±0.0002 — `CONST.alphaem` may differ from 1/137.036 in the last digits). And

```
sigma(100 keV) = S(0.1)/0.1 * exp(-9.9314) * 1e-24
               = 0.0765155/0.1 * 4.8552e-5 * 1e-24 = 3.7150e-26 cm^2
```

so `sigma(100 keV) = 3.7150e-26 cm^2` (±1 in the last digit).

Also confirm the header line reads `d_d__t_p: H2 + H2 -> H3 + p` and
`mu = 937.8064 MeV, Z1*Z2 = 1`.

- [ ] **Step 6: Verify the neutron-channel warning fires**

Temporarily set `REACTION = "n_p__d_g"` and `MODE = "S"` in §3, run
`python <scratchpad>/run_nb.py`, and confirm the output contains
`WARNING: n_p__d_g has Z1*Z2 = 0`. Then restore `REACTION = "d_d__t_p"` and
re-run to confirm `NOTEBOOK OK`.

- [ ] **Step 7: Verify the bad-reaction error message**

Temporarily set `REACTION = "d_d__t_q"`, run `python <scratchpad>/run_nb.py`,
and confirm it exits 1 with a message containing
`is not a reaction primat knows` and `Did you mean one of:`. Restore
`REACTION = "d_d__t_p"` and re-run to confirm `NOTEBOOK OK`.

- [ ] **Step 8: Commit**

```bash
git add generate_rates/thermal_average.ipynb
git commit -m "generate_rates: add user-input cell and uniform S(E)/sigma(E) adapter"
```

---

### Task 3: §5 thermal-average kernel

**Files:**
- Modify: `generate_rates/thermal_average.ipynb`

**Interfaces:**
- Consumes: `sigma_of_E_cm2`, `MU_MEV`, `Z1Z2`, `MEV_PER_GK`, `N_A`, `cfg` (Tasks 1–2).
- Produces: `T9_WAGONER: ndarray` (60 values), `X_NODES: ndarray`, `X_WEIGHTS: ndarray`, and `thermal_average(T9, theta) -> ndarray` returning `N_A⟨σv⟩` in cm³ mol⁻¹ s⁻¹.

- [ ] **Step 1: Append the §5 markdown cell**

````markdown
## §5 — The thermal average

Substituting `x = E/kT` turns the integral into

$$N_A\langle\sigma v\rangle
  = N_A\, c \sqrt{\frac{8}{\pi\,\mu c^2}}\;\sqrt{k_BT}
    \int_0^\infty \sigma(k_BT\,x)\; x\, e^{-x}\,{\rm d}x $$

which is dimensionally transparent: `σ` in cm², the square roots in MeV^-1/2
and MeV^1/2, and `c` in cm/s, giving cm³ s⁻¹.

The integrand is sharply peaked for charged particles — the Gamow peak sits at
`x₀ = (b/2)^{2/3} (k_BT)^{-1/3}` with `b = 2π α Z₁Z₂ √(μc²/2)` — so the domain
is split into **logarithmically spaced panels**, each integrated with a
20-point Gauss–Legendre rule. This resolves the peak at every temperature,
which a single uniform rule would not.

The upper cutoff `x_max` defaults to 200, reproducing the Mathematica
notebook's `v < 20 v_th`, and is widened automatically if the Gamow peak of a
highly-charged channel would otherwise sit near the edge.

**Temperature grid.** The output uses the 60-point Wagoner grid transcribed
from `PRIMAT-Main.m:2136`, matching the Mathematica notebook. Note this grid
is coarse below T9 = 0.01 compared with primat's 1000-point master grid, onto
which `network_data.py` resamples at load time; widen `T9_WAGONER` if your
reaction matters down there.
````

- [ ] **Step 2: Append the §5 code cell**

```python
from scipy.special import roots_legendre

# The 60-point Wagoner temperature grid, transcribed verbatim from
# generate_rates/PRIMAT-Main.m:2136 (where it is given in K, i.e. x 1e9).
T9_WAGONER = np.array([
    0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01,
    0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.018, 0.02, 0.025, 0.03,
    0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13,
    0.14, 0.15, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.25, 1.5, 1.75, 2.,
    2.5, 3., 3.5, 4., 5., 6., 7., 8., 9., 10.,
])
assert len(T9_WAGONER) == 60, len(T9_WAGONER)

# --- Quadrature nodes in x = E/kT -------------------------------------------
# X_MIN is not 0: the integrand vanishes there (as x e^{-x}, and far faster
# still for charged particles), and a strictly positive lower limit lets the
# panels be spaced logarithmically.  1e-8 is ~8 orders of magnitude below the
# thermal peak at x ~ 1, so nothing measurable is lost.
X_MIN = 1.0e-8
_X_MAX_DEFAULT = 200.0     # matches the Mathematica cutoff v < 20 v_th
_N_PANELS = 24             # log-spaced panels spanning [X_MIN, x_max]
_N_GL = 20                 # Gauss-Legendre nodes per panel


def _build_x_quadrature(T9_grid):
    """Build the Gauss-Legendre nodes and weights in x = E/kT.

    A single rule cannot serve both the thermal peak (x ~ 1) and the Gamow
    peak (x0 can reach several tens at the lowest temperatures), so the domain
    is cut into logarithmically spaced panels, each carrying its own 20-point
    Gauss-Legendre rule.  The upper cutoff is enlarged beyond the Mathematica
    default of 200 if a highly-charged channel would put the Gamow peak within
    a few peak-widths of it.

    Args:
        T9_grid: temperatures [1e9 K] the quadrature must serve; only the
            smallest matters, since x0 grows as (kT)^(-1/3).

    Returns:
        ``(nodes, weights)``, both 1-D arrays of length _N_PANELS * _N_GL,
        such that ``sum(weights * f(nodes))`` approximates the x-integral.

    Example:
        >>> nodes, weights = _build_x_quadrature(np.array([0.001, 10.0]))
        >>> abs(np.sum(weights * nodes * np.exp(-nodes)) - 1.0) < 1e-12
        True
    """
    x_max = _X_MAX_DEFAULT
    if Z1Z2 > 0:
        # Gamow peak position, standard saddle-point result:
        #   E0 = (b kT / 2)^{2/3},  b = 2 pi alpha Z1 Z2 sqrt(mu c^2 / 2)
        # so in the scaled variable  x0 = E0/kT = (b/2)^{2/3} (kT)^{-1/3}.
        b = 2.0 * np.pi * ALPHA_FS * Z1Z2 * np.sqrt(MU_MEV / 2.0)
        kT_min = MEV_PER_GK * np.min(T9_grid)
        x0 = (b / 2.0) ** (2.0 / 3.0) * kT_min ** (-1.0 / 3.0)
        # Keep at least a factor 3 of headroom above the peak.
        x_max = max(x_max, 3.0 * x0)

    edges = np.geomspace(X_MIN, x_max, _N_PANELS + 1)
    gl_x, gl_w = roots_legendre(_N_GL)          # nodes/weights on [-1, 1]
    nodes, weights = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        half = 0.5 * (hi - lo)
        nodes.append(0.5 * (hi + lo) + half * gl_x)
        weights.append(half * gl_w)
    return np.concatenate(nodes), np.concatenate(weights)


X_NODES, X_WEIGHTS = _build_x_quadrature(T9_WAGONER)


def thermal_average(T9, theta):
    """Maxwell-Boltzmann-averaged reaction rate N_A <sigma v>.

    Evaluates

        N_A <sigma v> = N_A c sqrt(8 / (pi mu c^2)) sqrt(kT)
                        * Integral_0^inf sigma(kT x) x exp(-x) dx

    the energy-space form of the velocity integral in the Mathematica notebook
    ``Thermal-Average.nb``.  The x-integral is done with the log-panelled
    Gauss-Legendre rule built above, and the whole thing is vectorized over
    temperature, so a full 60-point table costs one array evaluation of the
    user's cross-section.

    Args:
        T9: temperature(s) in 1e9 K; scalar or 1-D array.
        theta: parameter vector forwarded to the user's ``cross_section``.

    Returns:
        N_A <sigma v> in cm^3 mol^-1 s^-1, shaped like ``T9``.

    Example:
        >>> float(thermal_average(1.0, THETA0))    # d+d -> t+p at T9 = 1
        4.0e+05...
    """
    T9 = np.atleast_1d(np.asarray(T9, dtype=float))
    kT = MEV_PER_GK * T9                                  # [MeV]

    # sigma is evaluated on the outer product (temperature x quadrature node).
    E = kT[:, None] * X_NODES[None, :]                    # [MeV]
    integrand = sigma_of_E_cm2(E, theta) * X_NODES[None, :] * np.exp(-X_NODES)[None, :]
    integral = integrand @ X_WEIGHTS                      # [cm^2]

    prefactor = cfg.clight * np.sqrt(8.0 / (np.pi * MU_MEV)) * np.sqrt(kT)
    return N_A * prefactor * integral


print("thermal_average(T9=1) =", f"{float(thermal_average(1.0, THETA0)):.6e}",
      "cm^3 mol^-1 s^-1")
```

- [ ] **Step 3: Append the §5 self-test markdown cell**

```markdown
### §5a — Quadrature self-tests

Three checks, run every time the notebook executes:

1. **Constant σ.** For σ(E) = σ₀ the average is exactly the mean
   Maxwell–Boltzmann speed, `⟨σv⟩ = σ₀ √(8kT/πμ)`.
2. **1/v law.** For σ(E) = σ₀√(E₀/E) the rate is exactly
   `σ₀ c √(2E₀/μc²)`, independent of temperature — the reason thermal
   neutron-capture rates are quoted as constants.
3. **Peaked Gamow integrand.** The panelled rule is checked against an
   adaptive `scipy.integrate.quad` on the same integrand, at four
   temperatures spanning the grid.
```

- [ ] **Step 4: Append the §5 self-test code cell**

```python
from scipy.integrate import quad

_saved = (MODE, cross_section)          # the tests swap the user's input out

# --- Test 1: constant cross-section -> mean MB speed ------------------------
MODE = "sigma"
cross_section = lambda E, th: np.full_like(np.asarray(E, dtype=float), th[0])
_sigma0_barn = 1.0
_T9 = np.array([0.001, 0.01, 0.1, 1.0, 10.0])
_got = thermal_average(_T9, np.array([_sigma0_barn]))
_kT = MEV_PER_GK * _T9
_want = N_A * _sigma0_barn * BARN_CM2 * cfg.clight * np.sqrt(8.0 * _kT / (np.pi * MU_MEV))
_err1 = np.max(np.abs(_got / _want - 1.0))
print(f"test 1 (constant sigma):  max rel. error = {_err1:.3e}")
assert _err1 < 1e-12, _err1

# --- Test 2: 1/v law -> temperature-independent rate ------------------------
_E0 = 1.0e-8                            # reference energy [MeV] (arbitrary)
cross_section = lambda E, th: th[0] * np.sqrt(_E0 / np.asarray(E, dtype=float))
_got = thermal_average(_T9, np.array([_sigma0_barn]))
_want = (N_A * _sigma0_barn * BARN_CM2 * cfg.clight
         * np.sqrt(2.0 * _E0 / MU_MEV))
_err2 = np.max(np.abs(_got / _want - 1.0))
print(f"test 2 (1/v law):         max rel. error = {_err2:.3e}")
assert _err2 < 1e-10, _err2

# --- Test 3: peaked Gamow integrand vs adaptive quadrature ------------------
MODE, cross_section = _saved
_errs = []
for _t9 in (0.01, 0.1, 1.0, 10.0):
    _kt = MEV_PER_GK * _t9
    # Integrate in u = log x, where the Gamow-peaked integrand is smooth and
    # quad's adaptive subdivision has no trouble finding it.
    _f = lambda u: float(sigma_of_E_cm2(np.array([_kt * np.exp(u)]), THETA0)[0]
                         * np.exp(u) * np.exp(-np.exp(u)) * np.exp(u))
    _ref_int, _ = quad(_f, np.log(X_MIN), np.log(X_NODES[-1]), limit=400)
    _ref = (N_A * cfg.clight * np.sqrt(8.0 / (np.pi * MU_MEV)) * np.sqrt(_kt)
            * _ref_int)
    _errs.append(abs(float(thermal_average(_t9, THETA0)) / _ref - 1.0))
print(f"test 3 (Gamow vs quad):   max rel. error = {max(_errs):.3e}")
assert max(_errs) < 1e-8, _errs

print("all quadrature self-tests passed")
```

- [ ] **Step 5: Run and verify all three tests pass**

Run: `python <scratchpad>/run_nb.py`

Expected:
```
test 1 (constant sigma):  max rel. error = <below 1e-12>
test 2 (1/v law):         max rel. error = <below 1e-10>
test 3 (Gamow vs quad):   max rel. error = <below 1e-8>
all quadrature self-tests passed
NOTEBOOK OK
```

If test 3 fails, raise `_N_PANELS` (24 → 40) rather than loosening the
assertion; the panelled rule, not the reference, is the thing under test.

- [ ] **Step 6: Commit**

```bash
git add generate_rates/thermal_average.ipynb
git commit -m "generate_rates: add thermal-average kernel with quadrature self-tests"
```

---

### Task 4: §6 Monte-Carlo uncertainty

**Files:**
- Modify: `generate_rates/thermal_average.ipynb`

**Interfaces:**
- Consumes: `thermal_average`, `T9_WAGONER`, `THETA0`, `COV`, `SAMPLE_THETA`, `N_MC`, `SEED`, `NEGATIVE_SIGMA_COUNT` (Tasks 2–3).
- Produces: `RATE_CENTRAL: ndarray` (60,), `RATE_ERROR: ndarray` (60,), `RATE_SAMPLES: ndarray | None` of shape `(N_MC, 60)`.

- [ ] **Step 1: Append the §6 markdown cell**

```markdown
## §6 — Monte-Carlo uncertainty

The central rate is `thermal_average(T9, THETA0)`. The uncertainty is
obtained by resampling `theta` and re-integrating: `N_MC` draws from
`N(THETA0, COV)` (or from the user's `SAMPLE_THETA`), each giving a full
60-point rate curve.

The error column written to the table is the **multiplicative 1σ envelope**

```
error(T9) = sqrt( p84(T9) / p16(T9) )
```

with p16/p84 the 15.865th and 84.135th percentiles of the sampled rate. This
is the convention documented in the shipped
`d_d__*_parthenope3.0.txt` headers, and it is what primat's own rate-variation
model expects: setting `p_<reaction> = p` samples the rate at
`median * error**p`.

With `COV = None` and `SAMPLE_THETA = None` the column is filled with 1.0.
```

- [ ] **Step 2: Append the §6 code cell**

```python
RATE_CENTRAL = thermal_average(T9_WAGONER, THETA0)

if COV is None and SAMPLE_THETA is None:
    # No uncertainty requested: a multiplicative envelope of exactly 1.
    RATE_SAMPLES = None
    RATE_ERROR = np.ones_like(RATE_CENTRAL)
    print("no uncertainty requested (COV and SAMPLE_THETA are both None); "
          "error column set to 1.0")
else:
    rng = np.random.default_rng(SEED)
    if SAMPLE_THETA is not None:
        thetas = np.asarray(SAMPLE_THETA(rng, N_MC), dtype=float)
        if thetas.shape != (N_MC, len(THETA0)):
            raise ValueError(
                f"SAMPLE_THETA returned shape {thetas.shape}, "
                f"expected {(N_MC, len(THETA0))}"
            )
    else:
        cov = np.atleast_2d(np.asarray(COV, dtype=float))
        if cov.shape != (len(THETA0), len(THETA0)):
            raise ValueError(
                f"COV has shape {cov.shape}, expected "
                f"{(len(THETA0), len(THETA0))} to match THETA0"
            )
        try:
            # Cholesky is both the sampler and the positive-definiteness test:
            # it raises before any expensive integration if COV is not PSD.
            chol = np.linalg.cholesky(cov)
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                "COV is not positive definite, so it cannot be a covariance "
                f"matrix: {exc}.  Eigenvalues: {np.linalg.eigvalsh(cov)}"
            ) from exc
        thetas = THETA0 + rng.standard_normal((N_MC, len(THETA0))) @ chol.T

    RATE_SAMPLES = np.array([thermal_average(T9_WAGONER, th) for th in thetas])

    p16, p84 = np.percentile(RATE_SAMPLES, [15.865, 84.135], axis=0)
    if np.any(p16 <= 0.0):
        bad = T9_WAGONER[p16 <= 0.0]
        raise ValueError(
            "the 16th-percentile rate is non-positive at T9 = "
            f"{bad}, so a multiplicative error factor is undefined.  Your "
            "parameter uncertainties are large enough to drive the "
            "cross-section to zero there; narrow COV or use SAMPLE_THETA "
            "with a physically bounded prior."
        )
    RATE_ERROR = np.sqrt(p84 / p16)
    print(f"Monte Carlo: {N_MC} samples, "
          f"error factor ranges over [{RATE_ERROR.min():.4f}, "
          f"{RATE_ERROR.max():.4f}]")

if NEGATIVE_SIGMA_COUNT[0]:
    print(f"NOTE: clipped {NEGATIVE_SIGMA_COUNT[0]} negative cross-section "
          "evaluations to zero.  This means your parametrisation goes "
          "unphysical somewhere in the sampled energy/parameter range — check "
          "the S(E) plot in §8 before trusting the result.")
```

- [ ] **Step 3: Run and check the Monte Carlo**

Run: `python <scratchpad>/run_nb.py`

Expected: `NOTEBOOK OK`, and a line of the form
`Monte Carlo: 300 samples, error factor ranges over [1.0…, 1.0…]`. Because
the worked example applies a 2% uncertainty to each coefficient, every entry
of `RATE_ERROR` must lie in `(1.0, 1.1)` — a 2% input cannot produce a
double-digit-percent output envelope. If it does, the sampling is wrong.

- [ ] **Step 4: Verify the non-PSD covariance guardrail**

Temporarily set `COV = np.array([[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]])`
in §3 (symmetric but indefinite), run `python <scratchpad>/run_nb.py`, and
confirm it exits 1 with `COV is not positive definite`. Restore the original
`COV` and re-run to confirm `NOTEBOOK OK`.

- [ ] **Step 5: Commit**

```bash
git add generate_rates/thermal_average.ipynb
git commit -m "generate_rates: add Monte-Carlo uncertainty propagation"
```

---

### Task 5: §7 header construction and file writing

**Files:**
- Modify: `generate_rates/thermal_average.ipynb`

**Interfaces:**
- Consumes: `REACTION`, `REF`, `REACTANTS`, `PRODUCTS`, `OUTDIR`, `OVERWRITE`, `T9_WAGONER`, `RATE_CENTRAL`, `RATE_ERROR`, `cfg`, `MEV_PER_GK` (Tasks 1–4).
- Produces: `OUT_PATH: pathlib.Path`, `ALPHA_DB`, `BETA_DB`, `GAMMA_DB`, `Q_MEV` (floats), and `write_rate_table() -> Path`.

- [ ] **Step 1: Append the §7 markdown cell**

````markdown
## §7 — Detailed-balance header and file output

primat stores only the *forward* rate; the reverse rate is reconstructed as

$$ {\rm backward}(T_9) = \alpha\, T_9^{\beta}\, e^{\gamma/T_9}\;{\rm forward}(T_9) $$

with α, β, γ derived from nuclide masses and spins alone by
`primat.network_data.compute_detailed_balance_coefficients`. Reusing that
function — rather than copying numbers out of an existing header — is what
guarantees the new table's reverse rate is consistent with the rest of the
network.

The file is written to `OUTDIR/<reaction>/<reaction>_<REF>.txt` — the same
per-reaction-folder layout the shipped tree uses, which is the mechanism
primat uses for multiple candidate tables per reaction
(`network_data.available_rate_tables()`).

`OUTDIR` defaults to the untracked `generate_rates/rate_tables_out/tables`, so
**nothing is written into `primat/data/`**. That directory doubles as a
ready-made `user_nuclear_dir` overlay; §8 shows how to point a run at it.
Overlay resolution is per-file, so an overlaid `d_d__t_p` does not shadow any
other table.
````

- [ ] **Step 2: Append the §7 code cell**

```python
from primat.network_data import compute_detailed_balance_coefficients

ALPHA_DB, BETA_DB, GAMMA_DB = compute_detailed_balance_coefficients(
    REACTANTS, PRODUCTS, cfg)

# gamma is defined as -Q/(k_B * 1e9 K), so inverting it recovers Q in MeV.
# Deriving Q this way rather than re-summing binding energies guarantees the
# header's Q and gamma can never disagree.
Q_MEV = -GAMMA_DB * MEV_PER_GK


def write_rate_table():
    """Write the computed rate to a primat-format table file.

    Reproduces byte-for-byte the layout of the shipped tables (see
    ``generate_rates/convert_ac2024_rates.py``'s ``write_reaction_file``): two
    ``#`` header lines carrying the reaction, the reference label and the
    detailed-balance coefficients, then a ``#`` column-name line, then three
    ``%.6e`` columns ``T9 / rate / error``.

    The reaction display string is derived from ``REACTION`` itself: primat's
    reaction names join species short names with ``_`` and separate reactants
    from products with ``__``, and no species short name contains ``_``, so
    the split is unambiguous.

    Returns:
        The ``pathlib.Path`` written.

    Raises:
        FileExistsError: if the target exists and ``OVERWRITE`` is False.

    Example:
        >>> write_rate_table()
        PosixPath('.../generate_rates/rate_tables_out/tables/d_d__t_p/d_d__t_p_Mathematica-Sddp.txt')
    """
    lhs, rhs = REACTION.split("__")
    display = f"{' + '.join(lhs.split('_'))} > {' + '.join(rhs.split('_'))}"

    header = (
        f"{display}   [{REACTION}]   ref={REF}\n"
        f"detailed balance: alpha={ALPHA_DB:.6g} beta={BETA_DB:.6g} "
        f"gamma={GAMMA_DB:.6g}  Q={Q_MEV:.6g}\n"
        f"T9                 rate                error"
    )

    outdir = Path(OUTDIR) / REACTION
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{REACTION}_{REF}.txt"
    if path.exists() and not OVERWRITE:
        raise FileExistsError(
            f"{path} already exists.  Set OVERWRITE = True in §3 to replace "
            "it, or change REF to write a new variant alongside it."
        )

    np.savetxt(path, np.column_stack([T9_WAGONER, RATE_CENTRAL, RATE_ERROR]),
               fmt=["%.6e", "%.6e", "%.6e"], delimiter="   ", header=header)
    return path


print(f"detailed balance: alpha={ALPHA_DB:.6g} beta={BETA_DB:.6g} "
      f"gamma={GAMMA_DB:.6g}  Q={Q_MEV:.6g} MeV")

OUT_PATH = write_rate_table()
print(f"wrote {OUT_PATH}")
print(OUT_PATH.read_text().split("\n")[0])
print(OUT_PATH.read_text().split("\n")[1])
print(OUT_PATH.read_text().split("\n")[3])
```

- [ ] **Step 3: Run and compare the header against the shipped table**

Run: `python <scratchpad>/run_nb.py`

Expected `NOTEBOOK OK`, with the detailed-balance line reading

```
detailed balance: alpha=1.73492 beta=0 gamma=-46.7971  Q=4.03266 MeV
```

which must match line 2 of the *shipped* table, since both are derived from
the same nuclide data. Confirm with:

```bash
sed -n 2p primat/data/nuclear/tables/d_d__t_p/d_d__t_p_primat.txt
sed -n 2p generate_rates/rate_tables_out/tables/d_d__t_p/d_d__t_p_Mathematica-Sddp.txt
```

Expected: the two lines are byte-identical.

- [ ] **Step 4: Verify the generated table is loadable and well-shaped**

Run:

```bash
python -c "
import sys; sys.path.insert(0, '.')
import numpy as np
from primat.config import PRIMATConfig
from primat.network_data import available_rate_tables
cfg = PRIMATConfig(params={'user_nuclear_dir': 'generate_rates/rate_tables_out'})
print(available_rate_tables('d_d__t_p', cfg))
d = np.loadtxt('generate_rates/rate_tables_out/tables/d_d__t_p/d_d__t_p_Mathematica-Sddp.txt')
print(d.shape, d[0], d[-1])
"
```

Expected: the listing is `['d_d__t_p_Mathematica-Sddp.txt']` and the array
shape is `(60, 3)`.

Note the listing shows *only* the overlay's file, not the shipped
`_primat.txt`: for a reaction folder that exists in the overlay,
`available_rate_tables` reports the overlay folder's contents. Overlay
resolution is still per-*file* at solve time, so every other reaction keeps
falling back to the shipped tree — this is a listing quirk, not a shadowing
bug, and §8 does not depend on it.

- [ ] **Step 5: Verify the overwrite guardrail**

Run `python <scratchpad>/run_nb.py` a second time with `OVERWRITE = False`
still set. Expected: exit code 1 with `already exists.  Set OVERWRITE = True`.
Then set `OVERWRITE = True` in §3 and re-run; expected `NOTEBOOK OK`. Leave
`OVERWRITE = True` in the committed notebook so it is re-runnable.

- [ ] **Step 6: Commit**

Only the notebook is committed — `generate_rates/rate_tables_out/` is a
generated artifact and is gitignored in Task 7.

```bash
git add generate_rates/thermal_average.ipynb
git commit -m "generate_rates: write primat-format rate tables with detailed-balance header"
```

---

### Task 6: §8 worked example, validation plots, and the second d+d channel

**Files:**
- Modify: `generate_rates/thermal_average.ipynb`

**Interfaces:**
- Consumes: everything from Tasks 1–5.
- Produces: nothing consumed downstream; the section is self-contained validation.

- [ ] **Step 1: Append the §8 markdown cell**

````markdown
## §8 — Validation against the shipped d+d tables

The §3 worked example is the `d + d → t + p` S-factor from the Mathematica
notebook. Here it is compared with primat's shipped
`d_d__t_p_primat.txt` (Gómez et al. 2017, a full R-matrix evaluation), and the
same is done for the second d+d channel with

```
S_ddn(E) = 0.05225 + 0.3655 E − 0.1799 E² + 0.05832 E³ − 0.007393 E⁴
```

A simple polynomial S-factor is *not* expected to reproduce an R-matrix
evaluation to better than roughly 10–20% across the whole grid; the point of
this comparison is to confirm the machinery — units, Gamow factor, quadrature,
detailed balance — not the fit. A ratio flat to within tens of percent
validates the notebook; a ratio off by orders of magnitude, or with the wrong
temperature slope, means something is broken.
````

- [ ] **Step 2: Append the §8 code cell**

```python
import matplotlib.pyplot as plt

# The two d+d S-factor polynomials from the Mathematica notebook
# Thermal-Average.nb, in MeV*barn with E in MeV.
_DD_EXAMPLES = {
    "d_d__t_p":   np.array([0.05520, 0.2151, -0.02555]),
    "d_d__He3_n": np.array([0.05225, 0.3655, -0.1799, 0.05832, -0.007393]),
}

fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex="col")

for col, (reac, coeffs) in enumerate(_DD_EXAMPLES.items()):
    # Rebind the module-level knobs the kernel closes over, so the same code
    # path that produced the output file is exercised for both channels.
    MODE = "S"
    cross_section = lambda E, th: np.polyval(th[::-1], E)   # th[0] + th[1]E + ...
    REACTANTS, PRODUCTS = reaction_species(reac)
    MU_MEV = reduced_mass_MeV(*REACTANTS)
    Z1Z2 = charge(REACTANTS[0]) * charge(REACTANTS[1])
    X_NODES, X_WEIGHTS = _build_x_quadrature(T9_WAGONER)

    ours = thermal_average(T9_WAGONER, coeffs)

    ref = np.loadtxt(REPO / f"primat/data/nuclear/tables/{reac}/{reac}_primat.txt")
    ref_rate = np.exp(np.interp(np.log(T9_WAGONER),
                                np.log(ref[:, 0]), np.log(ref[:, 1])))

    ax = axes[0, col]
    ax.loglog(T9_WAGONER, ours, label="this notebook (polynomial S)")
    ax.loglog(ref[:, 0], ref[:, 1], "--", label="primat shipped table")
    ax.set_title(reac)
    ax.set_ylabel(r"$N_A\langle\sigma v\rangle$  [cm$^3$ mol$^{-1}$ s$^{-1}$]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, col]
    ax.semilogx(T9_WAGONER, ours / ref_rate)
    ax.axhline(1.0, color="k", lw=0.8)
    ax.set_xlabel(r"$T_9$")
    ax.set_ylabel("ratio  (ours / primat)")
    ax.set_ylim(0.5, 1.5)
    ax.grid(alpha=0.3)

    # The polynomial fits are only meant to hold over the BBN-relevant range;
    # quote the agreement there rather than over the whole grid.
    sel = (T9_WAGONER >= 0.1) & (T9_WAGONER <= 10.0)
    dev = np.max(np.abs(ours[sel] / ref_rate[sel] - 1.0))
    print(f"{reac}: max deviation over 0.1 <= T9 <= 10 is {100 * dev:.1f}%")
    assert dev < 0.30, (
        f"{reac} deviates from the shipped table by {100 * dev:.0f}%, far more "
        "than a polynomial-vs-R-matrix difference explains -- the machinery is "
        "probably wrong (units, Gamow factor, or quadrature)."
    )

fig.tight_layout()
plt.show()
```

- [ ] **Step 3: Append the §8b markdown cell, then the §8b code cell**

The prose below (down to the `large.txt` sentence) goes in a **markdown**
cell; the Python block inside it goes in a **code** cell immediately after.

````markdown
### Using the new table in a primat run

Which table file a run uses is decided by the **network list file**, whose
lines are `<reaction>, <table filename>`. So pointing primat at the new rate
takes two things: a network file naming it, and `user_nuclear_dir` pointing at
the overlay this notebook wrote into.

The cell below writes such a network file next to the table — a copy of the
shipped `small` network with this one reaction's filename swapped — and runs
BBN with it. Every *other* reaction still resolves to its shipped table,
because overlay resolution is per-file.

```python
import shutil
from primat.backend import run_bbn

overlay = Path(OUTDIR).parent                     # .../rate_tables_out
(overlay / "networks").mkdir(parents=True, exist_ok=True)

base = REPO / "primat/data/nuclear/networks/small.txt"
lines = []
for line in base.read_text().splitlines():
    reac = line.split(",")[0].strip()
    if reac == REACTION:
        line = f"{REACTION}, {OUT_PATH.name}"
    lines.append(line)
net = overlay / "networks" / "custom_rate.txt"
net.write_text("\n".join(lines) + "\n")

res = run_bbn(params={"network": "custom_rate",
                      "user_nuclear_dir": str(overlay)})
print(f"with {OUT_PATH.name}:  D/H = {res['DoH']:.7e}, YP = {res['YPBBN']:.8f}")

ref = run_bbn(params={"network": "small"})
print(f"shipped small network: D/H = {ref['DoH']:.7e}, YP = {ref['YPBBN']:.8f}")
```

This only works if `REACTION` is one of the 12 small-network reactions; for
anything else, start from `large.txt` instead of `small.txt`.
````

- [ ] **Step 4: Run and check both channels agree with the shipped tables**

Run: `python <scratchpad>/run_nb.py`

Expected two lines of the form
```
d_d__t_p: max deviation over 0.1 <= T9 <= 10 is <N>%
d_d__He3_n: max deviation over 0.1 <= T9 <= 10 is <N>%
```
with both `<N>` below 30, and `NOTEBOOK OK`.

If a deviation exceeds 30%, do **not** relax the assertion. Debug in this
order: (a) confirm the §5 self-tests still pass — if they do the quadrature is
fine; (b) print `sigma_of_E_cm2` at E = 0.1 MeV and check it against the hand
value in Task 2 Step 5; (c) check the Gamow exponent sign.

- [ ] **Step 5: Check the §8b BBN run reproduces the shipped result**

§8b is a live code cell, so `run_nb.py` already executed it. Its two printed
lines must agree, because the worked example's polynomial S-factor is close to
— but not identical with — the shipped table:

```
with d_d__t_p_Mathematica-Sddp.txt:  D/H = 2.4...e-05, YP = 0.246...
shipped small network:               D/H = 2.4359107e-05, YP = 0.24699714
```

The **second** line is the hard check: it must read exactly
`D/H = 2.4359107e-05, YP = 0.24699714`, the shipped `network="small"` result
verified during planning. If it does not, the overlay is leaking into the
reference run — investigate before continuing.

The first line should differ from it by no more than a few percent in D/H
(`d_d__t_p` competes with `d_d__He3_n` for deuterium burning, so a ~10%
rate change moves D/H by ~1%). A first line *identical* to the second means
the overlay was not picked up at all — check that
`user_nuclear_dir` points at `rate_tables_out`, not at its `tables/`
subdirectory.

- [ ] **Step 6: Commit**

```bash
git add generate_rates/thermal_average.ipynb
git commit -m "generate_rates: validate thermal-average notebook against shipped d+d tables"
```

---

### Task 7: §9 README pointer and final end-to-end check

**Files:**
- Modify: `generate_rates/README.md:41` (after the `PRIMAT-Main.m` bullet)
- Modify: `generate_rates/thermal_average.ipynb` (final clean re-execution)

**Interfaces:**
- Consumes: the finished notebook.
- Produces: the committed deliverable.

- [ ] **Step 1: Add the README entry**

In `generate_rates/README.md`, insert this bullet into the "Pipeline map"
list, immediately **before** the `PRIMAT-Main.m` bullet (so the notebook sits
with the other runnable entry points rather than with the raw source data):

```markdown
- **`thermal_average.ipynb`** — the user-facing entry point for adding a rate
  from your *own* cross-section data, as opposed to the bulk regeneration
  above. Give it an astrophysical S-factor `S(E)` (or a cross-section
  `σ(E)`) plus a parameter covariance, and it computes the
  Maxwell–Boltzmann-averaged `N_A<σv>` on the 60-point Wagoner temperature
  grid, propagates the uncertainty by Monte Carlo into the table's
  multiplicative error column, and writes a primat-format table with the
  correct detailed-balance header. All nuclide data (masses, charges, spins,
  Q-values) is read from a live `PRIMATConfig`, so it cannot drift from the
  solver's own nuclear data.

  Output goes to `generate_rates/rate_tables_out/` (gitignored), laid out as a
  ready-made `user_nuclear_dir` overlay — the shipped `primat/data/` tree is
  never modified. The notebook's last section builds a matching network list
  file and runs BBN with the new rate to show the effect. This notebook
  supersedes the Mathematica `Thermal-Average.nb`, whose loader stub
  `Thermal-Average.m` is kept here for reference.
```

- [ ] **Step 2: Gitignore the generated overlay**

Append to `.gitignore`:

```
# Rate tables generated by generate_rates/thermal_average.ipynb
generate_rates/rate_tables_out/
```

- [ ] **Step 3: Verify the README claims are true**

Run:

```bash
grep -n "thermal_average.ipynb" generate_rates/README.md
ls generate_rates/thermal_average.ipynb
ls generate_rates/rate_tables_out/tables/d_d__t_p/ generate_rates/rate_tables_out/networks/
git check-ignore -v generate_rates/rate_tables_out/tables/d_d__t_p/d_d__t_p_Mathematica-Sddp.txt
```

Expected: the grep hits; the notebook exists; the overlay contains
`d_d__t_p_Mathematica-Sddp.txt` and `custom_rate.txt`; and `git check-ignore`
prints the `.gitignore` rule, confirming the generated tree stays untracked.

- [ ] **Step 4: Clear outputs and do a final clean run**

Clearing execution counts and outputs keeps the committed notebook's diff
readable and proves it runs from a cold start.

```bash
python -c "
import nbformat
nb = nbformat.read('generate_rates/thermal_average.ipynb', as_version=4)
for c in nb.cells:
    if c.cell_type == 'code':
        c.outputs = []
        c.execution_count = None
nbformat.write(nb, 'generate_rates/thermal_average.ipynb')
print('cleared')
"
python <scratchpad>/run_nb.py
```

Expected: `cleared`, then the full output sequence ending in `NOTEBOOK OK`,
with all three §5 self-tests passing and both §8 deviations under 30%.

- [ ] **Step 5: Confirm nothing outside `generate_rates/` was touched**

Run: `git status --short`

Expected: the only modified/added paths are
`generate_rates/thermal_average.ipynb`, `generate_rates/README.md`, and
`.gitignore`. Nothing under `primat/` (source *or* data) and nothing under
`primat-c/`. In particular `primat/data/nuclear/tables/` must be untouched —
if it is not, `OUTDIR` was left pointing at the shipped tree; revert those
files and fix §3.

- [ ] **Step 6: Commit**

```bash
git add generate_rates/README.md generate_rates/thermal_average.ipynb .gitignore
git commit -m "generate_rates: document thermal_average.ipynb in the pipeline map"
```

---

## Self-review notes

**Spec coverage.** §1 → Task 1 Step 3. §2 → Task 1 Step 4. §3 → Task 2 Steps
1–2. §4 → Task 2 Steps 3–4. §5 → Task 3. §6 → Task 4. §7 → Task 5. §8 → Task
6. §9 → Task 7. All four guardrails are exercised: bad `REACTION` (Task 2
Step 7), `Z₁Z₂ = 0` warning (Task 2 Step 6), non-PSD `COV` (Task 4 Step 4),
negative-σ clipping (reported in Task 4 Step 2's code and surfaced in §8's
S(E) inspection). The overwrite guardrail is checked in Task 5 Step 5.

**Output location.** The notebook writes to the untracked overlay
`generate_rates/rate_tables_out/`, never into `primat/data/`. Task 5 Step 6
commits only the notebook; Task 7 Step 2 gitignores the overlay and Step 5
verifies `primat/` is untouched.

**Verified during planning, not assumed.** Table-variant selection is by
*network list file* (`networks/<name>.txt` lines of the form
`<reaction>, <filename>`), not by a `params` key — an earlier draft of §8b
invented a `rate_table_<reaction>` parameter that does not exist. The overlay
route was checked end to end with a real run: `network="custom_rate"` plus
`user_nuclear_dir` pointing at the overlay reproduces the shipped
`network="small"` result exactly (D/H = 2.4359107e-05, YP = 0.24699714) when
the overlaid table is a copy of the shipped one. Task 6 Step 5 pins those two
numbers.

**Known deviation from the spec.** The spec's §3 sketch showed `COV = None`;
the plan ships a small non-zero `COV` in the worked example instead, so that
the Monte-Carlo path is exercised on every clean run rather than being dead
code in the delivered notebook. The `COV = None` branch is still implemented
and documented.

**Grid-size correction.** The spec initially said 59 Wagoner points; the
transcribed list has **60**. The spec has been corrected; 60 is authoritative.
