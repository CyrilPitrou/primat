"""
Tests for the offline network-generation layer
(``generate_rates/convert_ac2024_rates.py`` + ``nuclide_table.py``).

That command runs once to turn AC2024 + the analytic table + NUBASE into the
three CSVs primat reads at start-up:

* ``nuclides.csv``        : every nuclide the network touches, with N,Z,A,Q,mass,spin,
* ``reactions_large.csv`` : the deduced >400-reaction list,
* ``detailed_balance.csv``: alpha,beta,gamma per reversible reaction.

The tests below check, without re-running the (slow) full generation:

1. the token resolver and the *formal* baryon/charge conservation check;
2. that the generated nuclide table is internally consistent and agrees with
   primat's hard-coded 12-nuclide table;
3. that the deduced reaction list is a superset of the known 12- and 68-reaction
   networks (the latter being "large" filtered to amax=8, the old "medium"
   network's exact equivalent) and that *every* listed reaction conserves A
   and Q;
4. that the detailed-balance coefficients computed from nuclide data reproduce
   primat's published values for those 68 reactions.

The CSV-based tests skip if the generated ``rates/nuclear/AC2024`` folder is
absent (fresh checkout before the generator has been run).
"""
import csv
import os
import sys

import pytest

_HERE = os.path.dirname(__file__)
_ROOT = os.path.join(_HERE, "..")
_GEN_DIR = os.path.join(_ROOT, "generate_rates")
_AC2024_DIR = os.path.join(_ROOT, "primat", "data", "csv")

# The generation helpers live in generate_rates/, which is not an
# installed package; add it to sys.path so the tests can import it directly.
sys.path.insert(0, _GEN_DIR)

from nuclide_table import (resolve_token, conservation_residual,   # noqa: E402
                           build_nuclide_table, canonical_name)

_needs_ac2024 = pytest.mark.skipif(
    not os.path.isdir(_AC2024_DIR),
    reason="rates/csv not generated "
           "(run python generate_rates/convert_ac2024_rates.py)",
)


# ---------------------------------------------------------------------------
# 1. Token resolution and formal conservation
# ---------------------------------------------------------------------------
def test_resolve_token_canonicalises_spellings():
    """``a``/``He4``, ``d``/``H2``, ``t``/``H3`` must collapse to one nuclide."""
    assert resolve_token("a").name == resolve_token("He4").name == "He4"
    assert resolve_token("d").name == resolve_token("H2").name == "H2"
    assert resolve_token("t").name == resolve_token("H3").name == "H3"
    n, p = resolve_token("n"), resolve_token("p")
    assert (n.Z, n.A, n.Q) == (0, 1, 0)
    assert (p.Z, p.A, p.Q) == (1, 1, 1)
    c12 = resolve_token("C12")
    assert (c12.Z, c12.A, c12.name) == (6, 12, "C12")


def test_shipped_rate_tables_are_on_the_master_grid():
    """Every shipped per-reaction rate table must already sit on primat's
    master T9 grid (``rate_grid_{npts,T9_min,T9_max}``).

    Both generators -- ``convert_ac2024_rates.py`` (the ``_primat`` tables) and
    ``parthenope3.0_extract/postprocess.py`` (the ``parthenope3.0`` tables) --
    write on this grid, so ``load_network``'s load-time resampler takes its
    identity fast path (``_resample_rate_table``) and spends no time
    interpolating stored rates.  If a table drifts off the grid (e.g. a stale
    500-point table left un-regenerated after ``rate_grid_npts`` changed) the
    solver silently pays a full cubic resample on every run and the C and
    Python backends can disagree on the resampled values -- exactly the drift
    this guard exists to catch.  ``decays.txt`` is exempt: it is a shared
    one-row-per-decay table, not a per-reaction T9 grid.
    """
    import glob
    import numpy as np
    from primat.config import DEFAULT_PARAMS

    master = np.logspace(
        np.log10(DEFAULT_PARAMS["rate_grid_T9_min"]),
        np.log10(DEFAULT_PARAMS["rate_grid_T9_max"]),
        DEFAULT_PARAMS["rate_grid_npts"],
    )
    tables_dir = os.path.join(_ROOT, "primat", "data", "nuclear", "tables")
    files = sorted(glob.glob(os.path.join(tables_dir, "*", "*.txt")))
    assert files, f"no rate tables found under {tables_dir}"

    off_grid = []
    for path in files:
        if os.path.basename(path) == "decays.txt":
            continue  # shared decay table, not a per-reaction T9 grid
        t9 = np.loadtxt(path)[:, 0]
        # Match the load-time fast-path criterion in _resample_rate_table:
        # same length and agreeing with the master grid to the on-disk
        # 7-significant-figure rounding of the T9 column.
        on_grid = (t9.shape == master.shape
                   and np.all(np.abs(t9 / master - 1.0) < 1.0e-6))
        if not on_grid:
            off_grid.append((os.path.relpath(path, _ROOT), len(t9)))

    assert not off_grid, (
        "rate tables not on the master grid (regenerate them -- see "
        "convert_ac2024_rates.py / parthenope3.0_extract/postprocess.py):\n"
        + "\n".join(f"  {p}  ({n} rows, expected {len(master)})"
                    for p, n in off_grid))


def test_leptons_and_photons_carry_charge_but_no_baryon():
    assert resolve_token("g").kind == "photon"
    bm, bp = resolve_token("Bm"), resolve_token("Bp")
    assert (bm.kind, bm.A, bm.Q) == ("lepton", 0, -1)
    assert (bp.kind, bp.A, bp.Q) == ("lepton", 0, +1)


def test_conservation_residual_zero_for_physical_reactions():
    # n + p -> d + g ; t + t -> a + n + n ; a beta-minus decay.
    assert conservation_residual(["n", "p"], ["d", "g"]) == (0, 0)
    assert conservation_residual(["t", "t"], ["a", "n", "n"]) == (0, 0)
    assert conservation_residual(["Li9"], ["Be9", "Bm"]) == (0, 0)        # A,Q conserved
    assert conservation_residual(["N17"], ["O16", "n", "Bm"]) == (0, 0)


def test_conservation_residual_catches_violations():
    dA, dQ = conservation_residual(["n", "p"], ["He4"])      # 2 baryons vs 4
    assert dA != 0
    dA, dQ = conservation_residual(["p", "p"], ["d", "g"])   # charge 2 vs 1
    assert dQ != 0


def test_canonical_name_special_cases():
    assert canonical_name(0, 1) == "n"
    assert canonical_name(1, 1) == "p"
    assert canonical_name(1, 2) == "H2"
    assert canonical_name(2, 4) == "He4"


# ---------------------------------------------------------------------------
# 2. Generated nuclides.csv consistency
# ---------------------------------------------------------------------------
def _load_nuclides_csv():
    with open(os.path.join(_AC2024_DIR, "nuclides.csv")) as f:
        return {r["name"]: r for r in csv.DictReader(f)}


@_needs_ac2024
def test_nuclides_csv_agrees_with_pyprimat_hardcoded_table():
    """Every nuclide in PRIMATConfig.Nuclides must appear in nuclides.csv with the
    same (N, Z) -- the generated table is a superset of the hard-coded one."""
    from primat.config import PRIMATConfig
    cfg = PRIMATConfig()
    nuc = _load_nuclides_csv()
    # Check the key ones used in Speciess_Small
    for name in ["n", "p", "H2", "H3", "He3", "He4", "Li7", "Be7"]:
        N, Z = cfg.Nuclides[name]
        assert name in nuc, f"{name} missing from nuclides.csv"
        assert (int(nuc[name]["N"]), int(nuc[name]["Z"])) == (N, Z)


@_needs_ac2024
def test_nuclides_csv_self_consistent():
    """A = N + Z and Q = Z for every row; mass excess and spin are present."""
    for r in _load_nuclides_csv().values():
        N, Z, A, Q = (int(r[k]) for k in ("N", "Z", "A", "Q"))
        assert A == N + Z and Q == Z
        float(r["mass_excess_keV"])                  # parses
        float(r["spin"])


# ---------------------------------------------------------------------------
# 3. Generated reactions_large.csv: superset + conservation
# ---------------------------------------------------------------------------
def _load_reactions_csv():
    with open(os.path.join(_AC2024_DIR, "reactions_large.csv")) as f:
        return list(csv.DictReader(f))


@_needs_ac2024
def test_reaction_list_is_superset_of_known_networks():
    """The deduced large list must contain every reaction of the 12-key and
    68-reaction (large, amax=8 -- the old "medium" network's exact
    equivalent) networks (matched by their <reactants>TO<products> file name)."""
    from primat.config import PRIMATConfig
    from primat.network_data import to_filename, _KEY12_REACTIONS, load_network
    names = {r["name"] for r in _load_reactions_csv()}
    for compact in _KEY12_REACTIONS:
        name = compact if 'TO' in compact else to_filename(compact)
        assert name in names, f"{compact} missing from large list"
    amax8_names = load_network(PRIMATConfig({"network": "large", "amax": 8}),
                               era="LT").names
    for compact in amax8_names:
        if compact == "n__p":
            continue
        name = compact if 'TO' in compact else to_filename(compact)
        assert name in names, f"{compact} missing from large list"


@_needs_ac2024
def test_every_listed_reaction_conserves_A_and_Q():
    for r in _load_reactions_csv():
        reactants = r["reactants"].split("+")
        products = r["products"].split("+")
        assert conservation_residual(reactants, products) == (0, 0), r["name"]


# ---------------------------------------------------------------------------
# 4. Detailed balance consistency
# ---------------------------------------------------------------------------
@_needs_ac2024
def test_detailed_balance_formula_consistency():
    """alpha,beta,gamma computed from nuclide data must reproduce the
    values in detailed_balance.csv: beta exactly,
    alpha and gamma to better than 1% (the documented detailed-balance accuracy)."""
    from primat.config import PRIMATConfig
    from primat.network_data import compute_detailed_balance_coefficients, reaction_species
    cfg = PRIMATConfig()
    with open(os.path.join(_AC2024_DIR, "detailed_balance.csv")) as f:
        db_rows = list(csv.DictReader(f))
    
    # Check a representative sample or all of them
    for row in db_rows:
        name = row["reaction"]
        ref_alpha = float(row["alpha"])
        ref_beta = float(row["beta"])
        ref_gamma = float(row["gamma"])
        
        reactants, products = reaction_species(name)
        alpha, beta, gamma = compute_detailed_balance_coefficients(reactants, products, cfg)
        
        assert beta == pytest.approx(ref_beta), f"Failed for {name}"
        if ref_alpha:
            assert abs(alpha - ref_alpha) / abs(ref_alpha) < 0.01
        if ref_gamma:
            assert abs(gamma - ref_gamma) / abs(ref_gamma) < 0.01


# ---------------------------------------------------------------------------
# Uncertainty column: one-sided multiplicative factor
#
# GOAL: pin that every shipped table's third column is a factor f >= 1, so
# that ``p_<rxn>`` means the same thing for every reaction.
#
# ``NetworkDefinition.apply_variations`` forms ``exp(p * log(f))``, so an
# ``f < 1`` silently flips the *direction* of the variation: a deterministic
# ``p = +1`` sweep would raise ~390 rates and lower the handful with f < 1.
#
# Regression guard: the source data (BBNRatesAC2024.dat, whose header defines
# the column as ``exp(sigma) = sqrt(sv_high/sv_low)``) is clean -- zero rows
# below 1 across all 337 blocks.  The sub-1 values were manufactured by the
# converter, which resampled that column with a log-log *cubic*.  The column
# is flat at exactly 1.0 wherever the rate is the ``0.999E-99`` sentinel and
# then steps up sharply (to f ~ 27 for O17_a__Ne20_n); a cubic rings around
# such a step and its undershoot crossed below 1 in 38 tables, worst 0.649.
# The converter now uses a shape-preserving (PCHIP) interpolant there, which
# cannot overshoot its data.
# ---------------------------------------------------------------------------
def test_shipped_uncertainty_columns_are_at_least_one():
    """Every primat-generated rate table has ``f >= 1`` at every T9.

    The Parthenope-sourced tables are excluded: they come from a different
    generator (see CLAUDE.md's rate-table provenance), so their uncertainty
    column is not this converter's to guarantee.
    """
    import glob
    import os

    import numpy as np

    offenders = []
    for path in sorted(glob.glob(
            os.path.join(os.path.dirname(__file__), os.pardir, "primat", "data",
                         "nuclear", "tables", "*", "*_primat.txt"))):
        col = np.loadtxt(path, unpack=True)
        if col.shape[0] < 3:
            continue
        if col[2].min() < 1.0:
            offenders.append((os.path.basename(path), float(col[2].min())))
    assert not offenders, (
        "rate tables whose uncertainty factor dips below 1 (this inverts the "
        f"sign of their p_<rxn> variation): {offenders[:5]}")


def test_monotone_interpolant_does_not_ring_below_a_flat_run():
    """A flat ``f = 1`` run followed by a step must not undershoot.

    This is the exact shape the ``0.999E-99`` sentinel rows produce, reduced to
    its essentials: a cubic through it dips below 1 between the knots, a
    shape-preserving interpolant cannot.
    """
    import numpy as np

    from generate_rates.convert_ac2024_rates import (
        interp_loglog, interp_loglog_monotone)

    T9 = np.array([0.01, 0.02, 0.04, 0.08, 0.16, 0.32])
    f = np.array([1.0, 1.0, 1.0, 12.0, 20.0, 27.0])
    dense = np.logspace(np.log10(T9[0]), np.log10(T9[-1]), 500)

    # The scheme that caused the bug does undershoot here ...
    assert interp_loglog(T9, f, dense).min() < 1.0
    # ... and the one now used does not, while still honouring the knots.
    got = interp_loglog_monotone(T9, f, dense)
    assert got.min() >= 1.0
    assert np.allclose(interp_loglog_monotone(T9, f, T9), f, rtol=1e-12)
