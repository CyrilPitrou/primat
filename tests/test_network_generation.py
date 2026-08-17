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

The CSV-based tests skip if ``primat/data/csv/`` is absent. That directory is
git-tracked and ships in the wheel, so the guard is a safety net for a
hand-assembled ``data_dir``, not a routine skip -- if these tests ever *do*
skip in CI, something has gone wrong with packaging, not with the generator.
"""
import csv
import os
import subprocess
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
    reason="primat/data/csv/ not present (it ships with the package; "
           "regenerate with python generate_rates/convert_ac2024_rates.py)",
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
    """Photons and beta-decay leptons resolve with A = 0 but the right charge.

    This is what lets check_conservation apply one uniform dZ = 0 rule to
    weak and beta-decay channels instead of special-casing them."""
    assert resolve_token("g").kind == "photon"
    bm, bp = resolve_token("Bm"), resolve_token("Bp")
    assert (bm.kind, bm.A, bm.Q) == ("lepton", 0, -1)
    assert (bp.kind, bp.A, bp.Q) == ("lepton", 0, +1)


def test_conservation_residual_zero_for_physical_reactions():
    """Real reactions -- radiative capture, a 3-body breakup, and two
    beta-minus decays -- have zero (dA, dQ) residual."""
    # n + p -> d + g ; t + t -> a + n + n ; a beta-minus decay.
    assert conservation_residual(["n", "p"], ["d", "g"]) == (0, 0)
    assert conservation_residual(["t", "t"], ["a", "n", "n"]) == (0, 0)
    assert conservation_residual(["Li9"], ["Be9", "Bm"]) == (0, 0)        # A,Q conserved
    assert conservation_residual(["N17"], ["O16", "n", "Bm"]) == (0, 0)


def test_conservation_residual_catches_violations():
    """The residual is non-zero for a baryon-violating and a charge-violating
    reaction: the check can actually fail, not just pass."""
    dA, dQ = conservation_residual(["n", "p"], ["He4"])      # 2 baryons vs 4
    assert dA != 0
    dA, dQ = conservation_residual(["p", "p"], ["d", "g"])   # charge 2 vs 1
    assert dQ != 0


def test_canonical_name_special_cases():
    """canonical_name maps (Z, A) to primat's spelling, including the two
    irregular cases n and p (which are not 'H1'-style names)."""
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
    """Every reaction in the generated catalog conserves baryon number and
    charge -- the formal gate on the whole generated network."""
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

    Only the ``*_primat.txt`` tables are globbed, i.e. those written by
    ``generate_rates/convert_ac2024_rates.py``. The ``*_parthenope3.0.txt``
    tables (used by the ``small_parthenope`` network) come from a different
    generator, ``generate_rates/parthenope3.0_extract/postprocess.py``, so
    their uncertainty column is not this converter's to guarantee.
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


# ---------------------------------------------------------------------------
# 5. The generation *commands* themselves (end-to-end)
#
# The tests above exercise the helper layer; these run the two scripts the way
# the README tells a maintainer to. Both scripts write into the source tree by
# default, which is exactly why they rot unnoticed: nothing had ever executed
# them, so `generate_qed_tables.py` spent months writing its tables to a
# directory (primat/rates/plasma) that had not existed since the data tree
# moved -- creating it on demand, printing "Done.", and leaving the shipped
# tables untouched.
# ---------------------------------------------------------------------------
_NUBASE = os.path.join(_GEN_DIR, "nubase_4.mas20.txt")
_AC2024_DAT = os.path.join(_GEN_DIR, "BBNRatesAC2024.dat")

_needs_nubase = pytest.mark.skipif(
    not os.path.isfile(_NUBASE),
    reason="generate_rates/nubase_4.mas20.txt not present",
)


@_needs_nubase
def test_nubase_halflives_are_read_at_the_documented_column_offsets():
    """GOAL: pin the fixed-width offsets of the NUBASE half-life field.

    NUBASE2020's format block numbers its columns from 1 (``70: 78   T #``),
    so the value is ``line[69:78]`` and the unit ``line[78:80]``. Slicing the
    value one column late still parses -- it just drops the leading digit of
    any half-life wide enough to fill the nine-character field, which is what
    made the generator report a spurious factor-2.5 disagreement for Ne18.
    Only a value whose field is *full* can catch that, hence Ne18 below
    ("1664.20  ms"); the two shorter ones guard the unit column and the
    year conversion.
    """
    from nuclide_table import load_nubase_halflives

    t12 = load_nubase_halflives(_NUBASE)

    # Ne18 -> F18: the field is full ("1664.20  ms"), so an off-by-one slice
    # reads 664.20 ms = 0.6642 s instead.
    assert t12[(10, 18)] == pytest.approx(1.6642, rel=1e-6)
    # Free neutron: a plain "s"-unit value (guards the unit column).  609.8 s
    # is NUBASE2020's half-life, i.e. ln(2) x 879.6 s -- consistent with the
    # tau_n = 878.4 s primat carries in DEFAULT_PARAMS to within the two
    # evaluations' spread.
    assert t12[(0, 1)] == pytest.approx(609.8, rel=1e-3)
    # C14: 5700 y, exercising the Julian-year conversion (365.2422 d).
    assert t12[(6, 14)] == pytest.approx(5700 * 86400 * 365.2422, rel=1e-3)
    # Stable nuclides carry no half-life.
    assert t12[(1, 1)] is None and t12[(2, 4)] is None


@_needs_nubase
def test_nubase_halflife_limits_are_not_reported_as_measurements():
    """GOAL: a bound (``>912.4 ys``) must not be served as a half-life.

    The marker sits in the field's first column -- exactly the character the
    old off-by-one slice discarded, silently turning limits into
    "measurements" that a cross-check would then compare against.
    """
    from nuclide_table import load_nubase_halflives

    assert load_nubase_halflives(_NUBASE)[(5, 20)] is None   # B20: ">912.4 ys"


def test_generate_qed_tables_default_output_dir_is_the_shipped_one():
    """GOAL: a plain ``generate_qed_tables.py`` run must overwrite the tables
    the solver actually loads, not a stray directory.

    The script creates its output directory on demand, so a stale default path
    fails silently -- it writes, reports success, and changes nothing. Pinning
    the constant against the shipped files is what makes that loud.
    """
    from generate_rates.generate_qed_tables import DEFAULT_PLASMA_DIR

    assert os.path.isdir(DEFAULT_PLASMA_DIR), (
        f"{DEFAULT_PLASMA_DIR} does not exist -- has the data tree moved? "
        f"A regeneration would silently write tables nothing reads.")
    for name in ("QED_pressure_correction_e2.txt",
                 "QED_pressure_correction_e3.txt"):
        assert os.path.isfile(os.path.join(DEFAULT_PLASMA_DIR, name))


@pytest.mark.slow
def test_generate_qed_tables_writes_fingerprinted_tables(tmp_path):
    """GOAL: the QED regeneration command still runs end to end.

    Uses a deliberately coarse grid (``--n-pts 8``): this is a "does the script
    work" check, not a physics check -- the *content* of the QED tables is
    pinned by ``test_qed_pressure.py`` and, cross-backend, by
    ``test_cache_parity.py``. What could rot here is the plumbing: the
    ``primat.qed_pressure`` entry points it calls, its argument parsing, and
    the fingerprint header ``save_qed_tables`` writes.
    """
    script = os.path.join(_GEN_DIR, "generate_qed_tables.py")
    out = tmp_path / "plasma"
    out.mkdir()
    result = subprocess.run(
        [sys.executable, script, "--n-pts", "8", "--output-dir", str(out)],
        cwd=tmp_path, capture_output=True, text=True, encoding="utf-8", timeout=300,
    )
    assert result.returncode == 0, result.stderr
    for name in ("QED_pressure_correction_e2.txt",
                 "QED_pressure_correction_e3.txt"):
        text = (out / name).read_text()
        assert "# fingerprint_hash:" in text
        assert '"n_pts":8' in text.replace(" ", "")
        # 8 grid rows plus the comment header.
        assert len([l for l in text.splitlines() if not l.startswith("#")]) == 8


@pytest.mark.slow
@_needs_ac2024
@_needs_nubase
def test_convert_ac2024_regenerates_the_shipped_tables_byte_for_byte(tmp_path):
    """GOAL: the shipped rate tables are exactly what the generator produces.

    ``primat/data/nuclear/tables/`` (428 tables + decays.txt), the three CSVs
    and ``networks/large.txt`` are all generated artifacts; nothing else checks
    that they can be *re*generated. A drift here means the committed data no
    longer corresponds to the script that claims to produce it -- so the next
    maintainer who reruns the generator gets an unexplained diff and has no way
    to tell which side is right.

    The script writes to paths relative to the current directory, so running it
    from ``tmp_path`` (with absolute input paths) regenerates a complete tree
    there without touching the repository. ~2 s.
    """
    import filecmp

    script = os.path.join(_GEN_DIR, "convert_ac2024_rates.py")
    result = subprocess.run(
        [sys.executable, script, "--input", _AC2024_DAT, "--nubase", _NUBASE],
        cwd=tmp_path, capture_output=True, text=True, encoding="utf-8", timeout=600,
    )
    assert result.returncode == 0, result.stderr
    # The generator's own checks must pass too -- they are the reason a bad
    # edit aborts instead of shipping.
    assert "formal check OK" in result.stdout
    assert "WARNING" not in result.stdout, result.stdout

    shipped_root = os.path.abspath(os.path.join(_ROOT, "primat", "data"))
    fresh_root = os.path.join(str(tmp_path), "primat", "data")

    compared = mismatched = 0
    for sub in (os.path.join("nuclear", "tables"),
                os.path.join("nuclear", "networks"),
                "csv"):
        fresh_dir = os.path.join(fresh_root, sub)
        for dirpath, _dirnames, filenames in os.walk(fresh_dir):
            for fname in filenames:
                fresh = os.path.join(dirpath, fname)
                shipped = os.path.join(
                    shipped_root, sub, os.path.relpath(fresh, fresh_dir))
                assert os.path.isfile(shipped), (
                    f"generator produced {fname}, which is not shipped")
                compared += 1
                if not filecmp.cmp(fresh, shipped, shallow=False):
                    mismatched += 1
    # 390 per-reaction tables (428 reactions less the 38 decays, which share
    # decays.txt) + decays.txt + 3 CSVs + networks/large.txt = 395.
    assert compared >= 390, f"only {compared} generated files found"
    assert mismatched == 0, (
        f"{mismatched} of {compared} regenerated files differ from the "
        f"shipped ones -- the committed data no longer matches the generator")
