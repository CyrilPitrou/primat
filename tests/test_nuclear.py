"""
Tests for the auto-derivation fallback in ``reaction_stoichiometry`` and the
duplicate-entry check in ``load_network``.

``reaction_stoichiometry`` normally locates the reactant/product split via the
detailed-balance exponent ``beta`` looked up in ``detailed_balance.csv``.  For a
reaction added directly to a network file that has no ``detailed_balance.csv``
entry yet, the literal ``"TO"`` token in the compact name marks the split
directly.  This fallback path is exercised here with synthetic reaction names
that are not present in ``detailed_balance.csv``:

1. ``ppTOdBp`` (p + p -> d + e+, the pp-fusion reaction): a balanced synthetic
   name.  The fallback must derive ``({"p": 2}, {"H2": 1, "Bp": 1})`` and pass
   the same A/Z conservation check as ``check_conservation``.
2. ``ppTOd`` (p + p -> d, missing the positron): an unbalanced synthetic name
   (charge mismatch).  The fallback must raise ``ValueError`` naming the
   reaction and the A/Z imbalance, instead of silently returning bad
   stoichiometry or raising a cryptic ``KeyError``.

``load_network`` also gained a check that rejects a network reaction list
containing the same entry twice (most likely a copy-paste mistake in a network
file), raising ``ValueError`` instead of silently dropping or double-counting
the repeat.
"""
import pytest

from primat.config import PRIMATConfig
from primat.network_data import (
    SPECIES_SMALL,
    load_network,
    reaction_stoichiometry,
)


def test_auto_derived_stoichiometry_for_unknown_reaction():
    """A synthetic name absent from detailed_balance.csv is split at the
    literal "TO" token, with H2/H3/He4 aliases (d/t/a) resolved and Bp/Bm
    bookkeeping tokens kept as-is."""
    react, prod = reaction_stoichiometry("ppTOdBp")
    assert react == {"p": 2}
    assert prod == {"H2": 1, "Bp": 1}


def test_auto_derived_stoichiometry_conserves_A_and_Z():
    """The A/Z totals of the auto-derived reactants and products agree, using
    the same nuclide (N, Z) data as check_conservation: p+p (A=2, Z=2) vs
    d + e+ (A=2, Z=1+1=2)."""
    from primat.network_data import _reaction_catalog, _default_data_dir, _LEPTON_Z
    _, _, _, nuc_NZ, _, _ = _reaction_catalog(_default_data_dir())

    def totals(counts):
        A = Z = 0
        for tok, mult in counts.items():
            if tok in _LEPTON_Z:
                Z += _LEPTON_Z[tok] * mult
                continue
            n, z = nuc_NZ[tok]
            A += (n + z) * mult
            Z += z * mult
        return A, Z

    react, prod = reaction_stoichiometry("ppTOdBp")
    assert totals(react) == totals(prod) == (2, 2)


def test_unbalanced_synthetic_reaction_raises_value_error():
    """``ppTOd`` (p + p -> d) drops the positron, leaving Z unbalanced
    (reactants Z=2, products Z=1).  The fallback must raise ValueError naming
    the reaction and the imbalance, not silently return bad stoichiometry."""
    with pytest.raises(ValueError, match="ppTOd.*conserve"):
        reaction_stoichiometry("ppTOd")


def test_reaction_with_no_TO_token_raises_value_error():
    """A name with no detailed_balance.csv entry and no "TO" separator cannot
    be split into reactants/products at all."""
    with pytest.raises(ValueError, match="cannot be derived"):
        reaction_stoichiometry("ppdtg")


def test_duplicate_reaction_entry_raises_value_error():
    """A network reaction list containing the same entry twice raises
    ValueError naming the duplicated entry, instead of silently dropping or
    double-counting it."""
    cfg = PRIMATConfig({"network": "small", "verbose": False})
    with pytest.raises(ValueError, match="n_p__d_g.*already present"):
        load_network(cfg, era="LT",
                      reaction_names=["n_p__d_g", "n_p__d_g", "d_d__He3_n"])


def test_small_network_reports_exactly_its_eight_nuclides():
    """The `small` network evolves exactly the 8 SPECIES_SMALL nuclides
    (n, p, H2, H3, He3, He4, Li7, Be7) and must *report* exactly those 8 in
    ``Y_final`` -- no phantom extras.

    Regression guard for a past bug: ``_solve_LT`` used to zero-fill
    ``Y_final`` with all of ``SPECIES_MT`` (SPECIES_SMALL + He6/Li8/Li6/B8),
    padding the small network up to a spurious 12-nuclide ``Y_final`` even
    though its ODE state vector (``abundance_names``) is only the 8 above.
    That mismatch also leaked into the MC nuclide set (``nuclide_names =
    list(Y_final.keys())``) and broke C-vs-Python backend parity (the C
    backend correctly reports 8).  The padding was narrowed to SPECIES_SMALL,
    so the invariant below must now hold: ``Y_final`` keys are *exactly* the
    evolved species, never a superset.
    """
    from primat import PRIMAT

    pr = PRIMAT({"network": "small", "verbose": False})
    pr.solve()

    # The ODE state vector and the reported abundances must be the same set of
    # nuclides -- Y_final must not carry any species the network never evolved.
    assert set(pr.nuclear.Y_final) == set(pr.nuclear.abundance_names)
    # And for `small` that set is exactly the 8 SPECIES_SMALL members.
    assert set(pr.nuclear.Y_final) == set(SPECIES_SMALL)
    assert len(pr.nuclear.Y_final) == 8
    # None of the four SPECIES_MT-only extras should appear.
    for phantom in ("He6", "Li8", "Li6", "B8"):
        assert phantom not in pr.nuclear.Y_final


# ---------------------------------------------------------------------------
# DT-era decay matrix: baryon-number conservation
#
# GOAL: pin the convention that ``Y`` is the *number* abundance per baryon
# (``Y_s = n_s/n_B``, ``sum_s A_s Y_s = 1``), which is what the LT/MT
# right-hand side uses, and which fixes the form of the DT-era decay matrix:
# each decay's product gain is the bare stoichiometric multiplicity, with no
# ``A_P/A_X`` mass weighting.
#
# Regression guard: ``_build_decay_matrix`` used to multiply the gain term by
# ``A_P/A_X`` *in addition to* the multiplicity, on the mistaken premise that
# Y was a mass fraction.  That silently destroyed baryon number for every
# decay whose products differ in mass number from the parent -- ``Li8 -> a+a``
# produced one alpha instead of two, ``C9 -> a+a+p`` lost 4/9 of the baryons.
# The 33 ordinary beta decays have ``A_P == A_X``, so they were unaffected and
# the error went unnoticed.  The column-sum identity below is exactly the
# check that catches it.
# ---------------------------------------------------------------------------
def _decay_matrix(network="large", amax=None):
    """Build the DT-era decay matrix without running a full BBN solve."""
    import numpy as np
    from primat.network_data import UpdateNuclearRates
    from primat.nuclear_network import NuclearNetwork

    cfg = PRIMATConfig({"network": network, "amax": amax, "verbose": False})
    nucl = UpdateNuclearRates(cfg)
    nn = NuclearNetwork.__new__(NuclearNetwork)   # only cfg is needed
    nn.cfg = cfg
    net = nucl._lt_net
    D = NuclearNetwork._build_decay_matrix(nn, net)
    A = (net.N + net.Z).astype(float)
    return net, np.asarray(A), D


def test_decay_matrix_conserves_baryon_number():
    """Every parent column of D must satisfy ``sum_s A_s D[s, X] = 0``.

    Leptons (Bm/Bp) and photons carry A = 0, so they remove no baryon number,
    and every decays.txt reaction balances A between parent and products.  The
    identity is therefore exact, not approximate.
    """
    import numpy as np

    net, A, D = _decay_matrix()
    colsum = A @ D
    # Scale the tolerance by each column's own decay rate, so a fast decay
    # (rate ~ 70 s^-1) is not held to the same absolute bound as a slow one.
    scale = np.maximum(np.abs(np.diag(D)), 1.0)
    worst = int(np.argmax(np.abs(colsum) / scale))
    assert np.allclose(colsum, 0.0, atol=1e-9 * scale, rtol=0), (
        f"baryon number not conserved in column {net.species[worst]!r}: "
        f"sum_s A_s D[s,X] = {colsum[worst]:.6e} (diagonal {D[worst, worst]:.6e})")


def test_decay_matrix_multi_fragment_decay_yields_both_alphas():
    """``Li8 -> a + a + Bm`` must put *two* alphas into the He4 row.

    This is the decay the ``A_P/A_X`` bug halved (``A_He4/A_Li8 = 4/8``
    exactly cancelled the multiplicity 2), so it is the sharpest single-number
    guard available.  ``B8 -> a + a + Bp`` is the mirror case.
    """
    net, A, D = _decay_matrix()
    for parent in ("Li8", "B8"):
        X = net.species.index(parent)
        P = net.species.index("He4")
        lam = -D[X, X]                       # the parent's own loss rate
        assert lam > 0.0, f"{parent} has no decay rate"
        assert D[P, X] == pytest.approx(2.0 * lam, rel=1e-12), (
            f"{parent} -> a + a should give D[He4,{parent}] = 2*lambda")


def test_decay_matrix_beta_decay_is_one_to_one():
    """A mass-preserving beta decay (``C14 -> N14``) gives a unit gain.

    The complement of the test above: where ``A_P == A_X`` the old and new
    formulas agree, so this pins that the fix did not disturb the 33 ordinary
    beta decays.
    """
    net, A, D = _decay_matrix()
    X = net.species.index("C14")
    P = net.species.index("N14")
    lam = -D[X, X]
    assert D[P, X] == pytest.approx(lam, rel=1e-12)


# ---------------------------------------------------------------------------
# Network-list parsing and rate-table resampling
#
# GOAL: two silent-corruption paths in load_network's front end -- a repeated
# reaction being integrated twice, and a user-supplied rate table narrower than
# the master grid being extrapolated without bound.
# ---------------------------------------------------------------------------
def test_duplicate_entry_detected_across_alternate_table_suffix():
    """A repeat differing only in its ``", filename"`` column is still a repeat.

    Regression guard: the duplicate check used to compare *raw entries*, so
    ``"n_p__d_g"`` alongside ``"n_p__d_g, n_p__d_g_primat.txt"`` slipped
    through and the reaction was compiled as two identical rows -- doubling its
    flux, with no warning.  The C backend has always keyed this check on the
    bare name (``cpr_load_network_list``), so this is a parity guard too.
    """
    cfg = PRIMATConfig({"network": "small", "verbose": False})
    with pytest.raises(ValueError, match="n_p__d_g.*already present"):
        load_network(cfg, era="LT",
                      reaction_names=["n_p__d_g",
                                      "n_p__d_g, n_p__d_g_primat.txt",
                                      "d_p__He3_g"])


def test_shipped_tables_resample_to_an_exact_identity():
    """A table already on the master grid must come back byte-identical.

    This is the fast path every shipped table takes; it is what keeps the two
    backends bit-for-bit equal on rate loading, so it must not regress when the
    out-of-range branches change.
    """
    import numpy as np
    from primat.network_data import _resample_rate_table

    grid = np.logspace(-3, 1, 1000)
    rate = np.exp(-1.0 / grid)                    # any smooth positive function
    out = _resample_rate_table(grid, rate, grid, label="test")
    assert np.array_equal(out, rate)


def test_short_rate_table_extrapolates_by_end_slope_and_warns():
    """A table narrower than the master grid warns, and stays bounded.

    Regression guard: out-of-range points used to be produced by evaluating the
    **cubic** spline outside its data, which is unbounded -- the case below
    came out a factor 3.2 low at the bottom of the grid, silently.  Continuing
    the table's end slope in log-log keeps the error to tens of percent and the
    warning makes it visible.
    """
    import numpy as np
    from primat.network_data import _resample_rate_table

    # A pure power law in log-log: the end slope is the exact continuation, so
    # any departure from it is the interpolator's own doing.
    grid = np.logspace(-3, 1, 200)
    T9_src = np.logspace(np.log10(0.05), np.log10(5.0), 100)
    out = pytest.warns(UserWarning, match=r"extrapolated by continuing")
    with out:
        got = _resample_rate_table(T9_src, T9_src ** 2.5, grid, label="x_y__z_g")

    # Exact for a power law, everywhere including the extrapolated ends.
    assert np.allclose(got, grid ** 2.5, rtol=1e-8)
    # And bounded: no runaway at either edge.
    assert np.isfinite(got).all() and (got > 0).all()


def test_full_coverage_rate_table_does_not_warn():
    """A table spanning the whole master grid must resample silently."""
    import numpy as np
    import warnings as _warnings
    from primat.network_data import _resample_rate_table

    grid = np.logspace(-3, 1, 200)
    T9_src = np.logspace(-3, 1, 137)             # different npts, same span
    with _warnings.catch_warnings():
        _warnings.simplefilter("error")           # any warning fails the test
        _resample_rate_table(T9_src, T9_src ** 2.5, grid, label="x_y__z_g")
