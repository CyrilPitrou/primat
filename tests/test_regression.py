"""
Regression tests on the final abundances, against the published reference
values in tests/README.md's "Validation reference" section.

Three layers:

* Default-precision sanity checks (via the ``solved_small`` / ``solved_large``
  fixtures) with loose tolerances — cheap, catch gross regressions.
* The **per-nuclide** default-precision check
  (``test_per_nuclide_abundances_match_the_reference_table``), which pins every
  cell of the published per-nuclide table at 1e-4 relative. This is the live
  half of that table's guard; the static half (README text vs. the constants)
  is ``tests/test_docs_consistency.py``'s
  ``test_per_nuclide_reference_table_matches_reference_constants``.
* High-precision *reference* checks (``reference`` marker) that rerun at the
  exact settings used to produce the published numbers
  (numerical_precision=1e-10, sampling_temperature_per_decade=2000,
  sampling_nTOp_per_decade=125, T_start_cosmo=100 MeV) and pin them to the
  tight published tolerances (YP +/-1e-5, D/H +/-3e-9).  These take ~60 s
  total and are the real guard for changes to the nuclear network.

All numbers come from tests/reference_values.py, which is the single source
shared with tests/README.md, test_cli.py, test_gui.py and test_runfiles.py.
"""
import pytest

# Single source for the reference observables (also parsed by
# tests/test_docs_consistency.py and quoted in tests/README.md's "Validation
# reference" tables) -- keep all three in sync via tests/reference_values.py.
from tests.reference_values import (REF_SMALL_YPBBN, REF_SMALL_DOH,
                                    REF_LARGE8_YPBBN, REF_LARGE8_DOH,
                                    NUCLIDE_REFERENCE, NUCLIDE_COLUMNS,
                                    NUCLIDE_REL_TOL)

pytestmark = pytest.mark.slow


# ---------------------------------------------------------------------------
# Default-precision sanity checks (loose)
# ---------------------------------------------------------------------------


@pytest.mark.solve
def test_small_network_YPBBN(solved_small):
    """Gross-regression guard on the small network's YP (loose, abs=1e-4)."""
    assert solved_small.results["YPBBN"] == pytest.approx(0.2469971, abs=1e-4)


@pytest.mark.solve
def test_small_network_DoH(solved_small):
    """Gross-regression guard on the small network's D/H (loose, rel=2e-3)."""
    assert solved_small.results["DoH"] == pytest.approx(2.43590e-5, rel=2e-3)


@pytest.mark.solve
def test_large_network_YPBBN(solved_large):
    """Gross-regression guard on the full large network's YP (loose, abs=1e-4)."""
    assert solved_large.results["YPBBN"] == pytest.approx(0.2470005, abs=1e-4)


@pytest.mark.solve
def test_large_network_DoH(solved_large):
    """Gross-regression guard on the full large network's D/H (loose, rel=2e-3)."""
    assert solved_large.results["DoH"] == pytest.approx(2.43658e-5, rel=2e-3)


# ---------------------------------------------------------------------------
# Per-nuclide table (the live half of its guard)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def solved_large8():
    """A solved PRIMAT instance (``large``, ``amax=8``) reused across the session.

    conftest.py already provides ``solved_small`` and ``solved_large``; this
    completes the trio of networks tests/README.md's per-nuclide table covers,
    so that table's check costs one extra solve rather than three.
    """
    from primat.main import PRIMAT
    r = PRIMAT({"network": "large", "amax": 8, "verbose": False})
    r.solve()
    return r


@pytest.mark.solve
@pytest.mark.parametrize("column,fixture_name", list(zip(
    NUCLIDE_COLUMNS, ("solved_small", "solved_large8", "solved_large"))))
def test_per_nuclide_abundances_match_the_reference_table(
        column, fixture_name, request):
    """Every cell of the published per-nuclide table must survive a real solve.

    GOAL: pin tests/README.md's "Per-nuclide final abundances" table against
    actual abundances. Only ``p`` and ``He4`` used to be pinned anywhere (via
    reference_values' scalars, at an ``abs=1e-4`` that is 1.6e-3 *relative* for
    He4 and vacuous for a ~4e-16 species), so 19 of the 21 cells were free to
    drift -- and by 2026-08-05 every row was stale in its 5th significant
    figure.

    ``n`` matters most here: a regression of the reverse-rate clamp in
    ``primat/network_data.py`` inflates it by a factor ~1750 (to ~7e-13), which
    this catches at ``NUCLIDE_REL_TOL`` = 1e-4 with ~4 orders of margin.
    ``tests/test_large_network.py`` guards the same property from the other
    direction (large vs. large+amax=8).

    The table was snapshotted on the auto (C) backend while these fixtures use
    ``primat.main.PRIMAT`` (pure Python), deliberately: the check must hold
    whether or not the C extension is built. The two backends agree on these
    abundances to <=2.2e-05, comfortably inside the 1e-4 bound -- which is
    itself why the table is documented as good to 5 significant figures, not
    the 7 it prints.
    """
    Y = request.getfixturevalue(fixture_name).nuclear.Y_final

    idx = NUCLIDE_COLUMNS.index(column)
    for nuclide, expected in NUCLIDE_REFERENCE.items():
        assert Y[nuclide] == pytest.approx(expected[idx], rel=NUCLIDE_REL_TOL), (
            f"{nuclide} in the {column!r} column: solved {Y[nuclide]:.7e}, "
            f"table says {expected[idx]:.7e}")


@pytest.mark.solve
def test_Neff_close_to_standard(solved_small):
    """Neff should be close to 3.044 for the standard model."""
    assert solved_small.results["Neff"] == pytest.approx(3.044, abs=0.005)


@pytest.mark.solve
def test_Born_mode_lowers_YP(solved_small):
    """Born-only n<->p rates (radiative/finite-mass corrections off) give lower YP."""
    from primat.main import PRIMAT
    r_born = PRIMAT({"radiative_corrections": False,
                   "finite_mass_corrections": False,
                   "network": "small"})
    r_born.solve()
    assert r_born.results["YPBBN"] < solved_small.results["YPBBN"] - 0.001


@pytest.mark.solve
def test_thermal_corrections_lower_YP(solved_small):
    """CCRTh (finite-temperature radiative corrections) must lower YP, by the
    measured amount.

    GOAL: give ``thermal_corrections`` the physical-effect test every other
    n<->p correction flag already had (``radiative_corrections`` /
    ``finite_mass_corrections`` via ``test_Born_mode_lowers_YP`` above,
    ``spectral_distortions`` via tests/test_spectral_distortions.py,
    ``tau_n_normalization`` via tests/test_mc.py). CCRTh is the most intricate
    of them -- a two-dimensional bremsstrahlung integral evaluated by vegas --
    yet nothing asserted that switching it on changed anything at all: the
    four modules that mention the flag all set it to ``False``, to avoid it.
    Its cached table's *contents* are deliberately excluded from
    tests/test_cache_parity.py too (independent Monte-Carlo streams), so only
    its fingerprint was pinned.

    Measured 2026-08-05 at numerical_precision=1e-8, network="small":

        YPBBN  0.2470107713 (off) -> 0.2469986676 (on)   -1.2104e-05 abs
        D/H    2.435933632e-05    -> 2.435868102e-05     -2.6901e-05 rel

    The sign is the physics (Brown & Sawyer 2001: the thermal photon bath
    slightly suppresses the n->p rate less than p->n, lowering the freeze-out
    n/p ratio); the +/-20% band catches both "silently not wired up"
    (difference exactly 0) and a wrong-magnitude regression, while tolerating
    ordinary vegas/solver noise -- the shift is 1.2x the +-1e-5 YP reference
    tolerance, so it cannot be pinned much more tightly at this precision.

    Read-only: save_nTOp*/=False keeps the fingerprint-mismatched recompute
    from writing a new cache file next to the shipped ones.
    """
    from primat.main import PRIMAT
    common = dict(network="small", numerical_precision=1e-8, verbose=False,
                  debug=False, save_nTOp=False, save_nTOp_thermal=False)
    on = PRIMAT(dict(common, thermal_corrections=True)).primat_results()
    off = PRIMAT(dict(common, thermal_corrections=False)).primat_results()

    assert on["YPBBN"] < off["YPBBN"]
    assert on["YPBBN"] - off["YPBBN"] == pytest.approx(-1.2104e-05, rel=0.2)
    assert (on["DoH"] - off["DoH"]) / off["DoH"] == pytest.approx(-2.6901e-05,
                                                                  rel=0.2)


@pytest.mark.solve
def test_Li7oH_order_of_magnitude(solved_small):
    """Li7/H should be in the range 1e-10 to 1e-9."""
    Li7 = solved_small.results["Li7oH"]
    assert 1e-10 < Li7 < 1e-9


@pytest.mark.solve
def test_He3oH_order_of_magnitude(solved_small):
    """He3/H should be in the range 1e-6 to 1e-4."""
    He3 = solved_small.results["He3oH"]
    assert 1e-6 < He3 < 1e-4


# ---------------------------------------------------------------------------
# High-precision reference checks (tight) — reproduce the published numbers
# ---------------------------------------------------------------------------
# Settings used to produce the published reference values: every solver-facing
# knob runfiles/primat_reference_run.py sets, so this tier checks the
# configuration tests/README.md actually documents. The last four were once
# missing here, which silently made this tier a different run from the
# published one -- worth 2.0e-08 in large+amax=8's D/H, 6.6x the bound below.
# Both weak-rate caches for this configuration are shipped, so it stays fast.
_REF_PARAMS = dict(numerical_precision=1e-10, sampling_temperature_per_decade=2000,
                   sampling_nTOp_per_decade=125, T_start_cosmo_MeV=100.0,
                   Omegabh2=0.02242, verbose=False, debug=False,
                   rate_grid_npts=4000, sampling_nTOp_thermal_per_decade=25,
                   vegas_n_eval=100000, vegas_n_itn=50)


@pytest.fixture(scope="session")
def ref_small():
    from primat.main import PRIMAT
    return PRIMAT({**_REF_PARAMS, "network": "small"}).primat_results()


@pytest.fixture(scope="session")
def ref_large():
    from primat.main import PRIMAT
    return PRIMAT({**_REF_PARAMS, "network": "large", "amax": 8}).primat_results()


@pytest.mark.reference
def test_reference_small_YPBBN(ref_small):
    """High-precision small-network YP reproduces the published value (±1e-5)."""
    assert ref_small["YPBBN"] == pytest.approx(REF_SMALL_YPBBN, abs=1e-5)


@pytest.mark.reference
def test_reference_small_DoH(ref_small):
    """High-precision small-network D/H reproduces the published value (±3e-9).

    This is the tightest pin in the suite -- ~1.2e-4 *relative* -- and the
    reason the reference tier exists: at the routine 1e-7 precision the same
    solve carries ~1e-8 of adaptive-step jitter, which would swamp it."""
    assert ref_small["DoH"] == pytest.approx(REF_SMALL_DOH, abs=3e-9)


@pytest.mark.reference
def test_reference_large_YPBBN(ref_large):
    """High-precision large/amax=8 YP reproduces the published value (±1e-5)."""
    assert ref_large["YPBBN"] == pytest.approx(REF_LARGE8_YPBBN, abs=1e-5)


@pytest.mark.reference
def test_reference_large_DoH(ref_large):
    """High-precision large/amax=8 D/H reproduces the published value (±3e-9)."""
    assert ref_large["DoH"] == pytest.approx(REF_LARGE8_DOH, abs=3e-9)


# ---------------------------------------------------------------------------
# No-numba full solve: pure-Python kernels must agree with the JIT path
# ---------------------------------------------------------------------------

@pytest.mark.solve
def test_no_numba_small_matches_numba(solved_small):
    """Pure-Python (use_numba=False) must agree with the JIT path to 1e-4."""
    from primat.main import PRIMAT
    r_nn = PRIMAT({"use_numba": False, "network": "small"}).primat_results()
    assert r_nn["YPBBN"] == pytest.approx(solved_small.results["YPBBN"], rel=1e-4)
    assert r_nn["DoH"]   == pytest.approx(solved_small.results["DoH"],   rel=1e-4)


@pytest.mark.solve
def test_no_numba_large_amax8_smoke():
    """Pure-Python large/amax=8 network solve completes and YP is physically
    reasonable (the old "medium" network's exact 68-reaction equivalent)."""
    from primat.main import PRIMAT
    r = PRIMAT({"use_numba": False, "network": "large", "amax": 8}).primat_results()
    assert 0.24 < r["YPBBN"] < 0.25
    assert 2.0e-5 < r["DoH"] < 3.0e-5


# ---------------------------------------------------------------------------
# amax cutoff: large network filtered to A <= 20 matches the full large
# network to ~1e-3
# ---------------------------------------------------------------------------

@pytest.mark.solve
def test_amax_filter_light_elements_match_large(solved_large):
    """With amax=20, heavy reactions (A>20) are dropped; light elements match
    the full large network."""
    from primat.main import PRIMAT
    r = PRIMAT({"network": "large", "amax": 20}).primat_results()
    # Light elements should still match the full large-network result to
    # ~1e-3 relative.
    assert r["YPBBN"] == pytest.approx(solved_large.results["YPBBN"], rel=1e-3)
    assert r["DoH"]   == pytest.approx(solved_large.results["DoH"],   rel=1e-3)


@pytest.mark.solve
def test_small_amax2_collapses_to_deuterium_channel():
    """``network="small", amax=2`` must collapse both MT and LT to just the
    n<->p weak rate + n_p__d_g.

    Regression guard for the MT-branch amax-ordering fix: the MT-era
    intersection used to be taken over the *unfiltered* bare reaction names,
    so an amax-violating reaction could still run in the MT era even though
    the LT era correctly dropped it."""
    from primat.main import PRIMAT
    from primat.config import PRIMATConfig
    from primat.network_data import load_network
    cfg = PRIMATConfig({"network": "small", "amax": 2})
    mt_names = load_network(cfg, era="MT").names
    lt_names = load_network(cfg, era="LT").names
    assert mt_names == ["n__p", "n_p__d_g"]
    assert lt_names == ["n__p", "n_p__d_g"]

    r = PRIMAT({"network": "small", "amax": 2}).primat_results()
    assert r["YPBBN"] == 0.0
    assert r["DoH"] > 0.0


def test_every_forward_rate_is_extrapolated_past_the_top_of_its_table():
    """The MT era starts above the master T9 grid, so every run extrapolates.

    `rate_grid_T9_max = 10` GK sits below the mid-temperature era's start at
    T_weak = 11.6045 GK, and the rate buffer continues the last cell's slope
    linearly rather than refusing to answer. The distance is what makes that
    matter, so it is pinned here: widening the grid, moving the era boundary or
    changing `rate_grid_npts` all move it, and none of them would otherwise
    show up anywhere. What it costs is in `docs/performance.md`; the guard that
    the extrapolation cannot go negative is in `tests/test_invariants.py`.
    """
    import numpy as np
    from primat.config import PRIMATConfig
    cfg = PRIMATConfig(params={"network": "small"})
    grid = np.logspace(np.log10(cfg.rate_grid_T9_min),
                       np.log10(cfg.rate_grid_T9_max), cfg.rate_grid_npts)
    T9_weak = cfg.T_weak / 1e9          # cfg.T_weak is in Kelvin, the grid in GK
    assert T9_weak > grid[-1], (
        "the master grid now covers the MT era's start, so nothing is "
        "extrapolated -- this test and the performance note it points at are "
        "both stale")
    cells = (T9_weak - grid[-1]) / (grid[-1] - grid[-2])
    assert cells == pytest.approx(17.48, rel=0.02), (
        f"the forward rates are extrapolated {cells:.2f} cells past the end of "
        "every table, not 17.48; re-measure what that costs before re-pinning")


@pytest.mark.slow
@pytest.mark.solve
def test_default_rate_grid_leaves_a_known_error_in_Li7():
    """What `rate_grid_npts = 1000` costs, once the ODE tolerance is converged.

    The master T9 grid every rate table is resampled onto is a second-order
    accuracy floor sitting under `numerical_precision`: refining it 4x still
    moves `Li7/H` by ~1e-4 and `D/H` by ~6e-6, where one more decade of
    tolerance moves them by ~5e-9. This is the number `docs/performance.md`'s
    "What the default grids cost" quotes, and the reason `Li7/H`'s last two
    reported decimals are grid artefacts. Both backends carry it identically,
    so no parity test can see it.
    """
    from primat.backend import run_bbn
    base = {"network": "small", "numerical_precision": 1e-10}
    coarse = run_bbn(dict(base, rate_grid_npts=1000))
    fine = run_bbn(dict(base, rate_grid_npts=4000))
    for key, expected in (("DoH", 5.9e-6), ("Li7oH", 9.1e-5)):
        shift = abs(fine[key] / coarse[key] - 1.0)
        assert shift == pytest.approx(expected, rel=0.5), f"{key}: {shift:.3e}"
    # YPBBN is unaffected at the level anything is pinned to.
    assert abs(fine["YPBBN"] / coarse["YPBBN"] - 1.0) < 1e-7
