"""Tests for mc_uncertainty, MCResult, and MCQuantityResult."""
import os
import pytest
import numpy as np
from primat.main import mc_uncertainty, MCResult, MCQuantityResult, _mc_run_batch

# Every test in this module runs at least one mc_uncertainty() loop, i.e.
# several full PRIMAT().solve() calls -- squarely in the "solve" tier.
pytestmark = [pytest.mark.slow, pytest.mark.solve]

# mc_uncertainty defaults to n_jobs=-1, which needs joblib -- an *optional*
# dependency (the ``mc``/``recommended`` extras). On a lean core install these
# tests have nothing to test, and must skip rather than fail the suite.
pytest.importorskip("joblib", reason="parallel Monte-Carlo needs joblib "
                                     '(pip install "primat[mc]")')

_BASE = {"network": "small"}
_NUM_MC = 8


@pytest.fixture(scope="module")
def mc_single():
    return mc_uncertainty(_NUM_MC, "YPBBN", params=_BASE, seed=0)


@pytest.fixture(scope="module")
def mc_multi():
    return mc_uncertainty(_NUM_MC, ["YPBBN", "DoH", "Li7oH"],
                          params=_BASE, seed=0)


# --- MCResult structure ---

def test_mc_single_returns_MCResult(mc_single):
    """A single-quantity request still returns a full MCResult container."""
    assert isinstance(mc_single, MCResult)


def test_mc_multi_returns_MCResult(mc_multi):
    """A multi-quantity request returns one MCResult covering them all."""
    assert isinstance(mc_multi, MCResult)


def test_mc_single_has_expected_key(mc_single):
    """Iterating an MCResult yields the quantity names that were requested."""
    assert "YPBBN" in list(mc_single)


def test_mc_multi_has_all_keys(mc_multi):
    """Every requested quantity is present -- none silently dropped."""
    for key in ("YPBBN", "DoH", "Li7oH"):
        assert key in list(mc_multi)


# --- MCQuantityResult attributes ---

def test_central_is_float(mc_single):
    """``central`` (the unvaried solve) is a plain float, not a numpy scalar --
    it crosses into JSON output and the Cobaya wrapper."""
    assert isinstance(mc_single["YPBBN"].central, float)


def test_mean_is_float(mc_single):
    """``mean`` is a plain float (see test_central_is_float)."""
    assert isinstance(mc_single["YPBBN"].mean, float)


def test_std_is_float(mc_single):
    """``std`` is a plain float (see test_central_is_float)."""
    assert isinstance(mc_single["YPBBN"].std, float)


def test_values_shape(mc_single):
    """``values`` holds exactly one entry per requested MC sample."""
    assert mc_single["YPBBN"].values.shape == (_NUM_MC,)


def test_mean_consistent_with_values(mc_single):
    """``mean`` is the mean of ``values`` -- the summary cannot drift from the
    samples it summarises."""
    q = mc_single["YPBBN"]
    assert q.mean == pytest.approx(np.mean(q.values), rel=1e-10)


def test_std_consistent_with_values(mc_single):
    """``std`` is the *sample* (ddof=1) standard deviation of ``values``.

    Pinning ddof is what makes diag(cov) == std**2 hold (the cov()/corr()
    matrices use ddof=1 too); a silent switch to the population estimator
    would break that identity at small N."""
    q = mc_single["YPBBN"]
    # std is the sample (ddof=1) standard deviation -- unbiased and consistent
    # with the ddof=1 cov()/corr() matrices, so diag(cov) == std**2 (see
    # test_cov_diag_equals_std_squared).
    assert q.std == pytest.approx(np.std(q.values, ddof=1), rel=1e-10)


def test_central_close_to_nominal(mc_single):
    """Central value should match a plain solve at nominal rates."""
    assert mc_single["YPBBN"].central == pytest.approx(0.2469, abs=1e-3)


# --- std > 0 (rates actually vary) ---

def test_std_positive(mc_single):
    """Varying the nuclear rates produces a non-zero spread -- guards against
    an MC loop that silently solves the same nominal rates every sample."""
    assert mc_single["YPBBN"].std > 0


def test_std_positive_multi(mc_multi):
    """Every quantity in a multi-quantity run has a non-zero spread."""
    for key in ("YPBBN", "DoH", "Li7oH"):
        assert mc_multi[key].std > 0


# --- Covariance / correlation matrices (MCResult.cov / MCResult.corr) ---
#
# These use a synthetic MCResult built from fixed sample arrays (no BBN solve
# needed) so the linear-algebra invariants can be checked deterministically and
# fast, plus one structural check on the real mc_multi fixture.

def _synthetic_mc(seed=0, n=200):
    """A 3-quantity MCResult with known structure: A ~ N(0,1); B strongly
    correlated with A; C constant (zero variance, to exercise the guard)."""
    rng = np.random.default_rng(seed)
    a = rng.standard_normal(n)
    b = 0.7 * a + 0.3 * rng.standard_normal(n)   # correlated with A
    c = np.full(n, 3.0)                            # identical in every sample
    data = {
        "A": MCQuantityResult(0.0, a),
        "B": MCQuantityResult(0.0, b),
        "C": MCQuantityResult(3.0, c),
    }
    return MCResult(data, seed=seed)


def test_cov_matrix_shape_and_symmetry():
    """cov() is square over the tracked quantities and symmetric."""
    mc = _synthetic_mc()
    C = mc.cov()
    assert C.shape == (3, 3)
    assert np.allclose(C, C.T)


def test_cov_diag_equals_std_squared():
    """diag(cov) == std**2 -- the reason MCQuantityResult.std uses ddof=1 too."""
    mc = _synthetic_mc()
    C = mc.cov()
    for i, q in enumerate(mc.quantity_names()):
        assert C[i, i] == pytest.approx(mc[q].std ** 2, rel=1e-12)


def test_corr_unit_diagonal():
    """corr() has a unit diagonal for every quantity, including the
    zero-variance one (whose 0/0 is defined to 1 on the diagonal)."""
    mc = _synthetic_mc()
    R = mc.corr()
    # Unit diagonal for every quantity, including the zero-variance C.
    assert np.allclose(np.diag(R), 1.0)


def test_corr_matrix_symmetry():
    """corr() is symmetric, NaN entries included (the zero-variance rows)."""
    mc = _synthetic_mc()
    R = mc.corr()
    assert np.allclose(R, R.T, equal_nan=True)


def test_corr_matches_cov_over_std():
    """R[i,j] == C[i,j] / (std_i std_j) for the varying quantities."""
    mc = _synthetic_mc()
    C, R = mc.cov(), mc.corr()
    names = mc.quantity_names()
    ia, ib = names.index("A"), names.index("B")
    expected = C[ia, ib] / (mc["A"].std * mc["B"].std)
    assert R[ia, ib] == pytest.approx(expected, rel=1e-12)
    # Strongly (positively) correlated by construction.
    assert 0.8 < R[ia, ib] < 1.0


def test_scalar_cov_corr_equal_matrix_entry():
    """The two-name scalar forms cov(a,b)/corr(a,b) return exactly the matrix
    entry, so callers can use either without a discrepancy."""
    mc = _synthetic_mc()
    C, R = mc.cov(), mc.corr()
    names = mc.quantity_names()
    ia, ib = names.index("A"), names.index("B")
    assert mc.cov("A", "B") == pytest.approx(C[ia, ib], rel=1e-12)
    assert mc.corr("A", "B") == pytest.approx(R[ia, ib], rel=1e-12)
    # Scalar of a quantity with itself: cov == var == std**2, corr == 1.
    assert mc.cov("A", "A") == pytest.approx(mc["A"].std ** 2, rel=1e-12)
    assert mc.corr("A", "A") == 1.0


def test_zero_variance_guard_no_warning():
    """A quantity identical in every sample gets NaN off-diagonal correlation
    (both matrix and scalar forms) and zero covariance, with NO RuntimeWarning
    storm from the 0/0 division."""
    mc = _synthetic_mc()
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")          # any RuntimeWarning -> failure
        R = mc.corr()
        scalar = mc.corr("A", "C")
    names = mc.quantity_names()
    ic = names.index("C")
    ia = names.index("A")
    assert np.isnan(R[ia, ic]) and np.isnan(R[ic, ia])
    assert R[ic, ic] == 1.0                     # unit diagonal even so
    assert np.isnan(scalar)
    # Covariance of the constant with itself is exactly 0 (ddof=1 variance).
    assert mc.cov("C", "C") == pytest.approx(0.0, abs=1e-30)


def test_cov_corr_unknown_name_raises_keyerror():
    """An unknown quantity name is a KeyError, not a silent NaN."""
    mc = _synthetic_mc()
    with pytest.raises(KeyError):
        mc.cov("not_a_quantity", "A")
    with pytest.raises(KeyError):
        mc.corr("A", "not_a_quantity")


def test_cov_corr_one_name_raises_typeerror():
    """Exactly zero or exactly two names -- one name is a usage error."""
    mc = _synthetic_mc()
    with pytest.raises(TypeError):
        mc.cov("A")
    with pytest.raises(TypeError):
        mc.corr("A")


def test_cov_corr_files_roundtrip():
    """dump_mc_covariance/dump_mc_correlation write the two-header-line TSV;
    reparsing the matrix body reproduces mc.cov()/mc.corr()."""
    from primat.backend import dump_mc_covariance, dump_mc_correlation
    mc = _synthetic_mc()

    def _parse(text):
        lines = text.splitlines()
        assert lines[0].startswith("#")                    # line 1: comment
        header = lines[1].split("\t")
        assert header[0] == "quantity"
        names = header[1:]
        rows = []
        for ln in lines[2:]:
            cells = ln.split("\t")
            assert cells[0] in names                       # row label
            rows.append([float(x) for x in cells[1:]])
        return names, np.array(rows)

    names_c, M_cov = _parse(dump_mc_covariance(mc))
    assert names_c == mc.quantity_names()
    assert np.allclose(M_cov, mc.cov(), equal_nan=True)

    names_r, M_corr = _parse(dump_mc_correlation(mc))
    assert names_r == mc.quantity_names()
    assert np.allclose(M_corr, mc.corr(), equal_nan=True)


def test_mc_multi_cov_shape_matches_quantities(mc_multi):
    """On a real MC result the covariance is square over all quantity_names."""
    C = mc_multi.cov()
    nq = len(mc_multi.quantity_names())
    assert C.shape == (nq, nq)
    assert np.allclose(C, C.T)


# --- Reproducibility ---

def test_same_seed_same_result():
    """A fixed seed reproduces the sample values exactly -- the property every
    MC-based uncertainty or covariance estimate is quoted against."""
    mc_a = mc_uncertainty(4, "YPBBN", params=_BASE, seed=42)
    mc_b = mc_uncertainty(4, "YPBBN", params=_BASE, seed=42)
    np.testing.assert_array_equal(mc_a["YPBBN"].values, mc_b["YPBBN"].values)


def test_different_seed_different_result():
    """Different seeds give different samples (the RNG is actually seeded from
    the argument, not from a constant)."""
    mc_a = mc_uncertainty(4, "YPBBN", params=_BASE, seed=0)
    mc_b = mc_uncertainty(4, "YPBBN", params=_BASE, seed=99)
    assert not np.allclose(mc_a["YPBBN"].values, mc_b["YPBBN"].values)


# --- Incremental reuse (prev=) ---

def test_extend_matches_full_run():
    """Extending an N-sample result to M>N must give *exactly* the same M
    samples as computing M from scratch -- the whole point of the ``prev``
    reuse is that the first N samples are seed-deterministic and untouched."""
    full = mc_uncertainty(6, ["YPBBN", "DoH"], params=_BASE, seed=0)
    part = mc_uncertainty(3, ["YPBBN", "DoH"], params=_BASE, seed=0)
    ext  = mc_uncertainty(6, ["YPBBN", "DoH"], params=_BASE, seed=0,
                          prev=part)
    for q in ("YPBBN", "DoH"):
        np.testing.assert_array_equal(full[q].values, ext[q].values)
        assert full[q].central == ext[q].central


def test_extend_truncates_when_fewer_requested():
    """Requesting fewer samples than ``prev`` truncates without solving."""
    big   = mc_uncertainty(6, "YPBBN", params=_BASE, seed=0)
    small = mc_uncertainty(4, "YPBBN", params=_BASE, seed=0, prev=big)
    np.testing.assert_array_equal(big["YPBBN"].values[:4], small["YPBBN"].values)


def test_truncating_reuse_truncates_nuclides_too():
    """GOAL: a truncating ``prev`` reuse must truncate EVERY column, not just
    the explicitly requested quantities.

    The nuclide block (always merged in, see _DEFAULT_MC_OBSERVABLES) used to
    keep all len(prev) rows while the quantity block was cut to num_mc, which
    (a) reported each nuclide's mean/std over the wrong sample count and
    (b) made samples_array()/cov()/corr() -- hence dump_mc_samples -- raise
    ValueError on the ragged columns. The C path (backend.run_mc) always
    truncated uniformly, so this was also a backend-parity bug."""
    big   = mc_uncertainty(6, "YPBBN", params=_BASE, seed=0)
    small = mc_uncertainty(4, "YPBBN", params=_BASE, seed=0, prev=big)

    # Every column -- observables AND nuclides -- carries exactly num_mc rows.
    assert {len(small[q].values) for q in small} == {4}
    np.testing.assert_array_equal(big["He4"].values[:4], small["He4"].values)
    # The aggregate accessors must work, not raise.
    assert small.samples_array().shape[0] == 4
    assert small.cov().shape[0] == small.samples_array().shape[1]
    # ... and the per-nuclide sigma must be the 4-sample one, not big's.
    assert small["He4"].std == pytest.approx(
        float(np.std(big["He4"].values[:4], ddof=1)), rel=1e-12)


def test_prev_ignored_when_seed_differs():
    """An incompatible ``prev`` (different seed) is silently ignored, giving a
    full recompute at the requested seed rather than reusing stale samples."""
    prev = mc_uncertainty(3, "YPBBN", params=_BASE, seed=0)
    ref  = mc_uncertainty(3, "YPBBN", params=_BASE, seed=5)
    got  = mc_uncertainty(3, "YPBBN", params=_BASE, seed=5, prev=prev)
    np.testing.assert_array_equal(ref["YPBBN"].values, got["YPBBN"].values)


def test_result_records_seed():
    """MCResult.seed is stored so callers (e.g. the GUI) can decide whether a
    cached result is reusable as ``prev``."""
    mc = mc_uncertainty(2, "YPBBN", params=_BASE, seed=7)
    assert mc.seed == 7


# --- nuclide name as quantity ---

def test_nuclide_quantity_works():
    """A bare nuclide name works as an MC quantity, not just a result-dict key."""
    mc = mc_uncertainty(4, "He4", params=_BASE, seed=0)
    assert isinstance(mc, MCResult)
    assert mc["He4"].central > 0
    assert mc["He4"].std > 0


# --- Large network variation ---

def test_mc_large_network_varies_heavy_elements():
    """Verify that MC on the large network varies species only present there."""
    # We choose B10, which is only produced in the large network (or at least
    # its variation depends on large-network-only reactions).
    # Using a tiny sample size for speed.
    mc = mc_uncertainty(4, ["DoH", "B10"], params={"network": "large"}, seed=0)
    assert mc["DoH"].std > 0
    assert mc["B10"].std > 0


# ---------------------------------------------------------------------------
# tau_n variation (Item 14)
# ---------------------------------------------------------------------------

def test_tau_n_alone_gives_nonzero_spread_in_YPBBN():
    """With no nuclear-rate offsets (rate_keys=[]), the only randomness left
    is tau_n_sample = tau_n_central + std_tau_n * randn() (one extra draw per
    sample, see _mc_run_batch).  Since YPBBN depends on the n<->p weak-rate
    normalisation 1/(Fn*tau_n), its spread across samples must be non-zero and
    of plausible magnitude (a fraction of a percent, comparable to the
    rate-driven spread in test_std_positive)."""
    res = np.array(_mc_run_batch({"network": "small", "verbose": False},
                                  rate_keys=[], quantities=["YPBBN"],
                                  seeds=list(range(8))))
    std = res[:, 0].std()
    assert 0 < std < 1e-3


def test_tau_n_normalization_false_disables_tau_n_effect():
    """With cfg.tau_n_normalization=False, tau_n does not enter background.NormWeakRates
    (see StandardBackground._setup_weak_rates), so the extra per-sample tau_n
    draw must be a no-op:
    with no rate offsets either, every sample reproduces the central value."""
    res = np.array(_mc_run_batch(
        {"network": "small", "verbose": False, "tau_n_normalization": False},
        rate_keys=[], quantities=["YPBBN"], seeds=list(range(8))))
    assert np.all(res[:, 0] == res[0, 0])


# ---------------------------------------------------------------------------
# custom_network support in mc_uncertainty / _mc_run_batch
# ---------------------------------------------------------------------------

import primat
_TABLES_DIR = os.path.join(os.path.dirname(primat.__file__),
                            "data", "nuclear", "tables", "d_d__He3_n")


def _table_text(T9, rate, err):
    """Build a 3-column rate-table text buffer (T9, rate, err), one row per
    sample point -- the format expected by custom_network["replaced"]."""
    lines = [f"{t:.6e} {r:.6e} {e:.6e}" for t, r, e in zip(T9, rate, err)]
    return "\n".join(lines) + "\n"


def test_removed_reaction_changes_central_value():
    """Removing a reaction alters the solved network, so DoH's central value
    (computed with custom_network) must differ from the default-network one."""
    default = mc_uncertainty(2, "DoH", params=_BASE, seed=0)
    removed = mc_uncertainty(2, "DoH", params=_BASE, seed=0,
                              custom_network={"removed": ["d_d__t_p"]})
    assert removed["DoH"].central != pytest.approx(default["DoH"].central)


def test_custom_error_column_drives_spread():
    """The core of the user's question: a reaction's MC spread should track
    its *custom* error column, not the shipped default. Same median rate,
    two different uncertainty factors -- low spread vs high spread."""
    T9, rate, _err = np.loadtxt(
        os.path.join(_TABLES_DIR, "d_d__He3_n_primat.txt"), unpack=True)

    noerr_table  = _table_text(T9, rate, np.full_like(rate, 1.0))
    bigerr_table = _table_text(T9, rate, np.full_like(rate, 3.0))

    base_params = {"network": "small", "verbose": False, "debug": False}
    seeds = list(range(8))

    res_noerr = np.array(_mc_run_batch(
        base_params, rate_keys=["p_d_d__He3_n"], quantities=["DoH"], seeds=seeds,
        custom_network={"replaced": {"d_d__He3_n": noerr_table}}))
    res_bigerr = np.array(_mc_run_batch(
        base_params, rate_keys=["p_d_d__He3_n"], quantities=["DoH"], seeds=seeds,
        custom_network={"replaced": {"d_d__He3_n": bigerr_table}}))

    std_noerr  = res_noerr[:, 0].std()
    std_bigerr = res_bigerr[:, 0].std()
    # expsigma=1 means p_d_d__He3_n no longer perturbs the rate (median *
    # exp(p*log(1)) = median); the residual std_noerr is from the unrelated
    # per-sample tau_n draw (_mc_run_batch), so it should be tiny compared to
    # the big-error case rather than exactly zero.
    assert std_bigerr > 100 * std_noerr


def test_replaced_table_std_via_public_api():
    """Same as above, but through the public mc_uncertainty() entry point,
    proving the custom_network plumbing works end-to-end (not just via the
    internal _mc_run_batch worker)."""
    T9, rate, _err = np.loadtxt(
        os.path.join(_TABLES_DIR, "d_d__He3_n_primat.txt"), unpack=True)
    bigerr_table = _table_text(T9, rate, np.full_like(rate, 5.0))

    default = mc_uncertainty(8, "DoH", params=_BASE, seed=0)
    replaced = mc_uncertainty(
        8, "DoH", params=_BASE, seed=0,
        custom_network={"replaced": {"d_d__He3_n": bigerr_table}})

    assert replaced["DoH"].std > default["DoH"].std


def test_added_reaction_is_varied():
    """GOAL: a brand-new reaction from custom_network["added"] must have its
    rate uncertainty propagated, like every other reaction in the network.

    ``_mc_resolve_rate_keys`` derives the varied set from the network *file*,
    which by construction cannot list an added reaction -- so added reactions
    used to be integrated but never varied, while the C sampler
    (primat-c/src/mc.c, which iterates the solved network) did vary them: the
    two backends propagated different uncertainty sets for the same custom
    network. Cheap check: no solve, just the key resolution + the solved
    network's reaction list."""
    from primat import PRIMAT
    from primat.main import _mc_resolve_rate_keys

    T9, rate, _err = np.loadtxt(
        os.path.join(_TABLES_DIR, "d_d__He3_n_primat.txt"), unpack=True)
    added_name = "t_t__He4_n_n"
    custom = {"added": {added_name: _table_text(T9, rate,
                                                np.full_like(rate, 1.2))}}
    params = {"network": "small", "verbose": False, "debug": False}

    keys = _mc_resolve_rate_keys(params, custom)
    solved = PRIMAT(params=params, custom_network=custom).nucl._lt_net.names[1:]

    assert added_name in solved, "added reaction should reach the solved network"
    assert f"p_{added_name}" in keys
    # Appended last, mirroring UpdateNuclearRates' `_selected_names +=
    # added_names`, so a run without custom_network keeps its RNG stream (and
    # therefore its sample values) unchanged.
    assert keys[-1] == f"p_{added_name}"
    assert keys[:-1] == _mc_resolve_rate_keys(params, None)


def test_prev_ignored_when_custom_network_differs():
    """A prev computed under one custom_network must not be silently reused
    for a different one -- mirrors test_prev_ignored_when_seed_differs."""
    T9, rate, _err = np.loadtxt(
        os.path.join(_TABLES_DIR, "d_d__He3_n_primat.txt"), unpack=True)
    bigerr_table = _table_text(T9, rate, np.full_like(rate, 5.0))
    cn = {"replaced": {"d_d__He3_n": bigerr_table}}

    prev = mc_uncertainty(3, "DoH", params=_BASE, seed=0)
    ref  = mc_uncertainty(3, "DoH", params=_BASE, seed=0, custom_network=cn)
    got  = mc_uncertainty(3, "DoH", params=_BASE, seed=0, custom_network=cn,
                          prev=prev)
    np.testing.assert_array_equal(ref["DoH"].values, got["DoH"].values)


def test_prev_ignored_when_params_differ():
    """A prev computed under different params (here: network) must not be
    silently reused -- closes the pre-existing blind spot in the reuse guard."""
    prev = mc_uncertainty(3, "DoH", params={"network": "small"}, seed=0)
    large_amax8 = {"network": "large", "amax": 8}
    ref  = mc_uncertainty(3, "DoH", params=large_amax8, seed=0)
    got  = mc_uncertainty(3, "DoH", params=large_amax8, seed=0,
                          prev=prev)
    np.testing.assert_array_equal(ref["DoH"].values, got["DoH"].values)
