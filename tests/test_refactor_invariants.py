"""
Invariants introduced by the performance/cleanup refactor.

These pin down the behaviour that the refactor relies on for correctness:

* MC results are independent of ``n_jobs`` (the per-worker reuse of the
  background + weak rates must not change the numbers).
* ``eta0b`` is recomputed whenever ``Omegabh2`` is reassigned.
* The electron-thermo tabulation reproduces the exact integrals.
* The fast ``_LinearRate`` evaluator matches ``interp1d(kind='linear')``.
"""
import numpy as np
import pytest

from primat.config import PRIMATConfig


# ---------------------------------------------------------------------------
# 1a — MC reuse must be independent of n_jobs
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.solve
def test_mc_njobs_independence():
    """Same seeds give numerically equivalent samples regardless of n_jobs.

    With n_jobs=1 all seeds run in the main process; with n_jobs>1 joblib
    spawns worker subprocesses whose JIT/floating-point environment may
    differ from the main process by up to the ODE-solver tolerance
    (numerical_precision=1e-7). assert_allclose at rtol=1e-6 verifies
    the samples agree well within that tolerance while tolerating the
    tiny process-environment differences.
    """
    from primat.main import mc_uncertainty
    base = {"network": "small"}
    mc1 = mc_uncertainty(6, ["YPBBN", "DoH"], params=base, n_jobs=1, seed=0)
    mcP = mc_uncertainty(6, ["YPBBN", "DoH"], params=base, n_jobs=3, seed=0)
    np.testing.assert_allclose(mc1["YPBBN"].values, mcP["YPBBN"].values, rtol=1e-6)
    np.testing.assert_allclose(mc1["DoH"].values,   mcP["DoH"].values,   rtol=1e-6)


# ---------------------------------------------------------------------------
# 2e — eta0b tracks Omegabh2
# ---------------------------------------------------------------------------

def test_eta0b_tracks_omegabh2_attribute():
    """Reassigning cfg.Omegabh2 recomputes eta0b, rather than leaving the
    baryon-to-photon ratio at its construction-time value."""
    cfg = PRIMATConfig({"Omegabh2": 0.02242})
    e0 = cfg.eta0b
    cfg.Omegabh2 = 0.024
    assert cfg.eta0b == pytest.approx(e0 * 0.024 / 0.02242, rel=1e-12)


def test_eta0b_tracks_omegabh2_setitem():
    """The same recomputation happens through the cfg[...] item interface."""
    cfg = PRIMATConfig({"Omegabh2": 0.02242})
    e0 = cfg.eta0b
    cfg["Omegabh2"] = 0.024
    assert cfg.eta0b == pytest.approx(e0 * 0.024 / 0.02242, rel=1e-12)


def test_gn_and_taun_come_from_defaults():
    """GN and tau_n must be present and overridable (single source of truth).

    ``cfg.GN`` is stored in SI units [m^3 kg^-1 s^-2]; ``cfg.Mpl`` (the
    natural-units Planck mass used by the Friedmann equation) is derived
    from it via ``CONST.GN_SI_to_MeV2``.
    """
    from primat.constants import CONST
    gn_si = 1.234e-10
    cfg = PRIMATConfig({"GN": gn_si, "tau_n": 880.0})
    assert cfg.GN == gn_si
    assert cfg.tau_n == 880.0
    gn_natural = gn_si * CONST.GN_SI_to_MeV2
    assert cfg.Mpl == pytest.approx(1.0 / np.sqrt(gn_natural), rel=1e-12)


# ---------------------------------------------------------------------------
# §6.1 — pluggable extra energy density (extra_rho)
# ---------------------------------------------------------------------------

def test_extra_rho_is_additive_in_hubble():
    """Each ``extra_rho`` callable adds ``rho(Tg)`` to ``rho_tot`` in
    ``background.Hubble``.

    ``background.Hubble`` returns
    ``H = MeV_to_secm1 * sqrt(rho_tot * 8*pi/(3*Mpl^2))``, so adding a
    constant extra energy density ``extra`` [MeV^4] through the
    ``extra_rho`` plug-in must increase ``H^2`` by exactly
    ``extra * 8*pi/(3*Mpl^2)``, independently of everything else
    `background.Hubble` computes.
    """
    from primat.main import PRIMAT
    base = {"network": "small", "verbose": False}
    p0 = PRIMAT(base)
    extra = 1.e-2  # MeV^4, an arbitrary but sizeable extra radiation density
    p1 = PRIMAT(base, extra_rho=[lambda Tg: extra])

    Tg = 1.0  # MeV
    H0 = p0.background.Hubble(Tg, Tg, Tg, Tg)
    H1 = p1.background.Hubble(Tg, Tg, Tg, Tg)
    assert H1 > H0

    dH2 = (H1 / p0.cfg.MeV_to_secm1)**2 - (H0 / p0.cfg.MeV_to_secm1)**2
    assert dH2 == pytest.approx(extra * 8. * np.pi / (3. * p0.cfg.Mpl**2), rel=1e-12)


def test_ede_is_appended_to_extra_rho():
    """``fEDE > 0`` appends exactly one ``rho_EDE`` callable to
    ``background.extra_rho``, via the same generic plug-in mechanism that
    ``extra_rho=`` callers use.

    Since the ΛCDM setup (``_setup_LCDM``) always pre-populates ``extra_rho``
    with two callables — ``rho_CDM(T)`` and ``rho_Lambda`` — the no-EDE
    baseline has exactly 2 entries.  EDE adds one more, giving 3 total.
    """
    from primat.main import PRIMAT
    p_no_ede = PRIMAT({"network": "small", "verbose": False})
    # 2 ΛCDM entries (CDM + cosmological constant) always present
    assert len(p_no_ede.background.extra_rho) == 2

    p_ede = PRIMAT({"network": "small", "verbose": False, "fEDE": 0.05})
    # EDE appends one more callable on top of the 2 ΛCDM ones
    assert len(p_ede.background.extra_rho) == 3


# ---------------------------------------------------------------------------
# 1b — electron-thermo tabulation reproduces the exact integrals
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("T", [0.05, 0.2, 0.5, 1.0, 5.0])
def test_tabulated_electron_thermo_matches_exact(T):
    """The cubic-interpolant table (always used) reproduces the exact quad
    integrals (``_*_exact``) to within the interpolation tolerance."""
    from primat.plasma import Plasma
    p = Plasma(PRIMATConfig())
    tab   = (p.rho_e(T), p.p_e(T), p.drho_e_dT(T), p.dp_e_dT(T))
    exact = (p._rho_e_exact(T), p._p_e_exact(T),
             p._drho_e_dT_exact(T), p._dp_e_dT_exact(T))
    for e, t in zip(exact, tab):
        assert t == pytest.approx(e, rel=1e-5)


# ---------------------------------------------------------------------------
# speedup — _LinearRate matches interp1d(kind='linear', fill_value='extrapolate')
# ---------------------------------------------------------------------------

def test_linear_rate_matches_interp1d():
    """_LinearRate reproduces scipy interp1d(kind='linear',
    fill_value='extrapolate') to 1e-12, including above the grid.

    _LinearRate exists purely as a faster stand-in for that call in the
    rate-table hot path, so it has to be numerically indistinguishable from
    it -- extrapolation behaviour included, since the master T9 grid is
    queried beyond its top edge."""
    from scipy.interpolate import interp1d
    from primat.network_data import _LinearRate
    rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0.001, 10.0, 50))
    y = rng.uniform(1e-3, 1e3, 50)
    fast = _LinearRate(x, y)
    ref = interp1d(x, y, kind="linear", bounds_error=False,
                   fill_value="extrapolate")
    # over the queried T9 range, including extrapolation above the grid (->11.6)
    xq = np.linspace(0.02, 11.6, 2000)
    np.testing.assert_allclose(fast(xq), ref(xq), rtol=1e-12, atol=0.0)


def test_linear_rate_scalar_and_array():
    """_LinearRate accepts both a scalar and an array, like interp1d."""
    from primat.network_data import _LinearRate
    f = _LinearRate(np.array([1.0, 2.0, 3.0]), np.array([10.0, 20.0, 30.0]))
    assert f(2.0) == pytest.approx(20.0)
    np.testing.assert_allclose(f(np.array([1.5, 2.5])), [15.0, 25.0])


# ---------------------------------------------------------------------------
# speedup — the background's scalar interpolant fast path is bit-exact
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("clamped", [False, True])
def test_scalar_linear_eval_reproduces_scipy_exactly(clamped):
    """_scalar_linear_eval matches interp1d(kind='linear') to the last bit.

    It replaces the scipy call on the T_of_t / a_of_t hot path, so "close
    enough" is not enough: any difference at all would move the BDF step
    sequence and every observable with it. Checked on the nodes, the
    midpoints between them and both extrapolation sides, for both fill_value
    forms the background uses.
    """
    from scipy.interpolate import interp1d
    from primat.background import _scalar_linear_eval
    rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0.0, 10.0, 400))
    y = rng.uniform(1e-3, 1e3, 400)
    # The two forms the background builds: T_of_t extrapolates, a_of_t clamps
    # to its end values.
    fill_value = (y[0], y[-1]) if clamped else "extrapolate"
    ref = interp1d(x, y, kind="linear", bounds_error=False, fill_value=fill_value)
    fast = _scalar_linear_eval(ref)
    assert fast is not None
    probe = np.concatenate([x, 0.5 * (x[:-1] + x[1:]), [-1.0, 11.0]])
    for q in probe:
        assert float(fast(q)) == float(ref(q))


def test_scalar_linear_eval_declines_what_it_cannot_reproduce():
    """_scalar_linear_eval returns None for a non-linear interp1d.

    The fast path is opt-in per interpolant: whatever it cannot reproduce
    exactly must fall back to scipy rather than silently interpolate by a
    different scheme.
    """
    from scipy.interpolate import interp1d
    from primat.background import _scalar_linear_eval
    x = np.linspace(0.0, 10.0, 40)
    cubic = interp1d(x, np.exp(-x), kind="cubic", bounds_error=False,
                     fill_value="extrapolate")
    assert _scalar_linear_eval(cubic) is None


# ---------------------------------------------------------------------------
# speedup — BDF's direct-LAPACK dense LU changes no digit
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.solve
def test_direct_lapack_bdf_matches_scipy_bdf():
    """The MT/LT solves give identical observables through scipy's own BDF.

    nuclear_network._BDF hands BDF getrf/getrs directly instead of going
    through scipy.linalg.lu_factor/lu_solve. Same LAPACK routines, same
    arguments -- so a run forced back onto plain method="BDF" must reproduce
    every digit, not merely agree within a tolerance.
    """
    import primat.nuclear_network as nn
    from primat.backend import run_bbn
    keys = ("Neff", "YPBBN", "DoH", "He3oH", "Li7oH")

    fast = run_bbn(params={"network": "small"}, force_backend="python")
    saved = nn._BDF
    try:
        nn._BDF = "BDF"
        plain = run_bbn(params={"network": "small"}, force_backend="python")
    finally:
        nn._BDF = saved
    assert [fast[k] for k in keys] == [plain[k] for k in keys]
