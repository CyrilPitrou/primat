"""Tests for the PRIMAT public API."""
import pytest
import numpy as np
from primat.main import PRIMAT


def test_A_N_Z_dicts():
    """A/N/Z expose each nuclide's mass, neutron and proton number."""
    r = PRIMAT()
    assert r.A["He4"] == 4
    assert r.Z["He4"] == 2
    assert r.N["He4"] == 2
    assert r.A["H2"] == 2
    assert r.A["Li7"] == 7
    assert r.A["n"] == 1
    assert r.Z["n"] == 0


def test_getitem_returns_callable(solved_small):
    """``primat["He4"]`` returns an abundance interpolator, not an array."""
    fn = solved_small["He4"]
    assert callable(fn)


def test_getitem_returns_positive_values(solved_small):
    """The interpolator is array-valued and non-negative above the noise floor."""
    t = np.logspace(0, 5, 20)
    vals = solved_small["He4"](t)
    assert vals.shape == (20,)
    # Before He4 forms its abundance is physically ~0; the stiff BDF solver
    # (and the linear interpolation between its output points) can leave
    # machine-noise-level negative excursions there.  Require non-negativity
    # only above that noise floor (final He4 ~ 0.06).
    assert np.all(vals >= -1e-12)


def test_getitem_scalar_input(solved_small):
    """A scalar time gives a plain positive float, not a 0-d array."""
    val = solved_small["He4"](100.0)
    assert isinstance(val, float)
    assert val > 0


def test_getitem_unknown_species_raises(solved_small):
    """An unknown species name is a KeyError, not a silent zero."""
    with pytest.raises(KeyError):
        solved_small["Unobtainium"]


def test_getitem_all_small_network_species(solved_small):
    """Every small-network nuclide has a usable, non-negative interpolator."""
    for sp in ["n", "p", "H2", "H3", "He3", "He4", "Li7", "Be7"]:
        fn = solved_small[sp]
        assert callable(fn)
        assert fn(100.0) >= 0


def test_T_of_t_and_t_of_T_are_inverses(solved_small):
    """t(T) and T(t) invert each other -- the background's two time coordinates
    must describe one history."""
    T_test = 0.5   # MeV
    t_val = float(solved_small.t_of_T(T_test))
    T_back = float(solved_small.T_of_t(t_val))
    assert T_back == pytest.approx(T_test, rel=1e-4)


def test_a_of_T_and_T_of_a_are_inverses(solved_small):
    """a(T) and T(a) invert each other (see t_of_T/T_of_t above)."""
    T_test = 0.5   # MeV
    a_val = float(solved_small.a_of_T(T_test))
    T_back = float(solved_small.T_of_a(a_val))
    assert T_back == pytest.approx(T_test, rel=1e-4)


def test_a_of_t_and_t_of_a_are_inverses(solved_small):
    """a(t) and t(a) invert each other (see t_of_T/T_of_t above)."""
    t_test = float(solved_small.t_of_T(0.5))   # s, at T_gamma = 0.5 MeV
    a_val = float(solved_small.a_of_t(t_test))
    t_back = float(solved_small.t_of_a(a_val))
    assert t_back == pytest.approx(t_test, rel=1e-4)


def test_a_of_T_consistent_with_a_of_t_and_T_of_t(solved_small):
    """a(T) and a(t(T)) must agree, since both trace the same a(T) solution."""
    T_test = 0.5   # MeV
    t_val = float(solved_small.t_of_T(T_test))
    assert float(solved_small.a_of_t(t_val)) == pytest.approx(
        float(solved_small.a_of_T(T_test)), rel=1e-4)


def test_get_quantity_result_key(solved_small):
    """get_quantity() resolves a result-dict key to that exact value."""
    assert solved_small.get_quantity("YPBBN") == pytest.approx(
        solved_small.results["YPBBN"], rel=1e-12)


def test_get_quantity_nuclide_name(solved_small):
    """get_quantity() also resolves a bare nuclide name to its final abundance,
    so MC/sensitivity callers can name either kind of quantity uniformly."""
    val = solved_small.get_quantity("He4")
    assert val > 0
    assert val == pytest.approx(solved_small.nuclear.Y_final["He4"], rel=1e-12)


def test_get_quantity_unknown_raises(solved_small):
    """A name that is neither a result key nor a nuclide raises ValueError."""
    with pytest.raises(ValueError):
        solved_small.get_quantity("not_a_thing")


def test_lazy_solve_triggers_on_accessor():
    """Accessing a result without calling solve() should auto-trigger it."""
    r = PRIMAT({"network": "small"})
    assert r.results is None
    yp = r.YPBBN()
    assert r.results is not None
    assert yp > 0


def test_solve_cached():
    """Calling solve() twice returns identical results (no re-computation)."""
    r = PRIMAT({"network": "small"})
    res1 = r.solve()
    res2 = r.solve()
    assert res1["YPBBN"] == res2["YPBBN"]


def test_primat_results_returns_dict(solved_small):
    """primat_results() returns the documented result dict, with the headline
    observables present."""
    res = solved_small.primat_results()
    assert isinstance(res, dict)
    for key in ("YPBBN", "YPCMB", "DoH", "He3oH", "Li7oH", "Neff"):
        assert key in res


def test_result_values_physical(solved_small):
    """Every observable lands in its physically sensible range -- the crudest
    possible guard, which catches a catastrophically broken solve."""
    res = solved_small.results
    assert 0.20 < res["YPBBN"] < 0.30
    assert 1e-5 < res["DoH"]   < 5e-5
    assert 1e-6 < res["He3oH"] < 1e-4
    assert 1e-10 < res["Li7oH"] < 1e-9
    assert 2.5 < res["Neff"] < 3.5


# ---------------------------------------------------------------------------
# Solver-failure reporting
#
# scipy's solve_ivp signals a step failure via sol.success/sol.status rather
# than raising, and still returns the partial trajectory in sol.y.  Reading
# sol.y[..., -1] unchecked would silently yield wrong abundances instead of an
# error (an LSODA convergence failure forced on the stiff MT era gives
# YP = 0.434 instead of 0.247, D/H 60% off, with no warning).  The C backend
# already fails loudly (nuclear_network.c checks cpr_ode_bdf's return code),
# so nuclear_network._check_solver keeps the Python backend's error behaviour
# in parity with it.  BDF converges for every supported configuration, so these
# guards never fire in normal use -- hence the failure is injected here.
# ---------------------------------------------------------------------------

def _fake_failed_sol(n_species):
    """A minimal stand-in for a failed scipy OdeResult (success=False)."""
    class _Sol:
        success = False
        status = -1
        message = "Required step size is less than spacing between numbers."
        y = np.ones((n_species, 2))
    return _Sol()


def test_check_solver_passes_on_success():
    """A converged solve must not raise."""
    from primat.nuclear_network import _check_solver

    class _OK:
        success = True
        message = "The solver successfully reached the end of the integration interval."

    _check_solver(_OK(), "LT", "small network, 8 nuclides")   # must not raise


def test_check_solver_raises_on_failure():
    """A failed solve must raise RuntimeError quoting the era and scipy's message."""
    from primat.nuclear_network import _check_solver

    with pytest.raises(RuntimeError) as exc:
        _check_solver(_fake_failed_sol(8), "MT", "small network, 8 species")
    msg = str(exc.value)
    assert "[MT]" in msg                       # era is identified
    assert "small network, 8 species" in msg   # run context is carried
    assert "Required step size" in msg         # scipy's own diagnosis is quoted


@pytest.mark.parametrize("failing_era, tag", [(0, "[HT]"), (1, "[MT]"), (2, "[LT]")])
def test_solve_raises_when_integrator_fails(monkeypatch, failing_era, tag):
    """Each of the three eras must surface an integrator failure as RuntimeError.

    Rather than returning fabricated abundances: the n-th solve_ivp call is made
    to report failure, and PRIMAT.solve() must raise instead of completing.
    """
    import primat.nuclear_network as nn
    real_solve_ivp = nn.solve_ivp
    calls = {"n": 0}

    def flaky(fun, tspan, y0, **kw):
        i = calls["n"]
        calls["n"] += 1
        if i == failing_era:
            return _fake_failed_sol(len(y0))
        return real_solve_ivp(fun, tspan, y0, **kw)

    monkeypatch.setattr(nn, "solve_ivp", flaky)
    with pytest.raises(RuntimeError) as exc:
        PRIMAT().solve()
    assert tag in str(exc.value)


def test_banner_falls_back_to_ascii_on_a_legacy_console(monkeypatch):
    """The startup banner must not crash a console that cannot encode it.

    GOAL: guard the Windows default console encoding (cp1252), which has no
    box-drawing or block-element characters. ``print(_banner())`` raised
    ``UnicodeEncodeError`` there and aborted every verbose run, including
    ``runfiles/primat_run.py``.

    The banner is written to stderr, so stderr's codec is what must decide
    which rendering it gets. Redirecting stdout alone must not change it.
    """
    import io
    import sys

    from primat.main import _banner, console_encodable

    def wrapper(encoding):
        return io.TextIOWrapper(io.BytesIO(), encoding=encoding)

    legacy = wrapper("cp1252")
    monkeypatch.setattr(sys, "stderr", legacy)
    monkeypatch.setattr(sys, "stdout", wrapper("utf-8"))
    assert not console_encodable("┏", sys.stderr)
    banner = _banner()
    banner.encode("cp1252")            # the point: this must not raise
    assert "PRIMAT" in banner

    monkeypatch.setattr(sys, "stderr", wrapper("utf-8"))
    monkeypatch.setattr(sys, "stdout", wrapper("cp1252"))
    assert console_encodable("┏", sys.stderr)
    assert "┏" in _banner()


def test_console_encodable_asks_the_stream_it_is_given(monkeypatch):
    """``console_encodable`` must judge the stream named, not always stdout.

    GOAL: stdout and stderr are redirected independently, so one stream's codec
    says nothing about the other's. Before this was fixed, the check always
    read ``sys.stdout``, which left ``configure_console`` reconfiguring stderr
    on the strength of stdout's encoding -- so a Windows run with stdout piped
    to a file and stderr on a cp1252 console kept the unprotected stderr the
    function exists to protect.
    """
    import io
    import sys

    from primat.main import console_encodable

    utf8 = io.TextIOWrapper(io.BytesIO(), encoding="utf-8")
    legacy = io.TextIOWrapper(io.BytesIO(), encoding="cp1252")

    assert console_encodable("┏", utf8)
    assert not console_encodable("┏", legacy)

    # Omitting the stream still means stdout, the documented default.
    monkeypatch.setattr(sys, "stdout", legacy)
    assert not console_encodable("┏")
    monkeypatch.setattr(sys, "stdout", utf8)
    assert console_encodable("┏")
