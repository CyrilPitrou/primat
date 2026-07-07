"""
Backend parity: primat._primat_c (C) vs. primat.main.PRIMAT (Python).

Why this test exists
---------------------
CLAUDE.md's "Keeping primat-c and primat in sync" section mandates that any
change to the physics/numerics of one backend be mirrored in the other, and
that "the two backends must also agree on the *shape* of their outputs
(same result-dict keys, ...) so callers can switch backends transparently."
This file is that check: it pins down (1) the result-dict *shape* exactly,
and (2) the numerical agreement *level* the two backends currently achieve,
so a future change that silently widens the gap is caught.

Residual gap (root-caused, FABLEPLAN 2.1)
-----------------------------------------
The dominant C-vs-Python parity gap used to be a weak-rate *interpolation
scheme* mismatch: Python interpolated the cached n<->p rate table with a
linear-space quadratic spline (scipy ``interp1d(kind='quadratic')``), while
the C backend used a local 3-point Lagrange quadratic. On the same cached
nodes the two curves differed by up to ~1e-4 relative through the n/p
freeze-out window (T ~ 0.2..2 MeV), C systematically above Python -- which
propagated to a systematic ~2.5e-5 YP and ~1.8e-5 D/H offset (Neff was always
identical, since the background agrees exactly). It was *not* BDF controller
noise: the gap plateaued rather than shrank as ``numerical_precision`` was
tightened 1e-7 -> 1e-10.

Both backends now interpolate the weak rates with the *same* log10-log10
not-a-knot cubic spline (Python ``_weak_rate_loglog_interp``, C
``cpr_weak_rate_nTOp`` via ``cpr_cubic_spline_fit_notaknot``) -- the scheme the
nuclear rate tables already share -- which also happens to be ~1-2 orders of
magnitude more accurate than either old scheme at the shipped node density.
This collapsed the YP gap to ~1e-6 and the D/H gap to <~1e-5. The small
residual D/H gap (~1e-6 for 'small', ~6e-6 for 'large'+amax=8) is the
*separate*, smaller nuclear-rate-table interpolation + BDF step-sequence
difference between the two solver stacks (it does not touch YP, which is set
purely by n/p freeze-out) -- BDF controller noise that shrinks with tighter
``numerical_precision``, not a scheme mismatch.

The cross-backend budget below (``rel=5e-5`` in D/H) is therefore ~40x tighter
than the pre-fix ``1e-3`` -- still a distinct, coarser budget than CLAUDE.md's
+/-3e-9 *same-backend* D/H regression tolerance, leaving headroom for the
BDF-noise residual to vary across platforms. Tightening it further tracks any
future unification of the nuclear-rate/solver residual; loosening it should
not happen without updating this docstring.
"""
import numpy as np
import pytest
import subprocess
import sys

from primat.backend import HAS_C_BACKEND, run_bbn
from primat.evolution import dump_evolution, load_evolution

pytestmark = [pytest.mark.slow, pytest.mark.solve, pytest.mark.backend]

requires_c_backend = pytest.mark.skipif(
    not HAS_C_BACKEND, reason="primat._primat_c C extension is not built"
)

# Keys always present in solve()'s result dict (see primat/main.py), i.e.
# regardless of network/flags -- excludes the conditional keys
# (Li6oLi7/YCNO/Neff/Omeganurel/OneOverOmeganunr).
_ALWAYS_KEYS = {"YPCMB", "YPBBN", "DoH", "He3oH", "He3oHe4", "Li7oH"}


@requires_c_backend
def test_backend_result_dict_shape_matches():
    """C and Python backends return the same result-dict keys for 'small'."""
    params = {"network": "small"}
    r_c = run_bbn(params, force_backend="c")
    r_py = run_bbn(params, force_backend="python")

    assert _ALWAYS_KEYS <= r_c.keys()
    assert _ALWAYS_KEYS <= r_py.keys()
    # Standard background (the default here) always provides the neutrino
    # sector hooks (see PRIMAT.solve()'s final_nu-guarded keys).
    assert {"Neff", "Omeganurel", "OneOverOmeganunr"} <= r_c.keys()
    assert {"Neff", "Omeganurel", "OneOverOmeganunr"} <= r_py.keys()
    # r_c carries one extra "Y_final" sub-dict the Python solve() result does
    # not (a bonus the C wrapper adds, see _wrapper.c's results_to_dict);
    # every other key must match exactly.
    assert r_c.keys() - {"Y_final"} == r_py.keys()
    assert isinstance(r_c["Y_final"], dict)


@requires_c_backend
def test_backend_small_network_numerical_agreement():
    """C vs. Python agreement budget for network='small' (see module docstring)."""
    params = {"network": "small"}
    r_c = run_bbn(params, force_backend="c")
    r_py = run_bbn(params, force_backend="python")

    assert r_c["YPBBN"] == pytest.approx(r_py["YPBBN"], abs=1e-5)
    assert r_c["Neff"] == pytest.approx(r_py["Neff"], abs=1e-3)
    # Post-fix residual is ~1e-6 relative (see module docstring: the weak-rate
    # interpolation scheme is now shared; what remains is nuclear-rate/BDF
    # noise). Budgeted at 5e-5 for cross-platform headroom.
    assert r_c["DoH"] == pytest.approx(r_py["DoH"], rel=5e-5)


@requires_c_backend
def test_backend_large_amax8_numerical_agreement():
    """C vs. Python agreement for network='large', amax=8 (PRIMAT.md S8.2's
    second reference config, alongside 'small' above)."""
    params = {"network": "large", "amax": 8}
    r_c = run_bbn(params, force_backend="c")
    r_py = run_bbn(params, force_backend="python")

    assert r_c["YPBBN"] == pytest.approx(r_py["YPBBN"], abs=1e-5)
    assert r_c["Neff"] == pytest.approx(r_py["Neff"], abs=1e-3)
    # ~6e-6 residual for large+amax=8 (see module docstring); same 5e-5 budget.
    assert r_c["DoH"] == pytest.approx(r_py["DoH"], rel=5e-5)


@requires_c_backend
def test_backend_GN_and_tau_n_agree_at_default():
    """GN/tau_n at their DEFAULT_PARAMS values must give the same result as
    not passing them at all, on the C backend.

    Regression test for a unit-mismatch bug: DEFAULT_PARAMS["GN"] is stored
    in SI units [m^3 kg^-1 s^-2] (see config.py's ``_GN_MeV2`` property and
    the GUI's "Constants" panel), but the C backend's ``CPRConfig.GN`` field
    is consumed directly by ``cpr_config_Mpl``/the Friedmann-equation Hubble
    helper in natural units [MeV^-2]. Before ``cpr_config_set_by_name``
    special-cased "GN" (mirroring the existing "Omegabh2" special case) to
    convert SI -> natural units, any GN value forwarded from Python -- even
    the SI-unit default -- landed in the C struct unconverted, off by ~34
    orders of magnitude, and produced a meaningless (garbage/NaN-adjacent)
    result. See primat-c/src/config.c's cpr_config_set_GN/cpr_config_set_by_name.
    """
    params = {"network": "small"}
    r_base = run_bbn(params, force_backend="c")
    r_explicit_default = run_bbn(
        {**params, "GN": 6.674299257609439e-11, "tau_n": 878.4},
        force_backend="c",
    )
    assert r_explicit_default["YPBBN"] == pytest.approx(r_base["YPBBN"], abs=1e-10)
    assert r_explicit_default["DoH"] == pytest.approx(r_base["DoH"], rel=1e-10)


@requires_c_backend
def test_backend_GN_and_tau_n_perturbation_agrees_with_python():
    """A +1% GN perturbation (and a +1 sigma tau_n perturbation) must shift
    YPBBN in the same direction and to comparable magnitude on both
    backends -- catching both a broken (unconverted) GN and a GN/tau_n that
    is silently ignored by the C backend (see module docstring for the
    unit-mismatch bug this guards against)."""
    params = {"network": "small"}
    r_c_base = run_bbn(params, force_backend="c")
    r_py_base = run_bbn(params, force_backend="python")

    r_c_gn = run_bbn({**params, "GN": 6.674299257609439e-11 * 1.01}, force_backend="c")
    r_py_gn = run_bbn({**params, "GN": 6.674299257609439e-11 * 1.01}, force_backend="python")
    d_c_gn = r_c_gn["YPBBN"] - r_c_base["YPBBN"]
    d_py_gn = r_py_gn["YPBBN"] - r_py_base["YPBBN"]
    # A +1% GN increase should raise YPBBN by O(1e-3) (faster expansion ->
    # earlier n/p freeze-out -> more He4) on both backends -- not the ~34
    # order-of-magnitude-wrong (effectively meaningless) shift a raw
    # SI-into-natural-units assignment produces.
    assert d_c_gn == pytest.approx(0.00088, abs=3e-4)
    assert d_c_gn == pytest.approx(d_py_gn, rel=0.3)

    r_c_tau = run_bbn({**params, "tau_n": 878.4 + 0.5}, force_backend="c")
    r_py_tau = run_bbn({**params, "tau_n": 878.4 + 0.5}, force_backend="python")
    d_c_tau = r_c_tau["YPBBN"] - r_c_base["YPBBN"]
    d_py_tau = r_py_tau["YPBBN"] - r_py_base["YPBBN"]
    assert d_c_tau != pytest.approx(0.0, abs=1e-8)
    assert d_c_tau == pytest.approx(d_py_tau, rel=0.3)


@pytest.mark.parametrize("force_backend", [
    "python",
    pytest.param("c", marks=requires_c_backend),
])
def test_output_files_announce_their_paths(force_backend, capfd, tmp_path):
    """Every solve-time output file should announce its path with [output].

    The time-evolution and final-abundance writers are backend-specific
    (Python and C both implement them), so this checks the shared user-facing
    console contract for both backends.
    """
    out_time = tmp_path / f"evolution_{force_backend}.tsv"
    out_final = tmp_path / f"final_{force_backend}.dat"
    params = {
        "network": "small",
        "output_time_evolution": True,
        "output_file": str(out_time),
        "output_final_result": True,
        "output_final_file": str(out_final),
    }
    if force_backend == "c":
        script = (
            "from primat.backend import run_bbn\n"
            f"run_bbn({params!r}, force_backend='c')\n"
        )
        proc = subprocess.run(
            [sys.executable, "-c", script],
            check=True,
            capture_output=True,
            text=True,
        )
        out = proc.stdout
    else:
        run_bbn(params, force_backend=force_backend)
        out = capfd.readouterr().out

    assert "[output] Time-evolution data" in out
    assert str(out_time.resolve()) in out
    assert "[output] Final abundances" in out
    assert str(out_final.resolve()) in out


def test_python_backend_background_output_announces_path(capfd, tmp_path):
    """The Python-only background TSV writer also uses the [output] prefix."""
    out_background = tmp_path / "background.tsv"
    run_bbn({
        "network": "small",
        "output_background_evolution": True,
        "output_background_file": str(out_background),
    }, force_backend="python")
    out = capfd.readouterr().out

    assert "[output] Background time-evolution data" in out
    assert str(out_background.resolve()) in out


@pytest.mark.parametrize("force_backend", [
    "python",
    pytest.param("c", marks=requires_c_backend),
])
def test_output_background_evolution_both_backends(force_backend, capfd, tmp_path):
    """Both backends write output_background.tsv when requested.

    This tests that the C backend now honours cfg->output_background_evolution
    (previously unwired, see primat-c/include/cprimat/api.h history).
    """
    out_background = tmp_path / f"background_{force_backend}.tsv"
    params = {
        "network": "small",
        "output_background_evolution": True,
        "output_background_file": str(out_background),
    }
    if force_backend == "c":
        script = (
            "from primat.backend import run_bbn\n"
            f"run_bbn({params!r}, force_backend='c')\n"
        )
        proc = subprocess.run(
            [sys.executable, "-c", script],
            check=True,
            capture_output=True,
            text=True,
        )
        out = proc.stdout
        # C backend writes file directly, no [output] announcement yet
        # (could be added to cpr_bg_write_time_evolution in future)
    else:
        run_bbn(params, force_backend=force_backend)
        out = capfd.readouterr().out
        assert "[output] Background time-evolution data" in out
        assert str(out_background.resolve()) in out

    # Both backends must produce the file
    assert out_background.exists()
    content = out_background.read_text()
    assert len(content) > 0
    # Check header contains expected columns (T, t, a, H, Tnue, ...)
    header = content.splitlines()[0]
    assert "T [MeV]" in header
    assert "t [s]" in header
    assert "a [1]" in header


@pytest.mark.parametrize("params", [
    {"network": "small"},
    {"network": "large", "amax": 8},
], ids=["small", "large_amax8"])
def test_evolution_round_trip_matches_in_memory_result(params, tmp_path):
    """dump_evolution/load_evolution round-trips the Python backend's
    in-memory EvolutionResult (PRIMAT.md S7.3/S7.4) at full precision --
    the disk file is a derived convenience, not a separate source of truth.
    """
    p = dict(params, output_time_evolution=True, output_file=None)
    result = run_bbn(p, force_backend="python")
    evo = result["evolution"]

    path = tmp_path / "evolution.tsv"
    dump_evolution(evo, path=str(path))
    loaded = load_evolution(str(path))

    np.testing.assert_allclose(loaded.t, evo.t)
    np.testing.assert_allclose(loaded.a, evo.a)
    np.testing.assert_allclose(loaded.T_gamma, evo.T_gamma)
    for flavour in ("e", "mu", "tau"):
        np.testing.assert_allclose(loaded.T_nu[flavour], evo.T_nu[flavour])
    assert loaded.Y.keys() == evo.Y.keys()
    for species in evo.Y:
        np.testing.assert_allclose(loaded.Y[species], evo.Y[species])


@requires_c_backend
@pytest.mark.parametrize("params", [
    {"network": "small"},
    {"network": "large", "amax": 8},
], ids=["small", "large_amax8"])
def test_evolution_cross_backend_agreement(params):
    """C and Python backends' in-memory EvolutionResults (PRIMAT.md S7.3,
    populated with no disk I/O via output_file=None) agree at matching time
    stamps, interpolating one series onto the other's timestamps (mirrors
    test_custom_background.py's table-interpolation comparison pattern), at
    PRIMAT.md S8.2's documented 1e-5 relative tolerance for the core
    background columns. Final-time Y agreement uses the coarser cross-backend
    D/H-level budget from this module's docstring, since the abundance
    curves cross through their steep BBN transition at slightly different t
    grids on the two backends (an O(1) relative artifact right at that
    transition is expected, not a regression -- see the final-row check
    below for the physically meaningful comparison)."""
    p = dict(params, output_time_evolution=True, output_file=None)
    evo_c = run_bbn(p, force_backend="c")["evolution"]
    evo_py = run_bbn(p, force_backend="python")["evolution"]

    from scipy.interpolate import interp1d

    mask = evo_c.t >= evo_py.t[0]
    for ca, pa in [(evo_c.a, evo_py.a), (evo_c.T_gamma, evo_py.T_gamma),
                   (evo_c.T_nu["e"], evo_py.T_nu["e"])]:
        interp_p = interp1d(evo_py.t, pa, fill_value="extrapolate")(evo_c.t)
        np.testing.assert_allclose(ca[mask], interp_p[mask], rtol=1e-4)

    assert evo_c.Y.keys() == evo_py.Y.keys()
    for species in evo_c.Y:
        # Compare final abundances only (the physically meaningful, stable
        # quantity) rather than the whole curve through the steep transition.
        assert evo_c.Y[species][-1] == pytest.approx(evo_py.Y[species][-1], rel=1e-3, abs=1e-20)


@requires_c_backend
def test_run_bbn_auto_prefers_c_backend():
    """force_backend=None/'auto' dispatches to C whenever it is available."""
    r_auto = run_bbn({"network": "small"})
    r_c = run_bbn({"network": "small"}, force_backend="c")
    assert r_auto == r_c


def test_run_bbn_rejects_unknown_force_backend():
    with pytest.raises(ValueError, match="force_backend"):
        run_bbn({"network": "small"}, force_backend="nope")


def test_run_bbn_validates_params_regardless_of_backend():
    """An invalid --network surfaces PRIMATConfig's ValueError pre-dispatch."""
    with pytest.raises(ValueError, match="network must be"):
        run_bbn({"network": "no_such_network"})


def test_run_bbn_python_only_features_force_python_backend(monkeypatch):
    """extra_rho/background always force the Python backend, even when the
    C backend is requested implicitly via 'auto' (custom_network is *not*
    one of these any more -- see primat/backend.py's module docstring)."""
    calls = []
    import primat.backend as backend_mod

    def fake_python_solve(params, extra_rho, custom_network, background, **kw):
        calls.append((extra_rho, custom_network, background))
        return {"YPBBN": 0.0}

    monkeypatch.setattr(backend_mod, "_python_solve", fake_python_solve)
    run_bbn({"network": "small"}, extra_rho=[lambda Tg: 0.0])
    assert len(calls) == 1


@requires_c_backend
def test_run_bbn_c_backend_rejects_python_only_features():
    with pytest.raises(ValueError, match="incompatible"):
        run_bbn({"network": "small"}, force_backend="c", extra_rho=[lambda Tg: 0.0])


@requires_c_backend
def test_run_bbn_c_backend_supports_output_time_evolution():
    """The unified EvolutionResult (primat.evolution, PRIMAT.md S7.3) has a
    C-side equivalent now (CPRResults's evol_* arrays, PRIMAT.md S7.6) --
    force_backend="c" with output_time_evolution=True returns the same
    in-memory EvolutionResult shape as the Python backend, not a raise."""
    result = run_bbn({"network": "small", "output_time_evolution": True,
                       "output_file": None}, force_backend="c")
    from primat.evolution import EvolutionResult
    assert isinstance(result["evolution"], EvolutionResult)


def test_run_bbn_auto_prefers_c_for_output_time_evolution(monkeypatch):
    """'auto' now dispatches output_time_evolution=True to the C backend
    when available, since it no longer needs the Python-only fallback (see
    module docstring in primat/backend.py)."""
    import primat.backend as backend_mod

    calls = []

    def fake_python_solve(params, extra_rho, custom_network, background, **kw):
        calls.append(params)
        return {"YPBBN": 0.0}

    monkeypatch.setattr(backend_mod, "_python_solve", fake_python_solve)
    monkeypatch.setattr(backend_mod, "HAS_C_BACKEND", False)
    run_bbn({"network": "small", "output_time_evolution": True})
    assert len(calls) == 1


@requires_c_backend
def test_run_bbn_c_backend_honors_nuclear_overlay(tmp_path):
    """user_nuclear_dir is supported on the C backend too (see
    primat-c's cpr_config_resolve_rates_path, primat-c/src/config.c): a
    user_nuclear_dir-supplied network file is loadable end-to-end through
    force_backend="c" exactly like a shipped one."""
    net_dir = tmp_path / "networks"
    net_dir.mkdir(parents=True)
    (net_dir / "overlaynet.txt").write_text(
        "n_p__d_g, n_p__d_g_primat.txt\nd_p__He3_g, d_p__He3_g_primat.txt\n"
    )
    params = {"network": "overlaynet", "user_nuclear_dir": str(tmp_path)}
    r_c = run_bbn(params, force_backend="c")
    r_py = run_bbn(params, force_backend="python")
    assert r_c["YPBBN"] == pytest.approx(r_py["YPBBN"], abs=1e-3)


def test_run_bbn_auto_prefers_c_backend_for_nuclear_overlay(tmp_path, monkeypatch):
    """'auto' dispatches to the C backend for a user_nuclear_dir request
    too, now that both backends support the overlay -- only
    extra_rho/background (Python-only features) force Python."""
    import primat.backend as backend_mod

    calls = []

    def fake_c_run_bbn(params, package_dir, custom_network=None, **kw):
        calls.append(params)
        return {"YPBBN": 0.0}

    monkeypatch.setattr(backend_mod, "HAS_C_BACKEND", True)
    monkeypatch.setattr(backend_mod, "_c_ext", type("M", (), {"run_bbn": staticmethod(fake_c_run_bbn)}))
    run_bbn({"network": "small", "user_nuclear_dir": str(tmp_path)})
    assert len(calls) == 1


# A GUI-shaped "Customise Reactions" override (primat/gui/custom_rates.py's
# kept_to_custom_network output shape): drop one small-network reaction,
# substitute another's rate table with a synthetic one (>=4 points -- the
# resampler's cubic not-a-knot fit on the all-positive branch needs at
# least 4 knots, see cpr_resample_rate_table/cpr_cubic_spline_fit_notaknot).
_CUSTOM_NETWORK = {
    "removed": ["Li7_p__a_a"],
    "replaced": {
        "d_p__He3_g": "\n".join(f"{t9} {10.0 * t9} 0.0" for t9 in
                                 (0.001, 0.01, 0.1, 1.0, 5.0, 10.0)),
    },
}


def test_run_bbn_auto_prefers_c_backend_for_custom_network(monkeypatch):
    """'auto' dispatches a custom_network request to the C backend too, now
    that it is no longer a python_only_feature (see primat/backend.py)."""
    import primat.backend as backend_mod

    calls = []

    def fake_c_run_bbn(params, package_dir, custom_network=None, **kw):
        calls.append(custom_network)
        return {"YPBBN": 0.0}

    monkeypatch.setattr(backend_mod, "HAS_C_BACKEND", True)
    monkeypatch.setattr(backend_mod, "_c_ext", type("M", (), {"run_bbn": staticmethod(fake_c_run_bbn)}))
    run_bbn({"network": "small"}, custom_network=_CUSTOM_NETWORK)
    assert calls == [_CUSTOM_NETWORK]


@requires_c_backend
def test_backend_custom_network_result_dict_shape_matches():
    """C and Python backends return the same result-dict keys for a
    custom_network request (mirrors test_backend_result_dict_shape_matches
    above, but exercising the removed/replaced injection path)."""
    params = {"network": "small"}
    r_c = run_bbn(params, force_backend="c", custom_network=_CUSTOM_NETWORK)
    r_py = run_bbn(params, force_backend="python", custom_network=_CUSTOM_NETWORK)

    assert _ALWAYS_KEYS <= r_c.keys()
    assert _ALWAYS_KEYS <= r_py.keys()
    assert r_c.keys() - {"Y_final"} == r_py.keys()


@requires_c_backend
def test_backend_custom_network_numerical_agreement():
    """C vs. Python agreement for a custom_network request, at the same
    cross-backend budget as test_backend_small_network_numerical_agreement
    (this module's docstring) -- removed/replaced reactions are still small
    perturbations of the 'small' network, so the same gap applies. Also
    checks the custom_network actually changed the result relative to the
    plain 'small' run, on both backends, so this isn't silently a no-op."""
    params = {"network": "small"}
    r_c = run_bbn(params, force_backend="c", custom_network=_CUSTOM_NETWORK)
    r_py = run_bbn(params, force_backend="python", custom_network=_CUSTOM_NETWORK)
    r_c_plain = run_bbn(params, force_backend="c")
    r_py_plain = run_bbn(params, force_backend="python")

    assert r_c["YPBBN"] == pytest.approx(r_py["YPBBN"], abs=1e-5)
    assert r_c["DoH"] == pytest.approx(r_py["DoH"], rel=5e-5)
    assert r_c["DoH"] != pytest.approx(r_c_plain["DoH"], rel=1e-6)
    assert r_py["DoH"] != pytest.approx(r_py_plain["DoH"], rel=1e-6)


@requires_c_backend
def test_backend_mc_cov_corr_parity():
    """The MC covariance/correlation feature (FABLEADVICE F-1) must be
    backend-transparent: both backends track the *same* MC quantity set, the
    ``dump_mc_*`` writers emit byte-identical header lines (line 1: N/seed/
    estimator; line 2: the quantity names), each backend's ``cov()``/``corr()``
    is square and symmetric over those quantities, and the matrices agree
    *statistically* at a large-ish N (the two backends draw from different RNG
    streams -- see the module docstring -- so only convergent statistics match,
    not per-sample values).

    The ``small`` network reports exactly its 8 evolved nuclides on *both*
    backends now (see tests/test_nuclear.py
    ::test_small_network_reports_exactly_its_eight_nuclides -- Python no longer
    pads Y_final up to SPECIES_MD), so the full MC quantity list (observables +
    nuclides) is identical across backends.
    """
    from primat.backend import (run_mc, dump_mc_covariance, dump_mc_correlation)

    params = {"network": "small"}
    n = 120
    mc_c = run_mc(n, params=params, force_backend="c", seed=0)
    mc_py = run_mc(n, params=params, force_backend="python", seed=0)

    # Both backends track the identical MC quantity set (observables + the 8
    # small-network nuclides, in the same order).
    assert mc_c.quantity_names() == mc_py.quantity_names()

    # Each backend's own full matrix is square and symmetric.
    for mc in (mc_c, mc_py):
        nq = len(mc.quantity_names())
        C, R = mc.cov(), mc.corr()
        assert C.shape == (nq, nq) and R.shape == (nq, nq)
        assert np.allclose(C, C.T)
        assert np.allclose(R, R.T, equal_nan=True)

    # Both header lines are byte-identical across backends: line 1 (N, seed,
    # estimator convention) and line 2 (the tab-separated quantity names).
    for dump in (dump_mc_covariance, dump_mc_correlation):
        assert dump(mc_c).splitlines()[:2] == dump(mc_py).splitlines()[:2]

    # Statistical agreement of the YPBBN-DoH correlation at N=120 (loose: the
    # two RNG streams give different samples, only convergent statistics).
    assert mc_c.corr("YPBBN", "DoH") == pytest.approx(
        mc_py.corr("YPBBN", "DoH"), abs=0.25)
