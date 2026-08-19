"""Guards for parameter values that are accepted but dangerous.

Round-3 pass 22 moved every ``DEFAULT_PARAMS`` key off its default on both
backends. These are the regressions for what it found: a configuration whose
weak rates come out NaN, a ``GN`` override that must reach the baryon-to-photon
ratio, the bool/int strictness the two parameter setters have to share, and the
two caches that were re-keyed by values that cannot change their contents.
"""
import os
import warnings

import numpy as np
import pytest

from primat.backend import HAS_C_BACKEND, run_bbn
from primat.cache_utils import fingerprint_hash
from primat.config import PRIMATConfig
from primat.weak_rates.api import (_weak_rate_loglog_interp,
                                   validate_weak_rates_finite)
from primat.weak_rates.cache import _thermal_fingerprint, _weak_rate_fingerprint

requires_c_backend = pytest.mark.skipif(
    not HAS_C_BACKEND, reason="primat._primat_c C extension is not built"
)

# Q = mn - mp = 0.355 MeV, below me: the rate integrands' sqrt(E^2 - me^2) has
# no real branch, so every entry of the n<->p table comes out NaN. That
# configuration is now rejected by PRIMATConfig before any integration runs
# (test_config.py's cross-field checks), which is what keeps the C backend's
# adaptive quadrature away from an integrand it cannot converge on. The guards
# below are the second line: a NaN table can still arrive from a cache file
# written before that check existed, or by any other route, and must stop the
# run rather than be reported as an abundance.
_MP_BELOW_ME = 939.2103602481599


def _poison_weak_cache(cache_dir):
    """Write an all-NaN n<->p table under the default config's own cache name.

    Returns the config the file belongs to. The loader keys on the filename,
    so this is exactly what a run that cached a NaN table leaves behind.
    """
    cfg = PRIMATConfig({"cache_dir": str(cache_dir)})
    fname = "nTOp_" + fingerprint_hash(_weak_rate_fingerprint(cfg)) + ".txt"
    weak = os.path.join(str(cache_dir), "weak")
    os.makedirs(weak, exist_ok=True)
    T = np.logspace(np.log10(1.16e7), np.log10(1.16e11), 320)
    with open(os.path.join(weak, fname), "w") as f:
        f.write("# T[K] Gamma_nTOp[1/tau_n] Gamma_pTOn[1/tau_n]\n")
        for t in T:
            f.write(f"{t:.17e} nan nan\n")
    return cfg


def test_non_finite_weak_rates_are_rejected_before_they_are_cached():
    """The validator that stands between a NaN table and the cache: without it
    the table was saved, reloaded by every later run of that configuration, and
    reported as YP = 0.98.

    Driven with a default config and an injected NaN table, because the mass
    combination that used to produce one no longer reaches this point -- see
    _MP_BELOW_ME above. The validator is generic: it judges the table it is
    handed, whatever produced it.
    """
    cfg = PRIMATConfig()
    T = np.logspace(7, 11, 8)
    nan = np.full_like(T, np.nan)
    with pytest.raises(ValueError, match="are not finite"):
        validate_weak_rates_finite(T, nan, nan, cfg, "computed")
    # A healthy table passes, including the zero prefix the p->n rate has.
    ok = np.linspace(1.0, 2.0, 8)
    zeros = np.concatenate([np.zeros(3), ok[3:]])
    validate_weak_rates_finite(T, ok, zeros, cfg, "computed")


def test_a_cached_nan_table_is_rejected_on_the_python_backend(tmp_path):
    """A poisoned cache file must stop the run, not flow into the abundances."""
    _poison_weak_cache(tmp_path)
    with pytest.raises(ValueError, match="are not finite"):
        run_bbn(params={"cache_dir": str(tmp_path), "show_progress": False},
                force_backend="python")


@requires_c_backend
def test_a_cached_nan_table_does_not_crash_the_c_backend(tmp_path):
    """Same file, C backend: this used to be a heap-buffer-overflow read in
    weak_interp_build and a SIGSEGV that killed the host process."""
    _poison_weak_cache(tmp_path)
    with pytest.raises(RuntimeError, match="are not finite"):
        run_bbn(params={"cache_dir": str(tmp_path), "show_progress": False},
                force_backend="c")


def test_an_all_zero_rate_column_is_rejected_not_read_past_the_end():
    """The empty-suffix case behind the same crash: no positive point to
    interpolate, which used to index one past the table on both backends."""
    T = np.logspace(7, 11, 50)
    with pytest.raises(ValueError, match="fewer than two positive"):
        _weak_rate_loglog_interp(T, np.zeros_like(T))


@requires_c_backend
@pytest.mark.slow
@pytest.mark.solve
def test_GN_override_reaches_eta0b_on_both_backends(tmp_path):
    """eta0b goes as 1/G, so a GN override must move D/H by the same amount on
    both backends -- and independently of where GN sits in the params dict."""
    GN, Omegabh2 = 6.6743e-11 * 1.01, 0.02242
    common = {"cache_dir": str(tmp_path), "show_progress": False}
    py = run_bbn(params=dict(common, GN=GN), force_backend="python")
    c_alone = run_bbn(params=dict(common, GN=GN), force_backend="c")
    # Omegabh2 at its own default value, set after GN: this used to be the only
    # spelling that gave C the right answer.
    c_after = run_bbn(params={**common, "GN": GN, "Omegabh2": Omegabh2},
                      force_backend="c")
    assert c_alone["DoH"] == pytest.approx(py["DoH"], rel=5e-5)
    assert c_alone["DoH"] == pytest.approx(c_after["DoH"], rel=1e-9)


def test_bool_and_int_are_not_interchangeable():
    """True is not the integer 1: the Python validator rejects both directions,
    and primat-c's cpr_config_set_by_name must agree (a bool reaching an int
    field made the standalone CLI print a D/H 5.8 % low at exit status 0)."""
    with pytest.raises(TypeError):
        PRIMATConfig({"sampling_nTOp_per_decade": True})
    with pytest.raises(TypeError):
        PRIMATConfig({"verbose": 1})


@requires_c_backend
@pytest.mark.slow
def test_c_cli_rejects_a_bool_for_an_int_parameter():
    """The standalone C CLI is the door with no Python validator behind it."""
    import subprocess
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    binary = os.path.join(root, "primat-c", "build", "primat-c")
    if not os.path.exists(binary):
        pytest.skip("primat-c CLI not built")
    out = subprocess.run([binary, "--set", "sampling_nTOp_per_decade=True"],
                         capture_output=True, text=True, cwd=root)
    assert out.returncode != 0
    assert ("sampling_nTOp_per_decade=True has the wrong type: expected int, "
            "got bool") in (out.stdout + out.stderr)


def test_n_electron_table_needs_four_knots():
    """The electron-thermo tables are fitted with a not-a-knot cubic; 1..3 used
    to pass the range check and die inside the spline fitter."""
    for n in (1, 2, 3):
        with pytest.raises(ValueError, match="not-a-knot"):
            PRIMATConfig({"n_electron_table": n})
    PRIMATConfig({"n_electron_table": 4})


def test_DeltaNeff_below_minus_three_is_rejected():
    """Below -3 the total neutrino energy density is negative; the failure used
    to surface as a NaN initial state from inside the ODE."""
    with pytest.raises(ValueError, match="must be >= -3"):
        PRIMATConfig({"DeltaNeff": -10.0})
    PRIMATConfig({"DeltaNeff": -3.0})


def _hashes(cfg):
    return (fingerprint_hash(_weak_rate_fingerprint(cfg)),
            fingerprint_hash(_thermal_fingerprint(cfg)))


def test_inert_overrides_do_not_re_key_the_caches():
    """Values that cannot change a cached table must not change its hash: the
    analytic-distortion amplitudes in tabulated-distortion mode, and a NEVO
    override naming the file the default already selects."""
    base = _hashes(PRIMATConfig({}))
    for params in ({"y_SZ": 0.01}, {"y_gray": 1e6},
                   {"nevo_file": "NEVOPRIMAT_col_1_7.csv"},
                   {"nevo_spectral_file": "NEVOPRIMAT.csv"},
                   {"nevo_grid_file": "NEVOGrid.csv"}):
        assert _hashes(PRIMATConfig(params)) == base, params


def test_shipped_weak_cache_hash_is_unchanged():
    """The default fingerprint must keep naming the shipped cache file: this is
    what makes the normalisation above free rather than a mass re-key."""
    cfg = PRIMATConfig({})
    fname = "nTOp_" + fingerprint_hash(_weak_rate_fingerprint(cfg)) + ".txt"
    shipped = os.path.join(cfg._resolved_data_dir, "cache_plasma_weak", "weak", fname)
    assert os.path.exists(shipped), fname


def test_analytic_distortion_amplitudes_still_re_key():
    """...while the same amplitudes in the mode that reads them must."""
    analytic = {"analytic_distortions": True, "incomplete_decoupling": False}
    ref = _hashes(PRIMATConfig(analytic))
    assert _hashes(PRIMATConfig(dict(analytic, y_SZ=0.01))) != ref


@pytest.mark.parametrize("params, expected", [
    ({"amax": 2}, "YPBBN = 0"),
    ({"munuOverTnu": 0.05}, "not self-consistent"),
    ({"decay_era": True}, "no effect"),
    ({"numerical_precision": 0.1}, "not converged"),
    ({"sampling_temperature_per_decade": 1}, "not converged"),
])
def test_accepted_but_dangerous_values_warn(params, expected):
    """Each of these runs to completion and reports a number that is not what
    the user is likely to think it is, so each says so."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        PRIMATConfig(params)
    assert any(expected in str(w.message) for w in caught), \
        [str(w.message) for w in caught]


def test_default_config_warns_about_nothing():
    """The guards above must stay silent on a default run."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        PRIMATConfig({})
    assert [str(w.message) for w in caught] == []
