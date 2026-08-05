"""Tests for plasma thermodynamics functions."""
import pytest
import numpy as np
from primat.config import PRIMATConfig
from primat.plasma import Plasma, rho_g, rho_nu


@pytest.fixture(scope="module")
def thermo():
    return Plasma(PRIMATConfig())


@pytest.mark.parametrize("T", [0.1, 1.0, 10.0])
def test_rho_g_positive_and_scales(T):
    """The photon energy density is positive and scales as T^4."""
    assert rho_g(T) > 0
    assert rho_g(2 * T) == pytest.approx(rho_g(T) * 16, rel=1e-6)


@pytest.mark.parametrize("T", [0.5, 1.0, 5.0])
def test_rho_e_positive(thermo, T):
    """The e+/e- energy density is positive wherever the pairs are present."""
    assert thermo.rho_e(T) > 0


def test_rho_e_vanishes_at_low_T(thermo):
    """Below the e+/e- cutoff the pair density is exactly 0, not a tiny
    Boltzmann residue -- the cutoff is what keeps the integrands finite."""
    assert thermo.rho_e(1e-5) == 0.0


@pytest.mark.parametrize("T", [0.5, 1.0, 5.0])
def test_p_e_positive(thermo, T):
    """The e+/e- pressure is positive wherever the pairs are present."""
    assert thermo.p_e(T) > 0


def test_spl_and_dspl_dT_consistent_with_standalone(thermo):
    """spl_and_dspl_dT must return the same values as spl and dspl_dT separately."""
    for T in [0.2, 0.5, 1.0, 5.0]:
        s_combined, ds_combined = thermo.spl_and_dspl_dT(T)
        assert s_combined  == pytest.approx(thermo.spl(T),     rel=1e-10)
        assert ds_combined == pytest.approx(thermo.dspl_dT(T), rel=1e-10)


def test_dspl_dT_finite_difference(thermo):
    """dspl_dT should agree with a finite-difference estimate of d(spl)/dT."""
    T = 1.0
    dT = 1e-4
    fd = (thermo.spl(T + dT) - thermo.spl(T - dT)) / (2 * dT)
    assert thermo.dspl_dT(T) == pytest.approx(fd, rel=1e-4)


def test_T_nu_decoupling_high_T_limit(thermo):
    """At high T >> me, entropy is dominated by photons+e±, so T_nu → T_γ."""
    T = 100.0
    assert thermo.T_nu_decoupling(T) == pytest.approx(T, rel=1e-3)


def test_T_nu_decoupling_low_T_limit(thermo):
    """At low T << me, only photon entropy survives, so T_nu → T_γ*(4/11)^(1/3)."""
    T = 0.001
    expected = T * (4.0 / 11.0) ** (1.0 / 3.0)
    assert thermo.T_nu_decoupling(T) == pytest.approx(expected, rel=1e-3)


def test_rho_nu_scaling():
    """rho_nu should scale as T^4."""
    T = 2.0
    assert rho_nu(2 * T) == pytest.approx(rho_nu(T) * 16, rel=1e-6)


def test_spl_positive(thermo):
    """The plasma entropy density is positive across the BBN temperature range."""
    for T in [0.1, 1.0, 10.0]:
        assert thermo.spl(T) > 0


def _expected_electron_thermo_name(cfg):
    """The electron-thermo cache filename `cfg` should produce.

    Rebuilt here from the documented fingerprint fields rather than imported
    from plasma.py, so the test genuinely pins the naming contract: if a field
    is added to (or dropped from) the fingerprint without updating this helper,
    the tests below fail rather than silently following along.
    """
    from primat.cache_utils import constants_hash, fingerprint_hash
    from primat.plasma import ELECTRON_THERMO_FORMAT_VERSION

    fp_hash = fingerprint_hash({
        "format_version":    ELECTRON_THERMO_FORMAT_VERSION,
        "n_electron_table":  cfg.n_electron_table,
        "T_start_cosmo_MeV": cfg.T_start_cosmo_MeV,
        "constants_hash":    constants_hash(),
    })
    return f"electron_thermo_{fp_hash}.txt"


# Grid sizes no shipped cache was ever built with. The cache_dir overlay falls
# back to the shipped tree on a read miss, so a test that used a *plausible*
# fingerprint could silently hit a shipped (or previously-written) file and
# never exercise the write path at all -- which is exactly how the first
# version of test_electron_thermo_cache_is_hash_named passed or failed
# depending on the order tests ran in. Deliberately odd values guarantee a miss
# and keep these tests order-independent. They only set the table's resolution,
# so the physics is unaffected; small values also keep the recompute quick.
_UNSHIPPED_N_TABLE_A = 307
_UNSHIPPED_N_TABLE_B = 311


def test_electron_thermo_cache_is_hash_named(tmp_path):
    """Each configuration gets its OWN electron_thermo_<hash>.txt.

    The cache carries its fingerprint in the *filename*, not just the header,
    so two configurations coexist instead of evicting one another. This is the
    property that fixes F6.12: before it, any run whose fingerprint differed
    from the file on disk overwrote the shipped, git-tracked copy, so a full
    test-suite run left the working tree dirty.

    Everything is written into a throwaway ``cache_dir``; reads still fall back
    to the shipped tree, so nothing here can touch it.
    """
    from primat.plasma import Plasma

    cfg_a = PRIMATConfig({"cache_dir": str(tmp_path),
                          "n_electron_table": _UNSHIPPED_N_TABLE_A})
    cfg_b = PRIMATConfig({"cache_dir": str(tmp_path),
                          "n_electron_table": _UNSHIPPED_N_TABLE_B})

    name_a = _expected_electron_thermo_name(cfg_a)
    name_b = _expected_electron_thermo_name(cfg_b)
    assert name_a != name_b, "different n_electron_table must give different names"

    Plasma(cfg_a)
    Plasma(cfg_b)

    # Both files exist side by side -- the second run did not evict the first.
    written = sorted(p.name for p in (tmp_path / "plasma").glob("electron_thermo_*.txt"))
    assert written == sorted([name_a, name_b])


def test_electron_thermo_cache_header_matches_filename(tmp_path):
    """The fingerprint in the header agrees with the one in the filename.

    The loader checks both (the header check is the guard against a
    hand-edited or truncated file being trusted on its name alone), so they
    must never disagree.
    """
    from primat.cache_utils import read_cache_fingerprint_hash
    from primat.plasma import Plasma

    cfg = PRIMATConfig({"cache_dir": str(tmp_path),
                        "n_electron_table": _UNSHIPPED_N_TABLE_A})
    Plasma(cfg)

    name = _expected_electron_thermo_name(cfg)
    path = tmp_path / "plasma" / name
    assert path.exists(), f"expected {name} to be written"
    # "electron_thermo_<hash>.txt" -> <hash>
    hash_from_name = name[len("electron_thermo_"):-len(".txt")]
    assert read_cache_fingerprint_hash(str(path)) == hash_from_name


def test_electron_thermo_cache_is_reused_on_second_build(tmp_path):
    """A second Plasma with the same config loads the cache instead of rewriting.

    Guards the hit path: with the hash in the filename it would be easy for a
    naming mismatch between writer and reader to go unnoticed, since every run
    would simply recompute (correct results, silently ~0.7 s slower each time).
    Comparing the file's mtime and bytes across two builds catches that.
    """
    from primat.plasma import Plasma

    cfg = PRIMATConfig({"cache_dir": str(tmp_path),
                        "n_electron_table": _UNSHIPPED_N_TABLE_A})
    Plasma(cfg)
    path = tmp_path / "plasma" / _expected_electron_thermo_name(cfg)
    first_bytes = path.read_bytes()
    first_mtime = path.stat().st_mtime_ns

    Plasma(PRIMATConfig({"cache_dir": str(tmp_path),
                         "n_electron_table": _UNSHIPPED_N_TABLE_A}))
    assert path.stat().st_mtime_ns == first_mtime, "cache was rewritten, so it was not reused"
    assert path.read_bytes() == first_bytes


# ---------------------------------------------------------------------------
# C backend plasma tests (require primat._primat_c to be built)
# ---------------------------------------------------------------------------

from primat.backend import HAS_C_BACKEND, run_bbn

requires_c_backend = pytest.mark.skipif(
    not HAS_C_BACKEND,
    reason="primat._primat_c C extension is not built"
)


@requires_c_backend
@pytest.mark.slow
@pytest.mark.backend
def test_c_backend_plasma_without_cache(tmp_path):
    """C backend can compute electron-thermo tables from scratch (cache miss).

    This verifies the fix for the NaN issue where the C backend's electron
    integrands could return NaN when the adaptive quadrature evaluated them
    slightly below E=x (the lower integration bound), causing sqrt(negative).

    The miss is arranged by asking for an ``n_electron_table`` no shipped file
    was built with, so no cache anywhere in the overlay matches and the tables
    must be computed. That replaces the old delete-and-restore dance on the
    shipped file: with hash-named caches a miss is a matter of choosing a
    fingerprint, not of removing a file, so the test no longer has to mutate
    (and hope to restore) the git-tracked tree. ``cache_dir`` catches the
    resulting write.
    """
    result = run_bbn({"network": "small", "n_electron_table": 401,
                      "cache_dir": str(tmp_path)}, force_backend="c")

    # The miss really happened: the C backend wrote its own cache file.
    written = list((tmp_path / "plasma").glob("electron_thermo_*.txt"))
    assert len(written) == 1, f"expected one freshly written cache, got {written}"

    # Verify we got reasonable results (not NaN or obviously wrong)
    assert np.isfinite(result["YPBBN"])
    assert np.isfinite(result["DoH"])
    assert result["YPBBN"] > 0.24
    assert result["YPBBN"] < 0.25
    assert result["DoH"] > 2e-5
    assert result["DoH"] < 3e-5


@requires_c_backend
@pytest.mark.slow
@pytest.mark.backend
def test_c_backend_plasma_with_cache():
    """C backend can read electron-thermo cache written by Python backend."""
    # Ensure cache exists (should be there from normal usage)
    result = run_bbn({"network": "small"}, force_backend="c")
    
    # Verify we got reasonable results
    assert np.isfinite(result["YPBBN"])
    assert np.isfinite(result["DoH"])
    assert result["YPBBN"] > 0.24
    assert result["YPBBN"] < 0.25
    assert result["DoH"] > 2e-5
    assert result["DoH"] < 3e-5


@requires_c_backend
@pytest.mark.slow
@pytest.mark.backend
def test_c_backend_plasma_recompute(tmp_path):
    """C backend can recompute electron-thermo cache when forced.

    The recompute is still redirected to a throwaway ``cache_dir``, but for a
    weaker reason than it used to be. ``recompute_electron_thermo=True`` writes
    the file named by the *current* fingerprint, which for a default config is
    the shipped one — so without the redirect this would still rewrite a
    git-tracked file (with numerically equivalent contents; the ~1e-4
    Python-vs-C gap that made that dangerous was closed at the root by F6.3,
    and the two now agree to ~1e-11).

    What is gone is the eviction hazard: a *non-default* config no longer
    overwrites anything, because its fingerprint names a different file. Reads
    still fall back to the shipped copies, so nothing else is recomputed here.
    """
    result = run_bbn({"network": "small", "recompute_electron_thermo": True,
                      "cache_dir": str(tmp_path)},
                    force_backend="c")


    # Verify we got reasonable results
    assert np.isfinite(result["YPBBN"])
    assert np.isfinite(result["DoH"])
    assert result["YPBBN"] > 0.24
    assert result["YPBBN"] < 0.25
    assert result["DoH"] > 2e-5
    assert result["DoH"] < 3e-5
