# -*- coding: utf-8 -*-
"""Cross-backend cache parity: both backends must produce the *same* cache files.

Why this module exists
----------------------
The C and Python backends **share** every on-disk cache: whichever one computes
a table first, the other reads it. That is deliberate, and it is not the same
thing as the two backends merely being compatible.

The alternative -- putting the backend identity into each fingerprint, so the
two never touch each other's files -- was considered and rejected. It converts
a parity *bug* into a parity *blind spot*. The shipped
``electron_thermo_cache.txt`` coming back modified after a test run is precisely
the observation that exposed a 1.0e-4 disagreement between the two backends'
electron thermodynamics, since closed at the root (1.0e-4 -> 8.9e-12). Had the
caches been segregated, that gap would still be in ``plasma.c`` today,
silently, each backend happily reading its own copy.
Segregation would also double the shipped cache tree (~1 MB -> ~2 MB in a wheel
served to Streamlit Community Cloud) and make every user who exercises both
backends pay the multi-minute vegas build twice.

Sharing is therefore the right default -- but it is only *safe* if something
actively checks that the two backends agree. That is this module. It is the
mechanism that makes the shared cache defensible.

What is checked
---------------
For each deterministic cache, both backends are driven with their own
``cache_dir`` on a coarse grid, then:

1. **Hash identity** -- both backends must emit the *same filename*. This pins
   ``cpr_weak_rate_fingerprint`` == ``weak_rate_fingerprint`` and
   ``cpr_constants_hash`` == ``constants_hash(cache)`` field for field, per
   cache -- the two sides declare the same constant subsets. A single
   field present on one side only, a float formatted differently, a key sorted
   differently: all show up here as a different hash.
2. **Column agreement** -- at tolerances held in named module constants whose
   comments carry the measured value and the date measured.

The coarse grid (low ``sampling_nTOp_per_decade``, ``weak_rate_cache=False``)
keeps the module inside the default suite. It does not weaken detection: the
divergence it is guarding against was a per-point quadrature-tolerance floor,
not a grid-resolution effect, so a coarse grid catches it just as well. The grid
density is itself a fingerprint field, so both backends still hash the same
coarse config and the hash-identity assertion is unaffected.

Deliberate gap: the CCRTh thermal cache (``nTOp_thermal_<hash>.txt``) is
**excluded** from the column comparison. Both backends compute it by
Monte-Carlo (vegas) with independent RNG streams, so they agree only to their
own noise floor -- each now reproduces *itself* exactly (both seed
deterministically: ``corrections._vegas_rng`` and ``weak_rates.c``'s
``th_vegas_seed``), but not the other. Forcing a recompute would also cost
minutes per backend, far too slow for the default suite. Its *fingerprint* is
still pinned, indirectly but exactly, by
:func:`test_thermal_cache_fingerprints_agree`. This exclusion is stated here so
the gap is visible rather than accidental.

Runtime: ~13 s. Each of the three module-scoped fixtures drives both backends
once and is shared by the tests that need it; giving every test its own pair of
runs cost ~30 s, which is the point at which this would have had to move behind
an opt-in marker instead of running by default.
"""
import os

import numpy as np
import pytest

from primat.constants import CONST
from primat.backend import HAS_C_BACKEND, run_bbn

pytestmark = pytest.mark.skipif(
    not HAS_C_BACKEND,
    reason="cross-backend cache parity needs the primat._primat_c extension",
)


# ---------------------------------------------------------------------------
# Tolerances. Each carries the value actually measured and when, so a future
# reader can tell a deliberate margin from a number someone made up.
# ---------------------------------------------------------------------------

# nTOp_<hash>.txt, both rate columns, max relative difference.
# Measured 2.5e-10 on the default grid; the pin is ~40x looser to absorb the
# coarse grid used here and platform libm variation.
NTOP_RTOL = 1e-8

# electron_thermo_<hash>.txt, all four columns (rho_e, p_e, drho_e_dT,
# dp_e_dT), max relative difference. Measured 1.3e-10: both backends use a
# quadrature tolerance relative to the integrand's own magnitude.
#
# Applied only above ELECTRON_THERMO_TAIL_FRAC of each column's peak. The
# columns fall ~25 decades into the Boltzmann-suppressed tail, where one
# backend can underflow to exact 0.0 against the other's ~1e-17 and a
# pointwise-relative metric reports 1.0 -- seen in CI on one platform at
# T = 1.70e-02 MeV, 1.1e-17 against a column peak of 5.8e+08. The tail is
# covered by ELECTRON_THERMO_SCALE_RTOL instead, on the same "could this
# difference move anything?" argument as QED_SCALE_RTOL below.
ELECTRON_THERMO_RTOL = 1e-9
ELECTRON_THERMO_TAIL_FRAC = 1e-15
ELECTRON_THERMO_SCALE_RTOL = 1e-11   # measured 1.0e-12

# QED_pressure_correction_e{2,3}.txt: max |difference| normalised by the
# COLUMN'S PEAK MAGNITUDE, not pointwise-relative. Measured 7.3e-20.
#
# Pointwise relative difference is the wrong metric for these two tables, and
# using it would force a meaninglessly loose pin. delta_P falls off the bottom
# of the double-precision range in the Boltzmann-suppressed low-T tail: the
# worst pointwise disagreement, 5.7e-3, sits at T = 1.5e-2 MeV where the value
# is -3.7e-24 MeV^4 against a column peak of order 1e4 MeV^4 -- pure quadrature
# noise at ~1e-28 of the physical scale. 468 of 500 rows already agree to within
# one ulp of the file's own "%.6E" format. Normalising by the column scale
# reports what actually matters: whether the two backends' QED correction
# differs anywhere it could move an observable.
QED_SCALE_RTOL = 1e-15

# Coarse but physically sane run. sampling_nTOp_per_decade is what makes this
# affordable; it is a fingerprint field, so both backends still agree on the
# hash. thermal_corrections=False keeps vegas out of the default suite (see the
# module docstring's "deliberate gap").
_COARSE = {
    "sampling_nTOp_per_decade": 8,
    "sampling_temperature_per_decade": 40,
    # Deliberately a grid size nothing ships a cache for. The cache_dir overlay
    # falls back to the shipped tree on a read miss, so with the DEFAULT
    # n_electron_table both backends would simply load the shipped
    # electron_thermo_<hash>.txt and write nothing -- and the comparison below
    # would be vacuous (or, worse, would compare a file with itself). An
    # unshipped value forces both backends to actually compute the table, which
    # is the only way this module can compare what each one produced. It sets
    # the table's resolution only, so no physics changes; 307 also keeps the
    # recompute quick.
    "n_electron_table": 307,
    "thermal_corrections": False,
    "weak_rate_cache": True,
    "verbose": False,
    "show_progress": False,
    "output_time_evolution": False,
    "output_final_result": False,
    "output_background_evolution": False,
}


def _run_both(tmp_path, extra=None):
    """Run both backends into separate cache_dirs; return the two cache roots.

    Each backend gets its own subtree so nothing is shared *within the test* --
    which is the point: if they agreed only because one read the other's file,
    the comparison would be vacuous.
    """
    dirs = {}
    for backend in ("c", "python"):
        d = tmp_path / backend
        params = dict(_COARSE, cache_dir=str(d))
        if extra:
            params.update(extra)
        run_bbn(params, force_backend=backend)
        dirs[backend] = d
    return dirs


# Each of the three fixtures below drives BOTH backends once and is shared by
# every test that needs that particular pair of runs. Doing it per-test instead
# cost ~4 s x 7 = 30 s; sharing brings the module to ~13 s, which is what keeps
# it in the default suite rather than behind an opt-in marker. Module scope is
# safe here because every test only READS the resulting files.

@pytest.fixture(scope="module")
def coarse_dirs(tmp_path_factory):
    """Both backends on the coarse grid: weak-rate + electron-thermo caches."""
    return _run_both(tmp_path_factory.mktemp("coarse"))


@pytest.fixture(scope="module")
def qed_dirs(tmp_path_factory):
    """Both backends with the QED tables force-recomputed.

    Without ``recompute_qed_corrections`` the shipped pair is a fingerprint hit
    for both backends, so each would read the same file and the comparison
    would be between a file and itself.
    """
    return _run_both(tmp_path_factory.mktemp("qed"),
                     extra={"recompute_qed_corrections": True})


@pytest.fixture(scope="module")
def thermal_dirs(tmp_path_factory):
    """Both backends with the CCRTh thermal correction enabled (cache hit)."""
    return _run_both(tmp_path_factory.mktemp("thermal"),
                     extra={"thermal_corrections": True})


def _names(d, subdir, prefix):
    """Sorted basenames of <d>/<subdir>/<prefix>*.txt."""
    p = d / subdir
    if not p.is_dir():
        return []
    return sorted(f.name for f in p.glob(f"{prefix}*.txt"))


def _max_rel(a, b):
    """Max pointwise relative difference between two equal-shaped arrays."""
    den = np.maximum(np.abs(b), np.finfo(float).tiny)
    return float(np.max(np.abs(a - b) / den))


def _max_scale_rel(a, b):
    """Max |a-b| normalised by b's peak magnitude, column by column.

    The right metric for a column that decays over many orders of magnitude
    (see QED_SCALE_RTOL): it asks "could this difference move anything?"
    rather than "how big is it next to a number that is itself noise?".
    """
    worst = 0.0
    for j in range(b.shape[1]):
        peak = np.max(np.abs(b[:, j]))
        if peak == 0.0:
            assert np.max(np.abs(a[:, j])) == 0.0, "column is zero on one backend only"
            continue
        worst = max(worst, float(np.max(np.abs(a[:, j] - b[:, j])) / peak))
    return worst


def _max_rel_above_scale(a, b, frac):
    """``_max_rel`` restricted to entries above ``frac`` of their column peak.

    A column that decays over ~25 decades ends in a Boltzmann-suppressed tail
    where one backend's quadrature can underflow to exact ``0.0`` while the
    other returns a denormal-scale value. Pointwise-relative scores that as
    1.0 -- a hard failure over a difference that cannot move any observable.
    Callers pair this with ``_max_scale_rel``, which covers the excluded tail.
    """
    worst = 0.0
    for j in range(b.shape[1]):
        peak = np.max(np.abs(b[:, j]))
        keep = np.abs(b[:, j]) > frac * peak
        if peak == 0.0 or not keep.any():
            continue
        worst = max(worst, _max_rel(a[keep, j], b[keep, j]))
    return worst


# ---------------------------------------------------------------------------
# 1. Hash identity -- the assertion that pins the two fingerprint ports.
# ---------------------------------------------------------------------------

def test_weak_rate_cache_filename_identical(coarse_dirs):
    """Both backends name the n<->p weak-rate cache identically.

    A mismatch means cpr_weak_rate_fingerprint and weak_rate_fingerprint have
    drifted -- a field on one side only, a differently formatted float, a
    different key order -- and the two backends would each silently maintain
    their own copy of a cache meant to be shared.
    """
    dirs = coarse_dirs
    c_names = _names(dirs["c"], "weak", "nTOp_")
    py_names = _names(dirs["python"], "weak", "nTOp_")
    assert c_names, "C backend wrote no weak-rate cache"
    assert c_names == py_names


def test_electron_thermo_cache_filename_identical(coarse_dirs):
    """Both backends name the e+- thermodynamic cache identically.

    Pins the electron-thermo fingerprint port, including constants_hash: the
    hash is in the filename, so this is a direct comparison of what each
    backend computed the fingerprint to be.
    """
    dirs = coarse_dirs
    c_names = _names(dirs["c"], "plasma", "electron_thermo_")
    py_names = _names(dirs["python"], "plasma", "electron_thermo_")
    assert c_names, "C backend wrote no electron-thermo cache"
    assert c_names == py_names


def test_thermal_cache_fingerprints_agree(thermal_dirs):
    """Both backends resolve the CCRTh thermal cache to the same shipped file.

    The thermal table's *contents* are deliberately not compared -- both
    backends compute it by Monte-Carlo with independent RNG streams, and
    forcing a recompute costs minutes per backend (see the module docstring).
    Its *fingerprint* is still worth pinning, and can be, for free.

    The trick is to assert a cache **hit on both sides**: each backend runs with
    its own empty ``cache_dir``, so the only thermal table either can find is
    the shipped one, reached through the overlay. If the two implementations of
    the thermal fingerprint ever disagreed, the one that drifted would name a
    file that does not exist, miss, recompute, and write a *differently named*
    table into its own ``cache_dir`` -- which is exactly what this asserts is
    empty.

    Note the failure mode: a genuine divergence makes this test **slow** (one
    backend pays a full vegas recompute) before it fails. That is the price of
    pinning this fingerprint without paying the recompute on every green run,
    and it is the right trade -- the slow path only happens when something is
    already broken.
    """
    dirs = thermal_dirs
    c_names = _names(dirs["c"], "weak", "nTOp_thermal_")
    py_names = _names(dirs["python"], "weak", "nTOp_thermal_")
    assert c_names == [], (
        "C backend missed the shipped thermal cache and recomputed "
        f"{c_names} -- cpr_thermal_fingerprint has drifted from "
        "thermal_fingerprint")
    assert py_names == [], (
        "Python backend missed the shipped thermal cache and recomputed "
        f"{py_names} -- the shipped table's fingerprint no longer matches "
        "thermal_fingerprint (was it re-keyed?)")


def test_constants_hash_identical_across_backends():
    """Python's constants_hash(cache) equals the C backend's per-cache hash.

    Each cache hashes only the constants it reads (cache_utils.CACHE_CONSTANTS,
    mirrored by cpr_constants_hash), in the same canonical-JSON form. If this
    ever fails, the cause is either the two subsets drifting apart -- which
    makes every shared cache file a cross-backend miss -- or a build flag
    (-ffast-math, FMA contraction) changing a constant's last bit between
    CPython and the C compiler. Both are worth surfacing loudly.

    Checked indirectly but exactly: the value appears inside the fingerprint
    JSON header that each backend writes, and the filenames compared above are
    derived from it. Here the Python side's own contract is pinned -- 16 hex
    digits, stable across calls, and distinct per cache.
    """
    from primat.cache_utils import CACHE_CONSTANTS, constants_hash
    hashes = {c: constants_hash(c) for c in CACHE_CONSTANTS}
    for c, h in hashes.items():
        assert len(h) == 16, c
        assert h == constants_hash(c), c
    assert len(set(hashes.values())) == len(hashes), hashes


# ---------------------------------------------------------------------------
# 2. Column agreement.
# ---------------------------------------------------------------------------

def test_weak_rate_cache_columns_agree(coarse_dirs):
    """The two backends' n<->p rate tables agree column by column."""
    dirs = coarse_dirs
    names = _names(dirs["c"], "weak", "nTOp_")
    name = [n for n in names if not n.startswith("nTOp_thermal_")][0]
    c = np.loadtxt(dirs["c"] / "weak" / name)
    p = np.loadtxt(dirs["python"] / "weak" / name)
    assert c.shape == p.shape
    assert _max_rel(c, p) < NTOP_RTOL


def test_electron_thermo_cache_columns_agree(coarse_dirs):
    """The two backends' e+- thermodynamic tables agree in all four columns.

    This is the comparison that catches a gap like the 1.0e-4 one a
    quadrature tolerance floor once opened in rho_e/p_e.
    """
    dirs = coarse_dirs
    name = _names(dirs["c"], "plasma", "electron_thermo_")[0]
    c = np.loadtxt(dirs["c"] / "plasma" / name)
    p = np.loadtxt(dirs["python"] / "plasma" / name)
    assert c.shape == p.shape
    assert _max_rel_above_scale(c, p, ELECTRON_THERMO_TAIL_FRAC) < ELECTRON_THERMO_RTOL
    assert _max_scale_rel(c, p) < ELECTRON_THERMO_SCALE_RTOL


def test_qed_pressure_tables_agree(qed_dirs):
    """The two backends' QED pressure-correction tables agree.

    Forced recompute on both sides (the shipped tables would otherwise be a hit
    for both, comparing a file with itself). Compared on the column-scale
    metric -- see QED_SCALE_RTOL for why pointwise-relative is not the right
    question for these two tables.
    """
    dirs = qed_dirs
    for name in ("QED_pressure_correction_e2.txt",
                 "QED_pressure_correction_e3.txt"):
        c = np.loadtxt(dirs["c"] / "plasma" / name)
        p = np.loadtxt(dirs["python"] / "plasma" / name)
        assert c.shape == p.shape, name
        # The T column must match exactly: it is a pure logspace grid, not an
        # integration result, so any difference is a grid bug rather than noise.
        assert np.array_equal(c[:, 0], p[:, 0]), f"{name}: T grids differ"
        assert _max_scale_rel(c, p) < QED_SCALE_RTOL, name


def test_qed_fingerprint_header_identical(qed_dirs):
    """Both backends write byte-identical QED fingerprint headers.

    The QED tables keep fixed filenames, so unlike the two hash-named families
    their fingerprint is only visible in the header -- which makes comparing
    the header the only way to pin cpr_qed_fingerprint against
    qed_pressure.qed_fingerprint.
    """
    dirs = qed_dirs

    def fp_lines(path):
        out = []
        with open(path) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                if line.startswith("# fingerprint"):
                    out.append(line.rstrip("\n"))
        return out

    for name in ("QED_pressure_correction_e2.txt",
                 "QED_pressure_correction_e3.txt"):
        c = fp_lines(dirs["c"] / "plasma" / name)
        p = fp_lines(dirs["python"] / "plasma" / name)
        assert c, f"{name}: C backend wrote no fingerprint header"
        assert c == p, name


@pytest.fixture(scope="module")
def perturbed_constant_dirs(tmp_path_factory):
    """Both backends on the coarse grid with me raised 1 %.

    me is the constant with the widest reach -- the e+- thermodynamics, the
    QED pressure tables and the n<->p phase space all read it -- so it
    exercises every fingerprint at once.
    """
    return _run_both(tmp_path_factory.mktemp("perturbed"),
                     extra={"me": CONST.me * 1.01})


def test_overridden_constant_rekeys_both_backends_identically(
        coarse_dirs, perturbed_constant_dirs):
    """A constant override changes the cache filenames, the same way on both.

    Two failure modes at once. If the names did not *change*, a run with an
    overridden constant would load the default-constant table and report its
    physics with no warning. If they changed *differently*, the backends would
    stop sharing a cache they are meant to share -- and the ordinary parity
    tests could not see it, since both would still agree with themselves.
    """
    for sub, prefix in (("weak", "nTOp_"), ("plasma", "electron_thermo_")):
        base_c = _names(coarse_dirs["c"], sub, prefix)
        bumped_c = _names(perturbed_constant_dirs["c"], sub, prefix)
        bumped_py = _names(perturbed_constant_dirs["python"], sub, prefix)
        assert bumped_c, f"C backend wrote no {sub} cache with me overridden"
        assert bumped_c != base_c, (
            f"{sub}/{prefix}* is unchanged by a 1 % shift in me: the override "
            "would silently reuse the default-constant table")
        assert bumped_c == bumped_py
