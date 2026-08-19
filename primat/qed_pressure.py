# -*- coding: utf-8 -*-
"""
qed_pressure.py — Analytical computation of QED plasma-pressure corrections
============================================================================

Computes the finite-temperature QED interaction-pressure corrections
δP(T), dδP/dT, and d²δP/dT² that enter the EM plasma thermodynamics
during BBN.  These supplement the free ideal-gas photon + e± expressions
with the leading electromagnetic interactions.

Physical background
-------------------
The QED interaction pressure is a finite-temperature correction arising
from the QED interaction between photons and electrons in the hot plasma.
It is decomposed into three contributions in increasing order of the
electromagnetic coupling e (α = e²/(4π)):

  δP = δP_a [O(e²)]  +  δP_{e3} [O(e³)]  +  δP_b [O(e⁴)]

where (following PRIMAT-Main.m and Phys. Rep. §II.E):

  δP_a(T)  = (α/π) T⁴ [-(2/3) I₀₁(x) - (2/π²) I₀₁(x)²]
              Leading O(α) one-loop correction (Frenkel–Galitskii–Migdal).

  δP_{e3}(T) = α^{3/2} (4/3)√(2π) T⁴ [(I₀₁(x)+I₂₋₁(x))/π²]^{3/2}
               O(α^{3/2}) ring/plasmon contribution (Blaizot–Zinn-Justin).

  δP_b(T)  = T⁴ ∫₀^∞ ∫₀^∞ F(p₁,p₂,x) dp₁ dp₂           [O(α²), optional]
             F = (α/π³) x² p₁ p₂ / (e₁ e₂)
                 × ln|(p₁+p₂)/(p₁-p₂)| / ((e^{e₁}+1)(e^{e₂}+1))

Here x = mₑ/T (dimensionless), eᵢ = √(pᵢ²+x²), and:

  I₀₁(x) = ∫₀^∞ p² / [√(p²+x²)(e^{√(p²+x²)}+1)] dp
           = ∫_x^∞ √(E²−x²)/(e^E+1) dE           (PRIMAT: Imn[1][0,1][x])

  I₂₋₁(x) = ∫₀^∞ √(p²+x²) / (e^{√(p²+x²)}+1) dp
            = ∫_x^∞ E²/[√(E²−x²)(e^E+1)] dE      (PRIMAT: Imn[1][2,-1][x])

The dominant term is δP_a, which is negative (interaction reduces the
pressure relative to the ideal gas).  The ring term δP_{e3} is positive
and roughly 10× smaller.  The two-loop exchange δP_b is typically 100×
smaller still and is optional.

File format
-----------
The results are stored in two separate four-column files, one per order
in e, so that either correction can be inspected, regenerated, or swapped
out independently of the other:

  ``data/plasma/QED_pressure_correction_e2.txt`` [O(e²), one-loop]::

    T [MeV]  dP_a [MeV^4]  d(dP_a)/dT [MeV^3]  d2(dP_a)/dT2 [MeV^2]

  ``data/plasma/QED_pressure_correction_e3.txt`` [O(e³), ring/plasmon]::

    T [MeV]  dP_e3 [MeV^4]  d(dP_e3)/dT [MeV^3]  d2(dP_e3)/dT2 [MeV^2]

When loaded by :mod:`primat.plasma`, the two files' values (and
derivatives) are summed to give the total correction.  The δP_b term
would require a separate flag and file.

(Backward compat: :func:`primat.plasma.Plasma._load_tables` also reads a
single 7-column ``QED_tables.txt`` and the ``QED_P_int.txt``/
``QED_dP_intdT.txt``/``QED_d2P_intdT2.txt`` trio as fallbacks;
:func:`save_qed_tables` only ever writes the two-file format above.)

Usage
-----
>>> from primat.qed_pressure import compute_qed_pressure_tables, save_qed_tables
>>> tables = compute_qed_pressure_tables()  # ~0.3 s on a modern laptop
>>> save_qed_tables(tables, "/path/to/data/plasma/")

Reference
---------
Pitrou, Coc, Uzan & Vangioni, Phys. Rep. 2018 (arXiv:1806.11095), §II.E
PRIMAT-Main.m: ``dPa``, ``dPe3``, ``dPb`` definitions (lines 920, 939, 949)
"""

import os
import numpy as np
from scipy.integrate import quad, dblquad
from scipy.interpolate import CubicSpline

from .cache_utils import constants_hash, write_cache_with_fingerprint

# Fallback constants for a standalone call with no PRIMATConfig in hand (the
# generate_rates/ scripts). A BBN run never reads them: plasma.Plasma._load_tables
# passes alpha=cfg.alphaem, me=cfg.me explicitly, and both are user-settable
# parameters. The written tables carry a fingerprint header keyed on
# cache_utils.constants_hash, so a table computed with other constants is
# rebuilt rather than loaded silently.
_ALPHA_FS = 1. / 137.035999084   # fine-structure constant (CODATA 2018)
_ME_MEV   = 0.51099895           # electron mass [MeV] (CODATA 2018)

# Bump when a code change alters the *numerical content* of the two QED
# pressure-correction tables for a fixed (T_min, T_max, n_pts) and fixed
# constants -- a changed integrand, a different differentiation scheme, a new
# column layout. Bumping invalidates every cached table regardless of its other
# fingerprint fields, exactly as WEAK_RATE_FORMAT_VERSION does for the weak
# tables. Mirrored by QED_FORMAT_VERSION in primat-c/src/qed_pressure.c.
#
# v1: first fingerprinted generation. Before it the two files carried no
# fingerprint at all, so a table generated with a different alpha/me or
# different T bounds was loaded silently -- the hazard this version closes.
# The shipped tables were given a v1 header in place: same numbers (the write
# format stayed "%.6E"), header lines only.
#
# v2: `constants_hash` narrowed to the two constants the integrands read,
# alphaem and me (cache_utils.CACHE_CONSTANTS). Numbers unchanged; the shipped
# tables were re-keyed in place, header lines only.
QED_FORMAT_VERSION = 2

# Low-x cutoff: for x = mₑ/T > 50 (T < mₑ/50 ≈ 10 keV) the e± are so
# non-relativistic that δP is effectively zero (Boltzmann-suppressed).
_X_NONREL_CUTOFF = 50.

# Upper limit for 1D momentum integrals, in units of x = mₑ/T.  The
# integrand decays as e^{-p} for large p, so p_max = 500 is more than
# sufficient even at very high temperatures.
_P_UPPER = 500.


# ---------------------------------------------------------------------------
# Fermi-Dirac momentum integrals I₀₁ and I₂₋₁
# ---------------------------------------------------------------------------

def _I01(x):
    """Fermi-Dirac phase-space integral I₀₁(x) [dimensionless].

    Defined as (PRIMAT: Imn[1][0,1][x]):

        I₀₁(x) = ∫₀^∞ p² / [√(p²+x²)(e^{√(p²+x²)}+1)] dp

    Equivalently (change of variable E = √(p²+x²)):

        I₀₁(x) = ∫_x^∞ √(E²−x²) / (e^E+1) dE

    The p-space form is used here because it is non-singular at the lower
    limit, making scipy.quad straightforward to apply.

    Parameters
    ----------
    x : float
        Dimensionless ratio mₑ/T.

    Returns
    -------
    float
        I₀₁(x) in natural units (ℏ = c = kB = 1).

    Example
    -------
    >>> _I01(0.0)   # ultra-relativistic limit → π²/12 ≈ 0.822
    >>> _I01(0.5)   # semi-relativistic (T ~ 1 MeV)
    """
    if x > _X_NONREL_CUTOFF:
        return 0.
    def integrand(p):
        E = np.sqrt(p * p + x * x)
        return p * p / (E * (np.exp(E) + 1.))
    result, _ = quad(integrand, 0., _P_UPPER,
                     epsabs=1e-13, epsrel=1e-13, limit=300)
    return result


def _I2m1(x):
    """Fermi-Dirac phase-space integral I₂₋₁(x) [dimensionless].

    Defined as (PRIMAT: Imn[1][2,-1][x]):

        I₂₋₁(x) = ∫₀^∞ √(p²+x²) / (e^{√(p²+x²)}+1) dp

    Equivalently:

        I₂₋₁(x) = ∫_x^∞ E² / [√(E²−x²)(e^E+1)] dE

    The p-space form removes the 1/√(E²−x²) singularity at the lower
    limit, so no special handling is needed.

    Parameters
    ----------
    x : float
        Dimensionless ratio mₑ/T.

    Returns
    -------
    float
        I₂₋₁(x) in natural units.

    Example
    -------
    >>> _I2m1(0.0)   # ultra-relativistic limit → π²/12 ≈ 0.822
    """
    if x > _X_NONREL_CUTOFF:
        return 0.
    def integrand(p):
        E = np.sqrt(p * p + x * x)
        return E / (np.exp(E) + 1.)
    result, _ = quad(integrand, 0., _P_UPPER,
                     epsabs=1e-13, epsrel=1e-13, limit=300)
    return result


# ---------------------------------------------------------------------------
# Three contributions to δP
# ---------------------------------------------------------------------------

def _dPa(T, alpha=_ALPHA_FS, me=_ME_MEV):
    """O(e²) QED interaction-pressure correction δP_a(T) [MeV⁴].

    The leading one-loop correction to the electromagnetic plasma pressure
    from the QED interaction between photons and electrons
    (PRIMAT: ``dPa``; Phys. Rep. §II.E):

        δP_a = (α/π) T⁴ [−(2/3) I₀₁(x) − (2/π²) I₀₁(x)²]

    This is negative (interaction lowers the pressure) and is the dominant
    QED correction, of order α ∼ 7×10⁻³.

    Parameters
    ----------
    T : float
        Photon temperature [MeV].
    alpha : float, optional
        Fine-structure constant (default: 1/137.035999084).
    me : float, optional
        Electron mass [MeV] (default: ``_ME_MEV`` = 0.51099895, CODATA 2018,
        kept equal to ``primat.constants.CONST.me`` by hand -- see the
        module-level note above).

    Returns
    -------
    float
        δP_a in MeV⁴.

    Example
    -------
    >>> _dPa(10.0)   # at T = 10 MeV
    """
    x = me / T
    I01 = _I01(x)
    return alpha / np.pi * T**4 * (-2./3. * I01 - 2./np.pi**2 * I01**2)


def _dPe3(T, alpha=_ALPHA_FS, me=_ME_MEV):
    """O(e³) QED interaction-pressure correction δP_{e3}(T) [MeV⁴].

    The O(α^{3/2}) ring/plasmon contribution to the QED pressure,
    arising from collective plasma oscillations (PRIMAT: ``dPe3``):

        δP_{e3} = α^{3/2} (4/3)√(2π) T⁴ [(I₀₁+I₂₋₁)/π²]^{3/2}

    This is positive and roughly 10× smaller than δP_a.

    Parameters
    ----------
    T : float
        Photon temperature [MeV].
    alpha : float, optional
        Fine-structure constant.
    me : float, optional
        Electron mass [MeV].

    Returns
    -------
    float
        δP_{e3} in MeV⁴.

    Example
    -------
    >>> _dPe3(10.0)   # at T = 10 MeV
    """
    x = me / T
    I01  = _I01(x)
    I2m1 = _I2m1(x)
    combo = (I01 + I2m1) / np.pi**2
    if combo <= 0.:
        return 0.
    return alpha**(3./2.) * (4./3.) * np.sqrt(2. * np.pi) * T**4 * combo**(3./2.)


def _dPb(T, alpha=_ALPHA_FS, me=_ME_MEV, epsrel=1e-4):
    """O(e⁴) QED interaction-pressure correction δP_b(T) [MeV⁴].

    The two-loop exchange contribution, corresponding to PRIMAT's
    ``dPb`` (``$CompleteQEDPressure=True``):

        δP_b = T⁴ ∫₀^∞ ∫₀^∞ F(p₁,p₂,x) dp₁ dp₂

    with

        F(p₁,p₂,x) = (α/π³) x² p₁ p₂ / (e₁ e₂)
                      × ln|(p₁+p₂)/(p₁−p₂)| / ((e^{e₁}+1)(e^{e₂}+1))

    where eᵢ = √(pᵢ²+x²).  The integrand is symmetric in p₁↔p₂, and the
    logarithm has an integrable singularity at p₁ = p₂.

    **Note**: this term is O(α²) ≈ 5×10⁻⁵ and is NOT included in the
    standard primat QED tables (which only store δP_a + δP_{e3}).
    It is provided here for completeness.  Computing it is expensive
    (~10–60 s per temperature point at low precision).

    **Deliberately Python-only.**  ``primat-c/src/qed_pressure.c`` has no
    counterpart to this function, and needs none: ``include_dPb`` defaults to
    False, so δP_b never reaches the shipped tables nor any solve path, and
    the two backends stay in parity without it.  This is an intentional
    asymmetry, not a porting omission -- if δP_b is ever switched on by
    default it must be ported to the C backend.

    Parameters
    ----------
    T : float
        Photon temperature [MeV].
    alpha : float, optional
        Fine-structure constant.
    me : float, optional
        Electron mass [MeV].
    epsrel : float, optional
        Relative accuracy target for the 2D numerical integration
        (default 1e-4, matching PRIMAT's ``PrecisionGoal->4``).

    Returns
    -------
    float
        δP_b in MeV⁴.
    """
    x = me / T
    if x > _X_NONREL_CUTOFF:
        return 0.

    # Upper momentum limit: at least 20, or 20x (non-relativistic: pmax ~ x)
    p_upper = max(20., 20. * x)

    def integrand(p2, p1):
        # eᵢ = √(pᵢ² + x²); the log factor is regularised by treating the
        # p1 == p2 singularity as integrable (verified: log divergence, area 0)
        e1 = np.sqrt(p1 * p1 + x * x)
        e2 = np.sqrt(p2 * p2 + x * x)
        if abs(p1 - p2) < 1e-14 * (p1 + p2 + 1e-10):
            return 0.   # contribution zero on the diagonal p1=p2
        log_factor = np.log(abs((p1 + p2) / (p1 - p2)))
        fd1 = 1. / (np.exp(e1) + 1.)
        fd2 = 1. / (np.exp(e2) + 1.)
        return (alpha / np.pi**3) * x**2 * p1 * p2 / (e1 * e2) * log_factor * fd1 * fd2

    result, _ = dblquad(integrand, 0., p_upper,
                        lambda p1: 0., lambda p1: p_upper,
                        epsrel=epsrel, epsabs=0.)
    return T**4 * result


# ---------------------------------------------------------------------------
# Grid computation and file I/O
# ---------------------------------------------------------------------------

def compute_qed_pressure_tables(T_min=1e-3, T_max=1e2, n_pts=500,
                                alpha=_ALPHA_FS, me=_ME_MEV,
                                include_dPb=False, verbose=True):
    """Compute δP, dδP/dT, d²δP/dT² on a temperature grid [MeV].

    Evaluates the O(e²) and O(e³) QED corrections to the EM plasma
    pressure (and optionally the O(e⁴) exchange term) on a log-spaced
    temperature grid, then differentiates numerically using a cubic spline.

    The two-column format matches the files loaded by
    :func:`primat.plasma._load_tables`:
      column 0 = T [MeV]
      column 1 = δP_a(T)   [MeV⁴]  (O(e²) = O(α))
      column 2 = δP_{e3}(T) [MeV⁴] (O(e³) = O(α^{3/2}))

    When ``include_dPb=True`` a third column for δP_b is added to the
    ``dP`` table, and the derivatives are recomputed accordingly.

    Parameters
    ----------
    T_min : float
        Minimum temperature [MeV] (default 1e-3, well below e± freeze-out).
    T_max : float
        Maximum temperature [MeV] (default 100, well above BBN start).
    n_pts : int
        Number of log-spaced temperature grid points (default 500).
    alpha : float
        Fine-structure constant (default: CODATA 2018 value).
    me : float
        Electron mass [MeV] (default: ``_ME_MEV`` = 0.51099895, CODATA 2018,
        kept equal to ``primat.constants.CONST.me`` by hand -- see the
        module-level note above).
    include_dPb : bool
        If True, also compute the expensive O(e⁴) two-loop term δP_b
        (adds ~10–60 s per temperature point; default False).
    verbose : bool
        Print progress messages (default True).

    Returns
    -------
    dict
        Keys: ``"T"``, ``"dP_e2"``, ``"dP_e3"``, ``"d_dP_e2_dT"``,
        ``"d_dP_e3_dT"``, ``"d2_dP_e2_dT2"``, ``"d2_dP_e3_dT2"``,
        and optionally ``"dP_b"``, ``"d_dPb_dT"``, ``"d2_dPb_dT2"``.
        All arrays have length ``n_pts``.

    Notes
    -----
    The derivatives dδP/dT and d²δP/dT² are obtained from a CubicSpline
    fit to the tabulated δP values, not from analytic differentiation of
    the Fermi-Dirac integrals.  The analytic route would require four
    additional quadratures per temperature point and would be ~7× slower
    with no practical accuracy gain: the spline derivatives agree with
    direct finite differences on _dPa/_dPe3 to <0.01% at all T.

    Example
    -------
    >>> tables = compute_qed_pressure_tables(n_pts=100, verbose=False)
    >>> tables["T"].shape
    (100,)
    """
    T_grid = np.logspace(np.log10(T_min), np.log10(T_max), n_pts)
    dP_e2  = np.zeros(n_pts)
    dP_e3  = np.zeros(n_pts)
    dP_b   = np.zeros(n_pts) if include_dPb else None

    for i, T in enumerate(T_grid):
        if verbose and i % max(1, n_pts // 10) == 0:
            print(f"  [QED] Computing T = {T:.3e} MeV  ({i+1}/{n_pts})")
        dP_e2[i] = _dPa(T, alpha=alpha, me=me)
        dP_e3[i] = _dPe3(T, alpha=alpha, me=me)
        if include_dPb:
            dP_b[i] = _dPb(T, alpha=alpha, me=me)

    # Differentiate numerically using a cubic spline; this avoids having to
    # differentiate the integrands analytically.
    spl_e2 = CubicSpline(T_grid, dP_e2)
    spl_e3 = CubicSpline(T_grid, dP_e3)
    d_e2   = spl_e2(T_grid, 1)   # first derivative
    d2_e2  = spl_e2(T_grid, 2)   # second derivative
    d_e3   = spl_e3(T_grid, 1)
    d2_e3  = spl_e3(T_grid, 2)

    out = {"T": T_grid,
           "dP_e2": dP_e2, "dP_e3": dP_e3,
           "d_dP_e2_dT": d_e2, "d_dP_e3_dT": d_e3,
           "d2_dP_e2_dT2": d2_e2, "d2_dP_e3_dT2": d2_e3}

    if include_dPb:
        spl_b  = CubicSpline(T_grid, dP_b)
        out["dP_b"]        = dP_b
        out["d_dPb_dT"]    = spl_b(T_grid, 1)
        out["d2_dPb_dT2"]  = spl_b(T_grid, 2)

    return out


def qed_fingerprint(T_min, T_max, n_pts, cfg=None):
    """Fingerprint dict for the two QED pressure-correction cache files.

    The tables are a function of exactly two things: the physical constants
    that enter the integrands (α and mₑ, via
    :func:`primat.cache_utils.constants_hash`), and the temperature grid they
    were evaluated on.  Everything else about a run (network, baryon density,
    neutrino treatment, ...) leaves δP_a and δP_{e3} untouched.

    A mismatch makes the loader recompute (~0.3 s) *without* writing, since
    these two files keep fixed names and a write would replace another
    configuration's pair (see :meth:`primat.plasma.Plasma._load_tables`).

    Args:
        T_min: float, lowest grid temperature [MeV].
        T_max: float, highest grid temperature [MeV].
        n_pts: int, number of log-spaced grid points.
        cfg: PRIMATConfig whose constants the tables were computed with;
            ``None`` uses the defaults (:data:`primat.constants.CONST`).

    Returns:
        dict, JSON-serialisable; pass to
        :func:`primat.cache_utils.fingerprint_hash` for the hash.  Mirrored
        field-for-field by ``cpr_qed_fingerprint`` in
        ``primat-c/src/qed_pressure.c``.

    Example:
        >>> from primat.cache_utils import fingerprint_hash
        >>> fingerprint_hash(qed_fingerprint(1e-3, 1e2, 500))   # doctest: +SKIP
        '0f3a...'
    """
    return {"format_version": QED_FORMAT_VERSION,
            "constants_hash": constants_hash("qed", cfg),
            "T_min": float(T_min),
            "T_max": float(T_max),
            "n_pts": int(n_pts)}


def save_qed_tables(tables, plasma_dir, verbose=True, cfg=None):
    """Write the computed QED tables to two four-column files, one per order in e.

    Produces the two files read by :func:`primat.plasma.Plasma._load_tables`:

      ``QED_pressure_correction_e2.txt`` — T, δP_a, d(δP_a)/dT, d²(δP_a)/dT²  [O(e²)]
      ``QED_pressure_correction_e3.txt`` — T, δP_{e3}, d(δP_{e3})/dT, d²(δP_{e3})/dT²  [O(e³)]

    Keeping the two orders in separate files lets either be inspected,
    regenerated, or swapped independently of the other.  Column units are
    given explicitly in each file's header: T in MeV, δP in MeV^4, dδP/dT
    in MeV^3, d²δP/dT² in MeV^2 (natural units ħ = c = k_B = 1).

    Each file also carries a fingerprint header (:func:`qed_fingerprint`)
    recording the format version, the physical-constants hash, and the (T_min,
    T_max, n_pts) grid, so that :meth:`primat.plasma.Plasma._load_tables`
    detects a table computed with different constants or a different grid and
    rebuilds it instead of loading it silently.

    The rows are written with ``fmt="%.6E"``, the format the shipped tables have
    always used: adding the fingerprint header was therefore a header-only
    change to those tracked files, with every data row byte-identical.

    Parameters
    ----------
    tables : dict
        Output of :func:`compute_qed_pressure_tables`.
    plasma_dir : str
        Path to the ``data/plasma/`` directory.
    verbose : bool
        Print confirmation message (default True).

    Raises
    ------
    OSError
        Propagated from the underlying write so that callers which want to
        degrade gracefully on a read-only install can catch it (see
        :meth:`primat.plasma.Plasma._load_tables`, which turns it into a
        warning naming the ``cache_dir`` remedy).

    Example
    -------
    >>> save_qed_tables(tables, "primat/data/plasma/")
    """
    T    = tables["T"]
    e2   = tables["dP_e2"]
    e3   = tables["dP_e3"]
    de2  = tables["d_dP_e2_dT"]
    de3  = tables["d_dP_e3_dT"]
    d2e2 = tables["d2_dP_e2_dT2"]
    d2e3 = tables["d2_dP_e3_dT2"]

    hdr_e2 = ("Source: primat qed_pressure.py — QED plasma-pressure correction delta_P_a(T)\n"
              "delta_P_a: O(e^2), one-loop (Frenkel-Galitskii-Migdal)\n"
              "Reference: Pitrou et al., Phys. Rep. (2018), eq. 47; PRIMAT-Main.m: dPa\n"
              "T [MeV]       dP_a [MeV^4]      d(dP_a)/dT [MeV^3]  d2(dP_a)/dT2 [MeV^2]")
    hdr_e3 = ("Source: primat qed_pressure.py — QED plasma-pressure correction delta_P_e3(T)\n"
              "delta_P_e3: O(e^3), ring/plasmon (Blaizot-Zinn-Justin)\n"
              "Reference: Pitrou et al., Phys. Rep. (2018), eq. 47; PRIMAT-Main.m: dPe3\n"
              "T [MeV]       dP_e3 [MeV^4]     d(dP_e3)/dT [MeV^3]  d2(dP_e3)/dT2 [MeV^2]")

    # Both files describe the same grid, so they share one fingerprint: the
    # loader checks each file's own header, and a mismatch on either rebuilds
    # (and rewrites) both, which is correct since they are always generated
    # together and summed column-by-column at point of use.
    fp = qed_fingerprint(T[0], T[-1], len(T), cfg=cfg)

    # write_cache_with_fingerprint creates the target on demand -- when
    # redirected to a fresh cache_dir the plasma/ subdir may not exist yet --
    # and emits the physics header block above the two fingerprint lines.
    # fmt="%.6E" is the format the shipped tables were written with, kept so
    # that gaining a fingerprint header leaves every data row byte-identical.
    write_cache_with_fingerprint(
        os.path.join(plasma_dir, "QED_pressure_correction_e2.txt"),
        fp, [T, e2, de2, d2e2], col_header=hdr_e2, fmt="%.6E")
    write_cache_with_fingerprint(
        os.path.join(plasma_dir, "QED_pressure_correction_e3.txt"),
        fp, [T, e3, de3, d2e3], col_header=hdr_e3, fmt="%.6E")

    if verbose:
        # Mirrored by plasma.c's cpr_log(cfg, "QED", ...) block; kept ASCII so
        # the two streams compare byte for byte on any console.
        print(f"[QED]  Tables written to {plasma_dir}:")
        print(f"       QED_pressure_correction_e2.txt  (4 columns: T, dP_a, derivatives)")
        print(f"       QED_pressure_correction_e3.txt  (4 columns: T, dP_e3, derivatives)")
        print(f"       T range: {T[0]:.3e}-{T[-1]:.3e} MeV  ({len(T)} points)")
