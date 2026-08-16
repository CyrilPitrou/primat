# -*- coding: utf-8 -*-
"""
weak_rates.api — public n<->p weak-rate entry points
===============================================================================

Top-level functions used by the rest of primat (mainly
primat.background): :func:`ComputeWeakRates` (raw computation),
:func:`InterpolateWeakRates` (load the fingerprinted on-disk cache without
recomputing), and :func:`RecomputeWeakRates` (the usual entry point: load the
cache if the fingerprint matches, else recompute and write it).

Normalisation: the stored rates are in units of 1/tau_n rather than absolute
s^-1 (Phys. Rep. Eqs. 89-91) — the caller multiplies by 1/cfg.tau_n (or by
the analytically computed K/Fn) after loading; see
:func:`weak_rates.corrections.ComputeFn` for the free-neutron-decay
phase-space normalisation Fn.
"""

import os

import numpy as np
from scipy.interpolate import interp1d

from .cache import n_points_per_decade, _weak_rate_fingerprint
from .corrections import _build_rate_context, _correction_terms, \
    _thermal_correction_interpolants, ComputeFn
from ..cache_utils import (fingerprint_hash, write_cache_with_fingerprint,
                           resolve_cache_file, cache_write_dir)

__all__ = ['ComputeWeakRates', 'InterpolateWeakRates', 'RecomputeWeakRates']


def _weak_rate_loglog_interp(T, rate):
    """Interpolant for a tabulated n<->p weak rate, in log10-log10 space.

    Physics/numerics rationale
    --------------------------
    The forward (n->p) and backward (p->n) rates span ~5 orders of magnitude
    across the BBN temperature grid and are smooth, nearly power-law functions
    of T, so a cubic spline in ``(log10 T, log10 Gamma)`` space is far more
    accurate than one in linear (T, Gamma) space at the same node density.
    This is the scheme the shipped nuclear rate tables use
    (:func:`network_data._resample_rate_table`) and the one the C backend uses
    for the weak rates (``cpr_weak_rate_nTOp``). ``interp1d(kind='cubic')`` is
    scipy's not-a-knot cubic, matching C's ``cpr_cubic_spline_fit_notaknot``
    boundary condition — the two backends must not differ here, since a
    weak-rate scheme mismatch dominates every other parity term (see
    ``tests/test_backend_parity.py``).

    Non-positive prefix (backward rate)
    -----------------------------------
    The backward p->n rate is ``exp(-Q/T)``-suppressed and is clamped to exact
    0 (below the ~1e-28 noise floor of ``ComputeWeakRates``) in a low-temperature
    prefix of the grid, where ``log10`` is undefined. The spline is therefore
    built only over the contiguous all-positive *suffix* of the table, and the
    rate is returned as 0 for ``T`` below that suffix -- physically the rate is
    negligible there and the caller
    (:meth:`background.Background.weak_nTOp_bkwrd`) clamps it to 0 anyway. The
    forward rate is strictly positive over the whole grid, so its suffix is the
    full table.

    Parameters
    ----------
    T    : array of temperatures [K], strictly increasing.
    rate : array of rates [1/tau_n] at those temperatures (>= 0).

    Returns
    -------
    callable ``f(T_K) -> rate`` accepting a scalar or array of temperatures [K]
    and returning the rate in units of 1/tau_n, log10-log10 cubic-interpolated
    (log10-log10 linear if the positive suffix has fewer than the 4 knots a
    not-a-knot cubic needs -- unreachable with the default grids, and matching
    the C backend's identical fallback).
    """
    T = np.asarray(T, dtype=float)
    rate = np.asarray(rate, dtype=float)
    # Largest index carrying a non-positive value: everything strictly above it
    # is the contiguous positive suffix we can take log10 of (the zeros sit in
    # a low-T prefix for the backward rate; the forward rate has none).
    nonpos = np.nonzero(rate <= 0.0)[0]
    i0 = int(nonpos[-1]) + 1 if nonpos.size else 0
    logT = np.log10(T[i0:])
    logR = np.log10(rate[i0:])
    # Pin the rate to 0 only below a genuine non-positive (clamped-to-zero)
    # prefix (the backward p->n rate). When the whole grid is positive (the
    # forward n->p rate, i0 == 0) there is no such region and we extrapolate
    # the spline below the first node, matching the old fill_value="extrapolate".
    T_zero_below = T[i0] if i0 > 0 else -np.inf
    kind = 'cubic' if logT.size >= 4 else 'linear'
    spline = interp1d(logT, logR, kind=kind, bounds_error=False,
                      fill_value='extrapolate')

    def _eval(T_K):
        Tq = np.asarray(T_K, dtype=float)
        # log10 is undefined at/below 0; clamp the argument (values there are
        # masked back to 0 below anyway) to avoid a spurious warning.
        val = np.power(10.0, spline(np.log10(np.maximum(Tq, 1e-300))))
        val = np.where(Tq >= T_zero_below, val, 0.0)
        return val if val.ndim else float(val)

    # Expose the fitted spline and the zero-mask threshold so the fast scalar
    # evaluator (weak_rates.fast_eval.FastWeakRate) can extract the PPoly and
    # JIT the hot scalar path. Does not change _eval's behaviour.
    _eval.spline = spline
    _eval.T_zero_below = T_zero_below
    return _eval

def ComputeWeakRates(Tvec, cfg, dFDneu_func=None, dFDneu_moments=None):
    """Compute the non-thermal n↔p weak rate tables over the BBN T range.

    Evaluates the forward rate Γ_{n→p}(T) and backward rate Γ_{p→n}(T) on the
    photon-temperature grid Tg_vec, by summing the additive correction terms
    (controlled by the cfg flags) returned by :func:`_correction_terms`:

    Γ_{n→p} = K × [ (CCR or Born) (+ FMCCR/FMNoCCR) (+ SD_CCR/SD) (+ SD_FMCCR/SD_FM) ]

    which flag selects which term is :func:`_correction_terms`'s business; the
    terms themselves are::

      _L_BORN    — Born rate ∫ p² [χ₊(E)+χ₊(−E)] dp  (Phys. Rep. Eqs. 77–78).
      _L_CCR     — Born integrand × FermiCoulomb × RadCorrResum (T=0 Coulomb
                   + resummed radiative corrections; Phys. Rep. Eq. 101).
      _L_FMCCR   — Finite-nucleon-mass correction × Coulomb × radiative
                   (Fokker-Planck expansion; Phys. Rep. §III.G).
      _L_FMNoCCR — Finite-nucleon-mass correction alone.
      _L_SD      — Born-level spectral-distortion correction (deviation of f_ν
                   from Fermi–Dirac, passed in via dFDneu_func).
      _L_SD_CCR  — Same with Coulomb × radiative factor.
      _L_SD_FMNoCCR — Finite-nucleon-mass correction applied to the distortion
                   instead of the plain Fermi-Dirac
                   (generate_rates/PRIMAT-Main-gray.m's δχFM).
      _L_SD_FMCCR — Same with Coulomb × radiative factor.

    The finite-temperature CCRTh correction (Brown & Sawyer 2001; Phys. Rep.
    §III.F, Eqs. 107–113) is deliberately NOT summed here — see
    :func:`_correction_terms` for why. It is cached separately
    (``nTOp_thermal_<hash>.txt``) and recombined in :func:`RecomputeWeakRates`.

    The overall rate constant K is normalised via the neutron lifetime::

        K = 1 / (τ_n × Fn)     (Phys. Rep. Eq. 89–91)

    where Fn = ComputeFn(cfg) is the free-decay phase-space integral.

    Parameters
    ----------
    Tvec       : list [Tg_vec, Tnu_vec], both float arrays in Kelvin (as output
                 by PRIMAT._setup_background_and_cosmo).
    cfg        : PRIMATConfig instance.
    dFDneu_func: callable or None.
        If provided, adds the spectral-distortion correction _L_SD.  Signature::

            dFDneu_func(en, x, znu, sgnq) → float

        where en = E/mₑ, x = mₑ/(kB Tγ), znu = mₑ/(kB Tν), sgnq = ±1.
        Must encode the sign convention for blocking factors (en < 0), as
        described in PRIMAT._setup_background_and_cosmo.
    dFDneu_moments: dict or None.
        If provided (only in analytic-distortion mode -- see
        neutrino_history.AnalyticDistortion.dFDneu_moments), adds the SD-FM
        correction (_L_SD_FMCCR/_L_SD_FMNoCCR) on top of the SD correction
        above, when cfg.finite_mass_corrections is also True.

    Returns
    -------
    [T_all, frwrd, bkwrd] : list
        T_all  — 1-D float array, photon temperatures in Kelvin.
        frwrd  — 1-D float array, non-thermal Γ_{n→p}(T) in units of 1/τ_n.
        bkwrd  — 1-D float array, non-thermal Γ_{p→n}(T) in units of 1/τ_n.

    Example:
        >>> rates = ComputeWeakRates([Tg_vec, Tnu_vec], cfg)
        >>> T_K, lam_nTOp, lam_pTOn = rates
    """
    ctx = _build_rate_context(Tvec, cfg)

    # Single grid spanning the whole BBN temperature range (T_end -> T_start).
    # cfg.sampling_nTOp_per_decade is points per decade of T (formerly
    # sampling_nTOp was the *total* point count, and before that the
    # per-era count when the network used three separate HT/MT/LT grids).
    n_pts = n_points_per_decade(cfg.sampling_nTOp_per_decade, cfg.T_end, cfg.T_start)
    T_all = np.logspace(np.log10(cfg.T_end), np.log10(cfg.T_start), n_pts)

    # Each correction term is already vectorised over T_all, so
    # the forward / backward rates are just the element-wise sum of the term
    # arrays -- no Python loop over the grid.  The finite-temperature CCRTh
    # term is intentionally absent here (see _correction_terms): it is stored
    # separately and recombined in RecomputeWeakRates.
    def nTOp_rate_(sgnq):
        return sum(value for _, value in
                   _correction_terms(ctx, T_all, sgnq, dFDneu_func, dFDneu_moments))

    frwrd = nTOp_rate_(+1)
    bkwrd = nTOp_rate_(-1)

    # Normalise by the neutron-decay phase-space integral Fn so that the
    # returned values are in units of 1/tau_n (multiply by 1/tau_n to get
    # the actual rate in s^-1).  Values below 1e-28 are purely numerical
    # noise (the p->n rate at very low T is exp(-Q/T)-suppressed to
    # ~1e-40 and cancellation makes it alternate sign); clamp them to 0.
    # The threshold is applied to the RAW integral, before the 1/Fn division,
    # so in stored (1/tau_n) units it sits at 1e-28/Fn ~ 5.7e-29 -- irrelevant
    # at 28 orders of magnitude below the rates that matter, but worth stating
    # since everything else on this path is quoted in 1/tau_n. The C backend
    # clamps identically (weak_rates.c, `(f < 1e-28) ? 0.0 : f / Fn`).
    Fn    = ComputeFn(cfg)
    frwrd = np.where(frwrd < 1e-28, 0.0, frwrd / Fn)
    bkwrd = np.where(bkwrd < 1e-28, 0.0, bkwrd / Fn)

    return [T_all, frwrd, bkwrd]


# ---------------------------------------------------------------------------
# Load / dispatch interface
# ---------------------------------------------------------------------------

def InterpolateWeakRates(cfg):
    """Load n↔p weak rates from the hash-named cache file and return interpolants.

    Reads ``rates/weak/nTOp_<hash>.txt`` (three columns: T in Kelvin,
    Gamma_{n→p} in units of 1/tau_n, Gamma_{p→n} in units of 1/tau_n) where
    ``<hash>`` is the 16-hex fingerprint of the current configuration.
    Raises FileNotFoundError if the file does not exist (it has not been
    computed yet for this configuration).

    Args:
        cfg : PRIMATConfig instance (provides _resolved_data_dir).

    Returns:
        [frwrd, bkwrd] : two scipy interp1d objects (extrapolating), each mapping
                         T in Kelvin → rate in units of 1/tau_n.
    """
    fp_hash = fingerprint_hash(_weak_rate_fingerprint(cfg))
    # Overlay read: cache_dir (if set) first, then the shipped package copy.
    path    = resolve_cache_file(cfg, "weak", "nTOp_" + fp_hash + ".txt")
    tab     = np.loadtxt(path)
    # log10-log10 cubic (matching the C backend and the nuclear rate tables);
    # see _weak_rate_loglog_interp for why linear-space quadratic was replaced.
    frwrd   = _weak_rate_loglog_interp(tab[:, 0], tab[:, 1])
    bkwrd   = _weak_rate_loglog_interp(tab[:, 0], tab[:, 2])
    return [frwrd, bkwrd]


def RecomputeWeakRates(Tvec, cfg, dFDneu_func=None, dFDneu_moments=None):
    """Load the n<->p weak-rate tables from the fingerprinted cache, or recompute.

    Implements the cache-loading logic for the *non-thermal* rate (the
    thermal CCRTh correction is always handled separately, see step 5):

    1. Compute the fingerprint hash of the current configuration
       (:func:`_weak_rate_fingerprint`).
    2. If `cfg.weak_rate_cache` is True and `rates/weak/nTOp_<hash>.txt`
       exists, load and interpolate it directly (cheap: no integration at
       all).  The fingerprint is enforced by the *filename*: a different
       configuration hashes to a different name and therefore misses.  The
       `# fingerprint_hash:`/`# fingerprint:` header lines inside the file are
       written for human/migration use and are never re-read by the loader.
    3. Otherwise call :func:`ComputeWeakRates` to recompute from scratch
       (~2 s).  If `cfg.save_nTOp` is True, save the new data and the current
       fingerprint header to that file.
    4. **Forced recompute**: if `cfg.spectral_distortions and
       cfg.analytic_distortions`, the cache is bypassed entirely (never
       loaded, never written).  Analytic distortions are continuous knobs
       (`y_SZ`, `y_gray`) typically scanned point-by-point in an MCMC;
       caching them would write one file per parameter point and pollute
       rates/weak/.  The same rule applies to any future user-supplied
       `dFDneu_func` that cannot be fingerprinted.
    5. Build the finite-temperature CCRTh correction with
       :func:`_thermal_correction_interpolants` (its own
       ``nTOp_thermal_<hash>.txt`` cache) and add it to the non-thermal
       interpolant from step 2/3.  The returned rate is the sum of the two.

    ``cfg.tau_n_normalization``/``cfg.tau_n`` do not enter the fingerprint:
    the stored rates are in units of 1/τ_n (Fn already applied inside
    ComputeWeakRates), so they need only multiplying by 1/tau_n after loading
    — the cached values themselves are tau_n-independent.

    The hash is embedded in the filename (``nTOp_<hash>.txt``), so different
    configurations coexist in ``rates/weak/`` without overwriting each other.
    ``cfg.save_nTOp`` defaults to True: every newly computed configuration is
    saved automatically so subsequent runs reuse it without recomputing.

    Parameters
    ----------
    Tvec        : [Tg_vec, Tnu_vec]  (arrays in MeV)
    cfg         : PRIMATConfig
    dFDneu_func : callable or None — spectral-distortion correction function;
                  forwarded to ComputeWeakRates on a cache miss.
    dFDneu_moments : dict or None — energy-moment functions of dFDneu_func
                  (analytic-distortion mode only); forwarded to
                  ComputeWeakRates on a cache miss to enable the SD-FM term.

    Returns
    -------
    [frwrd, bkwrd] : two callables (n->p and p->n) T[K] -> rate in units of
    1/τ_n, covering the whole BBN temperature range.  Each is the sum of the
    stored non-thermal interpolant (Born+FM+CCR+SD) and the separately-cached
    finite-temperature correction (CCRTh), so multiplying by 1/τ_n gives the
    full physical rate in s⁻¹.
    """
    forced_recompute = cfg.spectral_distortions and cfg.analytic_distortions

    fp       = _weak_rate_fingerprint(cfg)
    fp_hash  = fingerprint_hash(fp)
    fname    = "nTOp_" + fp_hash + ".txt"
    # Overlay read path (cache_dir first, else shipped copy); the write path
    # below always targets the writable base (cache_dir if set, else package).
    path     = resolve_cache_file(cfg, "weak", fname)
    write_path = os.path.join(cache_write_dir(cfg, "weak"), fname)

    # ---- Non-thermal rate (Born+FM+CCR+SD): load from cache or recompute. ----
    nonthermal = None
    if not forced_recompute and cfg.weak_rate_cache and os.path.exists(path):
        nonthermal = InterpolateWeakRates(cfg)
        if cfg.verbose:
            print("[weak-py] background n<->p weak rates: loaded from cache.")

    if nonthermal is None:
        if cfg.verbose and not forced_recompute and cfg.weak_rate_cache:
            print("[weak-py] Recomputing n<->p weak rates (no cache for this configuration).")
        T_all, frwrd, bkwrd = ComputeWeakRates(Tvec, cfg, dFDneu_func=dFDneu_func,
                                                dFDneu_moments=dFDneu_moments)

        if not forced_recompute and cfg.save_nTOp:
            # write_cache_with_fingerprint creates the dir and degrades to a
            # warning (returning False) on a read-only install -- never fatal.
            write_cache_with_fingerprint(
                write_path, fp, [T_all, frwrd, bkwrd],
                col_header="T[K] Gamma_nTOp[1/tau_n] Gamma_pTOn[1/tau_n]")

        # Same log10-log10 cubic scheme as the cache-load path above (and the
        # C backend), so a freshly recomputed run and a cache-reloaded run of
        # the same configuration interpolate identically.
        nonthermal = [_weak_rate_loglog_interp(T_all, frwrd),
                      _weak_rate_loglog_interp(T_all, bkwrd)]

    # ---- Thermal correction (CCRTh): cached separately, recombined here. ----
    # Both halves are in units of 1/τ_n, so the physical rate the solver uses
    # is simply their sum.  Returning thin closures (rather than re-tabulating)
    # keeps the two caches independent: a config that differs only in a
    # non-thermal flag reuses the same thermal table, and vice versa.
    nTOp_thermal, pTOn_thermal = _thermal_correction_interpolants(Tvec, cfg)
    frwrd_nt, bkwrd_nt = nonthermal

    # Original array-safe closures (kept as the fallback / array-input path).
    frwrd_orig = lambda T: frwrd_nt(T) + nTOp_thermal(T)
    bkwrd_orig = lambda T: bkwrd_nt(T) + pTOn_thermal(T)

    # Wrap in JIT scalar evaluators for the hot solver path. FastWeakRate
    # self-checks against the closures above at build time and transparently
    # falls back to them (per channel) if the JIT path ever disagrees, so this
    # can only speed things up, never change the numbers.
    from .fast_eval import FastWeakRate
    from .corrections import _T_CCRTH_MIN
    frwrd = FastWeakRate(frwrd_nt, nTOp_thermal, _T_CCRTH_MIN, frwrd_orig)
    bkwrd = FastWeakRate(bkwrd_nt, pTOn_thermal, _T_CCRTH_MIN, bkwrd_orig)
    return [frwrd, bkwrd]
