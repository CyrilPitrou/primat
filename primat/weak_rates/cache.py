# -*- coding: utf-8 -*-
"""
weak_rates.cache — fingerprinted on-disk cache for the n<->p weak-rate tables
===============================================================================

Fingerprint dicts and the log-spaced grid-density helper shared by
weak_rates.api (ComputeWeakRates/InterpolateWeakRates/RecomputeWeakRates) and
weak_rates.corrections (the CCRTh thermal-correction cache).  See
weak_rates.api for the cache-file layout and invalidation policy.
"""

import os

import numpy as np

from ..cache_utils import constants_hash, fingerprint_hash, resolve_cache_file

__all__ = ['WEAK_RATE_FORMAT_VERSION', '_WEAK_RATE_BG_FIELDS', '_THERMAL_BG_FIELDS',
           'n_points_per_decade', '_thermal_fingerprint', '_weak_rate_fingerprint',
           'thermal_cache_exists']

# ---------------------------------------------------------------------------
# Fingerprinted cache for the n<->p weak-rate tables
# ---------------------------------------------------------------------------
# Bump this whenever a code change alters the *numerical content* of the
# cached files for a fixed configuration (new physics term, changed formula,
# different file layout, ...).  Bumping it invalidates every existing cache
# file regardless of its fingerprint.
#
# v1: forward and backward rates stored together in nTOp_<hash>.txt (hash in
# filename, rates in units of 1/tau_n, clamped below 1e-28 to zero).
# Fingerprints simplified: thermal uses only the T range, incomplete-decoupling
# flag, and NEVO file selection; weak-rate drops sampling_temperature_per_decade,
# nevo_grid_file, external_scale_factor, thermal_corrections and
# thermal_fingerprint_hash.  tau_n_flag renamed to tau_n_normalization.
# v2: sampling_nTOp/sampling_nTOp_thermal (total grid points) replaced by
# sampling_nTOp_per_decade/sampling_nTOp_thermal_per_decade (points per decade
# of T), so the grid density now stays fixed when T_end_MeV changes the span.
# v3: the n<->p rate table nTOp_<hash>.txt now stores ONLY the non-thermal
# rate (Born + finite-mass + CCR + spectral-distortion); the finite-temperature
# radiative correction (CCRTh) is kept in its own nTOp_thermal_<hash>.txt and
# recombined at point of use (RecomputeWeakRates), matching the fingerprint
# which never included thermal_corrections.  The CCRTh table content also
# changed: the n->p direction is now clamped to 0 below ~10^8.2 K (see
# _L_CCRTh_compute) to remove a spurious infrared-divergent bremsstrahlung
# residual.  Both changes invalidate every v2 cache file.
# v4: the v2 and v3 content changes above were documented but the constant was
# never actually bumped past 1 -- so neither of them ever invalidated anything,
# and a pre-v3 file (whose nTOp_<hash>.txt still *included* CCRTh, and whose
# nTOp_thermal_<hash>.txt was unclamped below 10^8.2 K) was silently reloaded
# whenever its other fingerprint fields happened to match.  v4 pays that debt
# and bundles it with three new fingerprint fields:
#   * munuOverTnu (effective xi_e) added to the THERMAL fingerprint -- the
#     thermal integrand's Chitilde carries exp(-sgnq*xi_nu), so the CCRTh table
#     was being shared between runs with different neutrino degeneracies
#     (measured: xi_e = 0.3 moves CCRTh by ~4e-3 of the base rate at 1e10 K);
#   * nevo_grid_file added to the weak-rate fingerprint, which already had
#     nevo_spectral_file -- the two jointly define the tabulated distortion fed
#     to _L_SD, so keying on only one of them was inconsistent;
#   * sampling_temperature_per_decade added to the weak-rate fingerprint: it
#     sets the density of the Tg grid behind the LINEAR T_nu(T_gamma)
#     interpolant every rate integrand reads (corrections._build_rate_context),
#     and coarsening it from the default 600 to 40 points/decade moves the
#     rates by up to ~1e-3 (to ~1e-5 at the default) -- far above the D/H
#     regression tolerance, so it cannot be treated as cache-neutral.
# Every v1..v3 cache file is invalidated by this bump; the shipped tables were
# re-keyed in place (same numbers, new hash-named filenames and headers).
# v5: no *physical constant* was keyed by any fingerprint. Both tables here are
# explicit functions of primat.constants.CONST -- the n<->p rate integrands read
# me, alphaem, mn, mp, gA, Vud, radproton, kappa_n/p and GF, and the CCRTh
# correction is itself O(alphaem) -- yet none of those 26 values appeared in the
# hashed dict. Editing one of them in constants.py (or in the C mirror,
# primat-c/src/constants.c) therefore reloaded rates computed with the OLD
# value: not a slow run, a silently wrong answer, with the stale file still
# advertising a matching fingerprint. v5 adds a single `constants_hash` field --
# the hash of the whole constants struct (cache_utils.constants_hash; the C
# backend computes the identical value via cpr_constants_hash) -- to BOTH
# fingerprints below, so any constants edit invalidates every weak-rate and
# thermal cache file at once. The same field is added in the same change to the
# other two fingerprinted caches (electron thermo and, newly fingerprinted, the
# QED pressure tables), so that one re-key covers all of them.
# As in v4, the shipped tables were re-keyed in place: numbers untouched, only
# the hash-named filenames and the header lines changed.
WEAK_RATE_FORMAT_VERSION = 5

# Config fields entering the weak-rate fingerprint (nTOp_<hash>.txt).
# DeltaNeff is deliberately NOT listed: it only shifts the time-temperature
# relation Tg(t) and does not affect the rate integrand at fixed Tg (in decoupling approximation).
# In principle if we consider a DeltaNeff with incomplete decoupling we must also consider the associated NEVO file.
# We need to review the interplay between NEVO and primat.
# Note  that spectral distortions and incomplete decoupling effects are expected to have a small effect on weak rates.
_WEAK_RATE_BG_FIELDS = [
    "radiative_corrections",
    "finite_mass_corrections",
    # NOTE: the neutrino chemical potential enters the weak-rate fingerprint via
    # the EFFECTIVE ξ_e (cfg.xi_nu_e), added under the historical key name
    # "munuOverTnu" in _weak_rate_fingerprint below -- NOT listed here. This keeps
    # the default-run hash byte-identical to the previous single-xi fingerprint
    # (xi_nu_e == munuOverTnu whenever munuOverTnu_e is unset), so the shipped
    # data/cache_plasma_weak/weak/ caches stay valid, while a per-flavour
    # munuOverTnu_e override correctly
    # produces a distinct cache. ξ_μ/ξ_τ gravitate only and never enter here.
    "QED_corrections",
    "incomplete_decoupling",
    "spectral_distortions",
    "analytic_distortions",
    "y_SZ",
    "y_gray",
    "T_start_cosmo_MeV",
    "T_end_MeV",
    "sampling_nTOp_per_decade",
    # Density of the Tg grid the caller hands to ComputeWeakRates. It does not
    # set the rate table's own grid (that is sampling_nTOp_per_decade), but it
    # does set the node spacing of the LINEAR T_nu(T_gamma) interpolant built in
    # corrections._build_rate_context, which every rate integrand evaluates.
    # Measured on the default NEVO history (rel. change of Gamma_{n->p} against
    # a 2000-points/decade reference): 2e-3 at 10 points/decade, 1e-3 at 40,
    # 1.4e-4 at 200, 1.4e-5 at the default 600. That is orders of magnitude
    # above the +-3e-9 D/H regression tolerance, so a run that changes this
    # knob must NOT reuse another run's table.
    "sampling_temperature_per_decade",
    "nevo_file",
    "nevo_spectral_file",
    # nevo_grid_file completes the pair with nevo_spectral_file: the spectral
    # table supplies the distortion columns and the grid file the energy nodes
    # they are sampled on (neutrino_history.NEVOTable), so both are needed to
    # pin down the dFDneu_func that _L_SD/_L_SD_CCR integrate. Keying on only
    # one of them let a custom grid of the same length silently reuse the
    # shipped grid's rates.
    "nevo_grid_file",
    "nevo_file_prefix",
]

# Config fields entering the thermal-correction fingerprint
# (nTOp_thermal_<hash>.txt).  The temperature range and sampling, the neutrino
# decoupling mode (with or without QED corrections), the NEVO thermo table
# selection, and the electron-neutrino degeneracy matter for
# the double-integral over (E, k) that defines the finite-temperature
# radiative correction.
# When improving the interpolay with NEVO this could be improved.
#
# The effective xi_e is added under the historical "munuOverTnu" key by
# _thermal_fingerprint below (not listed here), exactly as in
# _weak_rate_fingerprint -- see that function for why the key name and the
# per-flavour fallback are what they are. It belongs in this fingerprint
# because the thermal integrands' Chitilde carries an explicit
# exp(znu*(en - sgnq*q) - sgnq*xi_nu) neutrino occupation
# (corrections._ccrth_IPENCCRT / _ccrth_IPENCCRDiffBremsstrahlung /
# _ccrth_C2dE1dE2): at xi_e = 0.3 the CCRTh term moves by ~4e-3 of the base
# CCR rate at 1e10 K (n->p) -- so before it was keyed here, every degenerate-BBN
# run silently reused the xi_e = 0 table (and, on a cold cache, wrote its own
# xi_e-specific numbers under the filename standard runs then load).
#
# PROVENANCE CAVEAT (thermal cache only). The CCRTh table is a Monte-Carlo
# (vegas) estimate, and the shipped nTOp_thermal_<hash>.txt files were produced
# by a Python run whose vegas was UNSEEDED. Python's vegas is now seeded
# deterministically (weak_rates/corrections.py's _vegas_rng, mirroring
# primat-c/src/weak_rates.c's th_vegas_seed), so from here on a Python
# recompute reproduces *itself* exactly.
#
# It does not reproduce the *shipped file* exactly, and cannot: the shipped
# numbers came from a random stream that no longer exists. So a cache hit and a
# fresh recompute of the same configuration -- same hash, same fingerprint --
# can differ by the MC noise floor. That is inherent to caching a Monte-Carlo
# integral, not a defect in the fingerprint: the fingerprint's job is to say
# "this file describes THIS configuration", which remains true.
#
# Measured (2026-08-04), recomputing the default shipped table
# (nTOp_thermal_c7da75afa7c0bf3b.txt, 57 rows) with the seeded vegas:
#
#   * two successive recomputes are now BYTE-IDENTICAL (max column difference
#     exactly 0, identical D/H to all digits) -- the property seeding buys;
#   * against the shipped (unseeded) table the CCRTh term itself moves by up to
#     3.8e-2 relative on L_nTOpCCRTh and 1.2e-2 on L_pTOnCCRTh, with medians of
#     9.7e-4 and 3.6e-4. Those percentages look large only because they are
#     relative to CCRTh, which is itself a ~1e-3 correction to the total n<->p
#     rate -- the worst point is ~4e-5 of the rate;
#   * the resulting shift in the headline observable is D/H 2.43589845e-05 ->
#     2.43589505e-05, i.e. |dD/H| = 3.4e-11, a factor 88 inside the +-3e-9 D/H
#     regression tolerance.
#
# So regenerating the shipped thermal tables is NOT required for correctness,
# which is why this change deliberately leaves them alone: doing so would cost
# a multi-minute vegas run per shipped configuration and churn several tracked
# data files, to move an observable by 1% of its own regression tolerance.
#
# The `provenance:` header line (cache_utils.write_cache_with_fingerprint)
# exists to record which backend and algorithm produced the numbers on disk
# right now; it is deliberately NOT part of the hash, precisely so this
# MC-noise difference never gates a cache hit/miss decision.
#
# T_end_MeV is deliberately NOT listed: the integral is clamped to exactly 0
# below ~10**8.2 K regardless of cfg.T_end (see corrections._T_CCRTH_MIN /
# _L_CCRTh_compute), and the cache grid is now built down to that fixed floor
# rather than down to cfg.T_end (corrections._L_CCRTh_interpolants), so the
# integral never depends on T_end_MeV in the first place. Including it here
# only forced spurious cache misses -- and the multi-minute vegas recompute
# that goes with them -- whenever a run changed T_end_MeV alone.
#
# sampling_temperature_per_decade is likewise deliberately NOT listed, even
# though it IS in _WEAK_RATE_BG_FIELDS: it acts on both tables through the same
# T_nu(T_gamma) interpolant, but CCRTh is itself only a ~1e-3 correction to the
# rate, so the ~1e-5 interpolant error at the default density propagates to
# ~1e-8 of the total rate -- far below any tolerance, and not worth forcing a
# multi-minute vegas recompute for. nevo_spectral_file/nevo_grid_file are
# absent for a stronger reason: the thermal integrands use the plain
# Fermi-Dirac occupation, never the distortion.
_THERMAL_BG_FIELDS = [
    "T_start_cosmo_MeV",
    "sampling_nTOp_thermal_per_decade",
    "QED_corrections",
    "incomplete_decoupling",
    "nevo_file",
    "nevo_file_prefix"
]


def n_points_per_decade(per_decade, T_lo, T_hi):
    """Number of log-spaced grid points spanning [T_lo, T_hi] at a fixed
    density of ``per_decade`` points per decade of T.

    Used so that ``sampling_nTOp_per_decade``/``sampling_nTOp_thermal_per_decade``
    keep a constant grid resolution even if ``T_end_MeV`` (and hence the
    number of decades spanned) changes, unlike the old total-point-count
    parametrisation.

    Args:
        per_decade: float, desired points per decade of T.
        T_lo, T_hi: float, grid endpoints [K], T_hi > T_lo.

    Returns:
        int, number of points (at least 2).
    """
    decades = np.log10(T_hi / T_lo)
    return max(2, int(round(per_decade * decades)))


def _thermal_fingerprint(cfg):
    """Fingerprint dict for the thermal radiative-correction cache file
    ``nTOp_thermal_<hash>.txt``.

    Only the fields that actually affect the finite-temperature double
    integral (Brown & Sawyer 2001) are included: the temperature integration
    range, the neutrino-to-photon temperature ratio T_ν(T_γ) (fixed by the
    NEVO table selection), the electron-neutrino degeneracy ξ_e (which enters
    the integrands' neutrino occupation directly), and the thermal-correction
    grid density.  See :data:`_THERMAL_BG_FIELDS` for what is deliberately
    left out and why.

    Args:
        cfg: PRIMATConfig instance.

    Returns:
        dict, JSON-serialisable.
    """
    fp = {"format_version": WEAK_RATE_FORMAT_VERSION,
          "sampling_nTOp_thermal_per_decade": cfg.sampling_nTOp_thermal_per_decade,
          # Physical constants (v5). The finite-temperature radiative
          # correction is itself O(alphaem) and its integrands carry me, so the
          # table is a direct function of the constants struct; hashing the whole
          # struct (cache_utils.constants_hash) rather than just those two is
          # deliberate -- see that function's docstring for why over-invalidating
          # is the safe side of this trade.
          "constants_hash": constants_hash()}
    for key in _THERMAL_BG_FIELDS:
        fp[key] = getattr(cfg, key)
    # Effective ξ_e under the same historical "munuOverTnu" key as
    # _weak_rate_fingerprint, so the two fingerprints name it identically
    # (ξ_μ/ξ_τ gravitate only and never reach these integrands).
    fp["munuOverTnu"] = cfg.xi_nu_e
    return fp


def thermal_cache_exists(cfg):
    """Whether ``cfg``'s CCRTh thermal-correction cache file is already on disk.

    ``True`` means :func:`weak_rates.corrections._L_CCRTh_interpolants` will
    load the fingerprinted ``nTOp_thermal_<hash>.txt`` file instead of
    running its (multi-minute, vegas-based) Monte-Carlo integration --
    callers that only need to know "is this about to be slow" (e.g. the GUI's
    progress message in ``gui/app.py``) should check this rather than just
    ``cfg.thermal_corrections``, since that flag alone says nothing about
    whether a cache hit is coming.

    Args:
        cfg: PRIMATConfig instance.

    Returns:
        bool.
    """
    # Overlay read: resolve_cache_file returns an existing cache_dir/shipped
    # copy if present, else the (non-existent) write path -- so os.path.exists
    # correctly reports a hit only when the thermal cache is actually reachable.
    path = resolve_cache_file(
        cfg, "weak",
        "nTOp_thermal_" + fingerprint_hash(_thermal_fingerprint(cfg)) + ".txt")
    return os.path.exists(path)


def _weak_rate_fingerprint(cfg):
    """Fingerprint dict for the n<->p weak-rate cache file ``nTOp_<hash>.txt``.

    ``cfg.tau_n_normalization``/``cfg.tau_n`` are deliberately excluded: the
    stored rates are in units of 1/τ_n (Fn already applied inside
    :func:`ComputeWeakRates`), so only 1/tau_n needs multiplying after
    loading — the cached values themselves are tau_n-independent.

    The thermal-correction cache has its own hash-named file and is not
    folded in here: the two caches are independent, and ``thermal_corrections``
    itself does not affect the stored non-thermal rates.

    Args:
        cfg: PRIMATConfig instance.

    Returns:
        dict, JSON-serialisable; pass to :func:`fingerprint_hash` for the hash.
    """
    fp = {"format_version":          WEAK_RATE_FORMAT_VERSION,
          "sampling_nTOp_per_decade": cfg.sampling_nTOp_per_decade,
          "radiative_corrections":   cfg.radiative_corrections,
          "finite_mass_corrections": cfg.finite_mass_corrections,
          # Physical constants (v5). The rate integrands read me, alphaem, mn,
          # mp, gA, Vud, radproton, kappa_n/p and GF directly, so the stored
          # rates are a function of the constants struct; before v5 none of them
          # was keyed, and editing one silently reused the old table. Hashing the
          # whole struct rather than that curated list is deliberate -- see
          # cache_utils.constants_hash.
          "constants_hash":          constants_hash()}
    for key in _WEAK_RATE_BG_FIELDS:
        fp[key] = getattr(cfg, key)
    # Neutrino chemical potential in the weak rates: only the electron-neutrino
    # ξ_e matters (n <-> p + e + nu_e). Store the EFFECTIVE ξ_e (per-flavour
    # override munuOverTnu_e, else the common munuOverTnu) under the historical
    # "munuOverTnu" key so a default run (munuOverTnu_e unset) hashes exactly as
    # before and keeps hitting the shipped data/cache_plasma_weak/weak/ caches;
    # the C backend
    # mirrors this (cpr_weak_rate_fingerprint, cache.c). ξ_μ/ξ_τ are omitted:
    # they gravitate only and do not touch the weak rates.
    fp["munuOverTnu"] = cfg.xi_nu_e
    # custom_background mode takes the Tg grid behind the T_nu(T_gamma)
    # interpolant from the *table file's own* T range (CustomBackground.
    # _setup_neutrino_history), NOT from T_start_cosmo_MeV/T_end_MeV -- so none
    # of the range/density fields above distinguish one custom table from
    # another, and two different custom backgrounds silently shared a single
    # cached nTOp table. Keyed on the path here, exactly as nevo_file is (an
    # in-place edit of the same path is still not caught -- same caveat as
    # there).
    #
    # Added CONDITIONALLY so a run without custom_background hashes exactly as
    # before and keeps hitting the shipped data/cache_plasma_weak/weak/ caches
    # -- the same trick as munuOverTnu above. Mirrored in C by
    # cpr_weak_rate_fingerprint (cache.c).
    if getattr(cfg, "custom_background", None) is not None:
        fp["custom_background"] = cfg.custom_background
    return fp

