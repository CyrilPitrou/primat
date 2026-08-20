# -*- coding: utf-8 -*-
"""
background.py
==============
Cosmological-background component for primat ("Class 1" of the
``PRIMAT`` split).

A *background* encapsulates everything the nuclear-network integration
(:class:`primat.nuclear_network.NuclearNetwork`, "Class 2") needs about the
expanding Universe.  The interface deliberately is **minimal**: only the
``T_of_t``/``t_of_T`` time<->temperature relations, the baryon mass density
``rhoB_BBN(t)`` *as a function of cosmic time*, and the (already normalised)
n<->p weak rates ``weak_nTOp_frwrd(T)``/``weak_nTOp_bkwrd(T)`` are compulsory.
Everything else (scale factor ``a``, the Hubble rate, individual neutrino
temperatures, ``N_eff``, ``Omega_nu``, the NEVO heating function, ...) is
optional, with safe ``None``/``NotImplementedError`` defaults in
:class:`Background`, so a minimal custom background (e.g. an externally
supplied ``T(t)``/``rho_B(t)`` table with no neutrino-sector model at all) can
still drive :class:`primat.nuclear_network.NuclearNetwork`.

:class:`Background` is the interface; :class:`StandardBackground` is today's
(and so far only) full implementation: NEVO non-instantaneous decoupling or
instantaneous decoupling, selected via ``cfg.incomplete_decoupling``, with the
neutrino sector itself delegated to
:func:`primat.neutrino_history.make_neutrino_history` (already a pluggable
seam).  Following the style of :mod:`primat.neutrino_history`, this is a
plain base class whose compulsory interface methods raise
``NotImplementedError`` -- *not* ``abc.ABC``/``@abstractmethod`` -- so a user
can subclass :class:`Background` for a custom cosmology (e.g. non-standard
expansion history) and hand an instance to a future
``PRIMAT(..., background=...)`` hook.
"""

import io
import os
import sys
import time
from collections import OrderedDict

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import zeta

from . import weak_rates as primat_weak_rates
from .neutrino_history import make_neutrino_history

__all__ = ["Background", "StandardBackground", "CustomBackground"]

# How many sub-intervals each grid interval is split into in the
# external_scale_factor mode -- both for the T(a) inverse and for the t(a)
# output grid t_of_T/T_of_t are built on.  a(T) is an algebraic read of the
# NEVO `x` column in that mode, so refining costs one table read per node and
# no extra ODE work, while both interpolants' error is O(h^2): 8 buys a factor
# ~64.  Mirrored by CPR_EXTERNAL_A_REFINE in primat-c/src/background.c.
_EXTERNAL_A_REFINE = 8

# Noise floor for the n<->p weak rates, in raw (1/tau_n) units -- i.e. the
# units returned by weak_nTOp_{frwrd,bkwrd}_raw, *before* multiplying by
# _norm_weak_rates.  Matches the threshold already used inside
# ComputeWeakRates (weak_rates.py) when it builds the cached tables: below
# this scale the rate is exp(-Q/T)-suppressed and the sum of phase-space
# correction terms has lost all significant digits to cancellation,
# alternating sign around zero.  Re-applied here (rather than trusted to the
# cache alone) because the quadratic interpolant built by
# InterpolateWeakRates can overshoot slightly negative between cached nodes
# even when every node value is >= 0.
_WEAK_RATE_FLOOR = 1e-28


def _clamp_raw_weak_rate(rate):
    """Zero out negative/sub-noise raw n<->p weak-rate values.

    Args:
        rate: array-like, raw rate in 1/tau_n units (as returned by
            weak_nTOp_{frwrd,bkwrd}_raw).

    Returns:
        Same shape as ``rate``, with any value below :data:`_WEAK_RATE_FLOOR`
        (this includes all negative values) replaced by 0.0.
    """
    if np.ndim(rate) == 0:      # hot path: one scalar per RHS call
        return 0.0 if rate < _WEAK_RATE_FLOOR else rate
    return np.where(rate < _WEAK_RATE_FLOOR, 0.0, rate)


def _scalar_linear_eval(interp):
    """Fast scalar evaluator for a linear ``interp1d``, or ``None``.

    ``interp1d.__call__`` costs ~10 µs of input validation per query, and
    ``T_of_t``/``a_of_t`` are queried once per nuclear-network RHS call.  This
    re-expresses the *evaluation only*, with the same two-node formula scipy
    uses -- the convex combination of ``_call_linear`` when the interpolant
    extrapolates, ``np.interp`` when it clamps to a constant ``fill_value``.
    Which one applies is settled by a build-time check against ``interp``
    itself on the nodes, the interval midpoints and both extrapolation sides;
    unless every digit matches -- or that check cannot even be run -- ``None``
    is returned and the caller keeps the scipy path.

    Args:
        interp: a ``scipy.interpolate.interp1d`` built with ``kind='linear'``.

    Returns:
        Callable float -> float, or ``None`` if neither form reproduces
        ``interp`` exactly.
    """
    try:
        xs = np.ascontiguousarray(interp.x, dtype=float)
        ys = np.ascontiguousarray(np.asarray(interp._y, dtype=float).ravel())
    except Exception:
        return None
    n = xs.size
    if n < 2 or ys.size != n:
        return None
    searchsorted = np.searchsorted

    def extrapolating(q):
        i = int(searchsorted(xs, q))
        if i < 1:
            i = 1
        elif i > n - 1:
            i = n - 1
        x_lo, x_hi = xs[i - 1], xs[i]
        dx = x_hi - x_lo
        return ((q - x_lo) / dx) * ys[i] + ((x_hi - q) / dx) * ys[i - 1]

    def clamping(q):
        return float(np.interp(q, xs, ys))

    # Probe the nodes (subsampled -- a background grid has thousands),
    # the midpoints between them, and one point off each end.
    step = max(1, n // 128)
    nodes = xs[::step]
    probe = np.concatenate([nodes, 0.5 * (nodes[:-1] + nodes[1:]),
                            [xs[0] - (xs[1] - xs[0]),
                             xs[-1] + (xs[-1] - xs[-2])]])
    try:
        reference = np.asarray(interp(probe), dtype=float).ravel()
    except Exception:
        # The probe is the check; a scipy whose evaluator has moved far enough
        # to raise on it has to be declined, not propagated to the caller.
        return None
    for candidate in (extrapolating, clamping):
        if all(float(candidate(q)) == r for q, r in zip(probe, reference)):
            return candidate
    return None


def _linear_with_powerlaw_tail(x, y):
    """Linear interpolant on ``(x, y)`` whose upper tail is a power law.

    Linear in ``x`` inside the grid and below it, exactly as
    ``interp1d(kind='linear', fill_value="extrapolate")``; above the last node
    it continues the power law through the last two nodes,
    ``y = y[-1] (x/x[-1])^p`` with ``p = dln y / dln x`` there.  Linear
    extrapolation of a decaying ``y`` crosses zero a finite distance past the
    end -- for ``T(a)`` at about ``2 a_end`` -- and a negative temperature is
    worse than an inaccurate one.

    Both arrays must be strictly positive and ``x`` strictly increasing.
    Mirrored by ``cpr_bg_T_of_a`` in ``primat-c/src/background.c``.

    Args:
        x, y : 1-D arrays of the same length (e.g. scale factor, T in MeV).

    Returns:
        callable ``f(x)`` accepting a scalar or array.

    Example:
        >>> f = _linear_with_powerlaw_tail(np.array([1., 2.]), np.array([2., 1.]))
        >>> float(f(4.0))          # 1 * (4/2)^-1, not the linear 0.0
        0.5
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    lin = interp1d(x, y, bounds_error=False, fill_value="extrapolate",
                   kind='linear')
    p = np.log(y[-1] / y[-2]) / np.log(x[-1] / x[-2])

    def _eval(q):
        qa = np.asarray(q, dtype=float)
        out = np.where(qa > x[-1],
                       y[-1] * np.power(np.maximum(qa, 1e-300) / x[-1], p),
                       lin(qa))
        return out if out.ndim else float(out)

    return _eval


def _loglog_interp1d(x, y):
    """Linear interpolant in (log x, log y), evaluated back in linear space.

    For the background's time/temperature relations both variables are
    strictly positive and vary as near-perfect power laws over many decades
    (t ∝ T⁻² in radiation domination), so a straight line in log-log space
    tracks them far better than one in linear space on the *same* nodes.
    Measured on the default ``sampling_temperature_per_decade = 600`` grid,
    switching ``t_of_T`` from linear to log-log cut its worst-case error from
    1.35e-05 to 8.8e-07, and its typical error to ~1e-09 -- gains of 10³-10⁵×
    at most temperatures. The residual peaks around 0.1-0.5 MeV, the e⁺e⁻
    annihilation era, which is exactly where t(T) departs most from a power
    law; nothing is lost there relative to the linear interpolant.

    This matters beyond t_of_T: ``T_of_t`` is evaluated on *every* nuclear
    network right-hand-side call, so its interpolation error enters the ODE
    solution directly.

    ``a_of_T`` already used this scheme (``np.exp(interp1d(log T))``) and was
    correspondingly accurate; this helper simply makes that treatment shared
    and explicit. Kept in lockstep with the C backend's log-log branch in
    ``cpr_bg_t_of_T`` / ``cpr_bg_T_of_t``.

    Args:
        x: array-like, strictly positive, the abscissa nodes.
        y: array-like, strictly positive, the ordinate nodes.

    Returns:
        Callable mapping x -> y, array-safe, extrapolating as a straight line
        in log-log space (i.e. continuing the end power law, which is the
        physically right behaviour here rather than the linear interpolant's
        straight line in (x, y) -- see :meth:`Background.t_of_T`'s note on
        out-of-range queries returning negative times).

    Example:
        >>> f = _loglog_interp1d([1.0, 2.0, 4.0], [1.0, 4.0, 16.0])  # y = x^2
        >>> float(f(3.0))                                # doctest: +ELLIPSIS
        9.0...
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    ln_interp = interp1d(np.log(x[order]), np.log(y[order]),
                         bounds_error=False, fill_value="extrapolate",
                         kind='linear')
    fast = _scalar_linear_eval(ln_interp)      # None if scipy's internals moved
    log, exp, ndim = np.log, np.exp, np.ndim

    def _f(q):
        if fast is not None and ndim(q) == 0:
            return exp(fast(log(q)))
        return exp(ln_interp(log(q)))

    return _f


class Background(object):
    """Minimal interface for the cosmological background consumed by
    ``NuclearNetwork``.

    ``__init__`` stores the constructor arguments common to every
    implementation and sets :attr:`has_scale_factor` to ``False`` (a
    subclass that builds ``a<->t<->T`` relations sets it to ``True``).

    Parameters
    ----------
    cfg : PRIMATConfig
        Run-time configuration (physical constants + flags).
    plasma : primat.plasma.Plasma
        Per-instance QED/electron-thermodynamics tables.  Always present
        (even a minimal background needs ``self.plasma.rho_g`` for
        :meth:`N_eff`).
    extra_rho : list of callable, optional
        Extra contributions to the total energy density entering the
        Friedmann equation.  Each element is a function
        ``rho(Tg) -> MeV^4`` of the photon temperature ``Tg`` [MeV].
        Stored (as a list, copied) on :attr:`extra_rho`; subclasses may
        append further components (e.g. Early Dark Energy) during their own
        setup.  A minimal background that has no Friedmann/Hubble machinery
        may simply ignore this list.

    Compulsory interface (raise ``NotImplementedError`` here)
    -----------------------------------------------------------
    * :meth:`T_of_t`, :meth:`t_of_T`  -- T_gamma(t) <-> t(T_gamma)
    * :meth:`rhoB_BBN`                -- rho_B(t) [g/cm^3], the prefactor for
      nuclear reaction rates (rate ~ rho_B)
    * :meth:`weak_nTOp_frwrd`, :meth:`weak_nTOp_bkwrd` -- normalised n<->p
      weak rates [s^-1] at photon temperature T [Kelvin]

    Optional interface (concrete defaults below)
    -----------------------------------------------
    * :meth:`a_of_T`/:meth:`T_of_a`/:meth:`a_of_t`/:meth:`t_of_a` -- scale
      factor relations; raise ``NotImplementedError`` unless
      :attr:`has_scale_factor` is ``True``.
    * :meth:`rho_nu_total_final` -- ``None`` unless the background tracks a
      neutrino sector.
    * :meth:`Omeganuh2_relnu`/:meth:`Omeganuh2_nrnu` -- ``None`` unless the
      background tracks a relic neutrino background.
    * :meth:`N_eff` -- concrete, generic formula (uses :attr:`plasma` only).
    * :meth:`write_time_evolution`/:meth:`_background_columns` -- concrete;
      the minimal background writes a two-column (``T``, ``t``) TSV.
    """

    def __init__(self, cfg, plasma, extra_rho=None):
        self.cfg = cfg
        self.plasma = plasma
        # Pluggable extra energy-density components (Early Dark Energy etc.);
        # copied so the caller's list is not mutated by subclass setup.
        self.extra_rho = list(extra_rho) if extra_rho is not None else []
        # Set True by a subclass that builds a<->t<->T relations
        # (StandardBackground); a minimal background has no scale factor.
        self.has_scale_factor = False

    # ======================================================================
    # Time <-> temperature (compulsory)
    # ======================================================================

    def T_of_t(self, t):
        """Photon temperature T_gamma(t) [MeV] at cosmic time ``t`` [s]
        (array-safe).

        Defined on the solved range only -- see :meth:`t_of_T` for what
        implementations do outside it."""
        raise NotImplementedError

    def t_of_T(self, T):
        """Cosmic time t(T_gamma) [s] at photon temperature ``T`` [MeV]
        (array-safe).

        **Valid range.** Implementations are only meaningful between
        ``cfg.T_end`` and ``cfg.T_start_cosmo``, the span over which the
        background was actually solved.  :class:`StandardBackground` builds
        this as a linear ``interp1d`` with ``fill_value="extrapolate"``, so a
        query outside that span is a straight-line continuation in ``(T, t)``
        -- not the ``t ∝ 1/T²`` of radiation domination, and far enough out it
        returns a *negative* time (with the default 0.001-40 MeV grid,
        ``t_of_T(80 MeV) ≈ -4.7e-04 s``).  Nothing inside primat queries
        outside the span -- the nuclear network starts at ``cfg.T_start``,
        comfortably inside -- but callers driving a background directly should
        stay in range or clamp.  Note this is a different question from the
        radiation-domination extrapolation *below* the NEVO table's lower edge
        (0.0141 MeV), which is inside the solved span and is handled correctly
        by both ``external_scale_factor`` modes."""
        raise NotImplementedError

    # ======================================================================
    # Scale factor (optional)
    # ======================================================================

    def a_of_T(self, T):
        """Scale factor a(T_gamma) (array-safe); requires
        :attr:`has_scale_factor`."""
        raise NotImplementedError

    def T_of_a(self, a):
        """Photon temperature T_gamma(a) [MeV] (array-safe); requires
        :attr:`has_scale_factor`."""
        raise NotImplementedError

    def a_of_t(self, t):
        """Scale factor a(t) (array-safe); requires
        :attr:`has_scale_factor`."""
        raise NotImplementedError

    def t_of_a(self, a):
        """Cosmic time t(a) [s] (array-safe); requires
        :attr:`has_scale_factor`."""
        raise NotImplementedError

    # ======================================================================
    # Baryon sector (compulsory)
    # ======================================================================

    def rhoB_BBN(self, t):
        """Baryon mass density rho_B(t) [g cm^-3] at cosmic time ``t`` [s]
        -- the prefactor for nuclear reaction rates (rate ~ rho_B)."""
        raise NotImplementedError

    # ======================================================================
    # n <-> p weak rates (normalised, compulsory)
    # ======================================================================

    def weak_nTOp_frwrd(self, T_K):
        """Normalised n -> p weak rate [s^-1] at photon temperature ``T_K``
        [Kelvin]."""
        raise NotImplementedError

    def weak_nTOp_bkwrd(self, T_K):
        """Normalised p -> n weak rate [s^-1] at photon temperature ``T_K``
        [Kelvin]."""
        raise NotImplementedError

    # ======================================================================
    # Derived cosmology (optional)
    # ======================================================================

    def Tnu_of_t(self, t):
        """Per-flavour neutrino temperature [MeV] at cosmic time ``t`` [s]
        (array-safe), or ``None`` if this background tracks no neutrino
        sector (the minimal background).

        A concrete subclass that tracks a neutrino sector returns
        ``{"e": Tnue, "mu": Tnumu, "tau": Tnutau}`` -- used by
        :meth:`_background_columns`/:meth:`time_evolution_text` and by
        :class:`primat.nuclear_network.NuclearNetwork`'s unified
        ``EvolutionResult``.
        """
        return None

    def rho_nu_total_final(self):
        """Final-time neutrino sector summary, or ``None``.

        A concrete subclass that tracks a neutrino sector returns
        ``(Tg_final, rho_nu_tot_final)`` -- the photon temperature [MeV] and
        the *total* neutrino energy density [MeV^4] (summed over flavours,
        plus any extra/spectral-distortion contributions) at the end of the
        integration.  Together with :meth:`N_eff` this gives
        ``Neff = N_eff(Tg_final, rho_nu_tot_final)``.  Returns ``None`` when
        no such information is available (the minimal background).
        """
        return None

    def Omeganuh2_relnu(self):
        """Omega_nu h^2 **per neutrino flavour** (nu + nubar), treating the
        relic neutrinos as relativistic today (massless-neutrino convention),
        or ``None`` if not tracked.

        Per flavour, not summed over the three: multiply by 3 to compare with
        the usual quoted total (Omega_nu h^2 ~ 1.7e-5 for Neff ~ 3.044).  The
        per-flavour convention matches :meth:`Omeganuh2_nrnu` below, whose
        reciprocal is the standard ~93 eV constant.
        """
        return None

    def Omeganuh2_nrnu(self):
        """Omega_nu h^2 **per neutrino flavour**, per unit neutrino mass
        [MeV^-1], treating the relic neutrinos as non-relativistic today
        (massive-neutrino convention), or ``None`` if not tracked.

        This is ``n_nu / rho_crit,100``, so ``1 / (Omeganuh2_nrnu * 1e-6)``
        recovers the familiar ``sum(m_nu) / 93.1 eV`` normalisation -- which is
        itself a per-flavour number density, hence the convention here.
        """
        return None

    def N_eff(self, Tg, rho_nu_tot):
        """Effective number of relativistic neutrino species.

        Generic formula, valid for any background that can supply the total
        neutrino energy density ``rho_nu_tot`` [MeV^4] at photon temperature
        ``Tg`` [MeV]:

            Neff = rho_nu_tot / rho_g(Tg) / ((7/8) (4/11)^(4/3))

        where ``rho_g(Tg)`` (from :attr:`plasma`, always present) is the
        photon energy density.  The ``(7/8)(4/11)^(4/3)`` factor converts one
        massless fermionic neutrino species (at its own temperature, energy
        density ``(7/8) rho_g(Tnu)``) into photon-temperature units via the
        standard instantaneous-decoupling ratio ``Tnu/Tg = (4/11)^(1/3)``, so
        that ``Neff = 3`` for three such species with no extra entropy
        injection.

        ``rho_nu_tot`` is whatever :meth:`rho_nu_total_final` (or an
        equivalent subclass-specific calculation) returns -- for
        :class:`StandardBackground` it already bundles the three
        flavour-dependent neutrino temperatures, any extra ΔNeff/extra_rho
        contributions to the neutrino sector, and the spectral-distortion
        extra energy density.
        """
        return rho_nu_tot / self.plasma.rho_g(Tg) / ((7. / 8.) * (4. / 11.) ** (4. / 3.))

    # ======================================================================
    # Background time-evolution output (optional, concrete default)
    # ======================================================================

    def _background_columns(self, t_out):
        """Return an ``OrderedDict`` of output columns evaluated on ``t_out``.

        The minimal background provides only ``T`` and ``t``; subclasses
        (e.g. :class:`StandardBackground`) extend this with ``a``, ``H``,
        individual neutrino temperatures, energy densities, etc.
        """
        return OrderedDict([("T [MeV]", self.T_of_t(t_out)), ("t [s]", t_out)])

    def time_evolution_text(self, n_points):
        """Return this background's time-evolution TSV as text, with no
        disk I/O.

        Same output grid and columns as :meth:`write_time_evolution`
        (``n_points`` log-spaced cosmic times from ``t_of_T(T_start_cosmo)``
        to ``t_of_T(T_end)``) -- factored out so a caller that wants the
        data in memory (e.g. ``primat-gui``'s download buttons) never
        needs a temporary file.
        """
        cfg = self.cfg
        T_start_cosmo = cfg.T_start_cosmo / cfg.MeV_to_Kelvin   # [MeV]
        T_end_MeV     = cfg.T_end         / cfg.MeV_to_Kelvin   # [MeV]
        t_lo = self.t_of_T(T_start_cosmo)
        t_hi = self.t_of_T(T_end_MeV)
        t_out = np.logspace(np.log10(t_lo), np.log10(t_hi), n_points)

        cols = self._background_columns(t_out)

        buf = io.StringIO()
        out_data = np.column_stack(list(cols.values()))
        out_header = "\t".join(cols.keys())
        np.savetxt(buf, out_data, delimiter='\t', header=out_header, comments='')
        return buf.getvalue()

    def write_time_evolution(self, path, n_points):
        """Write the background time evolution to ``path`` as a TSV.

        Columns are whatever :meth:`_background_columns` returns for this
        background -- at least ``T`` [MeV] and ``t`` [s];
        :class:`StandardBackground` adds ``a``, ``H``, the three flavour
        neutrino temperatures, the NEVO heating function (if available), and
        the plasma/neutrino/extra/total energy densities [MeV^4]. See
        :meth:`time_evolution_text` for the in-memory equivalent.

        Enabled by ``cfg.output_background_evolution=True``; the destination
        is ``cfg.output_background_file`` (relative paths resolve against the
        current working directory, like ``cfg.output_file``).
        """
        text = self.time_evolution_text(n_points)

        out_path = os.path.abspath(path)
        out_dir = os.path.dirname(out_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(out_path, 'w') as f:
            f.write(text)

        # Always announce: written only on explicit request
        # (output_background_evolution=True), like the nuclear-network TSV.
        n_rows = text.count("\n") - 1   # minus the header line
        print(f"[output] Background time-evolution data ({n_rows} rows) "
              f"written to {out_path}", file=sys.stderr)


class StandardBackground(Background):
    """The standard primat cosmological background.

    Builds the ``a <-> t <-> T`` relations and n<->p weak rates exactly as
    the pre-split ``PRIMAT`` did, under either of the two decoupling regimes
    selected by ``cfg.incomplete_decoupling`` (NEVO non-instantaneous
    decoupling, or instantaneous decoupling via EM entropy conservation), with
    the neutrino sector supplied by
    :func:`primat.neutrino_history.make_neutrino_history`.

    Parameters
    ----------
    cfg : PRIMATConfig
    plasma : primat.plasma.Plasma
    extra_rho : list of callable, optional
        See :class:`Background`.  Early Dark Energy (``cfg.fEDE > 0``) is
        appended automatically as the first such plug-in (see
        :meth:`_setup_EDE`); callers do not need to include it.
    """

    _a_of_t_scalar = None       # set by _store_background_arrays; see rhoB_BBN

    def __init__(self, cfg, plasma, extra_rho=None):
        super().__init__(cfg, plasma, extra_rho)
        self._setup_LCDM()
        self._setup_EDE()
        if cfg.verbose:
            print("[bg-py] Solving cosmological background: a(T), t(T) relations ...")
            _t_bg0 = time.time()
        self._setup_background_and_cosmo()
        if cfg.verbose:
            print(f"[bg-py] Background a(T), t(T) ready in {time.time()-_t_bg0:.2f} s")
        self._setup_weak_rates()

    # ======================================================================
    # ΛCDM energy-density components (CDM + cosmological constant)
    # ======================================================================

    def _setup_LCDM(self):
        r"""Append the ΛCDM cold dark matter and cosmological constant
        energy-density contributions to ``self.extra_rho``.

        During standard BBN (T ~ 1–0.001 MeV, t ~ 1–10^6 s) the CDM and Λ
        densities are completely negligible compared to the radiation density
        (by many orders of magnitude), so these terms have no effect on the
        standard run.  They become relevant only when ``T_end_MeV`` is reduced
        far below the standard 0.001 MeV, extending the integration towards
        matter-radiation equality or beyond.

        The cosmological constant is fixed by the flatness condition:

            Ω_Λ h² = h² − Ω_b h² − Ω_c h²

        so the total present-day energy density sums to the critical density,
        i.e. Ω_tot = 1 (flat universe).  If Ω_Λ < 0 (extreme non-standard
        parameters) a warning is printed but no exception is raised.

        The two extra densities as functions of the photon temperature Tg [MeV]
        are:

            ρ_CDM(Tg) = Ω_c h² · ρ_{crit,100} / a(Tg)³

        where ``a(Tg)`` is not yet known at call time.  The *comoving*
        amplitude ``rhocdm_a3 = Ω_c h² · ρ_{crit,100}`` is therefore stored and
        wrapped in a closure calling ``self.a_of_T(Tg)`` at evaluation time.
        ``_setup_background_and_cosmo`` publishes ``a_of_T`` at the end of its
        step 2 (the entropy ODE, itself ``a``-independent), strictly before
        step 3's time integration — the first thing to evaluate ``extra_rho``
        at all — so the closure always sees the exact ``a(T)``. Evaluating it
        earlier raises ``NotImplementedError`` from the base stub rather than
        silently substituting an approximation.

            ρ_Λ = (h² − Ω_b h² − Ω_c h²) · ρ_{crit,100}

        is a true constant (independent of a or T), appended as a constant
        callable.

        Both terms use ``cfg.rhocOverh2 = 3/(8π G_N) · H_{100}²``, the critical
        density at H = 100 km/s/Mpc in [MeV^4].  This follows the convention in
        ``PRIMATConfig.rhocOverh2`` and ``PRIMATConfig._update_derived``.

        References
        ----------
        See ``primat.config.PRIMATConfig.rhocOverh2`` for the unit convention.
        Planck 2018 + BAO fiducial values (Aghanim et al. 2020, A&A 641, A6,
        TT,TE,EE+lowE+lensing+BAO column): Ω_b h² = 0.02242, Ω_c h² = 0.11933,
        h = 0.6766.

        """
        cfg = self.cfg

        # Omegach2/h always exist and are always floats: PRIMATConfig type-
        # validates both (an explicit None raises "expected float"), so there
        # is no "None disables this contribution" case to handle. The C port
        # notes the same thing at background.c's setup_lcdm.
        Omegach2 = cfg.Omegach2
        h        = cfg.h

        # ρ_{crit,100} [MeV^4]: critical density at H = 100 km/s/Mpc.
        # cfg.rhocOverh2 = 3/(8π G_N) × H_{100}², so rhocrit100 = cfg.rhocOverh2.
        rhocrit100 = cfg.rhocOverh2

        # Comoving CDM amplitude: ρ_CDM(a) = rhocdm_a3 / a^3.
        # rhocdm_a3 = Ω_c h² × ρ_{crit,100}  [MeV^4]
        # Note: ρ_{crit,today} = h² × ρ_{crit,100}, so Ω_c × ρ_{crit,today}
        #       = Ω_c h² × ρ_{crit,100}.
        rhocdm_a3 = Omegach2 * rhocrit100   # [MeV^4]

        def rho_CDM(Tg):
            """CDM energy density [MeV^4] from the exact a(T) [MeV^4].

            Evaluates ρ_CDM = rhocdm_a3 / a(Tγ)³ using ``self.a_of_T``, which
            ``_setup_background_and_cosmo`` publishes as soon as step 2 has
            solved the entropy ODE -- before anything calls :meth:`Hubble`.
            Calling this earlier raises ``NotImplementedError`` from the base
            ``a_of_T`` stub, which is the intended loud failure.

            This used to be preceded by a ``rho_CDM_approx`` bootstrap that
            anchored ``a ≈ T0CMB_MeV / Tγ`` at today (a = 1, T = T0CMB) and was
            swapped out only *after* the background solve, so both ODEs ran on
            the approximation. That anchor is wrong in the BBN era by the e⁺e⁻
            reheating factor: a(T)·T is not constant across annihilation, it
            rises by (11/4)^{1/3}, so ``(a_approx/a_exact)³`` was 2.73 at
            Tγ = 10 MeV and 2.69 at 1 MeV, falling to 1.000 only below
            ~0.01 MeV. (The old comment justified it via matter-radiation
            equality, T ≫ T_eq -- the wrong mechanism, and one that made the
            approximation look best exactly where it was worst.) With
            ρ_CDM/ρ_γ ≈ 1.1e-6 at 1 MeV the 2.7× error biased ρ_tot by ~4e-7
            and H by ~2e-7 -- below the ±3e-9 D/H pin, but a needless
            divergence from primat-c, whose cpr_bg_Hubble always had the exact
            ``a`` available (it carries x = ln(a·T) in the ODE state).
            """
            a = self.a_of_T(Tg)
            return rhocdm_a3 / a**3

        # Cosmological constant: flatness condition Ω_tot = 1.
        # ρ_Λ = (h² − Ω_b h² − Ω_c h²) × ρ_{crit,100}
        Omegalambdah2 = h**2 - cfg.Omegabh2 - Omegach2   # = Ω_Λ h²
        rholambda = Omegalambdah2 * rhocrit100             # [MeV^4]

        if Omegalambdah2 < 0:
            import warnings
            warnings.warn(
                f"_setup_LCDM: Omega_Lambda h^2 = {Omegalambdah2:.4g} < 0 "
                f"(h={h}, Omegabh2={cfg.Omegabh2}, Omegach2={Omegach2}).  "
                "Cosmological constant is negative -- non-standard cosmology.",
                stacklevel=3,
            )

        def rho_Lambda(_Tg):
            """Cosmological constant energy density [MeV^4] (T-independent)."""
            return rholambda

        # Both are appended in their final form -- no bootstrap/swap step any
        # more (see rho_CDM's docstring). Nothing between here and the point
        # where _setup_background_and_cosmo publishes self.a_of_T evaluates
        # extra_rho: _setup_EDE only reads rho_g, step 1 builds the neutrino
        # history from tables/spl, and step 2's entropy ODE is a-independent
        # (it is pure spl/N_NEVO, never Hubble). The first extra_rho consumer
        # is the t(a) integration in step 3-4, by which time a_of_T is set.
        self.extra_rho.append(rho_CDM)
        self.extra_rho.append(rho_Lambda)

    # ======================================================================
    # Early Dark Energy setup
    # ======================================================================

    def _setup_EDE(self):
        """Build the EDE energy-density function from cfg.fEDE/zcEDE/wnEDE.

        If fEDE > 0, appends a ``rho_EDE(Tg) -> MeV^4`` callable to
        ``self.extra_rho`` (the generic extra-energy-density plug-in list);
        otherwise a no-op.  Must be called after ``self.plasma`` and
        ``self.extra_rho`` are set, since it evaluates rho_g and appends to
        that list.
        """
        cfg = self.cfg
        if cfg.fEDE == 0.:
            return

        thermo  = self.plasma
        acEDE   = 1. / (1. + cfg.zcEDE)
        amaxEDE = acEDE * (4. / (3. * cfg.wnEDE - 1.))**(1. / (3. * cfg.wnEDE + 3.))
        TmaxEDE = cfg.T0CMB / amaxEDE / cfg.MeV_to_Kelvin   # [MeV]
        TcEDE   = cfg.T0CMB / acEDE   / cfg.MeV_to_Kelvin   # [MeV]

        # The final Neff value in the standard case (cfg.Neff_SM) enters here.
        rhocEDEac = (cfg.fEDE / (1. - cfg.fEDE)
                     * thermo.rho_g(TmaxEDE)
                     * (1. + cfg.Neff_SM * 7./8. * (4./11.)**(4./3.))
                     / 2.
                     * (1. + 4. / (3. * cfg.wnEDE - 1.)))

        def rho_EDE(T):
            return 2. * rhocEDEac / (1. + (TcEDE / T)**(3. * cfg.wnEDE + 3.))

        self.extra_rho.append(rho_EDE)

    # ======================================================================
    # Friedmann expansion rate
    # ======================================================================

    def Hubble(self, Tg, Tnue, Tnumu, Tnutau):
        """Friedmann expansion rate H [s^-1] at photon temperature ``Tg``
        [MeV] and the three flavour neutrino temperatures
        ``Tnue``/``Tnumu``/``Tnutau`` [MeV] (array-safe).

        Internal to :class:`StandardBackground` (drives the a(T) and t(a)
        ODEs in :meth:`_setup_background_and_cosmo`, and the ``H`` column of
        :meth:`_background_columns`) -- not part of the minimal
        :class:`Background` interface.
        """
        cfg     = self.cfg
        thermo  = self.plasma
        rho_pl  = thermo.rho_g(Tg) + thermo.rho_e(Tg) - thermo.PQEDofT(Tg) + Tg * thermo.dPQEDdT(Tg)
        # Per-flavour neutrino energy density: each of ν_e, ν_μ, ν_τ
        # carries its own reduced chemical potential ξ (cfg.xi_nu_e/mu/tau,
        # each defaulting to cfg.munuOverTnu). All three gravitate; only ξ_e
        # also shifts the weak rates (handled in weak_rates, not here).
        rho_3nu = (thermo.rho_nu(Tnue,   cfg.xi_nu_e)
                   + thermo.rho_nu(Tnumu,  cfg.xi_nu_mu)
                   + thermo.rho_nu(Tnutau, cfg.xi_nu_tau))
        rho_tot = rho_pl + rho_3nu + thermo.rho_nu_extra(Tg)
        for rho_extra in self.extra_rho:
            rho_tot += rho_extra(Tg)
        # For analytic spectral distortions the neutrino phase-space distribution
        # is shifted from a perfect FD, adding extra energy density.  The NEVO
        # case needs no correction: the NEVO temperatures are defined as the
        # energy-equivalent FD temperature, so rho_nu already accounts for the
        # distortion.
        if self.rho_nu_SD is not None:
            # Average T_ν: use the energy-weighted mean of the three flavours.
            Tnu_avg = ((Tnue**4 + Tnumu**4 + Tnutau**4) / 3.)**0.25
            rho_tot += self.rho_nu_SD(Tnu_avg)
        return cfg.MeV_to_secm1 * (rho_tot * 8. * np.pi / (3. * cfg.Mpl**2))**0.5

    # ======================================================================
    # Background thermodynamics + cosmological setup
    # ======================================================================

    def _setup_background_and_cosmo(self):
        """Cosmological background thermodynamics.

        Two operating modes, selected by ``cfg.incomplete_decoupling``:

        *Incomplete decoupling* (``True``, default):
            Reads the pre-computed NEVO neutrino-decoupling table
            (``rates/NEVO/NEVOPRIMAT_col_1_7.csv``).  The three neutrino
            flavour temperatures T_νe, T_νμ, T_ντ are interpolated from the
            table, and the NEVO heating function N(T_γ) — representing the
            extra entropy injected into neutrinos during e+e- annihilations —
            drives the a(T_γ) ODE.

        *Instantaneous (complete) decoupling* (``False``):
            The NEVO table is **not** loaded.  All three neutrino flavours are
            assumed to have decoupled instantaneously from the EM plasma, so
            their common temperature is given by EM entropy conservation:

                T_ν = (spl(T_γ) / s_∞)^{1/3}

            where s_∞ = 11π²/45 is the high-T limit of spl(T)/T³ (photons +
            e+e- pairs at T >> m_e).  The NEVO heating function is set to
            N ≡ 0, reducing the a(T_γ) ODE to standard EM entropy
            conservation.  QED corrections to the plasma equation of state are
            still included via the spl/PQEDofT tables in both modes.

        In both modes the method builds all interpolants used by
        ``_setup_weak_rates``, the Omega_nu accessors, and the nuclear
        network's ``solve()``.

        Independently of the above, ``cfg.external_scale_factor`` selects how
        ``a(T_γ)`` itself is obtained (requires ``incomplete_decoupling=True``):

        *Minimal* (``False``, default):
            ``a(T_γ)`` is reconstructed by solving the entropy-conservation
            ODE ``d(ln a)/d(ln T) = -(3 s̄ + T ds̄/dT)/(N_NEVO + 3 s̄)`` driven
            by the NEVO heating function ``N_NEVO(T_γ)``.

        *External* (``True``):
            ``a(T_γ)`` is read directly from the NEVO table's ``x`` column
            (``x ∝ a`` by the NEVO convention), with radiation-domination
            extrapolation (``a ∝ 1/T_γ``) outside the table. No ODE is
            solved for ``a(T)``. ``t(a)`` is obtained the same way in both
            modes (Hubble integration below), since no NEVO file carries a
            cosmic-time column. The two modes agree to ~1e-5 on the
            observables: measured on the default config, ``Neff`` to 1.5e-16,
            ``YPBBN`` to 7.0e-07, ``D/H`` to 8.8e-06 and ``Li7/H`` to 2.0e-05
            relative.

        The five numbered steps below are delegated to
        :meth:`_setup_neutrino_sector` (1), :meth:`_build_a_of_T` (2), and
        :meth:`_invert_and_integrate_time` / :meth:`_store_background_arrays`
        (3-4 / 5); this method is the orchestrator threading their outputs
        together.
        """
        cfg    = self.cfg

        Tstartcosmo  = cfg.T_start_cosmo / cfg.MeV_to_Kelvin
        Tstart = cfg.T_start / cfg.MeV_to_Kelvin   # [MeV]
        Tend   = cfg.T_end   / cfg.MeV_to_Kelvin   # [MeV]

        # Step 1 - Neutrino-sector background (temperatures, heating,
        #          spectral distortion, extra neutrino energy density)
        nh, Tnue_of_Tg, Tnumu_of_Tg, Tnutau_of_Tg, N_NEVO_of_Tg = \
            self._setup_neutrino_sector()

        # Step 2 - Build a(T)
        a_of_T, T_sol, n_T_pts, dlnadlnT = self._build_a_of_T(
            nh, N_NEVO_of_Tg, Tend, Tstartcosmo)

        # Publish a(T) as soon as it exists. _setup_LCDM's rho_CDM closure
        # reads self.a_of_T, and steps 3-4 below are the first thing to
        # evaluate it (through Hubble); doing this here rather than in step 5
        # is what lets the CDM term be exact inside the time integration
        # instead of running on a radiation-domination approximation that is
        # 2.7x off in the BBN era -- see rho_CDM's docstring. Step 5 re-assigns
        # the same object, harmlessly.
        self.a_of_T = a_of_T

        # Steps 3-4 - Cosmic time.
        #
        # Minimal mode integrates dt/d(ln T) = dlnadlnT(T)/H(T) straight over
        # ln T (no a(T) inversion anywhere -- see _integrate_time_over_lnT for
        # the accuracy argument). Only external_scale_factor, which has no
        # closed-form d(ln a)/d(ln T), still needs the
        # invert-then-integrate-over-ln-a route.
        if dlnadlnT is not None:
            t_vec, Tg_vec, a_arr = self._integrate_time_over_lnT(
                a_of_T, dlnadlnT, T_sol, Tnue_of_Tg, Tnumu_of_Tg,
                Tnutau_of_Tg, Tstartcosmo)
        else:
            t_vec, Tg_vec, a_arr = self._invert_and_integrate_time(
                a_of_T, T_sol, n_T_pts, Tnue_of_Tg, Tnumu_of_Tg, Tnutau_of_Tg,
                Tstartcosmo, Tend)

        # Step 5 - Sample on the common time grid; set instance attributes
        self._store_background_arrays(
            t_vec, Tg_vec, a_arr, a_of_T, Tnue_of_Tg, Tnumu_of_Tg,
            Tnutau_of_Tg, N_NEVO_of_Tg)

    def _setup_neutrino_sector(self):
        """Step 1 of :meth:`_setup_background_and_cosmo`: the neutrino-sector
        background (temperatures, heating, spectral distortion, extra
        neutrino energy density).

        The neutrino sector is encapsulated in a NeutrinoHistory object
        (primat.neutrino_history): NEVOTable for incomplete
        decoupling, InstantaneousDecoupling otherwise, optionally decorated
        with the analytic μ+y spectral distortion (AnalyticDistortion).  It
        exposes the three flavour temperature functions, the NEVO heating
        N(T_γ) that drives the a(T_γ) ODE, the n<->p weak-rate distortion
        dFDneu_func, and the extra neutrino energy density rho_nu_SD.

        Sets ``self.dFDneu_func``/``self.dFDneu_moments``/``self.rho_nu_SD``
        (spectral-distortion hooks consumed by ``Hubble`` and
        ``_setup_weak_rates``; all ``None`` when there are no distortions).

        Returns
        -------
        (nh, Tnue_of_Tg, Tnumu_of_Tg, Tnutau_of_Tg, N_NEVO_of_Tg)
        """
        cfg    = self.cfg
        thermo = self.plasma
        nh = make_neutrino_history(cfg, thermo)

        Tnue_of_Tg   = nh.Tnue_of_Tg
        Tnumu_of_Tg  = nh.Tnumu_of_Tg
        Tnutau_of_Tg = nh.Tnutau_of_Tg
        N_NEVO_of_Tg = nh.N_NEVO_of_Tg

        self.dFDneu_func = nh.dFDneu_func   # None means "no spectral distortions"
        self.dFDneu_moments = nh.dFDneu_moments  # None unless analytic-distortion mode
        self.rho_nu_SD   = nh.rho_nu_SD     # None means "no extra energy density"

        return nh, Tnue_of_Tg, Tnumu_of_Tg, Tnutau_of_Tg, N_NEVO_of_Tg

    def _build_a_of_T(self, nh, N_NEVO_of_Tg, Tend, Tstartcosmo):
        """Step 2 of :meth:`_setup_background_and_cosmo`: build a(T_γ) (see
        that method's docstring for the *minimal* vs *external* mode
        details).

        Returns
        -------
        (a_of_T, T_sol, n_T_pts, dlnadlnT)
            ``a_of_T`` : the scale-factor callable (array-safe).
            ``T_sol``  : the temperature grid it was solved/sampled on.
            ``n_T_pts``: that grid's point count (reused by
            :meth:`_invert_and_integrate_time` for the log(a) sampling grid).
            ``dlnadlnT``: scalar callable ``T -> d(ln a)/d(ln T_γ)``, the
            entropy-conservation right-hand side, in *minimal* mode; ``None``
            in ``external_scale_factor`` mode (where a(T) is an algebraic
            table read, not an ODE solution, so no such closed form exists).
            :meth:`_integrate_time_over_lnT` reuses it as the chain-rule
            factor that lets t be integrated directly over ln T -- see that
            method for why this matters.
        """
        cfg    = self.cfg
        thermo = self.plasma

        def _sbar(T):
            return thermo.spl(T) / T**3   # dimensionless

        z0   = cfg.T0CMB / cfg.MeV_to_Kelvin   # [MeV]
        # Algebraic entropy-conservation boundary value a(Tend) = zend/Tend,
        # used as the ODE initial condition (minimal mode) and as the
        # table-normalisation anchor (external_scale_factor mode).
        # Requires no ODE: _sbar is the analytic
        # electron-thermo spline.
        zend = z0 / (_sbar(Tend) / cfg.s0bar) ** (1. / 3.)
        a_end = zend / Tend

        # Build the log-temperature grid directly with linspace and feed it
        # straight to t_eval.  Do NOT reconstruct it as log(logspace(...)):
        # the log->exp->log roundtrip can push an endpoint 1 ULP outside the
        # integration span [log(Tend), log(Tstartcosmo)], which makes
        # solve_ivp raise "Values in t_eval are not within t_span" for some
        # values of Tend (e.g. T_end = 2e-3 MeV fails while 1e-3 MeV happens
        # to roundtrip exactly).  np.linspace guarantees its endpoints equal
        # the span bounds exactly, so the check always passes.
        n_T_pts = primat_weak_rates.n_points_per_decade(cfg.sampling_temperature_per_decade, Tend, Tstartcosmo)
        lnT_sol = np.linspace(np.log(Tend), np.log(Tstartcosmo), n_T_pts)
        T_sol   = np.exp(lnT_sol)

        if cfg.external_scale_factor:
            # --------------------------------------------------------------
            # external_scale_factor=True: a(T) read directly from the NEVO
            # table's x column (a ∝ x by the NEVO convention),
            # normalised so a(Tend) matches the algebraic
            # a_end above -- the same boundary value the minimal-mode ODE
            # converges to.  No ODE solve.
            # --------------------------------------------------------------
            K = a_end / float(nh.x_of_Tg(Tend))

            def a_of_T(T):
                return K * nh.x_of_Tg(T)

            # No entropy-conservation ODE in this mode, so no closed-form
            # d(ln a)/d(ln T): the time integration falls back to the
            # inverse-interpolation path (see _setup_background_and_cosmo).
            dlnadlnT = None
        else:
            # --------------------------------------------------------------
            # minimal mode (default): solve the entropy-conservation ODE
            # d(ln a)/d(ln T) for the EM plasma.  The reduced entropy
            # sbar = s/T^3 and its T-derivative are both obtained
            # analytically from the electron-thermo spline via
            # thermo.spl_and_dspl_dT (a single evaluation returning s and
            # ds/dT): this is exact and fast, so no numerical-derivative
            # fallback (numdifftools / finite differences) is needed.
            # --------------------------------------------------------------
            def dlnadlnT(T):
                """d(ln a)/d(ln T_γ) from EM entropy conservation with NEVO
                heating: -(3 s̄ + T ds̄/dT) / (N_NEVO + 3 s̄), with s̄ = s/T³.

                Scalar (not solve_ivp-shaped) so it can be reused verbatim as
                the chain-rule factor in :meth:`_integrate_time_over_lnT`;
                ``_dlnadlnT_NEVO`` below is the list-returning wrapper the
                a(T) ODE needs.
                """
                s, ds_dT = thermo.spl_and_dspl_dT(T)
                sb     = s / T**3
                dsbdT  = ds_dT / T**3 - 3. * s / T**4   # d(s/T^3)/dT, chain rule
                N = float(N_NEVO_of_Tg(T))
                return -(3. * sb + T * dsbdT) / (N + 3. * sb)

            def _dlnadlnT_NEVO(lnT, y):
                return [dlnadlnT(np.exp(lnT))]

            lna_end = np.log(a_end)
            _t_nevo_a0 = time.time()
            sol_lna = solve_ivp(_dlnadlnT_NEVO,
                                [np.log(Tend), np.log(Tstartcosmo)],
                                [lna_end],
                                t_eval=lnT_sol,
                                method='LSODA', rtol=0.1*cfg.numerical_precision, atol=1e-10)
            if cfg.debug:
                print((f"[bckg]  Finished a(T) solve in {time.time()- _t_nevo_a0:.2f} s "
                       f"(status={sol_lna.status}, nfev={sol_lna.nfev})"), flush=True)
            _lnalnT = interp1d(sol_lna.t, sol_lna.y[0].flatten(),
                               bounds_error=False, fill_value="extrapolate")

            def a_of_T(T):
                return np.exp(_lnalnT(np.log(T)))

        return a_of_T, T_sol, n_T_pts, dlnadlnT

    def _integrate_time_over_lnT(self, a_of_T, dlnadlnT, T_sol, Tnue_of_Tg,
                                 Tnumu_of_Tg, Tnutau_of_Tg, Tstartcosmo):
        """Steps 3-4, *minimal* mode: integrate t directly over ln T_γ.

        Why not over ln a (the historical formulation, kept in
        :meth:`_invert_and_integrate_time` for ``external_scale_factor``):
        integrating ``dt/d(ln a) = 1/H(T(a))`` needs ``T(a)``, i.e. the
        *inverse* of a(T), which is only available as an interpolant over the
        a(T) solve's nodes.  That interpolant is a **linear** ``interp1d`` in
        (a, T) and carries a median relative error of 3.9e-06 at the default
        ``sampling_temperature_per_decade = 600`` -- and because it sits
        inside the ODE's right-hand side, that error caps t's accuracy no
        matter how tight ``cfg.numerical_precision`` is.  The stored ``Tg_vec``
        was then built by the *same* inverse interpolant, injecting it twice.
        Measured consequence: t(T) carried ~1e-05 relative error, which
        propagated to ~1e-06 in YP and D/H -- larger than the ODE tolerance
        implies, and invisible under the cross-backend test tolerance.

        Changing the independent variable to ln T removes the inversion
        entirely.  By the chain rule,

            dt/d(ln T) = dt/d(ln a) * d(ln a)/d(ln T) = dlnadlnT(T) / H(T)

        whose right-hand side depends on T *only* -- and T is the integration
        variable, so it is exact.  ``a`` enters only through H's CDM term, via
        the log-log-accurate ``a_of_T`` (published by
        :meth:`_setup_background_and_cosmo` before this runs), never through
        an inverse lookup.  The output grid is ``T_sol`` itself, so ``Tg_vec``
        is exact too.

        This mirrors ``combined_bg_rhs`` in
        ``primat-c/src/background.c``, which integrates the same chain-rule
        right-hand side over ``lnT``; the C co-integrates ``x = ln(aT)`` in the
        same pass to avoid any a(T) lookup, which Python does not need because
        step 2 has already produced an accurate ``a_of_T``.

        Args:
            a_of_T: scale-factor callable (array-safe), from step 2.
            dlnadlnT: scalar callable ``T -> d(ln a)/d(ln T)``, from step 2.
            T_sol: the temperature grid [MeV], ascending, from step 2.
            Tnue_of_Tg, Tnumu_of_Tg, Tnutau_of_Tg: flavour temperature
                interpolants [MeV].
            Tstartcosmo: cosmological start temperature [MeV], where the
                radiation-domination boundary condition is applied.

        Returns:
            ``(t_vec, Tg_vec, a_arr)``: cosmic time [s] ascending, the photon
            temperatures [MeV] it corresponds to (descending), and the scale
            factor at those temperatures -- the three arrays
            :meth:`_store_background_arrays` turns into the public
            interpolants.
        """
        cfg = self.cfg

        def Hubble_NEVO(Tg):
            return self.Hubble(Tg, Tnue_of_Tg(Tg), Tnumu_of_Tg(Tg), Tnutau_of_Tg(Tg))

        # Radiation-domination boundary condition t = 1/(2H), applied at
        # T_start_cosmo where the universe is deep in the relativistic era.
        # Same formula and same anchor point as the historical t(a) solve, so
        # this reformulation changes the discretisation, not the physics.
        t_ini = 1. / (2. * Hubble_NEVO(Tstartcosmo))

        def _dtdlnT(lnT, y):
            T = np.exp(lnT)
            return [dlnadlnT(T) / Hubble_NEVO(T)]

        # Integrate downward in temperature, from the boundary condition at
        # T_start_cosmo to T_end, reporting on step 2's own grid reversed to
        # descending order (t_eval must run in the integration direction).
        lnT_desc = np.log(T_sol)[::-1]
        _t_nevo_t0 = time.time()
        sol_t = solve_ivp(_dtdlnT,
                          [lnT_desc[0], lnT_desc[-1]],
                          [t_ini],
                          t_eval=lnT_desc,
                          method='LSODA', rtol=cfg.numerical_precision, atol=1e-12)
        if cfg.debug:
            print((f"[bckg]  Finished t(T) solve in {time.time()-_t_nevo_t0:.2f} s "
                   f"(status={sol_t.status}, nfev={sol_t.nfev})"), flush=True)

        Tg_vec = np.exp(sol_t.t)          # descending T, exactly step 2's grid
        t_vec  = sol_t.y[0].flatten()     # ascending t [s]
        a_arr  = a_of_T(Tg_vec)
        return t_vec, Tg_vec, a_arr

    def _invert_and_integrate_time(self, a_of_T, T_sol, n_T_pts, Tnue_of_Tg,
                                   Tnumu_of_Tg, Tnutau_of_Tg, Tstartcosmo, Tend):
        """Steps 3-4 of :meth:`_setup_background_and_cosmo`: invert a(T) ->
        T(a), then integrate dt/d(ln a) = 1/H(a).

        Used **only** for ``cfg.external_scale_factor=True``, where a(T) is an
        algebraic read of the NEVO ``x`` column rather than an ODE solution,
        so there is no closed-form ``d(ln a)/d(ln T)`` for the cheaper,
        inversion-free route :meth:`_integrate_time_over_lnT` takes.  See that
        method for what the ``T_of_a`` inversion below costs in accuracy
        (~1e-05 in t(T)); the same caveat applies here, and the C backend's
        ``dtdlna_rhs`` carries it in the same mode for the same reason.

        Returns ``(t_vec, Tg_vec, a_arr)`` -- the same triple
        :meth:`_integrate_time_over_lnT` returns, so
        :meth:`_store_background_arrays` is agnostic to which path produced it.
        """
        cfg = self.cfg

        # Refined against T_sol rather than reusing it. In this mode a(T) is a
        # closed-form table read (the NEVO `x` column), not an ODE solution, so
        # it can be sampled as finely as wanted for the price of the read --
        # while T_sol's density is set by sampling_temperature_per_decade,
        # which sizes the *ODE output*. Reusing it left this inverse carrying
        # an O(h^2) error of ~5e-06 in D/H at the default density, and the C
        # backend's cubic over the same nodes carried ~7e-06 of the opposite
        # sign, so the two backends sat 1.2e-05 apart at converged tolerance --
        # 12x the budget test_backend_parity's converged test allows. Both are
        # second order, so the factor below buys ~64x on each. Mirrored by
        # background.c, which refines the same way for the same reason.
        T_grid = np.logspace(np.log10(T_sol[0]), np.log10(T_sol[-1]),
                              (len(T_sol) - 1) * _EXTERNAL_A_REFINE + 1)
        a_grid = a_of_T(T_grid)                  # low a → high a (a_of_T is array-safe)

        # Linear, deliberately, and deliberately unlike the C backend, which
        # fits a not-a-knot cubic over these same nodes. That cubic is a
        # performance workaround for C's RK45 stepper, which otherwise rejects
        # steps on the curvature kinks at every node; LSODA here does not, so
        # there is nothing to buy. A cubic would in fact cost accuracy: in
        # external mode a(T) is a table read (the NEVO `x` column), so fitting
        # a cubic through these nodes manufactures curvature the data does not
        # contain. Each backend uses the interpolant suited to its integrator.
        T_of_a = interp1d(a_grid, T_grid, bounds_error=False, fill_value="extrapolate")

        a_ini = a_of_T(Tstartcosmo)
        a_fin = a_of_T(Tend)

        def Hubble_NEVO(Tg):
            return self.Hubble(Tg, Tnue_of_Tg(Tg), Tnumu_of_Tg(Tg), Tnutau_of_Tg(Tg))

        t_ini = 1. / (2. * Hubble_NEVO(Tstartcosmo))

        # Log-scale-factor grid, built directly in log space (see the note on
        # the a(T) solve above): feeding linspace endpoints straight to t_eval
        # avoids the log(logspace(...)) roundtrip that could land the last
        # point 1 ULP outside [log(a_ini), log(a_fin)].
        # Refined by the same factor as the inverse above, and for the same
        # reason: t_of_T/T_of_t are built directly on this output grid
        # (_store_background_arrays), so its density -- not the ODE tolerance
        # -- sets how accurately the nuclear network reads t(T) in this mode.
        # At the unrefined density the two backends' O(h^2) errors, opposite in
        # sign, left them 1.2e-05 apart in D/H at converged tolerance.
        lna_samp = np.linspace(np.log(a_ini), np.log(a_fin),
                                (n_T_pts - 1) * _EXTERNAL_A_REFINE + 1)

        def _dtdlna(lna, t):
            return [1. / Hubble_NEVO(T_of_a(np.exp(lna)))]

        _t_nevo_t0 = time.time()
        sol_t = solve_ivp(_dtdlna,
                          [np.log(a_ini), np.log(a_fin)],
                          [t_ini],
                          t_eval=lna_samp,
                          method='LSODA', rtol=cfg.numerical_precision, atol=1e-12)
        if cfg.debug:
            print((f"[bckg]  Finished t(a) solve in {time.time()-_t_nevo_t0:.2f} s "
                   f"(status={sol_t.status}, nfev={sol_t.nfev})"), flush=True)

        a_arr  = np.exp(sol_t.t)       # a values at ODE evaluation points
        t_vec  = sol_t.y[0].flatten()  # corresponding t [s]
        Tg_vec = T_of_a(a_arr)         # T_γ [MeV] (via the linear inverse)
        return t_vec, Tg_vec, a_arr

    def _store_background_arrays(self, t_vec, Tg_vec, a_arr, a_of_T,
                                  Tnue_of_Tg, Tnumu_of_Tg, Tnutau_of_Tg,
                                  N_NEVO_of_Tg):
        """Step 5 of :meth:`_setup_background_and_cosmo`: sample onto the
        common time grid and set the public interpolant/array instance
        attributes (``t_vec``, ``Tg_vec``, ``t_of_T``, ``T_of_t``, ``a_of_T``,
        ``a_of_t``, ...) consumed by ``_setup_weak_rates``, the Omega_nu
        accessors, and the nuclear network's ``solve()``.

        Takes the ``(t_vec, Tg_vec, a_arr)`` triple produced by whichever of
        :meth:`_integrate_time_over_lnT` (minimal mode) or
        :meth:`_invert_and_integrate_time` (``external_scale_factor``) ran,
        so this step is identical either way.
        """
        cfg = self.cfg

        Tnue_vec   = Tnue_of_Tg(Tg_vec)
        Tnumu_vec  = Tnumu_of_Tg(Tg_vec)
        Tnutau_vec = Tnutau_of_Tg(Tg_vec)
        # Energy-weighted average neutrino temperature (for weak rates / Omega_ν)
        Tnu_avg_vec = ((Tnue_vec**4 + Tnumu_vec**4 + Tnutau_vec**4) / 3.)**0.25

        self.t_vec      = t_vec
        self.Tg_vec     = Tg_vec
        self.Tnu_vec    = Tnu_avg_vec   # average, used by _Tnu_today and _setup_weak_rates
        self.Tnue_vec   = Tnue_vec
        self.Tnumu_vec  = Tnumu_vec
        self.Tnutau_vec = Tnutau_vec

        # Log-log rather than linear: both relations are near power laws over
        # ~4.6 decades, and T_of_t sits in every nuclear-network RHS call.
        # See _loglog_interp1d.
        self.t_of_T = _loglog_interp1d(Tg_vec, t_vec)
        self.T_of_t = _loglog_interp1d(t_vec, Tg_vec)
        self.TnuofT = interp1d(Tg_vec, Tnu_avg_vec, bounds_error=False,
                                fill_value="extrapolate", kind='linear')
        self.a_of_T = a_of_T   # already vectorised: np.exp(interp1d(log T))
        # Public T(a), rebuilt here from the *final* (a, T) pairs rather than
        # taken from the time integration: in minimal mode those pairs are
        # exact (Tg_vec is step 2's own grid, a_arr = a_of_T(Tg_vec) with the
        # log-log a_of_T). Linear, matching the C backend's cpr_bg_T_of_a.
        #
        # Past the last node the linear extrapolation reaches T = 0 at about
        # 2 a_end and goes negative beyond, so the tail is continued as the
        # power law through the last two nodes instead (T ~ 1/a in radiation
        # domination): positive for every a, and the same scheme a_of_T
        # already uses in the other direction.
        self.T_of_a = _linear_with_powerlaw_tail(a_arr, Tg_vec)
        self.a_of_t = interp1d(t_vec, a_arr, bounds_error=False,
                                fill_value=(a_arr[0], a_arr[-1]))
        self._a_of_t_scalar = _scalar_linear_eval(self.a_of_t)   # rhoB_BBN hot path
        self.t_of_a = interp1d(a_arr, t_vec, bounds_error=False,
                                fill_value=(t_vec[0], t_vec[-1]))
        self.has_scale_factor = True
        self.N_NEVO_of_Tg = N_NEVO_of_Tg

        # Whether N_NEVO_of_Tg is a *real* heating table (NEVOTable, read from
        # rates/NEVO/) or just the N=0 stub used by InstantaneousDecoupling to
        # close the a(T) ODE under EM entropy conservation.  Consumed by
        # _background_columns to decide whether the "Nheating" column carries
        # physical information or would just be a column of zeros.
        self.has_heating_table = cfg.incomplete_decoupling

    def _Tnu_today(self):
        """Relic neutrino temperature today [MeV], redshifted from the stored
        final Tnu/Tg ratio. Requires self.Tg_vec, self.Tnu_vec."""
        cfg = self.cfg
        return self.Tnu_vec[-1] / self.Tg_vec[-1] * cfg.T0CMB / cfg.MeV_to_Kelvin

    def rho_nu_total_final(self):
        """Final-time ``(Tg, rho_nu_tot)`` (see :meth:`Background.rho_nu_total_final`).

        ``rho_nu_tot`` sums the three flavour-dependent neutrino energy
        densities (at the final ``Tnue``/``Tnumu``/``Tnutau``), the extra
        neutrino energy density ``rho_nu_extra(Tg)`` (e.g. from
        ``cfg.DeltaNeff``), and -- if analytic spectral distortions are
        active (:attr:`rho_nu_SD` is not ``None``) -- the extra
        spectral-distortion energy density evaluated at the final
        energy-weighted average neutrino temperature.  This mirrors the
        ``rho_rad - rho_g`` numerator of the former ``N_eff`` closure, minus
        the ``rho_g`` term itself (now folded into the generic
        :meth:`Background.N_eff`).
        """
        thermo = self.plasma
        Tg_f     = self.Tg_vec[-1]
        Tnue_f   = self.Tnue_vec[-1]
        Tnumu_f  = self.Tnumu_vec[-1]
        Tnutau_f = self.Tnutau_vec[-1]

        cfg = self.cfg
        # Per-flavour ξ, matching Hubble(): each flavour contributes its
        # own chemical-potential excess to the final-time neutrino energy.
        rho_nu_tot_f = (thermo.rho_nu(Tnue_f,   cfg.xi_nu_e)
                        + thermo.rho_nu(Tnumu_f,  cfg.xi_nu_mu)
                        + thermo.rho_nu(Tnutau_f, cfg.xi_nu_tau)
                        + thermo.rho_nu_extra(Tg_f))

        if self.rho_nu_SD is not None:
            # Energy-weighted average T_ν at the final time (self.Tnu_vec is
            # exactly this average, see _setup_background_and_cosmo).
            Tnu_avg_f = self.Tnu_vec[-1]
            rho_nu_tot_f += self.rho_nu_SD(Tnu_avg_f)

        return Tg_f, rho_nu_tot_f

    # Both PER FLAVOUR (nu + nubar) -- there is deliberately no factor 3 here,
    # see the base-class docstrings. 7 pi^2/120 * T^4 is plasma.rho_nu(T), i.e.
    # one flavour; 3/2 zeta(3)/pi^2 * T^3 is one flavour's number density, the
    # normalisation behind the standard sum(m_nu)/93.1 eV relation.
    def Omeganuh2_relnu(self):
        """Omega_nu h^2 per flavour, relativistic-neutrino convention (see
        :meth:`Background.Omeganuh2_relnu` for the per-flavour caveat)."""
        return (7. * np.pi**2 / 120. * self._Tnu_today()**4) / self.cfg.rhocOverh2

    def Omeganuh2_nrnu(self):
        """Omega_nu h^2 per flavour per unit mass, non-relativistic-neutrino
        convention (see :meth:`Background.Omeganuh2_nrnu`)."""
        return (3. / 2. * zeta(3) / np.pi**2 * self._Tnu_today()**3) / self.cfg.rhocOverh2

    # ======================================================================
    # Background time-evolution output
    # ======================================================================

    def Tnu_of_t(self, t_out):
        """Per-flavour neutrino temperature [MeV] at cosmic time ``t_out``
        [s] (array-safe); see :meth:`Background.Tnu_of_t`.

        What is interpolated (linear, extrapolated, on ``self.t_vec``) is the
        dimensionless ratio T_nu/T_gamma, which is then multiplied by
        ``T_of_t``.  Interpolating the three temperatures directly would put
        them on a different scheme from ``T_of_t``'s log-log one, and between
        nodes the two schemes disagree enough to report T_nu > T_gamma above
        e+e- annihilation, where the two are physically equal.  Going through
        the ratio makes T_nu <= T_gamma hold by construction, since every
        tabulated ratio is <= 1 and linear interpolation cannot exceed its
        endpoints.  :class:`primat.neutrino_history.NEVOTable` interpolates its
        own T_nu(T_gamma) through the same ratio, for the same reason.
        """
        Tg_out = self.T_of_t(t_out)
        ratio_kw = dict(bounds_error=False, fill_value="extrapolate",
                        kind='linear')
        out = {}
        for key, Tnu_vec in (("e", self.Tnue_vec), ("mu", self.Tnumu_vec),
                             ("tau", self.Tnutau_vec)):
            ratio = interp1d(self.t_vec, Tnu_vec / self.Tg_vec, **ratio_kw)
            out[key] = ratio(t_out) * Tg_out
        return out

    def _background_columns(self, t_out):
        """Background output columns (see :meth:`Background._background_columns`).

        Extends the base ``T``, ``t`` columns with ``a``, ``H``, the three
        flavour neutrino temperatures, the NEVO heating function (if
        ``has_heating_table``), and the plasma/neutrino/extra/total energy
        densities [MeV^4].
        """
        cfg    = self.cfg
        thermo = self.plasma
        cols = super()._background_columns(t_out)
        T_out = cols["T [MeV]"]

        a_out = self.a_of_t(t_out)

        Tnu = self.Tnu_of_t(t_out)
        Tnue_out, Tnumu_out, Tnutau_out = Tnu["e"], Tnu["mu"], Tnu["tau"]

        H_out = np.array([
            self.Hubble(T_out[i], Tnue_out[i], Tnumu_out[i], Tnutau_out[i])
            for i in range(t_out.size)
        ])

        cols["a [1]"]          = a_out
        cols["H [s^-1]"]       = H_out
        cols["Tnue [MeV]"]     = Tnue_out
        cols["Tnumu [MeV]"]    = Tnumu_out
        cols["Tnutau [MeV]"]   = Tnutau_out

        # The "Nheating" column is the NEVO heating function N(T_gamma) that
        # drives the a(T_gamma) ODE (see _setup_background_and_cosmo).  It is
        # only physically meaningful when a real NEVO table was loaded
        # (cfg.incomplete_decoupling=True); under InstantaneousDecoupling,
        # N_NEVO_of_Tg is the N=0 stub used to close the ODE under plain EM
        # entropy conservation, so writing it out would just be a column of
        # zeros masquerading as data.
        if self.has_heating_table:
            cols["Nheating [1]"]  = self.N_NEVO_of_Tg(T_out)

        # Energy densities [MeV^4]: plasma (photons + e+- pairs, with QED
        # corrections), total neutrino sector (three flavours + extra +
        # spectral-distortion contribution), optional extra (Early Dark
        # Energy etc.), and their sum.
        #
        # thermo.rho_e/PQEDofT/dPQEDdT/rho_nu_extra are scalar-only (they
        # branch on ``Tg < me/_ELEC_THERMO_LOWT_RATIO`` or
        # ``cfg.DeltaNeff == 0.``), so evaluate them element-wise over T_out
        # -- same pattern as H_out above.
        rho_e_out      = np.array([thermo.rho_e(T)      for T in T_out])
        PQEDofT_out    = np.array([thermo.PQEDofT(T)    for T in T_out])
        dPQEDdT_out    = np.array([thermo.dPQEDdT(T)    for T in T_out])
        rho_nu_extra_out = np.array([thermo.rho_nu_extra(T) for T in T_out])

        rho_plasma = (thermo.rho_g(T_out) + rho_e_out
                      - PQEDofT_out + T_out * dPQEDdT_out)
        rho_nu_tot = (thermo.rho_nu(Tnue_out,   cfg.xi_nu_e)
                      + thermo.rho_nu(Tnumu_out,  cfg.xi_nu_mu)
                      + thermo.rho_nu(Tnutau_out, cfg.xi_nu_tau)
                      + rho_nu_extra_out)
        if self.rho_nu_SD is not None:
            Tnu_avg_out = ((Tnue_out**4 + Tnumu_out**4 + Tnutau_out**4) / 3.)**0.25
            rho_nu_tot = rho_nu_tot + self.rho_nu_SD(Tnu_avg_out)

        cols["rho_plasma [MeV^4]"] = rho_plasma
        cols["rho_nu_tot [MeV^4]"] = rho_nu_tot

        rho_tot = rho_plasma + rho_nu_tot
        if self.extra_rho:
            rho_extra = sum(rho_extra_fn(T_out) for rho_extra_fn in self.extra_rho)
            cols["rho_extra [MeV^4]"] = rho_extra
            rho_tot = rho_tot + rho_extra
        cols["rho_tot [MeV^4]"] = rho_tot

        return cols

    # ======================================================================
    # n <-> p weak rates
    # ======================================================================

    def _setup_weak_rates(self):
        cfg = self.cfg
        _t_weak0 = time.time()
        # Single forward and backward n<->p interpolant over the whole BBN
        # temperature range (the rate is continuous, so one grid suffices).
        self.weak_nTOp_frwrd_raw, self.weak_nTOp_bkwrd_raw = \
            primat_weak_rates.RecomputeWeakRates([self.Tg_vec, self.Tnue_vec], cfg,
                                        dFDneu_func=self.dFDneu_func,
                                        dFDneu_moments=self.dFDneu_moments)
        if cfg.debug:
            # Wording is generic on purpose: RecomputeWeakRates may have either
            # recomputed the rates (~2 s) or loaded them from a fingerprinted
            # cache file (~0 s) -- see primat.weak_rates.RecomputeWeakRates.
            print((f"[weak]  n <--> p weak rates ready in "
                   f"{time.time()-_t_weak0:.2f} s"), flush=True)

        # Normalisation factor: the stored rates are already in units of 1/tau_n
        # (ComputeFn was applied inside ComputeWeakRates), so multiplying by
        # 1/tau_n gives the actual rate in s^-1.  The absolute-normalisation
        # path (tau_n_normalization=False) still requires Fn to convert from
        # the 1/tau_n storage units to the GF-based normalisation.
        if cfg.tau_n_normalization:
            self._norm_weak_rates = 1. / cfg.tau_n   # [s^-1]
        else:
            Fn       = primat_weak_rates.ComputeFn(cfg)
            GFtilde2 = (cfg.GF * cfg.Vud)**2 * (1. + 3. * cfg.gA**2) / (2. * np.pi**3)
            self._norm_weak_rates = cfg.MeV_to_secm1 * (GFtilde2 * cfg.me**5) * Fn

    @property
    def NormWeakRates(self):
        """Normalisation factor applied to the raw n<->p weak-rate
        interpolants (settable -- the MC driver rescales it per-sample to
        propagate the neutron-lifetime uncertainty, see
        ``primat.main._mc_run_batch``)."""
        return self._norm_weak_rates

    @NormWeakRates.setter
    def NormWeakRates(self, value):
        self._norm_weak_rates = value

    def weak_nTOp_frwrd(self, T_K):
        """Normalised n -> p weak rate [s^-1] (see
        :meth:`Background.weak_nTOp_frwrd`)."""
        # Floor in raw (1/tau_n) units before normalising -- see
        # weak_nTOp_bkwrd below.  The quadratic interpolation used by
        # InterpolateWeakRates can overshoot slightly negative between
        # nodes even when the cached, already-clamped table values are all
        # >= 0, so the clamp is reapplied here rather than trusted to the
        # cache alone.
        return self._norm_weak_rates * _clamp_raw_weak_rate(self.weak_nTOp_frwrd_raw(T_K))

    def weak_nTOp_bkwrd(self, T_K):
        """Normalised p -> n weak rate [s^-1] (see
        :meth:`Background.weak_nTOp_bkwrd`)."""
        # Clamp to >= 0 and floor below the ~1e-28 (1/tau_n units) noise
        # level used in ComputeWeakRates (weak_rates.py).  At low T
        # (<~0.01 MeV) the tabulated/computed p->n rate is exp(-Q/T)-suppressed
        # down to ~1e-40 s^-1 and the sum of the phase-space correction terms
        # loses all its significant digits to cancellation, alternating sign
        # around zero (e.g. -9.3e-47 at 0.005 MeV); the quadratic
        # InterpolateWeakRates interpolant can also overshoot slightly
        # negative between cached nodes.  A rate is physically non-negative,
        # so this noise is replaced by 0 rather than left to feed a spurious
        # p->n *sink* of neutrons into the network.
        return self._norm_weak_rates * _clamp_raw_weak_rate(self.weak_nTOp_bkwrd_raw(T_K))

    # ======================================================================
    # Baryon sector
    # ======================================================================

    def _rhoB_of_a(self, a):
        """Baryon mass density at scale factor ``a`` [g cm^-3].

        ``n0CMB`` is the present-day CMB photon number density [MeV^3] and
        ``eta0b = n_B/n_gamma`` the baryon-to-photon ratio (from
        ``cfg.Omegabh2``); ``rho_B = m_B n_B`` with ``n_B = n0CMB eta0b / a^3``
        (comoving baryon number conservation).  Internal helper for
        :meth:`rhoB_BBN`.
        """
        cfg = self.cfg
        n0B = cfg.n0CMB * cfg.eta0b
        return cfg.ma * n0B * cfg.MeV4_to_gcmm3 / a**3  # [g cm^-3]

    def rhoB_BBN(self, t):
        """Baryon mass density rho_B(t) [g cm^-3] (see
        :meth:`Background.rhoB_BBN`)."""
        fast = self._a_of_t_scalar
        if fast is not None and np.ndim(t) == 0:
            return self._rhoB_of_a(fast(t))
        return self._rhoB_of_a(self.a_of_t(t))


class CustomBackground(Background):
    """Cosmological background driven by a user-supplied T(t), t, a(t) table.

    Reads the (T_γ, t, a) cosmological history from a file and uses the
    instantaneous-decoupling approximation for neutrino temperatures and n<->p
    weak rates.  This is the ``custom_background`` mode of
    :class:`~primat.config.PRIMATConfig`: ``cfg.incomplete_decoupling`` and
    ``cfg.spectral_distortions`` must both be ``False`` (enforced with warnings
    in ``PRIMATConfig.__init__`` if the caller forgot to set them).

    The file must be a tab- or comma-delimited text file with a header row
    naming the columns.  At minimum the columns **T** [MeV], **t** [s], and
    **a** (dimensionless scale factor, normalised so that ``a · T_γ → T0CMB``
    as ``T → 0``, i.e. ``a = 1`` today) must be present; extra columns are
    silently ignored.  Rows may be in any order; they are sorted internally by
    ascending cosmic time.

    **Neff estimation**: :meth:`rho_nu_total_final` estimates ``H`` at the
    final table point by numerical differentiation of ``ln a(t)``, then
    inverts the Friedmann equation ``H² = 8πG/3 · ρ_tot`` to obtain
    ``ρ_tot``.  The neutrino energy density follows as
    ``ρ_ν = ρ_tot − ρ_plasma``, and ``Neff = ρ_ν / ρ_γ / ((7/8)(4/11)^{4/3})``.
    A message is printed to inform the user that this indirect path was used.

    Parameters
    ----------
    cfg : PRIMATConfig
    plasma : primat.plasma.Plasma
    filename : str
        Path to the custom background file.

    Example
    -------
    The test suite in ``tests/test_custom_background.py`` writes a reference
    background to a temporary file (from a standard instantaneous-decoupling
    run) and checks that the custom-background run reproduces the same BBN
    observables to a relative error below ``1e-5``.
    """

    def __init__(self, cfg, plasma, filename):
        super().__init__(cfg, plasma)
        self._load_table(filename)
        self._setup_neutrino_history()
        self._setup_weak_rates()

    # ======================================================================
    # Table loading
    # ======================================================================

    def _load_table(self, filename):
        """Read T, t, a from the custom background file and build interpolants.

        The file is auto-detected as tab- or comma-delimited by inspecting the
        header line.  Rows are sorted by ascending cosmic time so that ``t`` is
        the natural independent variable for the interpolants.

        After this method, ``self.t_of_T``, ``self.T_of_t``, ``self.a_of_T``,
        ``self.T_of_a``, ``self.a_of_t``, and ``self.t_of_a`` are all set,
        and ``has_scale_factor`` is ``True``.  The raw sorted arrays are kept
        as ``_t_asc``, ``_T_by_t``, ``_a_by_t`` for :meth:`rho_nu_total_final`.

        Args:
            filename (str): path to the background file.

        Raises:
            ValueError: if required columns are missing or contain non-positive
                values.
        """
        # Auto-detect delimiter from the header.
        with open(filename) as fh:
            header_line = fh.readline()
        delimiter = '\t' if '\t' in header_line else ','

        raw = np.genfromtxt(filename, delimiter=delimiter, names=True, dtype=float)

        required = ('T', 't', 'a')
        missing  = [c for c in required if c not in raw.dtype.names]
        if missing:
            raise ValueError(
                f"custom_background file {filename!r} is missing required "
                f"columns: {missing}. Found: {list(raw.dtype.names)}"
            )

        T_raw = raw['T'].copy()
        t_raw = raw['t'].copy()
        a_raw = raw['a'].copy()

        if np.any(T_raw <= 0) or np.any(t_raw <= 0) or np.any(a_raw <= 0):
            raise ValueError(
                f"custom_background file {filename!r}: columns T, t, a "
                "must all be strictly positive."
            )

        # Sort by ascending t (earliest time = highest T first).
        idx_t  = np.argsort(t_raw)
        t_asc  = t_raw[idx_t]
        T_by_t = T_raw[idx_t]   # T_γ at each time (descending in T)
        a_by_t = a_raw[idx_t]   # a   at each time (ascending in a)

        # Sort by ascending T for T-based interpolants.
        idx_T  = np.argsort(T_raw)
        T_asc  = T_raw[idx_T]
        t_by_T = t_raw[idx_T]   # t at each T (decreasing t as T increases)
        a_by_T = a_raw[idx_T]   # a at each T (decreasing a as T increases)

        _kw = dict(bounds_error=False, fill_value='extrapolate', kind='linear')
        self.t_of_T = interp1d(T_asc,  t_by_T, **_kw)
        self.a_of_T = interp1d(T_asc,  a_by_T, **_kw)
        self.T_of_t = interp1d(t_asc,  T_by_t, **_kw)
        self.a_of_t = interp1d(t_asc,  a_by_t, **_kw)

        # a-based interpolants: a is ascending (a grows with time).
        idx_a  = np.argsort(a_by_t)
        a_sort = a_by_t[idx_a]
        self.T_of_a = interp1d(a_sort, T_by_t[idx_a], **_kw)
        self.t_of_a = interp1d(a_sort, t_asc[idx_a],  **_kw)

        self.has_scale_factor = True

        # Kept for rho_nu_total_final: time-ordered arrays.
        self._t_asc  = t_asc    # [s]
        self._T_by_t = T_by_t   # [MeV], decreasing
        self._a_by_t = a_by_t   # dimensionless, increasing

    # ======================================================================
    # Neutrino temperatures (instantaneous decoupling)
    # ======================================================================

    def _setup_neutrino_history(self):
        """Build instantaneous-decoupling T_ν(T_γ) and the grid for weak rates.

        Uses :class:`~primat.neutrino_history.InstantaneousDecoupling` (the
        same class used by ``StandardBackground`` when
        ``incomplete_decoupling=False``) to obtain T_νe(T_γ), T_νμ(T_γ),
        T_ντ(T_γ) from EM entropy conservation.

        The (Tg_vec, Tnue_vec) pair spans the T range of the loaded table and
        is passed to :func:`~primat.weak_rates.RecomputeWeakRates`.
        """
        from .neutrino_history import InstantaneousDecoupling
        cfg = self.cfg
        nh  = InstantaneousDecoupling(cfg, self.plasma)
        self._Tnue_of_Tg   = nh.Tnue_of_Tg
        self._Tnumu_of_Tg  = nh.Tnumu_of_Tg
        self._Tnutau_of_Tg = nh.Tnutau_of_Tg

        # Tg grid spanning the table's temperature range (low → high),
        # LOG-spaced -- as everywhere else in primat, and as the sizing helper
        # already assumes: n_points_per_decade returns
        # round(per_decade * log10(T_hi/T_lo)), a count of points *per decade*,
        # which only describes the grid if the points are laid out
        # logarithmically. This was a linear np.linspace, which over a table
        # spanning several decades starved the low-T end: for the standard
        # 40 MeV -> 0.001 MeV window (4.6 decades) it put 6 of its 2761 points
        # below 1e9 K, where log spacing puts 1161. That grid sets the node
        # spacing of the LINEAR T_nu(T_gamma) interpolant every n<->p rate
        # integrand reads -- the very quantity sampling_temperature_per_decade
        # is a weak-rate fingerprint field for (see weak_rates/cache.py's
        # _WEAK_RATE_BG_FIELDS: coarsening it moves the rates by ~1e-3).
        # Mirrored in primat-c/src/background.c's cpr_bg_init_custom.
        T_lo = float(np.min(self._T_by_t))
        T_hi = float(np.max(self._T_by_t))
        n_T_pts = primat_weak_rates.n_points_per_decade(cfg.sampling_temperature_per_decade, T_lo, T_hi)
        self.Tg_vec   = np.logspace(np.log10(T_lo), np.log10(T_hi), n_T_pts)
        self.Tnue_vec = self._Tnue_of_Tg(self.Tg_vec)

    # ======================================================================
    # n <-> p weak rates
    # ======================================================================

    def _setup_weak_rates(self):
        """Compute or load the n<->p weak rates for instantaneous decoupling.

        Delegates to :func:`~primat.weak_rates.RecomputeWeakRates` exactly
        as :meth:`StandardBackground._setup_weak_rates` does; ``dFDneu_func``
        is ``None`` because ``spectral_distortions=False`` is enforced by
        :class:`~primat.config.PRIMATConfig` when ``custom_background`` is set.
        """
        cfg  = self.cfg
        _t0  = time.time()
        self.weak_nTOp_frwrd_raw, self.weak_nTOp_bkwrd_raw = \
            primat_weak_rates.RecomputeWeakRates([self.Tg_vec, self.Tnue_vec], cfg,
                                        dFDneu_func=None)
        if cfg.debug:
            print(f"[weak-py] n <--> p weak rates ready in "
                  f"{time.time()-_t0:.2f} s", flush=True)

        if cfg.tau_n_normalization:
            self._norm_weak_rates = 1. / cfg.tau_n
        else:
            Fn       = primat_weak_rates.ComputeFn(cfg)
            GFtilde2 = (cfg.GF * cfg.Vud)**2 * (1. + 3. * cfg.gA**2) / (2. * np.pi**3)
            self._norm_weak_rates = cfg.MeV_to_secm1 * (GFtilde2 * cfg.me**5) * Fn

    @property
    def NormWeakRates(self):
        """Normalisation factor for n<->p rates (settable for MC tau_n scan).

        Mirrors :attr:`StandardBackground.NormWeakRates` so that
        :func:`~primat.main._mc_run_batch` can rescale rates per sample
        without recomputing the expensive weak-rate integrals.
        """
        return self._norm_weak_rates

    @NormWeakRates.setter
    def NormWeakRates(self, value):
        self._norm_weak_rates = value

    def weak_nTOp_frwrd(self, T_K):
        """Normalised n -> p weak rate [s^-1] (see :meth:`Background.weak_nTOp_frwrd`)."""
        # Floor in raw (1/tau_n) units before normalising: see the sibling
        # implementation above for why this is needed even on top of the
        # clamp already applied inside ComputeWeakRates.
        return self._norm_weak_rates * _clamp_raw_weak_rate(self.weak_nTOp_frwrd_raw(T_K))

    def weak_nTOp_bkwrd(self, T_K):
        """Normalised p -> n weak rate [s^-1] (see :meth:`Background.weak_nTOp_bkwrd`)."""
        # Clamp to >= 0 and floor below the numerical-noise level: see the
        # sibling implementation above -- the low-T p->n rate is dominated by
        # cancellation noise that can go negative, which is unphysical for a
        # rate.
        return self._norm_weak_rates * _clamp_raw_weak_rate(self.weak_nTOp_bkwrd_raw(T_K))

    # ======================================================================
    # Baryon sector
    # ======================================================================

    def rhoB_BBN(self, t):
        """Baryon mass density rho_B(t) [g cm^-3].

        Uses comoving baryon-number conservation: rho_B = m_B · n_B =
        m_B · n0CMB · eta0b / a(t)^3, with ``a(t)`` read from the supplied
        table.  The scale factor must be normalised so that ``a = 1`` today
        (``a · T_γ → T0CMB_MeV`` as ``T_γ → 0``), consistent with the
        standard primat convention (see
        :meth:`StandardBackground._rhoB_of_a`).

        Args:
            t (float or array): cosmic time [s].

        Returns:
            float or array: baryon mass density [g cm^-3].
        """
        cfg = self.cfg
        n0B = cfg.n0CMB * cfg.eta0b
        a   = self.a_of_t(t)
        return cfg.ma * n0B * cfg.MeV4_to_gcmm3 / a**3

    # ======================================================================
    # Neff via Friedmann equation
    # ======================================================================

    def rho_nu_total_final(self):
        """Estimate (T_γ_final, ρ_ν_final) from the Friedmann equation.

        The Hubble rate at the last row of the supplied table is obtained by
        numerical differentiation of ``ln a(t)``:

            H = d(ln a)/dt ≈ Δ(ln a) / Δt   (one-sided backward difference
                                               at the final point)

        The total energy density follows from the Friedmann equation
        ``H² = 8πG/3 · ρ_tot``, and the neutrino contribution is isolated as

            ρ_ν = ρ_tot − ρ_plasma(T_γ)

        where ``ρ_plasma`` is the photon + e⁺e⁻ + QED-correction energy
        density evaluated analytically at the final photon temperature.
        :meth:`Background.N_eff` then converts ``ρ_ν`` to ``Neff``.

        A message is always printed when this method is called so the user
        knows that the indirect Friedmann path was used.

        Returns:
            tuple: ``(Tg_final, rho_nu_total)`` — photon temperature [MeV]
                and total neutrino energy density [MeV^4] at the final time.
        """
        cfg    = self.cfg
        thermo = self.plasma

        t_arr = self._t_asc
        a_arr = self._a_by_t
        T_arr = self._T_by_t

        # Estimate H = d(ln a)/dt at the final table point.
        #
        # A simple endpoint finite difference on ln(a) vs t is only first-order
        # accurate and produces a ~1% bias in H that amplifies to a ~3% bias in
        # Neff (because rho_nu = rho_tot - rho_plasma with rho_plasma ~ rho_tot/2,
        # so a small relative error in rho_tot doubles into rho_nu).
        #
        # Instead we exploit that in the late radiation era a(t) follows a
        # power law: a ∝ t^p with p ≈ 0.5 exactly (radiation domination).
        # Fitting ln(a) = p·ln(t) + q over the last N_fit points gives a
        # slope p whose accuracy is limited only by the table range, not by
        # the local step size at the endpoint. Then H = p / t[-1].
        # N_fit = min(50, half the table) balances a wide baseline (reduces
        # noise) against non-stationarity (p evolves very slowly near T_end).
        N_fit   = min(50, len(t_arr) // 2)
        lna_fit = np.log(a_arr[-N_fit:])
        lnt_fit = np.log(t_arr[-N_fit:])
        # np.polyfit([1]): linear fit ln(a) = p·ln(t) + q
        p_slope = float(np.polyfit(lnt_fit, lna_fit, 1)[0])
        H_final = p_slope / t_arr[-1]   # [s^-1]
        Tg_f    = float(T_arr[-1])      # [MeV] — smallest T in the table

        # Friedmann: H [s^-1] = MeV_to_secm1 · sqrt(8π/(3 Mpl^2) · ρ_tot [MeV^4])
        # Inverted: ρ_tot = 3 Mpl^2 / (8π) · (H / MeV_to_secm1)^2
        H_MeV   = H_final / cfg.MeV_to_secm1
        rho_tot = 3. * cfg.Mpl**2 / (8. * np.pi) * H_MeV**2

        # Plasma energy density at Tg_f: photons + electrons + QED corrections.
        # (Eq. for rho_plasma: Phys. Rep. Eq. 19 / StandardBackground.Hubble.)
        rho_plasma = (thermo.rho_g(Tg_f) + thermo.rho_e(Tg_f)
                      - thermo.PQEDofT(Tg_f) + Tg_f * thermo.dPQEDdT(Tg_f))

        rho_nu_tot = rho_tot - rho_plasma

        print(
            f"[custom_background] Neff from Friedmann H²=8πG/3·ρ_tot "
            f"at Tγ = {Tg_f:.4e} MeV: "
            f"H = {H_final:.6e} s⁻¹, "
            f"ρ_tot = {rho_tot:.6e} MeV⁴, "
            f"ρ_plasma = {rho_plasma:.6e} MeV⁴, "
            f"ρ_ν = {rho_nu_tot:.6e} MeV⁴"
        )

        return Tg_f, rho_nu_tot
