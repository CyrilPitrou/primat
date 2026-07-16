# -*- coding: utf-8 -*-
"""
weak_rates.fast_eval — numba-JIT scalar evaluation of the n<->p weak rates
===============================================================================

The BDF solver evaluates the raw n->p and p->n weak rates (in units of 1/τ_n)
at a *scalar* photon temperature tens of thousands of times per BBN run, once
per ``NetworkDefinition.fill_buffer`` call.  Each evaluation is the sum of two
tabulated interpolants built in :mod:`weak_rates.api` /
:mod:`weak_rates.corrections`:

  * the non-thermal rate (Born+FM+CCR+SD), a **cubic** spline in log10-log10
    space (:func:`weak_rates.api._weak_rate_loglog_interp`), and
  * the finite-temperature CCRTh correction, a **quadratic** spline in linear
    space (:func:`weak_rates.corrections._thermal_correction_interpolants`),
    pinned to 0 below ``T_CCRTH_MIN``.

Evaluated through scipy's ``interp1d.__call__`` these two scalar lookups
dominated the pure-Python backend's solve time (each ``interp1d`` call pays a
fixed ~10 µs of ``asarray``/input-validation overhead before doing a trivial
spline evaluation).  This module removes that overhead on the hot *scalar* path
by extracting each interpolant's piecewise-polynomial (``PPoly``) breakpoints
and coefficients **once at setup** and evaluating them with a tiny
``@njit`` Horner loop -- the *same* spline, just without the per-call Python/
NumPy glue.

Numerical fidelity
------------------
The interpolants are still *fitted* by scipy exactly as before (this module
never changes the fitting scheme -- the log10-log10 not-a-knot cubic is a
cross-backend parity contract, see ``api._weak_rate_loglog_interp``); only the
*evaluation* is re-expressed.  Converting a fitted spline to ``PPoly`` and
evaluating it by Horner differs from scipy's B-spline evaluation only at the
~1e-15 (relative) rounding level, far below the ±3e-9 D/H same-backend
regression tolerance.  To be safe against scipy-version quirks, :class:`FastWeakRate`
runs a build-time self-check comparing the JIT path against the original scipy
closure over the whole temperature grid and **silently falls back** to the
scipy closure (``_ok = False``) if the agreement is worse than ``_SELFCHECK_RTOL``.
Array-valued queries (used by table dumps / time-evolution output, never on the
hot solver path) always take the original scipy closure.

If numba is not installed the ``@njit`` decorator degrades to a no-op and the
scalar path runs as plain Python (still correct, just without the speed-up);
callers are unaffected either way.
"""

import numpy as np

try:                                                    # numba is recommended, not required
    from numba import njit as _njit
    _HAS_NUMBA = True
except Exception:                                       # pragma: no cover - numba absent
    _HAS_NUMBA = False

    def _njit(*args, **kwargs):
        # Mimic numba.njit's two call forms: @_njit and @_njit(cache=True).
        if args and callable(args[0]):
            return args[0]
        return lambda f: f


# Relative tolerance for the build-time JIT-vs-scipy self-check.  Well above the
# ~1e-15 PPoly-vs-Bspline rounding difference, well below the 3e-9 D/H pin, so
# it only ever trips on a genuine extraction failure (e.g. a scipy change to the
# interp1d spline internals), in which case we fall back to scipy transparently.
_SELFCHECK_RTOL = 1e-11


@_njit(cache=True)
def _ppoly_scalar(xb, c, xq):
    """Evaluate a scipy ``PPoly`` at a single point ``xq`` (extrapolating).

    ``xb`` are the ``PPoly.x`` breakpoints (ascending, possibly with repeated
    boundary knots) and ``c`` the ``PPoly.c`` coefficient array of shape
    ``(order+1, nintervals)``; the local polynomial on interval ``i`` is
    ``sum_m c[m, i] * (xq - xb[i])**(order-m)``, evaluated here by Horner.
    Queries outside ``[xb[0], xb[-1]]`` reuse the first/last interval's
    polynomial (matching ``PPoly(..., extrapolate=True)`` and the interpolants'
    ``fill_value='extrapolate'``).
    """
    n = xb.shape[0]
    # First breakpoint strictly greater than xq, minus one -> bracketing
    # interval; clamp so out-of-range points extrapolate off the edge cell.
    i = int(np.searchsorted(xb, xq, side='right')) - 1
    if i < 0:
        i = 0
    elif i > n - 2:
        i = n - 2
    dx = xq - xb[i]
    r = 0.0
    for m in range(c.shape[0]):
        r = r * dx + c[m, i]
    return r


@_njit(cache=True)
def _weak_raw_scalar(cx, cc, T_zero_below, has_th, tx, tc, t_min, th_scale, Tq):
    """Full raw weak rate at scalar photon temperature ``Tq`` [K], in 1/τ_n.

    Reproduces ``lambda T: nonthermal(T) + thermal(T)`` from
    :func:`weak_rates.api.RecomputeWeakRates`:

      * non-thermal: ``10**cubic(log10 max(Tq, 1e-300))``, forced to 0 below
        ``T_zero_below`` (the backward rate's clamped-to-zero low-T prefix;
        ``-inf`` for the strictly-positive forward rate so the mask never fires);
      * thermal (only if ``has_th``): ``th_scale × quadratic(Tq)`` (the CCRTh
        correction divided by the neutron-decay factor Fn), pinned to 0 below
        ``t_min`` (``T_CCRTH_MIN``) exactly as the scipy closure does.
    """
    if Tq >= T_zero_below:
        lt = Tq if Tq > 1e-300 else 1e-300   # log10 undefined at/below 0
        val = 10.0 ** _ppoly_scalar(cx, cc, np.log10(lt))
    else:
        val = 0.0
    if has_th and Tq >= t_min:
        val += th_scale * _ppoly_scalar(tx, tc, Tq)
    return val


def _ppoly_from_interp1d(interp):
    """Extract ``(breakpoints, coeffs)`` of a scipy ``interp1d`` spline.

    Works for any spline ``kind`` (``'cubic'``, ``'quadratic'``, ...): reads the
    underlying B-spline (``interp._spline``), squeezes its coefficient array to
    1-D (interp1d stores it as ``(ncoef, 1)``), and converts to the equivalent
    piecewise polynomial via ``PPoly.from_spline``.  The conversion is exact up
    to floating-point rounding, so the returned arrays evaluate (by
    :func:`_ppoly_scalar`) to the same spline scipy would.
    """
    from scipy.interpolate import PPoly
    sp = interp._spline                      # BSpline backing the interp1d
    k = sp.k
    nc = sp.t.shape[0] - k - 1               # number of B-spline coefficients
    c1 = sp.c[:nc, 0] if sp.c.ndim == 2 else sp.c[:nc]
    pp = PPoly.from_spline((np.ascontiguousarray(sp.t),
                            np.ascontiguousarray(c1), k))
    return np.ascontiguousarray(pp.x), np.ascontiguousarray(pp.c)


class FastWeakRate:
    """Fast scalar evaluator for one raw n<->p weak-rate channel.

    Constructed from the same scipy interpolant objects the original closures
    captured (attached to those closures as ``.spline``/``.T_zero_below`` and
    ``.interp`` -- see :mod:`weak_rates.api`/:mod:`weak_rates.corrections`), plus
    the original combined closure ``orig`` used as the array-input / fallback
    path.  Calling the instance with a scalar temperature takes the JIT path;
    calling it with an array (or after a failed self-check) delegates to ``orig``.

    Parameters
    ----------
    nt_eval : the non-thermal closure, carrying ``.spline`` (interp1d) and
              ``.T_zero_below`` (float, ``-inf`` for the forward rate).
    th_eval : the thermal closure, carrying ``.interp`` (interp1d) and ``.scale``
              (the constant 1/Fn factor) when ``cfg.thermal_corrections`` is on,
              or neither attribute when it is off.
    t_min   : ``T_CCRTH_MIN`` [K] -- thermal floor below which the CCRTh term is 0.
    orig    : the original ``lambda T: nt_eval(T) + th_eval(T)`` closure.
    """

    def __init__(self, nt_eval, th_eval, t_min, orig):
        self._orig = orig
        self._ok = False
        try:
            self.cx, self.cc = _ppoly_from_interp1d(nt_eval.spline)
            self.T_zero_below = float(nt_eval.T_zero_below)
            th_interp = getattr(th_eval, "interp", None)
            if th_interp is not None:
                self.tx, self.tc = _ppoly_from_interp1d(th_interp)
                self.t_min = float(t_min)
                self.th_scale = float(th_eval.scale)
                self.has_th = True
            else:
                # Dummy arrays keep the njit signature monomorphic; never read
                # because has_th gates every thermal access.
                self.tx = np.zeros(2)
                self.tc = np.zeros((1, 1))
                self.t_min = 0.0
                self.th_scale = 0.0
                self.has_th = False
            self._ok = self._selfcheck()
        except Exception:
            # Any extraction failure -> stay on the scipy closure. Correctness
            # over speed; the fallback is exactly the previous behaviour.
            self._ok = False

    def _scalar(self, Tq):
        return _weak_raw_scalar(self.cx, self.cc, self.T_zero_below,
                                self.has_th, self.tx, self.tc, self.t_min,
                                self.th_scale, Tq)

    def _selfcheck(self):
        """Return True iff the JIT path matches the scipy closure to
        ``_SELFCHECK_RTOL`` over the fitted temperature domain.

        The check grid is built from the interpolants' own breakpoints so it is
        automatically in the right units (Kelvin): the non-thermal cubic lives
        in ``log10 T`` space, so its knots map to ``10**cx``; the thermal
        quadratic lives in linear ``T`` space (``tx``).  Interval midpoints are
        added to catch a between-nodes divergence, not just agreement at nodes.
        """
        T_nt = np.power(10.0, self.cx)
        pts = [T_nt, 0.5 * (T_nt[:-1] + T_nt[1:])]
        if self.has_th:
            pts.append(self.tx)
            pts.append(0.5 * (self.tx[:-1] + self.tx[1:]))
        T = np.concatenate(pts)
        T = T[np.isfinite(T) & (T > 0.0)]
        if T.size == 0:
            return False
        for Tq in T:
            ref = float(self._orig(Tq))
            got = self._scalar(float(Tq))
            denom = abs(ref) if ref != 0.0 else 1.0
            if abs(got - ref) > _SELFCHECK_RTOL * denom:
                return False
        return True

    def __call__(self, T):
        # Hot path: scalar query -> JIT. Array query (dumps/evolution output) or
        # a failed self-check -> original scipy closure.
        if self._ok and np.ndim(T) == 0:
            return self._scalar(float(T))
        return self._orig(T)
