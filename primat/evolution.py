# -*- coding: utf-8 -*-
"""
evolution.py
============
Unified time-evolution schema: the shared columns both
the Python and C backends populate, so a notebook can plot nuclide evolution
from *either* backend's output without caring which one ran.

This module docstring is the **authoritative contract** for that schema
(column names, order, and semantics): both backends' writers must conform to
it, and any new column must be added to both before it is considered part of
the schema.

``EvolutionResult`` is the in-memory primary artifact -- populated by
``NuclearNetwork.solve()`` as ``self.evolution`` (and surfaced as
``PRIMAT.solve()``'s returned dict's ``"evolution"`` key) whenever
``cfg.output_time_evolution=True``, with *no* disk I/O required to get it.
Disk output (the ``cfg.output_file`` TSV) is a derived convenience: it is
produced by calling :func:`dump_evolution` on that same object, not
something the solver opens a file for directly.

Columns (tab-separated, ``#``-free header line -- see :data:`_CORE_COLUMNS`):
``t_s  a  T_gamma_MeV  T_nue_MeV  T_numu_MeV  T_nutau_MeV  Y_<nuclide> ...``
The ``Y_<nuclide>`` block is network-dependent (small/large have different
nuclide lists), so the header line is the source of truth -- :func:`load_evolution`
reads it dynamically rather than assuming a fixed column count.

Per-reaction forward-rate columns (``<reaction>_frwrd``) are an **optional**
trailing block, appended after the ``Y_<nuclide>`` block when
``cfg.output_rates_time_evolution=True`` -- one column per LT reaction
except n<->p, which has no rate table (12 for
``small``/``small_parthenope``, 67 for ``large``+``amax=8``, 428 for full
``large``). They carry the forward
reaction-rate interpolant [same table units as the shipped nuclear rates]
evaluated at each row's photon temperature, and live in
:attr:`EvolutionResult.rates` (``None`` when the flag is off). Both backends
emit the identical column names in the identical (lexicographically sorted)
order, pinned by ``tests/test_backend_parity.py``. The ``a``/``T_nu`` columns
are ``np.nan`` when the active background has no scale-factor/neutrino-sector
tracking (e.g. a minimal custom background).
"""
from __future__ import annotations

import io
import os
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
from numpy.typing import NDArray

# Column order for the always-present, cross-backend-comparable block.  The
# Y_<nuclide> block (network-dependent) follows these six.
_CORE_COLUMNS = ("t_s", "a", "T_gamma_MeV", "T_nue_MeV", "T_numu_MeV", "T_nutau_MeV")
_Y_PREFIX = "Y_"
# Suffixes marking the optional trailing per-reaction rate columns. The
# legacy/current writers emit only forward-rate columns (``_frwrd``); ``_bkwrd``
# is reserved so a future backward-rate column would round-trip through
# ``load_evolution`` into ``rates`` without a schema change.
_RATE_SUFFIXES = ("_frwrd", "_bkwrd")


@dataclass
class EvolutionResult:
    """In-memory unified time-evolution result.

    Attributes
    ----------
    t : np.ndarray
        Cosmic time [s].
    a : np.ndarray
        Scale factor (``np.nan`` everywhere if the background has no
        scale-factor relation, e.g. a minimal custom background).
    T_gamma : np.ndarray
        Photon temperature [MeV].
    T_nu : dict of str -> np.ndarray
        Per-flavour neutrino temperature [MeV], keyed ``"e"``/``"mu"``/
        ``"tau"`` (``np.nan`` arrays if the background tracks no neutrino
        sector).
    Y : dict of str -> np.ndarray
        Per-nuclide abundance per baryon Y_i = n_i/n_b, keyed by nuclide name, in
        network order (``n``/``p`` first).
    rates : dict of str -> np.ndarray, optional
        Optional per-reaction forward-rate columns (populated when
        ``cfg.output_rates_time_evolution=True`` -- one per reaction in the
        active LT network; ``None`` otherwise). Keyed by column
        name ``<reaction>_frwrd`` (canonical rate syntax, e.g.
        ``n_p__d_g_frwrd``), each value an array aligned with :attr:`t`,
        holding the forward reaction-rate interpolant at that step's photon
        temperature. Serialised as the trailing column block after ``Y_`` (see
        the module docstring); both backends emit the identical names in the
        identical sorted order.
    """
    t: NDArray[np.float64]
    a: NDArray[np.float64]
    T_gamma: NDArray[np.float64]
    T_nu: dict[str, NDArray[np.float64]]
    Y: dict[str, NDArray[np.float64]] = field(default_factory=dict)
    rates: dict[str, NDArray[np.float64]] | None = None


def dump_evolution(result: EvolutionResult, path: str | None = None) -> str:
    """Serialise ``result`` to the shared TSV schema (module docstring).

    Always returns the TSV text; additionally writes it to ``path`` if
    given (relative paths resolve against the current working directory,
    like ``cfg.output_file`` elsewhere in this package). Called by:
    ``NuclearNetwork._write_time_evolution`` (the ``cfg.output_file``
    convenience path); and ``primat-gui``'s download buttons, which call
    this lazily on ``run.evolution`` to produce the file text for
    ``st.download_button(data=...)`` -- never via a tempfile.

    Parameters
    ----------
    result : EvolutionResult
    path : str, optional

    Returns
    -------
    str
        The TSV text (header line + one row per time step).
    """
    names = list(_CORE_COLUMNS) + [_Y_PREFIX + s for s in result.Y]
    columns = [result.t, result.a, result.T_gamma,
               result.T_nu["e"], result.T_nu["mu"], result.T_nu["tau"]] + list(result.Y.values())
    # Optional per-reaction rate block, appended AFTER the Y_ columns so
    # the default (rates=None) schema stays byte-identical. Insertion order of
    # the dict is the on-disk column order (the writer already sorts it).
    if result.rates:
        names += list(result.rates.keys())
        columns += list(result.rates.values())
    data = np.column_stack(columns)

    buf = io.StringIO()
    # comments='' (no leading "# ") matches the convention already used by
    # the rest of this package's TSV writers (e.g. background.py's
    # write_time_evolution), so a plain `header.split("\t")` recovers the
    # column names without stripping a comment marker first.
    np.savetxt(buf, data, delimiter='\t', header="\t".join(names), comments='')
    text = buf.getvalue()

    if path is not None:
        out_path = os.path.abspath(path)
        out_dir = os.path.dirname(out_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(out_path, 'w') as f:
            f.write(text)

    return text


def Y_interpolator(result: EvolutionResult, name: str) -> Callable[[float | NDArray[np.float64]], Any]:
    """Build a ``Y(t)`` callable for nuclide ``name`` from ``result.t``/
    ``result.Y[name]`` alone -- the backend-agnostic counterpart of
    ``primat.main.PRIMAT.__getitem__``'s live SciPy interpolator, usable on a
    plain :class:`EvolutionResult` from *either* backend's :func:`primat.backend.run_bbn`
    (no live ``PRIMAT``/``NuclearNetwork`` object required).

    Same convention as the live interpolator built by
    ``NuclearNetwork.solve()`` (``nuclear_network.py``'s ``self.Y_of_t``):
    piecewise-linear, with ``fill_value=(0, Y[-1])`` so a query before
    ``result.t[0]`` reads as zero abundance (not yet produced) and a query
    after ``result.t[-1]`` holds at the final value (no decay beyond the
    integrated era -- see ``run_bbn``'s ``decay_era`` gap, module docstring
    of ``primat.backend``, for the one case this does *not* cover).

    Args:
        result: EvolutionResult.
        name: str. A key of ``result.Y``.

    Returns:
        callable: ``t -> Y`` (accepts a scalar or array ``t`` [s]).
    """
    from scipy.interpolate import interp1d
    Y = result.Y[name]
    return interp1d(result.t, Y, bounds_error=False, fill_value=(0.0, Y[-1]))


def T_gamma_interpolator(result: EvolutionResult) -> Callable[[float | NDArray[np.float64]], Any]:
    """Build a ``T_gamma(t)`` [MeV] callable from ``result.t``/``result.T_gamma``
    alone -- the backend-agnostic counterpart of ``primat.main.PRIMAT.T_of_t``,
    usable on a plain :class:`EvolutionResult` from either backend.

    Piecewise-linear, clamped to the first/last sampled temperature outside
    ``[result.t[0], result.t[-1]]`` (``T_gamma`` decreases monotonically with
    ``t``, so this is the high-T/low-T extrapolation respectively).

    Args:
        result: EvolutionResult.

    Returns:
        callable: ``t -> T_gamma`` [MeV] (accepts a scalar or array ``t`` [s]).
    """
    from scipy.interpolate import interp1d
    T = result.T_gamma
    return interp1d(result.t, T, bounds_error=False, fill_value=(T[0], T[-1]))


def t_of_T_interpolator(result: EvolutionResult) -> Callable[[float | NDArray[np.float64]], Any]:
    """Build a ``T_gamma -> t`` [s] inverse-lookup callable, backend-agnostic.

    The counterpart of :func:`T_gamma_interpolator`, useful e.g. to add a
    secondary x-axis labelled in ``T_gamma`` on a plot whose primary axis is
    cosmic time ``t`` (see ``notebooks/AbundanceEvolution.ipynb``). Works from
    a loaded ``EvolutionResult``, so unlike ``primat.main.PRIMAT.t_of_T`` it
    needs neither a live instance nor the Python backend.

    ``T_gamma`` decreases monotonically with ``t`` over a BBN run, so the
    ``(T_gamma, t)`` pairs are reversed into ascending order and interpolated
    with :func:`numpy.interp` (clamped at the ends, matching the monotone
    extrapolation behaviour of :func:`T_gamma_interpolator`).

    Args:
        result: EvolutionResult.

    Returns:
        callable: ``T_gamma [MeV] -> t [s]`` (accepts a scalar or array).
    """
    T_asc = result.T_gamma[::-1]
    t_asc = result.t[::-1]
    return lambda T: np.interp(T, T_asc, t_asc)


def load_evolution(path: str) -> EvolutionResult:
    """Parse the shared TSV schema written by either backend's
    ``dump_evolution``/equivalent writer, returning the same
    :class:`EvolutionResult` structure as ``solve()``'s in-memory
    ``run.evolution`` -- for the case of reloading a previously-saved run
    without re-solving.

    Trailing per-reaction rate columns (``*_frwrd``/``*_bkwrd``, written when
    ``output_rates_time_evolution=True``) are collected into
    :attr:`EvolutionResult.rates` (``None`` if the file has none). Any other
    column beyond the core + ``Y_<nuclide>`` + rate blocks (e.g. a
    backend-specific bonus column) is simply ignored, so a file containing
    extra columns is still loadable.

    Parameters
    ----------
    path : str

    Returns
    -------
    EvolutionResult
    """
    with open(path) as f:
        header = f.readline().strip().split("\t")
        data = np.loadtxt(f)
    if data.ndim == 1:
        data = data[np.newaxis, :]

    col = {name: data[:, i] for i, name in enumerate(header)}
    Y = {name[len(_Y_PREFIX):]: arr for name, arr in col.items()
         if name.startswith(_Y_PREFIX)}
    # Trailing forward/backward per-reaction rate columns -> rates (in the
    # file's column order); None when the flag was off so callers can test it
    # as a plain truthiness/None flag.
    rates = {name: arr for name, arr in col.items()
             if name.endswith(_RATE_SUFFIXES)}

    return EvolutionResult(
        t=col["t_s"], a=col["a"], T_gamma=col["T_gamma_MeV"],
        T_nu={"e": col["T_nue_MeV"], "mu": col["T_numu_MeV"], "tau": col["T_nutau_MeV"]},
        Y=Y,
        rates=rates or None,
    )
