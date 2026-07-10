# -*- coding: utf-8 -*-
"""
primat.backend
================

Dispatch layer choosing between the compiled C extension
(``primat._primat_c``, wrapping ``primat-c``'s ``cprimat_run``) and the
pure-Python solver (``primat.main.PRIMAT``).

``HAS_C_BACKEND`` is probed once at import time (``True`` iff the extension
built successfully -- see ``setup.py``'s ``optional_build_ext``, which lets
``pip install`` succeed even without a C compiler). :func:`run_bbn` is the
single dispatch entry point; everything else in this module supports it.

Feature gap (the one remaining C-unsupported ``PRIMAT.__init__``
extension):

* ``background=`` (a custom :class:`primat.background.Background` object) --
  an inherently-Python extension point (arbitrary user Python subclassing the
  background), with no way to cross the C ABI. A non-``None`` ``background=``
  always forces the Python backend under ``force_backend in (None, "auto")``,
  and raises ``ValueError`` under ``force_backend="c"``.

The former ``extra_rho`` and ``decay_era`` gaps are now *closed* (O-8):

* ``extra_rho`` (extra Friedmann energy-density callables) is supported on
  the C backend via a tabulated handoff -- :func:`_tabulate_extra_rho`
  evaluates the summed ``rho(Tg)`` on a dense log-Tg grid and passes the
  ``(Tg[], rho[])`` arrays to the C extension, which splines them and adds
  ``rho(Tg)`` inside ``cpr_bg_Hubble`` (see ``primat-c/src/background.c`` and
  ``config.h``'s ``extra_rho_*`` fields). Both backends agree to the
  cross-backend tolerance.

* ``decay_era`` (the long-lived-isotope Decay-Time era past ``T_end``) is
  ported: ``cpr_nuclear_network_decay_era`` (``primat-c/src/nuclear_network.c``)
  mirrors ``_integrate_decay_era``'s matrix-exponential decay propagation
  (scaling-and-squaring Padé-13). It changes no result-dict observable on
  either backend (``Y_final`` is the end-of-LT state); its only output is the
  optional ``output_decay_evolution`` TSV, which both backends write in the
  identical schema.

Set ``PRIMAT_BACKEND_LOG=1`` in the environment (or call with
``log_backend=True``) to print, on every :func:`run_bbn`/:func:`run_mc` call,
which backend actually ran and why -- chiefly to catch a silent
``force_backend="auto"`` fallback to Python (e.g. because a C-unsupported
feature was requested, or the extension failed to build) during development.

``custom_network`` (the GUI "Customise Reactions" override: removed/replaced/
added reactions plus rate-table overrides) *is* supported on both backends:
``primat-c``'s ``cprimat_run``/``cpr_mc_uncertainty`` take an optional
``CPRCustomNetwork*`` (``primat-c/include/network_data.h``), and
``primat/_primat_c/_wrapper.c`` parses the same dict shape
(``UpdateNuclearRates``/``kept_to_custom_network``, see
``primat/network_data.py``/``primat/gui/custom_rates.py``) into one. It is no
longer part of ``python_only_feature`` below.

``output_time_evolution=True`` *is* supported on both backends: the C
extension's ``cprimat_run`` populates ``CPRResults``'s
``evol_*`` in-memory arrays (``primat-c/include/api.h``) and
``primat/_primat_c/_wrapper.c`` hands them back as an ``"evolution"`` dict
key (plain Python lists, no numpy C-API dependency in the extension); this
module assembles the same :class:`primat.evolution.EvolutionResult` shape
the Python backend produces, with no disk I/O on either backend's part.

``data_dir``/``user_nuclear_dir`` (see CLAUDE.md's "Rates directory
resolution" section) *are* supported on both backends: ``data_dir`` fully
replaces the shipped data tree; ``user_nuclear_dir`` is an additive overlay
for nuclear networks and rate tables.  They are ordinary ``params`` dict keys
applied generically via ``cpr_config_set_by_name`` on the C side, so no
special-casing is needed here — except that ``data_dir`` must also be
forwarded as the ``data_dir`` positional argument to ``_c_ext.run_bbn``/
``_c_ext.run_mc`` (the C extension's ``cpr_config_init_defaults`` takes the
data folder there rather than via ``cpr_config_set_by_name``).

:func:`run_mc` is the MC counterpart of :func:`run_bbn`: it dispatches between
``primat._primat_c``'s ``run_mc`` (wrapping ``primat-c/src/mc.c``'s threaded
``cpr_mc_uncertainty``) and ``primat.main.mc_uncertainty`` (joblib), returning
the same :class:`primat.main.MCResult` shape either way -- the "common
language" the two backends share for MC results (CLAUDE.md's backend-parity
mandate). The C path uses a pthread/xoshiro256** RNG, *not* NumPy's
``default_rng``, so individual samples are not bit-for-bit comparable across
backends (only statistically, mean/std convergence -- see ``mc.h``).

``prev`` (incremental sample reuse) *is* supported on the C path, mirroring
``cpr_mc_uncertainty``'s ``prev_centrals``/``prev_values`` parameters (see
``mc.h``): :func:`run_mc` checks the same reuse-guard ``mc_uncertainty`` does
internally (seed/quantities/params/custom_network all matching), plus one
more condition the C side cannot check for itself -- ``prev.backend`` must
equal the backend about to compute the extension, since the two backends'
RNG streams are not interchangeable. A ``prev`` that fails the guard (e.g.
computed by the other backend) is silently ignored, exactly like
``mc_uncertainty``'s own fallback -- never an error, and never a forced
backend switch. ``custom_network`` is supported on both backends, same as
:func:`run_bbn`.
"""
from __future__ import annotations

import os
import sys
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from .main import MCResult

__all__ = ["HAS_C_BACKEND", "run_bbn", "run_mc", "dump_mc_samples",
           "dump_mc_covariance", "dump_mc_correlation", "dump_final_with_sigma"]


def _log_backend(func_name: str, used: str, reason: str, log_backend: bool) -> None:
    """Print which backend ``func_name`` (``"run_bbn"``/``"run_mc"``) actually
    used, plus why, when asked to via ``log_backend=True`` or the
    ``PRIMAT_BACKEND_LOG`` environment variable (module docstring). Printed to
    stderr (not stdout) so it never pollutes a CLI's piped result output.
    """
    if log_backend or os.environ.get("PRIMAT_BACKEND_LOG"):
        print(f"[primat.backend] {func_name}: used {used} backend ({reason})",
              file=sys.stderr)

# Standard derived observables, unconditionally merged into every MC result
# (alongside every tracked nuclide's final Y -- see mc_uncertainty/_c_mc)
# regardless of what the caller explicitly requested via run_mc's
# `quantities` argument, so an MCResult is always complete enough to dump to
# disk via dump_mc_samples (the CLI's output_mc_samples/output_mc_file_prefix, or
# any programmatic caller writing a TSV) without the caller having to
# remember to ask for every ratio by name. Mirrors the GUI's
# primat.gui.panels._RATIO_FORMAT keys, which is where this set was
# originally curated; some entries (Li6oLi7, YCNO) only exist for networks
# that track Li6/CNO and are silently dropped when unavailable -- see
# _default_mc_quantities and _c_mc below.
_DEFAULT_MC_OBSERVABLES = ("Neff", "YPBBN", "YPCMB", "DoH", "He3oH", "He3oHe4",
                           "Li7oH", "Li6oLi7", "YCNO")

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
# The C extension's cpr_config_init_defaults() takes the data folder itself
# (containing NEVO/, weak/, plasma/, nuclear/, csv/), not its parent.
_C_DATA_DIR = os.path.join(_PACKAGE_DIR, "data")

_c_ext: Any = None
try:
    from . import _primat_c as _c_ext
    HAS_C_BACKEND = True
except ImportError:
    HAS_C_BACKEND = False


def _python_solve(params: dict[str, Any] | None, extra_rho: list | None,
                   custom_network: dict[str, Any] | None, background,
                   progress: bool = True) -> dict[str, Any]:
    """Run the pure-Python backend and return PRIMAT.solve()'s result dict."""
    from .main import PRIMAT
    return PRIMAT(params=params, extra_rho=extra_rho,
                  custom_network=custom_network, background=background).solve(progress=progress)


# Number of log-spaced Tg nodes used to tabulate extra_rho for the C backend
# (see _tabulate_extra_rho). Dense enough that a cubic spline over log10(Tg)
# reproduces any smooth rho(Tg) to well below the cross-backend tolerance;
# cheap since it is one array evaluation of the (already fast) callables.
_EXTRA_RHO_GRID_NPTS = 4000


def _tabulate_extra_rho(extra_rho: list, cfg) -> tuple[list[float], list[float]]:
    """Evaluate the *sum* of the ``extra_rho`` callables on a dense log-spaced
    Tg grid, for handoff to the C backend (O-8's tabulated interface -- see
    ``primat-c/include/config.h``'s ``extra_rho_*`` fields and
    ``primat/_primat_c/_wrapper.c``).

    Python's ``extra_rho`` is a list of ``rho(Tg) -> MeV^4`` callables summed
    inside ``StandardBackground.Hubble``; a live callable cannot cross the C
    ABI, so instead we sample the summed contribution once here and let the C
    background spline it. The grid spans the full temperature range the C
    background's Friedmann ODE queries (``[T_end_MeV, T_start_cosmo_MeV]``)
    with a generous half-decade margin each side, so the C-side cubic spline
    never has to *extrapolate* over the physical range -- it only interpolates,
    where a dense log grid is essentially exact for smooth ``rho(Tg)``.

    Args:
        extra_rho: list of callables ``Tg[MeV] -> rho[MeV^4]``.
        cfg: the :class:`primat.config.PRIMATConfig` for this run (read for
            ``T_end_MeV``/``T_start_cosmo_MeV`` to size the grid).

    Returns:
        ``(T_list, val_list)``: two equal-length lists of floats -- the Tg
        nodes [MeV] (strictly increasing) and the summed extra rho [MeV^4]
        at each node.

    Example:
        >>> T, v = _tabulate_extra_rho([lambda Tg: 1.0e-3], cfg)
        >>> all(abs(x - 1.0e-3) < 1e-15 for x in v)   # constant -> flat table
        True
    """
    import numpy as np
    # Half-decade (factor ~3.16) margin below T_end and above T_start so the
    # spline interpolates -- never extrapolates -- across the queried range.
    T_lo = cfg.T_end_MeV / (10.0 ** 0.5)
    T_hi = cfg.T_start_cosmo_MeV * (10.0 ** 0.5)
    T_grid = np.logspace(np.log10(T_lo), np.log10(T_hi), _EXTRA_RHO_GRID_NPTS)
    total = np.zeros_like(T_grid)
    for fn in extra_rho:
        # Each callable may be scalar-only; evaluate element-wise to be safe
        # (matches how StandardBackground calls them one Tg at a time).
        total += np.array([float(fn(T)) for T in T_grid])
    return T_grid.tolist(), total.tolist()


def run_bbn(params: dict[str, Any] | None = None, force_backend: str | None = None,
            extra_rho: list | None = None, custom_network: dict[str, Any] | None = None,
            background=None, log_backend: bool = False,
            progress: bool | None = None) -> dict[str, Any]:
    """Run one BBN computation, dispatching to the C or Python backend.

    This mirrors ``PRIMAT(params=params, ...).solve()``'s result dict (same
    keys: ``YPBBN``, ``DoH``, ``Neff``, ... -- see ``primat.main.PRIMAT.solve``
    and ``tests/test_backend_parity.py``), so callers can switch backends
    transparently.

    Args:
        params: dict, optional. Same ``PRIMATConfig`` overrides accepted by
            ``PRIMAT(params=...)``.
        force_backend: ``{None, "auto", "c", "python"}``. ``None``/``"auto"``
            (default) picks the C extension when it is available and the
            request has no C-unsupported feature (see module docstring),
            otherwise the Python backend. ``"c"``/``"python"`` force that
            backend, raising ``RuntimeError``/``ValueError`` respectively if
            the C backend is unavailable or the request uses a C-unsupported
            feature.
        extra_rho, custom_network, background: forwarded to ``PRIMAT.__init__``
            verbatim. ``extra_rho``/``background`` are Python-only (see module
            docstring), so any non-``None`` value forces the Python backend
            regardless of ``force_backend`` (except ``force_backend="c"``,
            which raises instead). ``custom_network`` is supported on both
            backends and never forces a fallback.
        log_backend: bool, default False. Print which backend actually ran
            and why (module docstring); also triggered by setting the
            ``PRIMAT_BACKEND_LOG`` environment variable.
        progress: bool, optional. ``None`` (default) defers to
            ``params['show_progress']`` (``DEFAULT_PARAMS`` default ``True``);
            pass an explicit ``True``/``False`` to override it for this call.
            Controls the compact ``[primat]  HT.  MT.  LT.  done.`` stderr
            phase markers on both backends (suppressed when ``verbose=True``).

    Returns:
        dict: the BBN result dict (``YPBBN``, ``DoH``, ``Neff``, ..., plus a
        ``Y_final`` sub-dict of every tracked nuclide's final mass fraction).

    Example:
        >>> run_bbn({"network": "small"})["YPBBN"]
        0.24700...
        >>> run_bbn({"network": "small"}, force_backend="python")["YPBBN"]
        0.24699...
    """
    if force_backend not in (None, "auto", "c", "python"):
        raise ValueError(f"force_backend must be one of None/'auto'/'c'/'python', "
                          f"got {force_backend!r}")

    params = params or {}

    # Validate params the same way regardless of backend (PRIMATConfig's
    # __init__ does all the checking -- e.g. an unknown --network name --
    # so a bad request raises the same ValueError whether or not the C
    # backend ends up being used; the resulting cfg itself is discarded for
    # the "c" path, which re-derives its own CPRConfig from params instead).
    from .config import PRIMATConfig
    cfg = PRIMATConfig(params)
    if progress is None:
        progress = cfg.show_progress

    # background= (a custom Background object) is an inherently-Python
    # extension point with no C-side equivalent (O-8 priority 3): it forces
    # the Python backend. extra_rho and decay_era are now BOTH supported on
    # the C backend -- extra_rho via the tabulated (Tg[], rho[]) handoff below
    # (_tabulate_extra_rho + the C spline), decay_era via cprimat_run's own
    # DT-era matrix-exponential propagation (mirrors _integrate_decay_era) --
    # so neither is a python_only_feature any more.
    python_only_feature = background is not None

    def _c_extra_rho_kwargs() -> dict[str, Any]:
        """Tabulate extra_rho for the C wrapper's extra_rho_T/extra_rho_val
        kwargs (empty when no extra_rho was given, so the call shape is
        unchanged for the common case)."""
        if extra_rho is None:
            return {}
        T_list, val_list = _tabulate_extra_rho(extra_rho, cfg)
        return {"extra_rho_T": T_list, "extra_rho_val": val_list}

    if force_backend == "python":
        _log_backend("run_bbn", "Python", "force_backend='python'", log_backend)
        return _python_solve(params, extra_rho, custom_network, background, progress=progress)

    if force_backend == "c":
        if not HAS_C_BACKEND:
            raise RuntimeError(
                "force_backend='c' requested but primat._primat_c is not "
                "available (the C extension failed to build or was not "
                "compiled -- see setup.py)."
            )
        if python_only_feature:
            raise ValueError(
                "force_backend='c' is incompatible with background= "
                "(a custom Background object is a Python-only extension point, "
                "no C-side equivalent)."
            )
        _log_backend("run_bbn", "C", "force_backend='c'", log_backend)
        _data_dir = (params or {}).get("data_dir") or _C_DATA_DIR
        return _assemble_c_result(_c_ext.run_bbn(params, _data_dir, custom_network,
                                                   show_progress=int(progress),
                                                   **_c_extra_rho_kwargs()))

    # force_backend in (None, "auto"): use the C backend opportunistically,
    # falling back to Python for anything it cannot express.
    if HAS_C_BACKEND and not python_only_feature:
        _log_backend("run_bbn", "C", "auto, no C-unsupported feature requested", log_backend)
        _data_dir = (params or {}).get("data_dir") or _C_DATA_DIR
        return _assemble_c_result(_c_ext.run_bbn(params, _data_dir, custom_network,
                                                   show_progress=int(progress),
                                                   **_c_extra_rho_kwargs()))
    reason = ("auto fallback: background= requested"
              if python_only_feature else "auto fallback: C extension unavailable")
    _log_backend("run_bbn", "Python", reason, log_backend)
    return _python_solve(params, extra_rho, custom_network, background, progress=progress)


def _assemble_c_result(result: dict[str, Any]) -> dict[str, Any]:
    """Replaces the C extension's plain-list ``"evolution"`` dict (see
    ``primat/_primat_c/_wrapper.c``'s ``evolution_to_dict``) with the same
    :class:`primat.evolution.EvolutionResult` the Python backend attaches
    under ``result["evolution"]`` -- so callers can switch backends
    transparently. No-op if ``output_time_evolution``
    wasn't requested (no ``"evolution"`` key at all)."""
    evo = result.get("evolution")
    if evo is None:
        return result
    import numpy as np
    from .evolution import EvolutionResult
    # Optional per-reaction rate columns: the C wrapper (evolution_to_dict)
    # only emits a "rates" key when it populated them
    # (output_rates_time_evolution on); its absence maps to
    # EvolutionResult.rates = None, matching the Python backend.
    rates = evo.get("rates")
    result["evolution"] = EvolutionResult(
        t=np.asarray(evo["t"]), a=np.asarray(evo["a"]), T_gamma=np.asarray(evo["T_gamma"]),
        T_nu={"e": np.asarray(evo["T_nue"]), "mu": np.asarray(evo["T_numu"]),
              "tau": np.asarray(evo["T_nutau"])},
        Y={name: np.asarray(arr) for name, arr in evo["Y"].items()},
        rates=({name: np.asarray(arr) for name, arr in rates.items()}
               if rates else None),
    )
    return result


def _default_mc_quantities(params: dict[str, Any] | None) -> list[str]:
    """Every tracked nuclide's final-Y name plus the standard observables.

    Resolved from one ordinary :func:`run_bbn` call (cheap relative to an
    ``num_mc``-sample MC run) rather than re-deriving the network's nuclide
    list from scratch, so this always matches exactly what the chosen
    ``network``/``amax``/``custom_network`` would track -- no duplicated
    network-introspection logic between here and ``NuclearNetwork``/
    ``cpr_nuclear_network``.
    """
    central = run_bbn(params)
    names = list(central["Y_final"].keys())
    names += [q for q in _DEFAULT_MC_OBSERVABLES if q in central]
    return names


def _assemble_c_mc_result(raw: dict[str, Any], quantities: list[str], seed: int | None,
                           params: dict[str, Any], custom_network: dict[str, Any] | None) -> "MCResult":
    """Converts the C extension's ``run_mc`` dict (``{name: {central, mean,
    std, values}}``, see ``_wrapper.c``) into the same
    :class:`primat.main.MCResult` :func:`primat.main.mc_uncertainty` returns,
    so callers can switch backends transparently. Mean/std are recomputed
    from ``values`` via :class:`primat.main.MCQuantityResult` (rather than
    trusting the C side's own mean/std fields) so both backends' MCResult
    objects are built by the exact same code, with only the sample source
    differing. ``backend="c"`` is recorded so a later ``prev=`` reuse-guard
    (here or in ``mc_uncertainty``) never mixes this result's xoshiro256**
    samples with the Python backend's NumPy samples.
    """
    from .main import MCQuantityResult, MCResult
    # Build MCResult from all keys in raw (includes both quantities and nuclides)
    data = {q: MCQuantityResult(raw[q]["central"], raw[q]["values"]) for q in raw}
    return MCResult(data, seed=seed, params=params, custom_network=custom_network, backend="c")


def _c_prev_reuse(prev: "MCResult | None", seed: int | None, quantities: list[str],
                   base_params: dict[str, Any], custom_network: dict[str, Any] | None) -> bool:
    """The C-path counterpart of ``mc_uncertainty``'s internal ``reuse``
    check (``primat/main.py``): delegates to the shared
    :func:`primat.main.mc_prev_is_reusable` guard with ``backend='c'``, so
    the two call sites (this one and ``mc_uncertainty``'s own) can never
    drift apart -- see that function's docstring for the exact conditions,
    in particular why a Python-origin ``prev`` must never be fed to the C
    side as if its samples were resumable (different, non-interchangeable
    RNG streams -- see this module's docstring).
    """
    from .main import mc_prev_is_reusable
    return mc_prev_is_reusable(prev, seed, quantities, base_params,
                                custom_network, backend='c')


def run_mc(num_mc: int, quantities: str | list[str] | None = None,
           params: dict[str, Any] | None = None, force_backend: str | None = None,
           seed: int | None = 0, n_jobs: int = -1, prev: "MCResult | None" = None,
           custom_network: dict[str, Any] | None = None, log_backend: bool = False,
           progress: bool | None = None) -> "MCResult":
    """Run an MC nuclear-rate/tau_n uncertainty propagation, dispatching to
    the C or Python backend (the MC counterpart of :func:`run_bbn`).

    This mirrors :func:`primat.main.mc_uncertainty`'s return value (an
    :class:`primat.main.MCResult`, indexed by quantity name -- same
    ``.central``/``.mean``/``.std``/``.values`` per quantity), so callers can
    switch backends transparently; see this module's docstring for the
    RNG caveat (C samples are statistically, not bit-for-bit, comparable to
    Python's).

    Args:
        num_mc: int. Number of MC samples.
        quantities: str or list of str, optional. A result-dict key
            (``'YPBBN'``, ``'DoH'``, ...) or nuclide name, or a list of
            either. ``None`` (default) uses every tracked nuclide's final Y
            plus the full ``_DEFAULT_MC_OBSERVABLES`` set (see
            :func:`_default_mc_quantities`). Regardless of what is passed
            here, the returned ``MCResult`` *always* additionally contains
            every tracked nuclide and every ``_DEFAULT_MC_OBSERVABLES`` entry
            this network/custom_network actually produces -- at no extra
            solving cost, since each MC sample already runs a full solve --
            so a TSV dump (:func:`dump_mc_samples`) is always complete even
            when ``quantities`` only asked for one or two values for display.
        params, seed, n_jobs: forwarded verbatim; see
            ``primat.main.mc_uncertainty``'s docstring.
        force_backend: ``{None, "auto", "c", "python"}``, same semantics as
            :func:`run_bbn`.
        prev: supported on both backends (see module docstring); a
            previously computed :class:`primat.main.MCResult` to *extend*
            rather than recompute from scratch. Reused only when it is
            sample-compatible (same seed/quantities/params/custom_network)
            *and* came from the same backend that will compute this call
            (``prev.backend``); otherwise silently ignored, mirroring
            ``mc_uncertainty``'s own fallback. Never forces a backend switch
            or raises.
        custom_network: supported on both backends (forwarded to
            ``cpr_mc_uncertainty``'s ``CPRCustomNetwork*``); never forces a
            fallback.
        log_backend: bool, default False. Print which backend actually ran
            and why (module docstring); also triggered by setting the
            ``PRIMAT_BACKEND_LOG`` environment variable.
        progress: bool, optional. ``None`` (default) defers to
            ``params['show_progress']`` (``DEFAULT_PARAMS`` default ``True``);
            pass an explicit ``True``/``False`` to override it for this call.
            Controls the ``[MC] Running N samples...`` banner and the
            ``N/total (XX%)`` counter on both backends.

    Returns:
        primat.main.MCResult

    Example:
        >>> run_mc(50, ['YPBBN', 'DoH'], params={'network': 'small'})['YPBBN'].std
        >>> run_mc(50, force_backend='python')['DoH'].mean
    """
    if force_backend not in (None, "auto", "c", "python"):
        raise ValueError(f"force_backend must be one of None/'auto'/'c'/'python', "
                          f"got {force_backend!r}")

    params = params or {}
    from .config import PRIMATConfig
    cfg = PRIMATConfig(params)  # validate params the same way regardless of backend
    if progress is None:
        progress = cfg.show_progress

    # Don't resolve quantities=None eagerly with a probe run_bbn -- each
    # backend resolves it from its own central solve (mc_uncertainty for
    # Python, the discovery run_bbn in _c_mc for C), so only one solve is
    # needed rather than two.
    quantities = [quantities] if isinstance(quantities, str) else quantities

    # mc_uncertainty() applies these same defaults to `base_params` before
    # storing it on the MCResult it returns (for its own reuse-guard) -- so
    # the C path's reuse-guard comparison below must use the identically
    # defaulted dict, or a Python-origin params dict would never compare
    # equal to itself.
    base_params = dict(params)
    base_params.setdefault('verbose', False)
    base_params.setdefault('debug', False)

    def _python_mc():
        from .main import mc_uncertainty
        return mc_uncertainty(num_mc, quantities, params=params, n_jobs=n_jobs,
                               seed=seed, prev=prev, custom_network=custom_network,
                               progress=progress)

    def _c_mc():
        # One ordinary (non-MC) solve to learn the nuclide list and which
        # optional derived observables (Li6oLi7/YCNO/Neff/...) this
        # network/config actually produces -- same role as
        # _default_mc_quantities's `central` for the Python backend's
        # quantities=None path, but needed here unconditionally since the
        # default-observable merge below always applies, not just when
        # quantities was omitted.
        # Probe solve: discovers which nuclides and optional observables
        # (Li6oLi7, YCNO, …) this network/config actually produces, so we
        # can build a complete quantities_with_nuclides list to pass to
        # cpr_mc_uncertainty.  progress=False suppresses phase markers here --
        # cpr_mc_uncertainty will print them for its own central solve.
        central = run_bbn(params, custom_network=custom_network,
                           force_backend="c", progress=False)
        all_nuclides = list(central["Y_final"].keys())
        # Always merge in the standard observables (filtered to those this
        # network actually has) and every nuclide, regardless of what the
        # caller explicitly requested -- so the returned MCResult is always
        # complete enough to dump to disk (see _DEFAULT_MC_OBSERVABLES).
        qty_set = set(quantities) if quantities is not None else set()
        extra_observables = [q for q in _DEFAULT_MC_OBSERVABLES
                              if q not in qty_set and q in central]
        quantities_plus_observables = (list(quantities) if quantities is not None
                                       else []) + extra_observables
        full_set = set(quantities_plus_observables)
        quantities_with_nuclides = (quantities_plus_observables
                                     + [nm for nm in all_nuclides if nm not in full_set])

        if _c_prev_reuse(prev, seed, quantities_with_nuclides, base_params, custom_network):
            n_prev = (min(len(prev[quantities_with_nuclides[0]].values), num_mc)
                      if quantities_with_nuclides else 0)
            prev_centrals = [prev[q].central for q in quantities_with_nuclides]
            prev_values = [list(prev[q].values[:n_prev]) for q in quantities_with_nuclides]
        else:
            prev_centrals = None
            prev_values = None
        _data_dir = (params or {}).get("data_dir") or _C_DATA_DIR
        raw = _c_ext.run_mc(params, _data_dir, num_mc, quantities_with_nuclides, seed, n_jobs,
                             custom_network, prev_centrals, prev_values,
                             progress=int(progress))
        return _assemble_c_mc_result(raw, quantities_with_nuclides, seed, base_params, custom_network)

    if force_backend == "python":
        _log_backend("run_mc", "Python", "force_backend='python'", log_backend)
        return _python_mc()

    if force_backend == "c":
        if not HAS_C_BACKEND:
            raise RuntimeError(
                "force_backend='c' requested but primat._primat_c is not "
                "available (the C extension failed to build or was not "
                "compiled -- see setup.py)."
            )
        _log_backend("run_mc", "C", "force_backend='c'", log_backend)
        return _c_mc()

    # force_backend in (None, "auto"): use the C backend opportunistically,
    # falling back to Python for anything it cannot express.
    if HAS_C_BACKEND:
        _log_backend("run_mc", "C", "auto, C extension available", log_backend)
        return _c_mc()
    _log_backend("run_mc", "Python", "auto fallback: C extension unavailable", log_backend)
    return _python_mc()


def dump_mc_samples(mc: "MCResult") -> str:
    """Serialise an :class:`primat.main.MCResult` to TSV text: one column per
    quantity (header = quantity names, in their original order), one row per
    MC sample -- the on-disk "common language" for MC results shared by both
    backends (CLAUDE.md's backend-parity mandate), and the same shape
    written to ``<output_mc_file_prefix>_samples.tsv`` when
    ``output_mc_samples=True``.

    Args:
        mc: primat.main.MCResult.

    Returns:
        str: TSV text, with a trailing newline.
    """
    names = mc.quantity_names()
    samples = mc.samples_array()
    lines = ["\t".join(names)]
    lines += ["\t".join(f"{v:.10e}" for v in row) for row in samples]
    return "\n".join(lines) + "\n"


def _mc_num_and_seed(mc: "MCResult") -> tuple[int, str]:
    """Helper for the covariance/correlation writers: return ``(N, seed_str)``
    where ``N`` is the MC sample count (rows of :meth:`~primat.main.MCResult.samples_array`)
    and ``seed_str`` is the base seed rendered for the file header (the integer,
    or ``"None"`` for a seedless result).  Kept in one place so both matrix
    headers -- and their byte-identical C-side counterparts (``primat-c``'s
    ``mc.c``, per CLAUDE.md's verbose/output-parity mandate) -- agree on the
    wording.
    """
    n = mc.samples_array().shape[0]
    seed = mc.seed
    return n, ("None" if seed is None else str(seed))


def dump_mc_covariance(mc: "MCResult") -> str:
    """Serialise an :class:`primat.main.MCResult`'s full sample **covariance**
    matrix (``mc.cov()``; ddof=1, all MC quantities in ``quantity_names``
    order) to the two-header-line TSV written to
    ``<output_mc_file_prefix>_covariance.tsv`` when ``output_mc_covariance=True``.

    The joint (off-diagonal) covariance -- e.g. between ``YPBBN`` and ``DoH``,
    which are driven by the same MC samples -- is exactly what a user needs to
    build a multi-observable likelihood; the per-observable variances live on
    the diagonal (``C[i, i] == mc[q_i].std**2``).

    File layout (author spec)::

        # Covariance matrix of the N=100 primat MC samples (seed=0): ...
        quantity	Neff	YPBBN	...
        Neff	<C[0,0]>	<C[0,1]>	...
        ...

    i.e. line 1 is a single ``#`` comment naming the file (with N, seed and the
    ddof=1 estimator convention); line 2 is the tab-separated quantity names
    labelling both columns and rows; then one row per quantity, its name first.

    Args:
        mc: primat.main.MCResult.

    Returns:
        str: TSV text, with a trailing newline.
    """
    import numpy as np
    names = mc.quantity_names()
    C = np.atleast_2d(mc.cov())
    n, seed = _mc_num_and_seed(mc)
    header = (f"# Covariance matrix of the N={n} primat MC samples "
              f"(seed={seed}): C[i,j] = sample covariance (ddof=1) of "
              f"quantities i and j.")
    lines = [header, "quantity\t" + "\t".join(names)]
    for i, nm in enumerate(names):
        lines.append(nm + "\t" + "\t".join(f"{C[i, j]:.10e}"
                                           for j in range(len(names))))
    return "\n".join(lines) + "\n"


def dump_mc_correlation(mc: "MCResult") -> str:
    """Serialise an :class:`primat.main.MCResult`'s full sample **correlation**
    matrix (``mc.corr()``; unit diagonal, ddof=1) to the two-header-line TSV
    written to ``<output_mc_file_prefix>_correlation.tsv`` when
    ``output_mc_correlation=True``.

    Same layout as :func:`dump_mc_covariance` (line 1 = a ``#`` comment; line 2
    = the quantity names; then one labelled row per quantity), with its own
    header wording and a unit diagonal.  A quantity that was identical in every
    sample (zero variance) has NaN off-diagonal entries -- see
    :meth:`primat.main.MCResult.corr`.

    Args:
        mc: primat.main.MCResult.

    Returns:
        str: TSV text, with a trailing newline.
    """
    import numpy as np
    names = mc.quantity_names()
    R = np.atleast_2d(mc.corr())
    n, seed = _mc_num_and_seed(mc)
    header = (f"# Correlation matrix of the N={n} primat MC samples "
              f"(seed={seed}): R[i,j] = Pearson correlation (ddof=1) of "
              f"quantities i and j; unit diagonal.")
    lines = [header, "quantity\t" + "\t".join(names)]
    for i, nm in enumerate(names):
        lines.append(nm + "\t" + "\t".join(f"{R[i, j]:.10e}"
                                           for j in range(len(names))))
    return "\n".join(lines) + "\n"


def dump_final_with_sigma(names: list[str], Y: dict[str, float],
                           sigma: dict[str, float] | None = None,
                           num_mc: int | None = None) -> str:
    """Render the ``output_final.dat``-format final-abundances text.

    Two columns (``nuclide  Y``) when ``sigma`` is ``None`` -- identical to
    the plain single-run format written by
    ``NuclearNetwork._write_final_result``. Three columns (``nuclide  Y
    sigma_N<num_mc>``) when an MC ``sigma`` dict is supplied, so the sample
    count backing the uncertainty estimate is recorded directly in the
    header rather than only in the (separate) MC-samples file.

    The header row uses the same column widths as the data rows so the
    column names sit directly above their respective values (no ``#``
    comment prefix that would shift the label two characters to the right).

    Args:
        names: list of str. Nuclide names, in the order to write them.
        Y: dict, name -> final mass-fraction abundance.
        sigma: dict, name -> 1-sigma MC uncertainty on ``Y[name]``, optional.
        num_mc: int, required when ``sigma`` is given -- the MC sample count,
            recorded in the header (e.g. ``sigma_N50``).

    Returns:
        str: the file text, with a trailing newline.
    """
    if sigma is None:
        # Header width (14) matches nuclide-name column in data rows so
        # "nuclide" sits directly above the names and "Y" above the values.
        lines = [f"{'nuclide':<14}Y"]
        lines += [f"{nm:<14}{Y[nm]:.6e}" for nm in names]
    else:
        if num_mc is None:
            raise ValueError("num_mc is required when sigma is given")
        lines = [f"{'nuclide':<14}{'Y':<14}sigma_N{num_mc}"]
        lines += [f"{nm:<14}{Y[nm]:<14.6e}{sigma[nm]:.6e}" for nm in names]
    return "\n".join(lines) + "\n"
