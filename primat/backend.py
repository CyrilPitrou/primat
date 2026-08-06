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

The former ``extra_rho`` and ``decay_era`` gaps are now *closed*:

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
``primat/_primat_c_src/_wrapper.c`` parses the same dict shape
(``UpdateNuclearRates``/``kept_to_custom_network``, see
``primat/network_data.py``/``primat/gui/custom_rates.py``) into one. It is no
longer part of ``python_only_feature`` below.

``output_time_evolution=True`` *is* supported on both backends: the C
extension's ``cprimat_run`` populates ``CPRResults``'s
``evol_*`` in-memory arrays (``primat-c/include/api.h``) and
``primat/_primat_c_src/_wrapper.c`` hands them back as an ``"evolution"`` dict
key (plain Python lists, no numpy C-API dependency in the extension); this
module assembles the same :class:`primat.evolution.EvolutionResult` shape
the Python backend produces, with no disk I/O on either backend's part.

``data_dir``/``user_nuclear_dir`` (see ``docs/howto/data-overlays.md``)
*are* supported on both backends: ``data_dir`` fully
replaces the shipped data tree; ``user_nuclear_dir`` is an additive overlay
for nuclear networks and rate tables.  They are ordinary ``params`` dict keys
applied generically via ``cpr_config_set_by_name`` on the C side, so no
special-casing is needed here — except that ``data_dir`` must *also* be
forwarded as the ``data_dir`` positional argument to ``_c_ext.run_bbn``/
``_c_ext.run_mc``, because the C extension's ``cpr_config_init_defaults``
loads ``csv/nuclides.csv`` from that argument before any ``params`` key is
applied. :func:`_c_data_dir` takes it from the *validated* config rather than
from the raw dict, so a ``~``-prefixed path reaches C already expanded (the
Python side expands it via ``config._PATH_PARAMS``).

:func:`run_mc` is the MC counterpart of :func:`run_bbn`: it dispatches between
``primat._primat_c``'s ``run_mc`` (wrapping ``primat-c/src/mc.c``'s threaded
``cpr_mc_uncertainty``) and ``primat.main.mc_uncertainty`` (joblib), returning
the same :class:`primat.main.MCResult` shape either way -- the "common
language" the two backends share for MC results. The C path uses a
pthread/xoshiro256** RNG, *not* NumPy's
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

import numbers
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
# that track Li6/CNO and are silently dropped when unavailable -- each
# backend filters them against its own central solve (mc_uncertainty for
# Python, _c_mc's probe run_bbn below for C).
_DEFAULT_MC_OBSERVABLES = ("Neff", "YPBBN", "YPCMB", "He4oH", "DoH", "He3oH",
                           "He3oHe4", "Li7oH", "Li6oLi7", "YCNO")

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
# The C extension's cpr_config_init_defaults() takes the data folder itself
# (containing NEVO/, weak/, plasma/, nuclear/, csv/), not its parent.
_C_DATA_DIR = os.path.join(_PACKAGE_DIR, "data")

_c_ext: Any = None
try:
    # `_primat_c` is the compiled C extension (built by setup.py from
    # primat/_primat_c_src/_wrapper.c); it is invisible to static analysis, so
    # mypy sees neither the `primat._primat_c` attribute nor that this rebinds
    # the `_c_ext: Any` declared just above -- both are intentional at runtime.
    from . import _primat_c as _c_ext  # type: ignore[attr-defined,no-redef]
    HAS_C_BACKEND = True
except ImportError:
    HAS_C_BACKEND = False

# Known limitation (accepted, not a bug to chase): on Windows an *editable*
# install (`pip install -e .`) does not build/expose the compiled `.pyd`, so
# HAS_C_BACKEND is False there and every run transparently uses the pure-Python
# backend. A normal *wheel* install on Windows DOES ship a working extension
# (see the green windows-latest leg of .github/workflows/wheels.yml), so end
# users are unaffected -- only Windows contributors developing from a source
# checkout. macOS/Linux build the extension in-place in editable installs as
# usual. Not worth pursuing, since research use of this project is macOS/Linux.


def _python_solve(params: dict[str, Any] | None, extra_rho: list | None,
                   custom_network: dict[str, Any] | None, background,
                   progress: bool = True) -> dict[str, Any]:
    """Run the pure-Python backend and return PRIMAT.solve()'s result dict.

    Backend-parity note: the C backend's ``run_bbn`` attaches a ``"Y_final"``
    sub-dict of every tracked nuclide's final mass fraction
    (``primat/_primat_c_src/_wrapper.c``). ``PRIMAT.solve()`` itself does not
    include it -- an in-process Python caller would query
    ``inst.get_quantity(...)`` / ``inst.nuclear.Y_final`` instead -- but
    :func:`run_bbn` returns only the result dict, with no instance to query,
    so we must mirror the C backend and attach ``"Y_final"`` here or callers
    (and ``tests/test_backend_parity.py``) see divergent result-dict keys
    across backends.
    """
    from .main import PRIMAT
    inst = PRIMAT(params=params, extra_rho=extra_rho,
                  custom_network=custom_network, background=background)
    result = inst.solve(progress=progress)
    result["Y_final"] = dict(inst.nuclear.Y_final)
    return result


# Number of log-spaced Tg nodes used to tabulate extra_rho for the C backend
# (see _tabulate_extra_rho). Dense enough that a cubic spline over log10(Tg)
# reproduces any smooth rho(Tg) to well below the cross-backend tolerance;
# cheap since it is one array evaluation of the (already fast) callables.
_EXTRA_RHO_GRID_NPTS = 4000


def _tabulate_extra_rho(extra_rho: list, cfg) -> tuple[list[float], list[float]]:
    """Evaluate the *sum* of the ``extra_rho`` callables on a dense log-spaced
    Tg grid, for handoff to the C backend's tabulated interface -- see
    ``primat-c/include/config.h``'s ``extra_rho_*`` fields and
    ``primat/_primat_c_src/_wrapper.c``).

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


def _c_data_dir(cfg) -> str:
    """Return the data root to hand the C extension as its positional
    ``data_dir`` argument: the run's own ``data_dir`` when set, else the
    package-shipped tree.

    Taken from the *validated* :class:`primat.config.PRIMATConfig` rather than
    from the raw ``params`` dict so the value is already ``~``-expanded (see
    ``config._PATH_PARAMS``); handing C the raw ``"~/mydata"`` would have it
    look for a literal ``./~/mydata`` tree.

    Args:
        cfg: the PRIMATConfig built for this call.

    Returns:
        str: the data root for this run -- an existing directory carrying at
        least ``csv/`` and ``nuclear/`` (validated by
        ``PRIMATConfig._validate_dir_field``), normally also ``NEVO/`` and the
        regenerable ``cache_plasma_weak/`` caches.

    Example:
        >>> _c_data_dir(PRIMATConfig({}))          # doctest: +SKIP
        '/.../site-packages/primat/data'
    """
    return cfg.data_dir or _C_DATA_DIR


def _c_params(params: dict[str, Any]) -> dict[str, Any]:
    """Return ``params`` without the keys ``PRIMATConfig`` has already reported
    as unknown, for handoff to the C extension.

    ``cpr_config_set_by_name`` rejects any name it does not recognise, and the
    wrapper turns that into a ``ValueError`` -- so a plain typo
    (``{"Omegab2h": 0.022}``) *raised* on the C path while the Python path
    warned "did you mean 'Omegabh2'?" and ran with the default cosmology. That
    contradicts ``strict_params``, whose documented default (``False``) is
    "warn and ignore"; with ``strict_params=True``, ``PRIMATConfig`` has
    already raised before we get here, so filtering can never hide a strict
    error.

    Only keys unknown to *both* sides are dropped: anything in
    ``DEFAULT_PARAMS`` (or a ``p_<rxn>``/``delta_<rxn>`` variation) is still
    forwarded, so a key Python accepts but the C field table lacks stays a hard
    error rather than being silently ignored -- that asymmetry is a parity bug
    worth failing on, and is exactly how the missing ``data_dir`` case was
    caught.

    Args:
        params: the numpy-unwrapped PRIMATConfig overrides for this run.

    Returns:
        dict: ``params`` itself when nothing needs dropping (the common path),
        otherwise a filtered copy.

    Example:
        >>> _c_params({"network": "small", "Omegab2h": 0.022})
        {'network': 'small'}
    """
    from .config import DEFAULT_PARAMS
    known = [k for k in params
             if k in DEFAULT_PARAMS or k.startswith(("p_", "delta_"))]
    if len(known) == len(params):
        return params
    return {k: params[k] for k in known}


def _validate_params(params: dict[str, Any]):
    """Build the throwaway :class:`primat.config.PRIMATConfig` that validates
    ``params`` identically for both backends, returning
    ``(cfg, warnings_list)``.

    Every request is type/range/choice-checked here (an unknown ``--network``
    name, ``Omegabh2="0.022"``, a ``p_<rxn>`` typo, ...) so a bad request fails
    the same way whether or not the C extension ends up servicing it — the
    resulting ``cfg`` is discarded on the C path, which re-derives its own
    ``CPRConfig`` from the same dict.

    Warnings raised during that construction (unknown keys, unmatched
    ``p_<rxn>``/``delta_<rxn>`` names) are *captured* rather than emitted,
    because the Python backend builds a second ``PRIMATConfig`` inside
    ``PRIMAT.__init__`` and would emit each of them a second time — one typo
    reading as two problems. The caller re-emits them via
    :func:`_reemit_warnings` on the C path only, where no second construction
    happens.

    Args:
        params: dict of PRIMATConfig overrides (already numpy-unwrapped).

    Returns:
        ``(cfg, caught)``: the validated config, and the list of
        :class:`warnings.WarningMessage` recorded while building it.

    Example:
        >>> cfg, caught = _validate_params({"network": "small"})
        >>> caught
        []
    """
    import warnings
    from .config import PRIMATConfig
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cfg = PRIMATConfig(params)
    return cfg, list(caught)


def _reemit_warnings(caught: list) -> None:
    """Re-emit warnings captured by :func:`_validate_params`, preserving each
    one's category and message.

    Called on the C path only: nothing else will construct a
    ``PRIMATConfig`` there, so without this a mistyped parameter key would be
    validated and then silently swallowed. ``stacklevel=3`` attributes the
    warning to the caller's ``run_bbn``/``run_mc`` call rather than to this
    helper (``warn`` -> ``_reemit_warnings`` -> ``run_bbn`` -> user).
    """
    import warnings
    for w in caught:
        warnings.warn(w.message, w.category, stacklevel=3)


def _unwrap_numpy_params(params: dict[str, Any]) -> dict[str, Any]:
    """Return ``params`` with every numpy scalar value replaced by its plain
    Python equivalent (``np.int64(80) -> 80``), leaving everything else
    untouched.

    Why both backends need this, at the dispatch boundary rather than in
    ``PRIMATConfig``: numpy scalars arrive naturally from a parameter scan
    (``for n in np.arange(40, 200, 20): run_bbn({"sampling_nTOp_per_decade": n})``)
    or from an external driver such as the Cobaya wrapper indexing a sampled
    array, and *neither* backend used to accept all of them. ``np.float64``
    passed (it subclasses ``float``) but ``np.int64``/``np.float32``/
    ``np.bool_`` failed -- the C path with a clear
    ``TypeError: unsupported parameter value type numpy.int64`` from the
    extension wrapper, the Python path with an opaque
    ``TypeError: Object of type int64 is not JSON serializable`` raised from
    inside the weak-rate cache fingerprint. Unwrapping here fixes both at once
    and keeps their accepted input identical, which is the point: a config a
    user can run on one backend must run on the other.

    Value-preserving and therefore fingerprint-preserving: ``.item()`` returns
    the exact same number, so no cache file is invalidated (see
    ``cache_utils._json_scalar``, the defensive second layer for callers that
    bypass this function by using ``PRIMAT(params=...)`` directly).

    Args:
        params: dict of PRIMATConfig overrides, possibly holding numpy scalars.

    Returns:
        dict: a new dict when anything was unwrapped, otherwise ``params``
        itself (so the common path allocates nothing).

    Example:
        >>> import numpy as np
        >>> _unwrap_numpy_params({"amax": np.int64(8)})
        {'amax': 8}
    """
    import numpy as np   # lazy, matching _tabulate_extra_rho's convention

    if not any(isinstance(v, np.generic) for v in params.values()):
        return params
    return {k: (v.item() if isinstance(v, np.generic) else v)
            for k, v in params.items()}


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
            verbatim. ``background`` is the one Python-only extension point
            (see module docstring), so a non-``None`` value forces the Python
            backend regardless of ``force_backend`` — except
            ``force_backend="c"``, which raises instead. ``extra_rho`` and
            ``custom_network`` are supported on both backends and never force
            a fallback (``extra_rho`` crosses to C as a tabulated
            ``(Tg[], rho[])`` handoff, see :func:`_tabulate_extra_rho`).
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

    params = _unwrap_numpy_params(params or {})

    # Validate params the same way regardless of backend; construction
    # warnings are held back and re-emitted only on the C path (see
    # _validate_params / _reemit_warnings).
    cfg, caught = _validate_params(params)
    if progress is None:
        progress = cfg.show_progress

    # background= (a custom Background object) is an inherently-Python
    # extension point with no C-side equivalent: it forces
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
        _reemit_warnings(caught)
        return _assemble_c_result(_c_ext.run_bbn(_c_params(params), _c_data_dir(cfg), custom_network,
                                                   show_progress=int(progress),
                                                   **_c_extra_rho_kwargs()))

    # force_backend in (None, "auto"): use the C backend opportunistically,
    # falling back to Python for anything it cannot express.
    if HAS_C_BACKEND and not python_only_feature:
        _log_backend("run_bbn", "C", "auto, no C-unsupported feature requested", log_backend)
        _reemit_warnings(caught)
        return _assemble_c_result(_c_ext.run_bbn(_c_params(params), _c_data_dir(cfg), custom_network,
                                                   show_progress=int(progress),
                                                   **_c_extra_rho_kwargs()))
    reason = ("auto fallback: background= requested"
              if python_only_feature else "auto fallback: C extension unavailable")
    _log_backend("run_bbn", "Python", reason, log_backend)
    return _python_solve(params, extra_rho, custom_network, background, progress=progress)


def _assemble_c_result(result: dict[str, Any]) -> dict[str, Any]:
    """Replaces the C extension's plain-list ``"evolution"`` dict (see
    ``primat/_primat_c_src/_wrapper.c``'s ``evolution_to_dict``) with the same
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
        num_mc: int. Number of MC samples; must be >= 1. Note that a sigma
            needs at least two samples, so ``num_mc=1`` legitimately reports
            ``std == 0``.
        quantities: str or list of str, optional. A result-dict key
            (``'YPBBN'``, ``'DoH'``, ...) or nuclide name, or a list of
            either. ``None`` (default) uses every tracked nuclide's final Y
            plus the full ``_DEFAULT_MC_OBSERVABLES`` set, as resolved by
            whichever backend runs (``mc_uncertainty``'s central solve for
            Python, ``_c_mc``'s probe solve for C). Regardless of what is passed
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

    # num_mc must be a positive count, checked here so BOTH backends reject the
    # same input: the C sampler allocates num_mc doubles per quantity, so a
    # negative value used to underflow size_t and abort the process with a
    # "primat: out of memory (18446744073709551576 bytes)" from mc.c, while the
    # Python path silently produced an empty sample set reported as
    # "value +/- 0". (cpr_mc_uncertainty now rejects it too, for callers that
    # reach the C API without passing through here.)
    if not isinstance(num_mc, numbers.Integral) or isinstance(num_mc, bool):
        raise TypeError(f"run_mc: num_mc must be an int, got {num_mc!r} of type "
                         f"{type(num_mc).__name__}.")
    if num_mc < 1:
        raise ValueError(f"run_mc: num_mc must be >= 1, got {num_mc}. (A sigma "
                          "needs at least 2 samples; num_mc=1 reports std=0.)")

    # Unwrap numpy scalars before anything else, exactly as run_bbn does, so an
    # MC scan driven by numpy arrays behaves identically on both backends and
    # the `prev` reuse-guard compares plain-Python params dicts (a np.float64
    # would compare equal to its float, but an np.int64 key would not survive
    # the C hand-off at all -- see _unwrap_numpy_params).
    params = _unwrap_numpy_params(params or {})
    # Validate params identically for both backends; construction warnings are
    # held back and re-emitted only on the C path, exactly as in run_bbn --
    # the Python path's own mc_uncertainty builds a PRIMATConfig and a PRIMAT
    # in this same process (main.py's _mc_rate_keys / central solve), which
    # emit them there.
    cfg, caught = _validate_params(params)
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
        _reemit_warnings(caught)
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
        raw = _c_ext.run_mc(_c_params(params), _c_data_dir(cfg), num_mc, quantities_with_nuclides, seed, n_jobs,
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
    backends, and the same shape
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
    ``mc.c``) -- agree on the wording.
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
