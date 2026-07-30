# -*- coding: utf-8 -*-
"""
main.py
=======
Main class for primat.

Design
------
* ``PRIMAT.__init__(params)`` accepts an optional dict of parameters,
  builds a ``PRIMATConfig``, loads all data files (thermodynamics tables and
  nuclear rate tables), and pre-computes the thermal background.
* ``PRIMAT.solve()`` runs the full nuclear network ODE integration and
  returns the BBN predictions.
* ``PRIMAT.primat_results()`` calls ``solve()`` and returns the result dict
  (for backwards compatibility).

"""

from __future__ import annotations

import re
import sys
import time
import warnings
from typing import Any, Callable, TYPE_CHECKING, cast
import numpy as np
from numpy.typing import NDArray

__all__ = ['PRIMAT', 'mc_uncertainty']

from .config       import PRIMATConfig
from . import plasma      as primat_thermo
from .background   import Background, StandardBackground, CustomBackground
from .nuclear_network import NuclearNetwork

if TYPE_CHECKING:
    from .evolution import EvolutionResult


# Column order for abundance interpolators; names match PRIMATConfig.Nuclides keys
_NUC_NAMES_SMALL = ["n", "p", "H2", "H3", "He3", "He4", "Li7", "Be7"]
_NUC_NAMES_FULL  = ["n", "p", "H2", "H3", "He3", "He4", "Li7", "Be7",
                    "He6", "Li8", "Li6", "B8"]

# Standard derived observables unconditionally merged into every MCResult
# returned by mc_uncertainty, on top of whatever `quantity` the caller
# explicitly requested and every tracked nuclide's final Y (see
# mc_uncertainty below) -- so the result is always complete enough to dump
# to disk (primat.backend.dump_mc_samples) without the caller having to
# remember to ask for every ratio by name. Mirrors
# primat.gui.panels._RATIO_FORMAT and primat.backend._DEFAULT_MC_OBSERVABLES
# (kept in sync by hand -- see primat.backend's copy of this constant for the
# C-backend counterpart). Li6oLi7/YCNO only exist for networks tracking
# Li6/CNO and are silently dropped when unavailable, exactly like nuclides
# that a custom_network removes.
_DEFAULT_MC_OBSERVABLES = ("Neff", "YPBBN", "YPCMB", "He4oH", "DoH", "He3oH",
                           "He3oHe4", "Li7oH", "Li6oLi7", "YCNO")

_BANNER_TEMPLATE = """
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                                               ┃
┃         ░█▀█░█▀▄░▀█▀░█▄█░█▀█░▀█▀              ┃
┃         ░█▀▀░█▀▄░░█░░█░█░█▀█░░█░              ┃
┃         ░▀░░░▀░▀░▀▀▀░▀░▀░▀░▀░░▀░              ┃
┃                                               ┃
┃  Welcome to PRIMAT (python backend) v{version}    ┃
┃                                               ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
"""


def _banner() -> str:
    """Render the startup banner with the installed package version.

    Reads ``primat.__version__`` lazily (rather than caching it at import
    time in a module-level constant) so the banner always reflects the
    version actually installed (``pyproject.toml``'s single source of
    truth), even if it changes within a long-lived process (e.g. tests
    reinstalling an editable checkout).
    """
    from . import __version__
    return _BANNER_TEMPLATE.format(version=__version__)


def _options_recap(cfg: PRIMATConfig, backend: str) -> str:
    """Render the verbose-mode "options recap" block (one line per item):
    backend, network/amax, numerical_precision, the five weak-rate flags,
    tau_n, Omegabh2/eta0b, and DeltaNeff. Printed right after the banner so
    a verbose run's header fully pins down the physics configuration before
    any solving starts.
    """
    lines = [
        f"[opts-py] backend              = {backend}",
        f"[opts-py] network              = {cfg.network!r} (amax={cfg.amax})",
        f"[opts-py] numerical_precision  = {cfg.numerical_precision:.3g}",
        f"[opts-py] radiative_corrections    = {cfg.radiative_corrections}",
        f"[opts-py] finite_mass_corrections  = {cfg.finite_mass_corrections}",
        f"[opts-py] thermal_corrections      = {cfg.thermal_corrections}",
        f"[opts-py] spectral_distortions     = {cfg.spectral_distortions}",
        f"[opts-py] tau_n_normalization      = {cfg.tau_n_normalization}",
        f"[opts-py] tau_n                = {cfg.tau_n:.4g} s",
        f"[opts-py] Omegabh2             = {cfg.Omegabh2:.8g} (eta0b={cfg.eta0b:.6g})",
        f"[opts-py] DeltaNeff            = {cfg.DeltaNeff:.8g}",
    ]
    return "\n".join(lines)


class PRIMAT:
    """
    Main primat class.

    Parameters
    ----------
    params : dict, optional
        Run-time parameters overriding defaults (see ``config.DEFAULT_PARAMS``).
    extra_rho : list of callable, optional
        Extra contributions to the total energy density entering the
        Friedmann equation.  Each element is a function
        ``rho(Tg) -> MeV^4`` of the photon temperature ``Tg`` [MeV],
        summed into ``rho_tot`` by :meth:`primat.background.StandardBackground.Hubble`.
        This is the generic plug-in point for "dark sector" components; Early
        Dark Energy (``cfg.fEDE > 0``) is implemented as the first such
        plug-in (see :meth:`primat.background.StandardBackground._setup_EDE`)
        and is appended automatically -- callers do not need to include it
        here.

        Example: a constant extra radiation density of dRho [MeV^4]::

            >>> PRIMAT({"network": "small"}, extra_rho=[lambda Tg: dRho])

        Ignored (with a warning) if ``background`` is supplied: a caller
        providing a full :class:`~primat.background.Background` instance
        is expected to fold any extra energy density into it directly.
    background : primat.background.Background, optional
        A pre-built background instance to drive the nuclear network with,
        in place of the standard ``cfg.custom_background``-driven dispatch
        below. This is the seam for a fully custom expansion history (e.g. a
        non-standard cosmology that needs more than ``extra_rho`` can
        express): subclass :class:`primat.background.Background` -- whose
        docstring lists the compulsory (``T_of_t``, ``t_of_T``, ``rhoB_BBN``,
        ``weak_nTOp_frwrd``/``weak_nTOp_bkwrd``) and optional methods -- and
        pass an instance here. Since ``Background.__init__`` already takes
        ``(cfg, plasma, extra_rho)``, ``PRIMAT`` takes ``self.cfg``/``self.plasma``
        *from the supplied instance* (``background.cfg``/``background.plasma``)
        rather than building its own, so the nuclear network and the
        background always agree on which config/plasma drove them -- build
        the instance with your own ``PRIMATConfig``/``Plasma`` first, then hand
        it to ``PRIMAT``. ``None`` (default) preserves today's
        ``cfg.custom_background``-based dispatch (:class:`StandardBackground`
        or :class:`CustomBackground`) using ``params``/``extra_rho`` as usual.
        Mutually exclusive with ``params``/``extra_rho``/``cfg.custom_background``,
        which only make sense for the default dispatch; supplying
        ``background`` together with ``params`` or ``extra_rho`` emits a
        warning, and the supplied ``background`` instance wins.

        .. warning::
           **The n<->p weak-rate cache cannot see your background.** It is
           keyed on the config alone
           (:func:`weak_rates.cache._weak_rate_fingerprint`), while the rates
           themselves are integrated against the ``[Tg_vec, Tnue_vec]`` grid
           *your* instance builds -- the grid behind the T_nu(T_gamma)
           interpolant every rate integrand reads. Two Background subclasses
           with different neutrino histories but the same ``cfg`` therefore
           hash to the same ``nTOp_<hash>.txt`` file, and whichever runs first
           defines the rates the other silently reuses. If your background
           departs from the standard T_nu(T_gamma) history, build its
           ``PRIMATConfig`` with ``weak_rate_cache=False`` (skip the load) and
           ``save_nTOp=False`` (skip the write-back). PRIMAT cannot detect
           this for you, and by the time it sees your instance the rates are
           already computed inside your constructor.

        Example: drive the network with a hand-built background::

            >>> from primat import Background
            >>> from primat.config import PRIMATConfig
            >>> from primat.plasma import Plasma
            >>> cfg = PRIMATConfig({"network": "small"})
            >>> plasma = Plasma(cfg)
            >>> class MyBackground(Background):
            ...     ...  # implement T_of_t, t_of_T, rhoB_BBN, weak_nTOp_*
            >>> PRIMAT(background=MyBackground(cfg, plasma))
    custom_network : dict, optional
        GUI/scripting "Customise Reactions" override, forwarded verbatim to
        :class:`primat.network_data.UpdateNuclearRates` (see its docstring
        for the ``{"removed": [...], "replaced": {...}, "added": {...}}``
        schema). ``None`` (default) uses the standard ``cfg.network`` reaction
        list unchanged. Not a ``PRIMATConfig`` field: it carries bulk table data
        rather than a fingerprintable scalar, so it does not participate in any
        rate cache fingerprint.

        Example: drop one reaction, override another's rate table, and add a
        brand-new reaction (its stoichiometry is read from the name)::

            >>> PRIMAT({"network": "small"}, custom_network={
            ...     "removed": ["d_d__t_p"],
            ...     "replaced": {"n_p__d_g": "0.001 1.2e3\\n10.0 4.5e1\\n"},
            ...     "added": {"t_t__He4_n_n": "0.001 1.0e2\\n10.0 1.0e2\\n"},
            ... })
    """

    def __init__(self, params: dict[str, Any] | None = None,
                 extra_rho: list[Callable[[float], float]] | None = None,
                 custom_network: dict[str, Any] | None = None,
                 background: Background | None = None) -> None:

        # ------------------------------------------------------------------
        # 1. Build configuration (+ thermodynamics, step 2) -- unless a
        #    pre-built `background` was supplied, in which case its own
        #    `.cfg`/`.plasma` are reused verbatim so the nuclear network
        #    and the background never disagree on which config/plasma drove
        #    them (see the `background` parameter docstring above).
        # ------------------------------------------------------------------
        if background is not None:
            if params:
                warnings.warn(
                    "PRIMAT: params is ignored when background= is supplied "
                    "(self.cfg is taken from background.cfg instead)."
                )
            if extra_rho is not None:
                warnings.warn(
                    "PRIMAT: extra_rho is ignored when background= is supplied "
                    "(fold any extra energy density into the background "
                    "instance directly before passing it in)."
                )
            self.cfg = background.cfg
            cfg = self.cfg
            if cfg.verbose:
                print(_banner())
                print(_options_recap(cfg, backend="python"))
                for msg in cfg._init_messages:
                    print(msg)
                self._t0 = time.time()
            # A per-instance Plasma object (rather than the module-level
            # default) so that several PRIMAT instances coexisting in the same
            # process (e.g. QED_corrections=True/False comparisons, MC
            # workers) each carry their own QED/electron-thermo tables
            # without overwriting one another's state.
            self.plasma = background.plasma
        else:
            self.cfg = PRIMATConfig(params or {})
            cfg = self.cfg
            # Print the banner/options-recap *before* building Plasma(cfg)
            # below: Plasma.__init__ emits its own "[init] ..." progress
            # messages (table loading/computing) gated on cfg.verbose, and
            # those must appear after the banner, not before it.
            if cfg.verbose:
                print(_banner())
                print(_options_recap(cfg, backend="python"))
                for msg in cfg._init_messages:
                    print(msg)
                self._t0 = time.time()
            self.plasma = primat_thermo.Plasma(cfg)

        self.N = {name: NZ[0]           for name, NZ in cfg.Nuclides.items()}
        self.Z = {name: NZ[1]           for name, NZ in cfg.Nuclides.items()}
        self.A = {name: NZ[0] + NZ[1]   for name, NZ in cfg.Nuclides.items()}

        # ------------------------------------------------------------------
        # 3. Initialize nuclear network (MT/LT eras)
        # ------------------------------------------------------------------
        from .network_data import UpdateNuclearRates
        self.nucl = UpdateNuclearRates(cfg, custom_network=custom_network)

        # ------------------------------------------------------------------
        # 4. Build the cosmological background (Class 1): a<->t<->T relations,
        #    rho_B(t), n<->p weak rates, Neff/Omega_nu -- everything the
        #    nuclear network needs about the expanding Universe.
        #
        #    Three modes:
        #    * Injected (background is not None): use the caller's instance
        #      as-is -- the seam for a fully custom expansion history (see
        #      the `background` parameter docstring above).
        #    * Standard (cfg.custom_background is None): StandardBackground
        #      solves the Friedmann / entropy-conservation ODEs, loading the
        #      NEVO non-instantaneous-decoupling table when available.
        #    * Custom (cfg.custom_background is a file path): CustomBackground
        #      reads T(t)/t/a(t) directly from that file and uses the
        #      instantaneous-decoupling approximation for neutrino temperatures
        #      and n<->p weak rates.  Neff is estimated via the Friedmann
        #      equation from the supplied a(t) (see primat.background).
        #
        #    Early Dark Energy (cfg.fEDE > 0) is only supported in the
        #    standard mode (appended to extra_rho by StandardBackground).
        # ------------------------------------------------------------------
        if background is not None:
            self.background = background
        elif cfg.custom_background is not None:
            self.background = CustomBackground(cfg, self.plasma, cfg.custom_background)
        else:
            self.background = StandardBackground(cfg, self.plasma, extra_rho)

        # ------------------------------------------------------------------
        # 5. Build the nuclear network (Class 2): the HT/MT/LT ODE
        #    integration, driven by the Background's T(t)/rho_B(t)/weak
        #    rates -- see primat.nuclear_network.NuclearNetwork.
        # ------------------------------------------------------------------
        self.nuclear = NuclearNetwork(cfg, self.nucl, self.background)

        # Populated by solve() with the BBN observables dict (Neff, YPBBN,
        # DoH, ...); None until solve() has run (see _ensure_solved).
        self.results: dict[str, Any] | None = None

        if cfg.verbose:
            print(f"[init-py]  Initialisation complete in {time.time()-self._t0:.1f} s")

    # ======================================================================
    # solve(): integrate nuclear network ODEs
    # ======================================================================

    def solve(self, progress: bool | None = None) -> dict[str, Any]:
        """
        Integrate the nuclear network over the three temperature eras and
        return a dict of BBN observables.

        Delegates the ODE integration to
        :meth:`primat.nuclear_network.NuclearNetwork.solve` (Class 2),
        which is driven by ``self.background`` (Class 1, see
        :mod:`primat.background`) and populates ``self.nuclear.Y_final``,
        ``self.nuclear.abundance_names`` and ``self.nuclear.Y_of_t``.  The
        "final observables" -- light-element ratios from ``Y_final``, plus
        ``Neff``/``Omeganurel``/``OneOverOmeganunr`` from the background's
        optional neutrino-sector hooks (:meth:`Background.rho_nu_total_final`,
        :meth:`Background.N_eff`, :meth:`Background.Omeganuh2_relnu`/
        :meth:`Background.Omeganuh2_nrnu`) -- are assembled here, into
        ``self.results``, which is what ``get_quantity``/``__getitem__``/
        ``Neff()``/``YPBBN()``/... and :meth:`primat_results` read.

        The neutrino-sector keys (``Neff``, ``Omeganurel``,
        ``OneOverOmeganunr``) are only added to the dict if the background
        actually provides that information (``None`` returned from the
        corresponding hook) -- a minimal background with no neutrino-sector
        model simply omits them.

        Args:
            progress: bool, optional. ``None`` (default) defers to
                ``cfg.show_progress`` (a ``DEFAULT_PARAMS`` key, default
                ``True``). Pass an explicit ``True``/``False`` to override
                the config for this call regardless of ``show_progress``
                (e.g. ``False`` inside MC workers to avoid per-sample spam).
        """
        if progress is None:
            progress = self.cfg.show_progress
        self.nuclear.solve(progress=progress)

        # For the large network, NuclearNetwork.solve() discovers nuclides
        # beyond the small/amax-restricted set (e.g. B10, C12, ...).  Extend the
        # N/Z/A maps (built in __init__ from cfg.Nuclides) so callers (e.g.
        # the AbundanceEvolution notebook, primat.gui.panels) can look up
        # A[name]/Z[name]/N[name] for every species returned by
        # abundance_names.
        if not self.cfg.is_small:
            for s, (N, Z) in self.nucl.large_NZ.items():
                self.N[s], self.Z[s], self.A[s] = N, Z, N + Z

        cfg  = self.cfg
        finL = self.nuclear.Y_final
        Yp_f, Yd_f, Yt_f, YHe3_f, Ya_f, YLi7_f, YBe7_f = (
            finL[k] for k in ("p", "H2", "H3", "He3", "He4", "Li7", "Be7"))

        # Primordial helium mass fraction (BBN definition: Y_p = 4 Y_He4,
        # i.e. mass of He4 over total baryon mass, since Y is the per-baryon
        # abundance normalised so sum_s A_s Y_s = 1).
        YPBBN = 4. * Ya_f
        # CMB-convention helium fraction n_He/(n_He+n_H) by mass (used by CMB
        # codes' "Y_He"): convert the BBN Y_p (mass fraction) via the He4/H
        # mass ratios cfg.He4Overma/cfg.HOverma.
        YPCMB = ((cfg.He4Overma / 4.) * YPBBN
                 / ((cfg.He4Overma / 4.) * YPBBN + cfg.HOverma * (1. - YPBBN)))

        # Custom networks may not produce some nuclide (e.g. He4 stripped out
        # of the reaction set), leaving its final abundance at exactly 0;
        # guard the corresponding ratios so we return inf/nan rather than
        # raising or silently producing a ZeroDivisionError-free but
        # misleading float (numpy floats just emit a RuntimeWarning and
        # produce inf/nan, but plain Python floats raise ZeroDivisionError).
        def _ratio(num, den):
            return num / den if den != 0.0 else (float("nan") if num == 0.0 else float("inf"))

        results = {
            "YPCMB":   YPCMB,
            "YPBBN":   YPBBN,
            # He4/H by *number* (not mass): the directly observed quantity in
            # e.g. extragalactic HII-region spectroscopy, related to the mass
            # fraction by Y_P = 4 Y_He4 and He4/H = Y_He4/Y_p, i.e.
            # He4/H = (Y_P/4)/(1-Y_P) to a very good approximation.
            "He4oH":   _ratio(Ya_f, Yp_f),
            "DoH":     _ratio(Yd_f, Yp_f),
            "He3oH":   _ratio(Yt_f + YHe3_f, Yp_f),
            "He3oHe4": _ratio(Yt_f + YHe3_f, Ya_f),
            "Li7oH":   _ratio(YLi7_f + YBe7_f, Yp_f),
        }

        # Li6/Li7: observable ratio after Be7→Li7 decay (large networks).
        # Only networks that actually track Li6 (large, large+amax>=6) carry it
        # in finL; the small network never evolves Li6, so `.get(..., 0.0)`
        # yields 0 there and the ratio is simply omitted (guarded by Y>0).
        if finL.get("Li6", 0.0) > 0:
            results["Li6oLi7"] = finL["Li6"] / (YLi7_f + YBe7_f)

        # YCNO (mass fraction): total baryon mass fraction in C, N, O isotopes
        # (large net). YCNO = sum_i A_i Y_i for all C (Z=6), N (Z=7), O (Z=8)
        # isotopes.
        cno = sum(self.A[s] * finL[s]
                  for s in finL
                  if len(s) >= 2 and s[0] in "CNO" and s[1:].isdigit())
        if cno > 0:
            results["YCNO"] = cno

        # Neff = rho_nu_tot / rho_g(Tg) / ((7/8)(4/11)^(4/3)) (generic formula,
        # Background.N_eff), evaluated at the final Tg and total neutrino
        # energy density -- only if the background tracks a neutrino sector.
        final_nu = self.background.rho_nu_total_final()
        if final_nu is not None:
            Tg_f, rho_nu_tot_f = final_nu
            results["Neff"] = self.background.N_eff(Tg_f, rho_nu_tot_f)

        relnu = self.background.Omeganuh2_relnu()
        if relnu is not None:
            results["Omeganurel"] = relnu * 1e+6

        nrnu = self.background.Omeganuh2_nrnu()
        if nrnu is not None:
            results["OneOverOmeganunr"] = 1. / (nrnu * 1e-6)

        if cfg.output_background_evolution:
            self.background.write_time_evolution(cfg.output_background_file,
                                                   cfg.output_n_points)

        # Unified time-evolution result (primat.evolution.EvolutionResult),
        # populated in memory by NuclearNetwork.solve()
        # above (no disk I/O required to get it) -- only when requested via
        # output_time_evolution=True.
        if self.nuclear.evolution is not None:
            results["evolution"] = self.nuclear.evolution

        self.results = results
        return results

    # ======================================================================
    # Public API
    # ======================================================================

    @property
    def T_of_t(self) -> Callable[[Any], Any]:
        """T_γ(t) interpolator [MeV], available after initialisation."""
        return self.background.T_of_t

    @property
    def t_of_T(self) -> Callable[[Any], Any]:
        """t(T_γ) interpolator [s], available after initialisation."""
        return self.background.t_of_T

    @property
    def a_of_T(self) -> Callable[[Any], Any]:
        """Scale factor a(T_γ), available after initialisation.

        ``a`` follows the same normalisation as the internal a(T) ODE
        (:meth:`primat.background.StandardBackground._setup_background_and_cosmo`):
        entropy conservation ``a^3 * spl(T) = const`` fixed so that
        ``a * T -> T0CMB`` [MeV] as ``T -> 0``, i.e. ``a = 1`` today up to the
        small entropy-injection correction from e+e- annihilation encoded in
        ``spl(T)``.

        Example
        -------
        >>> p.a_of_T(1.0)   # scale factor at T_gamma = 1 MeV
        """
        return self.background.a_of_T

    @property
    def T_of_a(self) -> Callable[[Any], Any]:
        """T_γ(a) interpolator [MeV], available after initialisation.

        Inverse of :attr:`a_of_T`; same normalisation convention for ``a``.
        """
        return self.background.T_of_a

    @property
    def a_of_t(self) -> Callable[[Any], Any]:
        """Scale factor a(t), available after initialisation.

        Same normalisation as :attr:`a_of_T`; ``t`` is the cosmic time [s]
        used by :attr:`T_of_t`/:attr:`t_of_T`.
        """
        return self.background.a_of_t

    @property
    def t_of_a(self) -> Callable[[Any], Any]:
        """Cosmic time t(a) [s], available after initialisation.

        Inverse of :attr:`a_of_t`; same normalisation convention for ``a``.
        """
        return self.background.t_of_a

    def __getitem__(self, species: str) -> Callable[[float | NDArray[np.float64]], Any]:
        """Return Y(t) for a species name (e.g. 'H2', 'He4', 'Li7').

        Calls solve() automatically if needed.
        """
        self._ensure_solved()
        names = self.nuclear.abundance_names
        if species not in names:
            raise KeyError(
                f"Unknown species '{species}'. Available: {names}"
            )
        idx = names.index(species)
        def fn(t):
            t_arr = np.atleast_1d(np.asarray(t, dtype=float))
            vals  = self.nuclear.Y_of_t(t_arr)[:, idx]
            return float(vals[0]) if np.ndim(t) == 0 else vals
        return fn

    def _ensure_solved(self) -> None:
        if self.results is None:
            self.solve()

    def primat_results(self) -> dict[str, Any]:
        """Return the BBN result dict, running ``solve()`` first if needed."""
        self._ensure_solved()
        return cast(dict, self.results)

    @property
    def abundance_names(self) -> list[str]:
        """Tracked nuclide names, in abundance-vector order (solves if needed).

        For the large network this is the full ~59-nuclide list; accessing it
        also guarantees ``self.A``/``N``/``Z`` cover every species (handy for
        plotting ``A_i Y_i`` for all nuclides)."""
        self._ensure_solved()
        return self.nuclear.abundance_names

    @property
    def evolution(self) -> "EvolutionResult | None":
        """The unified time-evolution result (``primat.evolution.EvolutionResult``,
        ``None`` unless ``cfg.output_time_evolution=True``), solving first if
        needed. Thin alias for ``self.nuclear.evolution`` so callers that
        don't care whether they hold a live ``PRIMAT`` or a backend-agnostic
        :class:`primat.gui.run_view.GuiRun` (which has no ``.nuclear`` at
        all) can read ``run.evolution`` uniformly either way."""
        self._ensure_solved()
        return self.nuclear.evolution

    # Convenience accessors
    def Neff(self) -> float:          self._ensure_solved(); return cast(dict, self.results)["Neff"]
    def Omeganurel(self) -> float:    self._ensure_solved(); return cast(dict, self.results)["Omeganurel"]
    def Omeganunonrel(self) -> float: self._ensure_solved(); return 1. / cast(dict, self.results)["OneOverOmeganunr"]
    def YPCMB(self) -> float:         self._ensure_solved(); return cast(dict, self.results)["YPCMB"]
    def YPBBN(self) -> float:         self._ensure_solved(); return cast(dict, self.results)["YPBBN"]
    def He4oH(self) -> float:         self._ensure_solved(); return cast(dict, self.results)["He4oH"]
    def DoH(self) -> float:           self._ensure_solved(); return cast(dict, self.results)["DoH"]
    def He3oH(self) -> float:         self._ensure_solved(); return cast(dict, self.results)["He3oH"]
    def Li7oH(self) -> float:         self._ensure_solved(); return cast(dict, self.results)["Li7oH"]

    def get_quantity(self, quantity: str) -> float:
        """Return a scalar BBN quantity by name.

        Accepts any key from the result dict ('YPBBN', 'DoH', 'He3oH',
        'Li7oH', 'Neff', 'YPCMB', ...) or a nuclide name from
        cfg.Nuclides ('H2', 'He4', 'Li7', ...) for the final mass fraction Y.
        """
        self._ensure_solved()
        results = cast(dict, self.results)
        Y_final = self.nuclear.Y_final
        if quantity in results:
            return results[quantity]
        if quantity in Y_final:
            return Y_final[quantity]
        raise ValueError(
            f"Unknown quantity '{quantity}'. "
            f"Valid result keys: {list(results.keys())}. "
            f"Valid nuclide names: {list(Y_final.keys())}."
        )


# ---------------------------------------------------------------------------
# MC result classes
# ---------------------------------------------------------------------------

class MCQuantityResult:
    """MC statistics for a single BBN quantity.

    Attributes
    ----------
    central : float
        Value at nominal rates (all p_* = 0).
    mean : float
        Mean of MC samples.
    std : float
        1σ **sample** standard deviation of the MC samples, i.e. the unbiased
        estimator ``sqrt( sum (x_i - mean)^2 / (N - 1) )`` (``np.std`` with
        ``ddof=1``).  The ``ddof=1`` normalisation is deliberate and matches
        :meth:`MCResult.cov`/:meth:`MCResult.corr` so that the diagonal of the
        covariance matrix equals ``std**2`` exactly (see those methods).  For a
        single sample (``N < 2``) the sample standard deviation is undefined
        (0/0); we report ``0.0`` there rather than NaN so that a degenerate
        ``--mc 1`` run still prints a finite ``value +/- 0`` line.
    values : np.ndarray, shape (num_mc,)
        All individual MC sample values.
    """
    __slots__ = ('central', 'mean', 'std', 'values')

    central: float
    values: NDArray[np.float64]
    mean: float
    std: float

    def __init__(self, central: float, samples: NDArray[np.float64] | list[float]) -> None:
        self.central = float(central)
        self.values  = np.asarray(samples)
        self.mean    = float(np.mean(self.values))
        # Sample (ddof=1) standard deviation -- unbiased, and consistent with
        # the ddof=1 sample covariance/correlation matrices (MCResult.cov/corr)
        # so diag(cov) == std**2. ddof=1 needs >= 2 samples; fall back to 0.0
        # for a single-sample run (population std of one point is 0 anyway).
        self.std     = (float(np.std(self.values, ddof=1))
                        if self.values.size >= 2 else 0.0)

    def __repr__(self) -> str:
        return (f"MCQuantityResult(central={self.central:.6g}, "
                f"mean={self.mean:.6g}, std={self.std:.6g}, "
                f"n={len(self.values)})")


class MCResult:
    """MC results for one or more BBN quantities, indexed by name.

    Usage::

        mc = mc_uncertainty(100, ['YPBBN', 'DoH'], params=...)
        mc['YPBBN'].mean
        mc['YPBBN'].std
        mc['YPBBN'].values
        mc['DoH'].central

    Attributes
    ----------
    seed : int or None
        The base random seed used to generate the samples (sample ``i`` used
        ``seed + i``).  Stored so that :func:`mc_uncertainty` can *extend* an
        existing result with more samples without recomputing the ones it
        already has: a previous result is only reused when its ``seed``, its
        set/order of quantities, its ``params`` and its ``custom_network`` all
        match the new request (see ``prev`` there).
    params : dict or None
        The ``base_params`` used to compute this result (after the
        ``verbose``/``debug`` defaults were applied), stored purely for the
        ``prev`` reuse-guard comparison above.
    custom_network : dict or None
        The "Customise Reactions" override used to compute this result (see
        :class:`PRIMAT`'s docstring), stored for the same reuse-guard comparison.
    backend : str or None
        Which backend produced this result (``"python"`` or ``"c"``), stored
        so a ``prev`` reuse-guard never mixes sample streams from the two
        RNGs (NumPy's ``default_rng`` vs. the C side's pthread/xoshiro256**
        -- see ``primat/backend.py``'s module docstring): a result is only
        reused as ``prev`` when its ``backend`` matches the backend about to
        compute the extension, in addition to the seed/quantities/params/
        custom_network checks above.
    """
    def __init__(self, data: dict[str, MCQuantityResult], seed: int | None = None,
                 params: dict[str, Any] | None = None,
                 custom_network: dict[str, Any] | None = None,
                 backend: str | None = None) -> None:
        self._data = data   # dict: str -> MCQuantityResult
        self.seed = seed
        self.params = params
        self.custom_network = custom_network
        self.backend = backend

    def __getitem__(self, quantity: str) -> MCQuantityResult:
        return self._data[quantity]

    def quantity_names(self) -> list[str]:
        """Quantity names in their original (insertion) order.

        Used by :func:`primat.backend.dump_mc_samples` as the column order
        of the MC-samples TSV, so the header always matches the order the
        caller requested -- regardless of whether this ``MCResult`` came
        from the Python (`mc_uncertainty`) or C (`backend.run_mc`) path.
        """
        return list(self._data)

    def samples_array(self) -> NDArray[np.float64]:
        """Stack every quantity's ``values`` into one ``(num_mc, n_quantity)``
        array, columns in :meth:`quantity_names` order.

        This is the backend-agnostic "common language" for MC samples: any
        ``MCResult`` (Python or C in origin) can be serialised the same way
        via :func:`primat.backend.dump_mc_samples`.
        """
        return np.column_stack([self._data[q].values for q in self.quantity_names()])

    def _quantity_index(self, name: str) -> int:
        """Column index of quantity ``name`` in :meth:`samples_array`, raising
        a helpful ``KeyError`` (listing the available names) on a typo -- used
        by the scalar two-name forms of :meth:`cov`/:meth:`corr`.
        """
        names = self.quantity_names()
        try:
            return names.index(name)
        except ValueError:
            raise KeyError(
                f"Unknown MC quantity {name!r}; available quantities are "
                f"{names}."
            ) from None

    def cov(self, a: str | None = None, b: str | None = None) -> "NDArray[np.float64] | float":
        """Sample covariance of the MC samples.

        Two call forms:

        * ``mc.cov()`` -> the full ``(n_q, n_q)`` covariance matrix over *all*
          MC quantities, rows/columns in :meth:`quantity_names` order.  This is
          the joint nuclear-rate + tau_n uncertainty a user needs when
          constraining cosmology with several abundances at once (e.g. YP and
          D/H together): the off-diagonal ``C[YPBBN, DoH]`` is the covariance
          the two share because they are driven by the *same* MC samples.
        * ``mc.cov("YPBBN", "DoH")`` -> the single scalar covariance between the
          two named quantities (a typo raises ``KeyError`` listing the valid
          names).

        Convention: the **sample** covariance with ``ddof=1`` (dividing by
        ``N - 1``), i.e. ``numpy.cov(..., ddof=1)``.  This matches
        :attr:`MCQuantityResult.std` (also ``ddof=1``), so the diagonal of the
        full matrix equals each quantity's ``std**2`` exactly.  With a single
        sample (``N < 2``) the estimator is undefined and entries are NaN.

        Example::

            mc = mc_uncertainty(500, ['YPBBN', 'DoH'], params=...)
            C = mc.cov()                 # 2x2 (or larger; all quantities)
            c_yp_dh = mc.cov('YPBBN', 'DoH')   # scalar, == C's off-diagonal

        Returns:
            np.ndarray (full form) or float (two-name form).
        """
        samples = self.samples_array()          # (num_mc, n_q)
        nq = samples.shape[1]
        # The ddof=1 sample covariance is undefined for a single sample (N < 2):
        # return NaN directly rather than calling np.cov (which would both emit
        # a "Degrees of freedom <= 0" RuntimeWarning and divide by zero).
        if a is None and b is None:
            if samples.shape[0] < 2:
                return np.full((nq, nq), np.nan)
            # rowvar=False: columns (= quantities) are the variables.
            return np.atleast_2d(np.cov(samples, rowvar=False, ddof=1))
        if a is None or b is None:
            raise TypeError(
                "cov() takes either no names (full matrix) or exactly two "
                "quantity names (scalar covariance)."
            )
        ia = self._quantity_index(a)
        ib = self._quantity_index(b)
        if samples.shape[0] < 2:
            return float("nan")
        # 2-variable sample covariance; [0, 1] is cov(a, b) (== var when a == b).
        pair = np.cov(samples[:, [ia, ib]], rowvar=False, ddof=1)
        return float(np.atleast_2d(pair)[0, 1])

    def corr(self, a: str | None = None, b: str | None = None) -> "NDArray[np.float64] | float":
        """Sample (Pearson) correlation of the MC samples.

        Same two call forms as :meth:`cov`:

        * ``mc.corr()`` -> the full ``(n_q, n_q)`` correlation matrix (unit
          diagonal), rows/columns in :meth:`quantity_names` order.
        * ``mc.corr("YPBBN", "DoH")`` -> the single scalar correlation
          coefficient between the two named quantities.

        Derived from :meth:`cov` (hence also ``ddof=1``) as
        ``R[i, j] = C[i, j] / (std_i * std_j)``.  A quantity that is identical
        in every sample (zero variance -- e.g. a nuclide untouched by the
        varied reactions) has an undefined correlation with anything else:
        those **off-diagonal** entries are NaN (computed under
        ``np.errstate`` so no divide-by-zero warning storm is emitted), while
        the diagonal is set to 1.0 by convention for every quantity.

        Example::

            mc.corr('YPBBN', 'DoH')   # e.g. -0.51: YP and D/H anti-correlate

        Returns:
            np.ndarray (full form) or float (two-name form).
        """
        if a is None and b is None:
            C = self.cov()
            d = np.sqrt(np.diag(C))
            # Suppress the 0/0 -> NaN warnings for any zero-variance quantity;
            # NaN off-diagonal is the documented, intended sentinel there.
            with np.errstate(invalid="ignore", divide="ignore"):
                R = C / np.outer(d, d)
            # Unit diagonal by convention (a quantity correlates perfectly with
            # itself), including zero-variance quantities whose row/column is
            # otherwise all-NaN.
            np.fill_diagonal(R, 1.0)
            return R
        if a is None or b is None:
            raise TypeError(
                "corr() takes either no names (full matrix) or exactly two "
                "quantity names (scalar correlation)."
            )
        ia = self._quantity_index(a)
        ib = self._quantity_index(b)
        if ia == ib:
            return 1.0
        samples = self.samples_array()
        va, vb = samples[:, ia], samples[:, ib]
        sa = np.std(va, ddof=1) if va.size >= 2 else 0.0
        sb = np.std(vb, ddof=1) if vb.size >= 2 else 0.0
        if sa == 0.0 or sb == 0.0:
            return float("nan")
        cab = np.atleast_2d(np.cov(np.column_stack([va, vb]),
                                   rowvar=False, ddof=1))[0, 1]
        return float(cab / (sa * sb))

    def to_flat_dict(self) -> dict[str, float]:
        """Flatten to a plain ``{name: value, ...}`` dict carrying both each
        quantity's nominal (``central``) value and its MC uncertainty, so
        callers get e.g. both ``YPBBN`` and ``sigma_YPBBN`` as ordinary dict
        entries instead of having to reach into :class:`MCQuantityResult`
        objects one quantity at a time.

        Returns:
            dict: for every quantity name ``q`` in :meth:`quantity_names`,
            ``flat[q] = self[q].central`` and ``flat[f"sigma_{q}"] = self[q].std``.
        """
        flat = {}
        for q in self.quantity_names():
            flat[q] = self._data[q].central
            flat[f"sigma_{q}"] = self._data[q].std
        return flat

    def __iter__(self):
        return iter(self._data)

    def __repr__(self) -> str:
        lines = [f"MCResult({len(self._data)} quantities):"]
        for k, v in self._data.items():
            lines.append(f"  {k}: central={v.central:.6g}, "
                         f"mean={v.mean:.6g}, std={v.std:.6g}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Module-level MC worker (must be at module level for joblib pickling)
# ---------------------------------------------------------------------------

def _mc_run_batch(base_params: dict[str, Any], rate_keys: list[str], quantities: list[str],
                  seeds: list[int], custom_network: dict[str, Any] | None = None,
                  include_nuclides: bool = False) -> list[list[float]]:
    """Run a batch of MC samples in one process, reusing a single PRIMAT.

    The cosmological background and n<->p weak rates depend only on
    ``base_params`` — *not* on the nuclear-rate offsets ``p_*`` — so they are
    computed once when the instance is built and then reused for every sample
    in the batch.  Each sample only re-draws the nuclear rates and re-solves the
    nuclear network (the cheap part), which is the whole point of the speed-up.

    Drawing the rate vector from ``default_rng(seed)`` per seed makes the result
    for a given seed independent of how the seeds are batched, so the output is
    identical (and reproducible) regardless of ``n_jobs``.

    In addition to the ``len(rate_keys)`` nuclear-rate offsets, each sample
    draws one further ``standard_normal()`` from the *same* per-sample
    ``Generator`` (after the rate offsets, so the RNG stream order does not
    depend on ``len(rate_keys)``) and uses it to perturb the neutron lifetime:
    ``tau_n_sample = tau_n_central + std_tau_n * randn()``.  When
    ``cfg.tau_n_normalization`` is True this rescales the weak-rate
    normalisation ``background.NormWeakRates = 1/tau_n``
    (``StandardBackground._setup_weak_rates``) without recomputing anything:
    the stored rates are in units of 1/tau_n, so only tau_n itself changes.
    When ``cfg.tau_n_normalization`` is False, ``tau_n`` does not enter the
    normalisation and the draw is a harmless no-op (kept so the RNG stream is
    identical either way).

    ``custom_network``, when given, is forwarded verbatim to ``PRIMAT`` (see its
    docstring): removed reactions are excluded from ``rate_keys`` by the
    caller (:func:`mc_uncertainty`), and replaced reactions are varied using
    the *custom* table's error column (``UpdateNuclearRates`` builds
    ``expsigma`` from it), so a custom rate's uncertainty flows through
    automatically.

    When ``include_nuclides=True``, each result row appends all nuclide
    abundances from Y_final after the requested quantities, in the order
    returned by the first sample's Y_final.keys().
    """
    inst = PRIMAT(params=base_params, custom_network=custom_network)
    cfg  = inst.cfg
    tau_n_central = cfg.tau_n
    # NormWeakRates is a StandardBackground/CustomBackground-only attribute
    # (module docstring above), not part of the minimal Background interface
    # -- both concrete backgrounds PRIMAT can build here define it, so the
    # cast just tells the type checker what runtime already guarantees.
    background = cast(Any, inst.background)
    # NormWeakRates = 1/tau_n (cfg.tau_n_normalization=True), so the product
    # is 1.0 -- scaling by tau_n_central and dividing by tau_n_sample gives
    # 1/tau_n_sample without any extra computation.
    norm_times_tau_n = background.NormWeakRates * tau_n_central
    results = []
    nuclide_names = None
    for seed in seeds:
        rng    = np.random.default_rng(seed)
        p_vals = rng.standard_normal(len(rate_keys))
        for k, v in zip(rate_keys, p_vals):
            setattr(cfg, k, float(v))
        tau_n_sample = tau_n_central + cfg.std_tau_n * rng.standard_normal()
        if cfg.tau_n_normalization:
            cfg.tau_n = tau_n_sample
            background.NormWeakRates = norm_times_tau_n / tau_n_sample
        inst.solve(progress=False)
        row = [inst.get_quantity(q) for q in quantities]
        if include_nuclides:
            # Capture nuclide names from the first sample for consistency
            if nuclide_names is None:
                nuclide_names = list(inst.nuclear.Y_final.keys())
            row.extend([inst.nuclear.Y_final[nm] for nm in nuclide_names])
        results.append(row)
    return results


def _mc_collect_samples(base_params: dict[str, Any], rate_keys: list[str], quantities: list[str],
                         seeds: list[int], n_jobs: int,
                         custom_network: dict[str, Any] | None = None,
                         include_nuclides: bool = False,
                         progress: bool = False) -> NDArray[np.float64]:
    """Run :func:`_mc_run_batch` for a list of seeds and stack the results.

    Splits ``seeds`` into one chunk per worker so the expensive cosmological
    background + n<->p weak-rate setup (which does *not* depend on the sampled
    nuclear rates) is paid once per worker instead of once per sample, then
    returns the ``(len(seeds), len(quantities))`` array of sampled quantity
    values (or ``(len(seeds), len(quantities) + n_nuclides)`` when
    ``include_nuclides=True``).  Because every sample draws its rate vector from
    ``default_rng(seed)`` (see ``_mc_run_batch``), the row for a given seed is
    independent of how the seeds are chunked -- so callers can safely build the
    full sample set incrementally by collecting disjoint seed ranges and
    stacking them (this is what the ``prev`` reuse in :func:`mc_uncertainty`
    relies on).

    An empty ``seeds`` list returns an empty ``(0, len(quantities))`` array so
    the result can always be ``np.vstack``-ed with an existing sample block.

    ``custom_network`` is forwarded to every :func:`_mc_run_batch` call; it is
    a plain JSON-serialisable dict, so it pickles fine for joblib workers.

    When ``include_nuclides=True``, nuclide abundances are appended to each
    row after the requested quantities.

    When ``progress=True``, print a ``\\r``-updated sample count to stderr as
    each worker chunk completes.  Seeds are split into several chunks per
    worker (not just one), so chunks finish at staggered times and the
    progress counter advances steadily instead of jumping straight from 0%
    to 100% -- without requiring per-sample IPC.
    """
    if not seeds:
        return np.empty((0, len(quantities)))

    # joblib is an *optional* dependency (the ``mc``/``recommended`` extras):
    # a lean ``pip install primat`` core can still run Monte-Carlo serially
    # (``n_jobs == 1``) with no joblib present, and the C backend never needs
    # it at all.  Only genuine parallelism (``n_jobs != 1``) reaches for
    # joblib, and a missing joblib there is a friendly, actionable error
    # rather than a bare ImportError surfacing from deep in the MC machinery.
    parallel = (n_jobs != 1)
    if parallel:
        try:
            from joblib import Parallel, delayed, effective_n_jobs
        except ImportError as exc:      # pragma: no cover - depends on env
            raise ImportError(
                f"Parallel Monte-Carlo (n_jobs={n_jobs!r}) requires joblib, "
                "which is not installed. Install it with "
                "`pip install \"primat[mc]\"` (or `pip install joblib`), or "
                "pass n_jobs=1 to run the samples serially without joblib."
            ) from exc
        n_workers = effective_n_jobs(n_jobs)
    else:
        n_workers = 1

    # Use several chunks per worker (not just one) so the progress generator
    # below yields updates steadily throughout the run instead of bunching
    # them all near the end: with exactly one chunk per worker, every chunk
    # takes about the same time and they all complete near-simultaneously.
    # In serial mode ``n_workers == 1``, so the chunking still gives the
    # progress counter (below) ~10 evenly-spaced updates over the run.
    _chunks_per_worker = 10
    n_chunks = max(1, min(len(seeds), n_workers * _chunks_per_worker))
    chunks   = [list(c) for c in np.array_split(seeds, n_chunks)]
    total    = len(seeds)
    done     = 0

    if progress:
        _w = len(str(total))
        print(f"\r[MC] {done:{_w}}/{total} samples ({0:3d}%)",
              end='', file=sys.stderr, flush=True)

    if parallel:
        # return_as='generator' yields results in submission order as each
        # worker completes its chunk, letting us update the progress counter
        # without waiting for all chunks to finish before seeing any output.
        gen = Parallel(n_jobs=n_jobs, return_as='generator')(
            delayed(_mc_run_batch)(base_params, rate_keys, quantities, chunk,
                                   custom_network=custom_network,
                                   include_nuclides=include_nuclides)
            for chunk in chunks
        )
    else:
        # Serial path: map over the chunks in-process, with no joblib import.
        # A plain generator keeps the progress-update loop below identical to
        # the parallel case (one chunk result at a time, in order).
        gen = (_mc_run_batch(base_params, rate_keys, quantities, chunk,
                             custom_network=custom_network,
                             include_nuclides=include_nuclides)
               for chunk in chunks)

    all_rows = []
    for chunk_result in gen:
        all_rows.extend(chunk_result)
        done += len(chunk_result)
        if progress:
            pct = int(100 * done / total)
            print(f"\r[MC] {done:{_w}}/{total} samples ({pct:3d}%)",
                  end='', file=sys.stderr, flush=True)

    if progress:
        print(file=sys.stderr)   # newline after final update
    return np.array(all_rows)


def mc_prev_is_reusable(prev: MCResult | None, seed: int | None, quantities: list[str],
                         params: dict[str, Any], custom_network: dict[str, Any] | None,
                         backend: str) -> bool:
    """Decide whether a previous MC result ``prev`` may be reused (extended)
    for a new :func:`mc_uncertainty`/:func:`primat.backend.run_mc` call,
    instead of being silently discarded and recomputed from scratch.

    Because MC sample ``i`` is fully determined by ``seed + i`` (see
    ``mc_uncertainty``'s docstring), ``prev`` is sample-for-sample
    compatible with a new call -- so its first ``min(len(prev), num_mc)``
    samples can be reused untouched -- iff *all* of the following hold:

    - ``prev`` was computed by the same ``backend`` ('python' or 'c') that
      will service this call.  The C and Python backends draw MC samples
      from different, non-interchangeable RNG streams (NumPy vs
      xoshiro256**), so mixing them would silently corrupt the sample
      statistics rather than raise -- this is the one condition that
      differs between the two call sites this helper replaces (the Python
      path also accepts a legacy ``prev.backend is None``, i.e. a
      :class:`MCResult` predating the ``backend`` attribute).
    - ``prev.seed == seed`` (same base seed -> same per-sample draws).
    - ``prev``'s stored quantity keys start with ``quantities`` in the same
      order (so the stacked sample columns line up column-for-column).
    - ``prev.params == params`` and ``prev.custom_network == custom_network``
      (a different network, physics flags, or reaction customisation must
      never silently reuse stale samples).

    Parameters
    ----------
    prev : MCResult or None
        The previous result to test, or ``None`` (never reusable).
    seed : int
        Base seed of the new call.
    quantities : list of str
        Quantity keys requested for the new call, in order.
    params : dict
        The new call's (already-defaulted) base parameters.
    custom_network : dict or None
        The new call's "Customise Reactions" override.
    backend : {'python', 'c'}
        Which backend will service the new call.

    Returns
    -------
    bool
        ``True`` iff ``prev`` may be reused as described above.
    """
    if prev is None:
        return False
    prev_backend = getattr(prev, 'backend', None)
    if backend == 'python':
        backend_ok = prev_backend in (None, 'python')
    else:
        backend_ok = prev_backend == 'c'
    return (backend_ok
            and getattr(prev, 'seed', None) == seed
            and list(prev)[:len(quantities)] == quantities
            and getattr(prev, 'params', None) == params
            and getattr(prev, 'custom_network', None) == custom_network)


def mc_uncertainty(num_mc: int, quantity: str | list[str] | None,
                    params: dict[str, Any] | None = None, n_jobs: int = -1,
                    seed: int | None = 0, prev: MCResult | None = None,
                    custom_network: dict[str, Any] | None = None,
                    progress: bool | None = None) -> MCResult:
    """Estimate nuclear-rate and neutron-lifetime uncertainties on BBN
    observables via Monte Carlo.

    Each MC sample draws all active nuclear rate offsets p_* independently from
    N(0,1), plus the neutron lifetime ``tau_n ~ N(cfg.tau_n, cfg.std_tau_n)``
    (which only reaches the weak-rate normalisation when
    ``cfg.tau_n_normalization=True``, the default), and runs a full primat
    solve.  By default all reactions in the selected network are varied.

    Parameters
    ----------
    num_mc : int
        Number of MC samples.
    quantity : str or list of str
        A key from the result dict ('YPBBN', 'DoH', 'He3oH',
        'Li7oH', 'Neff', 'YPCMB', ...) or a nuclide name ('H2', 'He4',
        'Li7', ...) for the final mass fraction Y.  Pass a list to evaluate
        multiple quantities in one MC pass (more efficient than separate calls).
        This only controls which quantities are *guaranteed* present and
        validated strictly (an unknown name raises); the returned
        :class:`MCResult` always additionally contains every tracked
        nuclide's final Y and every standard observable in
        ``_DEFAULT_MC_OBSERVABLES`` that this network/custom_network actually
        produces (``Neff``, ``YPBBN``, ``YPCMB``, ``He4oH``, ``DoH``,
        ``He3oH``, ``He3oHe4``, ``Li7oH``, ``Li6oLi7``, ``YCNO``), at no
        extra solving cost (each MC sample already runs a full solve). This
        keeps a TSV dump (``primat.backend.dump_mc_samples``) complete even when the
        caller only asked for one or two quantities for display purposes.
    params : dict, optional
        Base parameters for PRIMAT (e.g. Omegabh2, is_small, network).
    n_jobs : int
        Number of parallel workers passed to joblib.Parallel (-1 = all CPUs).
        ``n_jobs=1`` runs the samples serially *in-process* and does not import
        joblib at all, so a lean core install (``pip install primat``, no
        ``mc``/``recommended`` extra) can still do Monte-Carlo this way; any
        other value needs joblib and raises an actionable ImportError if it is
        missing.
    seed : int
        Base random seed; sample i uses seed + i for reproducibility.
        When evaluating on a parameter grid (e.g. scanning Ω_b h²), use the
        **same seed at every grid point** so that sample i draws the same rate
        vector p_* everywhere.  This correlates the MC noise across the grid,
        making any finite-sample bias cancel when comparing predictions at
        different parameter values.
    prev : MCResult, optional
        A previously computed result to *extend* rather than recompute from
        scratch.  Because sample ``i`` is fully determined by ``seed + i``, the
        first ``min(len(prev), num_mc)`` samples are identical to ``prev`` as
        long as the seed, the set/order of quantities, ``params`` and
        ``custom_network`` all match, so only the missing samples
        (``seed + n_prev .. seed + num_mc - 1``) are actually solved.  This
        makes it cheap to refine an estimate -- e.g. going from 30 to 50
        samples only runs the 20 new ones.  ``prev`` is silently ignored (full
        recompute) if its ``seed``, quantities, ``params``, ``custom_network``
        or ``backend`` (this function only reuses a ``prev`` whose ``backend``
        is ``"python"`` or unset -- a C-backend result has incompatible RNG
        samples, see ``primat.backend.run_mc``) differ from this call; if
        ``num_mc`` is *smaller* than ``len(prev)``, the result is just
        ``prev`` truncated to ``num_mc`` samples (nothing is solved).
    custom_network : dict, optional
        "Customise Reactions" override, forwarded to every ``PRIMAT`` instance
        built here (see :class:`PRIMAT`'s docstring for the
        ``{"removed": [...], "replaced": {...}, "added": {...}}`` schema).
        Reactions listed under ``"removed"`` are excluded from the set of
        varied rate offsets (``rate_keys``) below, since they no longer exist
        in the network; reactions listed under ``"replaced"`` stay in
        ``rate_keys`` and are varied using the *replacement* table's own error
        column (``UpdateNuclearRates`` builds ``expsigma`` from it), so a
        custom rate's uncertainty is honoured automatically; reactions listed
        under ``"added"`` are appended to ``rate_keys`` and varied the same
        way, since they are genuinely integrated (see
        :func:`_mc_resolve_rate_keys`).  ``None`` (default) uses the standard,
        uncustomised network.
    progress : bool, optional
        When truthy, print a running ``N/total (XX%)`` counter to stderr as
        samples complete, so the user can track advancement during long MC
        runs.  One update per parallel worker chunk -- for ``n_jobs=-1``
        (all CPUs) this gives roughly ``cpu_count`` updates over the full
        run.  ``None`` (default) defers to ``params['show_progress']``
        (``DEFAULT_PARAMS`` default ``True``); pass an explicit
        ``True``/``False`` to override it for this call.

    Returns
    -------
    MCResult
        Dict-like object indexed by quantity name.  Each value is an
        ``MCQuantityResult`` with attributes ``central``, ``mean``, ``std``,
        and ``values``.

    Example
    -------
    >>> mc = mc_uncertainty(100, ['YPBBN', 'DoH'], params={'Omegabh2': 0.022})
    >>> mc['YPBBN'].central
    >>> mc['YPBBN'].std
    >>> mc['DoH'].values   # full sample array
    """
    # quantity=None means "all standard observables + every tracked nuclide":
    # resolved below from the central solve, so no extra probe run is needed.
    explicit_quantities = ([] if quantity is None
                           else [quantity] if isinstance(quantity, str)
                           else list(quantity))

    base_params = dict(params or {})
    base_params.setdefault('verbose', False)
    base_params.setdefault('debug',   False)

    if progress is None:
        progress = base_params.get('show_progress', True)

    rate_keys = _mc_resolve_rate_keys(base_params, custom_network)

    # The standard observables (_DEFAULT_MC_OBSERVABLES) are always merged in
    # on top of whatever the caller explicitly requested, so the returned
    # MCResult is always complete enough to dump to disk (see that constant's
    # docstring) -- exactly like the nuclides merged in below. When reusing
    # ``prev`` we can only tell which of them are actually available for this
    # network/custom_network from ``prev`` itself (no solved instance yet);
    # when *not* reusing, availability is checked against the freshly solved
    # central_inst instead (a network without Li6/CNO simply lacks
    # Li6oLi7/YCNO, exactly as get_quantity would raise for them).
    prev_all_keys = list(prev) if prev is not None else []
    if prev is not None:
        extra_observables = [q for q in _DEFAULT_MC_OBSERVABLES
                              if q not in explicit_quantities and q in prev_all_keys]
        quantities = explicit_quantities + extra_observables
    else:
        quantities = explicit_quantities  # refined below once central_inst is solved

    # Reuse a previous result only when it is sample-for-sample compatible
    # with this call -- see mc_prev_is_reusable's docstring for the exact
    # conditions (same seed/quantities-order/params/custom_network, and a
    # Python-origin backend).
    reuse = mc_prev_is_reusable(prev, seed, quantities, base_params,
                                 custom_network, backend='python')

    (quantities, centrals, nuclide_names, nuclide_centrals,
     prev_samples, prev_nuclide_samples, n_prev) = _mc_resolve_centrals(
        reuse, prev, prev_all_keys, quantities, explicit_quantities,
        base_params, custom_network, num_mc)

    new_qty_samples, new_nucl_samples = _mc_run_new_samples(
        base_params, rate_keys, quantities, nuclide_names,
        n_prev, num_mc, seed, n_jobs, custom_network, progress)

    qty_samples = np.vstack([prev_samples, new_qty_samples])   # (num_mc, n_q)
    nucl_samples = np.vstack([prev_nuclide_samples, new_nucl_samples])  # (num_mc, n_nuclides)

    return _mc_assemble_result(quantities, centrals, qty_samples,
                               nuclide_names, nuclide_centrals, nucl_samples,
                               seed, base_params, custom_network)


def _mc_resolve_rate_keys(base_params: dict[str, Any],
                           custom_network: dict[str, Any] | None) -> list[str]:
    """Rate offsets to vary: all thermonuclear reactions in the selected network.

    Constructs a temporary config just to resolve the working directory and
    selected network filename correctly, then applies the same two
    ``custom_network`` adjustments ``UpdateNuclearRates.__init__`` applies to
    ``cfg.network``'s reaction list, so the varied set matches the set of
    reactions actually integrated:

    * ``custom_network["removed"]`` reactions are stripped (they no longer
      exist in the network, so a ``p_<rxn>`` for them would be inert);
    * ``custom_network["added"]`` reactions are appended -- in the same
      trailing position as ``UpdateNuclearRates``'s
      ``self._selected_names += added_names`` -- because a brand-new GUI-added
      reaction IS integrated and must have its rate uncertainty propagated
      like any other.  Appending (rather than interleaving) keeps every
      non-custom run's RNG stream, and hence its sample values, unchanged.
      This is also what the C sampler does, by iterating the *solved* network
      (``primat-c/src/mc.c``'s ``run_one_sample``): without the append the two
      backends would propagate different uncertainty sets for the same custom
      network.

    Note that for an ``amax``-restricted network the reaction *file* lists more
    reactions than the solved network carries (e.g. 428 vs 67 for
    ``large``+``amax=8``), so some returned keys are inert -- ``getattr(cfg,
    "p_<unused>")`` is harmless and every *live* reaction is still covered, so
    this only costs a few wasted ``standard_normal`` draws per sample.

    Returns the list of ``p_<reaction>`` parameter names to vary.
    """
    from .network_data import load_reaction_names
    tmp_cfg = PRIMATConfig(base_params)
    reactions = load_reaction_names(tmp_cfg)

    # Each entry is "bare_name" or "bare_name, filename.txt".
    # Extract only the bare_name for rate variation.
    removed = set(custom_network.get("removed", [])) if custom_network else set()
    bare_reactions = []
    for line in reactions:
        parts = re.split(r'[, ]+', line, maxsplit=1)
        bare_name = parts[0]
        if bare_name not in removed:
            bare_reactions.append(bare_name)

    # Brand-new reactions supplied by custom_network["added"], appended in the
    # same order UpdateNuclearRates appends them to the selected reaction list.
    if custom_network:
        for name in custom_network.get("added", {}):
            if name not in removed and name not in bare_reactions:
                bare_reactions.append(name)

    return [f'p_{rxn}' for rxn in bare_reactions]


def _mc_resolve_centrals(reuse, prev, prev_all_keys, quantities, explicit_quantities,
                          base_params, custom_network, num_mc):
    """Resolve the MC central values, finalised quantity list, and any
    reusable ``prev`` sample columns.

    When ``reuse`` is True (see :func:`mc_prev_is_reusable`), centrals and
    sample columns are pulled straight from ``prev`` instead of re-solving
    them (the central value, all p_*=0, does not depend on ``num_mc``).
    Otherwise runs one central :class:`PRIMAT` solve (all p_*=0) to get the
    centrals, and finalises ``quantities`` by merging in every
    ``_DEFAULT_MC_OBSERVABLES`` entry this network/``custom_network``
    actually produces (an explicitly requested but unknown quantity still
    raises via ``get_quantity``; a merged-in default observable that is
    unavailable is silently dropped).

    Returns
    -------
    (quantities, centrals, nuclide_names, nuclide_centrals,
     prev_samples, prev_nuclide_samples, n_prev)
    """
    if reuse:
        # The central value (all p_* = 0) does not depend on num_mc, so take it
        # straight from prev instead of re-solving it.
        centrals = [prev[q].central for q in quantities]
        # prev's per-quantity value arrays are the columns of the sample
        # matrix; transpose back to (n_prev, n_q) and keep at most num_mc rows
        # (a smaller num_mc simply truncates, solving nothing new).
        prev_samples = np.column_stack([prev[q].values for q in quantities])
        n_prev = min(prev_samples.shape[0], num_mc)
        prev_samples = prev_samples[:n_prev]
        # Nuclides: reuse from prev if available.  The [:n_prev] slice is
        # essential and must match the quantity block's above: when num_mc is
        # SMALLER than len(prev) (the documented "truncate, solve nothing"
        # path), an untruncated nuclide block would leave the returned
        # MCResult with two different sample counts -- every nuclide's
        # mean/std silently computed over len(prev) samples while the
        # requested quantities report num_mc, and samples_array()/cov()/
        # corr()/dump_mc_samples() raising ValueError on the ragged columns.
        nuclide_names = [q for q in prev_all_keys if q not in quantities]
        nuclide_centrals = [prev[q].central for q in nuclide_names] if nuclide_names else []
        prev_nuclide_samples = (
            np.column_stack([prev[q].values for q in nuclide_names])[:n_prev]
            if nuclide_names else np.empty((n_prev, 0)))
    else:
        # Central value (all p_* = 0).
        central_inst = PRIMAT(params=base_params, custom_network=custom_network)
        central_inst.solve()
        # The explicitly requested quantities must exist -- an unknown name
        # still raises (get_quantity), unchanged from before. The merged-in
        # default observables are silently dropped instead when unavailable
        # (e.g. Li6oLi7/YCNO on a network without Li6/CNO).
        explicit_centrals = [central_inst.get_quantity(q) for q in explicit_quantities]
        extra_observables, extra_centrals = [], []
        for q in _DEFAULT_MC_OBSERVABLES:
            if q in explicit_quantities:
                continue
            try:
                extra_centrals.append(central_inst.get_quantity(q))
                extra_observables.append(q)
            except ValueError:
                pass
        quantities = explicit_quantities + extra_observables
        centrals = explicit_centrals + extra_centrals
        nuclide_names = list(central_inst.nuclear.Y_final.keys())
        nuclide_centrals = [central_inst.nuclear.Y_final[nm] for nm in nuclide_names]
        prev_samples = np.empty((0, len(quantities)))
        prev_nuclide_samples = np.empty((0, len(nuclide_names)))
        n_prev = 0

    return (quantities, centrals, nuclide_names, nuclide_centrals,
            prev_samples, prev_nuclide_samples, n_prev)


def _mc_run_new_samples(base_params, rate_keys, quantities, nuclide_names,
                         n_prev, num_mc, seed, n_jobs, custom_network, progress):
    """Solve only the samples beyond any reused prefix (seeds ``seed+n_prev``
    through ``seed+num_mc-1``), printing an optional progress banner, and
    split the resulting column-stacked samples into quantity samples and
    nuclide samples.

    Returns ``(new_qty_samples, new_nucl_samples)``, shaped
    ``(num_mc - n_prev, n_q)`` and ``(num_mc - n_prev, n_nuclides)``.
    """
    new_seeds = [seed + i for i in range(n_prev, num_mc)]
    n_new     = len(new_seeds)

    if progress and n_new > 0:
        # Print a one-line banner before the progress counter starts so the
        # user knows the total workload, including any reused prefix.
        banner = f"[MC] Running {num_mc} sample{'s' if num_mc != 1 else ''}..."
        if n_prev > 0:
            banner += f" ({n_new} new, {n_prev} reused from prev)"
        print(banner, file=sys.stderr, flush=True)
    elif progress and n_new == 0:
        print(f"[MC] {num_mc} sample{'s' if num_mc != 1 else ''} "
              f"reused from prev (no new computation needed)",
              file=sys.stderr, flush=True)

    new_samples = _mc_collect_samples(base_params, rate_keys, quantities,
                                      new_seeds, n_jobs,
                                      custom_network=custom_network,
                                      include_nuclides=True,
                                      progress=progress and n_new > 0)
    # Parse results: first len(quantities) columns are quantities, rest are nuclides
    if new_samples.shape[0] > 0:
        new_qty_samples = new_samples[:, :len(quantities)]
        new_nucl_samples = new_samples[:, len(quantities):]
    else:
        new_qty_samples = np.empty((0, len(quantities)))
        new_nucl_samples = np.empty((0, len(nuclide_names)))

    return new_qty_samples, new_nucl_samples


def _mc_assemble_result(quantities, centrals, qty_samples, nuclide_names,
                        nuclide_centrals, nucl_samples, seed, base_params,
                        custom_network):
    """Build the :class:`MCResult` dict-like object from the (possibly
    reused+new-concatenated) quantity and nuclide sample matrices.

    Both matrices must carry the *same* number of rows (one row per MC
    sample): a ragged pair would give the nuclides a different sample count
    from the requested quantities, i.e. mean/std computed over the wrong set
    and a ValueError out of every aggregate accessor
    (:meth:`MCResult.samples_array`/:meth:`~MCResult.cov`/
    :meth:`~MCResult.corr`, hence also ``backend.dump_mc_samples``).  This
    once happened for real when a shrinking ``prev`` reuse truncated only the
    quantity block, so the invariant is asserted here rather than trusted.
    """
    if (nuclide_names and qty_samples.shape[0] != nucl_samples.shape[0]):
        raise AssertionError(
            "MC sample block mismatch: "
            f"{qty_samples.shape[0]} quantity rows vs "
            f"{nucl_samples.shape[0]} nuclide rows (must be equal)."
        )
    result_dict = {
        q: MCQuantityResult(centrals[j], qty_samples[:, j])
        for j, q in enumerate(quantities)
    }
    result_dict.update({
        nm: MCQuantityResult(nuclide_centrals[j], nucl_samples[:, j])
        for j, nm in enumerate(nuclide_names)
    })

    return MCResult(result_dict, seed=seed, params=base_params,
                    custom_network=custom_network, backend='python')
