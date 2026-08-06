# -*- coding: utf-8 -*-
"""
nuclear_network.py
==================
``NuclearNetwork`` (Class 2 of the PRIMAT split, see ``primat.background`` for
Class 1): the nuclear-reaction-network ODE integration across the HT/MT/LT
temperature eras.

Design
------
``NuclearNetwork`` is driven purely through the *minimal* public interface of
a ``primat.background.Background`` instance:

* ``T_of_t(t)`` / ``t_of_T(T)``  -- time <-> temperature
* ``rhoB_BBN(t)``                -- baryon mass density [g/cm^3] as a
  function of cosmic time (the prefactor for nuclear reaction rates)
* ``weak_nTOp_frwrd(T_K)`` / ``weak_nTOp_bkwrd(T_K)`` -- already-normalised
  n<->p weak rates [s^-1] at photon temperature ``T_K`` [Kelvin]

It knows nothing about *how* the background was constructed (NEVO table,
instantaneous decoupling, external background, scale factor, neutrino
sector, ...) -- this is exactly the seam that makes the background pluggable
(``primat.background.Background``).  In particular it does **not** use
``a_of_t``, ``Hubble``, the individual neutrino temperatures, or the NEVO
heating function: those are output-only quantities written directly by
``Background.write_time_evolution`` (see ``primat.background``), not
consumed by the nuclear solve.

``solve()`` integrates:

* **HT** (high temperature, T > T_weak ~ 1 MeV): n <-> p only.
* **MT** (mid temperature, T_weak -> T_nucl ~ 0.1 MeV): the fixed 18-reaction
  subset (n<->p + 17 reactions), regardless of network size.
* **LT** (low temperature, T_nucl -> T_end ~ 0.001 MeV): the chosen network
  (small/large, optionally amax-restricted).

and populates the public ``Y_final``, ``abundance_names`` and ``Y_of_t``
attributes consumed by ``PRIMAT``'s observable accessors (``get_quantity``,
``__getitem__``, ...) and by ``PRIMAT.solve()`` (which builds the BBN
observables dict -- ``Neff``, ``YPBBN``, ``YPCMB``, ``He4oH``, ``DoH``, ``He3oH``,
``He3oHe4``, ``Li7oH``, ``Omeganurel``, ``OneOverOmeganunr`` -- from
``Y_final`` and from ``background``'s optional neutrino-sector hooks).
"""

import os
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import zeta

from .evolution import EvolutionResult, dump_evolution


def _check_solver(sol, era, detail):
    """Raise if ``solve_ivp`` did not reach the end of the integration interval.

    Why this is needed
    ------------------
    scipy reports a step failure through ``sol.success`` / ``sol.status``, *not*
    by raising: on failure it still returns whatever partial trajectory it
    managed to compute in ``sol.y``.  Reading ``sol.y[..., -1]`` without
    checking therefore turns a solver failure into *silently wrong abundances*
    rather than an error -- the worst possible failure mode for a precision BBN
    code.  This is not hypothetical: forcing the (stiff) MT era onto scipy's
    LSODA makes it fail with "Repeated convergence failures", after which the
    unchecked final values give YP = 0.434 instead of 0.247 and a 60%-wrong
    D/H, with no warning of any kind.

    BDF converges for every supported configuration, so in normal use this
    helper never fires; it guards the cases the solver can still meet -- a
    pathological custom network, an extreme Monte-Carlo rate draw, or a
    user-tightened ``numerical_precision``.

    Backend parity
    --------------
    The C backend already fails loudly here (``primat-c/src/nuclear_network.c``
    checks ``cpr_ode_bdf``'s return code and propagates an error), and its
    Python bridge surfaces such failures as ``RuntimeError``
    (``_primat_c_src/_wrapper.c``'s ``PyErr_Format(PyExc_RuntimeError, ...)``).
    Raising ``RuntimeError`` here keeps the two backends' error behaviour
    identical.

    Parameters
    ----------
    sol : the ``OdeResult`` returned by ``scipy.integrate.solve_ivp``.
    era : str -- solver era tag used in the message ("HT", "MT" or "LT").
    detail : str -- extra context for the message (temperature range, network
        name and size), so a failure report identifies the run that produced it.

    Raises
    ------
    RuntimeError
        If ``sol.success`` is False, quoting scipy's own ``sol.message``
        (e.g. "Required step size is less than spacing between numbers.").

    Example
    -------
    >>> sol = solve_ivp(f, [t0, t1], y0, method="BDF")   # doctest: +SKIP
    >>> _check_solver(sol, "LT", "small network, 8 nuclides")   # doctest: +SKIP
    """
    if not sol.success:
        raise RuntimeError(f"[{era}] nuclear-network integration failed "
                           f"({detail}): {sol.message}")

__all__ = ["NuclearNetwork"]


class NuclearNetwork:
    """The nuclear reaction network (Class 2): HT/MT/LT ODE integration.

    Parameters
    ----------
    cfg : primat.config.PRIMATConfig
        Run-time configuration (network choice, temperature-era boundaries,
        numerical tolerances, rate-variation parameters p_*, ...).
    nucl : primat.network_data.UpdateNuclearRates
        Compiled MT/LT reaction-rate kernels (RHS + Jacobian) for the chosen
        network.
    background : primat.background.Background
        The cosmological background (Class 1) supplying ``T_of_t``/``t_of_T``,
        ``rhoB_BBN(t)``, and the normalised n<->p weak rates
        ``weak_nTOp_frwrd``/``weak_nTOp_bkwrd`` (see the module docstring for
        the full minimal interface).

    Attributes (populated by :meth:`solve`)
    ----------------------------------------
    Y_final : dict or None
        Final mass-fraction abundance ``Y`` of every nuclide in
        ``abundance_names``.
    abundance_names : list of str or None
        Tracked nuclide names, in abundance-vector order (LT species list).
    Y_of_t : scipy.interpolate.interp1d or None
        Abundance-vector interpolator ``Y(t)`` -> shape ``(len(abundance_names),)``,
        spanning HT+MT+LT.
    """

    def __init__(self, cfg, nucl, background):
        self.cfg = cfg
        self.nucl = nucl
        self.background = background
        self.Y_final = None
        self.abundance_names = None
        self.Y_of_t = None
        self._t_end = None   # cosmic time [s] at end of LT era; set by solve()
        self.evolution = None   # EvolutionResult; set by solve() iff output_time_evolution=True

    # ======================================================================
    # solve(): integrate nuclear network ODEs
    # ======================================================================

    def _era_boundaries(self):
        """Convert the four era temperature boundaries to cosmic times [s].

        Returns ``(t_start, t_weak, t_nucl, t_end)``, the HT/MT/LT era
        endpoints, via ``background.t_of_T`` (which expects a temperature in
        MeV, hence the ``/cfg.MeV_to_Kelvin`` conversion from the Kelvin
        values stored in ``cfg.T_start``/``T_weak``/``T_nucl``/``T_end``).
        """
        cfg = self.cfg
        t_of_T = self.background.t_of_T
        t_start = t_of_T(cfg.T_start / cfg.MeV_to_Kelvin)
        t_weak  = t_of_T(cfg.T_weak  / cfg.MeV_to_Kelvin)
        t_nucl  = t_of_T(cfg.T_nucl  / cfg.MeV_to_Kelvin)
        t_end   = t_of_T(cfg.T_end   / cfg.MeV_to_Kelvin)
        return t_start, t_weak, t_nucl, t_end

    def _compute_eta_b_weak(self, t_weak):
        """Baryon-to-photon ratio eta_b = n_B/n_gamma at T = T_weak.

        Evaluated once at T = T_weak from the two compulsory Background
        primitives ``rhoB_BBN(t)`` and ``t_of_T(T)``: ``_saha_YA`` (the
        Saha-equilibrium seed used by the MT era) is only ever called at
        T = cfg.T_weak, so this single value is exact -- no etab_of_T(T)
        interpolant is needed.
        """
        cfg = self.cfg
        T_weak_MeV  = cfg.T_weak / cfg.MeV_to_Kelvin
        nB_weak     = self.background.rhoB_BBN(t_weak) / (cfg.ma * cfg.MeV4_to_gcmm3)   # [MeV^3]
        ngamma_weak = (2. * zeta(3) / np.pi**2) * T_weak_MeV**3                          # [MeV^3]
        return nB_weak / ngamma_weak

    def _saha_YA(self, name, Yn, Yp, T, eta_b_weak):
        """Saha equilibrium mass-fraction abundance of nuclide `name`.

        At high temperature each nuclide is maintained in Nuclear Statistical
        Equilibrium (NSE) with free neutrons and protons via photo-dissociation.
        The Saha formula gives (Phys. Rep. §V.A):

            Y_A = g_A ζ(3)^{A-1} π^{(1-A)/2} 2^{(3A-5)/2}
                  × (M_A / mₙ^N mₚ^Z)^{3/2}
                  × (kB T)^{3(A-1)/2} η_b^{A-1}
                  × Yₙ^N Yₚ^Z exp(B_A / kB T)

        where A=N+Z is the mass number, g_A=2J+1 the spin degeneracy,
        B_A the binding energy (keV), and η_b = n_B/n_γ the baryon-to-photon
        ratio.  Used to seed the MT-era initial conditions at T = T_weak,
        where η_b = eta_b_weak (from :meth:`_compute_eta_b_weak`).

        Args:
            name : nuclide name string (key into cfg.Nuclides/NuclExcessMass).
            Yn   : free neutron mass fraction.
            Yp   : free proton mass fraction.
            T    : photon temperature in Kelvin (= cfg.T_weak at every
                   call site).
            eta_b_weak : baryon-to-photon ratio at T = T_weak (see
                   :meth:`_compute_eta_b_weak`).

        Returns:
            Y_A  : dimensionless mass fraction (≪ 1 at T ≫ BBN onset).
        """
        cfg = self.cfg
        x     = cfg.Nuclides[name]
        A     = x[0] + x[1]
        Z     = x[1]
        N     = A - Z
        Mass  = (A * cfg.ma * cfg.MeV
                 + cfg.keV * cfg.NuclExcessMass[name]
                 - Z * cfg.me * cfg.MeV)
        BindE = (N * cfg.NuclExcessMass["n"]
                 + Z * cfg.NuclExcessMass["p"]
                 - cfg.NuclExcessMass[name])
        # (M_A / mₙ^N mₚ^Z)^{3/2}: ratio of nuclear to free-nucleon masses
        NormYA = (Mass / ((cfg.mn * cfg.MeV)**(A - Z)
                          * (cfg.mp * cfg.MeV)**Z))**(3. / 2.)
        return ((2 * cfg.NuclSpin[name] + 1)
                * zeta(3)**(A - 1) * np.pi**((1 - A) / 2.)
                * 2**((3 * A - 5) / 2.)
                * NormYA
                * (cfg.kB * T)**(3. / 2. * (A - 1))
                * eta_b_weak**(A - 1)
                * Yp**Z * Yn**(A - Z)
                * np.exp(BindE * cfg.keV / (cfg.kB * T)))

    def _solve_HT(self, t_start, t_weak, _show):
        """Integrate the HT era (n <-> p only), T = T_start -> T_weak.

        Returns ``(sol_HT, Yn_HT_f, Yp_HT_f)``: the ``solve_ivp`` result and
        the final neutron/proton mass fractions, which seed the MT era.
        """
        cfg        = self.cfg
        background = self.background
        T_of_t     = background.T_of_t
        nTOp_frwrd = background.weak_nTOp_frwrd
        nTOp_bkwrd = background.weak_nTOp_bkwrd

        # Fixed era boundaries in MeV (10 / 1 MeV), used in verbose messages.
        T_start_MeV = cfg.T_start / cfg.MeV_to_Kelvin
        T_weak_MeV  = cfg.T_weak  / cfg.MeV_to_Kelvin
        if _show:
            print("[primat]  HT.", end='', file=sys.stderr, flush=True)
        if cfg.verbose:
            print(f"[nucl-py] Solving neutron decoupling at high temperature era"
                  f" (T = {T_start_MeV:.4g} -> {T_weak_MeV:.4g} MeV)")

        def Yn_i_func(T):
            b = nTOp_bkwrd(T)
            return b / (b + nTOp_frwrd(T))

        def Y_prime_HT(t, Y):
            T_K = T_of_t(t) * cfg.MeV_to_Kelvin
            f   = nTOp_frwrd(T_K)
            b   = nTOp_bkwrd(T_K)
            return b * Y[1] - f * Y[0], f * Y[0] - b * Y[1]

        Yn_i = Yn_i_func(cfg.T_start)
        Yp_i = 1. - Yn_i
        _t_ht0 = time.time()
        # HT integrator: LSODA here, Dormand-Prince RK45 in the C backend. A
        # KNOWN, accepted divergence -- recorded on both sides so it is not
        # "fixed" by accident. Same rtol (cfg.numerical_precision) and atol
        # (1e-10) on both, and the era is n <-> p only, so neither method is
        # more accurate: sweeping numerical_precision has LSODA, RK45 and BDF
        # converging to the same YPBBN. Aligning both backends on BDF was
        # tried and made cross-backend YP parity worse. Its contribution to
        # the cross-backend gap is measured by tests/backend_divergence.py.
        sol_HT = solve_ivp(Y_prime_HT, [t_start, t_weak], [Yn_i, Yp_i],
                           method='LSODA', rtol=cfg.numerical_precision, atol=1e-10)
        _check_solver(sol_HT, "HT",
                      f"T = {T_start_MeV:.4g} -> {T_weak_MeV:.4g} MeV")
        if cfg.verbose:
            print(f"[nucl-py] [HT] Finished solve_ivp in {time.time()-_t_ht0:.2f} s",
                  flush=True)
        if _show:
            print("  MT.", end='', file=sys.stderr, flush=True)
        Yn_HT_f, Yp_HT_f = sol_HT.y[0][-1], sol_HT.y[1][-1]
        return sol_HT, Yn_HT_f, Yp_HT_f

    def _solve_MT(self, t_weak, t_nucl, Yn_HT_f, Yp_HT_f, eta_b_weak, _show):
        """Integrate the MT era (fixed 18-reaction subset), T = T_weak -> T_nucl.

        Seeds every MT species except n/p (which come from the HT solution)
        at Saha (NSE) equilibrium via :meth:`_saha_YA`.  Returns
        ``(sol_MT, mt_final_raw)``: the ``solve_ivp`` result and a dict of
        final mass fractions by nuclide name, which seed the LT era.
        """
        cfg        = self.cfg
        background = self.background
        nucl       = self.nucl
        T_of_t     = background.T_of_t
        rhoB_BBN   = background.rhoB_BBN
        nTOp_frwrd = background.weak_nTOp_frwrd
        nTOp_bkwrd = background.weak_nTOp_bkwrd

        T_weak_MeV = cfg.T_weak / cfg.MeV_to_Kelvin
        T_nucl_MeV = cfg.T_nucl / cfg.MeV_to_Kelvin

        # One-slot memo on the cosmic time ``t``: BDF evaluates the RHS and the
        # (separately supplied) analytic Jacobian at the *same* ``t`` before
        # advancing the step, so ``rhoB_BBN(t)`` and ``T_of_t(t)`` -- each a
        # scipy interpolant call -- would otherwise be recomputed identically.
        # Keying on bit-identical ``t`` (no tolerance) keeps the result exact
        # while halving the background-interpolant traffic. Mirrors the same
        # pattern already used by ``NetworkDefinition.fill_buffer``.
        _bg = {"t": None, "rho": 0.0, "T_K": 0.0}
        def _bg_at(t):
            if t != _bg["t"]:
                _bg["rho"] = rhoB_BBN(t)
                _bg["T_K"] = T_of_t(t) * cfg.MeV_to_Kelvin
                _bg["t"]   = t
            return _bg["rho"], _bg["T_K"]

        def Y_prime_MT(t, Y):
            rho, T_K = _bg_at(t)
            return nucl.rhsMT(Y, T_K, rho, nTOp_frwrd, nTOp_bkwrd)

        def Jacobian_MT(t, Y):
            rho, T_K = _bg_at(t)
            return nucl.JacobianMT(Y, T_K, rho, nTOp_frwrd, nTOp_bkwrd)

        if cfg.verbose:
            print(f"[nucl-py] Solving nuclear network at mid temperature era"
                  f" (T = {T_weak_MeV:.4g} -> {T_nucl_MeV:.4g} MeV)")

        # Saha (NSE) seed for all MT species except n and p, which come from
        # the HT solution.  The MT network's species list is determined by the
        # NetworkDefinition, so this loop is independent of the network size.
        #
        # The seed is *added on top of* the HT solution's Yn + Yp = 1 rather
        # than renormalised, so the baryon budget sum_s A_s Y_s steps up by
        # sum_{A>=2} A_s Y_s^Saha(T_weak) at this handoff.  At T_weak = 1 MeV
        # every composite is still far below its BBN abundance, so the step is
        # ~1.6e-12 (measured end-to-end on `large`, amax=8: 1.0 at the start of
        # HT, 1.000000000001649 at the end of LT) -- ten orders of magnitude
        # below the last digit any observable is reported to.  Renormalising
        # here would perturb every reference number for no physical gain, so
        # the excess is documented rather than removed.
        mt_species = nucl._mt_net.species   # e.g. 8 for small, 12 for large/amax=8
        mt_saha = {"n": Yn_HT_f, "p": Yp_HT_f}
        for s in mt_species:
            if s not in mt_saha:
                mt_saha[s] = self._saha_YA(s, Yn_HT_f, Yp_HT_f, cfg.T_weak, eta_b_weak)
        Yi_MT = [mt_saha[s] for s in mt_species]

        _t_mt0 = time.time()
        sol_MT = solve_ivp(Y_prime_MT, [t_weak, t_nucl], Yi_MT,
                           method='BDF', jac=Jacobian_MT,
                           rtol=cfg.numerical_precision, atol=1e-15)
        _check_solver(sol_MT, "MT",
                      f"{cfg.network} network, {len(mt_species)} species, "
                      f"T = {T_weak_MeV:.4g} -> {T_nucl_MeV:.4g} MeV")
        if cfg.verbose:
            print(f"[nucl-py] [MT] Finished solve_ivp ({cfg.network} network, "
                  f"{len(mt_species)} species) in {time.time()-_t_mt0:.2f} s",
                  flush=True)
        if _show:
            print("  LT.", end='', file=sys.stderr, flush=True)
        # Extract MT final values by name — works for any network size.
        mt_final_raw = {s: sol_MT.y[i][-1] for i, s in enumerate(mt_species)}
        return sol_MT, mt_final_raw

    def _solve_LT(self, t_nucl, t_end, mt_final_raw, _show):
        """Integrate the LT era (chosen network), T = T_nucl -> T_end.

        Seeds the LT vector from ``mt_final_raw`` (filling any extra species
        absent from MT with 0).  Returns ``(sol_LT, finL)``: the ``solve_ivp``
        result and the final mass fractions by nuclide name (clamped to >= 0
        and zero-filled for any standard light species the chosen network
        does not track).
        """
        from .network_data import SPECIES_SMALL
        cfg        = self.cfg
        background = self.background
        nucl       = self.nucl
        T_of_t     = background.T_of_t
        rhoB_BBN   = background.rhoB_BBN
        nTOp_frwrd = background.weak_nTOp_frwrd
        nTOp_bkwrd = background.weak_nTOp_bkwrd

        T_nucl_MeV = cfg.T_nucl / cfg.MeV_to_Kelvin

        # One-slot memo on the cosmic time ``t`` (see the identical comment in
        # ``_solve_MT``): skips the duplicate ``rhoB_BBN``/``T_of_t`` interpolant
        # calls when BDF evaluates the RHS and Jacobian at the same ``t``.
        _bg = {"t": None, "rho": 0.0, "T_K": 0.0}
        def _bg_at(t):
            if t != _bg["t"]:
                _bg["rho"] = rhoB_BBN(t)
                _bg["T_K"] = T_of_t(t) * cfg.MeV_to_Kelvin
                _bg["t"]   = t
            return _bg["rho"], _bg["T_K"]

        def Y_prime_LT(t, Y):
            rho, T_K = _bg_at(t)
            return nucl.rhsLT(Y, T_K, rho, nTOp_frwrd, nTOp_bkwrd)

        def Jacobian_LT(t, Y):
            rho, T_K = _bg_at(t)
            return nucl.JacobianLT(Y, T_K, rho, nTOp_frwrd, nTOp_bkwrd)

        if cfg.verbose:
            print(f"[nucl-py] Solving nuclear network at low temperature era"
                  f" (T = {T_nucl_MeV:.4g} -> {cfg.T_end_MeV:.4g} MeV)")

        # Seed the LT vector from MT final values, filling any extra species
        # (present in the LT but absent in MT) with 0.  By looking up by name,
        # this works for any MT and LT network sizes without hardcoding.
        species_L = nucl.species_large
        Yi_LT = [mt_final_raw.get(s, 0.0) for s in species_L]

        _t_lt0 = time.time()
        # Universal LT absolute tolerance (cfg.atol_large_LT) for *every*
        # network, not just "large". Previously this was
        # `cfg.atol_large_LT if cfg.is_large else 1e-20`, i.e. keyed on the
        # literal network name -- which meant a custom network reproduced under
        # a renamed `user_nuclear_dir` overlay (is_large=False) silently used a
        # looser atol than the same network run as "large" in the GUI, breaking
        # bit-for-bit reproduction of the GUI's numbers (~1e-6). Using one atol
        # everywhere removes that name dependence. It only tightens `small`
        # (1e-20 -> 1e-26): ~8% slower, its abundances shift by ~1e-6 (a
        # tolerance artifact, not physics), and it never loosens `large`'s
        # heavy-nuclide tracking. Keep in lockstep with primat-c's
        # nuclear_network.c (bdf_opts_lt.atol).
        atol = cfg.atol_large_LT
        sol_LT = solve_ivp(Y_prime_LT, [t_nucl, t_end], Yi_LT,
                           method='BDF', jac=Jacobian_LT,
                           rtol=10.*cfg.numerical_precision, atol=atol)
        _check_solver(sol_LT, "LT",
                      f"{cfg.network} network, {len(species_L)} nuclides, "
                      f"T = {T_nucl_MeV:.4g} -> {cfg.T_end_MeV:.4g} MeV")
        if cfg.verbose:
            print(f"[nucl-py] [LT] Finished solve_ivp ({cfg.network} network, "
                  f"{len(species_L)} nuclides) in {time.time()-_t_lt0:.2f} s",
                  flush=True)
        if _show:
            print("  done.", file=sys.stderr)
        # Build LT final abundances by name.
        # Clamp to 0: the BDF solver can leave near-extinct nuclides at a
        # tiny negative value (numerical noise around zero), which is
        # unphysical for an abundance and breaks log-scale displays/ratios.
        finL = {s: max(sol_LT.y[i][-1], 0.0) for i, s in enumerate(species_L)}
        # Zero-fill only the eight *observable* light species (SPECIES_SMALL:
        # the n,p,H2,H3,He3,He4,Li7,Be7 that PRIMAT.solve() reads directly to
        # form YPBBN/D/H/He3/Li7).  This guards a custom LT network that drops
        # one of them (e.g. He4 stripped out) against a KeyError downstream,
        # *without* injecting phantom nuclides the network never evolved: the
        # small network's state vector is exactly SPECIES_SMALL, so this adds
        # nothing there and Y_final reports its true 8 nuclides -- matching
        # abundance_names (=species_L) and the C backend.  Networks that *do*
        # track the heavier He6/Li8/Li6/B8 (large, and large+amax>=8) keep them
        # because they are already in species_L above; we just no longer force
        # those four onto networks that don't (which was inflating `small` to a
        # spurious 12-nuclide Y_final).
        for s in SPECIES_SMALL:
            finL.setdefault(s, 0.0)

        if cfg.verbose:
            # Full list of every nuclide that was integrated numerically in the
            # LT era (species_L is exactly the LT solver's state vector).  The
            # list grows with the chosen network (8 / 12 / ~59 nuclides for
            # small / large, optionally amax-restricted).
            print("-" * 50)
            print(f"Primordial abundances ({len(species_L)} nuclides) at "
                  f"T = {cfg.T_end_MeV:.4g} MeV")
            print("-" * 50)
            for s in species_L:
                print(f"  Y{s:<5}= {finL[s]:.6e}")
        return sol_LT, finL

    def solve(self, progress=True):
        """
        Integrate the nuclear network over the three temperature eras.

        Populates ``self.Y_final``, ``self.abundance_names`` and
        ``self.Y_of_t`` and returns ``self.Y_final`` (the dict of final
        mass-fraction abundances by nuclide name).  The BBN observables dict
        (``Neff``, ``YPBBN``, ``DoH``, ...) is built by ``PRIMAT.solve()`` from
        ``self.Y_final`` and from ``background``'s optional neutrino-sector
        hooks -- it is no longer computed here.

        The three eras are delegated to :meth:`_solve_HT`, :meth:`_solve_MT`
        and :meth:`_solve_LT` (n<->p only / fixed 18-reaction subset / chosen
        network, see the module docstring); this method is the orchestrator
        that threads their outputs together, builds the combined abundance
        interpolator, and handles the optional outputs and DT (decay) era.

        Args:
            progress: bool, default True.  When True and ``cfg.verbose`` is
                False, print a compact one-line phase indicator to stderr
                (``[primat]  HT  MT  LT  done.``) so the user can see the
                solver is advancing without enabling full verbose output.
                Set to False inside MC workers (:func:`primat.main._mc_run_batch`)
                where per-sample dots would flood the terminal.
        """
        cfg  = self.cfg
        nucl = self.nucl

        # Refresh nuclear rates with current variation parameters (p_*, delta_*)
        nucl.apply_variations(cfg)

        # Quiet phase-progress: one compact stderr line when verbose=False so
        # the user can confirm the solver is advancing without full verbosity.
        # Suppressed when verbose=True (the verbose prints are more informative)
        # and when progress=False (set by _mc_run_batch to avoid per-sample spam).
        _show = progress and not cfg.verbose

        if cfg.verbose:
            _t0 = time.time()

        # ------------------------------------------------------------------
        # Temperature era boundaries [s]
        # ------------------------------------------------------------------
        t_start, t_weak, t_nucl, t_end = self._era_boundaries()
        self._t_end = t_end   # store for DT-era helpers and tests

        # ------------------------------------------------------------------
        # Baryon-to-photon ratio at T_weak, for the MT-era Saha (NSE) seed
        # ------------------------------------------------------------------
        eta_b_weak = self._compute_eta_b_weak(t_weak)

        # ------------------------------------------------------------------
        # HT -> MT -> LT era chain
        # ------------------------------------------------------------------
        sol_HT, Yn_HT_f, Yp_HT_f = self._solve_HT(t_start, t_weak, _show)
        sol_MT, mt_final_raw = self._solve_MT(
            t_weak, t_nucl, Yn_HT_f, Yp_HT_f, eta_b_weak, _show)
        sol_LT, finL = self._solve_LT(t_nucl, t_end, mt_final_raw, _show)

        # ------------------------------------------------------------------
        # Store final Y values for direct access (used by get_quantity)
        # ------------------------------------------------------------------
        # Use the LT species list as the canonical name list for any network.
        species_L = nucl.species_large
        self.abundance_names = species_L
        self.Y_final = dict(finL)

        # ------------------------------------------------------------------
        # Build abundance interpolator (always, so __getitem__ works)
        # ------------------------------------------------------------------
        # Each era integrates a different (growing) set of species; embed every
        # era's solution into the common abundance-vector columns *by name*, so
        # eras with fewer species (HT: n,p; MT: 12) line up with the LT layout.
        names = self.abundance_names
        col = {s: i for i, s in enumerate(names)}
        HT_names = ["n", "p"]
        MT_names = nucl._mt_net.species

        def _embed(sol_y, era_names):
            out = np.zeros((sol_y.shape[1], len(names)))
            for j, nm in enumerate(era_names):
                out[:, col[nm]] = sol_y[j]
            return out

        _t_nuc = np.concatenate((sol_HT.t, sol_MT.t[1:], sol_LT.t[1:]))
        _Y_nuc = np.vstack((_embed(sol_HT.y, HT_names),
                            _embed(sol_MT.y, MT_names)[1:, :],
                            _embed(sol_LT.y, names)[1:, :]))
        self.Y_of_t = interp1d(_t_nuc, _Y_nuc, axis=0, bounds_error=False,
                                fill_value=(0, _Y_nuc[-1]))

        # ------------------------------------------------------------------
        # Optional output: full time evolution of abundances + weak rates
        # ------------------------------------------------------------------
        if cfg.output_time_evolution:
            self._write_time_evolution(sol_HT, sol_LT, nucl)

        # ------------------------------------------------------------------
        # Optional output: two-column (nuclide, final abundance Y) table
        # ------------------------------------------------------------------
        if cfg.output_final_result:
            self._write_final_result()

        # ------------------------------------------------------------------
        # Decay Time (DT) era (optional; needs a network carrying decays)
        # ------------------------------------------------------------------
        # After BBN ends at t_end, long-lived radioactive isotopes (C14, Be10,
        # Na22, …) continue to decay on timescales of years to millions of
        # years.  The DT era propagates the abundance vector forward in time
        # using only the constant decay matrix (no Hubble expansion, no
        # thermal production), via matrix exponentiation:
        #   Y(t) = exp(D × Δt) × Y(t_end)
        # where D is the (constant) decay-rate matrix assembled from the decay
        # reactions in the LT network (see _build_decay_matrix).
        #
        # Gated on the network *actually carrying* a decay reaction, not on the
        # literal name "large" (cfg.is_large).  Name-keying was the same trap
        # the LT atol above was de-named to escape: a large-equivalent network
        # reproduced under a renamed user_nuclear_dir overlay carries exactly
        # the same decays.txt reactions, yet silently got no DT era at all.
        # weak_indices holds the n__p entry (index 0) plus every decay, so
        # "more than just n__p" is precisely "this network has decays".
        if cfg.decay_era and len(nucl._lt_net.weak_indices) > 1:
            Y0_DT = np.array([self.Y_final.get(s, 0.0) for s in self.abundance_names])
            D = self._build_decay_matrix(nucl._lt_net)
            t_decay_end = cfg.t_decay_end
            decay_n     = cfg.decay_n_points
            # Time grid log-spaced in the *elapsed* time Δt = t − t_end (not in
            # absolute t).  This is essential: the residual free neutron decays
            # with τ_n ≈ 880 s, a transient ~10 decades shorter than t_end
            # (~1.3×10⁶ s).  A grid log-spaced in absolute t would put its first
            # interior point ~10⁵ s past t_end, completely skipping the neutron
            # decay (linear interpolation between t_end and that point would
            # flatten it).  Spacing in Δt from Δt_min = 1 s gives dense sampling
            # immediately after t_end (resolving n, and any other fast residual)
            # while still reaching t_decay_end with coarse late-time sampling
            # for the slow decays (Na22, C14, Be10).
            t_DT = t_end + np.logspace(np.log10(1.0),
                                       np.log10(t_decay_end), decay_n)
            Y_DT = self._integrate_decay_era(D, Y0_DT, t_end, t_DT)
            if cfg.verbose:
                print(f"[nucl-py] [DT] Decay era: {decay_n} time points from "
                      f"t={t_end:.3g} s to t={t_end + t_decay_end:.3g} s")
                for i, s in enumerate(self.abundance_names[:12]):
                    if Y_DT[-1, i] > 0:
                        print(f"  Y{s:<5}= {Y_DT[-1, i]:.6e}")

            # Extend the public Y(t) interpolator across the DT era so that
            # callers (``run[species](t)``, ``get_quantity(..., t=...)``) see a
            # single seamless history t_start … t_end+t_decay_end, exactly like
            # the HT→MT→LT concatenation above.  t_DT[0] = t_end+1 > t_end, so
            # the appended grid stays strictly increasing; Y_DT is already in
            # ``abundance_names`` column order (it is built from Y0_DT, which is
            # itself indexed by ``abundance_names``), matching _Y_nuc's layout.
            # The right-hand fill_value becomes the late-time DT value (the
            # fully-decayed state) instead of the LT endpoint.
            _t_nuc = np.concatenate((_t_nuc, t_DT))
            _Y_nuc = np.vstack((_Y_nuc, Y_DT))
            self.Y_of_t = interp1d(_t_nuc, _Y_nuc, axis=0, bounds_error=False,
                                    fill_value=(0, _Y_nuc[-1]))

            if cfg.output_decay_evolution:
                self._write_decay_evolution(t_DT, Y_DT)

        return self.Y_final

    def _lt_t_end_s(self):
        """Return the cosmic time [s] at the end of the LT era.

        Populated by :meth:`solve`.  Used by DT-era helpers and tests to anchor
        the Δt = t − t_end offset for the decay matrix exponentiation.

        Returns
        -------
        float
            Cosmic time at T_end [s]; e.g. ~1.3×10^6 s (≈ 15 days) for the
            default T_end_MeV = 0.001.

        Raises
        ------
        RuntimeError
            If called before :meth:`solve`.
        """
        if self._t_end is None:
            raise RuntimeError("_lt_t_end_s() called before solve()")
        return self._t_end

    def _write_final_result(self):
        """Write a two-column ``nuclide  Y`` table of final abundances.

        Dumps every tracked nuclide of the active network and its final
        mass-fraction abundance ``Y`` at the end of BBN to
        ``cfg.output_final_file``.  ``Y`` is normalised so that
        ``sum_s A_s Y_s = 1`` (A = mass number), i.e. it is the per-baryon
        abundance weighted by A.  The rows are exactly the species of the
        chosen network: 8 for ``small``, ~59 for ``large`` (fewer with an
        ``amax`` cutoff, e.g. 12 for ``large, amax=8``),
        in abundance-vector order (``n`` and ``p`` first).

        Enabled by ``output_final_result=True``; the destination is
        ``output_final_file`` (relative paths resolve against the current
        working directory, like ``output_file``).  Typical use -- get the
        full nuclide vector of a
        single run without going through ``get_quantity`` for each name::

            PRIMAT(params={'output_final_result': True,
                              'output_final_file': 'results/run_final.dat',
                              'network': 'large'}).solve()

        produces a file whose first lines read::

            nuclide       Y
            n             4.032109e-16
            p             7.530243e-01
            H2            1.835287e-05
            ...
        """
        from .backend import dump_final_with_sigma

        cfg  = self.cfg
        # Resolve relative paths against the current working directory (the
        # universal convention), same rule as output_file.
        path = os.path.abspath(cfg.output_final_file)
        out_dir = os.path.dirname(path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        names = self.abundance_names
        with open(path, 'w') as f:
            f.write(dump_final_with_sigma(names, self.Y_final))
        # Always announce: this file is written only on explicit request
        # (output_final_result=True), so the user wants to know where it landed.
        print(f"[output] Final abundances ({len(names)} nuclides) written to {path}")

    def _write_time_evolution(self, sol_HT, sol_LT, nucl):
        """Build the unified ``EvolutionResult`` (see :mod:`primat.evolution`)
        and write it to ``cfg.output_file``.

        Enabled by ``output_time_evolution=True``.  Always sets
        ``self.evolution`` to the in-memory result (no disk I/O required to
        get it -- e.g. ``PRIMAT(...).solve()["evolution"]``); additionally
        writes ``cfg.output_file`` as a convenience via
        :func:`primat.evolution.dump_evolution`.  Works for all three
        networks (``small``/``large``, optionally ``amax``-restricted) --
        the ``Y_<species>`` columns are derived from
        ``self.abundance_names`` (8 / 12 / ~59 nuclides).

        Columns (see :mod:`primat.evolution` for the exact schema): cosmic
        time ``t_s``, scale factor ``a``, photon temperature
        ``T_gamma_MeV``, the three flavour neutrino temperatures, and one
        ``Y_<species>`` column per tracked nuclide (mass-fraction
        abundance). ``a``/the neutrino temperatures come from
        ``self.background`` (``np.nan`` if it tracks no scale factor/
        neutrino sector, e.g. a minimal custom background).

        Before the nuclear network starts integrating a given species (the
        HT era for everything but n/p, and the time before ``T_start_cosmo``
        for n/p too), its ``Y_<species>`` column is **exactly 0** -- the
        value ``_embed``/``Y_of_t`` produce there.  A previous version of
        this method filled that region with the Nuclear Statistical
        Equilibrium (Saha) prediction ``YA(name, Yn, Yp, T)`` for a smoother
        log-log plot; this was removed because the fill is *often wrong*:
        NSE need not hold for every nuclide all the way down to
        ``T_start_cosmo`` (e.g. for non-standard backgrounds with extra
        entropy injection or a non-thermal neutrino sector), and a
        silently-injected equilibrium value is worse than an honest 0 that a
        plotting tool can choose to mask. Consumers that want a smooth
        pre-MT curve can compute the Saha value themselves from the ``t_s``/
        ``T_gamma_MeV`` columns.

        Per-reaction forward-rate columns (``<reaction>_frwrd``) are appended
        after the ``Y_<species>`` block when
        ``cfg.output_rates_time_evolution=True`` -- one column per reaction in
        the active LT network (~12 for ``small``/``small_parthenope``, 68 for
        ``large``+``amax=8``, ~429 for full ``large``). Each column is the
        active forward reaction-rate interpolant (same units as the shipped
        rate tables) at the row's temperature, populated into
        ``EvolutionResult.rates`` and serialised by ``dump_evolution``. The
        n<->p weak rates are not duplicated on disk: recover them from
        ``run.background.weak_nTOp_frwrd``/``weak_nTOp_bkwrd`` evaluated at the
        ``T_gamma_MeV`` column. The C backend emits the identical rate columns.
        The richer background-only TSV (``H``, ``Nheating``, energy
        densities, ...) is still written separately by
        ``background.write_time_evolution``/``time_evolution_text`` when
        ``cfg.output_background_evolution=True`` (see
        :mod:`primat.background`).
        """
        cfg = self.cfg
        background = self.background
        # Derive column names from the actual abundance names so custom networks
        # (which may have fewer or different nuclides than the standard 8 or 12)
        # produce a result with the correct number of columns.
        names = self.abundance_names

        Y_of_t = self.Y_of_t

        # Uniform log-spaced output grid from T_start_cosmo to end of LT era
        t_cosmo = background.t_of_T(cfg.T_start_cosmo / cfg.MeV_to_Kelvin)
        t_end   = sol_LT.t[-1]
        t_out   = np.logspace(np.log10(t_cosmo), np.log10(t_end), cfg.output_n_points)

        T_out = background.T_of_t(t_out)
        a_out = (background.a_of_t(t_out) if background.has_scale_factor
                 else np.full_like(t_out, np.nan))
        Tnu = background.Tnu_of_t(t_out)
        if Tnu is None:
            nan_col = np.full_like(t_out, np.nan)
            Tnu = {"e": nan_col, "mu": nan_col, "tau": nan_col}

        # Abundances: zero before nuclear network starts (Y_of_t's fill_value)
        t_start = sol_HT.t[0]
        Y_out = np.zeros((len(t_out), len(names)))
        mask_nuc = t_out >= t_start
        Y_out[mask_nuc] = Y_of_t(t_out[mask_nuc])
        Y = {s: Y_out[:, j] for j, s in enumerate(names)}

        # Optional per-reaction forward-rate columns. One
        # <reaction>_frwrd column per reaction actually in the active LT
        # network (whatever the network / amax cutoff selects: ~12 for
        # small/small_parthenope, 68 for large+amax=8, ~429 for full large),
        # value = the active forward-rate interpolant at each output
        # temperature (plain rate, not a flux). Sorted by column name so the C
        # backend can emit the identical names in the identical order.
        # Computed directly from the LT rate table
        # (the same linear interpolation on the master T9 grid that the ODE
        # right-hand side and the C backend's cpr_network_fill_buffer use), so
        # it works for any reaction the active network carries.
        rates = None
        if cfg.output_rates_time_evolution:
            lt = nucl._lt_net
            g = lt.grid
            # T9 = T[K] / 1e9; searchsorted(g, T9) - 1 clamped to [0, len-2]
            # (edge-clamped interval index), vectorised over the output grid.
            T9 = T_out * cfg.MeV_to_Kelvin * 1e-9
            idx = np.clip(np.searchsorted(g, T9) - 1, 0, g.size - 2)
            w = (T9 - g[idx]) / (g[idx + 1] - g[idx])
            rates = {}
            # lt.names[0] is the prepended weak n__p (no fwd row); reaction i
            # (i >= 1) maps to forward-rate row i-1 of lt._fwd.
            for i in range(1, len(lt.names)):
                fwd_row = lt._fwd[i - 1]
                rates[f"{lt.names[i]}_frwrd"] = (
                    fwd_row[idx] * (1.0 - w) + fwd_row[idx + 1] * w)
            rates = {k: rates[k] for k in sorted(rates)}

        self.evolution = EvolutionResult(t=t_out, a=a_out, T_gamma=T_out,
                                         T_nu=Tnu, Y=Y, rates=rates)

        # cfg.output_file=None is the in-memory-only escape hatch (e.g.
        # primat-gui's _solve): self.evolution above is what
        # that caller actually wants, with no disk I/O at all -- this is the
        # only output_*=True flag in the package with that escape hatch,
        # since it is also the only one a hosted GUI needs to suppress.
        if cfg.output_file is None:
            return

        # Resolve relative paths against the current working directory (the
        # universal convention), not the installed-package directory.
        out_path = os.path.abspath(cfg.output_file)
        dump_evolution(self.evolution, out_path)

        # Always announce: written only on explicit request (output_time_evolution=True).
        print(f"[output] Time-evolution data ({len(t_out)} rows) written to {out_path}")

    # ======================================================================
    # Decay Time (DT) era helpers
    # ======================================================================

    def _build_decay_matrix(self, net):
        r"""Build the constant decay-rate matrix D for the DT era.

        In the DT era all thermonuclear reactions are frozen (T is too low for
        any thermal activation), so only radioactive decays remain.  The
        abundance vector Y evolves as:

            dY/dt = D · Y

        where D is a constant N×N matrix (N = number of nuclides in the LT
        network).  Each decay reaction ``X → P1 + P2 + ... + B±`` contributes:

            D[X_idx, X_idx] -= rate_X × mult_X        (loss term for parent X)
            D[P_idx, X_idx] += rate_X × mult_P        (gain term per product P)

        **Convention (important).**  ``Y`` is the *number* abundance per baryon,
        ``Y_s = n_s / n_B``, normalised so that ``Σ_s A_s Y_s = 1`` — it is
        **not** a mass fraction, despite the loose "mass fraction" wording used
        elsewhere in this file.  That is the convention the LT/MT right-hand
        side itself uses: ``network_builder._rhs_kernel`` applies the bare
        integer stoichiometry ``af_co = c_prod − c_react`` with no mass
        weighting, and ``network_builder.check_conservation`` verifies exactly
        ``Σ_s A_s ΔY_s = 0``.

        The gain term is therefore the bare multiplicity ``mult_P``, with **no**
        ``A_P / A_X`` factor.  A previous version carried such a factor (on the
        mistaken premise that Y was a mass fraction) *in addition to* ``mult_P``,
        which broke baryon conservation for every decay whose products differ in
        mass number from the parent: ``Li8 → α + α`` produced ``+λ`` of He4
        instead of ``+2λ``, i.e. half the alphas, and ``C9 → α + α + p`` lost 4/9
        of the baryon number.  The 33 ordinary β decays (``A_P = A_X``, e.g.
        ``C14 → N14``) were unaffected, which is why the error went unnoticed.

        Photons and leptons (Bm/Bp) are excluded from the ODE state vector; only
        nuclear species (those in ``net.species``) appear in D.

        In addition to the ``decays.txt`` reactions, the free-neutron β decay
        ``n → p`` is added explicitly with rate ``1/cfg.tau_n``: it is the
        T→0 limit of the thermal n↔p weak rate (which is handled by the
        background during HT/MT/LT, not stored as a decay table), so without it
        the residual free neutrons at ``t_end`` would never decay.

        Parameters
        ----------
        net : NetworkDefinition
            The LT network (``nucl._lt_net``); supplies species, N, Z,
            stoichiometry (``net.network``), decay-reaction flags
            (``net.weak_indices``), and rate tables (``net._fwd_median``).

        Returns
        -------
        D : np.ndarray, shape (N, N)
            Decay-rate matrix in [s^-1].  Off-diagonal entries are ≥ 0;
            diagonal entries are ≤ 0.

        Notes
        -----
        Baryon-number conservation: ``Σ_s A_s D[s, X] = 0`` for every parent
        column X.  This holds *exactly* (not approximately): the emitted
        leptons and photons carry A = 0, so they remove no baryon number, and
        every decay in ``decays.txt`` balances A between its parent and its
        nuclear products.  ``tests/test_nuclear.py`` pins it for the whole
        ``large`` network — it is the check that would have caught the
        ``A_P/A_X`` bug described above.

        Example
        -------
        ``C14 → N14 + Bm`` with rate λ (a mass-preserving β decay):

            D[C14, C14] = -λ      (C14 is lost)
            D[N14, C14] = +λ      (N14 is gained)

        ``Σ_s A_s D[s, C14] = 14×(-λ) + 14×(+λ) = 0`` ✓

        ``Li8 → α + α + Bm``, where the multiplicity — not any mass ratio — is
        what closes the budget:

            D[Li8,  Li8] = -λ     (one Li8 lost)
            D[He4,  Li8] = +2λ    (two alphas gained: mult_P = 2)

        ``Σ_s A_s D[s, Li8] = 8×(-λ) + 4×(+2λ) = 0`` ✓
        """
        N = len(net.species)
        D = np.zeros((N, N))

        # The rate tables (_fwd_median) are indexed without the n__p slot:
        # names[0] = "n__p", names[1:] = thermonuclear reactions.
        # _fwd_median[i] corresponds to names[i+1], so we need offset by 1.

        for rxn_idx in net.weak_indices:
            if rxn_idx == 0:
                continue   # n__p handled by the HT/MT/LT eras, not the DT era
            name = net.names[rxn_idx]
            # The decay rate is constant (T9-independent), stored as a
            # uniform array.  Read from _fwd_median at grid index 0.
            # rate_table_idx is rxn_idx - 1 because _fwd_median excludes n__p.
            rate = float(net._fwd_median[rxn_idx - 1, 0])   # [s^-1]
            if rate == 0.0:
                continue

            react, prod = net.network[rxn_idx]   # {species_idx: multiplicity}

            # Parent nuclide: the sole nuclear reactant (multiplicity 1 for all
            # beta/EC decays; multi-nucleon decays like Li8→α+α+Bm are handled
            # via the products dict below).
            for X_idx, X_mult in react.items():
                # Loss term for the parent X
                D[X_idx, X_idx] -= rate * X_mult

                # Gain terms for nuclear products (the lepton Bm/Bp and photons
                # are excluded from the ODE state and are already absent from
                # net.network's index-based stoichiometry).  Y is the number
                # abundance per baryon, so the gain is the bare multiplicity:
                # dY_P/dt = rate × mult_P × Y_X, with no A_P/A_X weighting (see
                # the "Convention" paragraph in this method's docstring).
                for P_idx, P_mult in prod.items():
                    D[P_idx, X_idx] += rate * P_mult

        # ------------------------------------------------------------------
        # Free-neutron β decay  n → p + e⁻ + ν̄
        # ------------------------------------------------------------------
        # The n→p transition is *not* a decays.txt entry: during BBN it is the
        # thermal weak rate (n__p, rxn_idx 0, T-dependent, computed by the
        # background) and is therefore skipped above.  In the DT era T→0, so
        # that thermal rate reduces to the vacuum decay constant λ_n = 1/τ_n
        # (τ_n = cfg.tau_n, the neutron lifetime).  Without this term the
        # residual free neutrons surviving at t_end (Y_n ~ 4×10⁻¹⁶) would be
        # frozen for all of cosmic time instead of decaying to protons within
        # ~minutes; including it lets the DT era track n correctly.  One
        # neutron makes one proton, so the gain is a bare +lam_n.
        if "n" in net.species and "p" in net.species:
            n_idx = list(net.species).index("n")
            p_idx = list(net.species).index("p")
            lam_n = 1.0 / self.cfg.tau_n   # [s^-1]; τ_n = neutron lifetime
            D[n_idx, n_idx] -= lam_n
            D[p_idx, n_idx] += lam_n

        return D

    def _integrate_decay_era(self, D, Y0, t_end, t_grid):
        r"""Propagate abundances through the DT era via matrix exponentiation.

        The DT (Decay Time) ODE ``dY/dt = D · Y`` with constant coefficient
        matrix D is solved exactly by:

            Y(t) = exp(D × (t − t_end)) · Y_0

        We form the dense matrix exponential with ``scipy.linalg.expm`` (Padé
        approximation with *scaling-and-squaring*) and apply it to Y0:

            Y(t_i) = expm(D × Δt_i) @ Y0

        where Δt_i = t_i − t_end is the elapsed time since BBN end.

        **Why not ``scipy.sparse.linalg.expm_multiply``?**  The decay matrix has
        a colossal eigenvalue spread: the fastest decay (B15, T½ ≈ 10 ms) gives
        an eigenvalue ~70 s⁻¹, while Δt reaches ~1 Gyr ≈ 3×10¹⁶ s, so
        ‖D·Δt‖ ~ 10¹⁸.  ``expm_multiply`` selects its number of internal
        matrix–vector products *linearly* in ‖D·Δt‖ (Al-Mohy & Higham 2011,
        Eq. 3.6), so for this norm it attempts ~10¹⁸ products and effectively
        hangs.  ``scipy.linalg.expm`` instead uses scaling-and-squaring whose
        cost grows only *logarithmically* in ‖D·Δt‖ (≈ log₂‖D·Δt‖ ~ 60
        squarings), so it handles the full 16-decade spread in milliseconds.
        Since D is small (N ≤ 60), forming the dense exp(D·Δt) is cheap
        (~3 ms per time point, ~0.1 s for the default 200-point grid).

        Parameters
        ----------
        D : np.ndarray, shape (N, N)
            Decay-rate matrix from :meth:`_build_decay_matrix` [s^-1].
        Y0 : np.ndarray, shape (N,)
            Initial abundance vector at t = t_end (end of LT era).
        t_end : float
            Cosmic time at end of BBN / start of DT era [s].
        t_grid : np.ndarray, shape (M,)
            Output times [s], all > t_end; log-spaced from solve().

        Returns
        -------
        Y_t : np.ndarray, shape (M, N)
            Abundance vectors at each output time.  Row i is Y(t_grid[i]).

        Notes
        -----
        D's eigenvalues are the negative decay constants (≤ 0), so exp(D·Δt)
        is a contraction and the result is numerically stable for any
        positive Δt.

        References
        ----------
        Al-Mohy & Higham (2009), "A New Scaling and Squaring Algorithm for the
        Matrix Exponential", SIAM J. Matrix Anal. Appl. 31, 970–989 (the
        algorithm behind ``scipy.linalg.expm``).
        """
        from scipy.linalg import expm

        N_t = len(t_grid)
        N   = len(Y0)
        Y_t = np.zeros((N_t, N))

        for k, t_k in enumerate(t_grid):
            dt = t_k - t_end   # elapsed time since end of BBN [s]
            # expm(D*dt) @ Y0 computes the exact solution Y(t_k) of dY/dt = D·Y.
            Y_t[k] = expm(D * dt) @ Y0
            # Clip small negative values that arise from floating-point
            # cancellation (the matrix exp may produce tiny negatives for
            # species whose abundance is near zero).
            np.clip(Y_t[k], 0.0, None, out=Y_t[k])

        return Y_t

    def _write_decay_evolution(self, t_grid, Y_t):
        """Write the DT-era abundance time series to a TSV file.

        Enabled by ``cfg.output_decay_evolution=True``; the destination is
        ``cfg.output_decay_file`` (relative paths resolve against the current
        working directory).

        Columns: ``t`` [s], then one ``Y<species>`` column per tracked
        nuclide in ``self.abundance_names`` (in abundance-vector order).

        Parameters
        ----------
        t_grid : np.ndarray, shape (M,)
            Output times [s] (log-spaced from t_end to t_end + t_decay_end).
        Y_t : np.ndarray, shape (M, N)
            Abundance vectors at each time, from :meth:`_integrate_decay_era`.
        """
        cfg  = self.cfg
        path = os.path.abspath(cfg.output_decay_file)
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

        nuc_cols = ["Y" + s for s in self.abundance_names]
        out_data = np.column_stack([t_grid, Y_t])
        out_header = "\t".join(["t"] + nuc_cols)
        np.savetxt(path, out_data, delimiter='\t', header=out_header, comments='')
        print(f"[output] Decay-era evolution ({len(t_grid)} rows) written to {path}")
