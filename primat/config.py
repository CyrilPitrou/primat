# -*- coding: utf-8 -*-
"""
config.py
=========
Central configuration for primat.

Physical constants and derived unit conversions are *fixed* and computed once
here.  All run-time flags and cosmological/nuclear parameters are carried in a
``PRIMATConfig`` instance and can be overridden by passing a parameter dictionary
to ``PRIMATConfig(params)``.

No file I/O happens here.  Nuclear rate data are loaded separately in
``nuclear_data.py``.
"""

import difflib
import numbers
import os
import re
import warnings
from typing import TYPE_CHECKING
import numpy as np

from .constants import CONST

__all__ = ['DEFAULT_PARAMS', 'PARAM_GROUPS', 'PRIMATConfig']

# ---------------------------------------------------------------------------
# Default parameter values exposed as a plain dict so callers can inspect them
# ---------------------------------------------------------------------------
DEFAULT_PARAMS: dict = {
    # ---- general behaviour and numerical settings ------------------------------------------------
    "verbose":               False, #If you want the messages from the code to be printed, set this to True.  This is separate from the debug, which controls the printing of extra messages for debugging purposes.
    "debug":                 False, #If you want the debug messages to be printed, set this to True.  This is separate from the verbose, which controls the printing of general messages from the code.
    "show_progress":         True,  # Set to False to hide the compact stderr progress indicators printed
    # when verbose=False: the "[primat]  HT.  MT.  LT.  done." phase markers from a single solve,
    # and the "[MC] Running N samples..." banner / "[MC] i/N samples (XX%)" counter from an MC run.
    "numerical_precision":        1.e-7, # for finite differences (solve_ivp). 1e-6 should be enough.
    "numba_installed":            True,  # will be re-checked at runtime. Allows just-in-time compilation for faster execution.
    "strict_params":              False, # How PRIMATConfig reacts to an unknown parameter key (a typo like "Omegab2h", or a key from a different code).
    # False (default): warn and ignore it, appending difflib "did you mean ...?" suggestions so the mistake is visible without being fatal. True: raise ValueError on the first unknown key 

    # ---- physics settings ------------------------------------------------
    # ---- neutrino decoupling ----------------------
    "incomplete_decoupling":      True, # True: non-instantaneous neutrino decoupling, read from the pre-computed NEVO table.
    # False: instantaneous decoupling (Tnu/Tgamma fixed by EM entropy conservation; see neutrino_history.InstantaneousDecoupling).
    # incomplete_decoupling=False with spectral_distortions=True (NEVO-based) is physically inconsistent and rejected; see PRIMATConfig validation.

    # ---- electromagnetic plasma -------------------
    "QED_corrections":            True,  # Whether to include QED interaction corrections to the EM plasma equation of state.
    "n_electron_table":           2000,  # number of log-spaced grid points for the electron-thermo (rho_e/p_e and derivatives) tables
    "recompute_electron_thermo":  False, # If False, load cache_plasma_weak/plasma/electron_thermo_cache.txt (via the cache_dir overlay) when its fingerprint matches; otherwise (or if True) recompute and overwrite it. See plasma.Plasma._build_electron_tables.
    "recompute_qed_corrections":  False, # True: always compute analytically and overwrite cache_plasma_weak/plasma/QED_*.txt (via the cache_dir overlay); False: load from files if present, otherwise compute on the fly without saving

    # ---- spectral distortions ---------------------
    "spectral_distortions":       True, # Corrections to n<->p weak rates from deviations of the neutrino phase-space distribution from a perfect Fermi-Dirac.
    # Two sub-modes, selected by analytic_distortions (see neutrino_history.py):
    #   False (default): read the distortion from the full NEVO spectrum file
    #     (86-column, not _col_1_7); requires incomplete_decoupling=True.
    #   True: analytic y-type (SZ) + gray-type distortion controlled by
    #     y_SZ/y_gray, also contributing rho_nuSD to the Friedmann equation.
    # NOTE: there is deliberately no mu-type (chemical-potential) spectral
    # distortion -- a neutrino chemical potential is not a spectral distortion;
    # use munuOverTnu instead (it shifts the weak rates AND the energy density).
    "analytic_distortions":       False, # False (default): read the distortion from the tabulated NEVO spectrum (requires incomplete_decoupling=True); True: use the analytic y_SZ/y_gray shapes (requires incomplete_decoupling=False).
    "y_SZ":                       0., # Amplitude of the y-type (Sunyaev-Zel'dovich-like, Compton) distortion; see neutrino_history.AnalyticDistortion.
    "y_gray":                     0., # Amplitude of the gray-type (gray-body temperature-rescaling) distortion: delta_f(y) = -fd(y) + fd(y/(1+y_gray))/(1+y_gray)**3.
    # ENERGY density shifts linearly, integral{y^3 delta_f dy} = y_gray * 7*pi**4/120 exactly -- a distinct, independent third distortion shape.

    # ---- custom NEVO tables ------------------------------------------------
    # Override the shipped rates/NEVO/ tables with custom ones (e.g. a
    # higher-resolution or non-standard neutrino-decoupling history).  Each is
    # a filename resolved relative to rates/NEVO/, or an absolute path; None
    # uses the shipped file selected by QED_corrections (see
    # neutrino_history.NEVOTable / resolve_nevo_path).
    "nevo_file":                  None, # 6/7-column thermo table (replaces NEVOPRIMAT[_NoQED]_col_1_7.csv)
    "nevo_spectral_file":         None, # 86-column spectral-distortion table (replaces NEVOPRIMAT[_NoQED].csv); only read when spectral_distortions=True and analytic_distortions=False
    "nevo_grid_file":             None, # y-grid for nevo_spectral_file (replaces NEVOGrid.csv); its length must match nevo_spectral_file's column count minus 6
    "nevo_file_prefix":           "NEVOPRIMAT", # base name for the *default* NEVO thermo/spectral
    # files: "<prefix>[_NoQED]_col_1_7.csv" (thermo) and "<prefix>[_NoQED].csv" (86-col
    # spectral). NEVOGrid.csv is NOT prefixed (shared y-grid). Ignored for any file
    # selected explicitly via nevo_file/nevo_spectral_file (those still win), and has no
    # effect when incomplete_decoupling=False (no NEVO file is read at all).

    # ---- data directory override and nuclear overlay -----------------------
    # See PRIMATConfig.resolve_rates_path and _resolved_data_dir. Both default
    # to None (shipped primat/data/ tree). When data_dir is set, it completely
    # replaces the shipped data tree (NEVO/, nuclear/, csv/, cache_plasma_weak/
    # must all be present under that directory; the regenerable weak-rate and
    # plasma caches live together under cache_plasma_weak/{weak,plasma}/, and
    # can be redirected to a writable location via cache_dir -- see below).
    # When user_nuclear_dir is set it is an additive overlay for nuclear networks and rate tables only
    # (checked before the shipped tree, so "small"/"large" remain available
    # even if only user_nuclear_dir is set and it doesn't contain them).
    # Overlay roots for user_nuclear_dir are treated as the equivalent of
    # primat/data/nuclear, so they should contain `networks/` and `tables/`
    # directly.
    "data_dir":          None, # Full-takeover data directory (must exist if set; replaces primat/data/)
    "user_nuclear_dir":  None, # Additive overlay for nuclear networks & rate tables (must exist if set)

    # ---- background mode ---------------------------------------------------
    "external_scale_factor":      False, # If True, read the scale factor a(T_gamma) directly
    # from the NEVO table's x column (a is proportional to x by the NEVO convention)
    # instead of solving the entropy-conservation ODE from the heating
    # function N_NEVO. t(a) is still obtained by Hubble integration (unchanged). Outside
    # the table's T range, both modes extrapolate assuming radiation domination
    # (a ~ 1/T, t ~ 1/T^2). Requires incomplete_decoupling=True.

    "custom_background":         None, # Path (str) to a user-supplied background file
    # containing at minimum the columns T [MeV], t [s], and a (scale factor, normalised
    # so that a·T_γ → T0CMB_MeV as T → 0, i.e. a = 1 today). The file must be
    # tab- or comma-delimited with a header row. Extra columns are ignored.  When set,
    # incomplete_decoupling and spectral_distortions are forced to False (with warnings
    # if they were True); the n<->p weak rates use the instantaneous-decoupling
    # approximation (T_ν(T_γ) from EM entropy conservation). Neff is estimated via the
    # Friedmann equation from the supplied a(t). Incompatible with external_scale_factor.

    # ---- fundamental constants (overridable for sensitivity studies) --------
    # CODATA-tabulated value, kept as the exact 5-significant-figure literal.
    # Do NOT replace with the result of converting some natural-units value
    # through CONST.GN_MeV2_to_SI -- that round-trips to a spurious 16-digit
    # float (6.674299257609439e-11, off from the tabulated constant at the
    # ~1e-7 relative level) which previously crept in here this way. Mirror
    # this literal digit-for-digit in primat-c/src/config.c's
    # cpr_config_set_GN default (see CLAUDE.md's primat/primat-c sync rule).
    "GN":                         6.6743e-11,   # Newton's constant, SI units [m^3 kg^-1 s^-2]

    # ---- background thermodynamics ----------------------------------------
    "T_start_cosmo_MeV":          40.0, # photon temperature at which the background integration starts [MeV]; must be > T_end_MeV. 40 MeV is well before any BBN-relevant weak or nuclear process (T_weak = 1 MeV), so the initial condition is pure radiation domination.
    "T_end_MeV":                  1.e-3,  # end temperature for nuclear integration [MeV]; default 0.001 MeV ≈ 11.6 MK
    "sampling_temperature_per_decade": 600,  # points per decade of T for the background a(T)/t(T) grid

    # ---- n <--> p weak rates ----------------------------------------------
    # cache_plasma_weak/weak/nTOp_*.txt carry a fingerprint header recording the config
    # fields that affect their content; RecomputeWeakRates loads the cache
    # only if its fingerprint matches, and otherwise recomputes from scratch
    # (~2 s).  See weak_rates.RecomputeWeakRates for the full cache logic.
    #
    # Four additive correction terms control which physical effects enter the
    # total n<->p rate (mirroring PRIMAT-Main.m §IV.B):
    #
    #   radiative_corrections   -- True: replace the Born chi function with the
    #                              Coulomb + T=0 resummed radiative correction
    #                              (CCR, Phys. Rep. Eq. 101; Czarnecki et al. 2004).
    #                              False: use the bare Born approximation.
    #   finite_mass_corrections -- True: add the Fokker-Planck finite-nucleon-mass
    #                              correction (FMCCR if radiative_corrections=True,
    #                              FMNoCCR otherwise; Phys. Rep. §III.G).
    #   thermal_corrections     -- True: add the finite-temperature radiative
    #                              correction (CCRTh; Brown & Sawyer 2001,
    #                              Phys. Rep. §III.F, Eq. 108 = Eqs. 109+112+113,
    #                              plus the bremsstrahlung correction Eq. 107).
    #   spectral_distortions    -- (controlled in the neutrino section above)
    #                              Corrections from non-FD neutrino distributions;
    #                              internally uses SD_CCR or SD_Born depending on
    #                              radiative_corrections.
    #
    # Born (crude) mode = radiative_corrections=False, finite_mass_corrections=False,
    #                     thermal_corrections=False.  All True = full PRIMAT rate.
    #
    # Two combinations are ACCEPTED but are not self-consistent orders of the
    # same expansion; they are allowed on purpose (they isolate one term at a
    # time, which is what the flags are for) but should not be read as "PRIMAT
    # minus one effect":
    #
    #   thermal_corrections=True with radiative_corrections=False
    #       CCRTh is itself an O(alpha) radiative correction, and its integrands
    #       keep the Coulomb factor unconditionally
    #       (weak_rates.corrections._fermi_stat inside _ccrth_IPENCCRT and
    #       _ccrth_IPENCCRDiffBremsstrahlung), so the result is a Born base rate
    #       carrying a finite-temperature radiative correction *and* a Coulomb
    #       factor the base rate itself does not have.  For a clean Born rate,
    #       set thermal_corrections=False too.
    #
    #   spectral_distortions=True with analytic_distortions=False and
    #   finite_mass_corrections=True
    #       The finite-mass correction to the distortion channel (SD-FM,
    #       PRIMAT-Main-gray.m's deltaChiFM) needs closed-form energy
    #       derivatives of the distortion, which only the analytic (y_SZ/y_gray)
    #       distortion provides -- the tabulated NEVO distortion has none.  So
    #       in the DEFAULT (table) distortion mode the SD term is included at
    #       Born/CCR level but its finite-mass companion is silently absent.
    #       See weak_rates.corrections._correction_terms.
    "radiative_corrections":      True,  # True: Coulomb + T=0 resummed radiative corrections (CCR); False: Born approximation.
    "finite_mass_corrections":    True,  # True: add Fokker-Planck finite-nucleon-mass correction (FMCCR or FMNoCCR).
    "thermal_corrections":        True,  # True: add finite-temperature radiative corrections (CCRTh; Brown & Sawyer 2001).

    ##################### caching/saving options
    "cache_dir": None, # single writable directory for ALL regenerable caches (n<->p weak-rate nTOp_*.txt AND the plasma electron-thermo/QED tables); None (default) = <data_dir>/cache_plasma_weak/ inside the (possibly installed) package, with weak/ and plasma/ subdirs. Set it when the install location is read-only (e.g. system-wide site-packages): caches are then WRITTEN to <cache_dir>/{weak,plasma}/ (created on demand) and READ from there first, falling back to the shipped caches in the package (overlay -- shipped caches are never shadowed). Not part of any fingerprint -- the cache LOCATION cannot affect physics.
    "weak_rate_cache":            True,  # If False, never load the cache (always recompute); save_nTOp still controls whether the result is written back.
    "save_nTOp":                  True,  # If True, the computed n<->p rates are saved to cache_plasma_weak/weak/ (or the cache_dir redirect) as nTOp_<hash>.txt (forward and backward columns together).
    "sampling_nTOp_per_decade":   80,    # points per decade of T (T_end -> T_start) in the single n<->p rate grid

    "save_nTOp_thermal":          True,  # If True, the computed thermal n<->p rates are saved to cache_plasma_weak/weak/ (or the cache_dir redirect) as nTOp_thermal_<hash>.txt (both directions in one file).
    "sampling_nTOp_thermal_per_decade": 20,   # points per decade of T (T_end -> T_start) for the thermal-correction table
    ##################### Normalization of weak rates
    "tau_n_normalization":        True,  # Use neutron lifetime to normalize weak rates (instead of absolute normalization from GF, Vud, gA, etc.)
    "tau_n":                      878.4,  # neutron lifetime [s]; overrides the class-level constant when tau_n_normalization=True
    "std_tau_n":                  0.5,    # 1σ uncertainty on tau_n [s], used for MC sampling

    # Accuracy knobs for the thermal n<->p radiative correction integral, used
    # only when the thermal-correction cache must be recomputed (see
    # weak_rates._L_CCRTh_interpolants).  Evaluated with the `vegas`
    # Monte-Carlo library when available, else scipy.integrate.dblquad.
    "vegas_n_eval":               20000,   # vegas: evaluations per iteration
    "vegas_n_itn":                20,      # vegas: number of iterations
    "epsrel_thermal":             1.e-2,   # dblquad fallback: relative tolerance
    
    # ---- Output options ------------------------------------------------------
    # Writes a TSV (cfg.output_file) with the time evolution of T, t, and of
    # every nuclide's abundance in the chosen network (8/~59 for small/large,
    # fewer for large with an amax cutoff) plus the n<->p weak rates; see
    # nuclear_network.NuclearNetwork._write_time_evolution.
    "output_time_evolution":      False,
    "output_rates_time_evolution": False, #whether to append per-reaction forward-rate columns (<reaction>_frwrd, e.g. n_p__d_g_frwrd) to the time-evolution output, after the Y_<nuclide> block. One column per reaction in the active LT network (~12 for small/small_parthenope, 68 for large+amax=8, ~429 for full large). Only useful to inspect the rate evolution; keep False otherwise to save disk space. Both backends emit the identical columns.
    "output_n_points":            500, # number of time samples written to the time-evolution output (log-spaced over the integration range)
    "output_file":                "results/output_tables.tsv", # destination of the time-evolution TSV; None still fills results["evolution"] in memory but writes nothing to disk
    # Two-column dump (nuclide name, final mass-fraction abundance Y) at the end of BBN.
    "output_final_result":        False,
    "output_final_file":          "results/output_final.dat", # destination of that two-column final-abundance dump

    # Writes a separate TSV (cfg.output_background_file) with the cosmological
    # background's own time evolution (T, t, and -- if available -- a, H,
    # individual neutrino temperatures, NEVO heating function, and
    # plasma/neutrino/extra/total energy densities); see
    # background.Background.write_time_evolution.
    "output_background_evolution": False,
    "output_background_file":     "results/output_background.tsv", # destination of the background time-evolution TSV

    # Monte-Carlo output files (backend.run_mc/main.mc_uncertainty only; a
    # plain solve() never writes them). All three share one filename stem set
    # by output_mc_file_prefix, each gated by its own boolean below:
    #   <prefix>_samples.tsv      (output_mc_samples)     -- every raw sample,
    #                              one column per quantity, one row per sample
    #                              (see primat.backend.dump_mc_samples)
    #   <prefix>_covariance.tsv   (output_mc_covariance)  -- the (n_q, n_q)
    #                              sample covariance matrix (ddof=1) over all MC
    #                              quantities (primat.backend.dump_mc_covariance)
    #   <prefix>_correlation.tsv  (output_mc_correlation) -- the matching
    #                              correlation matrix (primat.backend.dump_mc_correlation)
    # The covariance/correlation give the *joint* nuclear-rate uncertainty (the
    # off-diagonal YP-D/H covariance a user needs for a joint likelihood), not
    # just the per-observable sigmas in the samples file. Each flag is opt-in
    # plumbing (default False); the prefix is used verbatim (its directory is
    # created on demand). E.g. prefix "results/output_mc" writes
    # results/output_mc_samples.tsv etc.
    "output_mc_samples":           False,
    "output_mc_covariance":        False, # write the (n_q, n_q) sample covariance matrix (ddof=1) to <output_mc_file_prefix>_covariance.tsv
    "output_mc_correlation":       False, # write the matching correlation matrix to <output_mc_file_prefix>_correlation.tsv
    "output_mc_file_prefix":       "results/output_mc", # shared filename stem of the three MC outputs above (its directory is created on demand)


    # ---- nuclear network --------------------------------------------------
    # Master grid onto which every nuclear reaction rate table is resampled at
    # load time.  This makes load_network grid-agnostic: tables generated with
    # different grids (e.g. via --keep-source-grid in convert_ac2024_rates.py,
    # or from external sources) are all resampled onto this common grid so that
    # fill_buffer's single searchsorted path remains valid.
    "rate_grid_npts":             1000,       # number of points in the master T9 grid
    "rate_grid_T9_min":          1.0e-3,     # minimum T9 [GK] on the master grid
    # Maximum T9 [GK] on the master grid.  Deliberately *below* the MT era's
    # start (T_weak = 1 MeV = 11.6 GK): the shipped tables' own source data
    # stops at T9 = 10, so rates above it are extrapolated off the last grid
    # cell either way.  Verified numerically inert (<= 2e-6 on every
    # observable) -- see load_network's grid comment in network_data.py.
    "rate_grid_T9_max":          10.0,

    # Network selector.  "small" is the built-in ORDER_SMALL network.  Any other
    # value loads data/nuclear/networks/<network>.txt -- shipped options are
    # "small_parthenope" and "large"; any other name loads a custom network
    # file of the same form.
    "network":                    "small",

    # Maximum nuclide mass number A = N + Z to include, for *any* network
    # (not just "large" -- a network whose nuclides are all below the cutoff
    # simply sees no reaction dropped). Reactions involving any nuclide with
    # A > amax are dropped. None = no filter (keep all reactions). Must be a
    # positive integer when set.
    "amax":                       None,

    # Absolute solve_ivp tolerance for the LT era of EVERY network (not just
    # "large" -- despite the legacy name). It must be tight enough for the
    # large network's heavy nuclides, which reach very small abundances.
    "atol_large_LT":              1.e-26,
    "rescale_nuclear_rates":            False, #Use to vary some rates with a uniform factor to explore their impact.

    # Cap applied to the MC rate rescaling factor during Monte Carlo runs.
    # When a p_* parameter is non-zero, the effective variation factor is  variation = sigma^p + delta
    # which can grow very large for extreme draws of p.  This parameter clamps
    # the variation to [1/cap, cap] before multiplying the median rate.
    # A value of 30 means no more than a factor of 30 up or down. Reactions carrying a flat "uncertainty factor
    # f=10-100" placeholder (e.g. CF88 rates such as He3_t__a_d/He3_t__a_n_p)
    # can otherwise draw a >=3-sigma p and multiply their rate by up to 1000x,
    # which for non-trace species (He3/t, unlike the many trace heavy-nuclide
    # branches sharing the same placeholder error) dominates the MC variance
    # of D/H with an unphysically large single-sample outlier rather than a
    # smooth uncertainty estimate. Set to None to disable the cap entirely.
    "mc_rate_rescale_cap":         30,

    # QED correction to select radiative-capture nuclear rates (Pitrou & Pospelov 2020).
    # Applies a T9-dependent multiplicative rescaling to the forward rate tables of
    # n_p__d_g, d_p__He3_g, t_p__a_g, t_a__Li7_g, He3_a__Be7_g at load time.  When True the
    # corrected values become the new medians, so p_* and delta_* variations
    # work relative to the QED-corrected central value.
    "nuclear_qed_corrections":    True,

    # ---- cosmological inputs ----------------------------------------------
    "Omegabh2":                   0.02242,   # baryon density Omega_b h^2 = 0.02242 +/- 0.00014 (Planck 2018 + BAO, author decision 2026-07-08)
    "Omegach2":                   0.11933,  # cold dark matter density parameter Omega_c h^2 (Planck 2018)
    "h":                          0.6766,   # reduced Hubble constant h = H_0 / (100 km/s/Mpc) (Planck 2018)
    "DeltaNeff":                  0., # extra relativistic degrees of freedom beyond the three SM neutrinos, in units of one neutrino species: adds DeltaNeff * rho_nu(one flavour) to the Friedmann equation and shifts Neff by the same amount (0 = Standard Model).
    "munuOverTnu":                0., # Reduced chemical potential xi = mu/T of neutrinos (the COMMON default for all flavours, nu_e, nu_mu, nu_tau; antineutrinos carry -xi).
    # A genuine chemical potential: it shifts the n<->p weak rates (FD_nu3 in the
    # rate integrands) AND raises the neutrino energy density / Neff by
    # rho(xi) = T^4 (7pi^2/120 + xi^2/4 + xi^4/(8 pi^2)) per flavour
    # (plasma.rho_nu_chempot_excess). It is part of the weak-rate cache fingerprint.
    # munuOverTnu != 0 with incomplete_decoupling=True is physically inconsistent (the NEVO table assumes it vanishes); use incomplete_decoupling=False to explore non-zero values.
    # Per-flavour overrides : only xi_e enters the n<->p weak rates (nu_e appears in n <-> p + e + nu_e), while all three gravitate through the. Each of the three below defaults to None,
    # meaning "inherit munuOverTnu"; set a float to give that flavour its own xi. 
    "munuOverTnu_e":              None, # Reduced chemical potential xi_e of nu_e; None = inherit munuOverTnu. Enters BOTH the n<->p weak rates and the energy density.
    "munuOverTnu_mu":             None, # Reduced chemical potential xi_mu of nu_mu; None = inherit munuOverTnu. Gravitates only (energy density / Neff), no weak-rate effect.
    "munuOverTnu_tau":            None, # Reduced chemical potential xi_tau of nu_tau; None = inherit munuOverTnu. Gravitates only (energy density / Neff), no weak-rate effect.

    # ---- Decay-era options -------------------------------------------------
    # decay_reverse_rates: when True, compute detailed-balance reverse rates
    # for radioactive-decay reactions, instead of treating them as irreversible
    # (abg = (0, 0, 0)).  During standard BBN the forward decay rate is
    # negligible (e.g. C14 T1/2 = 5700 yr ≫ t_end ≈ 10^6 s), so the reverse
    # rate is likewise negligible; enabling this only matters when T_end_MeV
    # is extended far below the standard 0.001 MeV and thermal equilibrium of
    # long-lived isotopes becomes relevant.
    "decay_reverse_rates":        False,

    # decay_era: if True and network="large", run a fourth "Decay Time" (DT)
    # integration era after the LT era, propagating abundances forward in time
    # (at fixed comoving scale) purely under radioactive decay (no Hubble
    # expansion, no thermal production).  The DT era spans t ∈ [t_end, t_end +
    # t_decay_end], log-spaced on decay_n_points time points.
    "decay_era":                  False,
    "t_decay_end":                3.156e16,  # DT era duration [s] (default: 1 Gyr = 3.156e16 s)
    "decay_n_points":             200,        # log-spaced output points in the DT era
    "output_decay_evolution":     False,      # write TSV of DT-era abundance time evolution
    "output_decay_file":          "results/output_decay_evolution.tsv", # destination of that DT-era TSV

    # ---- Early Dark Energy ------------------------------------------------
    "fEDE":                       0.,     # EDE fraction at peak; 0 = disabled
    "zcEDE":                      1.e8,   # redshift of EDE peak
    "wnEDE":                      1.,     # EDE equation-of-state parameter
}


# ---------------------------------------------------------------------------
# Machine-readable grouping of DEFAULT_PARAMS, mirroring the "# ----" section
# comments above.  This exists so the GUI's full-parameter listing, the CLI's
# ``--list-params``, and the param-template generator (``generate_rates/
# gen_param_templates.py``) can all derive the same section headings/order
# from a single place instead of three independently hand-maintained copies
# (the standing chore CLAUDE.md calls out).  Every DEFAULT_PARAMS key appears
# in exactly one group -- test_config.py checks that this stays exhaustive
# and non-overlapping whenever a key is added, removed, or renamed.
# ---------------------------------------------------------------------------
PARAM_GROUPS: dict = {
    "General behaviour and numerical settings": (
        "verbose", "debug", "show_progress", "numerical_precision",
        "numba_installed", "strict_params",
    ),
    "Neutrino decoupling": (
        "incomplete_decoupling",
    ),
    "Electromagnetic plasma": (
        "QED_corrections", "n_electron_table", "recompute_electron_thermo",
        "recompute_qed_corrections",
    ),
    "Spectral distortions": (
        "spectral_distortions", "analytic_distortions", "y_SZ", "y_gray",
    ),
    "Custom NEVO tables": (
        "nevo_file", "nevo_spectral_file", "nevo_grid_file", "nevo_file_prefix",
    ),
    "Data directory override and nuclear overlay": (
        "data_dir", "user_nuclear_dir",
    ),
    "Background mode": (
        "external_scale_factor", "custom_background",
    ),
    "Fundamental constants": (
        "GN",
    ),
    "Background thermodynamics": (
        "T_start_cosmo_MeV", "T_end_MeV", "sampling_temperature_per_decade",
    ),
    "n <-> p weak rates": (
        "radiative_corrections", "finite_mass_corrections", "thermal_corrections",
    ),
    "Caching/saving options": (
        "cache_dir", "weak_rate_cache", "save_nTOp", "sampling_nTOp_per_decade",
        "save_nTOp_thermal", "sampling_nTOp_thermal_per_decade",
    ),
    "Normalization of weak rates": (
        "tau_n_normalization", "tau_n", "std_tau_n",
    ),
    "Thermal correction accuracy knobs": (
        "vegas_n_eval", "vegas_n_itn", "epsrel_thermal",
    ),
    "Output options": (
        "output_time_evolution", "output_rates_time_evolution",
        "output_n_points", "output_file", "output_final_result",
        "output_final_file", "output_background_evolution",
        "output_background_file", "output_mc_samples", "output_mc_covariance",
        "output_mc_correlation", "output_mc_file_prefix",
    ),
    "Nuclear network": (
        "rate_grid_npts", "rate_grid_T9_min", "rate_grid_T9_max", "network",
        "amax", "atol_large_LT", "rescale_nuclear_rates",
        "mc_rate_rescale_cap", "nuclear_qed_corrections",
    ),
    "Cosmological inputs": (
        "Omegabh2", "Omegach2", "h", "DeltaNeff", "munuOverTnu",
        "munuOverTnu_e", "munuOverTnu_mu", "munuOverTnu_tau",
    ),
    "Decay-era options": (
        "decay_reverse_rates", "decay_era", "t_decay_end", "decay_n_points",
        "output_decay_evolution", "output_decay_file",
    ),
    "Early Dark Energy": (
        "fEDE", "zcEDE", "wnEDE",
    ),
}


# String-valued config keys that represent filesystem paths.
# These are normalized with os.path.expanduser() so CLI users can pass
# quoted "~/" prefixes through --set and still get the expected home-dir
# expansion.  Must stay in step with primat-c/src/config.c's
# cpr_is_path_field(): a key expanded on one side but not the other makes the
# two backends read/write different directories from one config (cache_dir was
# missing here while C expanded it, so the same config sent the two backends'
# caches to different directories).
_PATH_PARAMS = {
    "nevo_file",
    "nevo_spectral_file",
    "nevo_grid_file",
    "custom_background",
    "data_dir",
    "user_nuclear_dir",
    "cache_dir",
    "output_file",
    "output_final_file",
    "output_background_file",
    "output_mc_file_prefix",
    "output_decay_file",
}


def _expanduser_path(value):
    """Expand a user-home prefix in a path-like config value.

    Parameters
    ----------
    value : str | os.PathLike | None
        Raw path value supplied by the caller. ``None`` is passed through
        unchanged so optional path parameters keep their sentinel value.

    Returns
    -------
    str | None
        The same path with a leading ``~`` resolved against the current
        user home directory, or ``None`` if that was the input.

    Example
    -------
        >>> _expanduser_path("~/Downloads/custom")
        '/home/user/Downloads/custom'
    """
    if value is None:
        return None
    return os.path.expanduser(os.fspath(value))


def _rates_overlay_notice(field: str, path: str) -> str:
    """Render the startup note for a custom data/nuclear overlay directory.

    Parameters
    ----------
    field : str
        Either ``"data_dir"`` (full-takeover data root) or
        ``"user_nuclear_dir"`` (additive nuclear overlay).
    path : str        Directory path already accepted by the config validator.

    Returns
    -------
    str
        Human-readable notice explaining the effect of the override.
    """
    if field == "data_dir":
        label = "full-takeover data directory"
        detail = "entire data tree (NEVO/, nuclear/, csv/, cache_plasma_weak/) replaced"
    else:
        label = "additive nuclear overlay"
        detail = "nuclear networks and rate tables"
    return (
        f"[init]  {field} {label} override: {detail} under "
        f"{os.path.abspath(os.path.expanduser(os.fspath(path)))!r}."
    )


def _overlay_candidates(base: str, relpath: str) -> list[str]:
    """Return overlay lookup candidates for a rates-relative path.

    The shipped tree is rooted at ``primat/data`` and therefore uses paths
    such as ``nuclear/networks/small.txt``.  Overlay directories are treated
    as the equivalent of ``primat/data/nuclear`` instead, so the primary
    lookup drops a leading ``nuclear/`` component when present and then
    falls back to the legacy nested layout for compatibility.
    """
    candidates = []
    # ``relpath`` is built by callers with os.path.join, so on Windows its
    # components are separated by "\\", not "/". Normalise to forward slashes
    # before detecting a leading "nuclear/" component, otherwise the
    # nuclear-stripped overlay candidate is never generated on Windows and a
    # user_nuclear_dir overlay's networks/rate tables become invisible there.
    norm = relpath.replace(os.sep, "/")
    if os.altsep:
        norm = norm.replace(os.altsep, "/")
    if norm.startswith("nuclear/"):
        stripped_parts = norm[len("nuclear/"):].split("/")
        candidates.append(os.path.join(base, *stripped_parts))
    candidates.append(os.path.join(base, relpath))
    return candidates

def _is_decorative_comment(text: str) -> bool:
    """True if a comment line is a section rule/heading rather than prose.

    ``DEFAULT_PARAMS`` is divided by decorative separators
    (``# ---- Output options ----------``,
    ``##################### caching/saving options``, bare ``# -----`` rules).
    Those describe a whole *group* of keys, so they must never be picked up as
    the description of whichever key happens to follow them -- see
    :func:`_default_params_comments`'s block fallback.

    Parameters
    ----------
    text : str
        One raw source line, already stripped of surrounding whitespace but
        *with* its leading ``#`` intact (the ``#####`` runs are only
        recognisable before stripping).

    Returns
    -------
    bool
        True for separators/headings, False for ordinary explanatory prose.

    Example
    -------
        >>> _is_decorative_comment("# ---- Output options ----")
        True
        >>> _is_decorative_comment("# Newton's constant, SI units")
        False
    """
    if text.startswith("####"):
        return True
    body = text.lstrip("#").strip()
    if not body:
        return True                      # a bare "#" spacer line
    if body.startswith("----") or body.startswith("===="):
        return True                      # "# ---- Section name ----"
    return not body.strip("-=# ")        # a pure rule, e.g. "# --------"


def _default_params_comments():
    """Parse this file's own source to extract a one-line description of every
    ``DEFAULT_PARAMS`` key (CLI discoverability: ``primat --list-params`` uses
    this to show users what every parameter means without duplicating the
    explanation in a second place that could drift out of sync).

    Two sources, in order of preference:

    1. The ``# ...`` text on the *same source line* as the key's value.
    2. Failing that, the **first sentence of the contiguous comment block
       immediately above the key**, with decorative section rules/headings
       skipped (:func:`_is_decorative_comment`).  A dozen keys are documented
       that way -- ``network``, ``amax``, ``atol_large_LT``, the ``output_*``
       group, the ``decay_*`` group -- and used to print with no description at
       all, even though ``--list-params`` is the documented way to discover
       parameters (``--set`` being hidden from ``--help``).

    Only the block's first *sentence* is taken, not the whole block: the rest is
    multi-paragraph prose for a human reading the dict, not a one-liner (a full
    explanation belongs in the source, which is where a curious user is
    pointed).

    Returns
    -------
    dict
        Maps each ``DEFAULT_PARAMS`` key to its one-line description
        (``""`` only if the key has neither a trailing comment nor a
        preceding comment block -- ``test_config.py`` asserts that never
        happens).

    Example
    -------
        >>> _default_params_comments()["network"]     # doctest: +SKIP
        'Network selector.  "small" is the built-in ORDER_SMALL network.'
    """
    import ast
    import tokenize

    path = __file__
    # Explicit UTF-8: config.py contains non-ASCII physics characters in its
    # comments (ν, ↔, →, …); the default locale encoding is cp1252 on Windows
    # and would raise UnicodeDecodeError re-reading this very file.
    with open(path, encoding="utf-8") as f:
        source = f.read()
    tree = ast.parse(source, filename=path)

    dict_node = None
    for node in ast.walk(tree):
        # DEFAULT_PARAMS is declared with a type annotation ("DEFAULT_PARAMS:
        # dict = {...}"), i.e. an ast.AnnAssign, not a plain ast.Assign.
        if (isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
                and node.target.id == "DEFAULT_PARAMS"):
            dict_node = node.value
            break
    if dict_node is None:
        return {}

    # Map each source line number to the comment token that starts on it, by
    # re-tokenizing the same file (comments aren't part of the AST).
    line_comments = {}
    # Explicit UTF-8: config.py contains non-ASCII physics characters in its
    # comments (ν, ↔, →, …); the default locale encoding is cp1252 on Windows
    # and would raise UnicodeDecodeError re-reading this very file.
    with open(path, encoding="utf-8") as f:
        for tok in tokenize.generate_tokens(f.readline):
            if tok.type == tokenize.COMMENT:
                line_comments[tok.start[0]] = tok.string.lstrip("#").strip()

    # Raw source lines (1-based indexing below, to match ast line numbers), so
    # the block fallback can tell a comment-ONLY line from a trailing comment.
    src_lines = source.splitlines()

    def _preceding_block_summary(lineno: int) -> str:
        """First sentence of the contiguous comment block ending just above
        source line ``lineno`` (1-based), or ``""`` if there is none.

        The block's lines are joined before cutting, so the result is a whole
        sentence rather than however much of it fitted on the first physical
        line ("Writes a TSV (cfg.output_file) with the time evolution of T, t,
        and of" was the pre-join behaviour).
        """
        block = []
        i = lineno - 1                       # 1-based line above the key
        while i >= 1:
            stripped = src_lines[i - 1].strip()
            if not stripped.startswith("#"):
                break
            block.append(stripped)
            i -= 1
        block.reverse()                      # top-to-bottom reading order
        prose = [raw.lstrip("#").strip() for raw in block
                 if not _is_decorative_comment(raw)]
        if not prose:
            return ""
        text = " ".join(prose)
        # Cut at the first sentence end. A ". " inside an identifier or an
        # abbreviation ("cfg.output_file", "e.g.") would truncate absurdly
        # early, so only accept a boundary once enough text has accumulated to
        # be a plausible sentence -- and require the next character to start a
        # new one (whitespace, then a capital or a quote/paren).
        _MIN_SENTENCE = 50      # chars; shorter "sentences" are abbreviations
        for match in re.finditer(r"\.(\s+|$)", text):
            end = match.start() + 1
            if end < _MIN_SENTENCE:
                continue
            rest = text[match.end():match.end() + 1]
            if rest and not (rest.isupper() or rest in "\"'(*`"):
                continue        # e.g. "... at 1e-3. see below" -> keep going
            return text[:end]
        # No sentence end at all (a bullet-style block): keep it to one line.
        _MAX_LEN = 200
        if len(text) <= _MAX_LEN:
            return text
        return text[:text.rfind(" ", 0, _MAX_LEN)] + " ..."

    comments = {}
    for key_node, value_node in zip(dict_node.keys, dict_node.values):
        key = ast.literal_eval(key_node)
        text = line_comments.get(value_node.end_lineno, "")
        if not text:
            text = _preceding_block_summary(key_node.lineno)
        comments[key] = text
    return comments


# ===========================================================================
# Parameter validation 
# ---------------------------------------------------------------------------
# A physicist's first contact with primat is often a typo in a params dict or a
# YAML file.  Without validation a wrong *type* surfaces much later as an opaque
# stack trace from deep inside the thermodynamics (e.g. Omegabh2="0.022" dying
# with "can't multiply sequence by non-int"), and a wrong *key* (Omegab2h)
# silently runs the default cosmology.  The helpers below validate every
# user-supplied DEFAULT_PARAMS override at construction time, raising a
# one-line, self-explanatory TypeError/ValueError that names the key, the
# received value, and the expected type/range.  The C backend mirrors the type
# checks in cpr_config_set_by_name and the range checks in cpr_config_validate
# (primat-c/src/config.c), per the CLAUDE.md parity mandate.
# ===========================================================================

# Kind tags used by the spec below.  Each maps to a predicate on a candidate
# value.  Numeric kinds use the ``numbers`` ABCs so numpy scalars (np.float64,
# np.int64 -- common in MCMC drivers) are accepted, while ``bool`` (a subclass
# of int) is explicitly excluded from the numeric kinds: passing True where a
# float is expected is a bug, not the number 1.0.
_KIND_CHECKS = {
    "bool":  lambda v: isinstance(v, bool),
    "int":   lambda v: isinstance(v, numbers.Integral) and not isinstance(v, bool),
    "float": lambda v: isinstance(v, numbers.Real) and not isinstance(v, bool),
    "str":   lambda v: isinstance(v, str),
    "none":  lambda v: v is None,
}
_KIND_ENGLISH = {"bool": "bool", "int": "int", "float": "float",
                 "str": "str", "none": "None"}

# Keys whose accepted kinds cannot be inferred from ``DEFAULT_PARAMS`` (their
# default is ``None``, so it carries no type) or that accept ``None`` in
# addition to their default's type (path parameters where ``None`` is a "skip"
# sentinel, and mc_rate_rescale_cap where ``None`` disables the cap).  Every
# other key infers its single accepted kind from the default value's type
# (see :func:`_param_kinds`).
_PARAM_TYPESPEC = {
    "amax":                  ("int", "none"),   # None = no A cutoff; else positive int
    "mc_rate_rescale_cap":   ("float", "none"),  # None = no cap; else positive number
    # Per-flavour neutrino chemical potentials: None = inherit munuOverTnu; else
    # any real xi (may be negative). See xi_nu_e/xi_nu_mu/xi_nu_tau properties.
    "munuOverTnu_e":         ("float", "none"),
    "munuOverTnu_mu":        ("float", "none"),
    "munuOverTnu_tau":       ("float", "none"),
    # Optional filesystem-path parameters: a str path or None ("use default" /
    # "skip this output").  Mirrors _PATH_PARAMS above.
    "nevo_file":             ("str", "none"),
    "nevo_spectral_file":    ("str", "none"),
    "nevo_grid_file":        ("str", "none"),
    "custom_background":      ("str", "none"),
    "data_dir":              ("str", "none"),
    "user_nuclear_dir":      ("str", "none"),
    "cache_dir":             ("str", "none"),
    "output_file":           ("str", "none"),
    "output_final_file":     ("str", "none"),
    "output_background_file": ("str", "none"),
    "output_mc_file_prefix": ("str", "none"),
    "output_decay_file":     ("str", "none"),
}

# Numeric range constraints, where the physics/numerics demand them.  Each
# entry is ``(predicate, human_text)``; the predicate is applied only to
# non-None values that already passed the type check.  fEDE (0<=fEDE<1) and
# amax (>=1) are validated by their own dedicated methods (_validate_fEDE /
# _validate_amax) with bespoke messages, so they are intentionally absent here.
_POSITIVE = (lambda v: v > 0, "must be > 0")
_POSITIVE_INT = (lambda v: v >= 1, "must be a positive integer (>= 1)")
_NON_NEGATIVE = (lambda v: v >= 0, "must be >= 0")
_PARAM_RANGE = {
    # ---- strictly positive floats (a physical scale, tolerance, or time) ----
    "numerical_precision": _POSITIVE,
    "GN":                  _POSITIVE,
    "T_start_cosmo_MeV":   _POSITIVE,
    "T_end_MeV":           _POSITIVE,
    "tau_n":               _POSITIVE,
    "Omegabh2":            _POSITIVE,
    "Omegach2":            _POSITIVE,
    "h":                   _POSITIVE,
    "atol_large_LT":       _POSITIVE,
    "epsrel_thermal":      _POSITIVE,
    "t_decay_end":         _POSITIVE,
    "zcEDE":               _POSITIVE,
    "rate_grid_T9_min":    _POSITIVE,
    "rate_grid_T9_max":    _POSITIVE,
    "mc_rate_rescale_cap": _POSITIVE,
    # ---- strictly positive counts (grid points, iterations, samplings) ------
    "n_electron_table":                  _POSITIVE_INT,
    "sampling_temperature_per_decade":   _POSITIVE_INT,
    "sampling_nTOp_per_decade":          _POSITIVE_INT,
    "sampling_nTOp_thermal_per_decade":  _POSITIVE_INT,
    "vegas_n_eval":                      _POSITIVE_INT,
    "vegas_n_itn":                       _POSITIVE_INT,
    "output_n_points":                   _POSITIVE_INT,
    "rate_grid_npts":                    _POSITIVE_INT,
    "decay_n_points":                    _POSITIVE_INT,
    # ---- non-negative (a 1-sigma width may legitimately be 0) ---------------
    "std_tau_n":           _NON_NEGATIVE,
}


def _param_kinds(key: str):
    """Return the tuple of accepted kind tags for a ``DEFAULT_PARAMS`` key.

    Uses the explicit :data:`_PARAM_TYPESPEC` entry when present (None-defaulted
    or None-able keys); otherwise infers a single kind from the default value's
    Python type.  ``bool`` is checked before ``int`` because ``bool`` is a
    subclass of ``int``.
    """
    if key in _PARAM_TYPESPEC:
        return _PARAM_TYPESPEC[key]
    default = DEFAULT_PARAMS[key]
    if isinstance(default, bool):
        return ("bool",)
    if isinstance(default, int):
        return ("int",)
    if isinstance(default, float):
        return ("float",)
    if isinstance(default, str):
        return ("str",)
    return None  # unreachable: every None-defaulted key is in _PARAM_TYPESPEC


_KIND_PYTYPE = {"bool": "bool", "int": "int", "float": "float",
                 "str": "str", "none": "None"}

# DEFAULT_PARAMS keys that must be skipped by the generator below because
# PRIMATConfig already declares them as a real, independently-typed
# ``@property`` further down in the class body (see Omegabh2's getter/
# setter): a bare TYPE_CHECKING annotation of the same name would be a
# duplicate definition as far as mypy is concerned.
_CONFIG_ANNOTATION_SKIP = {"Omegabh2"}


def _generate_config_type_annotations() -> str:
    """Generate the body of the ``if TYPE_CHECKING:`` attribute-annotation
    block inside :class:`PRIMATConfig` (between the ``BEGIN``/``END
    GENERATED PARAM ANNOTATIONS`` sentinel comments below).

    ``PRIMATConfig`` sets every ``DEFAULT_PARAMS`` key as a plain instance
    attribute at construction time (see ``_apply_user_overrides``), but does
    so dynamically from a loop over the dict, so neither an IDE nor mypy can
    see ``cfg.Omegabh2`` or ``cfg.network`` coming -- worse, ``__getattr__``
    (used for the unrelated ``p_<rxn>``/``delta_<rxn>`` dynamic pattern)
    makes static tools treat *any* misspelled attribute access as valid
    rather than flagging it. Bare ``name: type`` class-body annotations
    (no assignment) fix both: they cost nothing at runtime (``TYPE_CHECKING``
    is always ``False``) and are exactly what IDEs/mypy use for instance
    attribute completion and checking. Generated (rather than hand-typed) so
    it cannot drift from ``DEFAULT_PARAMS``/``_PARAM_TYPESPEC``;
    ``test_config.py`` fails if the block in the source is stale --
    regenerate it by running this function and pasting its output between
    the sentinels, or via ``python -m primat.tools.gen_docs --check``.
    """
    lines = []
    for key in DEFAULT_PARAMS:
        if key in _CONFIG_ANNOTATION_SKIP:
            continue
        kinds = _param_kinds(key)
        pytype = " | ".join(_KIND_PYTYPE[k] for k in kinds)
        lines.append(f"        {key}: {pytype}")
    return "\n".join(lines)


def _validate_param_value(key: str, value):
    """Type-, choice-, and range-check one user-supplied ``DEFAULT_PARAMS``
    override, raising an immediate, self-explanatory error on any mismatch.

    Parameters
    ----------
    key : str
        A key known to be in ``DEFAULT_PARAMS`` (p_<rxn>/delta_<rxn> dynamic
        keys are handled separately and never reach here).
    value : object
        The candidate override value.

    Raises
    ------
    TypeError
        If ``value``'s type does not match any accepted kind for ``key``.
    ValueError
        If ``value`` is outside the allowed numeric range
        (:data:`_PARAM_RANGE`).

    Example
    -------
        >>> _validate_param_value("Omegabh2", "0.022")   # doctest: +SKIP
        TypeError: PRIMATConfig: parameter 'Omegabh2' got '0.022' of type str;
        expected float.
        >>> _validate_param_value("Omegabh2", -0.1)      # doctest: +SKIP
        ValueError: PRIMATConfig: parameter 'Omegabh2' got -0.1, which is out
        of range: must be > 0.
    """
    kinds = _param_kinds(key)
    if kinds is None:
        return
    if not any(_KIND_CHECKS[k](value) for k in kinds):
        expected = " or ".join(_KIND_ENGLISH[k] for k in kinds)
        raise TypeError(
            f"PRIMATConfig: parameter {key!r} got {value!r} of type "
            f"{type(value).__name__}; expected {expected}."
        )
    if value is not None and key in _PARAM_RANGE:
        predicate, text = _PARAM_RANGE[key]
        if not predicate(value):
            raise ValueError(
                f"PRIMATConfig: parameter {key!r} got {value!r}, "
                f"which is out of range: {text}."
            )


class PRIMATConfig:
    """
    Immutable physical constants + mutable run-time parameters.

    Usage::

        cfg = PRIMATConfig()                    # all defaults
        cfg = PRIMATConfig({"Omegabh2": 0.022, "network": "large", "amax": 8})

    After construction every key in ``DEFAULT_PARAMS`` is an attribute, plus
    all physical constants listed below.
    """

    if TYPE_CHECKING:
        # BEGIN GENERATED PARAM ANNOTATIONS -- do not edit by hand; produced
        # by _generate_config_type_annotations() in this module (see its
        # docstring). test_config.py fails if this block goes stale.
        verbose: bool
        debug: bool
        show_progress: bool
        numerical_precision: float
        numba_installed: bool
        strict_params: bool
        incomplete_decoupling: bool
        QED_corrections: bool
        n_electron_table: int
        recompute_electron_thermo: bool
        recompute_qed_corrections: bool
        spectral_distortions: bool
        analytic_distortions: bool
        y_SZ: float
        y_gray: float
        nevo_file: str | None
        nevo_spectral_file: str | None
        nevo_grid_file: str | None
        nevo_file_prefix: str
        data_dir: str | None
        user_nuclear_dir: str | None
        external_scale_factor: bool
        custom_background: str | None
        GN: float
        T_start_cosmo_MeV: float
        T_end_MeV: float
        sampling_temperature_per_decade: int
        radiative_corrections: bool
        finite_mass_corrections: bool
        thermal_corrections: bool
        cache_dir: str | None
        weak_rate_cache: bool
        save_nTOp: bool
        sampling_nTOp_per_decade: int
        save_nTOp_thermal: bool
        sampling_nTOp_thermal_per_decade: int
        tau_n_normalization: bool
        tau_n: float
        std_tau_n: float
        vegas_n_eval: int
        vegas_n_itn: int
        epsrel_thermal: float
        output_time_evolution: bool
        output_rates_time_evolution: bool
        output_n_points: int
        output_file: str | None
        output_final_result: bool
        output_final_file: str | None
        output_background_evolution: bool
        output_background_file: str | None
        output_mc_samples: bool
        output_mc_covariance: bool
        output_mc_correlation: bool
        output_mc_file_prefix: str | None
        rate_grid_npts: int
        rate_grid_T9_min: float
        rate_grid_T9_max: float
        network: str
        amax: int | None
        atol_large_LT: float
        rescale_nuclear_rates: bool
        mc_rate_rescale_cap: float | None
        nuclear_qed_corrections: bool
        Omegach2: float
        h: float
        DeltaNeff: float
        munuOverTnu: float
        munuOverTnu_e: float | None
        munuOverTnu_mu: float | None
        munuOverTnu_tau: float | None
        decay_reverse_rates: bool
        decay_era: bool
        t_decay_end: float
        decay_n_points: int
        output_decay_evolution: bool
        output_decay_file: str | None
        fEDE: float
        zcEDE: float
        wnEDE: float
        # Omegabh2 is intentionally absent here: it is a real @property
        # below with its own getter/setter type annotations.
        # p_<reaction>/delta_<reaction> nuclear-rate-variation keys are a
        # dynamic, unbounded pattern (any reaction name) routed through
        # __getattr__/__setattr__ further below; not enumerable here.
        # END GENERATED PARAM ANNOTATIONS

    @property
    def is_small(self) -> bool:
        """True if using the 'small' network."""
        return self.network == "small"

    @property
    def is_large(self) -> bool:
        """True if using the 'large' network."""
        return self.network == "large"

    # ------------------------------------------------------------------
    # Effective per-flavour neutrino chemical potentials
    # ------------------------------------------------------------------
    # ``munuOverTnu_e/mu/tau`` default to None ("inherit the common
    # ``munuOverTnu``"); these three properties resolve None -> munuOverTnu so
    # the rest of the code (weak rates, Plasma.rho_nu, the Friedmann sum) always
    # reads a concrete float. Only ``xi_nu_e`` gates the n<->p weak rates (nu_e
    # is the flavour that appears in n <-> p + e + nu_e); all three enter the
    # neutrino energy density / Neff. When all three are left at their default
    # None with a single ``munuOverTnu``, xi_nu_e == xi_nu_mu == xi_nu_tau ==
    # munuOverTnu, reproducing the previous single-xi behaviour bit-for-bit.
    @property
    def xi_nu_e(self) -> float:
        """Effective reduced chemical potential xi_e = mu_{nu_e}/T of nu_e.

        Returns ``munuOverTnu_e`` if the user set it, else falls back to the
        common ``munuOverTnu``. This is the ONLY flavour that shifts the n<->p
        weak rates (via the FD_nu3 integrand in weak_rates), so it is also the
        one carried in the weak-rate cache fingerprint.
        """
        return self.munuOverTnu if self.munuOverTnu_e is None else self.munuOverTnu_e

    @property
    def xi_nu_mu(self) -> float:
        """Effective reduced chemical potential xi_mu of nu_mu (``munuOverTnu_mu``
        or, if None, ``munuOverTnu``). nu_mu only gravitates: it enters the
        neutrino energy density / Neff (Plasma.rho_nu) but not the weak rates.
        """
        return self.munuOverTnu if self.munuOverTnu_mu is None else self.munuOverTnu_mu

    @property
    def xi_nu_tau(self) -> float:
        """Effective reduced chemical potential xi_tau of nu_tau (``munuOverTnu_tau``
        or, if None, ``munuOverTnu``). Like nu_mu, nu_tau only gravitates
        (energy density / Neff), with no n<->p weak-rate effect.
        """
        return self.munuOverTnu if self.munuOverTnu_tau is None else self.munuOverTnu_tau

    # ------------------------------------------------------------------
    # Physical constants and unit-conversion factors
    # ------------------------------------------------------------------
    # All fixed PDG values, CGS<->natural-units conversion factors, and the
    # purely-constant derived quantities (sW2, s0bar, s0CMB, n0CMB, mB,
    # HubbleOverh, the fixed temperature eras T_start/T_weak/T_nucl/T_end,
    # ...) live in primat.constants.Constants (see that module for
    # definitions, formulas and citations). They are re-exposed here as
    # plain class attributes so existing code (cfg.me, cfg.MeV_to_Kelvin,
    # cfg.s0bar, ...) is unaffected; new physics code may instead import
    # CONST directly from primat.constants.
    Kelvin         = CONST.Kelvin
    second         = CONST.second
    cm             = CONST.cm
    gram           = CONST.gram
    erg            = CONST.erg
    kB             = CONST.kB
    clight         = CONST.clight
    hbar           = CONST.hbar
    Mpc            = CONST.Mpc
    MeV            = CONST.MeV
    keV            = CONST.keV
    alphaem        = CONST.alphaem
    GF             = CONST.GF
    mZ             = CONST.mZ
    me             = CONST.me
    mn             = CONST.mn
    mp             = CONST.mp
    T0CMB          = CONST.T0CMB
    MeV_to_Kelvin  = CONST.MeV_to_Kelvin
    MeV_to_secm1   = CONST.MeV_to_secm1
    MeV_to_g       = CONST.MeV_to_g
    MeV_to_cmm1    = CONST.MeV_to_cmm1
    MeV4_to_gcmm3  = CONST.MeV4_to_gcmm3
    T_start        = CONST.T_start
    T_weak         = CONST.T_weak
    T_nucl         = CONST.T_nucl
    sW2            = CONST.sW2
    geL            = CONST.geL
    geR            = CONST.geR
    gmuL           = CONST.gmuL
    gmuR           = CONST.gmuR
    gA             = CONST.gA
    kappa_p        = CONST.kappa_p
    kappa_n        = CONST.kappa_n
    deltakappa     = CONST.deltakappa
    Vud            = CONST.Vud
    radproton      = CONST.radproton
    s0bar          = CONST.s0bar
    s0CMB          = CONST.s0CMB
    n0CMB          = CONST.n0CMB
    ma             = CONST.ma
    He4Overma      = CONST.He4Overma
    HOverma        = CONST.HOverma
    Neff_SM        = CONST.Neff_SM
    mB             = CONST.mB
    maOvermB       = CONST.maOvermB
    HubbleOverh    = CONST.HubbleOverh

    # ------------------------------------------------------------------
    # Quantities depending on overridable parameters (GN, T_start_cosmo_MeV)
    # ------------------------------------------------------------------

    # Temperature era set by the overridable T_start_cosmo_MeV [K].
    @property
    def T_start_cosmo(self) -> float:
        return self.T_start_cosmo_MeV * self.MeV_to_Kelvin

    @property
    def T_end(self) -> float:
        """End temperature for nuclear integration [K].

        Set via ``T_end_MeV`` [MeV] in ``DEFAULT_PARAMS`` (default 0.001 MeV,
        i.e. the standard BBN endpoint at ~11.6 MK / ~1.3×10^6 s).  Making
        it configurable allows extending the integration into the Decay Time
        era (``decay_era=True``) or performing custom post-BBN analysis at
        lower temperatures.

        The default 0.001 MeV (≈ 11.6 MK, cosmic time ≈ 1.3×10⁶ s ≈ 15 days)
        is the standard end point of BBN integration.

        Example::

            # Extend BBN integration to 0.0001 MeV (10× lower than default):
            cfg = PRIMATConfig({"T_end_MeV": 1e-4})
        """
        return self.T_end_MeV * self.MeV_to_Kelvin

    # Gravity: GN is overridable, so it lives in DEFAULT_PARAMS only.
    # tau_n [s] is similarly overridable (DEFAULT_PARAMS), used by weak_rates.
    #
    # cfg.GN is stored in SI units [m^3 kg^-1 s^-2] (so it reads/edits like any
    # textbook value of Newton's constant), but the Friedmann equation below is
    # written in the natural-units (hbar=c=1) convention used throughout the
    # rest of the code, where G has dimension [energy]^-2. Convert once here via
    # CONST.GN_SI_to_MeV2 (see that property's docstring for the derivation).
    @property
    def _GN_MeV2(self) -> float:
        """Newton's constant in natural units [MeV^-2], converted from the
        SI-valued ``self.GN``."""
        return self.GN * CONST.GN_SI_to_MeV2

    @property
    def Mpl(self) -> float:
        return 1. / np.sqrt(self._GN_MeV2)

    @property
    def rhocOverh2(self) -> float:
        return 3. / (8. * np.pi * self._GN_MeV2) * self.HubbleOverh**2  # [MeV^4/h^2]

    # ------------------------------------------------------------------
    # Constructor: merge user params over defaults
    # ------------------------------------------------------------------
    def __init__(self, params: dict | None = None):
        self._init_defaults_and_nuclide_data(params)
        user_keys = self._apply_user_overrides(params)
        self._validate_fEDE()
        self._validate_custom_background()
        self._validate_data_dirs()
        # amax before _setup_rate_variation_defaults, which compares
        # reaction_category(...) > self.amax and so must not see an invalid one.
        self._validate_amax()
        self._validate_ranges()
        self._validate_network()
        valid_rxns = self._setup_rate_variation_defaults()
        self._warn_unknown_rate_variations(user_keys, valid_rxns)
        self._detect_optional_libraries()
        self._validate_nevo_files()
        self._validate_physics_flag_combos()

        # Derived cosmological quantity (depends on Omegabh2)
        self._update_derived()

    def _init_defaults_and_nuclide_data(self, params: dict | None):
        """Seed every ``DEFAULT_PARAMS`` key as an instance attribute, apply
        an early ``data_dir`` override, and load the nuclide CSV data.

        Bypasses ``__setattr__`` for the initial dict setup to avoid
        interference before the base dicts (``p_rxn``/``delta_rxn``) even
        exist. ``data_dir`` is applied -- and *validated* -- before
        ``_load_nuclide_data`` so ``nuclides.csv`` is read from the
        user-supplied root when one is provided, and so a typo'd root is
        reported as such instead of as a bare ``FileNotFoundError`` on
        ``<data_dir>/csv/nuclides.csv``.
        """
        for key, value in DEFAULT_PARAMS.items():
            # Deep copy dictionaries to avoid shared state between instances
            if isinstance(value, dict):
                object.__setattr__(self, key, value.copy())
            else:
                object.__setattr__(self, key, value)

        # Initialize nuclear rate variation dicts as empty for now.  They are
        # populated with the configured network's reactions *after* user
        # overrides are applied (self.network may itself be one of those
        # overrides), so that the per-reaction defaults match the network
        # actually requested by the caller -- see _setup_rate_variation_defaults.
        object.__setattr__(self, "p_rxn", {})
        object.__setattr__(self, "delta_rxn", {})
        # Created before the data_dir check below, which appends to it.
        object.__setattr__(self, "_init_messages", [])

        if params and "data_dir" in params:
            _validate_param_value("data_dir", params["data_dir"])
            object.__setattr__(self, "data_dir", _expanduser_path(params["data_dir"]))
            self._validate_dir_field("data_dir")

        self._load_nuclide_data()

    def _apply_user_overrides(self, params: dict | None) -> set:
        """Apply caller-supplied ``params`` on top of the defaults.

        Any key in ``DEFAULT_PARAMS``, or starting with the ``p_``/``delta_``
        rate-variation prefixes, is set via ``setattr`` (routed through the
        class's ``__setattr__``); anything else triggers an "unknown
        parameter" warning. Returns the raw set of user-supplied keys (used
        later by :meth:`_warn_unknown_rate_variations` to catch p_*/delta_*
        typos against the *finalised* network's reaction list).
        """
        user_keys = set(params.keys()) if params else set()
        if not params:
            return user_keys

        # strict_params governs how unknown keys are handled (warn vs. raise);
        # read its effective value up front -- from this very override if the
        # caller set it, else the default already seeded on self -- so the
        # decision is available while iterating.
        strict = bool(params.get("strict_params", self.strict_params))

        known_prefixes = ('p_', 'delta_')
        unknown = []
        for key, value in params.items():
            if key in DEFAULT_PARAMS:
                # Type/choice/range-check the value *before* storing it, so a
                # bad override never reaches the solver as a confusing later
                # error (see _validate_param_value's docstring).
                _validate_param_value(key, value)
                setattr(self, key, value)
            elif any(key.startswith(p) for p in known_prefixes):
                # p_<rxn>/delta_<rxn> rate variations: routed to the backing
                # dicts by __setattr__ (which float()s the value, raising on a
                # non-numeric one); validated against the network's reaction
                # list later by _warn_unknown_rate_variations.
                setattr(self, key, value)
            else:
                unknown.append(key)

        if unknown:
            self._report_unknown_keys(unknown, strict)
        return user_keys

    def _report_unknown_keys(self, unknown: list, strict: bool):
        """Report unknown parameter keys with "did you mean ...?" suggestions.

        Each unknown key is matched against the known ``DEFAULT_PARAMS`` names
        with :func:`difflib.get_close_matches`, so a typo like ``"Omegab2h"``
        is met with ``did you mean 'Omegabh2'?`` rather than a silent no-op.
        When ``strict`` (``strict_params=True``) the first batch of unknown
        keys raises ``ValueError``; otherwise it is a ``UserWarning`` (the
        back-compatible default) that also surfaces the suggestions.
        """
        details = []
        for key in sorted(unknown):
            matches = difflib.get_close_matches(key, DEFAULT_PARAMS.keys(), n=3)
            if matches:
                hint = " or ".join(repr(m) for m in matches)
                details.append(f"{key!r} (did you mean {hint}?)")
            else:
                details.append(repr(key))
        msg = "PRIMATConfig: unknown parameter key(s): " + ", ".join(details)
        if strict:
            raise ValueError(msg + " [strict_params=True]")
        warnings.warn(msg, stacklevel=3)

    def _validate_fEDE(self):
        """fEDE is a fraction of the total energy density at its peak, so it
        must satisfy 0 <= fEDE < 1.  The formula in background._setup_ede()
        has (1 - fEDE) in the denominator, which diverges at fEDE = 1.
        """
        if not (0. <= self.fEDE < 1.):
            raise ValueError(
                f"fEDE={self.fEDE!r} is out of range: must satisfy 0 ≤ fEDE < 1 "
                "(fEDE is the EDE fraction of the total energy density at its peak)."
            )

        # wnEDE is only consulted when EDE is switched on at all.
        if self.fEDE == 0.:
            return

        # background._setup_EDE locates the peak of the EDE fraction by solving
        #     d/du [ u^4 / (1 + u^q) ] = 0,   u = a/a_c,  q = 3(1 + wnEDE)
        # whose root is u^q = 4/(q - 4) = 4/(3 wnEDE - 1).  That root exists
        # only for q > 4, i.e. wnEDE > 1/3: for wnEDE <= 1/3 the EDE component
        # dilutes no faster than radiation, so its *fraction* of the total
        # never peaks during radiation domination and "fEDE at the peak" is
        # simply undefined in this parametrisation.  This is a genuine domain
        # limit, not a removable algebraic singularity.
        #
        # Caught here because the downstream failure is opaque: at wnEDE = 1/3
        # exactly, 4/(3 wnEDE - 1) raises ZeroDivisionError from inside
        # _setup_EDE; below 1/3 the base is negative and Python's float ** frac
        # silently returns a *complex* number, which then propagates through
        # amaxEDE/TmaxEDE/rhocEDEac and only surfaces hundreds of lines later
        # as solve_ivp's "`y0` is complex" -- naming neither EDE nor wnEDE.
        # (The C backend's pow() does not raise at all: it yields NaN.)
        #
        # Note the standard axion-like potential V ∝ (1 − cos φ)^n gives
        # wn = (n−1)/(n+1), so n = 1 -> 0 and n = 2 -> 1/3 both land in the
        # excluded region; n >= 3 (wn >= 1/2) is the usual EDE regime.
        if self.wnEDE <= 1. / 3.:
            raise ValueError(
                f"wnEDE={self.wnEDE!r} is out of range: must satisfy wnEDE > 1/3 "
                "when fEDE > 0. The EDE peak scale factor solves "
                "u^(3(1+wnEDE)) = 4/(3*wnEDE - 1), which has no solution for "
                "wnEDE ≤ 1/3 -- such a component dilutes no faster than "
                "radiation, so its energy fraction never peaks during radiation "
                "domination and fEDE (defined at that peak) is meaningless. "
                "For V ∝ (1 - cos φ)^n use wnEDE = (n-1)/(n+1) with n ≥ 3."
            )

    def _validate_custom_background(self):
        """custom_background: force instantaneous decoupling and no spectral
        distortions (the custom-background driver does not load NEVO tables
        and uses the analytic T_ν(T_γ) formula instead).  Must run before
        :meth:`_validate_physics_flag_combos` (external_scale_factor /
        spectral_distortions) so those see the corrected flag values.
        """
        if self.custom_background is not None:
            if self.external_scale_factor:
                raise ValueError(
                    "custom_background and external_scale_factor are mutually "
                    "exclusive: external_scale_factor reads a(T_γ) from the "
                    "NEVO table, which is not loaded in custom_background mode."
                )
            forced = []
            if self.incomplete_decoupling:
                forced.append("incomplete_decoupling=False")
                object.__setattr__(self, 'incomplete_decoupling', False)
            if self.spectral_distortions:
                forced.append("spectral_distortions=False")
                object.__setattr__(self, 'spectral_distortions', False)
            if forced:
                warnings.warn(
                    f"custom_background: forcing {', '.join(forced)} "
                    "(custom-background mode uses instantaneous-decoupling "
                    "weak rates; spectral distortions are not supported).",
                    stacklevel=2,
                )

    def _validate_dir_field(self, field: str):
        """Eagerly validate one directory-valued override and record its
        startup notice.

        Parameters
        ----------
        field : str
            ``"data_dir"`` (full-takeover data root) or ``"user_nuclear_dir"``
            (additive nuclear overlay).  ``None`` (the default for both) is a
            no-op.

        Raises
        ------
        ValueError
            If the path is not an existing directory, or -- for ``data_dir``,
            which *replaces* the shipped tree -- if it exists but is not a data
            tree at all (no ``csv/``, no ``nuclear/``).  The structural check
            matters because ``data_dir`` is consumed a few lines later by
            ``_load_nuclide_data``: without it the user sees a bare
            ``FileNotFoundError: <data_dir>/csv/nuclides.csv`` naming neither
            the parameter nor what the directory was supposed to contain.
            ``NEVO/`` and ``cache_plasma_weak/`` are deliberately *not*
            required: the former is unused when ``incomplete_decoupling=False``
            and its own resolver reports a clear "not found", and the latter
            holds only regenerable caches (a missing one costs recompute time,
            not correctness).

        Example
        -------
            >>> cfg._validate_dir_field("data_dir")   # doctest: +SKIP
            ValueError: data_dir='/tmp' exists but does not look like a primat
            data tree: missing subdirectories csv/, nuclear/ ...
        """
        value = getattr(self, field)
        if value is None:
            return
        if not os.path.isdir(value):
            raise ValueError(f"{field}={value!r} is not an existing directory")
        if field == "data_dir":
            missing = [sub for sub in ("csv", "nuclear")
                       if not os.path.isdir(os.path.join(value, sub))]
            if missing:
                raise ValueError(
                    f"data_dir={value!r} exists but does not look like a primat "
                    f"data tree: missing subdirectories "
                    f"{', '.join(sub + '/' for sub in missing)}. A data_dir "
                    "replaces the shipped primat/data/ tree entirely and must "
                    "carry csv/ (nuclides.csv, reactions_large.csv, "
                    "detailed_balance.csv) and nuclear/ (networks/, tables/), "
                    "plus NEVO/ unless incomplete_decoupling=False and "
                    "cache_plasma_weak/{weak,plasma}/ for the regenerable "
                    "caches. Use user_nuclear_dir for an additive overlay "
                    "instead of a full takeover."
                )
        self._init_messages.append(_rates_overlay_notice(field, value))

    def _validate_data_dirs(self):
        """Eagerly validate ``user_nuclear_dir`` (mirrors the ``nevo_file``
        pattern in :meth:`_validate_nevo_files`) so a typo'd override path
        fails fast at construction time rather than surfacing as a confusing
        "network not found" later.

        ``data_dir`` is *not* handled here: it is needed -- and therefore
        checked -- much earlier, in
        :meth:`_init_defaults_and_nuclide_data`, before ``nuclides.csv`` is
        read from it.  ``cache_dir`` is deliberately unvalidated: it is a
        write target created on demand (see its ``DEFAULT_PARAMS`` entry), so
        a not-yet-existing directory is a normal, supported input.
        """
        self._validate_dir_field("user_nuclear_dir")

    def _validate_ranges(self):
        """Cross-field numeric consistency checks -- the constraints that
        involve *two* parameters and so cannot live in the per-key
        :data:`_PARAM_RANGE` table.

        Each of these is accepted by the per-key checks yet leaves the solver
        with an impossible request, and each used to surface far downstream as
        an opaque integrator failure (or, worse, silently):

        - ``rate_grid_T9_min >= rate_grid_T9_max`` -- ``np.logspace`` then
          builds a *descending* master T9 grid, breaking ``fill_buffer``'s
          single ``searchsorted`` path; both backends die with a "step size
          underflowed"/"required step size is less than spacing between
          numbers" message from the MT era, naming the ODE solver rather than
          the two grid parameters.
        - ``T_end_MeV >= T_start_cosmo_MeV`` -- the background is integrated
          from ``T_start_cosmo_MeV`` *down* to ``T_end_MeV``, so an inverted
          pair asks for a zero-or-negative-length temperature range.
        - ``mc_rate_rescale_cap < 1`` -- the cap is applied as a clamp to
          ``[1/cap, cap]`` (``network_data.py``'s ``_apply_rate_variation``),
          whose bounds *cross* below 1: with ``cap = 0.5`` every sampled rate
          factor is pinned to 0.5, i.e. every MC sample divides every rate by
          two.  A cap of exactly 1 means "no variation at all", which is
          allowed (and ``None`` disables the cap).

        Mirrored in ``primat-c/src/config.c``'s ``cpr_config_validate`` per
        CLAUDE.md's parity mandate.

        Example
        -------
            >>> PRIMATConfig({"rate_grid_T9_min": 10., "rate_grid_T9_max": 1e-3})
            Traceback (most recent call last):
            ValueError: rate_grid_T9_min=10.0 must be < rate_grid_T9_max=0.001 ...
        """
        if self.rate_grid_T9_min >= self.rate_grid_T9_max:
            raise ValueError(
                f"rate_grid_T9_min={self.rate_grid_T9_min!r} must be < "
                f"rate_grid_T9_max={self.rate_grid_T9_max!r}: they bound the "
                "log-spaced master T9 grid every nuclear rate table is "
                "resampled onto, which must be increasing."
            )
        if self.T_end_MeV >= self.T_start_cosmo_MeV:
            raise ValueError(
                f"T_end_MeV={self.T_end_MeV!r} must be < "
                f"T_start_cosmo_MeV={self.T_start_cosmo_MeV!r}: the background "
                "and the nuclear network are integrated from T_start_cosmo_MeV "
                "down to T_end_MeV."
            )
        if self.mc_rate_rescale_cap is not None and self.mc_rate_rescale_cap < 1.:
            raise ValueError(
                f"mc_rate_rescale_cap={self.mc_rate_rescale_cap!r} must be >= 1 "
                "(or None to disable the cap): it clamps the MC rate-variation "
                "factor to [1/cap, cap], whose bounds cross below 1 -- a cap "
                "of 0.5 would pin every sampled rate to half its median."
            )

    def _validate_network(self):
        """Check that a non-``"small"`` ``network`` name resolves to an
        existing ``data/nuclear/networks/<name>.txt`` file (including any
        ``user_nuclear_dir`` overlay), raising a clear error listing every
        path searched otherwise.
        """
        if self.network != "small":
            path = self.resolve_rates_path("nuclear", "networks", f"{self.network}.txt")
            if not os.path.exists(path):
                searched = []
                if self.user_nuclear_dir is not None:
                    searched.extend(_overlay_candidates(
                        self.user_nuclear_dir,
                        os.path.join("nuclear", "networks", f"{self.network}.txt"),
                    ))
                searched.append(path)
                raise ValueError(
                    f"network must be 'small' or name an existing file in "
                    f"data/nuclear/networks; missing {path!r}"
                    + (f" (searched: {', '.join(repr(p) for p in searched)})" if searched else "")
                )

    def _setup_rate_variation_defaults(self) -> set:
        """Default every reaction of the *configured* network (self.network,
        finalised by :meth:`_apply_user_overrides`) to p_<rxn>=0 /
        delta_<rxn>=0, i.e. "no rate variation".  Uses ``setdefault`` so any
        p_<rxn>/delta_<rxn> override already applied is not clobbered.
        Returns the set of valid bare reaction names for this network (used
        by :meth:`_warn_unknown_rate_variations` to catch typos).
        """
        from .network_data import load_reaction_names, reaction_category
        reactions_with_tables = load_reaction_names(self, self.network)
        # Each entry is "bare_name" or "bare_name, filename.txt"; only the
        # bare reaction name is used as the p_<rxn>/delta_<rxn> key.
        # amax (now meaningful for any network, not just "large") must be
        # applied here too, so p_rxn/delta_rxn don't carry stale keys for
        # reactions load_network would have dropped.
        valid_rxns = set()
        for entry in reactions_with_tables:
            bare = re.split(r'[, ]+', entry, maxsplit=1)[0]
            if (self.amax is not None
                    and reaction_category(bare, self._resolved_data_dir) > self.amax):
                continue
            valid_rxns.add(bare)
        for rxn in valid_rxns:
            self.p_rxn.setdefault(rxn, 0.0)
            self.delta_rxn.setdefault(rxn, 0.0)
        return valid_rxns

    def _warn_unknown_rate_variations(self, user_keys: set, valid_rxns: set):
        """Catch p_<rxn>/delta_<rxn> typos in the constructor params. This
        has to happen after the network's reaction list is known (rather
        than in ``__setattr__`` at override time in
        :meth:`_apply_user_overrides`, which runs before ``self.network`` is
        finalised -- ``self.network`` may itself be one of the overrides),
        and by this point ``__setattr__``'s routing has already inserted the
        (possibly bogus) key into self.p_rxn/self.delta_rxn -- so we must
        check against ``valid_rxns`` (from :meth:`_setup_rate_variation_defaults`),
        not against those dicts.

        ``strict_params`` applies here too, exactly as it does to unknown
        ``DEFAULT_PARAMS`` keys in :meth:`_report_unknown_keys`: a mistyped
        reaction name (``p_n_p__d_gg``) is the *archetypal* silent no-op that
        strict mode exists to catch -- reaction names are long and
        underscore-heavy, so they are likelier to be mistyped than a config
        key, and the run otherwise proceeds with that rate unvaried. All bad
        keys are reported in one message rather than one exception per key.
        """
        unmatched = sorted(
            key for key in user_keys
            if (key.startswith('p_') and key[2:] not in valid_rxns)
            or (key.startswith('delta_') and key[6:] not in valid_rxns)
        )
        if not unmatched:
            return
        plural = "s" if len(unmatched) > 1 else ""
        msg = (f"PRIMATConfig: rate-variation key{plural} "
               f"{', '.join(repr(k) for k in unmatched)} "
               f"do{'' if plural else 'es'} not match any reaction in network "
               f"{self.network!r}; {'they have' if plural else 'it has'} no "
               "effect on the run.")
        if self.strict_params:
            raise ValueError(msg + " [strict_params=True]")
        warnings.warn(msg, stacklevel=2)

    def _detect_optional_libraries(self):
        """Detect optional libraries (numba) for flags not explicitly set by
        the caller.  Messages are stored in ``_init_messages`` for deferred
        printing (after the banner).
        """
        if self.numba_installed:
            try:
                import numba  # noqa: F401
                self._init_messages.append('[init]  numba detected: using it for JIT compilation.')
            except ImportError:
                self.numba_installed = False
                self._init_messages.append('[init]  numba not detected: running without JIT compilation.')

    def _validate_amax(self):
        """amax must be None or a positive integer."""
        # numbers.Integral (not bare int) so a numpy int (np.int64, common in
        # MCMC drivers) is accepted -- bool is an Integral subclass but is
        # already rejected upstream by _validate_param_value's type check.
        if self.amax is not None:
            if not (isinstance(self.amax, numbers.Integral) and self.amax >= 1):
                raise ValueError(
                    f"amax must be None or a positive integer (got {self.amax!r})."
                )

    def _validate_nevo_files(self):
        """Validate any custom NEVO table overrides: check each file exists
        and has the column/length count expected by
        ``neutrino_history.NEVOTable``, so a typo or malformed file is
        caught here with a clear message rather than as a confusing shape
        mismatch deep inside an interpolant. Covers ``nevo_file``,
        ``nevo_spectral_file``, ``nevo_grid_file``, and ``nevo_file_prefix``
        (which derives its filenames from the other three's naming
        convention and validates only the ones not already overridden
        individually).
        """
        from .neutrino_history import resolve_nevo_path
        if self.nevo_file is not None:
            path = resolve_nevo_path(self, self.nevo_file, "")
            if not os.path.exists(path):
                raise ValueError(f"nevo_file={self.nevo_file!r} not found "
                                  f"(resolved to {path!r})")
            ncols = np.loadtxt(path, delimiter=',', max_rows=1).size
            if ncols not in (6, 7):
                raise ValueError(f"nevo_file={self.nevo_file!r} ({path!r}) has "
                                  f"{ncols} columns; expected 6 or 7 (the NEVO "
                                  f"x,z,Tnue,Tnumu,Tnutau,N[,extra] thermo table)")

        # The shipped spectral table has 86 columns TOTAL = 6 thermo + 80
        # spectral, so the shipped NEVOGrid.csv has 80 y-nodes.  (Both numbers
        # used to be quoted as "86", which read as 86 spectral columns.)
        n_grid_nodes = 80  # shipped NEVOGrid.csv length; overridden below if a
                           # custom spectral table declares a different width
        spectral_overridden = self.nevo_spectral_file is not None
        if spectral_overridden:
            path = resolve_nevo_path(self, self.nevo_spectral_file, "")
            if not os.path.exists(path):
                raise ValueError(f"nevo_spectral_file={self.nevo_spectral_file!r} "
                                  f"not found (resolved to {path!r})")
            ncols = np.loadtxt(path, delimiter=',', max_rows=1).size
            if ncols <= 6:
                raise ValueError(f"nevo_spectral_file={self.nevo_spectral_file!r} "
                                  f"({path!r}) has {ncols} columns; expected "
                                  f"6 thermo columns plus at least one spectral "
                                  f"column (86 columns total in the shipped "
                                  f"tables: 6 thermo + 80 spectral)")
            n_grid_nodes = ncols - 6

        # Resolve the y-grid actually in play: the override if given, else the
        # shipped NEVOGrid.csv.  Checking the *shipped* grid too closes the
        # case where only nevo_spectral_file is overridden -- previously
        # n_grid_nodes was computed and then never compared against anything,
        # so a custom spectral table of the wrong width sailed past this method
        # and failed later inside NEVOTable's RegularGridInterpolator, which is
        # exactly the confusing deep-shape-mismatch this method exists to
        # prevent.
        if self.nevo_grid_file is not None:
            grid_path = resolve_nevo_path(self, self.nevo_grid_file, "")
            grid_desc = f"nevo_grid_file={self.nevo_grid_file!r}"
            if not os.path.exists(grid_path):
                raise ValueError(f"{grid_desc} not found "
                                  f"(resolved to {grid_path!r})")
        elif spectral_overridden:
            grid_path = resolve_nevo_path(self, None, "NEVOGrid.csv")
            grid_desc = "the shipped NEVOGrid.csv"
            if not os.path.exists(grid_path):
                grid_path = None   # nothing to check against
        else:
            grid_path = None       # neither overridden: shipped pair, known good

        if grid_path is not None:
            n_nodes = np.loadtxt(grid_path, delimiter=',').size
            if n_nodes != n_grid_nodes:
                raise ValueError(f"{grid_desc} ({grid_path!r}) has {n_nodes} "
                                  f"nodes; expected {n_grid_nodes} to match the "
                                  f"spectral table's {n_grid_nodes} y-columns. "
                                  f"Override nevo_grid_file to supply a matching "
                                  f"grid.")

        # Validate nevo_file_prefix: when not the shipped default, check that
        # the *derived* default filenames it implies exist and have the right
        # shape -- mirrors the nevo_file/nevo_spectral_file checks above, but
        # only for the files that aren't already overridden individually.
        if self.nevo_file_prefix != "NEVOPRIMAT" and self.incomplete_decoupling:
            prefix = self.nevo_file_prefix
            suffix = "" if self.QED_corrections else "_NoQED"

            if self.nevo_file is None:
                fname = f"{prefix}{suffix}_col_1_7.csv"
                path = resolve_nevo_path(self, None, fname)
                if not os.path.exists(path):
                    raise ValueError(f"nevo_file_prefix={prefix!r}: derived "
                                      f"thermo file {fname!r} not found "
                                      f"(resolved to {path!r})")
                ncols = np.loadtxt(path, delimiter=',', max_rows=1).size
                if ncols not in (6, 7):
                    raise ValueError(f"nevo_file_prefix={prefix!r}: "
                                      f"{fname!r} has {ncols} columns; "
                                      f"expected 6 or 7")

            if (self.spectral_distortions and not self.analytic_distortions
                    and self.nevo_spectral_file is None):
                fname = f"{prefix}{suffix}.csv"
                path = resolve_nevo_path(self, None, fname)
                if not os.path.exists(path):
                    raise ValueError(f"nevo_file_prefix={prefix!r}: derived "
                                      f"spectral file {fname!r} not found "
                                      f"(resolved to {path!r})")
                ncols = np.loadtxt(path, delimiter=',', max_rows=1).size
                if ncols <= 6:
                    raise ValueError(f"nevo_file_prefix={prefix!r}: "
                                      f"{fname!r} has {ncols} columns; "
                                      f"expected > 6")

    def _validate_physics_flag_combos(self):
        """Validate flag combinations that depend on ``incomplete_decoupling``:
        ``external_scale_factor`` (reads a(T) from the NEVO table, so
        requires it to be loaded) and ``spectral_distortions`` (analytic
        distortions require instantaneous decoupling; full-spectrum
        distortions require the NEVO table).
        """
        # external_scale_factor reads a(T) directly from the NEVO table's x
        # column (NEVOTable.x_of_Tg), so it requires the NEVO table to be
        # loaded in the first place.
        if self.external_scale_factor and not self.incomplete_decoupling:
            raise ValueError(
                "external_scale_factor=True requires incomplete_decoupling=True "
                "(a(T) is read from the NEVO table, which is only loaded by "
                "NEVOTable)."
            )

        # Validate spectral-distortion flag combination.
        if self.spectral_distortions:
            if self.analytic_distortions:
                if self.incomplete_decoupling:
                    raise ValueError(
                        "spectral_distortions=True with analytic_distortions=True "
                        "requires instantaneous decoupling (incomplete_decoupling=False)."
                    )
            else:
                if not self.incomplete_decoupling:
                    raise ValueError(
                        "spectral_distortions=True with analytic_distortions=False "
                        "requires incomplete_decoupling=True (the full NEVO spectrum "
                        "file is only available in the non-instantaneous decoupling mode)."
                    )

    def __getattr__(self, name: str):
        """Dynamic lookup for nuclear rate variations p_* and delta_*."""
        if name.startswith("p_"):
            return object.__getattribute__(self, 'p_rxn').get(name[2:], 0.0)
        if name.startswith("delta_"):
            return object.__getattribute__(self, 'delta_rxn').get(name[6:], 0.0)
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def __dir__(self):
        """Extend the default ``dir()``/tab-completion listing with the
        dynamic ``p_<rxn>``/``delta_<rxn>`` rate-variation attributes
        currently set on this instance (via ``__getattr__``/``__setattr__``
        above, so they never appear in ``object.__dir__`` on their own).
        The static ``DEFAULT_PARAMS`` keys are already real instance
        attributes by construction time and need no help from here; this
        only closes the gap for the unbounded p_*/delta_* pattern, which the
        ``TYPE_CHECKING`` annotation block above cannot enumerate either
        (any reaction name is valid, including ones from a custom network).
        """
        names = set(object.__dir__(self))
        names.update(f"p_{rxn}" for rxn in self.p_rxn)
        names.update(f"delta_{rxn}" for rxn in self.delta_rxn)
        return sorted(names)

    def __setattr__(self, name: str, value):
        """Dynamic routing for nuclear rate variations p_* and delta_*."""
        # The two backing dicts themselves start with those very prefixes, so
        # they must be assigned as plain attributes -- otherwise
        # `cfg.p_rxn = {...}` (e.g. round-tripping a saved dict) would be read
        # as a variation of a reaction literally named "rxn" and die in
        # float(dict).
        if name in ("p_rxn", "delta_rxn"):
            object.__setattr__(self, name, value)
        elif name.startswith("p_"):
            object.__getattribute__(self, 'p_rxn')[name[2:]] = float(value)
        elif name.startswith("delta_"):
            object.__getattribute__(self, 'delta_rxn')[name[6:]] = float(value)
        else:
            if name in _PATH_PARAMS:
                # Normalize "~" immediately so both direct assignment and
                # --set KEY=VALUE route through the same resolved path.
                value = _expanduser_path(value)
            object.__setattr__(self, name, value)
            if name == "Omegabh2":
                self._update_derived()

    # Omegabh2 is exposed as a property so that the derived baryon-to-photon
    # ratio eta0b is recomputed automatically whenever it is reassigned (by
    # attribute, by the constructor loop, or via __setitem__).
    @property
    def Omegabh2(self) -> float:
        return self._Omegabh2

    @Omegabh2.setter
    def Omegabh2(self, value: float):
        self._Omegabh2 = value
        self._update_derived()

    def _update_derived(self):
        """Recompute quantities that depend on mutable parameters."""
        self.Omegabh2_to_eta0b = (self.rhocOverh2 / self.n0CMB) / (self.ma / self.maOvermB)
        self.eta0b = self.Omegabh2_to_eta0b * self._Omegabh2

    # Convenience: allow dict-style access for backwards compat if needed
    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
        self._update_derived()

    # Class-level storage to avoid AttributeError if accessed before init
    Nuclides: dict[str, list[int]] = {}
    NuclExcessMass: dict[str, float] = {}
    NuclSpin: dict[str, float] = {}

    def _load_nuclide_data(self):
        """Load mass excesses, spins, and (N, Z) from nuclides.csv."""
        import csv
        path = os.path.join(self._resolved_data_dir, "csv", "nuclides.csv")
        
        self.Nuclides = {}
        self.NuclExcessMass = {}
        self.NuclSpin = {}
        
        with open(path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row['name']
                self.Nuclides[name] = [int(row['N']), int(row['Z'])]
                self.NuclExcessMass[name] = float(row['mass_excess_keV'])
                self.NuclSpin[name] = float(row['spin'])

    @property
    def _pkg_data_dir(self) -> str:
        """Package-shipped data root (``primat/data/``, contains NEVO/, nuclear/, weak/, plasma/, csv/).

        This is the fixed fallback used when ``data_dir`` param is ``None``.
        It always points to the installed package's own data tree regardless
        of any user override.
        """
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

    @property
    def _resolved_data_dir(self) -> str:
        """Resolved data root: the ``data_dir`` param when set, otherwise ``primat/data/``.

        Use this everywhere a data-root path is needed instead of the old
        ``cfg.data_dir + "/data"`` idiom.  Output paths are still resolved
        against the current working directory (see PRIMAT._write_time_evolution
        / _write_final_result).
        """
        return self.data_dir if self.data_dir else self._pkg_data_dir

    def resolve_rates_path(self, *parts: str) -> str:
        """Resolve a path inside the nuclear data tree through the overlay chain.

        Used by every caller that needs a nuclear network file or rate-table
        file, so a user's ``user_nuclear_dir`` additive overlay (or a full
        ``data_dir`` takeover — see those fields in ``DEFAULT_PARAMS``) is
        honoured without touching the installed ``primat`` package.

        Lookup order (first existing path wins):
          1. ``self.user_nuclear_dir`` (additive nuclear overlay), if set.
          2. ``self._resolved_data_dir`` (either the user-supplied ``data_dir``
             or the shipped ``primat/data/`` tree — always tried last so
             ``small``/``large`` and the default rate tables are never
             unreachable just because ``user_nuclear_dir`` is also configured).

        If the relative path is not found under any candidate base, the
        resolved-default path is returned anyway (not found), so callers get
        a "missing file" error that points at the expected default location
        rather than at whichever overlay happened to be checked last.

        Args:
            *parts: path components relative to a nuclear data root, e.g.
                ``"nuclear", "networks", "large.txt"``.

        Returns:
            str: an absolute path (existing, if found under any candidate
            base; otherwise the resolved-default path, for use in error
            messages).

        Example:
            >>> cfg.resolve_rates_path("nuclear", "networks", "large.txt")
            '/.../primat/data/nuclear/networks/large.txt'
        """
        relpath = os.path.join(*parts) if parts else ""
        bases = []
        if self.user_nuclear_dir:
            bases.append(self.user_nuclear_dir)
        bases.append(self._resolved_data_dir)  # shipped (or overridden) default, always last
        for base in bases:
            if relpath:
                for candidate in _overlay_candidates(base, relpath):
                    if os.path.exists(candidate):
                        return candidate
            else:
                if os.path.exists(base):
                    return base
        return os.path.join(bases[-1], relpath) if relpath else bases[-1]
