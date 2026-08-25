# -*- coding: utf-8 -*-
"""
GOAL: prove each fingerprinted cache declares exactly the physical constants
it reads.

Every cache is keyed on ``cache_utils.CACHE_CONSTANTS[<cache>]`` rather than on
the whole ``Constants`` struct, so an omission would silently serve one
configuration's table to another. This module perturbs every settable constant
and rebuilds each cache's DATA from scratch -- never through the cache files,
so the answer cannot depend on the fingerprint it is checking:

* **undeclared => inert.** A constant outside a cache's list must leave its
  data bit-identical, under four flag combinations that change which code path
  fills the cache. This is the safety direction: a failure here means the cache
  is under-keyed, i.e. wrong physics with no error.
* **declared => live.** A constant inside the list must move the data at the
  shipped default flags. This is the tightness direction: a failure means the
  cache re-keys on something it does not read, costing (for CCRTh) a
  multi-minute recompute.

Every step is deterministic -- quadrature, and a seeded vegas -- so "unchanged"
means exactly zero, not "within a tolerance"; ``test_rebuild_is_bit_identical``
pins that premise first.

The grids are deliberately coarse (a 3-row CCRTh table, 40-point electron
thermo): the dependency structure does not depend on resolution, and the
default grids would take hours.
"""
import numpy as np
import pytest

from primat.cache_utils import CACHE_CONSTANTS
from primat.config import PRIMATConfig
from primat.constants import OVERRIDABLE_CONSTANTS
from primat.plasma import Plasma
from primat.background import StandardBackground
from primat.qed_pressure import compute_qed_pressure_tables
from primat.weak_rates.api import ComputeWeakRates
from primat.weak_rates.corrections import (_build_rate_context, _L_CCRTh_compute,
                                           _ThermalIntegOpts, _T_CCRTH_MIN)

pytestmark = pytest.mark.slow

# Coarse everything: 16 constants x 4 caches is 17 full rebuilds per flag
# combination, so each one has to stay near a second.
_COARSE = dict(sampling_nTOp_per_decade=4, sampling_nTOp_thermal_per_decade=1,
               n_electron_table=40, vegas_n_eval=2000, vegas_n_itn=4,
               weak_rate_cache=False, save_nTOp=False, save_nTOp_thermal=False,
               recompute_electron_thermo=True)
_QED_NPTS = 12

# Big enough that a real dependence cannot hide in the last bits, small enough
# that Q = mn - mp stays physical when a nucleon mass is the one perturbed.
_PERTURBATION = 1e-4

# Flag combinations that route the caches through different code (analytic vs
# tabulated neutrino history, Born vs radiatively corrected rates, tabulated vs
# analytic spectral distortions, QED-corrected decoupling or not). The
# undeclared => inert direction must hold in all of them.
_FLAG_COMBOS = {
    "default": {},
    "instantaneous_decoupling": dict(incomplete_decoupling=False,
                                     spectral_distortions=False),
    "born_only": dict(radiative_corrections=False, finite_mass_corrections=False),
    "analytic_distortions": dict(analytic_distortions=True,
                                 incomplete_decoupling=False,
                                 y_SZ=1e-3, y_gray=1e-3),
    "no_QED_corrections": dict(QED_corrections=False),
}


def _cfg(flags, constant=None):
    params = dict(_COARSE, **flags)
    if constant is not None:
        params[constant] = getattr(PRIMATConfig(params), constant) * (1. + _PERTURBATION)
    return PRIMATConfig(params)


def _cache_data(cfg):
    """The four caches' data arrays, computed from scratch for this config."""
    plasma = Plasma(cfg)
    T_e = np.logspace(np.log10(cfg.me / 30.),
                      np.log10(max(cfg.T_start_cosmo_MeV, 100.) * 1.5),
                      cfg.n_electron_table)
    qed = compute_qed_pressure_tables(n_pts=_QED_NPTS, alpha=cfg.alphaem,
                                      me=cfg.me, verbose=False)

    bg = StandardBackground(cfg, plasma)
    Tvec = [bg.Tg_vec, bg.Tnue_vec]
    T_all, frwrd, bkwrd = ComputeWeakRates(Tvec, cfg,
                                           dFDneu_func=bg.dFDneu_func,
                                           dFDneu_moments=bg.dFDneu_moments)

    ctx = _build_rate_context(Tvec, cfg)
    opts = _ThermalIntegOpts(True, cfg.vegas_n_eval, cfg.vegas_n_itn,
                             cfg.epsrel_thermal)
    T_th = np.logspace(np.log10(_T_CCRTH_MIN), np.log10(cfg.T_start_nucl), 3)
    return {
        "weak": np.column_stack([T_all, frwrd, bkwrd]),
        "thermal": np.column_stack(
            [T_th,
             [_L_CCRTh_compute(ctx, t, +1, opts) for t in T_th],
             [_L_CCRTh_compute(ctx, t, -1, opts) for t in T_th]]),
        "electron_thermo": np.column_stack(
            [T_e, plasma._rho_e_tab(T_e), plasma._p_e_tab(T_e),
             plasma._drho_e_dT_tab(T_e), plasma._dp_e_dT_tab(T_e)]),
        "qed": np.column_stack([qed["T"], qed["dP_e2"], qed["dP_e3"]]),
    }


def _identical(a, b):
    return a.shape == b.shape and bool(np.array_equal(a, b))


@pytest.fixture(scope="module")
def reference():
    """Each flag combination's unperturbed cache data."""
    return {name: _cache_data(_cfg(flags)) for name, flags in _FLAG_COMBOS.items()}


def test_declared_caches_match_the_code(reference):
    """CACHE_CONSTANTS names the four caches, and only real constants."""
    assert set(CACHE_CONSTANTS) == set(reference["default"])
    for cache, names in CACHE_CONSTANTS.items():
        assert set(names) <= set(OVERRIDABLE_CONSTANTS), cache
        assert list(names) == sorted(names), cache


def test_rebuild_is_bit_identical(reference):
    """Rebuilding a cache without touching a constant reproduces it exactly.

    The premise the two tests below rest on: with a seeded vegas and
    deterministic quadrature there is no noise floor, so any difference they
    see is the perturbed constant and nothing else.
    """
    again = _cache_data(_cfg({}))
    for cache, data in reference["default"].items():
        assert _identical(data, again[cache]), cache


@pytest.mark.parametrize("combo", sorted(_FLAG_COMBOS))
def test_undeclared_constants_leave_every_cache_untouched(reference, combo):
    """A constant a cache does not declare cannot change its data.

    The safety direction. An omission from CACHE_CONSTANTS would let a run
    that changed this constant load -- and report -- another configuration's
    table, with no error anywhere.
    """
    ref = reference[combo]
    for name in OVERRIDABLE_CONSTANTS:
        undeclared = [c for c, names in CACHE_CONSTANTS.items() if name not in names]
        if not undeclared:
            continue
        data = _cache_data(_cfg(_FLAG_COMBOS[combo], name))
        for cache in undeclared:
            assert _identical(ref[cache], data[cache]), (
                f"{name} is not in CACHE_CONSTANTS[{cache!r}] but changes that "
                f"cache's data under flags {combo!r}: the cache is under-keyed")


def test_declared_constants_all_move_their_cache(reference):
    """Every declared constant really is read, at the shipped default flags.

    The tightness direction: a constant listed but not read re-keys the cache
    for nothing, which for the CCRTh table is a multi-minute recompute. Checked
    at the defaults only -- switching a correction off (`born_only`) legitimately
    removes constants from the live set, and over-declaration there is safe.
    """
    ref = reference["default"]
    for name in OVERRIDABLE_CONSTANTS:
        declared = [c for c, names in CACHE_CONSTANTS.items() if name in names]
        if not declared:
            continue
        data = _cache_data(_cfg({}, name))
        for cache in declared:
            assert not _identical(ref[cache], data[cache]), (
                f"{name} is declared in CACHE_CONSTANTS[{cache!r}] but does not "
                f"change that cache's data: the cache is over-keyed")
