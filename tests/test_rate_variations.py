"""
Tests for nuclear rate variation and MC uncertainty propagation.

The ``p_<reaction>`` mechanism shifts a reaction rate by ``exp(p × σ)``
relative to its median value, enabling MCMC sampling of nuclear-rate
uncertainties.  The ``delta_<reaction>`` mechanism adds a direct fractional
shift on top of any ``p_`` variation (delta=0.1 → +10% on the rate).
These tests verify that:
1. varying a rate (p_* or delta_*) actually changes the predicted abundances;
2. restoring p=0/delta=0 reproduces the baseline to floating-point precision;
3. delta_* works when passed via run_bbn's params dict (the original bug:
   delta was silently gated behind a separate flag and had no effect);
4. the MC runner propagates rate uncertainty to non-zero spread in observables.

The ``test_config_dynamic_attr`` test (attribute routing for p_* / delta_*)
lives in ``test_config.py`` where it logically belongs.
"""
import numpy as np
import pytest

from primat import PRIMAT, mc_uncertainty
from primat.backend import run_bbn


@pytest.mark.slow
@pytest.mark.solve
def test_solve_variation():
    """Varying p_n_p__d_g shifts D/H; reverting p=0 restores the baseline."""
    inst = PRIMAT(params={"network": "small", "verbose": False})
    res0 = inst.solve()
    dh0  = res0["DoH"]

    # Shift n_p__d_g by +1σ and re-solve
    inst.cfg.p_n_p__d_g = 1.0
    res1 = inst.solve()
    dh1  = res1["DoH"]
    assert dh1 != dh0, "Changing p_n_p__d_g should affect D/H"

    # Restore and verify exact match (deterministic ODE)
    inst.cfg.p_n_p__d_g = 0.0
    res2 = inst.solve()
    dh2  = res2["DoH"]
    assert np.isclose(dh2, dh0, rtol=1e-10), (
        f"Reverting p_n_p__d_g should match baseline: {dh2:.8e} vs {dh0:.8e}"
    )


@pytest.mark.slow
@pytest.mark.solve
def test_run_bbn_delta_variation():
    """delta_<rxn> passed via run_bbn params changes D/H on both backends.

    This is the regression test for the bug where delta_<rxn> was silently
    gated behind a separate enabling flag and had no effect when that flag
    was False (the default).  Passing delta_n_p__d_g=0.1 must shift D/H
    relative to the baseline without any extra flags.
    """
    base = {"network": "small", "verbose": False}
    res0 = run_bbn(base)
    dh0  = res0["DoH"]

    res1 = run_bbn({**base, "delta_n_p__d_g": 0.1})
    dh1  = res1["DoH"]
    assert dh1 != dh0, (
        "delta_n_p__d_g=0.1 via run_bbn should shift D/H "
        f"(got {dh1:.8e} == {dh0:.8e})"
    )

    # Confirm the Python backend also applies it
    res2 = run_bbn({**base, "delta_n_p__d_g": 0.1}, force_backend="python")
    dh2  = res2["DoH"]
    assert dh2 != dh0, (
        "delta_n_p__d_g=0.1 via run_bbn (Python backend) should shift D/H "
        f"(got {dh2:.8e} == {dh0:.8e})"
    )


@pytest.mark.slow
@pytest.mark.solve
def test_mc_large_network():
    """MC uncertainty spread is positive for D/H and B10 in the large network."""
    mc = mc_uncertainty(5, ["DoH", "B10"],
                        params={"network": "large"}, n_jobs=-1)
    assert mc["DoH"].std > 0, "D/H should have non-zero uncertainty"
    assert mc["B10"].std > 0, "B10 should have non-zero uncertainty in large network"


# ---------------------------------------------------------------------------
# Derived state that must be refreshed alongside the forward-rate table
#
# GOAL: pin that ``NetworkDefinition.apply_variations`` updates *everything*
# that depends on the active forward rates, not just ``_fwd`` itself.  Two
# pieces of derived state were previously left stale, so a variation was only
# half-applied -- silently, with no error and plausible-looking numbers.
# ---------------------------------------------------------------------------
def test_apply_variations_invalidates_fill_buffer_cache():
    """``fill_buffer`` memoises on ``(T_t, clamp)`` only, so it cannot see that
    the rate table changed underneath it.  ``apply_variations`` must drop that
    one-slot cache.

    Regression guard: without the invalidation, the first ``fill_buffer`` call
    after a variation returned the *previous* solve's rates whenever it landed
    on a bit-identical temperature -- the exact situation an MC loop reusing one
    ``UpdateNuclearRates`` across samples creates.  The C backend has always
    done this (``cpr_network_apply_variations``'s ``net->cache_valid = 0``), so
    this is also a backend-parity guard.
    """
    from primat.config import PRIMATConfig
    from primat.network_data import load_network

    cfg = PRIMATConfig({"network": "small", "verbose": False})
    net = load_network(cfg, era="LT")
    frwrd, bkwrd = (lambda T: 1.0), (lambda T: 2.0)

    T = 1.0e9
    before = net.fill_buffer(T, frwrd, bkwrd, clamp=True).copy()
    cfg.delta_n_p__d_g = 1.0                      # +100% on the forward rate
    net.apply_variations(cfg)
    after = net.fill_buffer(T, frwrd, bkwrd, clamp=True).copy()

    # r[2] is n_p__d_g's forward slot (r[0]/r[1] are the weak n<->p rates).
    assert after[2] == pytest.approx(2.0 * before[2], rel=1e-12), (
        "fill_buffer returned a stale rate after apply_variations")


def test_apply_variations_rescales_the_reverse_rate_cap():
    """The reverse-rate cap is ``bwd(T_nucl)``, hence proportional to the
    forward rate: it must scale with the variation, not stay at the median.

    Regression guard: ``_bwd_cap`` was computed once in ``load_network`` from
    the *median* table and never rebuilt.  Wherever the clamp binds, a varied
    reverse rate was therefore clamped back to its unvaried value -- breaking
    detailed balance by exactly the variation factor, in the one code path
    (MC uncertainty propagation) whose whole purpose is varying rates.

    ``B10_p__a_a_He3`` is one of only two reactions in the full ``large``
    network whose clamp actually binds inside the LT temperature range, which
    is what makes it the sharp test case.
    """
    from primat.config import PRIMATConfig
    from primat.network_data import load_network

    cfg = PRIMATConfig({"network": "large", "verbose": False})
    net = load_network(cfg, era="LT")
    i = net.names.index("B10_p__a_a_He3") - 1     # -1: names[0] is the weak n__p
    cap_median = float(net._bwd_cap[i])
    assert cap_median > 0.0

    cfg.delta_B10_p__a_a_He3 = 0.5                # +50% on the forward rate
    net.apply_variations(cfg)
    assert net._bwd_cap[i] == pytest.approx(1.5 * cap_median, rel=1e-12)

    # Restoring the baseline must restore the cap exactly, so repeated MC
    # samples do not drift.
    cfg.delta_B10_p__a_a_He3 = 0.0
    net.apply_variations(cfg)
    assert net._bwd_cap[i] == pytest.approx(cap_median, rel=1e-12)


if __name__ == "__main__":
    test_solve_variation()
    test_mc_large_network()
