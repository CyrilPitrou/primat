"""Pin a few of primat.constants.CONST's *derived* values.

CONST.alphaem, CONST.GF and CONST.mZ are primary inputs (PDG values, set
verbatim); CONST.sW2 (sin^2(theta_W)) and the effective electron/muon
couplings derived from it are computed from those three via the on-shell
relation. A typo in that formula (a stray factor of 2, a swapped GF/mZ) would
not raise any error -- it would just silently shift every weak rate that uses
sW2 by a few percent. These tests pin the derived numbers against an
independent hand-computation of the same formula, so such a typo fails loudly
here instead of showing up as an unexplained drift in Neff/YP.
"""
import dataclasses

import numpy as np
import pytest

from primat.cache_utils import constants_hash
from primat.config import DEFAULT_PARAMS, PRIMATConfig
from primat.constants import (CONST, Constants, FROZEN_CONSTANTS,
                              OVERRIDABLE_CONSTANTS)


def test_sW2_matches_onshell_relation():
    """sin^2(theta_W) = 1/2 * (1 - sqrt(1 - 2*sqrt(2)*pi*alphaem/(GF*mZ^2)))."""
    expected = 0.5 * (1. - np.sqrt(1. - 2. * np.sqrt(2.) * np.pi * CONST.alphaem
                                    / (CONST.GF * CONST.mZ**2)))
    assert CONST.sW2 == pytest.approx(expected, rel=1e-12)
    # Sanity check against the well-known PDG ballpark value (~0.223 in the
    # MSbar scheme; the on-shell scheme used here is close but not identical).
    assert 0.20 < CONST.sW2 < 0.24


def test_effective_couplings_consistent_with_sW2():
    """geL/geR/gmuL are simple offsets of sW2 (electron/muon neutral-current couplings)."""
    assert CONST.geL == pytest.approx(0.5 + CONST.sW2, rel=1e-12)
    assert CONST.geR == pytest.approx(CONST.sW2, rel=1e-12)
    assert CONST.gmuL == pytest.approx(-0.5 + CONST.sW2, rel=1e-12)


def test_MeV_to_Kelvin_round_trips_T_weak_and_T_nucl():
    """T_weak/T_nucl are MeV_to_Kelvin scaled by their defining MeV values."""
    assert CONST.T_weak == pytest.approx(1.0 * CONST.MeV_to_Kelvin, rel=1e-12)
    assert CONST.T_nucl == pytest.approx(0.11 * CONST.MeV_to_Kelvin, rel=1e-12)


# ---------------------------------------------------------------------------
# The 16 measured constants are ordinary parameters (pass 14). Goal of this
# group: an override must reach every consumer -- the derived quantities, the
# cache key, and both backends -- and the ten exact constants must stay
# unsettable, since no config can carry them across the C ABI.
# ---------------------------------------------------------------------------

def test_every_constants_field_is_either_overridable_or_frozen():
    """The two tuples partition the dataclass: no field is both, none is neither.

    A field in neither would be silently unsettable *and* unguarded -- poking
    it would be honoured by primat/ and ignored by primat-c/.
    """
    fields = {f.name for f in dataclasses.fields(Constants)}
    assert set(OVERRIDABLE_CONSTANTS).isdisjoint(FROZEN_CONSTANTS)
    assert set(OVERRIDABLE_CONSTANTS) | set(FROZEN_CONSTANTS) == fields


def test_overridable_constants_are_params_with_the_frozen_defaults():
    """Each of the 16 is a DEFAULT_PARAMS key defaulting to CONST's value."""
    for name in OVERRIDABLE_CONSTANTS:
        assert name in DEFAULT_PARAMS, name
        assert DEFAULT_PARAMS[name] == getattr(CONST, name), name


@pytest.mark.parametrize("name, derived", [
    ("alphaem", "sW2"), ("GF", "sW2"), ("mZ", "sW2"),
    ("kappa_p", "deltakappa"), ("kappa_n", "deltakappa"),
    ("T0CMB", "n0CMB"), ("T0CMB", "eta0b"),
    ("ma", "mB"), ("He4Overma", "mB"), ("HOverma", "mB"),
])
def test_derived_quantities_follow_an_overridden_constant(name, derived):
    """Overriding a measured constant moves everything computed from it.

    They were class attributes evaluated once at import, so an override used
    to leave e.g. cfg.sW2 reporting the default while cfg.alphaem reported
    the new value -- a config inconsistent with itself.
    """
    base = PRIMATConfig()
    bumped = PRIMATConfig({name: getattr(CONST, name) * 1.01})
    assert getattr(bumped, derived) != getattr(base, derived)


def test_frozen_constants_are_rejected_on_assignment():
    """The ten exact constants cannot be overridden, by params or by poking.

    They reach the C backend only through its own compiled-in copy, so an
    override honoured by primat/ and ignored by primat-c/ would make the same
    run give two answers depending on the backend.
    """
    cfg = PRIMATConfig()
    for name in FROZEN_CONSTANTS:
        with pytest.raises(ValueError, match="exact by definition"):
            setattr(cfg, name, getattr(CONST, name) * 2.0)
    # A params dict entry is not a DEFAULT_PARAMS key: warn-and-ignore, and
    # the value must not have changed.
    with pytest.warns(UserWarning, match="unknown parameter key"):
        cfg2 = PRIMATConfig({"kB": 1.0})
    assert cfg2.kB == CONST.kB


def test_object_setattr_bypass_is_caught_by_validate_frozen_constants():
    """validate_frozen_constants is the defence-in-depth re-check.

    __setattr__ rejects the assignment where it happens; this catches the
    routes that go around it (object.__setattr__, a subclass class attribute,
    an unpickled config).
    """
    cfg = PRIMATConfig()
    object.__setattr__(cfg, "hbar", CONST.hbar * 2.0)
    with pytest.raises(ValueError, match="hbar was overridden"):
        cfg.validate_frozen_constants()


def test_constants_hash_is_per_config_not_memoised_globally():
    """Two configs with different constants must hash differently.

    constants_hash was @lru_cache(maxsize=1) over a module singleton. Once
    constants became per-config that would have served the FIRST config's hash
    to every later one -- so a run overriding gA would load, and report, the
    default-gA weak rates with no warning anywhere.
    """
    default = constants_hash(PRIMATConfig())
    assert default == constants_hash()          # cfg=None means CONST
    for name in OVERRIDABLE_CONSTANTS:
        bumped = PRIMATConfig({name: getattr(CONST, name) * 1.01})
        assert constants_hash(bumped) != default, name


def test_perturbing_a_constant_changes_the_weak_cache_filename():
    """The cache key follows the constants, so no run can read another's table.

    This is the miniature of pass 15's proof test: the filename carries the
    fingerprint hash, so a changed constant must produce a different file
    rather than silently hitting the default one.
    """
    from primat.cache_utils import fingerprint_hash
    from primat.weak_rates.cache import _weak_rate_fingerprint

    def cache_name(cfg):
        return "nTOp_" + fingerprint_hash(_weak_rate_fingerprint(cfg)) + ".txt"

    base = cache_name(PRIMATConfig())
    for name in OVERRIDABLE_CONSTANTS:
        bumped = PRIMATConfig({name: getattr(CONST, name) * 1.01})
        assert cache_name(bumped) != base, name
