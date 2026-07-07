"""Tests for PRIMATConfig: defaults, overrides, derived quantities, nuclide data."""
import warnings
import pytest
from primat.config import PRIMATConfig, DEFAULT_PARAMS


def test_default_construction():
    cfg = PRIMATConfig()
    assert cfg.Omegabh2 > 0
    assert cfg.is_small is True
    assert cfg.numerical_precision > 0


def test_user_override():
    cfg = PRIMATConfig({"Omegabh2": 0.020})
    assert cfg.Omegabh2 == pytest.approx(0.020)


def test_unknown_key_warns():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        PRIMATConfig({"not_a_real_param": 42})
    assert any("not_a_real_param" in str(x.message) for x in w)


def test_unknown_key_does_not_raise():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        PRIMATConfig({"totally_unknown": 99})


def test_omegabh2_to_eta0b_positive():
    cfg = PRIMATConfig()
    assert cfg.Omegabh2_to_eta0b > 0


def test_nuclides_keys():
    cfg = PRIMATConfig()
    expected_subset = {"n", "p", "H2", "H3", "He3", "He4", "He6",
                       "Li6", "Li7", "Be7", "Li8", "B8"}
    assert expected_subset.issubset(set(cfg.Nuclides.keys()))


def test_p_rxn_typo_warns():
    """A p_<rxn> override whose reaction name isn't in the network must warn.

    Before this check existed, a typo'd reaction name (e.g. a stray
    underscore, or a name from a different network) was silently accepted: it
    became a no-op dict entry in cfg.p_rxn with no signal that the rate
    variation the caller asked for was never actually applied.
    """
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        PRIMATConfig({"p_not_a_real_reaction": 0.5})
    assert any("p_not_a_real_reaction" in str(x.message) for x in w)


def test_delta_rxn_typo_warns():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        PRIMATConfig({"delta_not_a_real_reaction": 0.5})
    assert any("delta_not_a_real_reaction" in str(x.message) for x in w)


def test_p_rxn_valid_reaction_does_not_warn():
    """A genuine reaction name (from the small network) must not warn."""
    cfg = PRIMATConfig()
    rxn = next(iter(cfg.p_rxn))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        cfg2 = PRIMATConfig({f"p_{rxn}": 0.3})
    assert not any(rxn in str(x.message) for x in w)
    assert getattr(cfg2, f"p_{rxn}") == pytest.approx(0.3)


def test_nuclides_NZ_values():
    cfg = PRIMATConfig()
    assert cfg.Nuclides["He4"] == [2, 2]
    assert cfg.Nuclides["H2"]  == [1, 1]
    assert cfg.Nuclides["Li7"] == [4, 3]
    assert cfg.Nuclides["n"]   == [1, 0]
    assert cfg.Nuclides["p"]   == [0, 1]


def test_p_rate_keys_count():
    """``p_rxn``/``delta_rxn`` carry one MCMC weight per *configured*
    network's reaction (small/large±amax), not always the full large set --
    see the "Corrected bug on MC" fix in ``PRIMATConfig.__init__``, which reads
    ``load_reaction_names(self._resolved_data_dir, self.network)`` rather than the
    hardcoded ``_REACTIONS_LARGE`` list."""
    from primat.network_data import load_network

    cfg = PRIMATConfig()  # default network="small" -> 12 reactions
    assert cfg.network == "small"
    assert len(cfg.p_rxn) == 12
    assert len(cfg.delta_rxn) == 12

    cfg_amax8 = PRIMATConfig({"network": "large", "amax": 8})
    net_amax8 = load_network(cfg_amax8, era="LT")
    n_thermonuclear = len(net_amax8.names) - 1  # exclude the prepended n__p
    assert len(cfg_amax8.p_rxn) == n_thermonuclear == 67
    assert len(cfg_amax8.delta_rxn) == n_thermonuclear == 67


def test_physical_constants_positive():
    cfg = PRIMATConfig()
    for attr in ("me", "mn", "mp", "Mpl", "kB", "MeV"):
        assert getattr(cfg, attr) > 0, f"cfg.{attr} should be positive"


def test_bad_type_raises_typeerror():
    """O-1: a value of the wrong *type* raises an immediate, self-explanatory
    TypeError naming the key, value, and expected type -- instead of dying much
    later inside the thermodynamics (Omegabh2="0.022" once produced a cryptic
    "can't multiply sequence by non-int" from deep in the solver)."""
    cases = [
        {"Omegabh2": "0.022"},   # str for a float field
        {"verbose": 1.5},         # non-bool for a bool field
        {"amax": "eight"},        # str for an int/None field
        {"network": 5},           # int for a str field
    ]
    for params in cases:
        with pytest.raises(TypeError) as exc:
            PRIMATConfig(params)
        key = next(iter(params))
        assert key in str(exc.value)


def test_bool_not_accepted_as_number():
    """True/False must not be silently taken as 1.0/0.0 for a numeric field:
    passing a bool where a float is expected is a bug, not the number one."""
    with pytest.raises(TypeError):
        PRIMATConfig({"Omegabh2": True})


def test_out_of_range_raises_valueerror():
    """O-1: physical/numerical range violations raise ValueError, always
    (independent of strict_params)."""
    cases = [
        {"Omegabh2": -0.1},
        {"tau_n": 0.0},
        {"numerical_precision": 0.0},
        {"rate_grid_npts": 0},
        {"h": -0.5},
        {"std_tau_n": -1.0},
    ]
    for params in cases:
        with pytest.raises(ValueError) as exc:
            PRIMATConfig(params)
        assert "out of range" in str(exc.value)


def test_std_tau_n_zero_allowed():
    """std_tau_n is a 1-sigma width that may legitimately be exactly 0."""
    cfg = PRIMATConfig({"std_tau_n": 0.0})
    assert cfg.std_tau_n == 0.0


def test_rate_interp_order_choice():
    """rate_interp_order is constrained to linear/quadratic/cubic."""
    for good in ("linear", "quadratic", "cubic"):
        assert PRIMATConfig({"rate_interp_order": good}).rate_interp_order == good
    with pytest.raises(ValueError) as exc:
        PRIMATConfig({"rate_interp_order": "spline"})
    assert "one of" in str(exc.value)


def test_numpy_scalars_accepted():
    """numpy scalar overrides (np.float64/np.int64, common in MCMC drivers)
    must pass the type check just like their Python counterparts."""
    import numpy as np
    cfg = PRIMATConfig({"Omegabh2": np.float64(0.022), "amax": np.int64(8)})
    assert cfg.Omegabh2 == pytest.approx(0.022)
    assert cfg.amax == 8


def test_nullable_params_accept_none():
    """None-able parameters (path sentinels, amax, the MC rate cap) accept
    None without a type error."""
    cfg = PRIMATConfig({"amax": None, "output_file": None,
                        "mc_rate_rescale_cap": None, "data_dir": None})
    assert cfg.amax is None
    assert cfg.output_file is None
    assert cfg.mc_rate_rescale_cap is None


def test_unknown_key_suggests_close_match():
    """O-1: an unknown key that is a near-miss of a real one gets a
    difflib 'did you mean ...?' suggestion in the warning."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        PRIMATConfig({"Omegab2h": 0.022})
    msgs = " ".join(str(x.message) for x in w)
    assert "Omegab2h" in msgs
    assert "Omegabh2" in msgs  # the suggested close match


def test_strict_params_raises_on_unknown_key():
    """O-1: strict_params=True upgrades the unknown-key warning to a
    ValueError (recommended in scripted/MCMC pipelines)."""
    with pytest.raises(ValueError) as exc:
        PRIMATConfig({"Omegab2h": 0.022, "strict_params": True})
    assert "Omegab2h" in str(exc.value)
    assert "Omegabh2" in str(exc.value)  # suggestion still surfaced


def test_strict_params_default_false_only_warns():
    """With strict_params at its default (False) an unknown key must warn,
    not raise (back-compatible behaviour)."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        PRIMATConfig({"Omegab2h": 0.022})
    assert any("Omegab2h" in str(x.message) for x in w)


def test_strict_params_is_a_default_param():
    """strict_params must be a real DEFAULT_PARAMS key (three-file sync)."""
    assert DEFAULT_PARAMS["strict_params"] is False


def test_config_dynamic_rate_attrs():
    """Dynamic ``p_*`` and ``delta_*`` attrs round-trip through the backing dicts."""
    cfg = PRIMATConfig()

    # p_<reaction> attribute routes to cfg.p_rxn dict
    cfg.p_n_p__d_g = 0.5
    assert cfg.p_rxn["n_p__d_g"] == 0.5
    assert cfg.p_n_p__d_g == 0.5

    # delta_<reaction> routes to cfg.delta_rxn dict
    cfg.delta_d_p__He3_g = 0.1
    assert cfg.delta_rxn["d_p__He3_g"] == 0.1
    assert cfg.delta_d_p__He3_g == 0.1

    # Unknown prefix falls through to object.__setattr__
    cfg.some_random_param = 42
    assert cfg.some_random_param == 42
