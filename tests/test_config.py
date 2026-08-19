"""Tests for PRIMATConfig: defaults, overrides, derived quantities, nuclide data."""
import os
import warnings
import pytest
from primat.config import (PRIMATConfig, DEFAULT_PARAMS, PARAM_GROUPS,
                            _default_params_comments)


def test_param_groups_covers_every_default_params_key_exactly_once():
    """PARAM_GROUPS is the single source of truth the GUI, CLI --list-params,
    and the param-template generator all derive their section headings from
    (see primat/tools/gen_param_templates.py, which derives both templates
    from it); it must stay
    exhaustive and non-overlapping as DEFAULT_PARAMS keys are added, removed,
    or renamed, or those consumers would silently drop or duplicate a key."""
    seen = []
    for group, keys in PARAM_GROUPS.items():
        seen.extend(keys)
    assert len(seen) == len(set(seen)), "a key appears in more than one PARAM_GROUPS group"
    assert set(seen) == set(DEFAULT_PARAMS), (
        set(DEFAULT_PARAMS) - set(seen), set(seen) - set(DEFAULT_PARAMS))


def test_config_type_annotation_block_is_up_to_date():
    """The TYPE_CHECKING attribute-annotation block inside PRIMATConfig
    (between the BEGIN/END GENERATED PARAM ANNOTATIONS sentinels) exists so
    IDEs/mypy see real completions for every DEFAULT_PARAMS key despite the
    dynamic __getattr__ used for p_*/delta_* rate variations. It is
    generated text, not hand-maintained -- this test fails loudly if a
    DEFAULT_PARAMS key was added/removed/renamed without regenerating it."""
    import inspect
    import primat.config as config_module

    source = inspect.getsource(config_module)
    begin = "# BEGIN GENERATED PARAM ANNOTATIONS"
    end = "# END GENERATED PARAM ANNOTATIONS"
    start_idx = source.index(begin) + len(begin)
    end_idx = source.index(end)
    block = source[start_idx:end_idx]

    expected = config_module._generate_config_type_annotations()
    for line in expected.splitlines():
        assert line in block, f"missing/stale annotation line: {line!r}"


def test_dir_includes_default_params_and_dynamic_rate_keys():
    """dir(cfg) should offer every DEFAULT_PARAMS key (ordinary instance
    attributes, already visible via object.__dir__) plus any p_*/delta_*
    rate-variation attribute actually set on the instance (only reachable
    via __getattr__/__setattr__, so object.__dir__ alone would miss it)."""
    cfg = PRIMATConfig({"p_n_p__d_g": 0.5, "delta_d_p__He3_g": 1.0})
    listing = dir(cfg)
    for key in DEFAULT_PARAMS:
        assert key in listing
    assert "p_n_p__d_g" in listing
    assert "delta_d_p__He3_g" in listing


def test_default_construction():
    """A no-argument PRIMATConfig is usable: the headline defaults are present
    and sane."""
    cfg = PRIMATConfig()
    assert cfg.Omegabh2 > 0
    assert cfg.is_small is True
    assert cfg.numerical_precision > 0


def test_user_override():
    """A user dict overrides the corresponding default."""
    cfg = PRIMATConfig({"Omegabh2": 0.020})
    assert cfg.Omegabh2 == pytest.approx(0.020)


def test_unknown_key_warns():
    """An unrecognised key warns, naming the offending key so a typo is
    findable."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        PRIMATConfig({"not_a_real_param": 42})
    assert any("not_a_real_param" in str(x.message) for x in w)


def test_unknown_key_does_not_raise():
    """...but it does not raise: strict_params defaults to False, i.e.
    'warn and ignore' (see test_strict_params_* for the opt-in strict mode)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        PRIMATConfig({"totally_unknown": 99})


def test_omegabh2_to_eta0b_positive():
    """The Omegabh2 -> eta0b conversion factor is positive and available."""
    cfg = PRIMATConfig()
    assert cfg.Omegabh2_to_eta0b > 0


def test_nuclides_keys():
    """The Nuclides table covers at least every species the shipped networks
    evolve."""
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
    """A delta_<rxn> override naming a reaction outside the network warns
    (the additive-variation twin of test_p_rxn_typo_warns above)."""
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
    """Nuclides maps each name to its [N, Z] pair -- spot-checked against the
    values every A/Z conservation check depends on."""
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
    """Every physical constant exposed on the config is positive -- catches a
    unit conversion or a default that got zeroed out."""
    cfg = PRIMATConfig()
    for attr in ("me", "mn", "mp", "Mpl", "kB", "MeV"):
        assert getattr(cfg, attr) > 0, f"cfg.{attr} should be positive"


def test_bad_type_raises_typeerror():
    """A value of the wrong *type* raises an immediate, self-explanatory
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
    """Physical/numerical range violations raise ValueError, always
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
    """An unknown key that is a near-miss of a real one gets a
    difflib 'did you mean ...?' suggestion in the warning."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        PRIMATConfig({"Omegab2h": 0.022})
    msgs = " ".join(str(x.message) for x in w)
    assert "Omegab2h" in msgs
    assert "Omegabh2" in msgs  # the suggested close match


def test_strict_params_raises_on_unknown_key():
    """strict_params=True upgrades the unknown-key warning to a
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


def test_strict_params_raises_on_unmatched_rate_variation():
    """strict_params=True upgrades the p_<rxn>/delta_<rxn> typo warning to a
    ValueError too.

    GOAL: a mistyped reaction name is the archetypal silent no-op strict mode
    exists to catch -- reaction names are long and underscore-heavy, and the
    run otherwise proceeds with that rate unvaried.
    """
    with pytest.raises(ValueError) as exc:
        PRIMATConfig({"p_n_p__d_gg": 1.0, "strict_params": True})
    assert "p_n_p__d_gg" in str(exc.value)
    assert "strict_params=True" in str(exc.value)


def test_cross_field_ranges_rejected():
    """Inconsistent *pairs* of parameters are rejected at construction.

    GOAL: each of these passes the per-key range checks yet leaves the solver
    with an impossible request, and each used to surface only as an opaque
    "step size underflowed" from inside the MT integration (or, for the MC cap,
    silently pinned every sampled rate to `cap`).
    """
    cases = [
        ({"rate_grid_T9_min": 10.0, "rate_grid_T9_max": 1e-3}, "rate_grid_T9_min"),
        ({"T_end_MeV": 100.0}, "T_end_MeV"),                  # > T_start_cosmo_MeV = 40
        ({"mc_rate_rescale_cap": 0.5}, "mc_rate_rescale_cap"),
        # Q = mn - mp at or below me: sqrt(E^2 - me^2) is imaginary over the
        # whole n<->p integration range. Rejecting it here is what keeps the C
        # backend out of a NaN-fed adaptive quadrature that recurses to its
        # full depth and takes hours per integral.
        ({"mp": 939.2101}, "mn - mp"),                        # Q = 0.355 MeV
        ({"mn": 938.6}, "mn - mp"),                           # Q < 0
        ({"me": 2.0}, "mn - mp"),                             # me above Q
    ]
    for params, expected in cases:
        with pytest.raises(ValueError) as exc:
            PRIMATConfig(params)
        assert expected in str(exc.value)
    # The boundary is Q = me: just below it the range collapses (and the tau_n
    # normalisation divides by zero), just above it the run is legal.
    d = PRIMATConfig()
    with pytest.raises(ValueError, match="mn - mp"):
        PRIMATConfig({"mp": d.mn - d.me * 0.999})
    assert PRIMATConfig({"mp": d.mn - d.me * 1.001}).mn - d.me * 1.001 > 0
    # A cap of exactly 1 ("no variation") and None ("no cap") stay legal.
    assert PRIMATConfig({"mc_rate_rescale_cap": 1.0}).mc_rate_rescale_cap == 1.0
    assert PRIMATConfig({"mc_rate_rescale_cap": None}).mc_rate_rescale_cap is None


def test_data_dir_validated_before_nuclides_are_read(tmp_path):
    """A bad data_dir is reported as such, not as a missing nuclides.csv.

    GOAL: data_dir is consumed by _load_nuclide_data during __init__, so the
    check has to run *before* it -- otherwise the user sees a bare
    FileNotFoundError on '<data_dir>/csv/nuclides.csv' naming neither the
    parameter nor what the directory should contain.
    """
    with pytest.raises(ValueError, match="not an existing directory"):
        PRIMATConfig({"data_dir": str(tmp_path / "nope")})

    # Exists, but is not a data tree: equally wrong for a *full takeover*.
    with pytest.raises(ValueError) as exc:
        PRIMATConfig({"data_dir": str(tmp_path)})
    assert "does not look like a primat data tree" in str(exc.value)
    assert "csv/" in str(exc.value) and "nuclear/" in str(exc.value)


def test_cache_dir_expands_tilde_like_the_c_backend(monkeypatch, tmp_path):
    """cache_dir is a path parameter: '~' must expand on the Python side too.

    GOAL: primat-c's cpr_is_path_field expands cache_dir, so a config with
    cache_dir='~/primat-cache' would otherwise have the two backends writing
    their caches to different places (a literal './~/primat-cache' on the
    Python side) and silently stop sharing them.
    """
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))   # Windows equivalent
    cfg = PRIMATConfig({"cache_dir": os.path.join("~", "primat-cache")})
    assert cfg.cache_dir == str(tmp_path / "primat-cache")


def test_every_default_param_has_a_one_line_description():
    """`primat --list-params` must describe every parameter.

    GOAL: --set is deliberately hidden from --help, which makes --list-params
    the documented way to discover parameters; a key printed as a bare
    "name = value" is undocumented as far as a CLI user is concerned. The
    description comes from config.py's own comments (trailing, else the first
    prose line of the comment block above the key), so this also enforces that
    a newly added parameter is commented at all.
    """
    comments = _default_params_comments()
    undocumented = sorted(k for k in DEFAULT_PARAMS if not comments.get(k))
    assert not undocumented, (
        "DEFAULT_PARAMS keys with no description for --list-params: "
        f"{undocumented}. Add a trailing '# ...' comment on the key's own line "
        "(or a comment block immediately above it) in primat/config.py."
    )


def test_config_dynamic_rate_attrs():
    """Dynamic ``p_*`` and ``delta_*`` attrs round-trip through the backing dicts."""
    cfg = PRIMATConfig()

    # The backing dicts themselves start with those very prefixes, so
    # assigning them must NOT be read as a variation of a reaction named
    # "rxn" (which would die in float(dict)).
    cfg.p_rxn = dict(cfg.p_rxn)
    cfg.delta_rxn = dict(cfg.delta_rxn)
    assert isinstance(cfg.p_rxn, dict) and isinstance(cfg.delta_rxn, dict)

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
