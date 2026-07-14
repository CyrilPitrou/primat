# -*- coding: utf-8 -*-
"""
test_evolution.py
==================
Unit tests for ``primat.evolution`` (the unified time-evolution schema,
``PRIMAT.md`` S7), independent of any actual BBN solve: round-trips a
synthetic ``EvolutionResult`` through :func:`dump_evolution`/
:func:`load_evolution`, both to a string and to a real file on disk.
"""
import numpy as np
import pytest

from primat.evolution import EvolutionResult, dump_evolution, load_evolution


def _make_result():
    t = np.array([1.0, 2.0, 3.0])
    return EvolutionResult(
        t=t,
        a=np.array([0.1, 0.2, 0.3]),
        T_gamma=np.array([10.0, 5.0, 1.0]),
        T_nu={"e": np.array([9.0, 4.5, 0.9]),
              "mu": np.array([9.1, 4.6, 1.0]),
              "tau": np.array([9.2, 4.7, 1.1])},
        Y={"n": np.array([1.0, 0.9, 0.8]), "p": np.array([0.0, 0.1, 0.2])},
    )


def test_dump_evolution_header_and_round_trip(tmp_path):
    result = _make_result()
    path = tmp_path / "evolution.tsv"

    text = dump_evolution(result, path=str(path))

    header = text.splitlines()[0].split("\t")
    assert header == ["t_s", "a", "T_gamma_MeV", "T_nue_MeV", "T_numu_MeV",
                       "T_nutau_MeV", "Y_n", "Y_p"]
    assert path.exists()
    assert path.read_text() == text

    loaded = load_evolution(str(path))
    np.testing.assert_allclose(loaded.t, result.t)
    np.testing.assert_allclose(loaded.a, result.a)
    np.testing.assert_allclose(loaded.T_gamma, result.T_gamma)
    for flavour in ("e", "mu", "tau"):
        np.testing.assert_allclose(loaded.T_nu[flavour], result.T_nu[flavour])
    for species in ("n", "p"):
        np.testing.assert_allclose(loaded.Y[species], result.Y[species])


def test_dump_evolution_without_path_writes_nothing(tmp_path):
    result = _make_result()
    text = dump_evolution(result)
    assert isinstance(text, str)
    assert list(tmp_path.iterdir()) == []


def test_dump_evolution_rates_columns_appended_after_Y(tmp_path):
    """EvolutionResult.rates (per-reaction forward-rate columns) are
    serialised AFTER the Y_ block and survive a synthetic dump/load round-trip
    (no solve needed -- exercises just the evolution.py schema plumbing)."""
    result = _make_result()
    result.rates = {"n_p__d_g_frwrd": np.array([1.0, 2.0, 3.0]),
                    "d_p__He3_g_frwrd": np.array([4.0, 5.0, 6.0])}
    path = tmp_path / "evo_rates.tsv"
    text = dump_evolution(result, path=str(path))

    header = text.splitlines()[0].split("\t")
    # core + Y block + rate block, in that order.
    assert header == ["t_s", "a", "T_gamma_MeV", "T_nue_MeV", "T_numu_MeV",
                      "T_nutau_MeV", "Y_n", "Y_p",
                      "n_p__d_g_frwrd", "d_p__He3_g_frwrd"]
    loaded = load_evolution(str(path))
    assert set(loaded.rates) == set(result.rates)
    for k in result.rates:
        np.testing.assert_allclose(loaded.rates[k], result.rates[k])


def test_load_evolution_rates_none_when_no_rate_columns(tmp_path):
    """A file with no *_frwrd/*_bkwrd columns loads with rates=None (the flag
    was off), so the default schema stays rate-free."""
    path = tmp_path / "evo.tsv"
    dump_evolution(_make_result(), path=str(path))
    loaded = load_evolution(str(path))
    assert loaded.rates is None


@pytest.mark.solve
def test_evolution_rates_columns_roundtrip(tmp_path):
    """output_rates_time_evolution=True adds <reaction>_frwrd columns after
    the Y_ block (small network), populated in EvolutionResult.rates, and they
    survive a dump/load round-trip."""
    from primat.backend import run_bbn
    out = tmp_path / "evo_rates.tsv"
    r = run_bbn({"network": "small", "output_time_evolution": True,
                 "output_rates_time_evolution": True,
                 "output_file": str(out)}, force_backend="python")
    evo = r["evolution"]
    assert evo.rates, "EvolutionResult.rates is empty with the flag on"
    assert any(k.endswith("_frwrd") for k in evo.rates)
    header = out.read_text().splitlines()[0].split("\t")
    y_cols = [c for c in header if c.startswith("Y_")]
    rate_cols = [c for c in header if c.endswith(("_frwrd", "_bkwrd"))]
    assert rate_cols and header.index(rate_cols[0]) > header.index(y_cols[-1])
    loaded = load_evolution(str(out))
    assert set(loaded.rates) == set(evo.rates)


@pytest.mark.solve
def test_evolution_schema_unchanged_when_flag_off(tmp_path):
    """The default schema must not grow columns (flag off = today's header)."""
    from primat.backend import run_bbn
    out = tmp_path / "evo.tsv"
    run_bbn({"network": "small", "output_time_evolution": True,
             "output_file": str(out)}, force_backend="python")
    header = out.read_text().splitlines()[0].split("\t")
    assert header[:6] == ["t_s", "a", "T_gamma_MeV", "T_nue_MeV",
                          "T_numu_MeV", "T_nutau_MeV"]
    assert not [c for c in header if c.endswith(("_frwrd", "_bkwrd"))]
