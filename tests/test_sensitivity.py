# -*- coding: utf-8 -*-
"""
Tests for :mod:`primat.sensitivity` — the logarithmic-sensitivity API.

Goal of this test group: pin down the public contract that the
``notebooks/Sensitivity.ipynb`` how-to (FABLEADVICE O-10) now relies on:

* :func:`sensitivity_table` returns a :class:`SensitivityTable` of the right
  shape with a shared central (fiducial) solve;
* each of the three target flavours (multiplicative parameter, nuclear-rate
  ``delta`` rescaling, additive knob) reproduces the exact symmetric
  finite-difference number a hand-written loop would produce — this is the
  "notebook output unchanged" acceptance criterion;
* the rendering helpers (``to_markdown``/``to_dataframe``) and the input
  validation behave as documented.
"""

import math

import numpy as np
import pytest

from primat.backend import run_bbn
from primat.config import PRIMATConfig
from primat.sensitivity import sensitivity_table, SensitivityTable, SensTarget

# A single small-network fiducial reused across tests. `small` keeps each solve
# sub-second so the ~9-solve tables here stay cheap.
BASE = {"network": "small", "Omegabh2": 0.02242}
DELTA = 0.01


def _log_sens(op, om, denom):
    """Reference symmetric logarithmic finite difference used by the notebook."""
    return (math.log(op) - math.log(om)) / denom


def test_shape_and_fiducial():
    """Matrix shape, labels, and the shared central solve are wired correctly."""
    tab = sensitivity_table(
        BASE,
        observables=["YPBBN", "DoH"],
        targets=["tau_n", "Omegabh2"],
    )
    assert isinstance(tab, SensitivityTable)
    assert tab.values.shape == (2, 2)
    assert tab.row_labels == ["tau_n", "Omegabh2"]
    # Default pretty column headers come from DEFAULT_OBS_LABELS.
    assert tab.obs_labels == ["$Y_P$", "D/H"]
    # The reported fiducial values must equal an independent unperturbed solve.
    r0 = run_bbn({**BASE, "verbose": False})
    assert tab.fiducial["YPBBN"] == pytest.approx(r0["YPBBN"], rel=0, abs=0)
    assert tab.fiducial["DoH"] == pytest.approx(r0["DoH"], rel=0, abs=0)


def test_multiplicative_matches_manual():
    """A multiplicative parameter (tau_n) reproduces a hand-rolled ±1% loop."""
    cfg = PRIMATConfig(BASE)
    fid = cfg.tau_n
    rp = run_bbn({**BASE, "tau_n": fid * (1 + DELTA), "verbose": False})
    rm = run_bbn({**BASE, "tau_n": fid * (1 - DELTA), "verbose": False})
    denom = 2 * math.log1p(DELTA)
    expect = _log_sens(rp["YPBBN"], rm["YPBBN"], denom)

    tab = sensitivity_table(BASE, observables=["YPBBN"], targets=["tau_n"])
    assert tab.values[0, 0] == pytest.approx(expect, rel=1e-12)
    # tau_n physically raises YP (more free neutrons survive to freeze-out).
    assert tab.values[0, 0] > 0


def test_rate_target_matches_delta_mechanism():
    """A bare reaction name is auto-classified and uses delta_<rxn> rescaling."""
    rxn = "n_p__d_g"
    rp = run_bbn({**BASE, "rescale_nuclear_rates": True,
                  f"delta_{rxn}": +DELTA, "verbose": False})
    rm = run_bbn({**BASE, "rescale_nuclear_rates": True,
                  f"delta_{rxn}": -DELTA, "verbose": False})
    denom = 2 * math.log1p(DELTA)
    expect = _log_sens(rp["DoH"], rm["DoH"], denom)

    tab = sensitivity_table(BASE, observables=["DoH"], targets=[rxn])
    assert tab.values[0, 0] == pytest.approx(expect, rel=1e-12)


def test_additive_target_matches_manual():
    """An additive knob (DeltaNeff, fiducial 0) varies by ±step, not ±p*step."""
    step = 1.0
    rp = run_bbn({**BASE, "DeltaNeff": +step, "verbose": False})
    rm = run_bbn({**BASE, "DeltaNeff": -step, "verbose": False})
    denom = 2 * math.log1p(DELTA)
    expect = _log_sens(rp["YPBBN"], rm["YPBBN"], denom)

    tab = sensitivity_table(
        BASE,
        observables=["YPBBN"],
        targets=[SensTarget("DeltaNeff", kind="additive", step=step)],
    )
    assert tab.values[0, 0] == pytest.approx(expect, rel=1e-12)
    # Extra relativistic species speed up expansion -> earlier freeze-out ->
    # more neutrons -> larger YP.
    assert tab.values[0, 0] > 0


def test_zero_fiducial_multiplicative_raises():
    """Multiplicative step on a zero-fiducial parameter is a clear error."""
    with pytest.raises(ValueError, match="fiducial value 0"):
        sensitivity_table(BASE, observables=["YPBBN"],
                          targets=[SensTarget("DeltaNeff", kind="param")])


def test_unknown_observable_raises():
    """Asking for a non-existent observable key fails fast and informatively."""
    with pytest.raises(KeyError, match="not in the result dict"):
        sensitivity_table(BASE, observables=["NoSuchObs"], targets=["tau_n"])


def test_to_markdown_and_dataframe():
    """The two rendering views agree with the underlying matrix."""
    tab = sensitivity_table(BASE, observables=["YPBBN", "DoH"],
                            targets=["tau_n", "n_p__d_g"])
    md = tab.to_markdown()
    # Header + separator + one row per target.
    lines = md.splitlines()
    assert lines[0].startswith("| Parameter |")
    assert len(lines) == 2 + len(tab.row_labels)
    assert "tau_n" in md and "n_p__d_g" in md

    df = tab.to_dataframe()
    assert list(df.index) == tab.row_labels
    assert df.shape == (2, 2)
    np.testing.assert_allclose(df.to_numpy(), tab.values)


def test_custom_obs_labels_length_check():
    """Mismatched obs_labels length is rejected."""
    with pytest.raises(ValueError, match="obs_labels"):
        sensitivity_table(BASE, observables=["YPBBN", "DoH"],
                          targets=["tau_n"], obs_labels=["only-one"])
