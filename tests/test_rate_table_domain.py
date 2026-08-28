# -*- coding: utf-8 -*-
"""A degenerate rate table read from disk must be refused, not integrated.

Goal: the library door agrees with the GUI door (``tests/test_gui_robustness.py``)
on what a valid rate table is, and both backends refuse the same tables with the
same message. The two failures this guards against: an unsorted ``T9`` column
supplied through a ``user_nuclear_dir`` overlay once ran to completion and
reported YP = 0.00000741 (Python) / 0.00000711 (C) against a true 0.24699907,
and an empty table file crashed the C backend with SIGSEGV.

The domain rules are pinned against ``validate_rate_table`` directly (fast); the
two end-to-end tests confirm they actually fire on the real loader path, on the
backend under test.
"""
import os

import numpy as np
import pytest

from primat.backend import HAS_C_BACKEND, run_bbn
from primat.config import PRIMATConfig
from primat.network_data import load_network, validate_rate_table


# A legitimate table: increasing T9, non-negative rates, three columns.
VALID = "1e-3 1e-3 1.0\n1e-2 2e-3 1.0\n1e-1 3e-3 1.0\n1.0 4e-3 1.0\n10.0 5e-3 1.0\n"


def _overlay(tmp_path, text, reaction="n_p__d_g"):
    """A ``user_nuclear_dir`` overlay whose only content is one rate table."""
    folder = tmp_path / "tables" / reaction
    folder.mkdir(parents=True)
    (folder / f"{reaction}_primat.txt").write_text(text)
    return str(tmp_path)


# ---------------------------------------------------------------------------
# The rules themselves
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,T9,rate,expected", [
    ("unsorted T9",   [1.0, 1e-3, 10.0], [1.0, 2.0, 3.0], "goes backwards"),
    ("duplicate T9",  [1e-3, 1e-3, 1.0], [1.0, 2.0, 3.0], "repeats"),
    ("zero T9",       [0.0, 1e-2, 1.0],  [1.0, 2.0, 3.0], "strictly positive"),
    ("negative T9",   [-1.0, 1e-2, 1.0], [1.0, 2.0, 3.0], "strictly positive"),
    ("negative rate", [1e-3, 1e-2, 1.0], [1.0, -2.0, 3.0], "must not be negative"),
    ("NaN rate",      [1e-3, 1e-2, 1.0], [1.0, np.nan, 3.0], "non-finite"),
    ("inf T9",        [1e-3, np.inf, 1.0], [1.0, 2.0, 3.0], "non-finite"),
    ("single row",    [1.0],             [1.0],            "at least two"),
    ("no rows",       [],                [],               "at least two"),
])
def test_degenerate_table_is_rejected(label, T9, rate, expected):
    """Each rule the log-log resampler depends on is enforced, naming the row."""
    with pytest.raises(ValueError, match=expected):
        validate_rate_table(np.array(T9, dtype=float),
                             np.array(rate, dtype=float))


def test_valid_table_passes():
    """The check rejects only degenerate input. A zero rate is legitimate."""
    validate_rate_table(np.array([1e-3, 1e-2, 1.0]),
                         np.array([0.0, 2.0, 3.0]),
                         np.array([1.0, 1.0, 1.0]))


def test_message_names_the_source_and_the_row():
    """The message carries the file it came from and the 1-based data row."""
    with pytest.raises(ValueError, match=r"my_table\.txt.*data row 3"):
        validate_rate_table(np.array([1e-3, 1e-2, 1.0]),
                             np.array([1.0, 2.0, -3.0]),
                             source="my_table.txt")


def test_every_shipped_table_passes(tmp_path):
    """All 402 shipped rate tables satisfy the rules the loader now enforces.

    Turning the check on without this would have broken the default run.
    """
    import primat
    root = os.path.join(os.path.dirname(primat.__file__),
                         "data", "nuclear", "tables")
    checked = 0
    for folder in sorted(os.listdir(root)):
        d = os.path.join(root, folder)
        if not os.path.isdir(d):
            continue
        for name in sorted(os.listdir(d)):
            data = np.loadtxt(os.path.join(d, name), unpack=True)
            validate_rate_table(data[0], data[1], data[2],
                                 source=f"{folder}/{name}")
            checked += 1
    assert checked >= 400, f"only {checked} tables found under {root}"


# ---------------------------------------------------------------------------
# The loader path
# ---------------------------------------------------------------------------

def test_two_column_table_loads(tmp_path):
    """A two-column table is the GUI's documented minimum and must load here too.

    It used to escape as ``IndexError: index 2 is out of bounds for axis 0 with
    size 2``, so a table exported from the GUI could not be re-used as a
    ``user_nuclear_dir`` overlay. The missing column reads as zero uncertainty.
    """
    overlay = _overlay(tmp_path, "1e-3 1e-3\n1e-2 2e-3\n1e-1 3e-3\n1.0 4e-3\n")
    cfg = PRIMATConfig({"network": "small", "user_nuclear_dir": overlay,
                         "verbose": False})
    net = load_network(cfg, era="LT")
    # The rate arrays skip the prepended weak ``n__p`` entry that ``names`` carries.
    i = net.names.index("n_p__d_g") - 1
    assert np.all(net._expsigma[i] == 0.0)
    assert net._fwd_median[i][0] == pytest.approx(1e-3, rel=1e-3)


@pytest.mark.parametrize("label,text,expected", [
    ("empty file", "", "no data rows"),
    ("comments only", "# nothing here\n", "no data rows"),
    ("four columns", "1e-3 1 1 9\n1e-2 2 1 9\n1e-1 3 1 9\n", "2 or 3 columns"),
    ("single row", "1.0 4e-3 1.0\n", "at least two"),
])
def test_malformed_table_file_names_itself(tmp_path, label, text, expected):
    """A malformed file is refused by a message naming it, not by an IndexError.

    Every one of these used to escape as a bare ``IndexError`` from the column
    indexing, which mentions neither the file nor the reaction.
    """
    overlay = _overlay(tmp_path, text)
    cfg = PRIMATConfig({"network": "small", "user_nuclear_dir": overlay,
                         "verbose": False})
    with pytest.raises(ValueError, match=expected):
        load_network(cfg, era="LT")


@pytest.mark.parametrize("backend", ["python",
    pytest.param("c", marks=pytest.mark.skipif(not HAS_C_BACKEND,
                                                reason="C extension not built"))])
def test_unsorted_overlay_table_is_refused_by_both_backends(tmp_path, backend):
    """The headline: this ran to completion on both backends, YP ~ 7e-06.

    The wording is shared, so the same table produces the same sentence
    whichever backend serves the run.
    """
    overlay = _overlay(tmp_path, "1.0 4e-3 1.0\n1e-3 1e-3 1.0\n10.0 5e-3 1.0\n")
    with pytest.raises((ValueError, RuntimeError),
                        match="T9 column must increase strictly"):
        run_bbn(params={"network": "small", "user_nuclear_dir": overlay,
                         "verbose": False}, force_backend=backend)


@pytest.mark.parametrize("backend", ["python",
    pytest.param("c", marks=pytest.mark.skipif(not HAS_C_BACKEND,
                                                reason="C extension not built"))])
def test_empty_overlay_table_is_refused_by_both_backends(tmp_path, backend):
    """An empty table file crashed the C backend with SIGSEGV.

    ``cpr_table_read`` treated "no data rows" as an error only when the column
    count was also unknown, and rate tables are read with a column hint, so the
    resampler indexed row 0 of an empty column.
    """
    overlay = _overlay(tmp_path, "")
    with pytest.raises((ValueError, RuntimeError), match="no data rows"):
        run_bbn(params={"network": "small", "user_nuclear_dir": overlay,
                         "verbose": False}, force_backend=backend)


# ---------------------------------------------------------------------------
# The loader's other user-supplied inputs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("backend", ["python",
    pytest.param("c", marks=pytest.mark.skipif(not HAS_C_BACKEND,
                                                reason="C extension not built"))])
def test_empty_network_file_is_refused_by_both_backends(tmp_path, backend):
    """A network with no reaction cannot nucleosynthesise and must say so.

    Python failed with ``IndexError: too many indices for array`` from the
    reverse-rate cap; C ran to completion and reported YP = 0.00000000 and
    He3/He4 = nan with exit status 0.
    """
    (tmp_path / "networks").mkdir()
    (tmp_path / "networks" / "empty_net.txt").write_text("\n# nothing here\n")
    with pytest.raises((ValueError, RuntimeError), match="lists no reactions"):
        run_bbn(params={"network": "empty_net",
                         "user_nuclear_dir": str(tmp_path), "verbose": False},
                 force_backend=backend)


@pytest.mark.parametrize("backend", ["python",
    pytest.param("c", marks=pytest.mark.skipif(not HAS_C_BACKEND,
                                                reason="C extension not built"))])
def test_amax_dropping_every_reaction_is_refused(tmp_path, backend):
    """The same empty network reached by filtering rather than by an empty file."""
    with pytest.raises((ValueError, RuntimeError), match="drops every reaction"):
        run_bbn(params={"network": "large", "amax": 1, "verbose": False},
                 force_backend=backend)


def test_malformed_decay_row_names_the_file_and_line(tmp_path):
    """A short or non-numeric ``decays.txt`` row used to raise a bare IndexError.

    The C loader has always named the file and line
    (``primat-c/src/network_data.c``'s ``cannot parse decay row``); this is
    Python catching up.
    """
    from primat.network_data import _load_decay_table
    (tmp_path / "decays.txt").write_text(
        "# name halflife_s rate_s uncertainty ref\nBm_H3 3.887e+08 1.784e-09\n")
    with pytest.raises(ValueError, match=r"decays\.txt:2: cannot parse decay row"):
        _load_decay_table(str(tmp_path))
