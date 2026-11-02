"""
Duplicate reaction channels: two reactions with the same reactants and
products (in any order or direction) add their rates, so loading such a
network warns. A name ending in ``___<label>`` declares a deliberate extra
channel: it is never reported, and it borrows its base reaction's
stoichiometry, detailed balance, nuclear-QED correction and table folder.
"""
import warnings

import numpy as np
import pytest

from primat.backend import HAS_C_BACKEND, run_bbn
from primat.config import PRIMATConfig
from primat.network_data import load_network

_BASE = ["n_p__d_g", "d_p__He3_g", "d_d__He3_n"]
_VARIANT = "d_p__He3_g___alt, d_p__He3_g_primat.txt"


def _duplicate_warnings(record):
    return [str(w.message) for w in record if "same channel" in str(w.message)]


def _table(T9):
    """A two-column rate table as GUI-upload text."""
    return np.column_stack([T9, 1e-30 * T9])


@pytest.mark.parametrize("network", ["small", "small_parthenope", "large"])
def test_shipped_networks_have_no_duplicate_channel(network):
    """No shipped network lists a channel twice."""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        load_network(PRIMATConfig({"network": network}))
    assert _duplicate_warnings(record) == []


def test_reversed_duplicate_warns_once_naming_both():
    """d+p->He3+g added again as He3+g->d+p warns once, naming both."""
    T9 = np.logspace(-3, 1, 50)
    tab = _table(T9)
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        load_network(PRIMATConfig({"network": "small"}),
                     reaction_names=_BASE + ["He3_g__d_p"],
                     custom_tables={"He3_g__d_p": (tab[:, 0], tab[:, 1],
                                                   np.full_like(T9, 1.1))})
    found = _duplicate_warnings(record)
    assert len(found) == 1
    assert "d_p__He3_g, He3_g__d_p" in found[0]
    assert "(H2 + p -> He3 + g)" in found[0]


def test_labelled_variant_is_silent_and_copies_its_base():
    """A ___<label> variant loads silently with its base reaction's data."""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        net = load_network(PRIMATConfig({"network": "small"}),
                           reaction_names=_BASE + [_VARIANT])
    assert _duplicate_warnings(record) == []
    i, j = net.names.index("d_p__He3_g"), net.names.index("d_p__He3_g___alt")
    assert net.network[i] == net.network[j]
    # _abg/_fwd_median rows skip n__p, hence the -1.
    np.testing.assert_array_equal(net._abg[i - 1], net._abg[j - 1])
    np.testing.assert_array_equal(net._fwd_median[i - 1], net._fwd_median[j - 1])


@pytest.mark.slow
@pytest.mark.solve
@pytest.mark.skipif(not HAS_C_BACKEND, reason="C extension not built")
def test_backends_agree_on_duplicates(tmp_path, capfd):
    """Both backends print the same warning and solve a variant alike."""
    T9 = np.logspace(-3, 1, 60)
    text = "\n".join(f"{t:.6e} {1e-30 * t:.6e}" for t in T9)
    custom = {"removed": [], "added": {"He3_g__d_p": text}}
    messages = {}
    for b in ("python", "c"):
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            run_bbn({"network": "small"}, force_backend=b, custom_network=custom)
        err = capfd.readouterr().err
        messages[b] = (_duplicate_warnings(record) if b == "python" else
                       [l[len("warning: "):] for l in err.splitlines()
                        if "same channel" in l])
    assert len(messages["python"]) == 1
    assert messages["python"] == messages["c"]

    # The shipped small network plus a second d+p->He3+g channel.
    small = PRIMATConfig({"network": "small"}).resolve_rates_path(
        "nuclear", "networks", "small.txt")
    with open(small) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    (tmp_path / "networks").mkdir()
    (tmp_path / "networks" / "variant.txt").write_text(
        "\n".join(lines + [_VARIANT]) + "\n")
    params = {"network": "variant", "user_nuclear_dir": str(tmp_path)}
    doh = {b: run_bbn(dict(params), force_backend=b)["DoH"]
           for b in ("python", "c")}
    assert doh["python"] == pytest.approx(doh["c"], rel=5e-5)
