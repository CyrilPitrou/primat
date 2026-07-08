# -*- coding: utf-8 -*-
"""
primat.gui.export_params
=========================

"Download params as .py / .ini" support for the Streamlit GUI (FABLEADVICE
S-15): serialise the params dict a GUI session has actually run into a
standalone Python script (mirroring ``runfiles/primat_run_explanatory.py``'s
layout) and into a ``primat-c`` ``.ini`` file (mirroring
``primat-c/examples/run_basic.ini``), so a configuration explored in the
browser can be reproduced from a script or the C CLI without the browser.

Only the keys the user actually changed are emitted (``params`` is already
that "changed subset", see ``params_form.render_sidebar_form``'s return-value
docstring) -- every other ``DEFAULT_PARAMS`` key is left at its default by
simply not appearing in the emitted ``cfg`` dict / ini file, exactly like
``primat.cli``'s "forward only what changed" behaviour. This keeps the
export self-contained: it does not need ``runfiles/primat_run_explanatory.py``
or ``primat-c/examples/run_basic.ini`` to exist on disk (neither ships in the
pip wheel, see ``pyproject.toml``'s ``packages`` list), so it works the same
whether the GUI is launched from a source checkout or ``pip install
".[gui]"``.
"""


def _py_literal(value):
    """Render ``value`` as a Python literal for a ``cfg = dict(...)`` line.

    ``repr`` already produces valid Python source for every value type a
    ``DEFAULT_PARAMS`` entry can hold (``bool``, ``int``, ``float``, ``str``,
    ``None``), so no per-type special-casing is needed here.
    """
    return repr(value)


def _ini_literal(value):
    """Render ``value`` as an ``.ini`` value, matching ``cpr_parse_literal``
    (``primat-c/src/config.c``): booleans as lower-case ``true``/``false``
    (``repr`` would give Python's capitalised ``True``/``False``, which
    ``cpr_parse_literal`` also accepts case-insensitively, but lower-case
    matches ``run_basic.ini``'s own style); strings double-quoted so a value
    containing spaces or starting with a digit round-trips unambiguously
    (``cpr_parse_literal`` strips matching quotes before falling through to
    its bool/int/float/string checks); everything else via ``str()``.
    """
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "none"
    if isinstance(value, str):
        return '"' + value.replace('"', '\\"') + '"'
    return str(value)


_PY_HEADER = '''# -*- coding: utf-8 -*-
"""
Standalone PRIMAT run, exported from the primat-gui sidebar.

Reproduces the configuration explored in the browser from a plain Python
script -- every parameter not listed below keeps its primat.config.
DEFAULT_PARAMS default (see that module, or
runfiles/primat_run_explanatory.py, for the full list of overridable keys
and their defaults).

Run with:

    python primat_gui_run.py
"""
from primat.backend import run_bbn  # , run_mc, dump_mc_samples (see bottom of file)

cfg = dict(
{body}
)
# force_backend: None/"auto" (default: C extension if built, else pure Python),
# "c", or "python" -- see primat/backend.py's module docstring for exactly
# which features always fall back to "python" regardless of this setting.
result = run_bbn(cfg, force_backend="auto")
print("Neff  =", result.get("Neff"))
print("YPBBN =", result["YPBBN"])
print("D/H   =", result["DoH"])
'''

_PY_MC_BLOCK = '''
# Monte-Carlo nuclear-rate/tau_n uncertainty propagation, mirroring the
# "Quick MC uncertainty" toggle active in the exported GUI session
# ({n_mc} samples):
from primat.backend import run_mc, dump_mc_samples, dump_mc_covariance, dump_mc_correlation

mc = run_mc({n_mc}, {quantities!r}, params=cfg, force_backend="auto")
for q in {quantities!r}:
    print(q, "=", mc[q].mean, "+/-", mc[q].std)
'''

_INI_HEADER = """# run_basic_from_gui.ini -- exported from the primat-gui sidebar.
#
# Reproduces the configuration explored in the browser as a primat-c
# KEY=VALUE ini file (see primat-c/examples/run_basic.ini for the full,
# heavily-commented list of overridable keys and their defaults). Every key
# not listed below keeps its DEFAULT_PARAMS default.
#
# Run with:
#   cd primat-c && make && ./build/primat-c --rates-dir <repo> --ini run_basic_from_gui.ini

"""

_CUSTOM_NETWORK_NOTE_PY = (
    "\n# NOTE: this GUI session used a customised reaction network (removed/"
    "\n# replaced/added reactions or rate-table overrides). That customisation"
    "\n# is not representable in this params dict -- export it separately from"
    "\n# the \"Reactions summary\" tab's \"Custom network\" download instead, and"
    "\n# pass it as run_bbn(cfg, custom_network=<the exported dict>).\n"
)

_CUSTOM_NETWORK_NOTE_INI = (
    "\n# NOTE: this GUI session used a customised reaction network (removed/"
    "\n# replaced/added reactions or rate-table overrides), which the .ini"
    "\n# format has no way to express -- export it separately from the"
    "\n# \"Reactions summary\" tab's \"Custom network\" download instead.\n"
)


def python_export_text(params, custom_network_active=False,
                        quick_mc=False, mc_samples=0, mc_quantities=None):
    """Return a standalone Python script reproducing ``params``.

    Parameters
    ----------
    params : dict
        The GUI's "changed from default" params dict (as returned by
        ``params_form.render_sidebar_form``, minus any ``"custom_network"``
        entry -- see ``custom_network_active``).
    custom_network_active : bool, optional
        Whether a customised network is active for this run; if so, a note
        is added pointing at the Reactions tab's own export instead of
        silently dropping the customisation.
    quick_mc : bool, optional
        Whether to append a ``run_mc`` block for the quantities the GUI's
        "Quick MC uncertainty" toggle covers.
    mc_samples : int, optional
        Sample count for the appended ``run_mc`` block.
    mc_quantities : list of str or None, optional
        Quantities passed to ``run_mc`` (only used when ``quick_mc``).

    Returns
    -------
    str
    """
    keys = [k for k in params if k != "custom_network"]
    if keys:
        body = "\n".join(
            f"    {k}={_py_literal(params[k])},"
            for k in sorted(keys)
        )
    else:
        body = "    # (every parameter left at its DEFAULT_PARAMS default)"
    text = _PY_HEADER.format(body=body)
    if custom_network_active:
        text += _CUSTOM_NETWORK_NOTE_PY
    if quick_mc and mc_quantities:
        text += _PY_MC_BLOCK.format(n_mc=mc_samples, quantities=list(mc_quantities))
    return text


def ini_export_text(params, custom_network_active=False):
    """Return a standalone ``primat-c`` ``.ini`` file reproducing ``params``.

    Parameters
    ----------
    params : dict
        Same "changed from default" dict as :func:`python_export_text`
        (``"custom_network"``, if present, is not representable in ``.ini``
        and is skipped -- see ``custom_network_active``).
    custom_network_active : bool, optional
        Whether a customised network is active for this run.

    Returns
    -------
    str
    """
    keys = [k for k in params if k != "custom_network"]
    text = _INI_HEADER
    if keys:
        text += "\n".join(
            f"{k} = {_ini_literal(params[k])}" for k in sorted(keys)
        )
        text += "\n"
    else:
        text += "# (every parameter left at its DEFAULT_PARAMS default)\n"
    if custom_network_active:
        text += _CUSTOM_NETWORK_NOTE_INI
    return text
