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


# The "Standard ratios" table of the Final abundances tab, as (key, format)
# pairs. Order and precision mirror primat.gui.panels._RATIO_FORMAT (the
# source of truth); duplicated here so this module has no import cycle with
# panels (panels imports this module). Precision follows CLAUDE.md.
_STANDARD_RATIOS = [
    ("Neff", ".8f"), ("YPBBN", ".8f"), ("YPCMB", ".8f"),
    ("DoH", ".7e"), ("He3oH", ".7e"), ("He3oHe4", ".6e"),
    ("Li7oH", ".6e"), ("Li6oLi7", ".6e"), ("YCNO", ".6e"),
]

_INI_HEADER = """# run_basic_from_gui.ini -- exported from the primat-gui "Final abundances" tab.
#
# Reproduces the tab's central run as a primat-c KEY=VALUE ini file. Every key
# not listed below keeps its DEFAULT_PARAMS default. Run from THIS file's own
# directory so a relative user_nuclear_dir resolves:
#   ./build/primat-c --ini run_basic_from_gui.ini
#
# For the +/- 1 sigma band, add:  --mc <N> --mc-seed 0

"""


def python_export_text(params, *, backend_used="auto", mc_quantities=None,
                       mc_samples=0, custom_network_name=None):
    """Return a standalone Python script reproducing the Final-abundances tab.

    Emits ``cfg = dict(...)`` from the GUI's "changed from default" params,
    calls ``run_bbn`` with the *pinned* backend that actually produced the
    displayed numbers, and prints every standard ratio the run produced from
    the deterministic ``run_bbn`` result (never an MC mean). When the session
    had "Quick MC uncertainty" on, a ``run_mc(seed=0, ...)`` block is appended
    that prints ONLY ``.std`` per quantity, matching ``app._quick_mc`` exactly
    so the +/- 1 sigma column reproduces bit-for-bit.

    Args:
        params: dict. The "changed from default" params (without a
            ``"custom_network"`` entry -- the caller strips it).
        backend_used: ``"c"``, ``"python"``, or ``"auto"``. Forwarded verbatim
            to ``run_bbn``/``run_mc`` as ``force_backend`` so the reproduction
            pins whichever backend the GUI used.
        mc_quantities: list[str] or None. Quantities for the appended
            ``run_mc`` std block; ``None``/empty appends no MC block.
        mc_samples: int. Sample count for the ``run_mc`` block.
        custom_network_name: str or None. When set, the script reproduces a
            customised network via the bundled ``nuclear/`` overlay: the base
            ``network`` is dropped and replaced by this name, and
            ``user_nuclear_dir="nuclear"`` is added.

    Returns:
        str. A runnable ``primat_gui_run.py``.

    Example:
        >>> python_export_text({"network": "small"}, backend_used="c")
    """
    cfg = {k: v for k, v in params.items() if k != "custom_network"}
    if custom_network_name is not None:
        # Reproduce the customised network from the self-contained overlay
        # packaged next to this script (nuclear/networks + nuclear/tables),
        # not the base network plus an inline custom_network dict.
        cfg.pop("network", None)
        cfg["network"] = custom_network_name
        cfg["user_nuclear_dir"] = "nuclear"

    body = ("\n".join(f"    {k}={_py_literal(cfg[k])}," for k in sorted(cfg))
            if cfg else
            "    # (every parameter left at its DEFAULT_PARAMS default)")
    ratios_src = "\n".join(f"    ({k!r}, {fmt!r})," for k, fmt in _STANDARD_RATIOS)

    lines = [
        "# -*- coding: utf-8 -*-",
        '"""',
        'Standalone PRIMAT run, exported from the primat-gui "Final abundances" tab.',
        "",
        "Reproduces that tab's standard ratios: central values come from run_bbn",
        "(not an MC mean); the +/- 1 sigma column, when present, comes only from",
        f"run_mc's .std. The backend is pinned to {backend_used!r} so the numbers",
        "reproduce bit-for-bit (C and Python differ at ~1e-5..1e-6).",
        "",
        "Run from THIS file's own directory:  python primat_gui_run.py",
        '"""',
        "from primat.backend import run_bbn",
        "",
        "cfg = dict(",
        body,
        ")",
        f"result = run_bbn(cfg, force_backend={backend_used!r})",
        "",
        "# Standard ratios, in the tab's order and precision.",
        "_RATIOS = [",
        ratios_src,
        "]",
        'print("Standard ratios (central values from run_bbn):")',
        "for _key, _fmt in _RATIOS:",
        "    if _key in result:",
        '        print(f"  {_key:8s} = {format(result[_key], _fmt)}")',
    ]
    if mc_quantities:
        q = list(mc_quantities)
        lines += [
            "",
            "# +/- 1 sigma from Monte-Carlo. Only .std is used; the central value",
            "# stays run_bbn's above. seed=0 + same backend => matches the tab.",
            "from primat.backend import run_mc",
            f"_mc = run_mc({mc_samples}, {q!r}, params=cfg, seed=0, "
            f"force_backend={backend_used!r})",
            f'print("\\n+/- 1 sigma (quick MC, {mc_samples} samples):")',
            "for _key, _fmt in _RATIOS:",
            f"    if _key in {q!r}:",
            '        print(f"  {_key:8s} +/- {format(_mc[_key].std, _fmt)}")',
        ]
    return "\n".join(lines) + "\n"


def ini_export_text(params, *, custom_network_name=None, mc_samples=0):
    """Return a standalone primat-c ``.ini`` reproducing the tab's central run.

    Same "changed from default" params as :func:`python_export_text`, emitted
    as ``KEY = VALUE`` for ``./build/primat-c --ini``. A customised network is
    reproduced via the bundled overlay: the base ``network`` is dropped and
    replaced by ``custom_network_name`` with ``user_nuclear_dir = "nuclear"``.
    The C CLI reproduces the central values; its own ``--mc N --mc-seed 0``
    flags reproduce the +/- 1 sigma band (see the bundle's README.txt).

    Args:
        params: dict. Changed-from-default params (no ``"custom_network"``).
        custom_network_name: str or None. See :func:`python_export_text`.
        mc_samples: int. Only used by the README's ``--mc`` hint, not emitted
            into the ini itself.

    Returns:
        str. A ``run_basic_from_gui.ini``.

    Example:
        >>> ini_export_text({"network": "large"})
    """
    cfg = {k: v for k, v in params.items() if k != "custom_network"}
    if custom_network_name is not None:
        cfg.pop("network", None)
        cfg["network"] = custom_network_name
        cfg["user_nuclear_dir"] = "nuclear"
    text = _INI_HEADER
    if cfg:
        text += "\n".join(f"{k} = {_ini_literal(cfg[k])}" for k in sorted(cfg))
        text += "\n"
    else:
        text += "# (every parameter left at its DEFAULT_PARAMS default)\n"
    return text
