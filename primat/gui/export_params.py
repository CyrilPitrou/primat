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
                       mc_samples=0, custom_network=None):
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
        custom_network: dict or None. The exact ``{"removed", "replaced",
            "added", ...}`` override the GUI passed to ``run_bbn``. When given,
            it is embedded verbatim and passed as ``run_bbn``'s / ``run_mc``'s
            ``custom_network`` argument, with the *base* ``network`` kept as-is
            -- i.e. the identical call the GUI made, so the reproduction is
            bit-for-bit. (The ``.ini`` reproduces the same run via the bundled
            overlay and is also bit-for-bit, since ``ORDER_MT`` is aligned with
            ``ORDER_SMALL`` -- see ``build_reproduction_zip``.)

    Returns:
        str. A runnable ``primat_gui_run.py``.

    Example:
        >>> python_export_text({"network": "small"}, backend_used="c")
    """
    cfg = {k: v for k, v in params.items() if k != "custom_network"}

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
    ]
    if custom_network is not None:
        # Reproduce a customised network via the EXACT override the GUI passed
        # to run_bbn (base network kept as-is), embedded verbatim -- this is
        # bit-for-bit even for a small-based network, whose MT-era reaction
        # ordering depends on the base name and so cannot be reproduced by the
        # renamed nuclear/ overlay alone (that overlay is for the .ini path).
        lines += [
            "# This GUI session used a customised reaction network. It is",
            "# reproduced via the exact custom_network override the GUI passed",
            "# to run_bbn, so the numbers match the tab bit-for-bit.",
            f"custom_network = {custom_network!r}",
            "",
        ]
    run_bbn_call = ("run_bbn(cfg, custom_network=custom_network, "
                    f"force_backend={backend_used!r})"
                    if custom_network is not None
                    else f"run_bbn(cfg, force_backend={backend_used!r})")
    lines += [
        "cfg = dict(",
        body,
        ")",
        f"result = {run_bbn_call}",
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
        run_mc_extra = ("custom_network=custom_network, "
                        if custom_network is not None else "")
        lines += [
            "",
            "# +/- 1 sigma from Monte-Carlo. Only .std is used; the central value",
            "# stays run_bbn's above. seed=0 + same backend => matches the tab.",
            "from primat.backend import run_mc",
            f"_mc = run_mc({mc_samples}, {q!r}, params=cfg, seed=0, "
            f"{run_mc_extra}force_backend={backend_used!r})",
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


def _readme_text(*, active, backend_used, mc_samples, network_name):
    """Return the bundle's README.txt: how to run each artifact + caveats.

    Explains that the .py reproduces the tab bit-for-bit (pinned backend; a
    custom network via the exact custom_network override), that the .ini uses
    the bundled nuclear/ overlay (with the small-based-network caveat below),
    and -- when MC was on -- how to get the +/- 1 sigma band from each
    artifact, including the RNG caveat that the C CLI's band is bit-identical
    only when the GUI itself ran on the C backend.
    """
    lines = [
        "primat reproduction bundle",
        "==========================",
        "",
        "Run everything from THIS directory (a relative user_nuclear_dir/overlay",
        "resolves against the current directory).",
        "",
        "1) Python (reproduces the Final-abundances tab exactly):",
        "     python primat_gui_run.py",
        "",
        "2) primat-c CLI (central values):",
        "     ./build/primat-c --ini run_basic_from_gui.ini",
    ]
    if active:
        lines += [
            "",
            f"This run used a CUSTOM network ({network_name!r}), reproduced",
            "bit-for-bit by both artifacts:",
            "  - primat_gui_run.py embeds the exact custom_network override the",
            "    GUI passed to run_bbn.",
            "  - run_basic_from_gui.ini loads the network from the bundled",
            "    nuclear/ overlay (networks/ + tables/) via user_nuclear_dir.",
        ]
    if mc_samples:
        lines += [
            "",
            f"+/- 1 sigma (Monte-Carlo, {mc_samples} samples, seed 0):",
            "  - primat_gui_run.py already prints it (via run_mc).",
            f"  - primat-c: add  --mc {mc_samples} --mc-seed 0  to the CLI call.",
            "",
            "Bit-exactness caveat: C and Python have independent RNG streams, so",
            "the band is bit-identical to the tab only on the backend it ran on",
            f"({backend_used!r}). On the other backend it is statistically",
            "equivalent but not bit-identical.",
        ]
    return "\n".join(lines) + "\n"


def build_reproduction_zip(params, *, backend_used="auto", mc=None, cfg=None,
                           custom_network=None, kept_names=None,
                           network_name=None):
    """Pack a self-contained reproduction bundle (.zip) for the tab's results.

    Always writes ``primat_gui_run.py`` (:func:`python_export_text`),
    ``run_basic_from_gui.ini`` (:func:`ini_export_text`) and ``README.txt``.
    When ``custom_network`` is given, the fully-resolved network + its rate
    tables (including uploads) are bundled under ``nuclear/`` by reusing
    :func:`primat.gui.custom_rates.export_zip`. The ``.py`` reproduces it via
    the exact ``custom_network`` dict the GUI passed to ``run_bbn`` (base
    network kept -- bit-for-bit); the ``.ini`` (C CLI, no custom_network dict)
    reproduces it via ``network=<network_name>`` + ``user_nuclear_dir=nuclear``
    reading that bundled ``nuclear/`` overlay.

    Args:
        params: dict. The "changed from default" params (a ``"custom_network"``
            entry, if present, is ignored -- pass the decoded dict as
            ``custom_network`` instead).
        backend_used: ``"c"``/``"python"``/``"auto"`` -- the pinned backend.
        mc: primat.main.MCResult or None. When given, its standard-ratio
            quantities and sample count drive the .py's run_mc std block.
        cfg: primat.config.PRIMATConfig. Required only when ``custom_network``
            is set (forwarded to ``export_zip`` for the overlay).
        custom_network: dict or None. ``{"removed", "replaced", "added"}``.
        kept_names: sequence[str] or None. Ordered kept reaction names for the
            overlay network file (required when ``custom_network`` is set).
        network_name: str or None. Basename for the overlay network file and
            the reproduced ``network=`` value (required when custom).

    Returns:
        bytes. The zip contents.

    Example:
        >>> build_reproduction_zip({"network": "small"}, backend_used="c")
    """
    import io
    import zipfile

    from primat.gui.custom_rates import export_zip

    active = custom_network is not None
    cn_name = network_name if active else None

    # MC quantities/sample count come straight from the cached MCResult, so the
    # reproduction's run_mc call matches app._quick_mc's quantities exactly.
    mc_quantities = None
    mc_samples = 0
    if mc is not None:
        mc_names = set(mc.quantity_names())
        mc_quantities = [k for k, _ in _STANDARD_RATIOS if k in mc_names]
        mc_samples = len(next(iter(mc._data.values())).values)

    # The .py reproduces a custom network via the EXACT custom_network dict the
    # GUI passed to run_bbn (base network kept). The .ini, which the C CLI reads
    # and which cannot carry a custom_network dict, uses the bundled nuclear/
    # overlay instead (network=<name> + user_nuclear_dir). Both are bit-for-bit:
    # ORDER_MT is aligned so ORDER_SMALL is a prefix subsequence of it (see
    # network_data.py), so a renamed overlay of a small-based network selects
    # the MT era in the same order as the base "small" run -- no ~1e-6 drift.
    py_text = python_export_text(
        params, backend_used=backend_used, mc_quantities=mc_quantities,
        mc_samples=mc_samples, custom_network=custom_network)
    ini_text = ini_export_text(
        params, custom_network_name=cn_name, mc_samples=mc_samples)
    readme = _readme_text(active=active, backend_used=backend_used,
                          mc_samples=mc_samples, network_name=cn_name)

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("primat_gui_run.py", py_text)
        zf.writestr("run_basic_from_gui.ini", ini_text)
        zf.writestr("README.txt", readme)
        if active:
            # Reuse the Reactions-tab exporter, which already writes a
            # self-contained networks/ + tables/ overlay (inlining uploaded and
            # edited tables). Re-nest it under nuclear/ so user_nuclear_dir
            # points at a single overlay root.
            inner = export_zip(cfg, custom_network, kept_names,
                               network_filename=network_name)
            with zipfile.ZipFile(io.BytesIO(inner)) as innerzf:
                for item in innerzf.namelist():
                    zf.writestr(f"nuclear/{item}", innerzf.read(item))
    return buf.getvalue()
