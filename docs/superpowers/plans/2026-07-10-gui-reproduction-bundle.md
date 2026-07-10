# GUI Reproduction Bundle Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the Output-tables `.py`/`.ini` reproduction downloads with a single self-contained "Download reproduction bundle (.zip)" button in the Final abundances tab that reproduces the tab's numbers exactly — including custom/uploaded networks.

**Architecture:** `primat/gui/export_params.py` gains a `build_reproduction_zip()` that emits `primat_gui_run.py` (run_bbn centrals + run_mc std-only, pinned backend), `run_basic_from_gui.ini`, `README.txt`, and — when a custom network is active — a bundled `nuclear/` overlay produced by reusing `custom_rates.export_zip()`. The reproduction reconstructs a custom network via `network=<name>` + `user_nuclear_dir="nuclear"`, an overlay both backends already consume. The button moves from `render_downloads_panel` to `render_results_panel`.

**Tech Stack:** Python 3.12, Streamlit (GUI, tested via `streamlit.testing.v1.AppTest`), pytest, `zipfile`/`io` (stdlib).

## Global Constraints

- Report BBN observables to CLAUDE.md precision: `Neff`/`YPBBN`/`YPCMB` `.8f`, `DoH`/`He3oH` `.7e`, `He3oHe4`/`Li7oH`/`Li6oLi7`/`YCNO` `.6e`.
- Heavy commenting: every new function needs a docstring (what/why, args w/ units, a usage example) and inline comments on non-obvious steps (scientific-code convention).
- Do NOT touch `primat-c/` — this is GUI/export plumbing only; the C backend already consumes `user_nuclear_dir` and offers `--mc N --mc-seed SEED`.
- Do NOT change the raw output-file downloads in the Output tables tab (`output_final.txt`, `output_time_evolution.tsv`, `output_background.tsv`, `nTOp_total.tsv`, `decays.txt`, `output_mc_*.tsv`) — only the two reproduction buttons move out.
- Backend pin values passed to `run_bbn`/`run_mc` are the literal strings `"c"` or `"python"` (never the display strings `"C"`/`"Python"`).
- MC reproduction must use `seed=0` and the same pinned backend, matching `app._quick_mc`, so `.std` is bit-identical.
- Run tests from the repo root: `python -m pytest tests/test_gui.py tests/test_gui_custom_network.py -q`.

---

## File Structure

- `primat/gui/export_params.py` — **modify**: rewrite `python_export_text`/`ini_export_text`; add `_STANDARD_RATIOS`, `_readme_text`, `build_reproduction_zip`; drop the `_CUSTOM_NETWORK_NOTE_*` "not representable" paths.
- `primat/gui/panels.py` — **modify**: add the reproduction-bundle button to `render_results_panel`; remove the two reproduction buttons + `run_params` param from `render_downloads_panel`.
- `primat/gui/app.py` — **modify**: pass `run_params`/`backend_used` into `render_results_panel`; drop `run_params` from the `render_downloads_panel` call.
- `tests/test_gui.py` — **modify**: repoint/rename the two export tests; add exporter-content tests.
- `tests/test_gui_custom_network.py` — **modify**: add zip-structure + overlay-vs-dict parity tests.

---

### Task 1: Rewrite `export_params.py` text exporters

**Files:**
- Modify: `primat/gui/export_params.py`
- Test: `tests/test_gui.py`

**Interfaces:**
- Produces:
  - `_STANDARD_RATIOS: list[tuple[str, str]]` — ordered `(key, format-spec)` pairs mirroring `panels._RATIO_FORMAT`.
  - `python_export_text(params, *, backend_used="auto", mc_quantities=None, mc_samples=0, custom_network_name=None) -> str`
  - `ini_export_text(params, *, custom_network_name=None, mc_samples=0) -> str`
- Consumes: `_py_literal`, `_ini_literal` (unchanged helpers already in the file).

- [ ] **Step 1: Write failing tests for the new exporter behaviour**

Replace the body of `tests/test_gui.py::test_gui_exported_python_script_reproduces_run` and `::test_export_params_ini_text_round_trips_values`, and add three new tests. Put this block in `tests/test_gui.py` (replacing the two existing functions named below):

```python
def test_gui_exported_python_script_reproduces_run():
    """The exported .py execs standalone and rebuilds the same cfg, pins the
    backend, and prints run_bbn central values (not an MC mean)."""
    import unittest.mock
    from primat.gui.export_params import python_export_text
    from primat.gui.session_keys import SessionKeys

    at = AppTest.from_file(APP_PATH)
    at.run(timeout=60)
    [omegabh2] = [ni for ni in at.sidebar.number_input if ni.key == "Omegabh2"]
    omegabh2.set_value(0.0225)
    at.run(timeout=60)
    _run_bbn(at)
    assert not at.exception

    stored_params = at.session_state[SessionKeys.params]
    script = python_export_text(
        {k: v for k, v in stored_params.items() if k != "custom_network"},
        backend_used="python")
    assert "Omegabh2=0.0225" in script
    assert "force_backend='python'" in script
    assert "run_bbn" in script and ".mean" not in script

    captured = {}

    def fake_run_bbn(cfg, force_backend="auto"):
        captured["cfg"] = cfg
        # Include every standard-ratio key so the print loop exercises them.
        return {k: 0.0 for k, _ in
                __import__("primat.gui.export_params", fromlist=["_STANDARD_RATIOS"])._STANDARD_RATIOS}

    with unittest.mock.patch("primat.backend.run_bbn", fake_run_bbn):
        exec(compile(script, "primat_gui_run.py", "exec"), {})
    assert captured["cfg"]["Omegabh2"] == 0.0225


def test_export_params_ini_text_round_trips_values():
    """ini_export_text keeps cpr_parse_literal value syntax and the empty note."""
    from primat.gui.export_params import ini_export_text

    text = ini_export_text({"verbose": True, "network": "large", "Omegabh2": 0.0225})
    assert "verbose = true" in text
    assert 'network = "large"' in text
    assert "Omegabh2 = 0.0225" in text
    assert "every parameter left at its DEFAULT_PARAMS default" in ini_export_text({})


def test_export_py_prints_full_standard_ratio_set():
    """The exported .py prints every standard ratio, in CLAUDE.md precision."""
    from primat.gui.export_params import python_export_text, _STANDARD_RATIOS

    script = python_export_text({"network": "small"}, backend_used="c")
    for key, _fmt in _STANDARD_RATIOS:
        assert repr(key) in script
    assert "force_backend='c'" in script


def test_export_py_mc_block_uses_std_only_and_seed_zero():
    """When MC quantities are given, the .py adds a run_mc(seed=0) std-only block."""
    from primat.gui.export_params import python_export_text

    script = python_export_text(
        {"network": "small"}, backend_used="python",
        mc_quantities=["YPBBN", "DoH"], mc_samples=30)
    assert "run_mc(30, ['YPBBN', 'DoH']" in script
    assert "seed=0" in script
    assert "force_backend='python'" in script
    assert ".std" in script
    assert ".mean" not in script


def test_export_custom_network_sets_user_nuclear_dir():
    """A custom-network export drops the base network and points at the overlay."""
    from primat.gui.export_params import python_export_text, ini_export_text

    py = python_export_text({"network": "small"}, backend_used="c",
                            custom_network_name="mynet")
    assert "network='mynet'" in py
    assert "user_nuclear_dir='nuclear'" in py

    ini = ini_export_text({"network": "small"}, custom_network_name="mynet")
    assert 'network = "mynet"' in ini
    assert 'user_nuclear_dir = "nuclear"' in ini
```

- [ ] **Step 2: Run the new tests to verify they fail**

Run: `python -m pytest tests/test_gui.py::test_export_py_prints_full_standard_ratio_set tests/test_gui.py::test_export_py_mc_block_uses_std_only_and_seed_zero tests/test_gui.py::test_export_custom_network_sets_user_nuclear_dir -q`
Expected: FAIL (`TypeError: python_export_text() got an unexpected keyword argument 'backend_used'` / missing `_STANDARD_RATIOS`).

- [ ] **Step 3: Implement the new exporters**

In `primat/gui/export_params.py`, add near the top (after the `_ini_literal` helper):

```python
# The "Standard ratios" table of the Final abundances tab, as (key, format)
# pairs. Order and precision mirror primat.gui.panels._RATIO_FORMAT (the
# source of truth); duplicated here so this module has no import cycle with
# panels (panels imports this module). Precision follows CLAUDE.md.
_STANDARD_RATIOS = [
    ("Neff", ".8f"), ("YPBBN", ".8f"), ("YPCMB", ".8f"),
    ("DoH", ".7e"), ("He3oH", ".7e"), ("He3oHe4", ".6e"),
    ("Li7oH", ".6e"), ("Li6oLi7", ".6e"), ("YCNO", ".6e"),
]
```

Replace the whole `python_export_text` function with:

```python
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
```

Replace the whole `ini_export_text` function with:

```python
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
```

Delete the now-unused `_CUSTOM_NETWORK_NOTE_PY` and `_CUSTOM_NETWORK_NOTE_INI` constants. Update `_PY_HEADER`/`_PY_MC_BLOCK`/`_INI_HEADER`: `_PY_HEADER` and `_PY_MC_BLOCK` are no longer referenced (the new `python_export_text` builds text inline) — delete them. Keep `_INI_HEADER` but update its "Run with" line to note running from the extracted dir:

```python
_INI_HEADER = """# run_basic_from_gui.ini -- exported from the primat-gui "Final abundances" tab.
#
# Reproduces the tab's central run as a primat-c KEY=VALUE ini file. Every key
# not listed below keeps its DEFAULT_PARAMS default. Run from THIS file's own
# directory so a relative user_nuclear_dir resolves:
#   ./build/primat-c --ini run_basic_from_gui.ini
#
# For the +/- 1 sigma band, add:  --mc <N> --mc-seed 0

"""
```

- [ ] **Step 4: Run the exporter tests to verify they pass**

Run: `python -m pytest tests/test_gui.py::test_export_py_prints_full_standard_ratio_set tests/test_gui.py::test_export_py_mc_block_uses_std_only_and_seed_zero tests/test_gui.py::test_export_custom_network_sets_user_nuclear_dir tests/test_gui.py::test_gui_exported_python_script_reproduces_run tests/test_gui.py::test_export_params_ini_text_round_trips_values -q`
Expected: PASS (5 passed).

- [ ] **Step 5: Commit**

```bash
git add primat/gui/export_params.py tests/test_gui.py
git commit -m "GUI export: full-ratio run_bbn centrals + MC std-only, pinned backend, overlay custom-network config"
```

---

### Task 2: Add `build_reproduction_zip` + README

**Files:**
- Modify: `primat/gui/export_params.py`
- Test: `tests/test_gui_custom_network.py`

**Interfaces:**
- Consumes: `python_export_text`, `ini_export_text`, `_STANDARD_RATIOS` (Task 1); `primat.gui.custom_rates.export_zip(cfg, custom_network, kept_names, network_filename=...) -> bytes`.
- Produces:
  - `_readme_text(*, active, backend_used, mc_samples, network_name) -> str`
  - `build_reproduction_zip(params, *, backend_used="auto", mc=None, cfg=None, custom_network=None, kept_names=None, network_name=None) -> bytes`

- [ ] **Step 1: Write failing tests for the zip builder**

Add to `tests/test_gui_custom_network.py`:

```python
def test_reproduction_zip_standard_run_has_three_files_no_overlay():
    """A non-custom bundle is exactly py + ini + README, no nuclear/ tree."""
    import io, zipfile
    from primat.gui.export_params import build_reproduction_zip

    data = build_reproduction_zip({"network": "small"}, backend_used="python")
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        names = set(zf.namelist())
    assert names == {"primat_gui_run.py", "run_basic_from_gui.ini", "README.txt"}


def test_reproduction_zip_custom_run_bundles_nuclear_overlay():
    """A custom bundle adds nuclear/networks + nuclear/tables from export_zip."""
    import io, zipfile
    from primat.config import PRIMATConfig
    from primat.network_data import UpdateNuclearRates
    from primat.gui.export_params import build_reproduction_zip

    cfg = PRIMATConfig({"network": "small"})
    nucl = UpdateNuclearRates(cfg)
    rows = [r for r in nucl.describe_reactions() if r[0] != "n__p"]
    kept_names = [r[0] for r in rows]
    custom_network = {"removed": [], "replaced": {}, "added": {}}

    data = build_reproduction_zip(
        {"network": "small"}, backend_used="python", cfg=cfg,
        custom_network=custom_network, kept_names=kept_names,
        network_name="mynet")
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        names = zf.namelist()
    assert "primat_gui_run.py" in names
    assert "nuclear/networks/mynet.txt" in names
    assert any(n.startswith("nuclear/tables/") for n in names)
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_gui_custom_network.py::test_reproduction_zip_standard_run_has_three_files_no_overlay tests/test_gui_custom_network.py::test_reproduction_zip_custom_run_bundles_nuclear_overlay -q`
Expected: FAIL (`ImportError: cannot import name 'build_reproduction_zip'`).

- [ ] **Step 3: Implement the builder + README**

Append to `primat/gui/export_params.py`:

```python
def _readme_text(*, active, backend_used, mc_samples, network_name):
    """Return the bundle's README.txt: how to run each artifact + caveats.

    Explains that the .py reproduces the tab bit-for-bit (pinned backend),
    that a custom network is reproduced via the bundled nuclear/ overlay, and
    -- when MC was on -- how to get the +/- 1 sigma band from each artifact,
    including the RNG caveat that the C CLI's band is bit-identical only when
    the GUI itself ran on the C backend.
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
            f"This run used a CUSTOM network ({network_name!r}). It is bundled",
            "under nuclear/ (networks/ + tables/) and loaded via",
            "user_nuclear_dir=nuclear -- no internet or extra files needed.",
        ]
    if mc_samples:
        lines += [
            "",
            f"+/- 1 sigma (Monte-Carlo, {mc_samples} samples, seed 0):",
            "  - primat_gui_run.py already prints it (via run_mc).",
            f"  - primat-c: add  --mc {mc_samples} --mc-seed 0  to the CLI call.",
            "",
            "Bit-exactness caveat: C and Python have independent RNG streams, so",
            f"the band is bit-identical to the tab only on the backend it ran on",
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
    :func:`primat.gui.custom_rates.export_zip`, and the emitted .py/.ini
    reproduce it via ``network=<network_name>`` + ``user_nuclear_dir=nuclear``.

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

    py_text = python_export_text(
        params, backend_used=backend_used, mc_quantities=mc_quantities,
        mc_samples=mc_samples, custom_network_name=cn_name)
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
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tests/test_gui_custom_network.py::test_reproduction_zip_standard_run_has_three_files_no_overlay tests/test_gui_custom_network.py::test_reproduction_zip_custom_run_bundles_nuclear_overlay -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add primat/gui/export_params.py tests/test_gui_custom_network.py
git commit -m "GUI export: build_reproduction_zip + README bundling nuclear/ overlay for custom networks"
```

---

### Task 3: Wire the button into the Final abundances tab, remove from Output tables

**Files:**
- Modify: `primat/gui/panels.py:70` (`render_results_panel`), `primat/gui/panels.py:358` (`render_downloads_panel`)
- Modify: `primat/gui/app.py:528-553` (`main`)
- Test: `tests/test_gui.py`

**Interfaces:**
- Consumes: `build_reproduction_zip` (Task 2); `run.nucl.describe_reactions()`; `SessionKeys.active_custom_network`.
- Produces: `render_results_panel(run, mc=None, run_params=None, backend_used=None)`; `render_downloads_panel(run, mc=None, background=None)` (no `run_params`).

- [ ] **Step 1: Write failing test for the button's new location**

Replace `tests/test_gui.py::test_downloads_panel_offers_gui_export` with:

```python
def test_results_panel_offers_reproduction_bundle():
    """The reproduction bundle download lives in the Final abundances tab,
    and is gone from the Output tables tab."""
    at = AppTest.from_file(APP_PATH)
    at.run(timeout=60)
    [omegabh2] = [ni for ni in at.sidebar.number_input if ni.key == "Omegabh2"]
    omegabh2.set_value(0.0225)
    at.run(timeout=60)
    _run_bbn(at)
    assert not at.exception

    # NOTE: _download_button matches on the button LABEL (returns None if
    # absent), so assert against the label text, not the file_name.
    assert _download_button(at, "Download reproduction bundle (.zip)") is not None
    # The old per-format reproduction buttons (whose labels WERE the filenames)
    # are gone.
    assert _download_button(at, "primat_gui_run.py") is None
    assert _download_button(at, "run_basic_from_gui.ini") is None
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_gui.py::test_results_panel_offers_reproduction_bundle -q`
Expected: FAIL (no `reproduction_bundle.zip` button yet; `primat_gui_run.py` still present).

- [ ] **Step 3: Add the button to `render_results_panel`**

In `primat/gui/panels.py`, change the signature and add the download block. Update the signature line:

```python
def render_results_panel(run, mc=None, run_params=None, backend_used=None):
```

At the end of `render_results_panel` (after the "Final abundances" per-nuclide `st.markdown("\n".join(lines))`), append:

```python
    # Reproduction bundle: the files that reproduce EXACTLY this tab's numbers.
    # Only shown once we know both the params that ran and which backend ran
    # (needed to pin force_backend for bit-for-bit reproduction).
    if run_params is not None and backend_used is not None:
        import json

        from primat.gui.export_params import build_reproduction_zip

        cn_json = run_params.get("custom_network")
        custom_network = json.loads(cn_json) if cn_json else None
        network_name = None
        kept_names = None
        if custom_network is not None:
            # Same title/kept-name derivation as _render_reaction_downloads,
            # so the overlay network file matches the Reactions-tab export.
            active = st.session_state.get(SessionKeys.active_custom_network)
            network_name = (active["title"] if active else
                            (run.cfg.network if run.cfg.network != "large" else "large"))
            kept_names = [name for name, equation, source, file
                          in run.nucl.describe_reactions() if name != "n__p"]
        params_only = {k: v for k, v in run_params.items() if k != "custom_network"}
        # _describe_backend_used returns display strings ("C"/"Python"/...);
        # force_backend wants "c"/"python".
        backend_pin = "c" if backend_used == "C" else "python"
        try:
            zip_bytes = build_reproduction_zip(
                params_only, backend_used=backend_pin, mc=mc, cfg=run.cfg,
                custom_network=custom_network, kept_names=kept_names,
                network_name=network_name)
        except Exception as exc:  # never let an export failure break the tab
            st.warning(f"Could not build the reproduction bundle: {exc}")
        else:
            st.markdown("**Reproduce these results**")
            st.download_button(
                "Download reproduction bundle (.zip)", data=zip_bytes,
                file_name="reproduction_bundle.zip", mime="application/zip",
                key="dl_reproduction",
                help="A self-contained script + primat-c .ini (+ the custom "
                     "network, if any) that reproduce exactly the values above.",
            )
```

Ensure `SessionKeys` is imported at the top of `panels.py` (it already imports from `primat.gui`; if `SessionKeys` is not imported, add `from primat.gui.session_keys import SessionKeys`).

Run: `grep -n "SessionKeys" primat/gui/panels.py | head -2`
If no import line exists, add `from primat.gui.session_keys import SessionKeys` alongside the other `primat.gui` imports.

- [ ] **Step 4: Remove the reproduction buttons from `render_downloads_panel`**

In `primat/gui/panels.py`, change the signature:

```python
def render_downloads_panel(run, mc=None, background=None):
```

Delete the entire trailing `if run_params is not None:` block (the two `_file_download` calls for `primat_gui_run.py` and `run_basic_from_gui.ini`, plus the `custom_network_active`/`params_only`/`quick_mc`/`mc_quantities`/`mc_samples` locals feeding them). Remove the `run_params` paragraph from the docstring.

- [ ] **Step 5: Update `app.main` wiring**

In `primat/gui/app.py`, in `main`, change the two render calls (currently around lines 543 and 552):

```python
        panels.render_results_panel(run, mc=mc, run_params=stored_params,
                                    backend_used=backend_used)
```
and
```python
        panels.render_downloads_panel(run, mc=mc, background=background)
```

- [ ] **Step 6: Run the wiring test + the earlier tests**

Run: `python -m pytest tests/test_gui.py::test_results_panel_offers_reproduction_bundle tests/test_gui.py -q`
Expected: PASS (the whole `test_gui.py` file, including the Task 1 exporter tests).

- [ ] **Step 7: Commit**

```bash
git add primat/gui/panels.py primat/gui/app.py tests/test_gui.py
git commit -m "GUI: move reproduction download to Final abundances tab as a single bundle"
```

---

### Task 4: Custom-network overlay-vs-dict parity (the must-verify)

**Files:**
- Test: `tests/test_gui_custom_network.py`

**Interfaces:**
- Consumes: `build_reproduction_zip` (Task 2); `primat.backend.run_bbn(params, custom_network=..., force_backend=...)`.

- [ ] **Step 1: Write the parity test**

Add to `tests/test_gui_custom_network.py`:

```python
def test_reproduction_overlay_matches_custom_network_dict_run(tmp_path):
    """The bundled nuclear/ overlay reproduces the custom_network-dict run
    bit-for-bit -- the core "reproduce exactly" guarantee for custom networks.

    Customisation: replace one small-network reaction's rate table with its
    own shipped text (a no-op edit that still forces the full overlay path --
    a genuinely different table would only widen the gap this test guards
    against). Both runs use the Python backend so they are deterministic and
    need no C build.
    """
    import io
    import zipfile

    from primat.backend import run_bbn
    from primat.config import PRIMATConfig
    from primat.network_data import UpdateNuclearRates
    from primat.gui.export_params import build_reproduction_zip

    cfg = PRIMATConfig({"network": "small"})
    nucl = UpdateNuclearRates(cfg)
    rows = [r for r in nucl.describe_reactions() if r[0] != "n__p"]
    kept_names = [r[0] for r in rows]

    # Pick the first non-weak reaction with a real per-reaction table file and
    # read its shipped text, to feed a no-op "replaced" entry.
    name, _eq, _src, fname = next(r for r in rows if r[3] not in (None, "", "--"))
    shipped_path = cfg.resolve_rates_path("nuclear", "tables", name, fname)
    with open(shipped_path) as fh:
        shipped_text = fh.read()
    custom_network = {"removed": [], "replaced": {name: shipped_text}, "added": {}}

    # A) the run the GUI actually makes: base network + custom_network dict.
    a = run_bbn({"network": "small"}, custom_network=custom_network,
                force_backend="python")

    # B) the reproduction: extract the bundle and run via the overlay.
    data = build_reproduction_zip(
        {"network": "small"}, backend_used="python", cfg=cfg,
        custom_network=custom_network, kept_names=kept_names,
        network_name="small")
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        zf.extractall(tmp_path)
    b = run_bbn(
        {"network": "small", "user_nuclear_dir": str(tmp_path / "nuclear")},
        force_backend="python")

    for key in ("YPBBN", "DoH", "He3oH", "Li7oH", "Neff"):
        assert abs(a[key] - b[key]) <= 1e-9 * abs(a[key]) + 1e-12, (
            f"{key}: overlay {b[key]!r} != dict-run {a[key]!r}")
```

- [ ] **Step 2: Run the parity test**

Run: `python -m pytest tests/test_gui_custom_network.py::test_reproduction_overlay_matches_custom_network_dict_run -q`
Expected: PASS. If it FAILS, the overlay and the custom_network-dict path genuinely diverge (e.g. an MT-era reaction whose table the overlay applies but the dict path does not, or vice-versa) — that is the bug the spec flagged; STOP and investigate before continuing (do not loosen the tolerance to pass).

- [ ] **Step 3: End-to-end drive via the actual GUI (verification, not just unit test)**

Confirm the button works in a live-ish AppTest custom-network session. Add:

```python
def test_gui_custom_network_offers_reproduction_bundle():
    """A GUI session with an active custom network offers the bundle in the
    Final abundances tab (companion to the direct build_reproduction_zip
    tests -- exercises the panel wiring end to end)."""
    at = _run_gui_with_custom_network()  # reuse this file's existing helper
    assert not at.exception
    assert _download_button(at, "Download reproduction bundle (.zip)") is not None
```

Check the helper name first: `grep -n "def _run_gui_with_custom_network\|custom.*network.*helper\|def _apply_custom" tests/test_gui_custom_network.py`. If no reusable helper exists, model the setup on this file's existing custom-network AppTest test (whichever test drives "Create custom network" / applies a `custom_network` param), then assert the button. If `_download_button`'s accessor differs, match the file's style.

- [ ] **Step 4: Run the full GUI suite**

Run: `python -m pytest tests/test_gui.py tests/test_gui_custom_network.py -q`
Expected: PASS (all).

- [ ] **Step 5: Commit**

```bash
git add tests/test_gui_custom_network.py
git commit -m "GUI export: parity + wiring tests for custom-network reproduction bundle"
```

---

## Self-Review Notes

- **Spec coverage:** always-zip (Task 2/3), Final-abundances placement + Output-tables removal (Task 3), full standard-ratio print (Task 1), run_bbn-central + MC-std-only (Task 1), pinned backend (Task 1/3), custom-network overlay via `export_zip` (Task 2), `.ini` `--mc` README hint + RNG caveat (Task 2), overlay-vs-dict parity (Task 4). All covered.
- **Backend string:** `_describe_backend_used` returns `"C"`/`"Python"`/`"Python (C extension unavailable)"`; Task 3 maps `"C" -> "c"`, everything else `-> "python"`. Consistent with Task 1's `backend_used`/`force_backend` literals.
- **No import cycle:** `export_params` does NOT import `panels`; `_STANDARD_RATIOS` is duplicated with a source-of-truth comment. `build_reproduction_zip` imports `custom_rates.export_zip` locally.
- **AppTest bytes caveat:** `download_button` protos don't expose bytes, so the content assertions run against `build_reproduction_zip`/`python_export_text` directly (Tasks 1/2/4); AppTest tests only assert button presence (Task 3/4 Step 3).
