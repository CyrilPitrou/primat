# GUI reproduction bundle (Final abundances tab)

**Date:** 2026-07-10
**Status:** Approved design, pre-implementation
**Area:** `primat/gui/` (Streamlit app)

## Problem

The GUI currently offers "GUI configuration as Python script" (`.py`) and
"GUI configuration as `.ini`" downloads in the **Output tables** tab
(`panels.render_downloads_panel`, via `export_params.py`). Three shortcomings:

1. **Custom networks don't work.** When a custom/uploaded network is active,
   the `.py`/`.ini` only emit a *note* ("customisation is not representable —
   export it separately from the Reactions tab") instead of a runnable
   reproduction. Uploaded rate tables and removed/added/replaced reactions are
   silently dropped from the reproduction.
2. **MC mean instead of the run_bbn central value.** The appended MC block
   prints `mc[q].mean ± mc[q].std`. The user wants the **central value from
   `run_bbn`** (matching the "Final abundances" tab), with the MC used **only**
   for the `± 1σ` (std).
3. **Too few observables.** The `.py` prints only `Neff`, `YPBBN`, `D/H`. It
   should reproduce the full **Standard ratios** table shown in the Final
   abundances tab.

Additionally, the reproduction download lives in the wrong tab: it belongs in
the **Final abundances** tab, so that the tab a user is looking at is the one
offering "download the files that reproduce exactly what you see here".

## Goals

- One **"Download reproduction bundle (.zip)"** button in the **Final
  abundances** tab, directly under the Standard-ratios table.
- The bundle reproduces the displayed numbers **exactly** (bit-for-bit on the
  same machine/backend), including custom networks with uploaded rate tables.
- The `.ini` becomes fully functional for custom networks (first time).

## Non-goals

- No change to the raw output-file downloads (`output_final.txt`,
  `output_time_evolution.tsv`, `output_background.tsv`, `nTOp_total.tsv`,
  `decays.txt`, `output_mc_*.tsv`) — they stay in the **Output tables** tab.
- No new physics/numerics; no `primat-c` numerical changes (this is GUI +
  export plumbing only). The C backend already consumes `user_nuclear_dir`.

## Decisions (resolved with the author)

- **Packaging:** *always a single `.zip`* (uniform UX), even for a standard run
  (then the zip is just `.py` + `.ini` + `README.txt`, no `nuclear/`).
- **Backend:** *pin the exact backend that ran* (`force_backend="c"` or
  `"python"`, from `_solve`'s `backend_used`) in both the `run_bbn` and
  `run_mc` calls, so central values and MC std reproduce exactly. Accepted
  tradeoff: a `force_backend="c"` script raises on a machine without the C
  extension built (rare — C is the default build) rather than silently using
  slightly different Python numbers.

## Bundle layout

```
reproduction_bundle.zip
  primat_gui_run.py          # run_bbn (+ run_mc for std only), pinned backend
  run_basic_from_gui.ini     # primat-c CLI equivalent (KEY=VALUE)
  README.txt                 # how to run each; CWD / user_nuclear_dir note
  nuclear/                   # PRESENT ONLY when a custom network is active
    networks/<name>.txt      # fully-resolved network (one reaction per line)
    tables/<name>/<file>     # every kept reaction's table, incl. uploads/edits
```

The `nuclear/` subtree is produced by the existing
`primat.gui.custom_rates.export_zip()` (it already writes exactly the
`networks/` + `tables/` layout that is the `user_nuclear_dir` overlay root, and
already inlines uploaded/edited tables, shipped-table matches, and inline decay
overrides). We reuse it verbatim (writing its entries under a `nuclear/`
prefix inside the combined zip).

## Reproduction semantics

### `primat_gui_run.py`

- Builds `cfg = dict(...)` from the GUI's "changed from default" params, with
  two custom-network adjustments (see below).
- `result = run_bbn(cfg, force_backend="<pinned>")`.
- Prints **every Standard ratio the run produced**, in the tab's order:
  `Neff, YPBBN, YPCMB, DoH, He3oH, He3oHe4, Li7oH, Li6oLi7, YCNO`
  (the `_RATIO_FORMAT` set in `panels.py`, filtered to keys present in
  `result`), each from `result[key]` — i.e. the deterministic central value,
  **not** an MC mean. Same decimal precision as the tab (`_RATIO_FORMAT`).
- If quick-MC was active in the session:
  `mc = run_mc(num_mc, quantities, params=cfg, seed=0, force_backend="<pinned>")`
  where `quantities = [q for q in _RATIO_FORMAT if q in result]` and
  `num_mc` is the session's sample count — matching `app._quick_mc` exactly
  (`seed=0`, same quantities, same backend) so `mc[q].std` is bit-identical to
  the tab. The script prints the `± 1σ` column using **only** `mc[q].std`
  (never `mc[q].mean`).

### `run_basic_from_gui.ini`

- Same changed params as KEY=VALUE. Consumed by
  `./build/primat-c --ini run_basic_from_gui.ini` (from the extracted dir).
- Reproduces the **central** values via the C CLI's normal output.
- The `± 1σ` band is reproducible natively in C via the CLI's own MC flags
  (`cli.c`: `--mc N --mc-seed SEED`, writing
  `<output_mc_file_prefix>_covariance.tsv`/`_samples.tsv`). When quick-MC was
  active the `README.txt` shows the concrete command:
  `./build/primat-c --ini run_basic_from_gui.ini --mc <num_mc> --mc-seed 0`.
- **Bit-exactness caveat (RNG streams differ per backend):** the C-CLI MC band
  is bit-identical to the tab **only when the GUI session itself ran on the C
  backend** (`backend_used == "c"`). If the session ran on Python, the C-CLI
  gives a statistically-equivalent (not bit-identical) band; the bit-exact
  reproduction in that case is the `.py` (which pins `force_backend="python"`).
  The `README.txt` states this explicitly, keyed on the pinned backend.

### Custom-network reproduction (AMENDED during implementation)

**Original plan** reproduced a custom network on *both* `.py` and `.ini` via a
renamed `user_nuclear_dir="nuclear"` overlay (`network="<custom_name>"`).
End-to-end verification proved this cannot be bit-exact for a **small-based**
custom network: the MT-era reaction *ordering* is keyed on the literal base
name (`_select_era_reactions`: `cfg.network == "small"` → `ORDER_SMALL`, else
`ORDER_MT`), and `load_reaction_names` hardcodes `small`'s 12 reactions. So a
renamed overlay flips the MT branch, permuting the stiff BDF solve and drifting
observables ~1e-6. (Large-based and table-only-small customisations are
unaffected.) Author decision: fix in the GUI only, no physics/ordering change
(which would also shift `small_parthenope`'s pinned references).

**Resolved design:**

- **`.py`** embeds the *exact* `custom_network` dict the GUI passed to
  `run_bbn` and keeps the **base** network — i.e. `run_bbn(cfg,
  custom_network=<dict>, force_backend=<pinned>)`, the identical call the GUI
  made — so it reproduces the tab bit-for-bit for every custom network,
  including small-based-with-removal. `run_mc` likewise gets
  `custom_network=<dict>`.
- **`.ini`** (C CLI, which cannot carry a `custom_network` dict) uses the
  bundled `nuclear/` overlay: `network="<custom_name>"` +
  `user_nuclear_dir="nuclear"`. Exact when the base network is **not**
  `small`; ~1e-6 for **any** small-based customisation (removed/added
  reactions *or* rate-table edits), because the overlay is loaded under a
  different name than `small` and the MT-era ordering is keyed on the base
  name (documented in `README.txt` and the ini header). The `.py` is always
  exact.
- The `nuclear/` overlay is still bundled (for the `.ini` and for re-import),
  so no rate tables are inlined except inside the `.py`'s embedded dict.

## Code changes

- **`primat/gui/export_params.py`**
  - `python_export_text(...)`: new signature accepting the pinned backend, the
    MC quantity list + sample count, and an optional `custom_network` **dict**
    (embedded verbatim, passed to `run_bbn`/`run_mc`, base network kept).
    Replaces the `mc.mean ± std` block with the `run_bbn`-central +
    `run_mc`-std pattern and the full ratio print. Drops the
    `_CUSTOM_NETWORK_NOTE_PY` path.
  - `ini_export_text(...)`: `custom_network_name` overlay override
    (`network=<name>` + `user_nuclear_dir="nuclear"`); replaces
    `_CUSTOM_NETWORK_NOTE_INI`.
  - New `build_reproduction_zip(params, *, results, backend_used, mc=None,
    cfg=None, custom_network=None, kept_names=None, network_name=None)`:
    orchestrates the combined zip — always writes `primat_gui_run.py`,
    `run_basic_from_gui.ini`, `README.txt`; when `custom_network` is set,
    calls `export_zip(...)` and re-packs its entries under `nuclear/`.

- **`primat/gui/panels.py`**
  - `render_results_panel(run, mc=None, run_params=None, backend_used=None)`:
    add the single "Download reproduction bundle (.zip)" button under the
    Standard-ratios table, built via `build_reproduction_zip`.
  - `render_downloads_panel(...)`: remove the two reproduction (`.py`/`.ini`)
    download buttons. `run_params` was consumed *only* by those buttons, so
    drop the `run_params` parameter from this function entirely.

- **`primat/gui/app.py` (`main`)**
  - Pass `stored_params` (as `run_params`) and `backend_used` into
    `render_results_panel`.
  - Stop passing `run_params` to `render_downloads_panel`.

## Edge cases

- **No custom network:** zip has no `nuclear/`; `.py`/`.ini` omit
  `user_nuclear_dir`/`network` override.
- **No quick-MC:** no `± 1σ` column, no `run_mc` block/import.
- **`export_zip` unavailable data:** it already degrades gracefully
  (skips a missing on-disk table) — no new failure mode introduced.
- **Large network custom:** the zip inlines every kept table (~429 for
  `large`) — large but correct and self-contained, matching today's
  Reactions-tab custom export.

## Testing / verification

- **Unit (`tests/test_gui.py`, `tests/test_gui_custom_network.py`):**
  - zip structure for standard vs custom runs (presence/absence of `nuclear/`).
  - `.py` content: pinned `force_backend`, full standard-ratio print from
    `run_bbn`, `run_mc(..., seed=0, ...)` producing a std-only column, no
    `mc.mean` usage.
  - `.ini` content: custom-network drops base `network`, sets
    `user_nuclear_dir=nuclear` + `network=<name>`.
- **End-to-end parity (must drive, not just unit-test):** run the emitted
  overlay config through `run_bbn` and confirm it reproduces the numbers the
  GUI displayed for the same custom network (base `network` + `custom_network`
  dict). They must agree, including any MT-era reaction whose rate table was
  overridden. If they diverge, that divergence is a bug to fix or explicitly
  document before this is considered complete.

## References

- `primat/gui/export_params.py` — current exporter.
- `primat/gui/custom_rates.py::export_zip` — self-contained overlay zip builder.
- `primat/gui/panels.py::render_results_panel`, `render_downloads_panel`,
  `final_abundances_text`, `_RATIO_FORMAT`, `_RATIO_LABELS`.
- `primat/gui/app.py::_solve` (`backend_used`), `_quick_mc` (`seed=0`).
- `CLAUDE.md` → "Rates directory resolution (overlay)" (`user_nuclear_dir`).
