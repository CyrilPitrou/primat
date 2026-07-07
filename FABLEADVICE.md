# FABLEADVICE.md — Improvement plan for primat (PyPI-grade, community-scale)

Audit date: 2026-07-07, at version 0.3.2 (branch `master`, clean tree).
Written by Fable 5 after a full review of `primat/`, `primat/gui/`,
`primat-c/`, `doc/`, `notebooks/`, `runfiles/`, `tests/`, and the packaging/CI
files. Goal: primat installed by hundreds of cosmologists via `pip install
primat`, robust on every platform, beautiful to read and use, and flexible
enough that physicists reach for it for *every* BBN use case.

**Overall verdict**: the codebase is unusually healthy — rich physics
docstrings, a layered 38-file test suite, backend parity discipline, a real
GUI, fingerprinted caches, an honest PyPI publishing guide. What is missing
is mostly *around* the code: CI that actually runs the Python tests, wheels
that are verified before shipping, an HTML documentation site, accurate
README examples, friendlier error messages, and the ecosystem glue
(Cobaya, CITATION.cff, changelog, conda-forge) that turns a good code into
the community-standard one.

---

## How to use this plan

- **Sonnet** gets tasks that are well-specified, mechanical, or high-volume
  editing with a clear acceptance test. It should not touch physics or
  numerics.
- **Opus** gets tasks that are cross-cutting, involve design decisions, new
  subsystems, or anything that must be mirrored between the Python and C
  backends (per CLAUDE.md's parity mandate).
- Work in phase order (**Phase 0 first** — author decision, 2026-07-07 —
  then Phase 1 before the PyPI release; later phases after).
- Items are numbered `F-n` (the Phase-0 feature), `S-n` (Sonnet) and `O-n`
  (Opus).

### Ground rules for BOTH models (non-negotiable, from CLAUDE.md)

1. Read `CLAUDE.md` first and obey it: reporting precision (Neff to 8
   decimals, D/H to 7, …), heavy physics commenting with paper citations,
   the C↔Python numerics parity rule, verbose-output parity, and the
   three-way sync of `DEFAULT_PARAMS` ↔ `runfiles/primat_run_explanatory.py`
   ↔ `primat-c/examples/run_basic.ini`.
2. After any change that could affect numbers, run
   `python runfiles/primat_run.py` and check against the CLAUDE.md
   "Validation before committing" table.
3. Run `graphify update .` after code changes (project rule).
4. **Never** perform the 🔴 steps of `PyPiGuide.md` (claiming the PyPI name,
   tagging/publishing a release). Those are human-only actions.
5. Documentation edits must be pinned by extending
   `tests/test_docs_consistency.py` wherever a number or claim can be
   machine-checked — that file exists precisely because README/CLAUDE.md
   have staled before.

---

## Phase 0 — MC covariance & correlation matrices (START HERE)

Requested by the author 2026-07-07. Motivation: when a user constrains
cosmology with several abundances at once (mostly YP and D/H), they need the
*joint* nuclear-rate uncertainty — the covariance between observables from
the same MC samples — not just per-observable sigmas. The samples already
exist in `MCResult` (`samples_array()` stacks them as `(num_mc, n_quantity)`),
so this is a presentation feature, but it must land uniformly across the
Python API, CLI, GUI, output files, both backends, and the docs.

Author decisions already made (do not re-litigate):
- Full matrices cover **all** MC quantities (every standard observable +
  every tracked nuclide's final Y — the same set as the samples file), in
  `quantity_names()` order.
- **Hard rename** `output_mc_file` → `output_mc_file_prefix` (no deprecated
  alias; primat is not on PyPI yet).
- CLI `--mc` prints the 4×4 covariance and correlation of the four main
  products: `YPBBN`, `DoH`, `He3oHe4`, `Li7oH`.
- Every behavior change is reflected in README.md **in the same PR**.

### F-1. Core implementation, both backends  «Opus»
**Python API** (`primat/main.py`, `MCResult` — build on `samples_array()`):
- `mc.cov()` → full `(n_q, n_q)` sample covariance matrix (`np.cov` on the
  stacked samples, `ddof=1` — document the convention; it must match how
  `MCResult.std` is computed, check and align if it uses ddof=0).
- `mc.corr()` → correlation matrix (`np.corrcoef`); guard zero-variance
  quantities (a nuclide identical in every sample) → NaN off-diagonal with a
  docstring note, not a RuntimeWarning storm.
- Scalar two-name form: `mc.cov("YPBBN", "DoH")` / `mc.corr("YPBBN", "DoH")`
  (raise KeyError with the available names on a typo).
- Works identically for C- and Python-backend `MCResult`s since both fill
  `values` — the API layer is backend-agnostic by construction.
**Config keys** (`DEFAULT_PARAMS` goes 74 → 76; full three-file sync per
CLAUDE.md — config.py, `primat_run_explanatory.py`, `run_basic.ini`, count
comments, plus the C `CPRConfig` field table / `cpr_config_set_by_name`):
- Rename `output_mc_file` → `output_mc_file_prefix`, default
  `"results/output_mc"`. Written files: `<prefix>_samples.tsv`,
  `<prefix>_covariance.tsv`, `<prefix>_correlation.tsv`
  (`output_mc_samples.tsv` as a default filename is gone).
- New booleans `output_mc_covariance`, `output_mc_correlation` (default
  False), siblings of the existing `output_mc_samples`.
**File format** (author spec — exactly two header lines):
```
# Covariance matrix of the N=100 primat MC samples (seed=0): C[i,j] = sample covariance (ddof=1) of quantities i and j.
quantity	Neff	YPBBN	YPCMB	DoH	...
Neff	...	...
```
i.e. line 1 = one `#` line saying what the file is (include N, seed, and the
estimator convention); line 2 = tab-separated quantity names labelling both
columns and rows; then one row per quantity, name first. Correlation file
identical with its own line-1 wording and unit diagonal.
**Writers**: shared helpers next to `dump_mc_samples` in `primat/backend.py`
(`dump_mc_covariance(mc)`, `dump_mc_correlation(mc)`), used by CLI and GUI.
The **standalone C CLI** (`primat-c` binary: `mc.c`/`cli.c`) must write the
same three files with byte-identical header wording (verbose-output-parity
spirit) — compute cov/corr in C from its own samples.
**CLI** (`primat/cli.py`): rename the `--output_mc_file` flag to
`--output_mc_file_prefix`; when `--mc N` runs, print after the
`value ± sigma` block the 4×4 **correlation** matrix (aligned, 3 decimals)
and the 4×4 **covariance** matrix (`%.3e`) of YPBBN/DoH/He3oHe4/Li7oH,
each with a one-line title; file writing gated by the three booleans.
**Tests**: extend `tests/test_mc.py` (cov symmetric, `diag(corr)==1`,
`diag(cov)==std**2` under the chosen ddof, scalar form == matrix entry,
zero-variance guard); a file round-trip parse test; backend parity — same
header lines and matrix shape from both backends, statistical agreement of
matrices at large-ish N; a C unit test injecting fixed samples and
cross-checking cov/corr values against numpy-computed references.
**README.md** (same PR): update the MC section — new key names, the three
files, a `mc.corr("YPBBN","DoH")` example, and the CLI printout sample.
**Accept**: `primat --mc 100 --set output_mc_samples=True --set
output_mc_covariance=True --set output_mc_correlation=True
--output_mc_file_prefix results/demo` prints the two 4×4 matrices and writes
the three files on both `--backend c` and `--backend python` with identical
headers; validation table (CLAUDE.md) unchanged; full suite green.

### F-2. Demos, notebook, GUI  «Sonnet» (after F-1 merges)
- **`runfiles/primat_mc.py`** (new, author-mandated): heavily-commented demo
  in the house style of `primat_run_explanatory.py` — run
  `run_mc(500, params={...})`, print each observable `value ± sigma` at
  CLAUDE.md precision, print the 4×4 correlation of the main products, show
  the scalar `mc.cov("YPBBN", "DoH")` / `mc.corr("YPBBN", "DoH")` access,
  and write the three `<prefix>_*.tsv` files. Add it to
  `tests/test_runfiles.py`'s script list (use a small `num_mc` via an env
  override or a `--quick` arg so the smoke test stays fast) and to the
  runfiles list in README/CLAUDE.md.
- **`notebooks/MonteCarloRates.ipynb`**: add a section computing and
  displaying the full correlation matrix via `mc.corr()` (pandas-styled
  heat-styled table for the 4 main products; replace any hand-rolled
  `np.corrcoef(values)` with the new API), with a sentence on *why* the
  YP–D/H correlation matters for joint likelihoods. Keep the corner plot.
- **GUI** (`primat/gui/panels.py`): alongside the existing MC-samples
  download, add "Download covariance (.tsv)" and "Download correlation
  (.tsv)" buttons sharing the radical name (`output_mc_samples` /
  `output_mc_covariance` / `output_mc_correlation` naming), fed by the F-1
  `dump_mc_*` helpers; follow the existing `_download_button` pattern and
  add AppTest coverage in `tests/test_gui.py`.
- **Docs-consistency**: pin README's quoted key names
  (`output_mc_file_prefix`, `output_mc_covariance`, `output_mc_correlation`)
  in `tests/test_docs_consistency.py`.
**Accept**: `python runfiles/primat_mc.py` runs standalone from the repo
root; `pytest -m "gui or notebook"` green; README examples reproduce.

---

## Phase 1 — Release blockers (do these before `pip install primat` is real)

### S-1. Add a Python test workflow to GitHub Actions  «Sonnet»
**Problem**: `.github/workflows/` contains only `wheels.yml`,
`c-sanitizers.yml`, and a legacy `build_linux.yml`. **pytest never runs in
CI.** Worse, `tests/README.md` and `pytest.ini` already *claim* CI lanes
("CI: every push, ~3 min", "CI: nightly") that do not exist.
**Do**: create `.github/workflows/tests.yml`:
- On push/PR: matrix `{ubuntu-latest, macos-14, windows-latest} ×
  {3.10, 3.13}` (full 3.10–3.13 on ubuntu only, to keep minutes sane),
  `pip install -e ".[recommended,gui,notebooks]"` + pytest with
  `-m "not slow or solve"` (the documented fast+solve lane).
- Nightly (schedule): full suite including `-m reference` and the `wheel`
  marker on ubuntu.
- Cache pip; upload the pytest junit report as an artifact.
**Accept**: green runs on all matrix legs; `tests/README.md`'s CI claims are
now true (update the wording to match the actual schedule).

### S-2. Smoke-test wheels inside cibuildwheel  «Sonnet»
**Problem**: `setup.py`'s `optional_build_ext` deliberately lets the C
extension fail without failing the build. That is right for source installs
and **wrong for wheels**: a wheel that silently shipped without
`primat._primat_c` would give every user the 25×-slower Python backend, and
nothing in `wheels.yml` would notice (there is no `CIBW_TEST_COMMAND`).
**Do**: in `wheels.yml` add
```yaml
CIBW_TEST_COMMAND: python -c "import primat.backend as b; assert b.HAS_C_BACKEND, 'wheel shipped without C extension'; r=b.run_bbn({'network':'small'}); assert abs(r['YPBBN']-0.247) < 2e-3, r['YPBBN']"
```
(one line, coarse tolerance — it is a smoke test, not a regression test; the
tight regression lives in pytest). Skip-list any platform where this is
expected to fail *only* with an explicit comment and a tracking issue.
**Accept**: a `workflow_dispatch` run of `wheels.yml` is green on all four
runners **including windows-latest** with the assert active. If the Windows
leg fails to build the extension, stop and hand the compiler error to Opus
(O-6) rather than weakening the assert.

### S-3. Fix the factually wrong / stale parts of README.md  «Sonnet»
All verified against the running code on 2026-07-07:
- **Quick start values are wrong** (README lines ~74–76): the snippet
  actually prints `YPBBN = 0.24699911`, `DoH = 2.4350167e-05` (verified by
  running it), not `~0.246915` / `~2.43647e-05`.
- **`--set` syntax is wrong** (line ~132): README shows `--set tau_n 880.1`
  but `primat/cli.py` requires `KEY=VALUE` (`--set tau_n=880.1`) and errors
  otherwise.
- **`run_mc` example uses a nonexistent kwarg** (line ~321):
  `run_mc(params=..., n_samples=100, ...)` — the signature is
  `run_mc(num_mc, quantities=None, params=None, ...)`. The second example
  (line ~369) is correct; fix the first to match.
- **Stale "Python-only features" list** (lines ~170–174): claims
  `output_time_evolution=True` forces the Python backend, but per
  `primat/backend.py`'s module docstring (and CLAUDE.md) it is supported on
  both backends now. The list should be: `extra_rho`, `background=`,
  `decay_era`, MC `prev` reuse across backends.
- **`output_file` default contradicts itself**: the parameter table
  (line ~226) says default `results/output_tables.tsv` (correct —
  `DEFAULT_PARAMS` confirms), but the Output section (line ~419) says
  "`output_file=None` (the default)". Reconcile.
- `He3oH` row has a typo: `((He3+T)/H`.
- **Cobaya section is out of date** ("will be available"): the wrapper
  already exists in the author's separate `primat_tools` GitHub repo —
  rewrite to present tense with the link (details in O-4).
(CLAUDE.md's stale "currently 73 keys", data paths, PDF location, and
version-check claim were already fixed directly by Fable on 2026-07-07 —
but F-1 bumps the count again to 76; keep it true.)
**Accept**: every number/example in README is reproduced by actually running
the shown command; add pins to `tests/test_docs_consistency.py` for the
`--set` syntax and the Python-only-features list (parse README, assert).

### S-4. Complete the PyPI metadata in pyproject.toml  «Sonnet»
Currently missing entirely: `classifiers`, `keywords`, `[project.urls]`, and
a machine-readable license expression. For a package physicists will find
via pypi.org search, this is the shop window.
**Do**:
- `license = "GPL-3.0-or-later"` (SPDX string, PEP 639; the LICENCE text is
  GPLv3) + `license-files = ["LICENCE"]`.
- Classifiers: Development Status, Intended Audience :: Science/Research,
  Topic :: Scientific/Engineering :: Astronomy and :: Physics, Programming
  Language :: Python :: 3.10–3.13, Programming Language :: C, Operating
  System :: OS Independent.
- `keywords = ["BBN", "big bang nucleosynthesis", "cosmology", "primordial
  abundances", "early universe", "helium", "deuterium"]`.
- `[project.urls]`: Homepage, Repository
  (`https://github.com/CyrilPitrou/primat`), Issues, Documentation (docs
  site once O-1 exists), Changelog.
- Give the maintainer (Julien Froustey) an email or drop the field.
**Accept**: `twine check dist/*` clean; `pip show primat` displays the URLs.

### S-5. Repository hygiene for the public repo  «Sonnet»
**Load-bearing exception — do NOT remove**: the committed wheel
`wheels/primat-0.3.2-…whl` and `.github/workflows/build_linux.yml` are the
deployment chain for **primat.streamlit.app**: Streamlit Community Cloud
installs from the repo-root `requirements.txt`, whose last line points at
`./wheels/primat-<version>-cp312-…linux_x86_64.whl`, and `build_linux.yml`
is how that wheel is produced. Instead of deleting:
- Add a `wheels/README.md` and a header comment in `build_linux.yml`
  explaining this chain, so no future cleanup (human or agent) breaks the
  website.
- Modernize `build_linux.yml` without changing its role: `python -m build
  --wheel` instead of the deprecated `python setup.py bdist_wheel`, and
  `actions/setup-python@v5`; pin Python 3.12 to match the wheel tag
  Streamlit Cloud needs.
- Add a check to `tests/test_docs_consistency.py`: the wheel filename in
  `requirements.txt` must exist in `wheels/` and its version must equal
  `pyproject.toml`'s (today they match at 0.3.2; a version bump must not
  silently leave the website on an old wheel — extend the CLAUDE.md
  version-bump checklist accordingly).
- Deduplicate `joblib` (listed twice) in `requirements.txt`.
- Generated outputs in `results/` that a default run overwrites
  (`output_background.tsv`, `output_final.dat`, `output_tables.tsv`): untrack
  them. Keep the two `PRIMAT_Yp_DH_ErrorMC_*_2026.dat` files only if
  something loads them (check `notebooks/` and `runfiles/` first; if
  StandardPlots.ipynb reads them, move them under `notebooks/results/` and
  say so in the notebooks README).
- Ensure `.gitignore` covers `build/`, `*.egg-info/`, `wheelhouse/`,
  `results/*.tsv` (minus any kept references), `__pycache__`.
- `MANIFEST.in`: verify the sdist also ships `primat-c/Makefile`,
  `primat-c/examples/*.ini`, and `tests/` (decide: shipping tests in the
  sdist is good practice for debian/conda packagers). Run
  `python -m build --sdist` and inspect the tarball listing.
**Accept**: fresh clone + `python -m build` + install of *both* artifacts
passes `tests/test_wheel_smoke.py`; `git status` stays clean after a default
run (no tracked file modified by running the code).

### S-6. Restructure wheels.yml for the testpypi→pypi switch  «Sonnet»
`publish` currently hard-points at TestPyPI with the real index commented
out — easy to forget to flip, and flipping requires a commit at the most
delicate moment (PyPiGuide step 6).
**Do**: two publish jobs — `publish-testpypi` (on `workflow_dispatch` with an
input flag, environment `testpypi`, `repository-url: https://test.pypi.org/legacy/`)
and `publish-pypi` (only on `release: published`, environment `pypi`, no
repository-url). Keep `id-token: write` per-job. Document in PyPiGuide.md
that the flip is now automatic and step 6's manual edit is obsolete.
**Accept**: `workflow_dispatch` publishes to TestPyPI only; the real-PyPI
job is demonstrably unreachable except from a published release.

### O-1. Config validation and error-message UX (both backends)  «Opus»
**Problem** (verified): `PRIMATConfig({"Omegabh2": "not_a_number"})` dies
much later with `TypeError: can't multiply sequence by non-int of type
'numpy.float64'` — a physicist's first contact with primat after a typo in a
yaml file will be a stack trace from deep inside the thermodynamics. Unknown
keys only produce a `UserWarning` ("unknown parameter keys ignored"), which
scrolls past unseen in an MCMC log; a typo like `Omegab2h` silently runs the
default cosmology.
**Do** (design carefully, this touches every entry path):
- In `PRIMATConfig.__init__`, validate each user-supplied key's *type*
  against the `DEFAULT_PARAMS` default's type (float-accepting-int, bool,
  str, None-able fields like `amax`/`output_file` need a small per-key
  spec), raising `ValueError`/`TypeError` immediately with the key name, the
  received value, and the expected type/range.
- Range checks where physics demands them (`Omegabh2 > 0`, `tau_n > 0`,
  `numerical_precision > 0`, `amax` positive int, …).
- Unknown keys: keep the warning by default for back-compat, but (a) append
  `difflib.get_close_matches` suggestions ("did you mean 'Omegabh2'?"), and
  (b) add a `strict_params` config key (default False now, consider True at
  v0.4) that upgrades it to a ValueError. `p_<rxn>`/`delta_<rxn>` dynamic
  keys must keep working.
- **Mirror on the C side** (`cpr_config_set_by_name` /
  `cpr_config_validate` in `primat-c/src/config.c` + INI loader), with the
  same messages, per the parity mandate. The wrapper
  (`primat/_primat_c/_wrapper.c`) must surface C-side validation failures as
  Python exceptions, not exit codes.
- Tests: extend `tests/test_config.py` and the C unit test `test_config.c`.
**Accept**: the two probe cases above produce one-line, self-explanatory
errors on both backends; full suite green; no reference value moves.

### O-2. Decide and execute the dependency diet  «Opus» (decision) → Sonnet (execution)
`plotly` and `joblib` are hard dependencies. plotly is a ~10–20 MB install
that only `primat/plotting.py` and the GUI/notebooks need; joblib is only
used by the pure-Python MC path. A cosmologist embedding primat in a
pipeline (or a cluster admin building a container) should be able to install
a lean core.
**Do**: Opus decides between (a) status quo (simplest, one install always
works) and (b) moving `plotly` to a `plots` extra (also included in `gui`
and `notebooks` extras) with a lazy import inside `plotting.py` that raises
a friendly "pip install primat[plots]" message, and making `joblib` a lazy
import with a clear error only when `mc_uncertainty` runs on the Python
backend with n_jobs≠1. Recommendation: **(b)** — but check first that
nothing in the default `run_bbn` path imports plotting, and grep the GUI and
notebooks for import-time assumptions. If (b), Sonnet executes: pyproject
extras, lazy imports, README install matrix, test that `pip install primat`
core in a clean venv can `run_bbn` without plotly present (extend
`test_wheel_smoke.py`).
**Accept**: documented decision in the commit message; if (b), core wheel
installs and solves with neither plotly nor streamlit present.

---

## Phase 2 — Documentation & first impressions (the "beautiful" part)

### O-3. Build the documentation site (Sphinx + Read the Docs)  «Opus» (architecture) + «Sonnet» (migration)
This is the single highest-leverage item in the whole plan. Today the docs
are: a 525-line README, a LaTeX PDF (`doc/primat_documentation_v0.3.1.pdf`),
EXTENDING.md, notebook READMEs, and superb docstrings *that no user can
browse*. Community-scale packages (astropy, CLASS, emcee) live or die by
their HTML docs.
**Opus does**: choose and scaffold the stack — recommended: Sphinx +
`furo` (or `pydata-sphinx-theme`), `myst-parser` for the existing .md
content, `myst-nb` to render the `notebooks/` gallery (executed in CI's
nightly lane only, cached), `sphinx.ext.napoleon` + `autodoc` for the API
reference, `sphinx-copybutton`. Decide the information architecture:
1. Landing page (what primat is, 10-line quick start, citation).
2. Installation (wheels, extras, from-source, conda-forge once it exists).
3. Tutorials (the notebooks, rendered).
4. How-to guides: choose a network / `amax`; custom networks (GUI zip and
   API); rate variations & MC; data_dir/user_nuclear_dir overlays; custom
   NEVO tables; custom backgrounds & `extra_rho`; CLASS/CAMB tables; CLI
   reference (auto-generated from argparse via `sphinx-argparse`).
5. Physics manual: link the PDF; long-term, decide whether to port the LaTeX
   chapters to Sphinx (do NOT do this now — just link it).
6. API reference (autodoc of `primat`, `primat.backend`, `primat.config`,
   `primat.evolution`, `primat.plotting`, `primat.weak_rates.api`).
7. Changelog and Citing pages (no Contributing page — author decision).
Also: wire Read the Docs (`.readthedocs.yaml`), and a `docs` extra in
pyproject.
**Sonnet then does**: the bulk migration — split README's deep sections
(weak-rate cache workflow, rate variation, NEVO overrides, overlay) into
how-to pages; convert docstring formatting glitches that break napoleon;
add intersphinx links to numpy/scipy; screenshot the GUI for its page.
**Accept**: `sphinx-build -W` clean in CI (add to tests.yml); RTD preview
renders; README slims to ~150 lines that link out (see S-7).

### S-7. README facelift  «Sonnet» (after O-3 skeleton exists)
- Badges: PyPI version, Python versions, CI status, docs status, license,
  arXiv:1801.08023 (the Phys. Rep.), later Zenodo DOI.
- A hero figure: the Schramm plot from `notebooks/plots/schramm_plot.jpg`
  (or regenerate as SVG/PNG at repo-friendly size).
- Keep: what it is, install, 10-line quick start, the four ways to run,
  backend table, citation, links to docs for everything else.
- One short "primat vs other BBN codes" positioning paragraph (PArthENoPE,
  PRyMordial, AlterBBN exist; state primat's differentiators — precision
  physics (CCRTh/SD/FM corrections), C speed + Python flexibility, MC
  uncertainty machinery, GUI — without disparaging others).
**Accept**: README under ~200 lines, zero deep parameter tables (those live
in docs), every claim still test-pinned where feasible.

### S-8. Notebook and template polish  «Sonnet»
- `notebooks/README.md` table is missing `AnimatedAbundances.ipynb` and
  `CompareSmallNetworks.ipynb` — add rows.
- **Replace the baryon-density band with the author's choice**: everywhere
  the notebooks (and notebooks/README.md) draw or quote the CMB baryon
  density band — currently "Planck Ω_b h² = 0.02285 ± 0.00016" — use
  **Ω_b h² = 0.022425 ± 0.000136** instead (author decision, 2026-07-07;
  note this is also the code's default `Omegabh2`). Update the band in
  every notebook that hard-codes it (StandardPlots, PosteriorBaryons, …),
  the README prose, and regenerate the affected plots in
  `notebooks/plots/`. Ask the author for the citation to print next to the
  number (Planck/ACT combination?) rather than inventing one, and stop
  labelling it bare "Planck" if the reference is a combination.
- Add "Open in Colab" badges to each notebook header cell (they must then
  `pip install primat` in the first cell, guarded by an
  `importlib.util.find_spec` check so local runs don't reinstall).
- Ensure all notebooks are committed with cleared outputs or intentionally
  kept outputs (pick one policy; recommend cleared + rendered on the docs
  site) and record the chosen policy in `notebooks/README.md`.
- `runfiles/primat_run_explanatory.py` and `run_basic.ini`: spot-check the
  74-key sync (CLAUDE.md mandates it) and fix the "73"/"74" count comments.
**Accept**: `pytest -m notebook` green; docs-consistency test extended to
assert the notebooks README lists every `*.ipynb` in the folder.

### S-9. Public API ergonomics  «Sonnet»
- **Canonical import is `from primat.backend import run_bbn`** (author
  decision): make every doc, docstring, notebook, and template use exactly
  this form — no top-level re-export is required. Optional nice-to-have,
  only if the author later approves: additionally alias `run_bbn`/`run_mc`/
  `HAS_C_BACKEND` in `primat/__init__.py` for discoverability (cheap,
  backwards-compatible; watch for import cycles since `backend` imports
  `main` lazily). Default to NOT doing it.
- Add `primat.__citation__` (BibTeX string) and surface it in `--credits`.
**Accept**: a grep for `from primat import run_bbn` over docs/notebooks
returns nothing; all examples use the `primat.backend` form and run.

### S-10. Type hints + `py.typed`  «Sonnet»
Public-facing modules have almost no annotations (0 annotated returns in
`backend.py`/`main.py`; no `py.typed` marker), so IDE users and Pylance get
nothing. Annotate the *public* surface only: `backend.run_bbn`/`run_mc`
signatures and return types (`dict[str, float | EvolutionResult]` — consider
a `TypedDict` for the result dict), `PRIMAT` methods, `PRIMATConfig.__init__`,
`evolution.py` (mostly typed already via dataclass), `plotting.py`,
`mc_uncertainty`/`MCResult`. Ship `primat/py.typed`, add it to
package-data. Do NOT chase full-package strictness; internal numerics can
stay untyped. Add a lenient `mypy` (or `pyright --level warning`) step to
tests.yml scoped to the annotated modules.
**Accept**: `mypy primat/backend.py primat/evolution.py primat/config.py`
clean under the chosen config; no runtime behavior change.

### S-11. CLI discoverability  «Sonnet»
- Add `primat --list-params`: print all 74 `DEFAULT_PARAMS` keys with their
  defaults and the one-line comments (parse them from `config.py` the same
  way `test_docs_consistency` does, or maintain a dict) — today the `--set`
  escape hatch is documented as "intentionally undocumented", which is
  clever but hostile to the exact power users primat wants.
- `--version` should also report backend availability
  (`primat 0.3.2 (C backend: available)`).
- Mention `--list-params` in the `--help` epilog next to `--set`.
**Accept**: `tests/test_cli.py` covers both; help text stays under one
screen.

---

## Phase 3 — Ecosystem & community

### O-4. Wire up the existing Cobaya wrapper (primat_tools)  «Opus» (small)
The Cobaya wrapper already exists — it lives in the author's separate
`primat_tools` GitHub repository, **not** in this repo, and is intentionally
staying there. Do not implement a wrapper here.
**Do**:
- Rewrite README's Cobaya section from future tense ("will be available")
  to present tense, linking `primat_tools` with a 3-line usage teaser.
- Add a docs-site how-to page (O-3) "Using primat in MCMC chains" that
  points to `primat_tools` and states the version-compatibility contract
  (which primat versions each primat_tools release supports).
- Add a cross-repo compatibility smoke test to the nightly CI lane (S-1):
  `pip install git+https://github.com/CyrilPitrou/primat_tools` alongside
  the current checkout and run its import/smoke entry point, so a breaking
  change in primat's API (e.g. a result-dict key rename) is caught here
  before users hit it in a chain. Skip gracefully if the repo is
  unreachable.
**Accept**: README/docs link works; nightly cross-repo smoke job green.

### O-5. Promote the CLASS/CAMB table generator to a first-class feature  «Opus» (design) → Sonnet
`runfiles/generate_table_CLASS_CAMB.py` (checkpointed MC grid over
Ω_b h²×ΔN_eff) is exactly what CLASS's `bbn` module and CAMB's
`bbn_table` consume, but it is hidden in runfiles and repo-only.
**Do**: move the logic into `primat/tables.py` with a CLI subcommand
(`primat table --format class --out primat_table.dat --num-mc 400 …`),
keep the checkpoint/fingerprint machinery, document grid/precision
trade-offs, and publish a pre-computed reference table as a release asset
(not in the wheel — it is big) with its generation command recorded in the
header. Sonnet wires CLI/tests/docs after Opus fixes the API.
**Accept**: a generated table loads in CLASS (documented recipe) and CAMB;
`primat table` appears in CLI docs; nightly CI generates a tiny 3×3 grid as
a smoke test.

### S-12. Citation & archival infrastructure  «Sonnet»
- `CITATION.cff` (the Phys. Rep. 754 (2018) 1 paper as preferred-citation,
  software authors Pitrou & Froustey, version + DOI fields ready for
  Zenodo).
- Enable Zenodo–GitHub integration note in PyPiGuide (human clicks the
  toggle; Sonnet documents the step and adds the badge placeholder).
- Consider a JOSS paper (`paper.md` skeleton) — leave as a stub PR for the
  author to decide; JOSS gives citable software credit reviewers respect.
**Accept**: `cffconvert --validate` passes; GitHub shows the "Cite this
repository" button.

### S-13. CHANGELOG.md + release process notes  «Sonnet»
Keep-a-Changelog format; back-fill 0.3.0→0.3.2 from `git log` at coarse
granularity (renames, C backend default, GUI custom networks, weak-rate
interpolation unification — the git history is descriptive enough). Add an
"Unreleased" section; the PR template's checklist (S-14) reminds
contributors to add a line. Link from pyproject `[project.urls]` and the
docs.
**Accept**: changelog exists, rendered in docs; PyPiGuide references it in
the release checklist.

### S-14. Issue/PR templates  «Sonnet»
**No CONTRIBUTING.md** (author decision, 2026-07-07 — do not create one;
`tests/README.md` and CLAUDE.md already carry the dev conventions).
Just the GitHub templates:
- Issue templates: bug report (ask for `primat --version`, backend used,
  platform, and the exact params dict), physics question, feature request.
- PR template with a short checklist: validation run against the CLAUDE.md
  reference table, C↔Python parity considered, DEFAULT_PARAMS three-file
  sync done if parameters changed, changelog line added.
**Accept**: templates appear when opening issues/PRs on GitHub.

### O-6. Windows & platform hardening  «Opus»
The MSVC risk is known (PyPiGuide documents the `<complex.h>` rewrite in
`primat-c/src/weak_rates.c`) but **unverified end-to-end**. Also unverified:
paths with spaces/UTF-8 (the data overlay does string path assembly in C —
`cpr_config_resolve_rates_path` with fixed buffers), long-path limits, and
CRLF in INI/table parsing.
**Do**: drive S-2's Windows wheel leg to green; audit `primat-c` for MSVC
warnings (`/W4`) and fixed-size path buffer overflows (add a unit test with
a long `user_nuclear_dir`); confirm the INI parser tolerates CRLF; run the
pytest suite on windows-latest in S-1's matrix and fix fallout (tmpdir,
path-sep, encoding). Anything numeric that differs on MSVC (x87/fp:precise
vs fast) gets documented tolerance treatment, not ad-hoc weakening.
**Accept**: Windows wheels build with the C extension, pass the smoke
assert, and the pytest suite is green on windows-latest.

### O-7. conda-forge feedstock  «Opus» (prep) + human
After the PyPI release: write the feedstock recipe (`meta.yaml` building
from the sdist; the optional-C-extension trick must become a *required*
build there, since conda-forge always has a compiler), submit to
staged-recipes, and document the maintenance flow. The human maintainer
must be listed and approve.
**Accept**: `conda install -c conda-forge primat` works with the C backend.

---

## Phase 4 — Flexibility & delight (the "best BBN code for all use cases" part)

### O-8. Close the C-backend feature gaps  «Opus»
The silent-fallback list (`extra_rho`, `background=`, `decay_era`, MC
`prev`) is small but each entry costs 25× performance exactly when users do
non-standard physics — primat's flagship use case.
Priority order:
1. **`extra_rho` on C**: design a tabulated interface — Python evaluates
   each callable on a dense T-grid once and passes `(T[], rho[])` arrays
   through the wrapper; C interpolates (cubic, same spline module) inside
   `cpr_bg_Hubble`. This preserves the flexible Python API while keeping the
   C solver fast. Parity test: `extra_rho=[lambda T: const]` must equal
   `DeltaNeff`-equivalent runs on both backends.
2. **`decay_era` on C**: port `_integrate_decay_era` (expm-based decay
   propagation; see `tests/` "C Cache I/O" decay-era tests for the physics
   pins).
3. `background=` (custom Background object) can stay Python-only — document
   it as such, it is an inherently-Python extension point.
**Accept**: fallback list shrinks in `backend.py`'s docstring, README, and
docs simultaneously; parity tests added to `test_backend_parity.py`.

### O-9. Per-flavour neutrino degeneracies  «Opus»
`munuOverTnu` applies one ξ to all flavours. The literature (and Parthenope,
PRyMordial) commonly scans ξ_e separately from ξ_μ,τ — ξ_e enters the n↔p
weak rates directly while ξ_μ,τ only gravitate, so the physics is genuinely
different and primat currently cannot express it.
**Do**: add `munuOverTnu_e`/`munuOverTnu_mu`/`munuOverTnu_tau` (defaulting
to `munuOverTnu` for back-compat), thread ξ_e through
`weak_rates/integrands.py` FD kernels and the fingerprint fields
(`_WEAK_RATE_BG_FIELDS`), the plasma ν energy densities, and mirror all of
it in `primat-c/src/weak_rates.c`/`plasma.c`. Extend
`tests/test_spectral_distortions.py`'s ξ tests (±ξ symmetry, YP slope sign)
per-flavour. This is deep physics+parity work — Opus only, heavy citation
duty (Phys. Rep. §5, Froustey et al. for consistency limits).
**Accept**: ξ_e-only run reproduces known YP sensitivity (ΔYP ≈ −0.25 ξ_e to
leading order); all-flavours-equal reproduces today's numbers exactly;
CLAUDE.md/templates/INI updated (74→77 keys, three-file sync!).

### O-10. Sensitivity API  «Opus» (small)
`notebooks/Sensitivity.ipynb` computes ∂ln(obs)/∂ln(param) tables by
finite-differencing solves — a thing referees ask for constantly. Promote it
to `primat.sensitivity.sensitivity_table(params, observables, targets,
rel_step=…)` returning a small dataclass with a `to_markdown()`/DataFrame
view, using the C backend and sharing the central solve. The notebook then
becomes a thin demo of the API.
**Accept**: notebook output unchanged (same numbers) while calling the new
API; documented how-to page.

### S-15. GUI polish  «Sonnet» (with GUI tests updated in the same PR)
The GUI is already strong (custom networks, zip import/export, backend
pinning). Gaps worth closing, all mechanical:
- **Export the full run**: a "Download params as .py / .ini" button so a
  GUI-explored configuration can be reproduced from a script or the C CLI
  (serialize the current params dict into the two template formats — reuse
  `runfiles/primat_run_explanatory.py`'s layout).
- Progress feedback during quick-MC (sample counter, mirrors CLI's `[MC]`).
- A "Copy citation" element next to the credits.
- Surface `HAS_C_BACKEND` and the active backend in the sidebar footer.
Follow the existing `tests/test_gui*.py` AppTest patterns for every added
widget.
**Accept**: `pytest -m gui` green; exported .py file runs standalone and
reproduces the GUI result.

### S-16. Performance & benchmark page  «Sonnet»
`primat-c/examples/baseline_timings.txt` and `make bench` exist but nothing
is published. Add a docs page: typical wall times per network/backend
(small ~x s, large ~y s, MC-100 ~z s on a stated reference machine),
the weak-rate cache warm/cold split (~1.8 s + vegas minutes), and guidance
(when to raise `numerical_precision`, `sampling_*`, what the reference-run
settings cost). Regenerate numbers with `make bench` + a small timing
script rather than copying stale ones.
**Accept**: docs page with a reproducible script in `runfiles/` that
regenerates the table.

---

## Suggested sequencing

| Order | Items | Why |
|-------|-------|-----|
| 0 | F-1 → F-2 | the MC covariance/correlation feature (author priority) |
| 1 | S-1, S-2, S-3, S-4, S-5, S-6 | CI + wheels + truthful README = safe to release |
| 2 | O-1, O-2 | error UX and dependency shape are hard to change after users exist |
| 3 | O-3 → S-7, S-8, S-9, S-10, S-11 | docs site, then everything that links into it |
| 4 | S-12, S-13, S-14, O-6 | community scaffolding + Windows confidence |
| 5 | O-4, O-5, O-7 | Cobaya/CLASS/conda — the adoption multipliers |
| 6 | O-8, O-9, O-10, S-15, S-16 | flexibility and delight |

## Explicitly out of scope for the models
- Any 🔴 PyPiGuide step (PyPI name claim, release tagging/publishing).
- Changing the license.
- Silently "fixing" physics numbers beyond what the author has explicitly
  decided (the Ω_b h² band value in S-8 is decided: 0.022425 ± 0.000136;
  its printed citation still needs the author's input).
- Porting the LaTeX physics manual to Sphinx (link it; porting is a separate
  project).
