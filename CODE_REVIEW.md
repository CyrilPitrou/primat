# PRIMAT — Codebase Review (2026-07-15)

An independent, whole-repo assessment: weaknesses, risks, and concrete
improvement suggestions, each tagged with the Claude model best suited to
execute it. Overall verdict first: **this is a healthy codebase.** The test
suite is layered and self-documenting, CI covers a 8-job OS×Python matrix plus
nightly regression and C sanitizers, the C code uses `snprintf` throughout (no
`strcpy`/`sprintf`), packaging is lean with well-reasoned extras, and the
docstring/comment culture demanded by CLAUDE.md is actually followed. The
findings below are about reducing friction and long-term maintenance cost, not
about anything being broken.

**How to read the model tags.** Each suggestion names the cheapest model that
can do the job safely:

- **Haiku** — mechanical, low-risk edits (docs, config, renames, deletions);
  verify with the existing test suite.
- **Sonnet** — well-scoped refactors and new tooling where the tests fully pin
  the behavior; no physics judgment needed.
- **Opus** — cross-cutting refactors, C-side changes, performance work; needs
  to hold both backends in mind but not make physics decisions.
- **Fable/Mythos-class (or you)** — anything that touches numerics, weak-rate
  fingerprints, backend parity tolerances, or architecture decisions with no
  test to hide behind.

---

## 1. `PRIMATConfig` as a god object

**Problem.** 191 graph edges — everything talks to the config. 80 flat keys
mix physical constants, solver tolerances, IO paths, GUI conveniences, and
per-reaction `p_*` variations routed through dynamic `__getattr__`. This is
the single biggest obstacle for a newcomer, and dynamic attribute routing
defeats IDE completion and mypy.

**Suggestions.**
- Don't break the flat dict (it's the C-backend contract and user API), but
  *document the taxonomy in code*: group DEFAULT_PARAMS into clearly titled
  sections (largely done via comments) and add a machine-readable
  `PARAM_GROUPS` dict so the GUI, CLI `--list-params`, and templates all
  derive their grouping from one place instead of three hand-maintained
  copies. — **Sonnet**
- Consider a `__dir__` implementation and a generated `.pyi` stub listing all
  80 keys + `p_*` pattern, so completion works despite `__getattr__`. —
  **Sonnet**
- A deeper split (config → PhysicsParams/SolverParams/IOParams) is *not*
  recommended: the churn across both backends and every test would outweigh
  the benefit. Decision to revisit only if the key count doubles. — **you**

## 2. The three hand-synchronized param templates

**Problem.** CLAUDE.md's biggest standing chore: every DEFAULT_PARAMS change
must be mirrored by hand in `runfiles/primat_run_explanatory.py` and
`primat-c/examples/run_basic.ini`, plus two key-count comments plus CLAUDE.md
itself. `test_docs_consistency.py` catches *omissions* but the work is still
manual and the counts ("currently 80 keys") rot.

**Suggestion.** Invert it: write a small generator
(`generate_rates/gen_param_templates.py` or `python -m primat.tools.templates`)
that emits both template files from `config.py`'s
`_default_params_comments` (which already exists and is imported by the CLI).
The templates become build products you regenerate in the same commit; the
existing test then diffs committed vs. freshly generated text, making drift
impossible rather than merely detected. Also drop the literal key-count
sentences in favor of "all DEFAULT_PARAMS keys" so nothing rots. — **Sonnet**
(generator + test), CLAUDE.md rewording — **Haiku**

## 3. Stale and overgrown project documentation

**Problems.**
- CLAUDE.md references `FINAL.md` as the current roadmap — **the file no
  longer exists in the tree** (only `CHANGELOG.md`, `CLAUDE.md`, `PyPiGuide.md`,
  `README.md` remain). Anyone (human or model) following that pointer hits a
  wall. Either restore the roadmap or delete the paragraph.
- CLAUDE.md is 509 lines and increasingly reads as an accreted changelog
  (root-caused bug histories, superseded key counts, long narrative
  parentheticals). Its job is *durable* guidance; history belongs in
  `CHANGELOG.md`/git. A 40% trim would make the remaining rules more likely to
  be followed.
- `docs/superpowers/{plans,specs}` — process artifacts (a 2026-07-10 plan and
  design doc) are committed inside the Sphinx docs tree. Move them to a
  `dev/` or `docs/dev/` folder excluded from the built docs, or delete them if
  the feature shipped.

**Suggestions.** Fix the FINAL.md pointer and relocate the plan/spec files —
**Haiku**. The CLAUDE.md trim requires judgment about which history is still
load-bearing — **you** (or Opus drafting, you approving).

## 4. Error handling

**Problem.** 12 `except Exception` sites. The 8 in `gui/` are defensible
(never crash the app; most log the exception). The ones to audit:
`network_data.py:1432` and `:1758` (swallowing inside rate-table loading — a
malformed table should fail loudly, not degrade silently),
`plasma.py:643/675` (cache-write warnings — fine if they warn, verify they
name the path), `network_builder.py:427` (numba-absent fallback — fine).

**Suggestion.** Narrow the two `network_data.py` handlers to the specific
expected exceptions (`OSError`, `ValueError`) and make sure each path either
re-raises or emits a warning naming the file. — **Sonnet**, with **you**
confirming which failures are legitimately recoverable.

## 5. Performance

**Observations.** `import primat` is 0.41 s with the C backend available —
fine. The C backend is the default fast path, so Python-side speed matters
mainly for the fallback, MC workers, and the GUI.

**Suggestions.**
- **Numba cache coverage:** only 3 `njit(cache=True)` sites
  (`network_builder.py`); the ~15 FD integrand kernels in
  `weak_rates/integrands.py` and the four electron-thermo integrands in
  `plasma.py` are JIT-compiled without `cache=True`, so every fresh process
  pays the compile again (relevant for `joblib` MC workers and the Streamlit
  app's process). Measure first (time a cold `ComputeWeakRates` with/without),
  then add `cache=True` where the closure pattern allows it. — **Opus**
  (measure + apply; closures over config values can silently defeat the cache,
  which is why this isn't a Haiku job)
- **Rate-table loading:** every `PRIMAT.__init__` re-reads and re-resamples
  every rate table onto the 1000-point master grid (~429 reactions for
  `large`). If profiling shows this matters for MC or GUI latency, add a
  fingerprinted compiled-network cache alongside the existing weak/plasma
  caches. Profile before building — the C backend may already make this moot.
  — **Opus**, go/no-go by **you**
- Don't chase Python-side solver speed further; the C backend exists precisely
  so the Python code can stay readable. Treat readability as the Python
  backend's performance metric.

## 6. C backend robustness

**Observations.** Genuinely good shape: no `strcpy`/`sprintf`, 163 `snprintf`
calls, a sanitizer CI lane, its own unit tests. Remaining nits:

- 258 `malloc/calloc` sites vs. ~93 NULL-comparison checks — some allocation
  failures are likely unchecked. For a scientific CLI, aborting on OOM is
  acceptable; silently dereferencing NULL is not. Add a tiny
  `cpr_xmalloc`/`cpr_xcalloc` (log + exit) and sweep call sites where no
  graceful path exists; keep explicit checks where a graceful degrade is
  intended (cache writers). — **Opus**
- `snprintf` truncation: path-building helpers
  (`cpr_config_resolve_rates_path` etc.) take `outsize` — verify each caller
  checks the return or that truncation is detected once centrally, and add a
  unit test with an absurdly long `data_dir`. — **Sonnet**

## 7. Repo hygiene

**Observations.** Pack size 32.5 MiB — acceptable, but three growth vectors:

- `wheels/` currently holds one 6.4 MB wheel (documented Streamlit-deploy
  requirement — keep), but each version bump *replaces* it in the worktree
  while git history keeps every old blob forever. This is the price of the
  Streamlit chain; note it in `wheels/README.md` and consider whether the
  Streamlit app could instead install from PyPI once releases are regular
  (that would end the committed-wheel chain entirely — a decision for **you**).
- `biblio/` PDFs (~4 papers) and `manual/*.pdf` figures: fine, they're the
  citation backbone CLAUDE.md depends on. No action.
- `runfiles/__pycache__/` exists on disk but is properly gitignored; nothing
  tracked. No action.

## 8. Smaller items (one Haiku afternoon)

- `primat/weak_rates/__init__.py` participates in a 3-file import cycle with
  `api.py`/`corrections.py` (the top-level `primat` cycle is already handled
  by documented lazy imports). Flatten by having `corrections.py` import the
  FD setup directly from `integrands.py` instead of via the package
  `__init__`. — **Haiku**, run the weak-rates tests after.
- `cli._build_parser` (217 lines, 28 `add_argument`s): generate the `--set`
  help and param listing from `_default_params_comments` instead of prose
  duplication where any remains. — **Sonnet**
- `CHANGELOG.md` exists — good; make sure the 0.3.2→next entry captures the
  parameter removals visible in recent commits (`rate_interp_order`). —
  **Haiku**

---

## Suggested execution order

1. **Docs truth restored** (item 3: FINAL.md pointer, plans/specs relocation) — Haiku, minutes.
2. **Template generator** (item 2) — Sonnet; kills the most annoying recurring chore.
3. **Exception narrowing** (item 4) — Sonnet.
4. **C-side xmalloc + truncation tests** (item 6) — Opus.
5. **Numba cache measurement** (item 5) — Opus.
6. **CLAUDE.md trim** (item 3) — you, with a model drafting.

Items deliberately *not* recommended: splitting `PRIMATConfig`, restructuring
`corrections.py`/`background.py`, any Python-solver micro-optimization, and
any change to tolerances or fingerprints without your explicit involvement.
