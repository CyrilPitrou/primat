# Cache parity guards and constants fingerprinting — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task, inline in the current session. Do **not** dispatch subagents. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the shared C/Python caches *safe* rather than segregating them —
key every cache on the physical constants it was computed from, give the
electron-thermo cache a hash-named file so configurations stop evicting one
another, fingerprint the QED tables, seed Python's vegas, and add a
cross-backend cache-parity test module.

**Spec:** `docs/superpowers/specs/2026-08-04-cache-parity-and-constants-fingerprint-design.md`
— read it first. It carries the audit, the rejected alternatives, and the
rationale for each choice below.

**Branch:** create off `review/c-physics-core`.

## Global constraints

- **CLAUDE.md's backend-parity rule applies to every task here.** Each
  fingerprint change, filename change and verbose message must land in *both*
  `primat/` and `primat-c/src/` in the same task.
- **No observable may move.** Tasks 2–5 change cache *keys*, not cache
  *contents*. Any shift in Neff / YP / D/H means something was recomputed
  differently and must be investigated, not re-pinned.
- **Comment heavily** (CLAUDE.md): docstrings say what is computed and why, with
  units; magic numbers are named and explained; the fingerprint comment blocks in
  `weak_rates/cache.py` are prose-heavy by convention — match that density when
  adding `constants_hash`.
- **Intermediate verification must use a scratch `cache_dir`.** From Task 2
  onward the shipped caches are stale until Task 5 re-keys them; running without
  a redirect would trigger multi-minute vegas rebuilds and dirty the tracked
  tree.
- `python -m primat.tools.gen_param_templates --check` must exit 0 throughout —
  no `DEFAULT_PARAMS` key is added by this plan.

## File structure

| File | Responsibility |
|---|---|
| `primat/cache_utils.py` | **Modify.** `constants_hash()` helper; extend cache listing/clearing to `plasma/`. |
| `primat/weak_rates/cache.py` | **Modify.** `constants_hash` into both fingerprints; bump `WEAK_RATE_FORMAT_VERSION` 4 → 5 with a documented rationale block. |
| `primat/plasma.py` | **Modify.** Electron-thermo hash-in-filename; QED fingerprint; pass `cfg.alphaem`/`cfg.me` at line 562. |
| `primat/qed_pressure.py` | **Modify.** Note that callers now supply α/mₑ; keep the module defaults as a documented fallback. |
| `primat/config.py` | **Modify.** Override guard for `me`/`alphaem`; update the `recompute_electron_thermo` comment. |
| `primat/weak_rates/corrections.py` | **Modify.** Seed vegas deterministically. |
| `primat-c/src/cache.c` | **Modify.** `constants_hash` over `g_const`; mirror both fingerprints. |
| `primat-c/src/plasma.c` | **Modify.** Electron-thermo hash-in-filename; QED fingerprint. |
| `primat-c/src/qed_pressure.c` | **Modify.** QED fingerprint header write. |
| `primat-c/include/cache.h` | **Modify.** Declare the new helper. |
| `tests/test_cache_parity.py` | **Create.** Cross-backend cache-file parity. |
| `tests/test_plasma.py` | **Modify.** 3 tests name the fixed electron-thermo filename. |
| `primat-c/tests/unit/test_plasma.c` | **Modify.** Same. |
| `runfiles/generate_weak_rate_caches.py` | **Modify.** `_expected_filenames` / `_thermal_cache_path` learn the new hashes. |
| `primat/data/cache_plasma_weak/**` | **Re-key.** 55 weak files renamed + headers; electron-thermo renamed; 2 QED files gain headers. |
| `README.md`, `CHANGELOG.md`, `.claude/review_findings.md`, `.claude/review_plan.md` | **Modify.** Task 8. |
| `<scratchpad>/rekey_caches.py` | Throwaway. Task 5's rename/rewrite script. Not committed. |

---

### Task 1: Baseline and recon

No code changes. Establish what must not move, and confirm two unknowns.

- [ ] **Step 1: Record the baseline.** From a clean tree, run
  `python runfiles/primat_run.py` and record Neff, YPBBN, He4/H, D/H, Li7/H at
  CLAUDE.md's precision. Run the full Python suite and the C suite; record
  pass/fail counts and `git status` afterwards (expect
  `electron_thermo_cache.txt` modified — that is the F6.12 behaviour this plan
  removes). Restore the tree with `git checkout` before proceeding.
- [ ] **Step 2: Confirm the packaging mechanism.** No `MANIFEST.in` or
  `package_data` entry was found for `primat/data/`. Determine how the data tree
  reaches the wheel (`setup.py`, `pyproject.toml`, or a setuptools default) and
  confirm a *renamed* file is still included. This gates Tasks 4 and 5.
- [ ] **Step 3: Enumerate touch sites.** Every reference to
  `electron_thermo_cache.txt` and to `QED_pressure_correction_e{2,3}.txt` across
  `primat/`, `primat-c/`, `tests/`, `runfiles/`, and the docs. Write the list
  into the task notes; Tasks 3 and 4 work from it.

**Verify:** baseline numbers and suite counts recorded; packaging mechanism
named; touch-site list complete.

---

### Task 2: `constants_hash` in all fingerprints (both backends)

**Files:** `primat/cache_utils.py`, `primat/weak_rates/cache.py`,
`primat/plasma.py`, `primat-c/src/cache.c`, `primat-c/include/cache.h`,
`primat-c/src/plasma.c`

- [ ] **Step 1: Python helper.** Add `constants_hash()` to `cache_utils.py`:
  `fingerprint_hash(dataclasses.asdict(CONST))` over all 26 fields, memoised.
  Docstring explains *why one broad hash rather than a curated per-cache list*
  (spec §3 Part 1: over-invalidation is safe, under-invalidation is a silent
  wrong answer) and notes that it deliberately over-covers.
- [ ] **Step 2: C helper.** Mirror in `cache.c` over `g_const`, using the
  existing `cpr_fingerprint_hash` / `CPR_DOUBLE` machinery — `cache.c:95-106`
  already serialises doubles as Python `repr` does, so the two hashes must come
  out equal. Declare in `cache.h`.
- [ ] **Step 3: Wire into all four fingerprints**, Python and C together:
  `_weak_rate_fingerprint`, `_thermal_fingerprint`, the electron-thermo
  fingerprint (`plasma.py:765` / `plasma.c:403`), and — created in Task 3 — the
  QED one.
- [ ] **Step 4: Bump `WEAK_RATE_FORMAT_VERSION` 4 → 5** with a v5 rationale
  paragraph in the existing comment block, matching the prose density of the
  v1–v4 entries: no physical constant was keyed, so a change to `me`, `alphaem`,
  `gA`, `mn`, `mp`, `Vud`, `radproton`, `kappa_n/p` or `GF` silently reused
  stale rates.
- [ ] **Step 5: Assert the two hashes agree.** A throwaway check that Python's
  `constants_hash()` and C's produce the same 16 hex digits. If they differ,
  stop — it means a derived double (`hbar`) is not bit-identical between the two
  builds (spec §5), which is a build-flag bug, not something to work around.
- [ ] **Step 6: Override guard.** `PRIMATConfig` raises if `cfg.me` or
  `cfg.alphaem` differs from `CONST`, with a message naming the deferred work
  (spec §6) as the supported route. Mirror the check on the C side where a
  config crosses the ABI.

**Verify:** both backends' `constants_hash` identical; a scratch-`cache_dir` run
of each backend produces the *same* new `nTOp_<hash>.txt` filename; suites green
(with the scratch redirect).

---

### Task 3: QED tables — cfg-sourced constants and a fingerprint header

**Files:** `primat/plasma.py`, `primat/qed_pressure.py`,
`primat-c/src/plasma.c`, `primat-c/src/qed_pressure.c`

- [ ] **Step 1:** `plasma.py:562` passes `alpha=cfg.alphaem, me=cfg.me`,
  matching `plasma.c:135`. Update `qed_pressure.py`'s comment block: the module
  constants remain as a standalone-use default, but the solver now supplies them
  from `cfg`, so the "kept identical by hand" hazard is gone for solver runs.
- [ ] **Step 2:** Give both files a fingerprint header
  `{format_version, constants_hash, T_min, T_max, n_pts}`, written via
  `write_cache_with_fingerprint` and mirrored in `qed_pressure.c`/`plasma.c`.
  A mismatch recomputes and overwrites, exactly as the electron-thermo cache
  already does (~0.3 s).
- [ ] **Step 3:** Keep the fixed filenames (spec §2) and keep the legacy
  `QED_tables.txt` / three-file fallbacks working — they are header-less, so
  they read as "unknown fingerprint" and must count as a miss, not a crash.
- [ ] **Step 4:** Leave `network_data`'s third `ALPHA`/`ME_MEV` copy untouched;
  add a one-line cross-reference comment pointing at `network_data.c:551`'s
  explanation so a future reader does not "fix" it.
- [ ] **Step 5: Verbose parity** — if a QED cache-miss message is added on one
  backend, add the matching `cpr_log` on the other.

**Verify:** both backends recompute the QED tables into a scratch `cache_dir`
and produce byte-identical fingerprint headers; measure and record the
cross-backend column agreement on `e2`/`e3` (this number is new — it becomes the
Task 7 pin).

---

### Task 4: Electron-thermo cache gets its hash in the filename

**Files:** `primat/plasma.py`, `primat-c/src/plasma.c`, `primat/config.py`,
`primat/cache_utils.py`, `tests/test_plasma.py`,
`primat-c/tests/unit/test_plasma.c`

- [ ] **Step 1:** `electron_thermo_cache.txt` → `electron_thermo_<hash>.txt`,
  built from the existing 3-field fingerprint (now 4 with `constants_hash`), on
  both backends. Drop the legacy fixed-name read path — a miss costs ~0.7 s
  (`plasma.py:724`), cheaper than compat code.
- [ ] **Step 2:** Extend `list_weak_cache_files` / `clear_weak_cache` (or add
  plasma-side siblings) so `--cache-info` and `--cache-clear` see the `plasma/`
  tree, since orphan hashes now accumulate there. Update the CLI help text.
  *This is scope added beyond the original request — see spec §3 Part 2.*
- [ ] **Step 3:** Update the three `tests/test_plasma.py` tests and
  `primat-c/tests/unit/test_plasma.c` that name the fixed filename. The
  backup/restore dances in those tests exist only because the shipped file was
  being overwritten; with hash-named files a non-default config writes a *new*
  file, so simplify rather than transliterate.
- [ ] **Step 4:** Update `config.py`'s `recompute_electron_thermo` comment
  (it names the old filename).

**Verify:** two configs with different `T_start_cosmo_MeV` produce two coexisting
files instead of evicting each other; suites green under a scratch redirect.

---

### Task 5: Re-key the shipped caches, once

**Files:** `primat/data/cache_plasma_weak/**`,
`runfiles/generate_weak_rate_caches.py`

The numbers do not change — this is a rename plus a header rewrite, so it must
be done by transformation, never by recomputation.

- [ ] **Step 1:** Write `<scratchpad>/rekey_caches.py`: for each of the 55
  `weak/` files, parse the existing `# fingerprint:` JSON, insert
  `constants_hash` and the bumped `format_version`, recompute the hash, rewrite
  the header, rename the file. Same for the electron-thermo file. Add fresh
  headers to the two QED files from the known `T_min=1e-3, T_max=1e2, n_pts=500`.
- [ ] **Step 2:** Run it, then **prove the data rows are untouched**: a
  column-wise diff of every old/new pair must be exactly zero. Record the check.
- [ ] **Step 3:** Update `runfiles/generate_weak_rate_caches.py`'s
  `_expected_filenames` / `_thermal_cache_path` to the new names.
- [ ] **Step 4:** Confirm a default run now hits every cache (no recompute) with
  **no** `cache_dir` redirect, and leaves `git status` clean.

**Verify:** `python runfiles/primat_run.py` reproduces Task 1's baseline
observables *exactly*; `git status` clean after a full suite run — the headline
success criterion; `git diff --stat` shows renames and header lines only.

---

### Task 6: Deterministic Python vegas

**Files:** `primat/weak_rates/corrections.py`, `primat/weak_rates/cache.py`

- [ ] **Step 1:** Seed Python's vegas the way `primat-c/src/weak_rates.c:1047`
  seeds C's, deriving the seed from the fingerprint so it is stable per
  configuration. Document the parallel in a comment on both sides.
- [ ] **Step 2:** Do **not** regenerate the 17 shipped thermal tables. Document
  the residual in `weak_rates/cache.py` next to the existing `provenance:` note:
  a Python recompute now reproduces *itself*, but not the shipped file, under
  the same hash — inherent to a Monte-Carlo cache.
- [ ] **Step 3: Measure the residual.** Recompute one shipped thermal table with
  the seeded vegas, report the max relative spread against the shipped values,
  and the resulting D/H shift against the ±3e-9 pin.

**Verify:** two successive Python recomputes of the same thermal table are now
bit-identical; the measured residual and its D/H response recorded.

---

### Task 7: `tests/test_cache_parity.py`

**Files:** `tests/test_cache_parity.py` (new), `tests/README.md`

- [ ] **Step 1:** New module, skipped wholesale when `HAS_C_BACKEND` is False.
  Each test drives both backends with `weak_rate_cache=False` and a per-backend
  `cache_dir=tmp_path/{c,py}` on a coarse grid (low
  `sampling_nTOp_per_decade`).
- [ ] **Step 2: Hash identity.** Both backends must emit the *same* filename for
  `nTOp_<hash>.txt` and `electron_thermo_<hash>.txt`. This is the assertion that
  pins the two fingerprint implementations — including `constants_hash` —
  field-for-field, replacing prose in `.claude/review_findings.md`.
- [ ] **Step 3: Column agreement**, tolerances as named module constants whose
  comments carry the measured value and the date:

  | Cache | measured | pin |
  |---|---|---|
  | `nTOp_<hash>.txt` | 2.5e-10 (pass 6) | 1e-8 |
  | electron-thermo, all four columns | 8.9e-12 (pass 6) | 1e-9 |
  | QED `e2`/`e3` | from Task 3 | set from measurement |

  CCRTh is **excluded** — vegas cost. Note that exclusion in the module
  docstring so the gap is deliberate and visible.
- [ ] **Step 4:** Module docstring explains *why the caches are shared rather
  than segregated by backend*, and that this module is what makes sharing safe
  (spec §1.1). Add the group to `tests/README.md`.
- [ ] **Step 5:** Measure the added wall-clock. If it exceeds ~30 s, move the
  slowest case behind an opt-in marker and say so.

**Verify:** the module passes; deliberately perturbing one backend's constant in
a scratch build makes it fail (confirm it actually detects what it claims to).

---

### Task 8: Documentation, findings, and final verification

**Files:** `README.md`, `CHANGELOG.md`, `.claude/review_findings.md`,
`.claude/review_plan.md`

- [ ] **Step 1:** Rewrite README's caching section ("Working from a git
  checkout") — the F6.12 caveat it documents is now fixed for the
  electron-thermo cache. Say what `constants_hash` does and why a constants
  change now invalidates everything.
- [ ] **Step 2:** `CHANGELOG.md` entry: format-version bump, the re-key, the
  electron-thermo rename, the QED fingerprint, seeded vegas.
- [ ] **Step 3:** `.claude/review_findings.md` — resolution note under pass 6
  recording that **F6.12** and **F6.13** are closed here, with every measured
  number from Tasks 3, 6 and 7.
- [ ] **Step 4:** `.claude/review_plan.md` — pass-6 status cell points at this
  work. **Do not reopen pass 6**; it is done and committed.
- [ ] **Step 5: Final verification.** Full Python suite + C suite green;
  `git status` clean afterwards; `runfiles/primat_run.py` matches Task 1's
  baseline at CLAUDE.md's precision; `gen_param_templates --check` exits 0;
  `graphify update .`.

**Verify:** every box above ticked, with the actual command output quoted — not
asserted.
