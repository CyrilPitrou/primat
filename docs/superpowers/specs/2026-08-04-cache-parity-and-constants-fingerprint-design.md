# Cache parity guards and constants fingerprinting — design

**Date:** 2026-08-04
**Status:** approved, ready for implementation
**Branch:** off `review/c-physics-core` (pass 6 stays closed; this is a follow-up)

---

## 1. Motivation

The C and Python backends **share** every on-disk cache: whichever backend
computes a table first, the other one reads it. This was raised as a possible
source of "strange behaviour", with the proposal to add the backend identity to
the cache fingerprints so the two never mix.

The audit below concludes that backend segregation is the wrong lever, and
identifies the two things that are actually broken: nothing keys the caches on
the *physical constants* they were computed from, and one cache file name cannot
represent more than one configuration at a time.

### 1.1 Why not put the backend in the fingerprint

For three of the four cached artifacts the two backends compute the *same*
mathematical object, and now agree to 1e-10–1e-12 (measured, pass 6). The shared
cache is not a defect there — it is what made finding **F5.14** possible: the
shipped `electron_thermo_cache.txt` coming back modified after a test run is the
observation that exposed a 1.0e-4 gap between the two backends' electron
thermodynamics, subsequently closed at the root by **F6.3** (1.0e-4 → 8.9e-12).
Had the caches been segregated by backend, that gap would still be in `plasma.c`
today, silently, each backend reading its own copy.

Backend segregation converts a parity *bug* into a parity *blind spot*. It also
costs: ~55 shipped weak-cache files and 1.0 MB become ~110 and 2.0 MB in a wheel
that is served to Streamlit Community Cloud, git churn doubles, and any user who
exercises both backends pays the multi-minute vegas build twice for any
configuration not shipped.

The valuable half of the original proposal — *a test that asserts both backends
produce the same cache file from the same config* — is adopted in full (§4).

### 1.2 What is actually broken

**(a) No physical constant is in any fingerprint.** `primat.constants.CONST` has
26 fields; `DEFAULT_PARAMS` has none of them in common. None appears in any of
the four fingerprints. Changing `me`, `alphaem`, `gA`, `mn`, `mp`, `Vud`,
`radproton`, `kappa_n/p` or `GF` therefore silently reuses cached weak rates
computed with the old value.

**(b) One cache file name, many configurations.** The weak caches carry their
fingerprint hash **in the filename**, so configurations coexist.
`electron_thermo_cache.txt` is a single fixed name with the hash only in its
**header**, so configurations evict one another. This is the direct cause of
**F6.12** (a full test-suite run leaves the git-tracked shipped file modified)
and of **F6.13**'s latent flakiness.

**(c) The QED pressure tables have no fingerprint at all**, and the Python path
does not even read the constants it should. `qed_pressure.py`'s `_dPa`/`_dPe3`
default to module-level `_ALPHA_FS`/`_ME_MEV`, and `plasma.py:562` calls
`compute_qed_pressure_tables(T_min=1e-3, T_max=1e2, n_pts=500, verbose=False)` —
α and mₑ are **not passed**. The C backend does the opposite: `plasma.c:135`
passes `g_const.alphaem, g_const.me` explicitly. The values are byte-equal today
(verified: `me = 0.51099895`, `alpha = 0.0072973525692838015`), so nothing is
numerically wrong — but the two backends are bound to two independent sources of
truth for the same constant, agreeing by convention only.

**(d) The Python CCRTh table is not reproducible even against itself** — Python's
vegas is unseeded, while C's is deterministically seeded
(`primat-c/src/weak_rates.c:1047`).

### 1.3 Audit tables

Copies of mₑ/α in the tree:

| Copy | Consumed by | Status |
|---|---|---|
| `CONST.me` / `CONST.alphaem` → `cfg.me` / `cfg.alphaem` | weak rates, plasma, background, NEVO, `network_data`, `nuclear_network` | nominal source of truth |
| `qed_pressure._ME_MEV` / `_ALPHA_FS` | QED pressure tables (Python) | "kept identical by hand" per its own comment |
| `network_data`'s local `ALPHA` / `ME_MEV` | QED nuclear-rate rescale | **deliberately** duplicated; `network_data.c:551` documents why (must reproduce Python bit-for-bit) |
| C `g_const.alphaem` / `g_const.me` | everything on the C side | process-wide global, ~30 read sites, no config path |

Which cache depends on which constants:

| Cache | Constants entering it |
|---|---|
| `weak/nTOp_<hash>.txt` | mₑ, α, `mn`, `mp`, `gA`, `Vud`, `radproton`, `kappa_n/p`, `GF` |
| `weak/nTOp_thermal_<hash>.txt` | α (the correction *is* O(α)), mₑ |
| `plasma/electron_thermo_cache.txt` | **mₑ only** — α does not enter; the QED piece is a separate table |
| `plasma/QED_pressure_correction_e{2,3}.txt` | mₑ, α |

Shipped cache inventory (2026-08-04): `weak/` holds 55 files (38 `nTOp_*`, 17
`nTOp_thermal_*`), 1.0 MB; `plasma/` holds 3 files, 304 KB.

---

## 2. Non-goals

- **No backend field in any fingerprint.** The caches stay shared.
- **No doubling of shipped cache data.**
- **`me`/`alphaem` do not become `DEFAULT_PARAMS` keys in this work.** Making
  them genuinely user-settable requires rewiring ~30 `g_const.*` read sites on
  the C side, several in functions with no `cfg` in scope (e.g. the `plasma.c`
  integrand statics), plus template regeneration and de-duplicating the three
  Python copies. That is a separate piece of work of comparable size — see §8.
- **The QED tables keep their fixed filenames.** Nothing user-facing varies for
  them, so they cannot proliferate and there is no eviction churn to fix.

---

## 3. Design

### Part 1 — `constants_hash` in all four fingerprints

Add **one** field, `constants_hash`, to each of the four fingerprint dicts. Its
value is the 16-hex-digit hash of the entire constants struct: on the Python
side `fingerprint_hash(dataclasses.asdict(CONST))` (26 fields, canonical
sorted-key JSON, already the repo's convention); on the C side the same
computation over `g_const` via the existing `cpr_fingerprint_hash` machinery.

Rationale for one broad field rather than a curated per-cache list:

- It covers mₑ and α *and* the seven other constants that affect the weak rates
  today and are unkeyed, with no list to maintain.
- It stays correct when a constant is added to `CONST` later — a curated list
  is exactly the kind of thing that goes stale.
- It over-invalidates: changing `Mpc` rebuilds the weak tables even though it
  cannot affect them. That is **safe** (caches are regenerable by construction),
  whereas under-invalidation is a silent wrong answer. The asymmetry decides it.

Cross-backend hashing of arbitrary doubles is already solved:
`primat-c/src/cache.c:95-106` implements a Python-`repr`-compatible
shortest-round-trip float serializer precisely so both backends hash floats
identically. Adding 26 constants is mechanical.

**Override guard.** Until `me`/`alphaem` are real params (§8), an attempt to
override them must fail loudly rather than diverge silently: today
`params={"me": …}` hits the unknown-key path and is warn-and-ignored on both
backends, so a user's override would be honoured by neither — but a direct
attribute poke (`cfg.me = x`) *is* honoured by most of the Python code and by
none of C. `PRIMATConfig` gains a check that raises if `cfg.me`/`cfg.alphaem`
differ from `CONST` at the point the backend is selected, naming §8 as the
supported route.

### Part 2 — `electron_thermo_cache.txt` → `electron_thermo_<hash>.txt`

Both backends already build an identical 3-field fingerprint
(`plasma.py:765`, `plasma.c:403`); the filename simply gains the hash, mirroring
the weak tree.

- Configurations coexist instead of evicting one another. A non-default run
  writes a *new* file rather than overwriting the shipped one — **this is what
  closes F6.12's git churn and F6.13's root cause.**
- The legacy fixed-name read path is dropped rather than kept for compatibility:
  a miss costs ~0.7 s (`plasma.py:724`), cheaper than carrying compat code.
- `--cache-info` / `--cache-clear` currently walk only `weak/`
  (`cache_utils.list_weak_cache_files`). Since orphan hashes will now accumulate
  under `plasma/` too, they are extended to cover it. *This is scope added
  beyond the original request; it exists so the new files remain cleanable.*

### Part 3 — QED tables: cfg-sourced constants and a fingerprint header

- `plasma.py:562` passes `alpha=cfg.alphaem, me=cfg.me`, matching
  `plasma.c:135`. This removes the second Python source of truth.
- The two files gain a fingerprint header:
  `{format_version, constants_hash, T_min, T_max, n_pts}`, mirrored in
  `qed_pressure.c` / `plasma.c`. A mismatch recomputes and overwrites, exactly
  as the electron-thermo cache already does; a miss costs ~0.3 s.
- Fixed filenames are kept (§2 non-goals). The file is rewritten only when a
  constant changes — precisely when a rewrite is wanted.
- Side effect: the two shipped tracked files gain header lines. Data rows
  unchanged.

`network_data`'s third copy of `ALPHA`/`ME_MEV` is **left alone**: it is a
deliberate, documented divergence (`network_data.c:551`) feeding a rate rescale
that is applied at load time and never cached in a fingerprinted file.

### Part 4 — deterministic Python vegas

Seed Python's vegas the way `weak_rates.c:1047` already seeds C's, so a Python
CCRTh recompute reproduces itself. The shipped thermal tables are **not**
regenerated, so a residual remains: a Python recompute reproduces itself but not
the shipped file, under the same hash. That residual is inherent to a
Monte-Carlo cache, is what the existing `provenance:` header line
(`cache_utils.py:136-148`) exists to record, and is documented in
`weak_rates/cache.py` alongside it. Its size is measured during verification
(§6).

### Part 5 — `tests/test_cache_parity.py`

For each deterministic cache: force **both** backends to recompute into
**separate `cache_dir`s** on a **coarse grid**, then assert

1. **Hash identity** — both backends emit the same `nTOp_<hash>.txt` /
   `electron_thermo_<hash>.txt` filename. This pins
   `cpr_weak_rate_fingerprint` ≡ `_weak_rate_fingerprint` field-for-field, and
   `constants_hash` equality on top; today that equality is asserted only by
   prose in `.claude/review_findings.md`.
2. **Column-wise agreement** at a tolerance held in a named module constant
   whose comment carries the measured value and the date measured.

| Cache | Pass-6 measured | Proposed pin |
|---|---|---|
| `nTOp_<hash>.txt` | 2.5e-10 | 1e-8 |
| `electron_thermo` (`rho_e`, `p_e`, `drho_e_dT`, `dp_e_dT`) | 8.9e-12 | 1e-9 |
| QED `e2`/`e3` | never measured | measure, then pin |
| CCRTh thermal | never measured | **excluded** (vegas cost) |

The coarse grid (low `sampling_nTOp_per_decade`, `weak_rate_cache=False`) keeps
the module inside the default suite. It does not weaken detection: the F5.14
divergence was a per-point quadrature-tolerance floor, not a grid-resolution
effect, so a coarse grid catches it just as well. The density is itself part of
the fingerprint, so both backends still hash the same coarse config and the
hash-identity assertion is unaffected. The module skips wholesale when
`HAS_C_BACKEND` is False.

### Part 6 — one shipped-cache re-key

Parts 1–3 change every fingerprint, so every shipped cache file must be
re-keyed: rename to the new hash and rewrite the header, **numbers unchanged**.
There is precedent — `WEAK_RATE_FORMAT_VERSION`'s v4 note records exactly such a
re-key ("the shipped tables were re-keyed in place"). Parts 1, 2 and 3 are
therefore implemented together so this happens **once**, not three times.

Scope: 55 files in `weak/`, `electron_thermo_cache.txt` →
`electron_thermo_<hash>.txt`, and header lines added to the two QED files.
`runfiles/generate_weak_rate_caches.py` (which knows the expected filenames,
`_expected_filenames`/`_thermal_cache_path`) is updated with them.

### Part 7 — documentation and review-plan bookkeeping

- `README.md`'s caching section ("Working from a git checkout") is rewritten:
  the F6.12 caveat it documents is now fixed for the electron-thermo cache.
- `CHANGELOG.md` entry.
- `.claude/review_findings.md`: a resolution note under the pass-6 section
  recording that **F6.12** (left by decision, documented only) and **F6.13**
  (out of pass-6 scope) are closed here, plus every newly measured tolerance.
- `.claude/review_plan.md`: pass-6 status cell points at this follow-up. **Pass 6
  is not reopened** — it is done and committed; this is downstream work.

---

## 4. Verification

- Python suite and C suite green.
- **`git status` clean after a full suite run** — the observable success
  criterion for Part 2. (Today a full run leaves
  `primat/data/cache_plasma_weak/plasma/electron_thermo_cache.txt` modified.)
- `python runfiles/primat_run.py` observables checked against
  `tests/README.md`'s validation table at the documented precision (Neff 8 dp,
  YP 8 dp, D/H 7 sf, Li7/H 6 sf). No observable may move: Parts 1–3 change
  cache *keys*, not cache *contents*.
- The shipped-vs-seeded-recompute CCRTh spread (Part 4) measured, together with
  its D/H response, against the ±3e-9 D/H pin.
- The new QED cross-backend tolerance measured and recorded.
- `python -m primat.tools.gen_param_templates --check` exits 0 (no params
  change, so this must stay clean).

---

## 5. Risks

| Risk | Mitigation |
|---|---|
| Derived constants (`hbar = 6.62607015/(2π)·1e-27`) must be bit-identical between Python and a C compiler for `constants_hash` to agree; `-ffast-math` or FMA contraction would break it | The Part 5 hash-identity assertion catches it immediately and loudly, rather than silently. Any divergence is a build-flag bug worth surfacing. |
| The re-key touches 55+ tracked data files in one commit | Numbers are unchanged; verify with a column-wise diff of old vs new content before committing, and state the check in the commit message. |
| Packaging: no `MANIFEST.in`/`package_data` glob was found, so how `primat/data/` is packaged must be confirmed before renaming shipped files | Confirm the packaging mechanism first; a rename must not drop files from the wheel. Covered by `tests/test_wheel_smoke.py`. |
| `tests/test_plasma.py` (3 tests) and `primat-c/tests/unit/test_plasma.c` name the fixed electron-thermo filename | Enumerated in the plan; updated with Part 2. |
| Part 5 adds C-backend recomputes to the default suite | Coarse grid keeps it fast; measure the added wall-clock and report it. If it exceeds ~30 s, move to an opt-in marker. |

---

## 6. Deferred: making `me`/`alphaem` user-settable

Recorded here so the follow-up does not have to re-derive it. To make the two
constants genuinely modifiable:

1. Add `me`, `alphaem` to `DEFAULT_PARAMS`, `PARAM_GROUPS`, and
   `gen_param_templates._TEMPLATE_DESCRIPTIONS`; regenerate
   `runfiles/primat_run_explanatory.py` and `primat-c/examples/run_basic.ini`
   (CLAUDE.md's sync rule).
2. Python: point `qed_pressure` at `cfg` (already done by Part 3) and decide
   what to do about `network_data`'s deliberate third copy.
3. C: either add both to `CPRConfig` and rewire ~30 `g_const.*` reads — several
   in functions with no `cfg` in scope — or assign `g_const.me = cfg->me` per
   run, which needs care because `g_const` is a process-wide global shared by
   coexisting configs.
4. Extend `cpr_config_set_by_name` and the INI round-trip.
5. `constants_hash` (Part 1) already makes the cache side correct, so no further
   fingerprint work is needed.
6. Remove the Part 1 override guard once the above lands.
