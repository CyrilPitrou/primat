# Contributing to primat

primat is a scientific code, read by physicists as much as by programmers, and
it ships two interchangeable backends that must stay in step. Most of what
follows exists to serve those two facts.

If you only read one section, read **"Decisions already made"** — it is the
list of things that have been tried, measured and settled, so that nobody
spends a week rediscovering them.

## Reporting numerical results

When you quote a BBN observable — in a commit message, an issue, a docstring,
a pull-request description — use at least this many decimals:

| Observable | Minimum decimals | Format specifier |
|------------|------------------|------------------|
| Neff       | 8                | `%.8f` |
| YP (BBN)   | 8                | `%.8f` |
| He4/H      | 7                | `%.7e` |
| D/H        | 7                | `%.7e` |
| Li7/H      | 6                | `%.6e` |

This is not pedantry. Flags such as `incomplete_decoupling` and
`QED_corrections` move Neff at the level of 1e-2 to 1e-3, and the two backends
agree on D/H only to about 1e-5 — a number quoted to three decimals cannot
show either.

No example values are given here on purpose: they would go stale. The live
ones are in `tests/reference_values.py`, the single source that both
`tests/README.md` and the tests read.

## Comments and docstrings

Comments exist to make the physics legible, which means they must be short
enough to actually be read.

**The 15-second rule.** A reader must absorb any docstring in under 15
seconds. Anything longer belongs in `docs/`, in a module header, or nowhere.

**A docstring has four parts and no fifth:**

1. One line: what it computes.
2. The physics, formula or convention, with its citation — two or three
   sentences at most.
3. Arguments and return value, with units.
4. One usage example, where a caller would otherwise have to guess.

**Inline comments explain *why*, not *what*.** A comment restating the line
below it is noise. Naming and explaining a magic number is a *why*.

**Citations.** The primary reference is
`biblio/Pitrou_etal_PhysReptArxivVersion.pdf`, cited by equation number;
secondary PDFs are in `biblio/`. `generate_rates/PRIMAT-Main.m` is a fallback
only — its equation numbers may be off by one from the published paper.

**No archaeology.** Comments describe the code as it is: not what it used to
do, not what a review found, not which alternative was rejected, not what a
deleted function did not do. That is what commit messages and this file's
"Decisions already made" section are for. A comment whose subject no longer
exists is deleted outright.

**No narrated measurements.** A number measured from a run belongs in the test
that pins it, not in a comment that goes stale. A comment may name the test;
it may not quote the number.

**No unverifiable claims about other files.** "Mirrored by `setup_ede` in
`primat-c/src/background.c`" cannot be checked while reading and decays
silently. State each cross-file requirement once — here or in
`tests/README.md` — and let a test enforce it.

**The test suite is documentation too.** Every test says what its goal is, and
`tests/README.md` explains every test file, under the same 15-second rule.

## What the tests enforce

These rules are not honour-system. Breaking one fails the suite, so you will
find out immediately — but knowing they exist saves you the surprise.

| Rule | Enforced by |
|---|---|
| Physics and numerics exist in **both** backends — every formula, correction, clamp, tolerance, cache-fingerprint field and default | `tests/test_backend_parity.py`; the contract is in `README.md`, "Backend parity contract" |
| Both backends narrate a run identically: same stages, same wording, same order | `tests/test_verbose_parity.py` |
| Both backends produce the same result-dict keys and the same time-evolution columns | `tests/test_backend_parity.py`; the schema contract is `primat/evolution.py`'s module docstring |
| Both backends compute the same cache fingerprints, so they share cache files instead of evicting each other's | `tests/test_cache_parity.py` |
| The two parameter templates match `DEFAULT_PARAMS` | `tests/test_docs_consistency.py`; regenerate with `python -m primat.tools.gen_param_templates` |
| `CPRIMAT_VERSION` matches `pyproject.toml`'s version | `tests/test_docs_consistency.py` |
| The documented validation numbers match a live run | `tests/test_runfiles.py`, `tests/test_regression.py` |
| A rate table is valid before it reaches the solver | `tests/test_rate_table_domain.py` |

**Adding a parameter** means three edits in one commit: `DEFAULT_PARAMS` and
`PARAM_GROUPS` in `primat/config.py`, a one-line description in
`gen_param_templates._TEMPLATE_DESCRIPTIONS`, then run the generator and commit
what it writes.

**Changing physics** means changing it twice — once in `primat/`, once in
`primat-c/src/`. Purely cosmetic changes (renames, comments, restructuring
with no numerical effect) need no mirroring. If you are unsure which kind you
have, treat it as numerical and mirror it.

## Decisions already made

Each of these was measured, not assumed. Re-opening one is fine — but start
from the measurement, not from scratch.

**The two backends do not agree to the last digit, and that is understood.**
`tests/README.md`'s "Known cross-backend divergences" is the authoritative
list: which term causes what, which were aligned, and which were deliberately
left. Two were *tried* and made agreement **worse** — matching the two
high-temperature integrators, and making Python use C's interpolant in
`external_scale_factor` mode. Read that table before "fixing" a divergence.

**No tolerance has ever been loosened to make a failure go away.** The ±3e-9
deuterium regression bound has been that tight since the first commit; only
its central value has moved, for a documented physics change each time. The
cross-backend budget has only ever been tightened. Keep it that way.

**`rate_grid_npts = 1000` is converged.** Measured against 8000: deuterium to
4e-6, He3/H to 8e-6, Li7/H to 8e-5. Neff does not move at 8 decimals. The one
caveat worth knowing is that Li7/H's sixth decimal is *not* grid-converged.
Dropping to 500 is visibly under-resolved.

**`--cache-clear` clears the shipped cache files too, not only a user
overlay.** This is deliberate: recomputing them reproduces the shipped values
to within the solver's own jitter, so "clear the cache" clearing everything is
the honest reading of the command.

**Exported network archives carry user rate tables verbatim, on their original
grid, at full precision — deliberately not resampled onto the master grid.**
That is what makes an export/import round trip reproduce the run exactly.
Resampling at export would round the values and stretch a coarse upload across
a wider grid, which re-importing would then resample a second time.

**The Streamlit deployment chain is not legacy and must not be tidied away.**
The public demo installs primat from the repo-root `requirements.txt`, whose
last line points at a committed wheel under `wheels/`, built by
`.github/workflows/build_linux.yml`. Deleting any of the three breaks the
website. When the version is bumped, rebuild the wheel, commit it, and update
the filename in `requirements.txt` in the same change.

## Before you commit

Run `python runfiles/primat_run.py` from the repository root and check its
numbers against `tests/README.md`'s "Validation reference".

Use the right tolerance column. That script runs at the default precision,
while the published reference numbers come from a high-precision run — the
tight bound belongs to the *reference* column and fails on a healthy tree if
applied to the routine one. `tests/README.md` spells out which is which, and
both are automated, so `pytest tests/` is the real answer.

The fast lane (`pytest tests/ -m "not slow"`) takes well under a minute; the
full suite takes around twenty. Run the full suite before anything that
touches physics.
