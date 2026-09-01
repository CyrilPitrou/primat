# Development notes

primat is a scientific code, read by physicists as much as by programmers, and
it ships two interchangeable backends that must stay in step. Most of what
follows exists to serve those two facts.

If you only read one section, read **"Decisions already made"** — it is the
list of things that have been tried, measured and settled, so that nobody
spends a week rediscovering them.

## Reporting numerical results

When a BBN observable is quoted — in a commit message, a docstring, a note —
use at least this many decimals:

| Observable | Minimum decimals | Format specifier |
|------------|------------------|------------------|
| Neff       | 8                | `%.8f` |
| YP (BBN)   | 8                | `%.8f` |
| He4/H      | 7                | `%.7e` |
| D/H        | 7                | `%.7e` |
| Li7/H      | 6                | `%.6e` |

This is not pedantry. Flags such as `incomplete_decoupling` and
`QED_corrections` move Neff at the level of 1e-2 to 1e-3, and the two backends
agree on D/H only to a budget of 5e-5 — a number quoted to three decimals
cannot show either.

No example values are given here on purpose: they would go stale. The live
ones are in `tests/reference_values.py`, the single source that both
`tests/README.md` and the tests read.

Mind the abundance convention when quoting a per-nuclide number: `Y_i` is the
abundance *per baryon*, `n_i/n_b`, so `sum(A_i Y_i) = 1`. The mass fraction
is `A_i Y_i` — `YP` is `4 Y_He4`, four times the `He4` entry of `Y_final`.

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

**No unverifiable claims about other files.** A comment asserting that
another file behaves some way cannot be checked while reading, and decays
silently. State each cross-file *requirement* once — here or in
`tests/README.md` — and let a test enforce it.

One exception, deliberate: a comment **may** name the C symbol that mirrors
the line it sits on, as in "Mirrored by `cpr_bg_T_of_a` in
`primat-c/src/background.c`". That is a signpost, not a requirement — someone
changing physics needs to find the other side, and a symbol name is the
fastest way to it. The rule above still governs the requirement itself: the
parity is what a test asserts, never the comment. Keep the pointer to a named
symbol in a named file, so a rename breaks it visibly at the next `git grep`
rather than quietly.

**The test suite is documentation too.** Every test says what its goal is, and
`tests/README.md` explains every test file, under the same 15-second rule.

**Shorthand is expanded at first use.** `docs/glossary.md` holds one line per
term this project uses on a reader — `HT`/`MT`/`LT`, `CCR`, `FM`, `SD`,
`CCRTh`, `NEVO`, `T9`, `YP`, `expsigma`, `p_<reaction>`, `amax` — with units
where there are units. The rule in every file, source or prose: expand a term
the first time that file uses it, or link to the glossary; a file that uses
one term many times links once rather than repeating the expansion. Add the
entry to the glossary before using a new abbreviation anywhere else.

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
| The version matches `pyproject.toml` everywhere it is repeated — `CPRIMAT_VERSION`, `CITATION.cff`, `manual/`, the Streamlit wheel | `tests/test_docs_consistency.py`; the full bump procedure is `PyPiGuide.md`'s Step 1 table |
| The documented validation numbers match a live run | `tests/test_runfiles.py`, `tests/test_regression.py` |
| A rate table is valid before it reaches the solver | `tests/test_rate_table_domain.py` |

**Adding a parameter** means three edits in one commit: `DEFAULT_PARAMS` and
`PARAM_GROUPS` in `primat/config.py`, a one-line description in
`gen_param_templates._TEMPLATE_DESCRIPTIONS`, then run the generator and commit
what it writes.

**Changing physics** means changing it twice — once in `primat/`, once in
`primat-c/src/`, where the mirror of a module carries the same name
(`background.py`/`background.c`, `weak_rates/`/`weak_rates.c`, and so on;
`primat-c/README.md`'s "Code structure" lists them). Purely cosmetic changes
(renames, comments, restructuring with no numerical effect) need no mirroring.
If you are unsure which kind you have, treat it as numerical and mirror it.

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

**`rate_grid_npts = 1000` is converged enough, and the residual is known.**
Past about `numerical_precision = 1e-8` the answer is set not by the solver
but by two fixed sampling grids, and both backends carry the same error, so no
cross-backend test can see it. `docs/performance.md`'s "What the default grids
cost" is the authoritative table — quote it rather than re-measuring. The
caveat worth carrying in your head: the last **two** of the six decimals
`Li7/H` is reported to are grid artefacts, not physics, so compare `Li7/H`
only between runs at the same grid settings. Dropping to 500 is visibly
under-resolved.

**`--cache-clear` clears the shipped cache files too, not only a user
overlay.** This is deliberate: recomputing them reproduces the shipped values
to within the solver's own jitter, so "clear the cache" clearing everything is
the honest reading of the command.

**Exported network archives carry user rate tables verbatim, on their original
grid, at full precision — deliberately not resampled onto the master grid.**
That is what makes an export/import round trip reproduce the run exactly.
Resampling at export would round the values and stretch a coarse upload across
a wider grid, which re-importing would then resample a second time.

**Parameter renames were a clean break, with no compatibility aliases.**
`numba_installed` became `use_numba`, `atol_large_LT` became `atol_LT`, and
`rescale_nuclear_rates` — accepted but read by nothing — was deleted. An old
spelling is now an unknown key: a warning normally, an error under
`strict_params`, with the usual "did you mean …?" hint. The judgement behind
that was that this package has virtually no outside users, so a shim would
cost more in confusion than the break costs in migration. Weigh the same
question the same way, or decide it differently on new evidence — but do not
add aliases on the assumption nobody thought about it.

**Do not merge `_ccrth_FD2_vec` into `integrands.FD2`.** They compute the same
thing and look like an obvious duplication. They are not interchangeable: the
two differ by about e^-300 in the far tail, so the merge cannot be promised to
leave every digit unchanged, and a rename alone buys too little to be worth
re-pinning anything.

**The README is deliberately self-contained**, and about two-thirds of it
restates nine `docs/` pages. Splitting it was proposed with the duplication
measured, and declined: a single file that answers a new reader end to end is
worth maintaining twice. The cost is real — keep both copies right when you
change either.

**Three large files stay single files**, and that was priced rather than
assumed: `network_data.py`, `config.py` and `primat/gui/params_form.py`. The
job on each is to make it readable where it stands — section headers, a module
header saying what lives where — not to split it.

`params_form.py` is the one that costs a reader: **it is two modules
interleaved in one file**, priced and left alone rather than split: the
sidebar form's metadata and helpers, then the whole custom-network dialog
layer, then the form's entry point and group renderers. What costs a contributor time is that adding one
parameter group means editing three sites on either side of that dialog layer
— `_FORM_METADATA`, then `GROUP_ORDER`/`_EXPANDED_GROUPS`/`_SUBHEADING`, then
`_render_curated_groups` — with nothing in the file saying so.

**The shipped NEVO table's `T_numu` and `T_nutau` are not identical**, and
that is an input, not a defect. They differ by up to 1.6e-05, alternating sign
at 14 of the table's 2761 nodes, which is NEVO solver noise; the two flavours
are physically degenerate below the muon mass. Both backends read the same
table and the difference is far below the sensitivity of any observable.

**Two off-default nucleon masses look like a hang and are not.** Setting `mn`
about 1 % high (or `mp` 1 % low) puts the neutron–proton mass difference far
from its physical value, which forces a full Monte-Carlo (`vegas`) rebuild of
the finite-temperature correction table: over 26 minutes on the pure-Python
backend before the first observable appears. The same configuration with
`thermal_corrections=False` finishes in seconds. If a parameter sweep appears
to stall, check whether it moved a nucleon mass.

**The Streamlit deployment chain is not legacy and must not be tidied away.**
The public demo installs primat from the repo-root `requirements.txt`, whose
last line points at a committed wheel under `wheels/`, built by
`.github/workflows/build_linux.yml`. Deleting any of the three breaks the
website. When the version is bumped, rebuild the wheel, commit it, and update
the filename in `requirements.txt` in the same change — this is one row of
`PyPiGuide.md`'s Step 1 table, which is the complete list.

**A numba BDF integrator was measured and declined.** The projection was about
3.5x further on the pure-Python backend for roughly 600 lines of new numba
following `primat-c/src/ode_bdf.c`. The objection that decided it is
structural, not numerical: numba is optional here, every njit kernel has a
pure-numpy fallback, so a numba BDF would have to keep scipy's BDF as its own
fallback — two integrators giving slightly different numbers, with the
regression pins holding on only one. The cheaper middle option, fusing the
whole RHS chain into one njit call and keeping scipy's BDF, buys about a third
of that and was not taken either, for the same reason.

**INI files are a C-side feature by decision.** `primat-c` reads a run
configuration from `--ini`; the Python CLI has no equivalent and is not missing
one — `--set KEY=VALUE` and a `params` dict cover the same ground from Python.
Do not describe `--ini` as available on both: two docs pages once promised it
and had to be corrected.

**The vegas-less thermal fallback stays slow.** Without `vegas`, the
finite-temperature correction falls back to `dblquad`, which is correct but
turns a cold thermal cache from minutes into hours. Making it fast means
reimplementing the integral; what was done instead is the warning, which now
states the cost and names both ways out. That is the honest fix, and it stands
until someone wants to write the integrator.

**The C unit suite is not built under MSVC.** `primat-c/Makefile` is
POSIX-only, so `make test` does not run on Windows. The same sources are
compiled under MSVC through the Python extension and exercised by the Python
suite there, so the code is covered; what is missing is only the standalone C
programs. Closing it means a second build file (CMake or nmake), which has not
been judged worth maintaining.

**`debug` is retired at the next version cut, and the other 94 settings stay.**
Every parameter, CLI flag and GUI control was audited once and eleven were put
up for retirement; ten were kept, each earning its place. `debug` is the
exception: `verbose` does the same job on both backends and is pinned by
`tests/test_verbose_parity.py`. Removing a key is a breaking change of the same
class as a rename, so it waits for the version cut rather than landing on its
own.

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
