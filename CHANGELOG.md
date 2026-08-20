# Changelog

All notable changes to `primat` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project has
not yet reached a stable public API (`0.x`), so backwards-incompatible
changes may still land in minor bumps.

Entries here are coarse-grained by design — one line per user-visible change
or theme, not a mirror of `git log`. For full detail on any entry, `git log`
in this repository is the authoritative source.

## [Unreleased]

### Fixed
- **A string parameter set to `None` no longer crashes the C backend.**
  `--set network=None`, or `network = none` in an `.ini` file, left the field
  empty and the run died on a segmentation fault (exit 139) printing nothing,
  where the pure-Python backend raises a `TypeError` naming the key. The C
  backend accepted `None` for every string parameter because most of them —
  the output and data-directory paths — legitimately take it; the two that do
  not now reject it in the same words Python uses.
- **A line longer than a reader's buffer is no longer parsed in pieces.**
  Reading a data file line by line handed an over-long line back in chunks,
  each treated as a line of its own, so the tail of a long `#` comment header
  could arrive without its `#` and be read as a data row — silently, and only
  for a file that happens to straddle the buffer. Comment and blank lines of
  any length are now skipped whole, matching what the Python backend's
  `numpy.loadtxt` does, and an over-long *data* line is an error naming the
  file, the line and the limit. Applies to rate tables, NEVO tables, cache
  files, network lists, `decays.txt`, the reaction catalogue CSVs,
  `nuclides.csv` and `.ini` files.
- **The two pure-Python fast paths now decline a scipy that has moved, instead
  of breaking.** Both are opt-in substitutions for a scipy internal and both
  document a fallback, but only one of the two checks was real: the BDF dense
  LU replacement guarded its import and not the three solver attributes it
  patches, so a scipy that keeps the class and restructures those internals
  failed inside the solve; and the linear-interpolant evaluator let an
  exception from its own build-time probe escape. Neither changes any number
  on a scipy where the fast path applies.
- **The test suite runs on a lean install again.** `joblib` and `pandas` are
  optional dependencies, but three test modules imported them unconditionally,
  so `pytest` on a core `pip install primat` failed rather than skipping.

- **A cache file that could not be written completely is no longer installed,
  and a damaged one no longer crashes the C backend.** Running out of disk
  space left a truncated cache carrying a valid header and no data rows: the
  writing run reported correct numbers and exited 0 in silence, and every
  later run of that configuration then aborted the process outright — taking
  the calling Python session with it when reached through the compiled
  extension. Short writes are now detected and the file discarded with a
  warning, and a cache that fails to parse despite a matching header is
  recomputed, as the pure-Python backend already did.
- **A failed cache write now says so without `--verbose`.** The C backend
  reported it only in verbose output, and on standard output, where it could
  corrupt a `--json` pipe; it now warns on standard error like the Python
  backend, naming the file and how to redirect the cache.
- **A truncated n↔p cache is reported the same way by both backends**, naming
  the file and the offending row, instead of a raw numpy message mentioning an
  argument the caller never passed.
- **Monte Carlo and concurrent runs no longer race over shared state.** The
  physical constants every configuration is built from were rewritten by each
  worker thread rather than read; the Gauss–Legendre integration nodes were
  built without synchronisation; two threads populating one cache shared a
  single temporary file; and the progress display's stop flag was read outside
  its mutex. All four are now safe, and ThreadSanitizer reports nothing where
  it previously reported 53 races.
- **Two configurations differing in `numba_installed` can be built from
  different threads.** The just-in-time kernel rebinding is process-wide and
  marked itself complete before it was, so a second thread could re-wrap an
  already-compiled function and die on a numba error naming nothing about
  primat.
- **Memory is released on the failure paths too.** An unreadable n↔p cache
  leaked the background and neutrino-decoupling tables already built (828 KB
  per attempt, which accumulates in a long-lived server), and a bad
  `--data_dir` leaked the configuration it had begun to fill.
- **The caches that outlive a run are bounded.** The reaction-catalog cache and
  the graphical interface's network-label fallback both grew without limit in a
  process handing out many data directories or network names.
- **A proton and neutron closer in mass than an electron no longer hangs the C
  backend.** With `mn - mp` at or below `me` every n↔p rate integrand is
  imaginary, and the C adaptive quadrature — whose stopping test is false for
  NaN — recursed to its full depth, running for hours at 100 % CPU with no
  output, reachable from the default backend. Both backends now reject the
  combination up front, in the same words.
- **`external_scale_factor` runs are ~80× more accurate, and the two backends
  now agree in that mode.** Both built the `T(a)` inverse and the `t(T)` the
  nuclear network reads on the background ODE's output grid, whose spacing
  carries a second-order error; the two errors had opposite signs, leaving the
  backends 1.2e-05 apart in D/H — 12× their converged budget — with the gap
  *growing* as the tolerance was tightened. Both now refine that grid, which
  costs nothing in this mode because `a(T)` is a table read rather than an ODE
  solution. Default-mode results are untouched.
- **The three thermal-correction accuracy knobs do something again.**
  `vegas_n_eval`, `vegas_n_itn` and `epsrel_thermal` were not part of the CCRTh
  cache's name, and that cache is consulted before they are read — so raising
  them changed nothing whenever a table existed, and whichever setting computed
  a configuration first served every later run of it. They are now part of the
  name on both backends; the shipped tables were renamed with their contents
  untouched, so no published number moves.
- **Both command-line tools now answer a bad setting identically.** Every
  rejected configuration exits 2 (the C tool returned 1 for everything its
  validator caught, against a convention its own source already stated), and
  all 34 rejected inputs tested produce the same message — the range and type
  templates, the explanatory notes, the number formatting, the "did you mean
  …?" suggestion and the paths in file-not-found errors were all divergent.
  An unwritable `--output_file` is now one `error:` line on the Python side
  instead of a raw traceback.
- **`primat --json` output is valid JSON again when a file is written
  alongside it.** All six "[output] … written to …" announcements — time
  evolution, final abundances, the background TSV, the decay-era table and the
  two Monte-Carlo matrices — went to stdout on both backends, ahead of the JSON
  document, so `primat --json --output_time_evolution … | jq` failed. They now
  go to stderr, where progress messages belong; stdout carries only the
  results.
- **The startup banner no longer aborts a run on a Windows console.** Its
  box-drawing characters cannot be encoded in cp1252, the Windows default, so
  `print(_banner())` raised `UnicodeEncodeError` before any physics — killing
  `runfiles/primat_run.py`, the documented validation entry point. Both
  backends fall back to an ASCII banner and separator when the console codec
  cannot carry the originals.
- **The published changelog no longer links to an unrelated commercial site.**
  A bare `README.md` in one entry was auto-linkified by MyST into
  `http://README.md` — `.md` being a live top-level domain — and rendered as an
  outbound link on the documentation site.
- **A run with a non-default `alphaem` or `me` no longer overwrites the
  shipped QED pressure tables.** Those two files keep fixed names, and a
  fingerprint mismatch used to rebuild *and rewrite* them — harmless while the
  constants were compile-time, a repository-dirtying footgun once they became
  parameters. They are now written only when `recompute_qed_corrections` asks,
  which is what that flag has always documented; a mismatched run recomputes
  the ~0.3 s of integrals in memory. Use `cache_dir` to cache a second
  configuration's pair.
- **A mistyped rate-variation key is no longer silent on the standalone C
  CLI.** `cprimat --set p_n_p__d_gg=1` ran to completion with that rate
  unvaried and no message, where the Python CLI warns (and, under
  `strict_params`, raises). Both now report the same text and the same exit
  status. Runs reaching the C solver through `primat/backend.py` were already
  covered, since `PRIMATConfig` validates the params dict first.
- **The `reference` test tier now runs the configuration it documents.** Its
  parameter set listed four of the eight solver settings
  `runfiles/primat_reference_run.py` uses, so the tier meant to reproduce the
  published validation reference was reproducing a different run — by 2.0e-08
  in `large, amax=8`'s D/H, 6.6x the ±3e-9 bound the same table advertises. It
  passed because the published constants had been produced by the tier's own
  configuration rather than the documented one. The settings are now mirrored,
  the whole "Validation reference" section is re-snapshotted against one tree,
  and the reference run's two weak-rate caches are shipped, taking that run
  from 1395 s to 23 s with bit-identical output.
- **Both backends now interpolate the CCRTh thermal correction the same way.**
  The finite-temperature n↔p correction is read from a cache the two backends
  share, but they drew different curves between its nodes: Python fitted
  scipy's global quadratic B-spline, the C backend a local 3-point Lagrange
  quadratic. That is the same class of mismatch that once dominated the
  *non-thermal* rate and was fixed there by adopting a shared cubic; the
  thermal channel was never converted, and its cache grid is ~8× coarser, so
  the two curves parted by up to 1e-05 of the rate between nodes while
  agreeing to 8e-13 at them. Both now fit a not-a-knot cubic in linear space
  (linear rather than log-log because that correction changes sign). At
  converged tolerance this removes ~94 % of the cross-backend YPBBN gap and
  ~85 % of the D/H gap; every published number stays inside its pinned
  tolerance. Found by the new cross-backend divergence harness.
- **The background's time coordinate no longer flows through a linear `T(a)`
  inverse.** Python's `t(a)` ODE took its right-hand side through
  `interp1d(a_grid, T_grid)`, a linear inversion carrying a median 3.9e-06
  error at the default node density. Sitting *inside* the RHS, it capped `t`'s
  accuracy regardless of the ODE tolerance, and the stored temperature grid was
  built by the same inverse, injecting the error twice — including into the
  grid the n↔p weak-rate tables are tabulated on. Python now integrates
  `dt/d(lnT)` directly (no `T(a)` inversion anywhere, verified against
  `scipy.quad` to 1.6e-09), both backends interpolate `t_of_T`/`T_of_t` in
  log-log, and the C samples its output arrays on the ODE's own `lnT` grid
  instead of recovering T from a log(a) grid through that same inverse.
  `t_of_T` self-convergence improves ~15× on both backends; the published
  validation reference still holds within its tolerances, so nothing was
  re-pinned. `external_scale_factor` keeps the inversion, having no closed-form
  `d(ln a)/d(ln T)`.
- **The C backend no longer runs on freed memory after an `--ini` type error.**
  A string field was released *before* its replacement was type-checked, so a
  rejected value left the field dangling; both C loaders then warned and
  continued, and `cpr_config_free` double-freed at exit. Proven under ASan with
  `network = 3`. Type errors are now fatal on the C side as they always were on
  the Python side, and a nonexistent `user_nuclear_dir` is rejected through
  every entry path rather than only the dedicated flag.
- **`--ini` parameters now reach the Monte-Carlo workers.** The MC override set
  was built from CLI flags alone, so `--ini Omegabh2=0.030 --mc 2` printed a
  single-run D/H for the requested cosmology while every MC central value and σ
  described the *default* one, with nothing to indicate the mismatch. A related
  aliasing bug is fixed with it: two string `--set` overrides shared one static
  scratch buffer, so the second silently overwrote the first.
- **The GUI's exported network zip carries the rate tables the network actually
  pins.** For a reaction the user had not customised, `export_zip` assumed a
  `<name>_primat.txt` filename instead of reading the base network file's own
  `name, filename` pairing. Downloading `small_parthenope` therefore shipped
  `small`'s tables under the Parthenope name; re-importing reproduced D/H
  −2.6 % and Li7/H +15.5 % off the run it claimed to reproduce, silently.
  `small` and `large` exports are unaffected.
- **The MT era's species set no longer depends on the network's name.** Python
  special-cased `network="small"`, so a GUI custom network adding, say,
  `d_a__Li6_g` integrated the MT era with one fewer nuclide than the C backend
  did. Both now intersect the same fixed 18-reaction MT set. Default runs are
  bit-identical. The MT BDF `atol` is aligned at 1e-15 on both backends too
  (C was 1e-16, worth 1.4e-06 in D/H on the Python side).
- **A failed `output_file` write is no longer silent.** The C CLI ignored the
  writer's return code and exited 0 having written nothing; it now exits 1 with
  a message, as the Python side raises `OSError`.
- **`run_bbn` releases the GIL** for the duration of the solve, so Ctrl-C works
  and a GUI thread is no longer blocked; `run_mc` no longer leaks a `PyFloat`
  per quantity per call.
- **A non-UTF-8 upload in the GUI reports a clear error** instead of an
  uncaught `UnicodeDecodeError` — the `.decode()` sat one line above the `try`
  meant to catch it. GUI memory is also bounded now: the heavy solve/preview
  caches were unbounded and process-global (~16 MB per configuration) and the
  export zip was rebuilt on every rerun rather than on demand, both of which
  matter on the ~1 GB public demo.
- **C backend hardening.** `cpr_load_network`'s four name matrices are
  heap-allocated (219 KB → 27 KB stack frame, against MC workers' 512 KiB
  thread stack); `rate_grid_npts` must be ≥ 2 on both backends (it was ≥ 1,
  giving C an out-of-bounds read and Python a `ZeroDivisionError`); JSON string
  escaping reserves for the `\u00xx` worst case; and per-solve scratch buffers
  that were malloc/free'd 10⁴–10⁵ times per LT solve are hoisted.
- **The C CLI's output now matches the Python CLI's exactly** — `--json` in
  both its shapes, the plain-text report (including the running-time line and
  the centred header), `--credits`, the overlay notice, and the verbose tag
  stream, all of which diff empty against Python.
- **`generate_qed_tables.py` writes where the solver reads.** Its output path
  had not existed since the data tree moved, so the documented regeneration
  command silently created a stray directory and left the shipped tables
  untouched. A missing target is now a hard error naming the cause.
- **The NUBASE half-life cross-check reads the right column**, per the file's
  own format header; the off-by-one mis-parsed 8 ground states and printed a
  spurious warning on every otherwise-clean run.
- **The docs build is green again** (`sphinx-build -W --keep-going`): one
  heading level was breaking it.
- **Cache fingerprints now include the physical constants they were computed
  from.** `primat.constants.CONST` has 26 fields and *none* of them appeared in
  any fingerprint, so editing `me`, `alphaem`, `gA`, `mn`, `mp`, `Vud`,
  `radproton`, `kappa_n/p` or `GF` silently reused cached weak rates, e±
  thermodynamics and QED tables computed with the old value — a wrong answer
  with no warning. All four fingerprints gain one `constants_hash` field,
  computed identically by both backends (`cache_utils.constants_hash` /
  `cpr_constants_hash`, verified equal: `672484f3068a1c59`). It hashes the whole
  struct rather than a per-cache list of the constants each table consumes:
  over-invalidating costs a recompute, under-invalidating is a wrong answer.
  `WEAK_RATE_FORMAT_VERSION` is bumped 4 → 5 and the shipped tables were
  re-keyed in place — renamed and re-headered, with every data row verified
  byte-identical (max column difference exactly 0). No observable moves.
  Overriding `me`/`alphaem` on a `PRIMATConfig` now raises instead of diverging
  silently, since the C backend reads its own compiled-in constants.
- **The e± thermodynamic cache carries its hash in the filename**
  (`electron_thermo_<hash>.txt`, was a fixed `electron_thermo_cache.txt` with
  the hash only in its header). Configurations coexist instead of evicting one
  another, so a full test-suite run no longer leaves the shipped, git-tracked
  copy modified, and alternating between two configurations stops paying the
  rebuild every time. `primat --cache-info` / `--cache-clear` now cover the
  `plasma/` tree as well as `weak/` on both backends, and newly generated
  files are `.gitignore`d like the weak-rate ones.
- **The QED plasma-pressure tables are fingerprinted, and take α and mₑ from the
  config.** The two `QED_pressure_correction_e{2,3}.txt` files had no
  fingerprint at all, so a table built with different constants or a different
  T grid was loaded silently; they now carry a header (`format_version`,
  `constants_hash`, `T_min`, `T_max`, `n_pts`) and are rebuilt on a mismatch.
  The Python path also read α and mₑ from `qed_pressure`'s own module-level
  copies rather than from `cfg` — the C path always passed `g_const` — so the
  two backends were bound to independent sources of truth for the same
  constant, agreeing by convention only. Both now source them from `cfg`.
  The shipped tables gained header lines only; data rows are unchanged.
- **Python's vegas is seeded deterministically**, mirroring the C backend's
  existing `th_vegas_seed`. A Python CCRTh recompute previously did not
  reproduce even itself, so the cache's promise that the hash identifies the
  contents held only up to Monte-Carlo noise. Two successive recomputes are now
  byte-identical. The shipped thermal tables are deliberately *not* regenerated:
  a seeded recompute differs from them by up to 3.8e-2 on the CCRTh term itself
  (a ~1e-3 correction to the rate), which moves D/H by 3.4e-11 — 88× inside the
  ±3e-9 regression tolerance.
- **Custom NEVO tables of non-shipped width now work on the C backend.** A
  6-column `nevo_file`, and a `nevo_spectral_file` of any width other than the
  shipped 86, passed `PRIMATConfig`'s validation and ran on the Python backend
  but aborted on the C one (`expected 7 columns, found 6`) — i.e. on the
  `force_backend="auto"` default. Both NEVO readers now take the width from the
  file. The same hard-coded 86 also left an out-of-bounds read when a
  `nevo_grid_file` declared more y-nodes than the spectral table had columns;
  that combination is now rejected with a clear message on every path,
  including the standalone `primat-c` CLI, which does not run `PRIMATConfig`.
- **e± thermodynamic tables are now accurate at the low-temperature edge of
  their grid.** Both backends integrated to an *absolute* tolerance (1e-12)
  that exceeded the integrand's own magnitude there (~e⁻³⁰ ≈ 9e-14), so the
  tabulated values carried no significant digits below T ≈ 0.024 MeV and the
  Python- and C-written tables disagreed by ~1e-4 while sharing one fingerprint
  — whichever backend last recomputed silently became the other's plasma input.
  Both now use a relative tolerance and agree to ~1e-11 across the whole grid.
  `ELECTRON_THERMO_FORMAT_VERSION` is bumped to 2, so any existing
  `electron_thermo_cache.txt` is rebuilt automatically on first use. Observable
  impact is nil: the affected region has ρ_e/ρ_γ ≤ 5e-7 and sits at or below
  the hard-zero cutoff.
- **The ΛCDM term no longer enters the background solve through a
  radiation-domination approximation.** `ρ_CDM` was evaluated with
  `a ≈ T0CMB/T_γ`, anchored at today, and swapped for the exact `a(T)` only
  *after* both background ODEs had run. That anchor is wrong in the BBN era by
  the e⁺e⁻ reheating factor cubed — 2.73× at 10 MeV, 2.69× at 1 MeV, exact only
  below ~0.01 MeV — biasing H by ~2e-7 and diverging needlessly from the C
  backend, which always had the exact `a`. `a(T)` is now published as soon as
  the entropy ODE has produced it, before anything evaluates the Hubble rate;
  the approximation and its swap step are gone. D/H moves by 6e-13 against a
  ±3e-9 regression tolerance.
- **`custom_background` runs sample their T grid logarithmically.** The grid was
  sized per *decade* but laid out linearly, so for a table spanning the usual
  40 MeV → 0.001 MeV window only 6 of its 2761 points fell below 1e9 K (1161
  when log-spaced). That grid sets the node spacing of the T_ν(T_γ) interpolant
  every n↔p rate integrand reads; the interpolation error it caused reached
  1.2e-4 on the p→n rate and is now ~280× smaller. **If you have run
  `custom_background` before and kept a `cache_dir`, delete its cached
  `nTOp_*.txt` files once** — the fingerprint cannot see a grid-layout change,
  so a stale table would otherwise be reused.
- **Two different `custom_background` tables no longer share one cached rate
  table.** In that mode the T grid comes from the supplied table's own range,
  not from `T_start_cosmo_MeV`/`T_end_MeV`, so nothing in the weak-rate
  fingerprint distinguished one custom background from another. It is now keyed
  on the table path (as `nevo_file` already was), emitted only when set — runs
  without `custom_background` hash exactly as before and keep hitting the
  shipped caches.
- **Additive sensitivity rows are now genuine elasticities** (`d ln O / d ln p`,
  the meaning of every other row and of the table as a whole). An additive
  target such as `DeltaNeff` was divided by the *multiplicative* rows'
  `2 ln(1+rel_step)`, so its cell scaled with a `rel_step` that plays no part in
  an absolute ±step perturbation — inflating it ~100x at the documented
  `step=1.0`, `rel_step=0.01` setting. `SensTarget` gains `ref`: the fiducial of
  the physical parameter the offset displaces (`ref="Neff"` reads the run's own
  central value), making the row `∂ln O/∂ln N_eff`. `Y_P` vs `N_eff` now reads
  `+0.1647` (was `+5.4994`), D/H `+0.4104` (was `+13.6513`) — the textbook
  `Y_P ∝ N_eff^0.16`, `D/H ∝ N_eff^0.41`. Without `ref` the fallback is the
  linear `2*step`, documented as per-unit and not comparable with other rows.
  `notebooks/Sensitivity.ipynb` and the *How-to → Sensitivity tables* page were
  updated and the notebook re-executed.
- **`primat --cache-info`/`--cache-clear` honour `--data_dir` and
  `--set cache_dir=…`** on the Python CLI, which previously always inspected and
  cleared the *default* data tree whatever the user pointed at (the C CLI
  already honoured `--data_dir` for both commands). `--cache-clear` still clears
  every cached file; both help texts now name the regeneration cost.
- **A test no longer overwrites a shipped weak-rate cache.**
  `test_recomputed_rates_match_cached` set `weak_rate_cache=False` but left
  `save_nTOp` at its default `True`, so it wrote its freshly integrated rates
  over the git-tracked `nTOp_<hash>.txt` for that fingerprint. The recompute
  differs by only ~1e-10, but the default `numerical_precision=1e-7` leaves ~1e-6
  of adaptive-step jitter in the observables, so merely running the test suite
  shifted the repo's validation reference (YP 0.24700086 → 0.24700060) while
  every pinned tolerance still passed.
- **Monte-Carlo `prev` reuse no longer corrupts a shrinking request.** When
  `num_mc < len(prev)`, only the requested quantities' samples were truncated
  while every nuclide kept all `len(prev)` rows, so per-nuclide `mean`/`std`
  were silently computed over the wrong sample count and
  `MCResult.samples_array()`/`cov()`/`corr()` (hence `dump_mc_samples`) raised
  `ValueError`. The C backend was already correct; this was a Python-only
  parity bug.
- **Monte-Carlo now varies reactions added through `custom_network["added"]`.**
  The varied set was derived from the network *file*, which cannot list a
  brand-new reaction, so a GUI-added reaction was integrated but never sampled —
  while the C backend, which iterates the solved network, did sample it.
- Numpy scalars (`np.int64`, `np.float32`, `np.bool_`) are accepted as parameter
  values on both backends; previously they aborted the Python backend with an
  opaque `TypeError: Object of type int64 is not JSON serializable` from the
  weak-rate cache fingerprint, and the C backend with a type error from the
  extension wrapper. Hash-preserving, so no cache file is invalidated.
- The C backend prints `[init-c] Initialisation complete in X s` and quotes the
  `network` name in its options recap, restoring line-for-line diffability of
  the two backends' verbose headers.
- `cache_utils.write_cache_with_fingerprint` can write a cache path with no
  directory component (`os.makedirs("")` used to raise, turning a perfectly
  writable target into a "could not write cache" warning — including for the
  form shown in its own docstring example).
- **QED plasma-pressure tables are now read with the same interpolant they are
  written with** (both backends). `Plasma._load_tables` documents three
  interchangeable modes (load from file / compute analytically / recompute),
  but the file path built a *linear* `interp1d` while the analytic path built a
  cubic spline. On the shipped 500-point log grid over [1e-3, 1e2] MeV
  (`d(lnT) = 0.023`, `δP ~ T⁴`) that cost ~8e-4 relative on δP, so — since
  `δP/ρ_pl ~ 4e-4` during BBN — Neff shifted in its 6th decimal according only
  to whether the cache files happened to exist. Both paths now go through
  `plasma._qed_spline` (and `cpr_cubic_spline_fit_notaknot` in C). Residual
  file-vs-analytic disagreement is 3.9e-06, set by the tables' own `%.6E` write
  precision rather than by interpolation, and pinned by
  `tests/test_qed_pressure.py::test_file_and_analytic_paths_agree`.
- **The NEVO heating function is clamped to N ≥ 0** (both backends). Heating is
  entropy flowing from the EM plasma into the neutrinos and cannot reverse, but
  74 of the 600 rows of the shipped `NEVOPRIMAT_col_1_7.csv` carry a negative
  residual from the NEVO solve itself (all within T_γ ∈ [0.0315, 0.0835] MeV,
  reaching −4.2e-06 against a peak N of 4.2e-03), which the a(T_γ) ODE
  integrated as a spurious reverse transfer.
- Together these two shift the default run by ~4e-6 relative in D/H and Li7/H
  (`D/H` 2.4358955e-05 → 2.4359049e-05 on the Python backend) — two orders of
  magnitude inside the ±3e-9 D/H regression tolerance, and cross-backend
  agreement stays at ~5e-6 against the documented 5e-5 bound.
- **`wnEDE ≤ 1/3` is now rejected** when `fEDE > 0`, on both backends. The EDE
  peak scale factor solves `u^(3(1+wnEDE)) = 4/(3·wnEDE − 1)`, which has no root
  for `wnEDE ≤ 1/3`: such a component dilutes no faster than radiation, so its
  energy *fraction* never peaks during radiation domination and `fEDE` (defined
  at that peak) is meaningless. Previously `wnEDE = 1/3` raised a bare
  `ZeroDivisionError` and `wnEDE = 0` silently produced a *complex* scale factor
  that surfaced hundreds of lines later as solve_ivp's "`y0` is complex", while
  the C backend produced a NaN background without complaint. Both are standard
  axion-like values (`wn = (n−1)/(n+1)` for n = 2 and n = 1).
- `AnalyticDistortion` no longer omits `x_of_Tg`, which made a documented
  `NeutrinoHistory` protocol attribute raise `AttributeError` instead of
  returning `None`. Latent only — `PRIMATConfig` currently forbids the flag
  combination that would reach it.
- Custom NEVO tables: overriding `nevo_spectral_file` *without* also overriding
  `nevo_grid_file` is now validated against the shipped `NEVOGrid.csv`, instead
  of being computed and then never compared — a width mismatch used to surface
  as a shape error deep inside `RegularGridInterpolator`.
- `Constants.erg` was missing a square on `second` (it read `gram·cm²/second`
  against its own docstring). Numerically inert under the natural-units
  convention, where all three base units are 1.
- `primat/qed_pressure.py`'s local `_ME_MEV` was the CODATA 2014 electron mass
  while `CONST.me` is CODATA 2018 — and the C backend already used the latter,
  so the two backends generated QED tables at different electron masses. Now
  identical (an 8e-9 relative change to freshly computed tables).
- **n↔p weak-rate cache keys** (`WEAK_RATE_FORMAT_VERSION` 1 → 4, both
  backends). Three configuration fields changed the rates but were absent from
  the fingerprint, so runs that differed only in one of them silently shared a
  cache file:
  - `munuOverTnu`/`munuOverTnu_e` was missing from the **thermal** (CCRTh)
    fingerprint, although the thermal integrands carry an explicit
    `exp(−sgnq·ξ_ν)` neutrino occupation. Degenerate-BBN runs were reusing the
    ξ=0 table — worth ~4e-3 of the base rate at ξ_e = 0.3, T = 1e10 K, i.e. far
    above anything YP tolerates — and, on a cold cache, writing their own
    ξ-specific numbers under the filename standard runs load.
  - `nevo_grid_file` was missing from the weak-rate fingerprint while its
    partner `nevo_spectral_file` was present; the two jointly define the
    tabulated distortion the SD term integrates.
  - `sampling_temperature_per_decade` was missing from the weak-rate
    fingerprint. It sets the node spacing of the linear T_ν(T_γ) interpolant
    every rate integrand reads: coarsening it moves the rates by ~1e-3 (40
    points/decade) down to ~1e-5 at the default 600.

  The version constant is also bumped past the v2/v3 generations that were
  documented in the changelog comment but never actually applied, so pre-v3
  cache files (whose `nTOp_*.txt` still included CCRTh, and whose thermal table
  was unclamped below 10^8.2 K) can no longer be loaded. The shipped tables were
  re-keyed in place — same numbers, new hash-named filenames — so default runs
  still hit the cache. **Existing `cache_dir` trees and editable installs:**
  stale files are simply never loaded again (delete them, or `primat
  --cache-clear`), but a compiled C extension built before this change computes
  the old hashes and will miss the re-keyed tables — rebuild it
  (`python setup.py build_ext --inplace`).
- No observable changes: with a rebuilt extension both backends reproduce the
  previous D/H, YP and Neff bit-for-bit.

### Removed
- **Four unused names dropped from the public surface.** `cfg.erg` and
  `cfg.s0CMB` (with their `primat.constants.Constants` properties and their
  `cpr_erg`/`cpr_s0CMB` counterparts) took part in no formula on either
  backend; `s0CMB` was additionally carried in `DERIVED_OVERRIDABLE`, so the
  per-config constants machinery was maintaining a quantity nothing read.
  `cpr_cubic_spline_fit_natural` had no caller among the C backend's twenty
  spline fits and no Python counterpart at all. Nothing else moved: the
  present-day photon *number* density `n0CMB` is unaffected, no cache is
  re-keyed (cache fingerprints declare `Constants` fields, never derived
  properties), and every observable is bit-identical on both backends.

### Changed
- **`numpy` and `scipy` now carry tested lower bounds** (`numpy>=1.26`,
  `scipy>=1.11`), where the dependency list previously declared none at either
  end. Continuous integration gained two lanes to keep them honest: one that
  pins exactly those versions on the oldest supported Python, and one that
  runs against the newest numpy and scipy including pre-releases, so an
  upstream change is seen here rather than in an install.
- **The pure-Python backend is ~1.6× faster, with every printed digit
  unchanged.** A warm `small` run goes from 0.750 s to 0.460 s and
  `large, amax=8` from 1.085 s to 0.684 s (C: 0.040 s / 0.134 s). Three
  changes, none of them physics: the background's `T_of_t`/`t_of_T`/`a_of_t`
  lookups now evaluate the same two-node linear formula scipy does without
  paying `interp1d`'s ~10 µs of per-query input validation; the MT/LT BDF
  solves call LAPACK's `getrf`/`getrs` directly instead of through
  `scipy.linalg.lu_factor`/`lu_solve`, whose batch-dispatch and finiteness
  wrappers cost more than the factorisation itself at ~15k calls per solve;
  and the n↔p rate floor takes a scalar branch instead of building a 0-d array
  per query. Both fast paths verify at build time that they reproduce what
  they replace, and fall back to scipy if they cannot. Observables are
  identical float-for-float on `small`, `small_parthenope` and
  `large, amax=8`, so the cross-backend gap is unmoved.
- **Each cache is now keyed on the constants it actually reads,** rather than
  on all 26 at once: eight for the n↔p rate table, five for the CCRTh thermal
  correction, two for the QED pressure tables, one for the electron
  thermodynamics (`cache_utils.CACHE_CONSTANTS`, mirrored in
  `primat-c/src/cache.c`). Eight of the sixteen measured constants — `T0CMB`,
  `GF`, `mZ`, `Vud`, `ma`, `He4Overma`, `HOverma`, `Neff_SM` — change no cached
  number, yet used to re-key every cache: `--T0CMB 2.7250` cost a 116 s
  Monte-Carlo rebuild of a bit-identical CCRTh table. The declared lists are
  proven, not asserted: `tests/test_cache_constant_deps.py` perturbs every
  settable constant, rebuilds each cache from scratch, and fails if the data
  moves when the constant is undeclared (under-keyed, silently wrong physics)
  or stays put when it is declared (over-keyed, a needless recompute). Cache
  format versions bumped accordingly (weak/thermal 5→6, electron thermo 2→3,
  QED 1→2) and the shipped files re-keyed with byte-identical data rows.

### Added
- **The 16 measured physical constants are now ordinary parameters.**
  `alphaem`, `GF`, `mZ`, `me`, `mn`, `mp`, `T0CMB`, `gA`, `Vud`, `kappa_p`,
  `kappa_n`, `radproton`, `ma`, `He4Overma`, `HOverma` and `Neff_SM` are
  `DEFAULT_PARAMS` keys, settable through a params dict, `--me 0.511`, a
  `me = 0.511` INI line, `cpr_config_set_by_name`, and a "Constants" group in
  the GUI. The dividing line is whether a number has an error bar: the other
  ten (`Kelvin`/`second`/`cm`/`gram`, `kB`/`clight`/`hbar`/`MeV`/`keV`, `Mpc`)
  are exact by the natural-units convention, the 2019 SI redefinition or an
  IAU definition, and stay frozen — overriding one is now rejected with an
  error rather than silently honoured by Python and ignored by the C backend.
  Every derived quantity follows an override (`sW2`, `deltakappa`, `mB`,
  `n0CMB` and the `eta0b` chain), and each cache is re-keyed on the constants
  it was computed with, so a run with a shifted `me` cannot load the
  default-`me` tables. At default values both backends are bit-identical to
  before and every shipped cache filename is unchanged. One documented
  exception: the Pitrou & Pospelov QED correction to the radiative-capture
  nuclear rates is a fit performed at the CODATA α and keeps its own literal,
  so `alphaem` does not reach it.
- **`--list-reactions` on both CLIs.** The `p_<reaction>`/`delta_<reaction>`
  rate-variation keys are a per-network, unbounded family, so they cannot
  appear in `--list-params`; a user previously had to already know a
  reaction's exact internal name to vary its rate. The new flag prints the
  reaction names of the selected `--network`/`--amax`, byte-identically on
  both backends.
- **`--list-params` and `--mc-jobs` on the C CLI.** The first enumerates the
  same field table the setter dispatches on, so the listing cannot drift from
  what is settable; the second exposes the MC job count, which was hard-coded
  to all cores.
- **Test coverage for what was documented but unpinned**: all 21 cells of the
  per-nuclide reference table (live and static); the evolution TSV header
  compared byte-for-byte across backends; `thermal_corrections` (CCRTh) given a
  physical-effect test, so every n↔p correction flag now has one; the free
  neutron added to the large-vs-`amax=8` comparison; and the first end-to-end
  coverage of `generate_rates/`, including a byte-for-byte regeneration check
  of all 395 shipped artifacts.
- **`tests/test_cache_parity.py`** — cross-backend *cache* parity, the companion
  to `test_backend_parity.py`. The two backends share every on-disk cache; this
  module is what makes that safe rather than merely convenient. It asserts that
  both emit the same cache filenames (pinning every fingerprint implementation
  field-for-field) and that their columns agree at documented, measured
  tolerances. Verified to detect what it claims: perturbing the C electron mass
  by one ulp fails 7 of its 8 tests.

### Documented
- **The two backends' known, deliberate divergences** are recorded in
  `tests/README.md`'s "Known cross-backend divergences", with `README.md`'s
  "Backend parity contract" stating what parity means and which tests enforce
  it. Two divergences are intentional — the HT-era integrator (`LSODA` vs
  RK45; aligning both on BDF was tried and *degraded* YP parity) and
  `external_scale_factor`'s interpolant (a C-side performance workaround that
  measures worse when mirrored) — and the residual D/H gap is open, structural
  rather than round-off, and downstream of the background.
- **`docs/` corrected against live runs**: the landing page's quick-start
  numbers (stale by ~2900× the regression pin), the MC how-to's worked example
  and its correlation prose (both fabricated), the evolution-TSV schema in the
  output how-to (which described columns that do not exist and omitted the
  background TSV entirely), and the custom-network export contract, which is
  verbatim/original-grid precisely so a round trip is bit-for-bit.
- **`docs/performance.md`** gained a profile of where the pure-Python backend's
  time goes, the optimisations already applied, and the conclusion that the
  remaining floor is scipy's pure-Python BDF stepper.
- **`Omeganurel` and `OneOverOmeganunr` are per neutrino flavour** (ν + ν̄), not
  summed over the three — now stated in `Background.Omeganuh2_relnu`/`_nrnu`,
  `docs/howto/output.md` and `primat-c/include/background.h`. The values are
  unchanged (multiply `Omeganurel` by 3 for the usual quoted total ≈ 17); the
  per-flavour convention is the natural one for `OneOverOmeganunr`, whose ≈ 93
  reproduces the standard Σm_ν / 93.1 eV normalisation.
- **The ΔNeff extra species deliberately uses a different "instantaneous
  decoupling" normalisation from the SM neutrinos** when `QED_corrections=True`:
  `T_nu_decoupling`'s free-gas σ_∞ rather than the QED-corrected `_sbar_ref`,
  leaving it ~0.31 % low in energy density during BBN. That is what makes the
  reported `Neff` come out as `Neff_SM + ΔNeff` to machine precision, i.e. makes
  the knob mean what its name says; `rho_nu_extra` now spells out the trade-off
  and the measured cost of the alternative, and `T_nu_decoupling` no longer
  claims to be the SM neutrino temperature in that mode.
- `plasma.rho_SM`/`p_SM` are labelled as the ξ = 0, no-spectral-distortion
  reference quantities they are — **not** the Friedmann source, which is
  `StandardBackground.Hubble` (it adds each flavour's own ξ and `rho_nu_SD`).
- `Background.t_of_T` documents its valid range: outside
  `[T_end, T_start_cosmo]` it extrapolates linearly and can return a negative
  time. This is distinct from the radiation-domination extrapolation below the
  NEVO table's edge, which is inside the solved span and correct in both
  `external_scale_factor` modes.
- The `external_scale_factor` True/False agreement is quoted at its measured
  ~1e-5 (per-observable figures given) instead of "~1e-6".
- Every constant in `primat/constants.py` now records its edition (SI 2019
  exact / CODATA 2018 / CODATA 2010 / PDG 2020 / PDG 2018 / AME2020 / AME2016 /
  Fixsen 2009); the set is deliberately not single-vintage, and `gA`/`Vud` are
  flagged as the group that feeds the n↔p rates.
- Assorted `plasma.py` docstring drift corrected: a reference to a
  `_setup_qed_pressure` method that does not exist, a recompute-mode paragraph
  naming the wrong output file, stale `rates/plasma/` paths, and a "no
  module-level mutable state" claim that overlooked the four numba integrand
  handles (harmless, but now described accurately).
- `primat/weak_rates/corrections.py` now cites the Phys. Rep. **equation**
  numbers for every correction term, instead of section numbers alone: the
  relativistic Fermi function (Eq. 100), the resummed radiative factor
  (Eq. B35, with g = B32, Spence L = B33, constants B31/B36), the finite-mass
  terms (Eqs. 114, 115a/115b, χ_FM from App. B.3) and the four CCRTh
  sub-integrands (Eqs. 107, 109, 112a, 113, with the F_± kernels of Eq. B41 and
  the B51a/B51b kernel). Three section references were also off against
  `biblio/Pitrou_etal_PhysReptArxivVersion.pdf` and are corrected: the T=0
  radiative corrections are §III.E (was §III.D) and the finite-temperature ones
  §III.F (was §III.H, which is "Weak magnetism").
- The `F_+`/`F_−` asymmetry in the bremsstrahlung soft subtraction is now
  documented as deliberate on both backends, citing Phys. Rep. Eq. B43 where it
  is printed, plus the measured consequence of "correcting" it (the CCRTh sum
  would grow to ~0.8% of the base rate at 3e10 K).
- Two accepted-but-not-self-consistent flag combinations are called out in
  `config.py`: `thermal_corrections=True` with `radiative_corrections=False`,
  and the absence of the SD-FM term outside analytic-distortion mode.
- `background=`: documented that the weak-rate cache is keyed on the config
  alone and cannot see a custom background's temperature grid — use
  `weak_rate_cache=False`/`save_nTOp=False` for a non-standard history.

### Changed
- **Test tolerances tightened, none loosened.** The only guard on the shipped
  n↔p weak-rate cache against a fresh integration asserted `rel=2e-3` where the
  measured agreement is 2.3e-08, so a 0.2 % error in any correction term passed
  while moving YP by ~2e-4; it is now `1e-6`. The baryon-conservation bound
  went 1e-6 → 1e-10 (measured 1.6e-12). The documented "run `primat_run.py`
  after any change" workflow gained a second, correct tolerance column — the
  ±3e-9 D/H bound belongs to the high-precision reference run and *fails* on a
  healthy tree when applied to the default-precision script — and the script's
  printed numbers are now checked by a test rather than by hand.
- Python backend: the ten weak-rate Fermi-Dirac integrand kernels
  (`weak_rates/integrands.py`) and the four e± electron-thermo integrands
  (`plasma.py`) are now numba-compiled with `cache=True`, so a fresh process
  (joblib MC worker, Streamlit server, re-run CLI) loads the compiled machine
  code from numba's on-disk cache instead of recompiling — ~2.3 s of cold-start
  JIT saved per process (measured). The plasma integrands were moved to module
  level and now take the electron mass `me` as an explicit argument (rather than
  closing over `cfg.me`), which is what makes the on-disk cache safe. No effect
  on any observable.
- C backend (`primat-c`) now routes its ~250 unrecoverable heap allocations
  (ODE work vectors, spline tables, loaded network, result arrays, …) through
  new checked helpers `cpr_xmalloc`/`cpr_xcalloc`/`cpr_xrealloc`
  (`include/xalloc.h`): a failed allocation now prints
  `primat: out of memory (<bytes>) at <file>:<line>` and exits, instead of
  dereferencing NULL and crashing anonymously. Sites that intentionally
  degrade gracefully (cache writers with their own NULL checks, `errmsg`-return
  paths) are unchanged. No effect on any observable — purely OOM-diagnostic
  robustness.

### Removed
- `rate_interp_order` `DEFAULT_PARAMS`/C config parameter: never consumed by
  any solver, resampler, or rate-lookup path on either backend (rate-table
  resampling always hardcoded log-log cubic, and the per-step master-grid
  lookup always used linear `searchsorted` regardless of its value), so
  setting it to `quadratic`/`cubic` changed neither observables nor runtime.
  Removed together with its C field/default/validation, both param
  templates' entries, both backends' enum tests, and the now-empty
  `_PARAM_CHOICES` machinery it was the last user of.

## [0.3.2] - 2026-07-11

### Added
- `output_rates_time_evolution` now writes per-reaction forward-rate columns
  (`<reaction>_frwrd`, small/small_parthenope networks), on both backends —
  previously a no-op. Populated in `EvolutionResult.rates` and round-tripped
  by `primat.evolution.load_evolution`.
- `primat.sensitivity.sensitivity_table` — a one-call API returning the
  logarithmic-sensitivity matrix ∂ln(observable)/∂ln(parameter) as a
  `SensitivityTable` dataclass (`.to_markdown()`/`.to_dataframe()` views),
  with `notebooks/Sensitivity.ipynb` reduced to a thin demo and a new
  *How-to → Sensitivity tables* docs page.
- `primat --list-params` and `--version` (with backend build status) for CLI
  discoverability.
- Type hints across the public API, `py.typed` marker, and a lenient mypy CI
  job.
- `primat.__citation__` (ready-made BibTeX string) and top-level
  `run_bbn`/`run_mc`/`HAS_C_BACKEND` aliases, so common entry points no
  longer require reaching into `primat.backend`.
- `CITATION.cff` and Zenodo archival instructions in `PyPiGuide.md`.
- This changelog.
- Config validation with clearer error messages, shared between both
  backends.
- MC covariance/correlation output (`MCResult.cov()`/`.corr()`) and a
  `runfiles/primat_mc.py` demo script.
- Sphinx documentation site (`docs/`, published to Read the Docs), migrating
  and superseding the old `README.md`/`EXTENDING.md` prose content.
- `primat-gui`'s Final abundances tab gained a single reproduction-bundle
  download (`.py`/`.ini` + README, `primat/gui/export_params.py`): prints the
  full standard-ratio `run_bbn` centrals plus a `run_mc(seed=0)` std-only
  block, pins `force_backend` to whichever backend actually ran, and for a
  custom network embeds the exact `custom_network` dict (Python) or a
  `nuclear/` overlay directory (`.ini`, via `user_nuclear_dir`) — including
  any uploaded/edited rate tables — so a downloaded bundle reproduces the
  GUI run bit-exactly.
- `notebooks/ReactionRates.ipynb` (⟨σv⟩(T9) of any reaction vs. the Hubble
  rate, with the master-grid reinterpolation overlaid) and
  `notebooks/AnimatedAbundances.ipynb` (animated GIFs of the small-network
  abundance evolution vs. ΔNeff/Ω_b h²), both wired into the docs tutorial
  gallery with a guard test against future gallery drift.

### Changed
- The C backend now supports `extra_rho` and `decay_era`, closing two of the
  three former Python-only feature gaps. `extra_rho` callables are
  sampled onto a dense temperature grid and splined into the C Friedmann
  equation; `decay_era`'s long-lived-isotope Decay-Time propagation is ported
  via a scaling-and-squaring Padé matrix exponential, writing an identical
  `output_decay_evolution` TSV. Only `background=` (a custom `Background`
  object) remains Python-only.
- Default `Omegabh2` changed to the Planck 2018 + BAO value (`0.02242`),
  from the previous default.
- `plotly` and `joblib` moved from hard dependencies to the `mc`/`plots`/
  `gui` extras — a plain `pip install primat` install is now lighter.
- Nuclear rate tables (`large` network) regenerated at 1000 points instead
  of being reinterpolated onto the master T9 grid from a coarser source
  grid, for both backends.
- Cache trees consolidated under `cache_plasma_weak/{weak,plasma}/` with an
  additive `cache_dir` overlay redirect for read-only installs, and
  non-fatal (warn, not crash) cache-write failures on both backends
  (80 → 81 `DEFAULT_PARAMS` keys).
- Per-flavour neutrino degeneracies `munuOverTnu_e/mu/tau`.

### Fixed
- Assorted release-blocker fixes ahead of the first PyPI publish: packaging
  metadata, wheel build matrix, and related polish.
- The Python backend's `run_bbn` now also exposes a `Y_final` sub-dict,
  matching the C backend and restoring result-dict parity.
- Windows editable-install C-extension shadowing, and assorted
  Windows-portability failures in the CI Tests matrix.
- MSVC POSIX-header/pthreads build failures on the Windows leg of
  `wheels.yml`, so Windows users now get binary wheels (with the fast C
  backend) from PyPI instead of falling back to a source build (64-bit only).
- `GN` (Newton's constant) default corrected to the exact CODATA literal
  `6.6743e-11`; the previous default (`6.674299257609439e-11`) was off at
  the ~1.1e-7 relative level.

### Documented
- Electron-thermo / QED-pressure table extrapolation behaviour beyond the
  tabulated range, with a runtime warning when a run goes past it.

## [0.3.1] - 2026-07-02

### Added
- `show_progress` config flag to control `[primat]`/`[MC]` stderr progress
  messages, wired through both the Python and C backends' CLIs
  (including `--flag`/`--no-flag` boolean parsing).
- `mc_rate_rescale_cap` parameter to cap Monte Carlo rate-rescaling factors
  (default lowered from `1e3` to `30` after further validation).
- Flat `sigma_<name>` fields in MC results on both backends.
- Ctrl-C abort support for a running `primat-c` Monte Carlo (`run_mc`)
  sample.
- Memory-leak checking (`make leak-test`) and an ASan/UBSan CI job for the
  C backend.
- A pre-computed `PRIMAT_Yp_DH_ErrorMC_1000_2026.dat` table for CLASS/CAMB
  consumption.

### Changed
- Renamed `rates_dir`/`user_rates_dir` config fields to `data_dir`/
  `user_nuclear_dir` (clearer overlay semantics — see
  `docs/howto/data-overlays.md`).
- Split the three separate QED plasma-pressure correction table files into
  one consolidated `QED_tables.txt`.
- Integrated three background-ODE performance branches: dense-output RK45,
  a combined 2D background ODE, and monotone spline lookups for rate
  interpolation.
- Smoothed MC progress reporting (was jumping straight from 0% to 100%).
- Moved C backend headers from `include/cprimat/` to `include/`.

### Fixed
- `--no-show_progress` being silently ignored by the `primat-c` CLI's
  `--mc` path.
- `delta_<rxn>` rate perturbations not applying when
  `rescale_nuclear_rates` was left at its default.
- A NaN issue in the C backend's `electron_thermo` plasma cache.
- Several gcc-14 warnings (`-Wformat-truncation`, `-Wmaybe-uninitialized`,
  unused-parameter) in the C backend.

## [0.3.0] - 2026-06-25

Initial PyPI-track release. By this point the project already had its
current two-backend architecture in place:

### Added
- Dual backend: a pure-Python implementation (`primat/`) and a fast C99
  port (`primat-c/`), exposed to Python via a compiled extension and
  dispatched through `primat.backend.run_bbn(force_backend={"auto","c","python"})`.
- `primat-gui`, a Streamlit application (four usage modes: Python API, CLI,
  GUI, notebooks), including a "Customise Reactions" flow to build and
  import/export custom nuclear networks.
- Monte Carlo uncertainty propagation (`run_mc`) with rate-key resolution
  and incremental sample reuse (`prev`) on the Python backend.
- Non-instantaneous neutrino decoupling via the NEVO tables, with
  overridable `nevo_file`/`nevo_spectral_file`/`nevo_grid_file`/
  `nevo_file_prefix` parameters.
- Analytic QED plasma-pressure corrections and per-reaction weak-rate
  correction terms (Born, CCR, finite-mass, thermal, spectral-distortion).
- Grid-agnostic nuclear rate loading (tables resampled onto a configurable
  master T9 grid at load time) and the `large`/`amax` network-filtering
  mechanism.
- Unified log-log cubic not-a-knot interpolation for n↔p weak rates on both
  backends, collapsing the former C-vs-Python D/H gap to ≲1e-5.
- Unified time-evolution TSV schema (`primat/evolution.py`,
  `EvolutionResult`/`load_evolution`), implemented identically by both
  backends.

[Unreleased]: https://github.com/CyrilPitrou/primat/compare/v0.3.2...HEAD
[0.3.2]: https://github.com/CyrilPitrou/primat/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/CyrilPitrou/primat/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/CyrilPitrou/primat/releases/tag/v0.3.0
