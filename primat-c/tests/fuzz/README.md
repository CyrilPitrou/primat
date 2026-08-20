# Fuzzing the C parsers

Every C function that reads a user-supplied file has a fuzz target here. The
point is the class of failure nobody thinks to write a test for: pass 21 found
that an empty or comments-only rate table segfaulted the C backend by trying
inputs one at a time, and a fuzzer finds that whole class without being told
what to look for.

## Running

```bash
make fuzz-build              # build the targets (ASan + UBSan + edge coverage)
make fuzz RUNS=200000        # run every target to RUNS executions
python3 tests/fuzz/run_fuzz.py --target fuzz_ini --runs 50000   # one target
make fuzz-coverage           # what each target actually reached
```

A failure leaves the reproducing input and the sanitizer report in
`artifacts/<kind>/` (untracked). Minimise it, fix the defect, and move the
minimised input to `regressions/<kind>/`, where `tests/unit/test_fuzz_regressions.c`
replays it — a fuzz finding with no committed reproducer is not closed.

## Why an engine lives here

Apple's Command Line Tools clang instruments with `-fsanitize-coverage` but
ships no libFuzzer runtime (`libclang_rt.fuzzer_osx.a` is absent), so
`fuzz_engine.c` supplies the loop: an edge map fed by
`__sanitizer_cov_trace_pc_guard`, a corpus that grows whenever an input reaches
a new edge, and a mutator carrying a dictionary of the tokens these
line-oriented numeric formats key on. Targets keep libFuzzer's
`LLVMFuzzerTestOneInput` signature, so they build against the real runtime
unchanged wherever it exists.

Crashes are not caught in-process — the sanitizers end the process — so the
engine records each input in `artifacts/.current` before running it and
`run_fuzz.py` turns that into a saved artifact after an abnormal exit.

## Proving the harness can fail

"No crashes in 200,000 executions" means nothing from a pipeline that cannot
report one. `fuzz_selftest.c` crashes on the input `BOOM`:

```bash
python3 tests/fuzz/run_fuzz.py --selftest --runs 300000
```

exits non-zero unless an artifact reproducing the fault comes back. It has
already earned its keep: the first version of that target allocated, copied and
freed a buffer whose result was unused, which clang deletes outright at `-O1`,
so it never crashed at all.

## Layout

| Path | Contents |
|---|---|
| `fuzz_engine.c` | the coverage-guided loop, corpus and mutator |
| `fuzz_<kind>.c` | one target per parser |
| `corpus/<kind>/` | seeds (shipped data + the degenerate shapes passes 20/21 collected) and the inputs the campaign evolved |
| `regressions/<kind>/` | minimised reproducers for defects found |
| `artifacts/<kind>/` | crash inputs and sanitizer logs from a local run (untracked) |
| `seed_corpus.py` | regenerates the seeds from the shipped data tree |
| `run_fuzz.py` | drives the targets, saves artifacts, `--selftest` |
| `fuzz_coverage.py` | line/region coverage each target reached |
