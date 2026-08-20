#!/usr/bin/env python3
"""Drive the fuzz targets and collect crashes.

The engine runs in-process, so a sanitizer report ends the process. This
script restarts the target from its persisted corpus, saves the input the
engine had recorded in ``artifacts/.current`` as a crash artifact, and moves
on -- the same job libFuzzer's ``-artifact_prefix`` does.

Usage:
    python3 tests/fuzz/run_fuzz.py --runs 100000 [--target fuzz_ini]
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent                      # primat-c/
BUILD = ROOT / "build" / "fuzz"
DATA_DIR = ROOT.parent / "primat" / "data"


def run_target(name: str, runs: int, timeout: int, seed: int,
               stop_at_first: bool = False) -> dict:
    binary = BUILD / name
    corpus = HERE / "corpus" / name.replace("fuzz_", "")
    artifacts = HERE / "artifacts" / name.replace("fuzz_", "")
    corpus.mkdir(parents=True, exist_ok=True)
    artifacts.mkdir(parents=True, exist_ok=True)
    current = artifacts / ".current"

    env = dict(os.environ)
    env["CPRIMAT_FUZZ_DATA_DIR"] = str(DATA_DIR)
    env["ASAN_OPTIONS"] = "abort_on_error=0:detect_stack_use_after_return=1"
    env["UBSAN_OPTIONS"] = "print_stacktrace=1:halt_on_error=1"

    done, crashes, restarts = 0, [], 0
    last_line = ""
    t0 = time.time()
    while done < runs:
        batch = runs - done
        cmd = [str(binary), str(corpus), f"-runs={batch}",
               f"-artifacts={artifacts}", f"-timeout={timeout}",
               f"-seed={seed + restarts}"]
        # errors="replace": the loaders echo the offending key/value back in
        # their diagnostics, so a fuzzed byte string comes out on stderr and
        # strict UTF-8 decoding would fail on the child's own report.
        proc = subprocess.run(cmd, env=env, cwd=ROOT, capture_output=True,
                              text=True, errors="replace")
        for line in proc.stdout.splitlines():
            if line.startswith("[fuzz]"):
                last_line = line
        if proc.returncode == 0:
            done = runs
            break
        # Abnormal exit: the input in .current is the reproducer.
        restarts += 1
        blob = current.read_bytes() if current.exists() else b""
        digest = hashlib.sha1(blob).hexdigest()[:12]
        kind = "timeout" if proc.returncode == 88 else "crash"
        art = artifacts / f"{kind}-{digest}"
        art.write_bytes(blob)
        (artifacts / f"{kind}-{digest}.log").write_text(proc.stderr[-20000:])
        crashes.append(str(art))
        # Drop it from the corpus so the restart does not replay it forever.
        for f in corpus.iterdir():
            if f.is_file() and f.read_bytes() == blob:
                f.unlink()
        current.unlink(missing_ok=True)
        done += 1
        if stop_at_first:
            break
        if restarts > 200:
            print(f"  {name}: 200 distinct failures, stopping early", file=sys.stderr)
            break
    return {"name": name, "summary": last_line, "crashes": crashes,
            "restarts": restarts, "seconds": time.time() - t0}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=int, default=20000)
    ap.add_argument("--timeout", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--target", action="append", default=None)
    ap.add_argument("--selftest", action="store_true",
                    help="run only fuzz_selftest and require it to crash")
    args = ap.parse_args()

    if args.selftest:
        res = run_target("fuzz_selftest", args.runs, args.timeout, args.seed,
                         stop_at_first=True)
        print(res["summary"] or "[fuzz] fuzz_selftest: no summary")
        if not res["crashes"]:
            print("SELFTEST FAILED: the harness reported no crash for a target "
                  "that crashes on 'BOOM'", file=sys.stderr)
            return 1
        blob = Path(res["crashes"][0]).read_bytes()
        print(f"selftest ok: crash reproduced by {blob!r}")
        return 0

    names = args.target or sorted(
        p.stem for p in HERE.glob("fuzz_*.c")
        if p.stem not in ("fuzz_engine", "fuzz_selftest"))
    rc = 0
    for name in names:
        if not (BUILD / name).exists():
            print(f"{name}: not built (make fuzz-build)", file=sys.stderr)
            rc = 1
            continue
        res = run_target(name, args.runs, args.timeout, args.seed)
        print(res["summary"] or f"[fuzz] {name}: no summary")
        if res["crashes"]:
            rc = 1
            print(f"  !! {len(res['crashes'])} failure(s):")
            for c in res["crashes"]:
                print(f"     {c}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
