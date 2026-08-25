# -*- coding: utf-8 -*-
"""
benchmark.py
============
Regenerates the timing table published in ``docs/performance.md``.
Measures wall-clock time for:

- a single BBN solve, small network, both backends
- a single BBN solve, large network restricted to A <= 8 (68 reactions,
  the old "medium" network's exact equivalent), both backends
- an MC uncertainty-propagation run (C backend, small network)

All timings are *warm-cache* (the shipped ``primat/data/cache_plasma_weak/weak/*.txt``
n<->p weak-rate cache and ``primat/data/cache_plasma_weak/plasma/QED_*.txt`` tables already
match the default configuration's fingerprint) -- this is the timing a user
sees on every run after the first. The cold-cache cost (first run ever, or
first run after changing a field in ``primat.weak_rates.cache``'s
fingerprint) is *not* measured here since it is dominated by a multi-minute
``vegas`` Monte-Carlo integration when ``thermal_corrections=True`` (the
default); see :doc:`docs/howto/weak-rate-cache` for that cost breakdown.

Run from the repo root so that ``primat/data/`` resolves correctly:

    python runfiles/benchmark.py

Pass ``--quick`` (or set ``PRIMAT_BENCHMARK_QUICK``) to shrink the MC sample
count for a fast smoke test (see ``tests/test_runfiles.py``); the numbers
printed under ``--quick`` are not meaningful as a benchmark, only as a
"did it run" check. Each solve/MC call is repeated ``--repeats`` times
(default 3) and the *minimum* wall time is reported, to reduce noise from
OS scheduling/thermal throttling on the first call of a warm process.

Prints one line per row and then the whole table in Markdown, ready to paste
over docs/performance.md's. Writes nothing. Takes about five seconds at the
default ``--repeats 3``.
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from primat.backend import HAS_C_BACKEND, run_bbn, run_mc

_quick = "--quick" in sys.argv or bool(os.environ.get("PRIMAT_BENCHMARK_QUICK"))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--quick", action="store_true", help="shrink MC sample count for a smoke test")
parser.add_argument("--repeats", type=int, default=3, help="repeat each solve this many times, report the min")
args, _unknown = parser.parse_known_args()

REPEATS = args.repeats
NUM_MC = 5 if _quick else 100

NETWORKS = [
    ("small", {"network": "small"}),
    ("large, amax=8", {"network": "large", "amax": 8}),
]
BACKENDS = ["c", "python"] if HAS_C_BACKEND else ["python"]


def _time_min(fn, repeats=REPEATS):
    """Call ``fn()`` ``repeats`` times and return the minimum wall time [s].

    The minimum (rather than the mean) is reported because the first call in
    a freshly-started Python process pays JIT/import warm-up costs (numba
    kernel compilation, module-level table loads) that later calls don't --
    taking the min isolates the steady-state per-solve cost this table is
    meant to document.
    """
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def main():
    rows = []
    for net_label, net_kwargs in NETWORKS:
        for backend_name in BACKENDS:
            params = {**net_kwargs, "Omegabh2": 0.02242}
            elapsed = _time_min(lambda: run_bbn(params, force_backend=backend_name, progress=False))
            rows.append((f"{net_label} ({backend_name})", elapsed))
            print(f"{net_label:16s} {backend_name:6s} : {elapsed:.4f} s (min of {REPEATS})")

    if HAS_C_BACKEND:
        mc_params = {"network": "small", "Omegabh2": 0.02242}
        t0 = time.perf_counter()
        run_mc(NUM_MC, params=mc_params, force_backend="c", seed=0, progress=False)
        mc_elapsed = time.perf_counter() - t0
        rows.append((f"MC-{NUM_MC} (small, c)", mc_elapsed))
        print(f"MC-{NUM_MC:<9d} c      : {mc_elapsed:.4f} s (single run, n_jobs=-1)")

    print()
    print("Markdown table (paste into docs/performance.md):")
    print()
    print("| Run | Wall time |")
    print("|-----|-----------|")
    for label, elapsed in rows:
        print(f"| {label} | {elapsed:.3f} s |")


if __name__ == "__main__":
    main()
