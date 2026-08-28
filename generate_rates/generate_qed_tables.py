#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_qed_tables.py
======================
Standalone script to recompute the QED plasma-pressure correction tables
and write them to ``primat/data/cache_plasma_weak/plasma/`` (the shipped
location -- the two regenerable cache trees, n<->p weak rates and plasma
electron-thermo/QED tables, live together under ``cache_plasma_weak/``).

These tables store δP(T), dδP/dT, and d²δP/dT² — the finite-temperature
QED corrections to the EM plasma pressure that enter the background
evolution of the BBN code.  They were originally computed with
PRIMAT-Main.m (Mathematica); this script provides the equivalent Python
computation so the files can be regenerated without Mathematica.

The computation uses :mod:`primat.qed_pressure` which implements the
analytic formulas from PRIMAT-Main.m:

    δP(T) = δP_a(T)  [O(α), leading]
           + δP_{e3}(T)  [O(α^{3/2}), ring/plasmon]

(The O(α²) two-loop exchange term δP_b is available via --include-dPb
but is not included in the standard files.)

Usage::

    # From the repository root:
    python generate_rates/generate_qed_tables.py

    # Higher-resolution grid:
    python generate_rates/generate_qed_tables.py --n-pts 1000

    # Also compute the O(e^4) two-loop exchange term (very slow):
    python generate_rates/generate_qed_tables.py --include-dPb

The output files are written to ``primat/data/cache_plasma_weak/plasma/``:
  - ``QED_pressure_correction_e2.txt`` — T, δP_a, d(δP_a)/dT, d²(δP_a)/dT²  [O(e²)]
  - ``QED_pressure_correction_e3.txt`` — T, δP_{e3}, d(δP_{e3})/dT, d²(δP_{e3})/dT²  [O(e³)]

Each file carries a fingerprint header (:func:`primat.qed_pressure.qed_fingerprint`:
format version, ``constants_hash``, and the T grid), so a table computed with
different constants or bounds is detected and rebuilt by the loader rather
than used silently.

Reproducibility note: the *shipped* tables were produced by an earlier build
and agree with a fresh run to ~1e-6 relative — the last digit of the ``%.6E``
write format, i.e. far below any BBN tolerance but **not** byte-identical.
The fingerprint deliberately keys on the constants and the grid, not on the
table's content, so it cannot (and need not) flag that difference. Overwriting
the shipped files with a fresh run is therefore physically a no-op; do it only
when the formulas or constants actually change.

Physical background
-------------------
The QED interaction pressure corrects the ideal-gas (photon + e±) EM
plasma equation of state.  It is decomposed into an O(e²) leading term
and an O(e³) ring/plasmon term following Phys. Rep. §II.E (PRIMAT
variables ``dPa``, ``dPe3``).  At T = 10 MeV:
  δP_a   ≈ −17  MeV⁴  (negative: reduces pressure)
  δP_{e3} ≈ +0.3 MeV⁴ (positive: ring contribution)
  total  ≈ −16.7 MeV⁴

Reference
---------
Pitrou, Coc, Uzan & Vangioni, Phys. Rep. 2018 (arXiv:1806.11095), §II.E
PRIMAT-Main.m: ``dPa``, ``dPe3``, ``dPb`` definitions
"""

import sys
import os
import argparse
import time

# Ensure the repo root is on sys.path so that primat is importable.
_repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from primat.qed_pressure import compute_qed_pressure_tables, save_qed_tables

# Where the shipped QED tables live, i.e. where a plain (no --output-dir) run
# writes.  A module-level constant rather than an expression buried in main()
# so tests/test_network_generation.py can assert it still points at a real
# directory: the previous value (primat/rates/plasma) had not existed since
# the data tree moved under primat/data/, and because the script created it on
# demand a regeneration silently wrote tables nothing would ever read.
# Mirrors primat.cache_utils.plasma_cache_dir's shipped branch.
DEFAULT_PLASMA_DIR = os.path.join(
    _repo_root, "primat", "data", "cache_plasma_weak", "plasma")


def main():
    parser = argparse.ArgumentParser(
        description="Recompute QED plasma-pressure correction tables.")
    parser.add_argument("--n-pts", type=int, default=500,
                        help="Number of log-spaced temperature grid points "
                             "(default: 500, matching the PRIMAT file).")
    parser.add_argument("--T-min", type=float, default=1e-3,
                        help="Minimum temperature [MeV] (default: 1e-3).")
    parser.add_argument("--T-max", type=float, default=100.,
                        help="Maximum temperature [MeV] (default: 100).")
    parser.add_argument("--include-dPb", action="store_true",
                        help="Also compute the O(e^4) two-loop exchange "
                             "correction δP_b (very slow: ~10 s per point).")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory for QED_*.txt files.  Default: "
                             "<repo_root>/primat/data/cache_plasma_weak/plasma/")
    args = parser.parse_args()

    plasma_dir = args.output_dir or DEFAULT_PLASMA_DIR
    if not os.path.isdir(plasma_dir):
        # A missing default target means the data tree moved again; say so
        # instead of silently creating a stray directory and writing tables
        # nothing will ever read (the failure mode this check replaces).
        if args.output_dir is None:
            raise SystemExit(
                f"error: the shipped plasma-table directory {plasma_dir} does "
                f"not exist -- has the data tree moved? Pass --output-dir "
                f"explicitly to write elsewhere.")
        os.makedirs(plasma_dir, exist_ok=True)

    print("Computing QED plasma-pressure tables:")
    print(f"  T grid: {args.T_min:.2e}–{args.T_max:.2e} MeV, {args.n_pts} points")
    print(f"  include δP_b (O(e^4)): {args.include_dPb}")
    print(f"  output: {plasma_dir}/")
    print()

    t0 = time.time()
    tables = compute_qed_pressure_tables(
        T_min=args.T_min,
        T_max=args.T_max,
        n_pts=args.n_pts,
        include_dPb=args.include_dPb,
        verbose=True,
    )
    dt = time.time() - t0
    print(f"\nComputation finished in {dt:.1f} s")

    save_qed_tables(tables, plasma_dir, verbose=True)
    print("\nDone.")


if __name__ == "__main__":
    main()
