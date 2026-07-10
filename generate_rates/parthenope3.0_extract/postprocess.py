#!/usr/bin/env python3
"""Turn the raw Parthenope generator output into primat rate tables.

For each reaction we read the raw ``out_<name>.dat`` (columns: T9, f, error)
produced by the compiled Parthenope-3.0 kernel, repair the unphysical low-T9
tail (where the published polynomial fits extrapolate to negative rates well
below the temperatures Parthenope ever evaluates them), and write a clean
3-column table ``<name>_parthenope3.0.txt`` with a documentary header.

Low-T9 repair: below the lowest T9 where the raw fit is positive and finite,
the rate is filled by *log-log linear extrapolation* of the two lowest valid
points (continuing the steep Gamow fall-off as T9 -> 0), and the multiplicative
error is held constant at its value on the validity floor.  This region is
dynamically irrelevant (the reactions are frozen out there) but the values must
stay positive and monotone so primat's log-log resampler can ingest them.

Grid: after repair, each table is resampled onto primat's master T9 grid
(``rate_grid_{npts,T9_min,T9_max}`` in ``primat.config.DEFAULT_PARAMS``,
currently 1000 log-spaced points from 1e-3 to 10 GK) using the *exact*
runtime resampler ``primat.network_data._resample_rate_table`` (log-log cubic).
Writing the tables already on the master grid makes ``load_network``'s
load-time resample a no-op, so the shipped file and what the solver sees are
identical -- and keeps this stream in lock-step with the AC2024 ``_primat``
tables, which ``convert_ac2024_rates.py`` also writes on the master grid.
``tests/test_rate_tables_on_grid.py`` enforces that every shipped table stays
on this grid.
"""
import os
import sys
import numpy as np

# Repo root is three levels up from generate_rates/parthenope3.0_extract/.
_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from primat.config import DEFAULT_PARAMS
from primat.network_data import _resample_rate_table

OUT_DIR = os.path.join(_REPO_ROOT, "primat", "data", "nuclear", "tables")

# primat's master T9 grid -- the single source of truth is DEFAULT_PARAMS, so
# this generator can never drift from the runtime grid the solver resamples to.
STD_GRID = np.logspace(
    np.log10(DEFAULT_PARAMS["rate_grid_T9_min"]),
    np.log10(DEFAULT_PARAMS["rate_grid_T9_max"]),
    DEFAULT_PARAMS["rate_grid_npts"],
)

# primat name -> (header arrow string, Parthenope source/branch)
META = {
    "npTOdg":     ("n + p > d + g",        "Serpico et al 2004"),
    "dpTOHe3g":   ("d + p > He3 + g",      "Pisanti et al 2020 (PIS2020, default)"),
    "ddTOHe3n":   ("d + d > He3 + n",      "Pisanti et al 2020 (PIS2020, default)"),
    "ddTOtp":     ("d + d > t + p",        "Pisanti et al 2020 (PIS2020, default)"),
    "tpTOag":     ("t + p > a + g",        "Caughlan-Fowler 1988"),
    "tdTOan":     ("t + d > a + n",        "Serpico et al 2004"),
    "taTOLi7g":   ("t + a > Li7 + g",      "Serpico et al 2004"),
    "He3nTOtp":   ("He3 + n > t + p",      "Serpico et al 2004"),
    "He3dTOap":   ("He3 + d > a + p",      "Serpico et al 2004"),
    "He3aTOBe7g": ("He3 + a > Be7 + g",    "Serpico et al 2004"),
    "Be7nTOLi7p": ("Be7 + n > Li7 + p",    "Serpico et al 2004"),
    "Li7pTOaa":   ("Li7 + p > a + a",      "Serpico et al 2004"),
}


def repair(t9, f, err):
    """Return (f, err) with the invalid low-T9 tail extrapolated log-log.

    A point is "valid" when its rate is finite and strictly positive and its
    error is finite and >= 1.  We take the contiguous valid block that reaches
    the high-T9 end of the grid; everything below its lower edge is refilled.
    """
    valid = np.isfinite(f) & (f > 0.0) & np.isfinite(err) & (err >= 1.0)
    # Lower edge = first index from which all remaining points are valid.
    # (BBN-relevant and high-T region; the bad tail is purely at the bottom.)
    k = 0
    for i in range(len(valid) - 1, -1, -1):
        if not valid[i]:
            k = i + 1
            break
    if k == 0:
        return f.copy(), err.copy()  # nothing to repair
    f = f.copy()
    err = err.copy()
    # Log-log slope from the two lowest valid points.
    x0, x1 = np.log(t9[k]), np.log(t9[k + 1])
    y0, y1 = np.log(f[k]),  np.log(f[k + 1])
    slope = (y1 - y0) / (x1 - x0)
    for i in range(k):
        f[i] = np.exp(y0 + slope * (np.log(t9[i]) - x0))
        err[i] = err[k]
    return f, err


def main():
    for name, (arrow, ref) in META.items():
        raw = np.loadtxt(os.path.join("/tmp/pgen", f"out_{name}.dat"))
        t9, f, err = raw[:, 0], raw[:, 1], raw[:, 2]
        nbad = int(np.sum(~((np.isfinite(f) & (f > 0) & np.isfinite(err) & (err >= 1.0)))))
        f, err = repair(t9, f, err)
        # Resample onto primat's master grid with the exact runtime resampler
        # (log-log cubic), so the shipped file needs no load-time interpolation.
        f = _resample_rate_table(t9, f, STD_GRID)
        err = _resample_rate_table(t9, err, STD_GRID)
        t9 = STD_GRID
        path = os.path.join(OUT_DIR, f"{name}_parthenope3.0.txt")
        with open(path, "w") as fh:
            fh.write(f"# {arrow}   [{name}]   ref=Parthenope3.0 ({ref})\n")
            fh.write("# Forward rate N_A<sigma v> [cm^3 mol^-1 s^-1] extracted from\n")
            fh.write("# parthenope3.0.f (default PIS2020 nuclear-rate selection).\n")
            fh.write("# error column = sqrt((1+drate_up)/(1+drate_lo)), the full Parthenope\n")
            fh.write("# 1-sigma multiplicative envelope (statistical fp/fm + systematic floor).\n")
            if nbad:
                fh.write(f"# NOTE: lowest {nbad} T9 points (below the fit validity floor, where\n")
                fh.write("#   the published polynomial extrapolates non-physically) are filled by\n")
                fh.write("#   log-log extrapolation; this region is dynamically frozen for BBN.\n")
            fh.write("# T9                 rate                error\n")
            for a, b, c in zip(t9, f, err):
                fh.write(f"{a:.6e} {b:.6e} {c:.6e}\n")
        print(f"{name:12s} -> {os.path.basename(path)}  (repaired {nbad} low-T9 pts)")


if __name__ == "__main__":
    main()
