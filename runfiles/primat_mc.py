# -*- coding: utf-8 -*-
"""
primat_mc.py
============
Heavily-commented demo of primat's Monte-Carlo nuclear-rate/tau_n
uncertainty propagation, in the house style of ``primat_run_explanatory.py``:
run ``run_mc()``, print each observable's ``value +/- sigma`` at the
decimal counts tests/README.md's validation tables use (Neff 8, YP 8,
D/H 7, Li7/H 6 -- enough to resolve the ~1e-3 flag-level effects), show
the joint (covariance and
correlation) uncertainty between abundances -- both the full matrix and the
scalar two-name form -- and write the three ``<prefix>_*.tsv`` files
(samples / covariance / correlation).

Run from the repo root so that the shipped ``data/`` data resolve correctly:

    python runfiles/primat_mc.py

Pass ``--quick`` (or set the ``PRIMAT_MC_QUICK`` environment variable) to use
a small sample count suitable for a fast smoke test (see
``tests/test_runfiles.py``); the default sample count is large enough for the
printed correlations/uncertainties to be meaningful.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from primat.backend import dump_mc_correlation, dump_mc_covariance, dump_mc_samples, run_mc

# `--quick`/PRIMAT_MC_QUICK trade statistical precision for speed -- useful
# for CI/smoke-testing this script, not for a real uncertainty estimate.
_quick = "--quick" in sys.argv or bool(os.environ.get("PRIMAT_MC_QUICK"))
NUM_MC = 20 if _quick else 500

# Same default cosmology as primat_run_explanatory.py: Omegabh2 = 0.02242 is
# the code's own DEFAULT_PARAMS value (Planck 2018 TT,TE,EE+lowE+lensing).
PARAMS = {"Omegabh2": 0.02242}

# `quantities=None` still returns every standard observable (Neff, YPBBN,
# YPCMB, He4oH, DoH, He3oH, He3oHe4, Li7oH, Li6oLi7, YCNO) plus every tracked
# nuclide's final Y, at no extra solving cost (see run_mc's docstring) --
# so the three dump_mc_* files below are always complete even though we
# only print a handful of quantities here.
mc = run_mc(NUM_MC, params=PARAMS, seed=0)

print(f"Monte-Carlo BBN uncertainty propagation, N = {NUM_MC} samples "
      f"(seed=0), Omegabh2 = {PARAMS['Omegabh2']}")
print()

# --- Per-observable value +/- sigma, at the validation tables' precision ---
print(f"Neff        = {mc['Neff'].mean:.8f} +/- {mc['Neff'].std:.8f}")
print(f"YP (BBN)    = {mc['YPBBN'].mean:.8f} +/- {mc['YPBBN'].std:.8f}")
print(f"D/H         = {mc['DoH'].mean:.7e} +/- {mc['DoH'].std:.7e}")
print(f"Li7/H       = {mc['Li7oH'].mean:.6e} +/- {mc['Li7oH'].std:.6e}")
print()

# --- Joint uncertainty: full 4x4 correlation of the main products ---------
# The off-diagonal terms are exactly what a joint (YPBBN, DoH, ...)
# likelihood needs: e.g. YP and D/H are driven by the same underlying
# nuclear-rate/tau_n samples and are therefore correlated, not independent.
main_products = ["YPBBN", "DoH", "He3oHe4", "Li7oH"]
names = mc.quantity_names()
idx = [names.index(q) for q in main_products]
R = mc.corr()[idx][:, idx]

print("Correlation matrix (YPBBN, DoH, He3oHe4, Li7oH):")
header = "".join(f"{q:>10}" for q in main_products)
print(f"{'':>10}{header}")
for i, q in enumerate(main_products):
    print(f"{q:>10}" + "".join(f"{R[i, j]:10.3f}" for j in range(len(main_products))))
print()

# --- Scalar two-name access: mc.cov()/mc.corr() also take two quantity
# names directly, without slicing the full matrix yourself. -----------------
print(f"Scalar access: mc.cov('YPBBN', 'DoH')  = {mc.cov('YPBBN', 'DoH'):.6e}")
print(f"Scalar access: mc.corr('YPBBN', 'DoH') = {mc.corr('YPBBN', 'DoH'):.6f}")
print()

# --- Write the three MC output files --------------------------------------
os.makedirs("results", exist_ok=True)
prefix = "results/output_mc"
with open(f"{prefix}_samples.tsv", "w") as f:
    f.write(dump_mc_samples(mc))
with open(f"{prefix}_covariance.tsv", "w") as f:
    f.write(dump_mc_covariance(mc))
with open(f"{prefix}_correlation.tsv", "w") as f:
    f.write(dump_mc_correlation(mc))
print(f"Wrote {prefix}_samples.tsv, {prefix}_covariance.tsv, "
      f"{prefix}_correlation.tsv")
