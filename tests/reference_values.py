"""Centralised reference BBN observables for the default (small-network) run.

Single source of truth for ``tests/test_cli.py``, ``tests/test_gui.py`` and
any other test asserting the default-config small-network result, so that a
routine default-parameter tweak (e.g. ``rate_grid_npts``, commit ``e00f062``)
only needs a tolerance check here, not a hunt for every literal pin scattered
across the suite.

The tolerances mirror CLAUDE.md's "Validation before committing" table
(``YP (BBN)`` and ``D/H``, small network): a result outside these bounds
indicates a *physics* regression, not test brittleness. ``NEFF_ABS_TOL`` is
not separately documented in CLAUDE.md but is given the same ``1e-5``
margin used for ``YP``, since both observables are driven by the same
n<->p weak-rate / background machinery.
"""

# Default small-network run (network="small", spectral_distortions=True,
# nuclear_qed_corrections=True -- the PRIMATConfig defaults), as produced by
# `primat.cli.main([])`, `primat-gui`'s default "Run BBN", and
# `runfiles/primat_run.py`. Re-snapshotted 2026-07-08 after the default
# `Omegabh2` was changed to 0.02242 (Planck 2018 + BAO); values are the
# auto-backend (C, when available) CLI output, matching what these tests
# actually invoke.
NEFF_REFERENCE  = 3.0439772986
YPBBN_REFERENCE = 0.24699701   # primat.cli.main() default (auto backend), Omegabh2=0.02242
DOH_REFERENCE   = 2.435908e-5  # primat.cli.main() default (auto backend), Omegabh2=0.02242

# Tolerances (CLAUDE.md: "A result outside these bounds indicates a regression").
NEFF_ABS_TOL  = 1e-5
YPBBN_ABS_TOL = 1e-5
DOH_ABS_TOL   = 3e-9

# Per-nuclide final mass fractions, small network (primat.cli.main() default, auto backend).
P_REFERENCE   = 7.529428e-01
HE4_REFERENCE = 6.174925e-02
NUCLIDE_ABS_TOL = 1e-4  # mirrors the table's own documented precision

# High-precision reference values (tests/README.md "Validation reference"),
# produced by runfiles/primat_reference_run.py (numerical_precision=1e-10,
# sampling_temperature_per_decade=2000, sampling_nTOp_per_decade=125,
# T_start_cosmo_MeV=100, rate_grid_npts=4000). Single source for
# tests/test_regression.py's reference tier AND the table parser in
# tests/test_docs_consistency.py -- update all three places together (the
# parser test fails if the README table and these constants drift).
REF_SMALL_YPBBN  = 0.24699814
REF_SMALL_DOH    = 2.43589e-5
REF_LARGE8_YPBBN = 0.24700149
REF_LARGE8_DOH   = 2.43660e-5
