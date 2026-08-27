"""Centralised reference BBN observables shared across the test suite.

Single source of truth for ``tests/test_cli.py``, ``tests/test_gui.py``,
``tests/test_regression.py`` and ``tests/test_docs_consistency.py``, so that a
routine default-parameter tweak (e.g. ``rate_grid_npts``, commit ``e00f062``)
only needs a tolerance check here, not a hunt for every literal pin scattered
across the suite.

Three groups live here, each with its own tolerances:

1. the **default-precision** small-network observables
   (``NEFF_REFERENCE``/``YPBBN_REFERENCE``/``DOH_REFERENCE``) — what
   ``primat.cli.main([])``, ``primat-gui``'s "Run BBN" and
   ``runfiles/primat_run.py`` actually print;
2. the **high-precision** observables (``REF_SMALL_*``/``REF_LARGE8_*``)
   produced by ``runfiles/primat_reference_run.py``, which
   ``tests/test_regression.py``'s ``reference`` tier reproduces and
   ``tests/README.md``'s "Validation reference" tables quote;
3. the **per-nuclide** default-precision final abundances
   (``NUCLIDE_REFERENCE``), quoted in the same README section.

The tolerances mirror the bounds published in ``tests/README.md``'s
"Validation reference" section: a result outside them indicates a *physics*
regression, not test brittleness. ``NEFF_ABS_TOL`` is not separately published
but is given the same ``1e-5`` margin used for ``YP``, since both observables
are driven by the same n<->p weak-rate / background machinery.
"""

# ---------------------------------------------------------------------------
# 1. Default-precision small-network run
# ---------------------------------------------------------------------------
# network="small", spectral_distortions=True, nuclear_qed_corrections=True --
# the PRIMATConfig defaults -- as produced by `primat.cli.main([])`,
# `primat-gui`'s default "Run BBN", and `runfiles/primat_run.py`.
#
# Snapshotted on the auto backend (C when available), which is what these tests
# actually invoke. The pure-Python backend's own values are PY_YPBBN_REFERENCE /
# PY_DOH_REFERENCE below -- 1.1e-07 and 1.6e-10 away, so both backends sit
# inside the bounds below.
NEFF_REFERENCE  = 3.0439772986
YPBBN_REFERENCE = 0.24699907   # primat.cli.main() default (auto backend), Omegabh2=0.02242
DOH_REFERENCE   = 2.4358767e-5  # primat.cli.main() default (auto backend), Omegabh2=0.02242

# Tolerances (tests/README.md: "A result outside these bounds indicates a
# regression").
NEFF_ABS_TOL  = 1e-5
YPBBN_ABS_TOL = 1e-5
DOH_ABS_TOL   = 3e-9

# The same run on the pure-Python backend. Quoted alongside the C values in
# docs/index.md's quick start, which test_docs_consistency.py checks against
# these constants; pinned to a live solve by test_regression.py's solved_small
# tier, which runs primat.main.PRIMAT (i.e. the Python backend) against the
# constants above.
PY_YPBBN_REFERENCE = 0.24699896
PY_DOH_REFERENCE   = 2.4358605e-5

# ---------------------------------------------------------------------------
# 2. High-precision reference observables
# ---------------------------------------------------------------------------
# tests/README.md "Validation reference", produced by
# runfiles/primat_reference_run.py. Its full settings live in that script and
# are mirrored by tests/test_regression.py's _REF_PARAMS, which must stay
# byte-identical to them: the four it once omitted (rate_grid_npts,
# sampling_nTOp_thermal_per_decade, vegas_n_eval, vegas_n_itn) moved
# large+amax=8's D/H by 2.0e-08, i.e. 6.6x the +/-3e-9 bound these constants
# are checked against. Single source for that reference tier AND the table
# parser in tests/test_docs_consistency.py -- update all three places together
# (the parser test fails if the README table and these constants drift).
REF_SMALL_YPBBN  = 0.24699844
REF_SMALL_DOH    = 2.4358977e-5
REF_LARGE8_YPBBN = 0.24700179
REF_LARGE8_DOH   = 2.4366098e-5

# Bound for a *routine* `runfiles/primat_run.py` check against the two
# constants above. That script runs at the default numerical_precision=1e-7,
# not the reference run's 1e-10, so the tight DOH_ABS_TOL=3e-9 above does NOT
# apply to it: the routine run lands 6.2e-10 (C backend) and 7.8e-10 (Python
# backend) below REF_LARGE8_DOH. 2e-8 leaves ample headroom over a margin that
# is a property of this platform and network, while staying ~120x tighter than
# the loose sanity tier in tests/test_regression.py. See tests/README.md's
# "Validation reference".
ROUTINE_RUN_DOH_ABS_TOL = 2e-8
ROUTINE_RUN_YPBBN_ABS_TOL = 1e-5

# ---------------------------------------------------------------------------
# 3. Per-nuclide final abundances (default precision, auto backend)
# ---------------------------------------------------------------------------
# Final abundances Y_s = n_s/n_b of the small-network nuclides at the end of
# BBN, for the three networks tests/README.md tabulates. These are abundances
# per baryon, not mass fractions: the mass fraction is A_s Y_s, so YPBBN is
# 4 * Y_He4. Quoted to 7
# significant figures; the two backends agree on them to <=2.2e-05 relative,
# so NUCLIDE_REL_TOL below is what a check should actually use -- the last two
# quoted digits are backend-dependent.
#
# Snapshotted on the auto backend, alongside the high-precision constants
# above; both move together whenever the solver's numerics change.
NUCLIDE_REFERENCE = {
    #        small          large, amax=8   large
    "n":   (3.997246e-16, 3.996325e-16, 3.996357e-16),
    "p":   (7.529408e-01, 7.529374e-01, 7.529374e-01),
    "H2":  (1.834071e-05, 1.834568e-05, 1.834577e-05),
    "H3":  (5.851941e-08, 5.838979e-08, 5.839007e-08),
    "He4": (6.174977e-02, 6.175059e-02, 6.175060e-02),
    "Li7": (2.181376e-11, 9.178228e-11, 9.178164e-11),
    "Be7": (3.966454e-10, 3.223689e-10, 3.223658e-10),
}

# Column order of NUCLIDE_REFERENCE's tuples, matching the README table.
NUCLIDE_COLUMNS = ("small", "large_amax8", "large")

# Relative bound for a live default-precision solve against NUCLIDE_REFERENCE.
# Measured cross-backend spread is <=2.2e-05 and the accumulated drift from
# successive numerics improvements <=6.2e-05, so 1e-4 is the smallest round
# bound that does not make such an improvement a test failure -- while still
# catching, say, the ~1750x jump in `n` that a regression of the reverse-rate
# clamp (primat/network_data.py, "exothermic blow-up") would produce.
NUCLIDE_REL_TOL = 1e-4

# Baryon-number conservation bound, sum_s A_s Y_s == 1. Measured 1.6e-12 for
# all three networks (2026-08-05); pinned at 1e-10 so a real stoichiometry
# leak fails rather than hiding under a decorative bound.
BARYON_ABS_TOL = 1e-10
