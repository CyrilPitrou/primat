"""
Smoke test for the example scripts in ``runfiles/``, plus the one numeric
check that backs the documented validation workflow.

``runfiles/primat_run.py`` is the canonical "run this to validate a change"
entry point (tests/README.md's "Validation reference"), and
``primat_compare.py``/``primat_run_explanatory.py`` are further worked
examples -- but nothing in the test suite used to execute them
(``tests/reference_values.py`` only mirrors their expected numbers, and
``tests/test_docs_consistency.py`` only string-checks parameter names in
``primat_reference_run.py``). This means an import-path bug or an API rename
in ``primat.backend``/``primat.main`` could break these scripts silently until
a human runs one by hand.

``test_primat_run_matches_the_validation_reference`` closes the other half of
that hole: it parses the script's *printed numbers* and checks them against
the published reference, so the workflow tests/README.md tells a contributor
to follow by hand is also run by CI. It deliberately uses
``ROUTINE_RUN_*_ABS_TOL`` rather than the tight ``DOH_ABS_TOL``: this script
runs at the default ``numerical_precision=1e-7``, not the reference run's
1e-10, so the +-3e-9 reference bound does not apply to it (the C backend lands
3.8e-9 below the reference D/H -- outside +-3e-9 and correctly inside +-2e-8).
See tests/README.md's "Which tolerance applies to which command".

Each script is run as a real subprocess (``python <script>``) rather than
imported, since none of them are wrapped in a ``main()``/``if __name__``
guard -- they are plain top-level scripts meant to be run directly. The
subprocess's *working directory* is a throwaway ``tmp_path``, not the repo
root: all three scripts write ``results/*.tsv``/``*.dat`` to a path relative
to the *current directory* (not ``__file__``), so running from ``tmp_path``
keeps every output out of the tracked (albeit gitignored) ``results/``
directory at the repo root.

Deliberately excluded:

* ``primat_reference_run.py`` -- several minutes by design (high-precision
  reference run for updating tests/README.md's benchmarks), out of scope for a
  smoke test.
* ``generate_weak_rate_caches.py`` -- (re)writes the fingerprinted
  ``rates/weak/nTOp_*.txt`` cache files that are force-added to git; not
  safe to run unattended even from a throwaway cwd (it resolves the cache
  directory relative to the installed package, not cwd).
* ``generate_table_CLASS_CAMB.py`` -- needs an external CLASS/CAMB
  installation and a multi-hour Monte Carlo table generation run.
"""
import re
import subprocess
import sys
from pathlib import Path

import pytest

from tests.reference_values import (REF_LARGE8_YPBBN, REF_LARGE8_DOH,
                                    ROUTINE_RUN_YPBBN_ABS_TOL,
                                    ROUTINE_RUN_DOH_ABS_TOL)

pytestmark = [pytest.mark.slow]

RUNFILES_DIR = Path(__file__).resolve().parents[1] / "runfiles"

FAST_RUNFILES = [
    "primat_run.py",
    "primat_run_explanatory.py",
    "primat_compare.py",
]

# primat_mc.py needs the "--quick" flag (20 samples instead of 500) to stay
# within this smoke test's time budget -- it is not one of FAST_RUNFILES
# since it takes an extra CLI argument the others don't.
MC_RUNFILE = "primat_mc.py"

# benchmark.py also needs "--quick" (5 MC samples
# instead of 100; still 3 repeats per network/backend solve) to stay within
# this smoke test's time budget -- its numbers under --quick are not
# meaningful as a benchmark, only as a "did it run" check (see its module
# docstring).
BENCHMARK_RUNFILE = "benchmark.py"


@pytest.mark.parametrize("name", FAST_RUNFILES)
def test_runfile_executes_cleanly(name, tmp_path):
    """Run an example script as a subprocess; fail on a nonzero exit code
    or a traceback, and sanity-check it printed the headline observables."""
    result = subprocess.run(
        [sys.executable, str(RUNFILES_DIR / name)],
        cwd=tmp_path,
        capture_output=True,
        text=True, encoding="utf-8",
        timeout=120,
    )

    assert result.returncode == 0, (
        f"{name} exited with code {result.returncode}\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )
    assert "Traceback" not in result.stderr, result.stderr

    # Every one of these scripts prints Neff and D/H (spelled "D/H") for at
    # least the small network -- a quick sanity check that it actually ran
    # the solver rather than exiting early/silently.
    assert "Neff" in result.stdout
    assert "D/H" in result.stdout


@pytest.mark.solve
def test_primat_run_matches_the_validation_reference(tmp_path):
    """primat_run.py's printed YP/(D/H) must match the published reference.

    GOAL: make the documented "after any modification, run primat_run.py and
    check the output against these tables" workflow (tests/README.md,
    "Validation reference") an automated check rather than an honour-system
    one. The script solves ``network="large", amax=8`` at the default
    ``numerical_precision=1e-7``, so it is compared against the *routine*
    tolerance column, not the reference one -- see this module's docstring.

    Parsing stdout rather than importing the script is deliberate: what a
    contributor actually reads is the printed line, so that is what is pinned.
    """
    result = subprocess.run(
        [sys.executable, str(RUNFILES_DIR / "primat_run.py")],
        cwd=tmp_path, capture_output=True, text=True, encoding="utf-8", timeout=120,
    )
    assert result.returncode == 0, result.stderr

    def _printed(label):
        # Lines look like " YP (BBN) =  0.24700262002081247".
        m = re.search(rf"^\s*{re.escape(label)}\s*=\s*([0-9.eE+-]+)\s*$",
                      result.stdout, re.M)
        assert m, f"{label!r} not found in primat_run.py output:\n{result.stdout}"
        return float(m.group(1))

    assert _printed("YP (BBN)") == pytest.approx(
        REF_LARGE8_YPBBN, abs=ROUTINE_RUN_YPBBN_ABS_TOL)
    assert _printed("D/H") == pytest.approx(
        REF_LARGE8_DOH, abs=ROUTINE_RUN_DOH_ABS_TOL)


def test_primat_mc_runfile_executes_cleanly(tmp_path):
    """Run primat_mc.py --quick (20 MC samples) as a subprocess; fail on a
    nonzero exit code or a traceback, and sanity-check it printed the
    headline observables and wrote the three MC output files."""
    result = subprocess.run(
        [sys.executable, str(RUNFILES_DIR / MC_RUNFILE), "--quick"],
        cwd=tmp_path,
        capture_output=True,
        text=True, encoding="utf-8",
        timeout=120,
    )

    assert result.returncode == 0, (
        f"{MC_RUNFILE} exited with code {result.returncode}\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )
    assert "Traceback" not in result.stderr, result.stderr
    assert "Neff" in result.stdout
    assert "D/H" in result.stdout
    assert "Correlation matrix" in result.stdout

    for suffix in ("samples", "covariance", "correlation"):
        assert (tmp_path / "results" / f"output_mc_{suffix}.tsv").exists()


def test_benchmark_runfile_executes_cleanly(tmp_path):
    """Run benchmark.py --quick (5 MC samples, 3 repeats/solve) as a
    subprocess; fail on a nonzero exit code or a traceback, and sanity-check
    it printed the regenerated Markdown timing table."""
    result = subprocess.run(
        [sys.executable, str(RUNFILES_DIR / BENCHMARK_RUNFILE), "--quick"],
        cwd=tmp_path,
        capture_output=True,
        text=True, encoding="utf-8",
        timeout=120,
    )

    assert result.returncode == 0, (
        f"{BENCHMARK_RUNFILE} exited with code {result.returncode}\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )
    assert "Traceback" not in result.stderr, result.stderr
    assert "Markdown table" in result.stdout
    assert "| Run | Wall time |" in result.stdout
