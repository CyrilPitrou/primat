# -*- coding: utf-8 -*-
"""The two backends must narrate a run identically.

GOAL: pin the one sync rule that had no test. A verbose run is how a user (and
a reviewer chasing a parity bug) sees what each backend did, so adding a stage
to one and not the other is a divergence like any other -- but nothing checked
it, and by the time this test was written the two streams had already drifted
in two places: C printed each era's nuclide list after *both* "network" lines
rather than after its own, and the LT finish line said "species" on C against
"nuclides" on Python.

The comparison is line for line on the message stream, after normalising the
three things that are *meant* to differ: the ``-py``/``-c`` tag suffix, the
name each backend reports for itself, and wall-clock timings.
"""
import re
import subprocess
import sys

import pytest

from primat.backend import HAS_C_BACKEND


pytestmark = pytest.mark.solve


# Python-only lines, excluded by content with a reason. numba is a Python-only
# optional dependency (the C backend is compiled, so it has no JIT to report).
PYTHON_ONLY = (
    "numba detected",
    "numba not detected",
)

_TAG = re.compile(r"^\[([a-z]+)(?:-py|-c)\]")
_SECONDS = re.compile(r"\d+\.\d+ s")
_BACKEND_LINE = re.compile(r"^\[opts\] backend\s*=.*$")


def _message_stream(backend):
    """The normalised ``[tag] message`` lines a verbose run prints."""
    proc = subprocess.run(
        [sys.executable, "-m", "primat.cli", "--backend", backend, "--verbose"],
        capture_output=True, text=True, encoding="utf-8", timeout=600,
    )
    assert proc.returncode == 0, proc.stderr
    lines = []
    for line in (proc.stdout + proc.stderr).splitlines():
        if not line.startswith("["):
            continue
        if any(marker in line for marker in PYTHON_ONLY):
            continue
        line = _TAG.sub(r"[\1]", line)          # [rates-py]/[rates-c] -> [rates]
        line = _SECONDS.sub("<t> s", line)      # wall-clock timings
        line = _BACKEND_LINE.sub("[opts] backend = <backend>", line)
        lines.append(line.rstrip())
    assert lines, f"{backend} backend printed no tagged messages"
    return lines


@pytest.mark.skipif(not HAS_C_BACKEND, reason="primat._primat_c is not built")
def test_backends_narrate_a_run_identically():
    """Same stages, same wording, same order on both backends."""
    # Warm the shared caches first. Both backends announce whether they loaded
    # the electron-thermo/QED tables or computed them, so a cold cache left by
    # an earlier test would make whichever backend ran first say something
    # different -- a failure about cache state, not about parity.
    _message_stream("python")

    py = _message_stream("python")
    c = _message_stream("c")

    if py != c:
        import difflib
        diff = "\n".join(difflib.unified_diff(py, c, "python", "c", lineterm=""))
        pytest.fail(
            "the two backends' verbose output has diverged. Add or mirror the "
            "matching cpr_log()/print() on the other side, or -- if the line is "
            "genuinely backend-specific -- add it to PYTHON_ONLY with a "
            f"reason.\n\n{diff}"
        )


@pytest.mark.skipif(not HAS_C_BACKEND, reason="primat._primat_c is not built")
def test_the_stream_is_long_enough_to_be_worth_comparing():
    """Guard the guard: a normalisation bug that empties the stream must fail.

    Both filters above drop lines, so an over-eager one could reduce the
    comparison to a handful of lines that trivially match.
    """
    assert len(_message_stream("python")) >= 20
