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
import os
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


# Off-default configurations the stream is compared over. Each one reaches
# stages the default point never does -- a bigger network's reaction listing, a
# different background mode, a decay era, a weak-rate table computed rather
# than loaded, a QED recompute that writes files. The default point alone left
# a real divergence unseen: the C backend never printed the four-line "[QED]
# Tables written to ..." block, and its ellipsis was ASCII where Python's was
# U+2026.
OFF_DEFAULT_CONFIGS = {
    "default": [],
    "large_amax8": ["--set", "network=large", "--set", "amax=8"],
    "small_parthenope": ["--set", "network=small_parthenope"],
    "external_scale_factor": ["--set", "external_scale_factor=True"],
    "no_QED": ["--set", "QED_corrections=False"],
    "born_only": ["--set", "radiative_corrections=False",
                  "--set", "finite_mass_corrections=False",
                  "--set", "thermal_corrections=False",
                  "--set", "spectral_distortions=False"],
    "decay_era": ["--set", "decay_era=True"],
    "recompute_qed": ["--set", "recompute_qed_corrections=True"],
}


# The child is read back as UTF-8, so it must write UTF-8. A Python child
# encodes a redirected stream in the locale encoding, cp1252 on Windows, where
# the decay listing's "s⁻¹" becomes the lone byte 0xb9: the parent's
# read then raises UnicodeDecodeError and loses that whole stream. The C
# backend writes its bytes straight through and is unaffected.
_UTF8_ENV = dict(os.environ, PYTHONIOENCODING="utf-8")


def _message_stream(backend, args=()):
    """The normalised ``[tag] message`` lines a verbose run prints."""
    proc = subprocess.run(
        [sys.executable, "-m", "primat.cli", "--backend", backend, "--verbose"]
        + list(args),
        capture_output=True, text=True, encoding="utf-8", timeout=600,
        env=_UTF8_ENV,
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
@pytest.mark.parametrize("config", sorted(OFF_DEFAULT_CONFIGS),
                          ids=sorted(OFF_DEFAULT_CONFIGS))
def test_backends_narrate_a_run_identically(config, tmp_path):
    """Same stages, same wording, same order on both backends.

    Run over OFF_DEFAULT_CONFIGS, not only the default point: the rule is that
    the two backends narrate *a run*, and a stage only one configuration
    reaches is exactly where a missing mirror hides.
    """
    args = ["--set", f"cache_dir={tmp_path}"] + OFF_DEFAULT_CONFIGS[config]
    # Warm the shared caches first. Both backends announce whether they loaded
    # the electron-thermo/QED tables or computed them, so a cold cache left by
    # an earlier test would make whichever backend ran first say something
    # different -- a failure about cache state, not about parity.
    _message_stream("python", args)
    _message_stream("c", args)

    py = _message_stream("python", args)
    c = _message_stream("c", args)

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
