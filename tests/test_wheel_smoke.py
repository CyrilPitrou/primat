"""
"pip install" smoke test: build a wheel, install it in a
clean virtual environment, and run a small BBN solve.

Why this test exists
---------------------
The editable install used during development (``pip install -e .`` / running
straight from a git checkout) resolves ``rates/`` relative to the source tree
no matter how the path is computed, so a bug in the package-data path
resolution (e.g. the stale ``../rates`` left over from the package
reorganisation) is invisible there.  It only shows up once the
package is installed as a *wheel* into ``site-packages`` -- a different
directory layout entirely.  Building the wheel also exercises the
``[tool.setuptools.package-data]`` declaration in ``pyproject.toml``: if a
required file under ``rates/`` were ever excluded, the import would still
succeed but ``PRIMAT(...).solve()`` would fail with a ``FileNotFoundError``
deep inside ``primat.network_data``/``primat.weak_rates``.

The venv is created with ``--system-site-packages`` so the already-installed
numpy/scipy/joblib (and any optional numba/vegas) are reused --
this test checks the *primat* packaging, not whether its dependencies can
be downloaded.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import venv
import zipfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _venv_python(venv_dir):
    """Interpreter inside a ``venv.create``d directory, on any platform.

    Windows puts it in ``Scripts/python.exe``, POSIX in ``bin/python``.
    """
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


@pytest.mark.slow
@pytest.mark.wheel
def test_wheel_install_smoke_solve():
    """Build a wheel, pip-install it in a clean venv, and run a small solve.

    Steps:
      1. ``pip wheel`` the repo root into a temporary directory (using the
         setuptools build backend already installed in this environment, so
         no network access is required).
      2. Create a fresh venv (``--system-site-packages`` to reuse numpy/scipy/
         joblib already present) and ``pip install --no-deps`` the wheel.
      3. In that venv, run a default-configuration small-network solve and
         check YP/D-H against the same loose tolerances
         ``tests/test_regression.py``'s default-precision checks use.

    A failure here most likely means ``rates/`` data files are missing from
    the wheel, or a path is computed relative to the source tree instead of
    the installed package (``primat.config.PRIMATConfig.data_dir``).
    """
    with tempfile.TemporaryDirectory(prefix="primat_wheel_") as tmp:
        tmp_path = Path(tmp)
        wheel_dir = tmp_path / "wheel"
        venv_dir  = tmp_path / "venv"

        # ------------------------------------------------------------
        # 1. Build the wheel (no build isolation: reuse the setuptools
        #    already installed here, avoiding any network access for the
        #    build itself).
        #
        #    ``--no-build-isolation`` requires setuptools to be importable in
        #    THIS interpreter, but PEP 668 / Python >= 3.12 venvs no longer
        #    ship setuptools by default (only the dev env has it incidentally),
        #    so the nightly CI leg failed with pip exit code 2. Guarantee the
        #    build requirement first -- a no-op (no network) when setuptools is
        #    already present, a one-time fetch on a bare 3.12 environment.
        # ------------------------------------------------------------
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "-q", "setuptools>=61"],
            check=True,
        )
        subprocess.run(
            [sys.executable, "-m", "pip", "wheel", str(REPO_ROOT),
             "-w", str(wheel_dir), "--no-deps", "--no-build-isolation", "-q"],
            check=True,
        )
        wheels = list(wheel_dir.glob("*.whl"))
        assert len(wheels) == 1, f"expected exactly one wheel, got {wheels}"

        # ------------------------------------------------------------
        # 2. Fresh venv + install the wheel (reusing system site-packages
        #    for numpy/scipy/joblib/numba/...).
        # ------------------------------------------------------------
        venv.create(venv_dir, with_pip=True, system_site_packages=True)
        venv_python = _venv_python(venv_dir)
        subprocess.run(
            [str(venv_python), "-m", "pip", "install", "--no-deps", "-q",
             str(wheels[0])],
            check=True,
        )

        # ------------------------------------------------------------
        # 3. Smoke solve: default config, small network.  save_nTOp
        #    defaults to False, so this does not write into site-packages.
        # ------------------------------------------------------------
        smoke_script = (
            # Ignore the development checkout's editable import hook, which
            # would otherwise serve its own primat submodules to this venv --
            # see test_install_falls_back_to_python_... below.
            "import sys\n"
            "sys.meta_path = [f for f in sys.meta_path if not\n"
            "                 getattr(f, '__module__', '').startswith('__editable__')]\n"
            "from primat import PRIMAT\n"
            "import primat\n"
            "assert 'site-packages' in primat.__file__, primat.__file__\n"
            "r = PRIMAT({'network': 'small', 'verbose': False, 'debug': False}).solve()\n"
            "print(r['YPBBN'], r['DoH'])\n"
        )
        result = subprocess.run(
            [str(venv_python), "-c", smoke_script],
            # cwd outside the checkout: `python -c` puts the working directory
            # first on sys.path, so running from the repo root would import the
            # source tree and never touch the wheel this test just installed.
            cwd=str(tmp_path), check=True, capture_output=True, text=True,
        )

    # Same loose tolerances as tests/test_regression.py::test_small_network_*
    yp_str, doh_str = result.stdout.split()
    yp, doh = float(yp_str), float(doh_str)
    assert yp  == pytest.approx(0.2469983, abs=1e-4)
    assert doh == pytest.approx(2.43490e-5, rel=2e-3)


# Script run in a subprocess for test_core_runs_without_plotly_or_joblib
# (below).  It installs a meta-path finder that makes ``import plotly`` and
# ``import joblib`` fail, *then* imports primat and exercises the two paths a
# lean core install must support without those two now-optional dependencies:
# a single ``run_bbn`` solve, and a serial
# (``n_jobs=1``) Monte-Carlo run.  A subprocess is used so the import blocker
# and any already-imported plotly/joblib in the test session cannot interfere.
_NO_OPTIONAL_DEPS_SCRIPT = r"""
import sys
import importlib.abc


class _Blocker(importlib.abc.MetaPathFinder):
    '''Refuse to import the named top-level packages (and their submodules).'''
    def __init__(self, blocked):
        self.blocked = set(blocked)

    def find_spec(self, name, path, target=None):
        if name.split('.')[0] in self.blocked:
            raise ImportError(f"{name} is blocked for this test")
        return None


# Drop anything already imported so the blocker actually bites, then install it.
for _m in list(sys.modules):
    if _m.split('.')[0] in {"plotly", "joblib"}:
        del sys.modules[_m]
sys.meta_path.insert(0, _Blocker({"plotly", "joblib"}))

# Importing primat.backend and running a solve must touch neither package.
import primat.backend as b
r = b.run_bbn({"network": "small", "verbose": False, "debug": False},
              force_backend="python")
assert abs(r["YPBBN"] - 0.247) < 2e-3, r["YPBBN"]

# Serial Monte-Carlo (n_jobs=1) must run without importing joblib.
mc = b.run_mc(2, ["YPBBN"], params={"network": "small"},
              force_backend="python", n_jobs=1, seed=0)
assert mc["YPBBN"].mean > 0.0

# Sanity: the blocker really is active (a real import would have raised).
assert "plotly" not in sys.modules and "joblib" not in sys.modules
print("OK")
"""


@pytest.mark.slow
@pytest.mark.solve
def test_core_runs_without_plotly_or_joblib():
    """A lean core install (no plotly, no joblib) can still ``run_bbn`` and run
    serial Monte-Carlo.

    plotly (GUI figures) and joblib (parallel MC) were moved out of the hard
    dependencies into extras.  This pins that promise: with both packages made
    un-importable, ``primat.backend.run_bbn`` and a serial ``run_mc(n_jobs=1)``
    on the pure-Python backend must both succeed -- proving neither the
    ``import primat`` path nor the core solve/serial-MC path secretly depends
    on the two now-optional packages.
    """
    result = subprocess.run(
        [sys.executable, "-c", _NO_OPTIONAL_DEPS_SCRIPT],
        cwd=str(REPO_ROOT), capture_output=True, text=True,
    )
    assert result.returncode == 0, (
        "core solve / serial MC failed without plotly+joblib:\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    assert "OK" in result.stdout


# ---------------------------------------------------------------------------
# The documented degradation path: no working compiler -> pure Python
# ---------------------------------------------------------------------------

# Files a wheel build needs. The whole checkout is ~500 MB (git history, the
# built docs); these are the build inputs, ~26 MB, which keeps the copy below
# affordable.
_BUILD_INPUTS = ("pyproject.toml", "setup.py", "MANIFEST.in", "README.md",
                 "LICENSE", "primat", "primat-c")


@pytest.mark.slow
@pytest.mark.wheel
def test_install_falls_back_to_python_when_the_extension_cannot_build():
    """`pip install` must succeed with a broken C toolchain, minus the extension.

    GOAL: pin the promise setup.py's ``optional_build_ext`` makes and
    README/docs repeat -- that a missing or broken compiler degrades to the
    pure-Python backend instead of failing the install. It rests on a broad
    ``except Exception`` around setuptools' ``build_ext``, which nothing
    exercised, so a change in how setuptools reports build failures would turn
    the graceful degradation into a hard install error unnoticed.

    The failure is forced by corrupting the wrapper's C source in a throwaway
    copy of the build inputs -- a real compiler error, on whichever compiler
    the platform uses, rather than a simulated one.
    """
    with tempfile.TemporaryDirectory(prefix="primat_nocc_") as tmp:
        tmp_path = Path(tmp)
        src = tmp_path / "src"
        src.mkdir()
        for name in _BUILD_INPUTS:
            origin = REPO_ROOT / name
            if not origin.exists():
                continue
            if origin.is_dir():
                # Skip build outputs: an editable checkout has the compiled
                # extension sitting in primat/, and copying it would leave a
                # working _primat_c next to the deliberately broken source.
                shutil.copytree(origin, src / name, symlinks=True,
                                ignore=shutil.ignore_patterns(
                                    "*.so", "*.pyd", "*.dll", "*.o",
                                    "__pycache__", "build", "*.egg-info"))
            else:
                shutil.copy(origin, src / name)

        wrapper = src / "primat" / "_primat_c_src" / "_wrapper.c"
        wrapper.write_text("#error primat test: forced extension build failure\n")

        wheel_dir = tmp_path / "wheel"
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "-q", "setuptools>=61"],
            check=True,
        )
        build = subprocess.run(
            [sys.executable, "-m", "pip", "wheel", str(src), "-w", str(wheel_dir),
             "--no-deps", "--no-build-isolation"],
            capture_output=True, text=True,
        )
        # The point of the test: a compile error is not an install error.
        assert build.returncode == 0, build.stdout + build.stderr
        wheels = list(wheel_dir.glob("*.whl"))
        assert len(wheels) == 1, f"expected exactly one wheel, got {wheels}"

        venv_dir = tmp_path / "venv"
        venv.create(venv_dir, with_pip=True, system_site_packages=True)
        venv_python = _venv_python(venv_dir)
        subprocess.run(
            [str(venv_python), "-m", "pip", "install", "--no-deps", "-q",
             str(wheels[0])],
            check=True,
        )

        # The wheel itself is the primary evidence: no compiled extension was
        # packaged, so the install cannot have one.
        with zipfile.ZipFile(wheels[0]) as zf:
            compiled = [n for n in zf.namelist() if n.endswith((".so", ".pyd", ".dll"))]
        assert compiled == [], f"wheel carries a compiled extension: {compiled}"

        # And it still solves. HAS_C_BACKEND is asserted only when `primat`
        # actually resolves inside the venv: the venv reuses system site
        # packages for numpy/scipy, which on a development machine also exposes
        # the editable checkout (extension included) through its import hook.
        probe = (
            # The venv reuses system site packages (for numpy/scipy), which on
            # a development machine also carries the editable checkout's
            # import hook -- and that hook resolves primat._primat_c to the
            # compiled extension in the source tree even when primat itself
            # comes from the venv. Drop it, or the assertion below tests the
            # developer's build rather than the wheel just installed.
            "import sys\n"
            "sys.meta_path = [f for f in sys.meta_path if not\n"
            "                 getattr(f, '__module__', '').startswith('__editable__')]\n"
            "import primat, primat.backend as b\n"
            "print(primat.__file__)\n"
            "print(b.HAS_C_BACKEND)\n"
            "r = b.run_bbn({'network': 'small'})\n"
            "print(r['YPBBN'], r['DoH'])\n"
        )
        result = subprocess.run(
            [str(venv_python), "-c", probe],
            cwd=str(tmp_path),      # see test_wheel_install_smoke_solve
            check=True, capture_output=True, text=True,
        )

    origin, has_c, values = result.stdout.strip().split("\n")
    assert "site-packages" in origin, f"probe imported {origin}, not the install"
    assert has_c == "False", "extension present in a wheel built without it"

    yp, doh = (float(v) for v in values.split())
    assert yp  == pytest.approx(0.2469983, abs=1e-4)
    assert doh == pytest.approx(2.43490e-5, rel=2e-3)
