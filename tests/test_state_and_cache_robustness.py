"""Behaviour when many runs share one process, and when a cache file is damaged.

These guard three promises that only break off the happy path:

* a damaged cache file is *recovered from* (electron thermo) or *reported*
  (n<->p weak rates) rather than trusted, identically on both backends;
* a cache that could not be written completely is never installed;
* the process-global state primat keeps -- the numba kernel rebinding, the
  reaction catalog, the GUI's network-label fallback -- stays bounded and
  safe when several configurations are alive at once.
"""
import glob
import os
import subprocess
import sys
import threading
import warnings

import pytest

from primat.backend import HAS_C_BACKEND, run_bbn

pytestmark = [pytest.mark.slow, pytest.mark.solve]

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "primat", "data")


def _cache_copy(tmp_path, subdir, mutate):
    """Build a cache_dir holding every shipped file of ``subdir``, mutated.

    Args:
        tmp_path : pytest tmp_path.
        subdir   : "weak" or "plasma".
        mutate   : text -> text, applied to each file's contents.

    Returns:
        str path usable as the ``cache_dir`` parameter.
    """
    d = tmp_path / "cache"
    (d / subdir).mkdir(parents=True, exist_ok=True)
    pattern = {"weak": "weak/nTOp_*.txt",
               "plasma": "plasma/electron_thermo_*.txt"}[subdir]
    for f in glob.glob(os.path.join(DATA_DIR, "cache_plasma_weak", pattern)):
        (d / subdir / os.path.basename(f)).write_text(mutate(open(f).read()))
    return str(d)


def _half(text):
    return text[: len(text) // 2]


# ---------------------------------------------------------------------------
# A damaged cache must not be trusted, and must not take the process with it
# ---------------------------------------------------------------------------

def test_truncated_electron_thermo_cache_is_recomputed_not_trusted(tmp_path):
    """A half-written electron-thermo cache is recovered from, not obeyed.

    The file keeps its fingerprint header, so it is selected and then fails to
    parse. Both backends must warn and recompute; the C backend used to free
    the fingerprint string twice on exactly this path and abort the process,
    which through the extension killed the caller's interpreter.
    """
    cache = _cache_copy(tmp_path, "plasma", _half)
    for backend in (["c", "python"] if HAS_C_BACKEND else ["python"]):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = run_bbn({"cache_dir": cache}, force_backend=backend)
        assert res["DoH"] == pytest.approx(2.436e-05, rel=1e-3), backend


def test_truncated_weak_rate_cache_reports_the_file_and_row_on_both_backends(tmp_path):
    """A truncated n<->p cache raises one sentence, the same on both backends.

    numpy's own message ("the number of columns changed from 3 to 1 at row
    157") names neither the file nor a cache, and points at a `usecols`
    argument the caller never passed.
    """
    cache = _cache_copy(tmp_path, "weak", _half)
    with pytest.raises(ValueError) as exc:
        run_bbn({"cache_dir": cache}, force_backend="python")
    msg = str(exc.value)
    assert "expected 3 columns, found" in msg
    assert "nTOp_" in msg and cache in msg
    assert "usecols" not in msg

    if not HAS_C_BACKEND:
        return
    with pytest.raises(Exception) as c_exc:
        run_bbn({"cache_dir": cache}, force_backend="c")
    assert msg in str(c_exc.value), "the two backends must report it identically"


@pytest.mark.skipif(sys.platform.startswith("win"), reason="needs ulimit -f")
def test_a_cache_that_cannot_be_written_completely_is_not_installed(tmp_path):
    """Running out of disk space must not leave a truncated cache behind.

    The writer renames a temporary file into place; without checking that the
    rows actually reached the disk, a short write installs a file that keeps
    its fingerprint header and loses its data rows, which every later run then
    trusts.
    """
    binary = os.path.join(os.path.dirname(__file__), "..", "primat-c", "build", "primat-c")
    if not os.path.exists(binary):
        pytest.skip("primat-c CLI not built")
    cache = tmp_path / "cache"
    # `trap '' XFSZ` because the file-size limit delivers SIGXFSZ where a real
    # full disk returns an error from write(2); ignoring the signal turns the
    # limit into the ENOSPC-shaped failure this guards.
    cmd = (f"ulimit -f 20; trap '' XFSZ; exec {binary} --data_dir {DATA_DIR} "
           f"--set cache_dir={cache} --set me=0.510996 "
           f"--set recompute_electron_thermo=True --set thermal_corrections=False")
    proc = subprocess.run(["/bin/sh", "-c", cmd], capture_output=True, text=True)
    installed = [p for p in glob.glob(str(cache / "**" / "*"), recursive=True)
                 if os.path.isfile(p)]
    assert installed == [], f"a short write installed {installed}"
    assert "could not write cache" in proc.stdout + proc.stderr


# ---------------------------------------------------------------------------
# Process-global state
# ---------------------------------------------------------------------------

def test_numba_rebinding_survives_concurrent_flipping():
    """Two configurations flipping ``use_numba`` must not kill a thread.

    The rebinding is process-wide; setting its "already done" flag before the
    names were actually rebound let a second thread re-wrap an already-jitted
    function, which numba rejects with a TypeError naming nothing about primat.
    """
    from primat import plasma
    from primat.weak_rates import integrands

    errors = []

    def flip(value, n):
        for _ in range(n):
            try:
                integrands._setup_fd_impls(value)
                plasma._setup_electron_integrands(value)
            except Exception as exc:                      # pragma: no cover
                errors.append(exc)

    threads = [threading.Thread(target=flip, args=(v, 40))
               for v in (True, False) for _ in range(3)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert errors == [], errors
    # Leave the process on the default so later tests are unaffected.
    integrands._setup_fd_impls(True)
    plasma._setup_electron_integrands(True)


def test_coexisting_configurations_reproduce_their_solo_results():
    """Configurations alive at once must each give what they give alone.

    Every instance shares the module-level plasma/weak-rate machinery, so a
    configuration that leaked state into another would show up here as a
    digit that depends on what else was built.
    """
    cases = [{}, {"gA": 1.2700}, {"network": "small_parthenope"},
             {"Omegabh2": 0.0200}]
    alone = [run_bbn(dict(c), force_backend="python")["DoH"] for c in cases]
    built = [dict(c) for c in cases]
    together = [run_bbn(c, force_backend="python")["DoH"] for c in reversed(built)]
    assert list(reversed(together)) == alone


def test_process_global_caches_are_bounded():
    """The two process-global caches must evict rather than grow for ever.

    Both outlive every run in the process: a server handing out temporary data
    trees or custom-network names would otherwise retain each one it ever saw.
    """
    from primat.network_data import _reaction_catalog
    assert _reaction_catalog.cache_info().maxsize is not None

    pytest.importorskip("streamlit")
    from primat.gui import params_form

    params_form._network_label_cache.clear()
    for i in range(params_form._NETWORK_LABEL_CACHE_MAX + 20):
        params_form._remember_network_label(f"net{i}", f"net{i} (3)")
    cache = params_form._network_label_cache
    assert len(cache) == params_form._NETWORK_LABEL_CACHE_MAX
    assert f"net{params_form._NETWORK_LABEL_CACHE_MAX + 19}" in cache
    assert "net0" not in cache
    cache.clear()
