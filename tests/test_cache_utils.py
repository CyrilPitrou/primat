"""Direct unit tests for primat.cache_utils' fingerprinted-cache helpers.

These are exercised indirectly by every test that builds a PRIMATConfig and
solves (the n<->p weak-rate cache and the electron-thermo cache both use this
module), but a direct round-trip test pins the contract precisely:
fingerprint_hash is order-independent, read_cache_fingerprint_hash recovers
exactly the hash write_cache_with_fingerprint wrote, and a corrupted or
missing file is treated as "unknown fingerprint" (None) rather than raising.
The last test is a regression check for the atomic-write fix: the cache file
is written via a per-process temp file + os.replace so concurrent writers
racing on a missing cache cannot tear it (see cache_utils.py's module
docstring for the motivating incident).
"""
import os

import numpy as np
import pytest

from primat.cache_utils import (
    fingerprint_hash,
    read_cache_fingerprint_hash,
    write_cache_with_fingerprint,
)


def test_fingerprint_hash_is_order_independent():
    """The fingerprint hashes the dict's *content*, not its insertion order, so
    two configs that differ only in how their fields were assembled share a
    cache file. Also pins the 16-character digest length the filenames use."""
    h1 = fingerprint_hash({"a": 1, "b": 2})
    h2 = fingerprint_hash({"b": 2, "a": 1})
    assert h1 == h2
    assert len(h1) == 16


def test_fingerprint_hash_distinguishes_values():
    """A changed field value changes the hash -- the whole point of the scheme."""
    assert fingerprint_hash({"a": 1}) != fingerprint_hash({"a": 2})


def test_write_then_read_round_trip(tmp_path):
    """A written cache reads back with the same fingerprint and intact data rows
    (the fingerprint header lines must stay '#' comments np.loadtxt skips)."""
    path = str(tmp_path / "cache.txt")
    fp = {"format_version": 1, "sampling_nTOp_per_decade": 80}
    write_cache_with_fingerprint(path, fp, [np.array([1., 2., 3.]),
                                             np.array([4., 5., 6.])],
                                  col_header="T[K] rate[1/s]")
    assert read_cache_fingerprint_hash(path) == fingerprint_hash(fp)
    # Data rows must be intact (the fingerprint header lines are '#' comments
    # that np.loadtxt ignores).
    data = np.loadtxt(path)
    assert data.shape == (3, 2)
    assert data[:, 0].tolist() == [1., 2., 3.]


def test_read_missing_file_returns_none(tmp_path):
    """A missing cache file is reported as 'unknown fingerprint' (None), not an
    exception: a cache miss is normal control flow."""
    assert read_cache_fingerprint_hash(str(tmp_path / "does_not_exist.txt")) is None


def test_read_truncated_header_returns_none(tmp_path):
    """A header-less or corrupt file is also an unknown fingerprint, so a
    hand-edited or legacy cache degrades to a recompute rather than crashing."""
    path = tmp_path / "corrupt.txt"
    # No '# fingerprint_hash:' line at all -- a header-less legacy file.
    path.write_text("# just a comment\n1.0 2.0\n")
    assert read_cache_fingerprint_hash(str(path)) is None


def test_read_mismatched_hash_is_detected(tmp_path):
    """A stored hash differs from a different config's hash -- what makes the
    caller's 'stale cache' check work."""
    path = str(tmp_path / "cache.txt")
    write_cache_with_fingerprint(path, {"a": 1}, [np.array([1.])])
    stored = read_cache_fingerprint_hash(path)
    assert stored != fingerprint_hash({"a": 2})  # caller's mismatch check


def test_write_is_atomic_no_leftover_tmp_file(tmp_path):
    """Regression test for the os.replace() atomic-write fix.

    write_cache_with_fingerprint must write to f"{path}.tmp.{pid}" and
    os.replace() it into place, never np.savetxt directly to `path` -- so a
    reader can never observe a half-written file, and no stray .tmp.<pid>
    file is left behind afterwards.
    """
    path = str(tmp_path / "cache.txt")
    write_cache_with_fingerprint(path, {"a": 1}, [np.array([1., 2.])])
    assert os.path.exists(path)
    leftover = [f for f in os.listdir(tmp_path) if ".tmp." in f]
    assert leftover == []


def test_write_overwrites_existing_file_atomically(tmp_path):
    """Rewriting an existing cache replaces both its fingerprint and its data,
    with no residue from the previous contents."""
    path = str(tmp_path / "cache.txt")
    write_cache_with_fingerprint(path, {"a": 1}, [np.array([1.])])
    write_cache_with_fingerprint(path, {"a": 2}, [np.array([99.])])
    assert read_cache_fingerprint_hash(path) == fingerprint_hash({"a": 2})
    assert np.loadtxt(path).item() == pytest.approx(99.)


# ---------------------------------------------------------------------------
# cache_dir redirect + cache_plasma_weak/ overlay. The two writable
# cache trees (weak/ + plasma/) live under primat/data/cache_plasma_weak/,
# and the cache_dir parameter redirects WRITES elsewhere while still READING
# the shipped caches through an overlay (never shadowing them). A failed
# cache write degrades to a UserWarning (naming cache_dir), never a crash.
# ---------------------------------------------------------------------------

def test_cache_dir_param_redirects_writes(tmp_path):
    """cache_dir=<dir> makes the WRITE dirs live under it (<dir>/weak,
    <dir>/plasma); unset falls back to <data_dir>/cache_plasma_weak/*."""
    from primat.config import PRIMATConfig
    from primat.cache_utils import cache_write_dir
    cfg = PRIMATConfig({"cache_dir": str(tmp_path)})
    assert cache_write_dir(cfg, "weak")   == os.path.join(str(tmp_path), "weak")
    assert cache_write_dir(cfg, "plasma") == os.path.join(str(tmp_path), "plasma")
    cfg_default = PRIMATConfig({})
    assert cache_write_dir(cfg_default, "weak").endswith(
        os.path.join("data", "cache_plasma_weak", "weak"))


def test_cache_dir_overlay_still_reads_shipped_caches(tmp_path):
    """OVERLAY semantics: with cache_dir set but a file absent there, the
    resolver falls back to the shipped <data_dir>/cache_plasma_weak/<sub>/
    copy (so shipped caches are never shadowed); a file present in cache_dir
    wins over the shipped one."""
    from primat.config import PRIMATConfig
    from primat.cache_utils import resolve_cache_file
    cfg = PRIMATConfig({"cache_dir": str(tmp_path)})
    # Absent in cache_dir -> resolves to a shipped nTOp_*.txt (there is at
    # least one shipped weak cache): the returned path must NOT be under
    # tmp_path and must exist.
    import glob, primat
    shipped = os.path.join(os.path.dirname(primat.__file__),
                           "data", "cache_plasma_weak", "weak")
    name = os.path.basename(sorted(glob.glob(os.path.join(shipped, "nTOp_*.txt")))[0])
    got = resolve_cache_file(cfg, "weak", name)
    assert got == os.path.join(shipped, name) and os.path.exists(got)
    # Present in cache_dir -> that copy wins.
    (tmp_path / "weak").mkdir()
    (tmp_path / "weak" / name).write_text("# local\n")
    assert resolve_cache_file(cfg, "weak", name) == os.path.join(
        str(tmp_path), "weak", name)


def test_cache_write_failure_warns_instead_of_raising(tmp_path):
    """An unwritable target directory must degrade to a UserWarning, not a
    PermissionError crash -- the freshly computed in-memory values are valid.
    The warning must both explain AND point at the cache_dir remedy (author
    decision 2026-07-09: read-only installs are what cache_dir exists for)."""
    import warnings
    # Force a portable write failure by pointing the cache at a path whose
    # parent is a regular file: os.makedirs() inside write_cache_with_fingerprint
    # then raises NotADirectoryError (an OSError) on every platform. A read-only
    # directory via chmod(0o500) is NOT portable -- on Windows os.chmod cannot
    # clear a directory's write permission, so the write would spuriously succeed.
    blocker = tmp_path / "blocker"
    blocker.write_text("not a directory")
    target = blocker / "sub" / "nTOp_test.txt"
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ok = write_cache_with_fingerprint(
            str(target), {"field": 1.0},
            [np.ones(3), np.zeros(3)], col_header="a b")
    assert ok is False
    assert any("could not write" in str(x.message)
               and "cache_dir" in str(x.message) for x in w)


def test_cache_dir_not_in_weak_rate_fingerprint():
    """Cache LOCATION must not invalidate caches (it cannot affect numbers)."""
    from primat.weak_rates.cache import WEAK_RATE_BG_FIELDS
    assert "cache_dir" not in WEAK_RATE_BG_FIELDS


def test_shipped_data_uses_cache_plasma_weak_layout():
    """The relocated cache tree exists and the old top-level dirs are gone
    (hard cutover -- no legacy fallback is supported)."""
    import primat
    data = os.path.join(os.path.dirname(primat.__file__), "data")
    assert os.path.isdir(os.path.join(data, "cache_plasma_weak", "weak"))
    assert os.path.isdir(os.path.join(data, "cache_plasma_weak", "plasma"))
    assert not os.path.exists(os.path.join(data, "weak"))
    assert not os.path.exists(os.path.join(data, "plasma"))


# --- numpy scalars in a fingerprint (must not abort a run, must not change
# --- any existing hash) ---

def test_fingerprint_hash_accepts_numpy_scalars():
    """GOAL: a numpy scalar config value must hash, not crash.

    Numpy scalars arrive naturally from a parameter scan built with
    np.arange/np.linspace, or from an external driver (the Cobaya wrapper)
    indexing a sampled array. np.float64 always worked because it subclasses
    float, but np.int64/np.float32/np.bool_ raised
    "TypeError: Object of type int64 is not JSON serializable" from deep inside
    the weak-rate cache, aborting the whole run.
    """
    for value, plain in ((np.int64(80), 80),
                         (np.float64(0.5), 0.5),
                         (np.bool_(True), True)):
        # Hash-preserving: unwrapping must not invalidate any cache file that a
        # plain-Python value already produced.
        assert fingerprint_hash({"x": value}) == fingerprint_hash({"x": plain})


def test_fingerprint_hash_still_rejects_unserialisable():
    """Only numpy scalars are unwrapped; a genuinely unserialisable value must
    still raise rather than hash to something arbitrary (a silently-wrong
    fingerprint would be far worse than a crash)."""
    with pytest.raises(TypeError):
        fingerprint_hash({"x": object()})


def test_write_accepts_bare_filename(tmp_path, monkeypatch):
    """GOAL: a cache path with no directory part must be writable.

    os.makedirs("") raises FileNotFoundError, which the writer's `except
    OSError` turned into a spurious "could not write cache" warning for a
    perfectly writable target -- including for the form used in
    write_cache_with_fingerprint's own docstring example.
    """
    monkeypatch.chdir(tmp_path)
    fp = {"format_version": 1}
    assert write_cache_with_fingerprint(
        "nTOp_frwrd.txt", fp, [np.array([1., 2.]), np.array([3., 4.])])
    assert read_cache_fingerprint_hash("nTOp_frwrd.txt") == fingerprint_hash(fp)
