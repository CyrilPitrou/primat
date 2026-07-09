# -*- coding: utf-8 -*-
"""
cache_utils.py — fingerprinted self-validating cache files
============================================================================

Several expensive precomputations (n<->p weak rates, their finite-temperature
radiative corrections, the e+- thermodynamic tables) are written to plain-text
``np.savetxt`` files under ``rates/`` and reloaded on the next run instead of
being recomputed.  Historically these caches were trusted unconditionally:
whatever was on disk was used, even if the configuration that produced it
(neutrino-decoupling treatment, spectral distortions, sampling density, ...)
no longer matches the current run.  This silently makes flags such as
``spectral_distortions`` a no-op.

The fix is a *fingerprint*: a dict of every configuration entry that affects
the cached numbers, serialised as canonical (sorted-key, whitespace-free) JSON
and hashed with sha256 (truncated to 16 hex digits -- short enough to read,
long enough that two different configurations colliding by accident is
astronomically unlikely).  The hash and the JSON dict are written as
``#``-comment header lines of the cache file:

    # fingerprint_hash: a3f9c1b2e4d5f607
    # fingerprint: {"format_version":1,"sampling_nTOp_per_decade":80,...}

``np.loadtxt`` ignores ``#`` lines by default, so the data rows are unaffected.
The JSON line is for humans ("with which flags was this produced?"); only the
hash line is compared by the loader.  A cache file with no header (or an
unparsable one) is reported as having an unknown fingerprint -- the caller
decides whether that counts as a cache hit or a miss.
"""

import hashlib
import json
import os
import warnings

import numpy as np


def fingerprint_hash(fingerprint: dict) -> str:
    """Return the sha256 hash (first 16 hex digits) of a fingerprint dict.

    The dict is serialised to canonical JSON first (``sort_keys=True`` and no
    extra whitespace) so that the hash depends only on the *values*, not on
    the order in which the caller happened to build the dict.

    Args:
        fingerprint: dict of config values that determine a cache file's
            content (e.g. ``{"format_version": 1, "sampling_nTOp_per_decade": 80, ...}``).

    Returns:
        16-hex-character hash string, e.g. ``"a3f9c1b2e4d5f607"``.
    """
    blob = json.dumps(fingerprint, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def read_cache_fingerprint_hash(path: str):
    """Return the fingerprint hash stored in a cache file's header, or None.

    Reads only the leading ``#``-comment lines of `path`, looking for a line
    of the form ``# fingerprint_hash: <hash>``.  Stops at the first
    non-comment line (the data rows are never parsed).

    Args:
        path: path to a file previously written by
            :func:`write_cache_with_fingerprint`, or a legacy file with no
            header.

    Returns:
        The hash string if found, otherwise ``None`` -- which covers a
        missing file, a header-less legacy file, and a corrupt header.
    """
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                if line.startswith("# fingerprint_hash:"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        return None
    return None


def write_cache_with_fingerprint(path: str, fingerprint: dict, columns, col_header: str = "",
                                  provenance: str = None):
    """Write a ``np.savetxt`` cache file with a fingerprint header.

    Args:
        path: output file path; parent directory must already exist.
        fingerprint: dict to hash and embed verbatim as JSON (see
            :func:`fingerprint_hash`).
        columns: sequence of equal-length 1-D arrays, written column-wise
            (``np.column_stack(columns)``).
        col_header: optional human-readable column-name line, written before
            the fingerprint lines (e.g. ``"T[K] rate[1/s]"``).
        provenance: optional human-readable string recording which backend
            and algorithm computed this file (e.g.
            ``"backend=python algorithm=vegas"``), written as its own
            ``# provenance: ...`` header line *after* the fingerprint lines.
            Deliberately NOT part of `fingerprint` / the hash: this cache
            file is shared between backends (whichever computes it first,
            the other just reads it -- see weak_rates/cache.py), and for the
            deterministic (non-thermal) caches both backends always agree
            to machine precision anyway. For the thermal (CCRTh) cache,
            where both backends use independent Monte-Carlo estimates
            (vegas) with their own noise floor, this is purely informational
            provenance ("who produced the number on disk right now"), not a
            cache key -- it must never gate a cache hit/miss decision.

    Returns:
        True if the file was written, False if the write failed on an
        ``OSError`` (e.g. a read-only install). A failure is NOT fatal: the
        freshly computed in-memory values are perfectly valid, only the on-disk
        cache is skipped, so the caller continues and the next run just
        recomputes. A ``UserWarning`` is emitted that names the ``cache_dir``
        parameter as the remedy (redirect the whole cache tree to a writable
        directory — see :func:`cache_write_dir`). This graceful degradation is
        why the two writable cache trees can live inside the installed package
        without crashing read-only installs (B-1).

    Example:
        >>> write_cache_with_fingerprint(
        ...     "nTOp_frwrd.txt",
        ...     {"format_version": 1, "sampling_nTOp_per_decade": 80},
        ...     [T_all, frwrd], col_header="T[K] rate[1/s]")
    """
    fp_hash = fingerprint_hash(fingerprint)
    fp_json = json.dumps(fingerprint, sort_keys=True, separators=(",", ":"))
    header_lines = []
    if col_header:
        header_lines.append(col_header)
    header_lines.append("fingerprint_hash: " + fp_hash)
    header_lines.append("fingerprint: " + fp_json)
    if provenance:
        header_lines.append("provenance: " + provenance)
    # Write to a per-process temp file then atomically rename into place
    # (os.replace), so concurrent MC workers racing to populate a missing
    # cache never observe a partially-written file. The whole write (parent-dir
    # creation included, so a fresh cache_dir subtree materialises on demand) is
    # guarded: a read-only install must degrade to a warning, never a crash.
    tmp_path = f"{path}.tmp.{os.getpid()}"
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        np.savetxt(tmp_path, np.column_stack(columns), header="\n".join(header_lines))
        os.replace(tmp_path, path)
    except OSError as e:
        # Best-effort cleanup of a partial temp file (ignore if it too fails).
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        warnings.warn(
            f"could not write cache to {path}: {e}; results are unaffected, "
            "but the next run will recompute. Set the cache_dir parameter to "
            "redirect the cache to a writable directory.", stacklevel=2)
        return False
    return True


# ---------------------------------------------------------------------------
# Plasma+weak cache tree location (B-1). Both regenerable cache trees -- the
# n<->p weak rates (nTOp_<hash>.txt / nTOp_thermal_<hash>.txt) AND the plasma
# electron-thermo / QED-pressure tables -- live together under
# ``<data_dir>/cache_plasma_weak/{weak,plasma}/``. The ``cache_dir`` parameter
# redirects that whole tree to a writable directory for read-only installs,
# with OVERLAY semantics: reads try ``cache_dir`` first then fall back to the
# shipped package copy (so shipped caches are never shadowed), writes go only
# to ``cache_dir``. This mirrors config.py:resolve_rates_path / the
# user_nuclear_dir overlay exactly. The cache LOCATION is never part of any
# fingerprint -- it cannot affect the cached numbers.
# ---------------------------------------------------------------------------

def _cache_bases(cfg):
    """Overlay bases for the plasma+weak cache tree, first-wins on READ.

    ``cfg.cache_dir`` (the writable redirect) first if set, then the package
    ``<resolved_data_dir>/cache_plasma_weak`` -- always last so the shipped
    caches stay reachable even with cache_dir set. Cache LOCATION only:
    never part of any fingerprint.
    """
    bases = []
    if getattr(cfg, "cache_dir", None):
        bases.append(cfg.cache_dir)
    bases.append(os.path.join(cfg._resolved_data_dir, "cache_plasma_weak"))
    return bases


def cache_write_dir(cfg, subdir: str) -> str:
    """Directory a cache file is WRITTEN into: ``<writable base>/<subdir>``
    (the first base -- cfg.cache_dir if set, else the package root). Created
    on demand by the writer. ``subdir`` is ``"weak"`` or ``"plasma"``."""
    return os.path.join(_cache_bases(cfg)[0], subdir)


def resolve_cache_file(cfg, subdir: str, filename: str) -> str:
    """Resolve a cache file for READING through the overlay: the first
    existing ``<base>/<subdir>/<filename>``; if none exists, the write path
    (so a miss points where the file WILL be written)."""
    for base in _cache_bases(cfg):
        cand = os.path.join(base, subdir, filename)
        if os.path.exists(cand):
            return cand
    return os.path.join(cache_write_dir(cfg, subdir), filename)


# ---------------------------------------------------------------------------
# Writable weak-rate cache directory: inspection / cleanup (`primat
# --cache-info` / `--cache-clear`). Every new
# PRIMATConfig fingerprint run with spectral_distortions/incomplete_decoupling
# etc. drops another nTOp_<hash>.txt / nTOp_thermal_<hash>.txt file under the
# weak/ subdir of the cache tree; these are regenerable on demand (a fresh run
# just recomputes and re-caches them), so it is always safe to delete them.
# ---------------------------------------------------------------------------

def weak_cache_dir(cfg) -> str:    # write dir; used by the CLI cleanup helpers
    """Return the ``weak/`` WRITE directory (``cache_dir`` if set, else the
    package's ``cache_plasma_weak/weak``)."""
    return cache_write_dir(cfg, "weak")


def plasma_cache_dir(cfg) -> str:
    """Return the ``plasma/`` WRITE directory (``cache_dir`` if set, else the
    package's ``cache_plasma_weak/plasma``)."""
    return cache_write_dir(cfg, "plasma")


def list_weak_cache_files(cfg):
    """Return the sorted list of ``nTOp_*.txt`` cache file paths on disk.

    Iterates over the overlay bases (the ``cache_dir`` redirect, if set, and
    the shipped package tree) so ``--cache-info``/``--cache-clear`` see every
    reachable cache file, deduplicated by basename with the redirect winning.
    """
    seen = {}
    for base in _cache_bases(cfg):
        d = os.path.join(base, "weak")
        if not os.path.isdir(d):
            continue
        for name in os.listdir(d):
            if name.startswith("nTOp_") and name.endswith(".txt"):
                seen.setdefault(name, os.path.join(d, name))
    return sorted(seen.values())


def clear_weak_cache(cfg) -> int:
    """Delete every cached ``nTOp_*.txt`` file. Returns the count removed.

    The cache is purely an optimisation (every entry is reproducible from
    ``cfg`` by recomputing), so removing all of it is always safe -- the
    next run simply pays the one-time recompute cost again per
    configuration touched.
    """
    paths = list_weak_cache_files(cfg)
    removed = 0
    for path in paths:
        try:
            os.remove(path)
            removed += 1
        except OSError:
            # A shipped cache on a read-only install cannot be removed; skip it
            # rather than crash the --cache-clear command.
            pass
    return removed
