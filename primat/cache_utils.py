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

import functools
import hashlib
import json
import os
import warnings

import numpy as np

from .constants import CONST, Constants


def _json_scalar(obj):
    """``json.dumps(default=...)`` hook: unwrap a numpy scalar to its Python
    equivalent, and let anything else raise as before.

    A config value can perfectly well arrive as a numpy scalar -- a parameter
    scan built with ``np.arange``/``np.linspace``, or an external driver such as
    the Cobaya wrapper handing over an element of a sampled array.  ``np.float64``
    happened to survive ``json.dumps`` (it subclasses ``float``) while
    ``np.int64``/``np.float32``/``np.bool_`` did not, so an ``np.int64`` for e.g.
    ``sampling_nTOp_per_decade`` aborted the whole run with an opaque
    ``TypeError: Object of type int64 is not JSON serializable`` raised from deep
    inside the weak-rate cache.

    ``.item()`` yields the exact Python scalar (``np.int64(80) -> 80``), whose
    canonical JSON is byte-identical to what a plain ``80`` produces -- so
    hardening this is hash-preserving: no existing cache file is invalidated.
    """
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _canonical_json(fingerprint: dict) -> str:
    """Canonical JSON of a fingerprint dict: ``sort_keys=True``, no padding
    whitespace, numpy scalars unwrapped (:func:`_json_scalar`).

    Single source of truth for both the hashed blob (:func:`fingerprint_hash`)
    and the human-readable ``# fingerprint:`` header line
    (:func:`write_cache_with_fingerprint`), so the two can never disagree.
    """
    return json.dumps(fingerprint, sort_keys=True, separators=(",", ":"),
                      default=_json_scalar)


def fingerprint_hash(fingerprint: dict) -> str:
    """Return the sha256 hash (first 16 hex digits) of a fingerprint dict.

    The dict is serialised to canonical JSON first (``sort_keys=True`` and no
    extra whitespace) so that the hash depends only on the *values*, not on
    the order in which the caller happened to build the dict.  Numpy scalars
    are unwrapped to their Python equivalents on the way (see
    :func:`_json_scalar`), so ``np.int64(80)`` and ``80`` hash identically.

    Args:
        fingerprint: dict of config values that determine a cache file's
            content (e.g. ``{"format_version": 1, "sampling_nTOp_per_decade": 80, ...}``).

    Returns:
        16-hex-character hash string, e.g. ``"a3f9c1b2e4d5f607"``.
    """
    return hashlib.sha256(_canonical_json(fingerprint).encode("utf-8")).hexdigest()[:16]


# Physical constants each fingerprinted cache actually reads -- the only ones
# whose value can change its numbers, and so the only ones it is keyed on.
# The ten frozen constants (kB, MeV, hbar, clight, ... ) are absent by
# construction: PRIMATConfig rejects an override of them, so they cannot vary
# within a process even though the integrands read them. Editing one in the
# source is a code change like any other, and is covered by bumping the
# affected caches' FORMAT_VERSION.
# tests/test_cache_constant_deps.py perturbs every settable constant and
# asserts each cache's data moves iff it is listed here.
CACHE_CONSTANTS = {
    # Stored in units of 1/tau_n, i.e. already divided by ComputeFn, which is
    # what pulls in gA and the anomalous moments alongside the integrands' own
    # me/mn/mp/alphaem/radproton.
    "weak":            ("alphaem", "gA", "kappa_n", "kappa_p",
                        "me", "mn", "mp", "radproton"),
    # Stores the RAW L_CCRTh (the 1/Fn division happens at point of use), so
    # only the integrands' own constants appear: the O(alphaem) prefactor,
    # Q = mn - mp, and FermiCoulomb's me/radproton.
    "thermal":         ("alphaem", "me", "mn", "mp", "radproton"),
    # e+- integrands and the grid's lower edge me/30.
    "electron_thermo": ("me",),
    # delta_P_a / delta_P_e3 integrands.
    "qed":             ("alphaem", "me"),
}


@functools.lru_cache(maxsize=64)
def _constants_hash_cached(consts: Constants, cache: str) -> str:
    return fingerprint_hash({k: getattr(consts, k) for k in CACHE_CONSTANTS[cache]})


def constants_hash(cache: str, cfg=None) -> str:
    """Return the 16-hex-digit hash of the constants one cache reads.

    Each of the four fingerprinted caches embeds this in place of a hash of
    the whole ``Constants`` struct, so a constant the cache does not read
    cannot invalidate it -- ``--T0CMB`` no longer costs a two-minute CCRTh
    recompute.  The declared sets are :data:`CACHE_CONSTANTS`.  Only dataclass
    *fields* are hashed, never the ``@property``-derived quantities: they add
    no information, and excluding them keeps the hash free of any float a C
    compiler might contract differently from CPython.

    The C backend computes the identical value per cache
    (``cpr_constants_hash``, ``primat-c/src/cache.c``), so both backends key a
    shared cache file the same way; ``tests/test_cache_parity.py`` asserts it.

    Args:
        cache: which cache's set to hash -- a key of :data:`CACHE_CONSTANTS`
            (``"weak"``, ``"thermal"``, ``"electron_thermo"``, ``"qed"``).
        cfg: a ``PRIMATConfig`` (or anything exposing ``.constants``).
            ``None`` uses the all-defaults :data:`primat.constants.CONST`.

    Returns:
        16-hex-character hash string, e.g. ``"6e0c1c4c95a2b6b0"``.

    Example:
        >>> constants_hash("qed", cfg)         # doctest: +SKIP
        '6e0c1c4c95a2b6b0'
    """
    consts = CONST if cfg is None else cfg.constants
    return _constants_hash_cached(consts, cache)


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
                                  provenance: str = None, fmt: str = None):
    """Write a ``np.savetxt`` cache file with a fingerprint header.

    Args:
        path: output file path; any missing parent directory is created on
            demand (a bare filename writes into the current working directory).
        fingerprint: dict to hash and embed verbatim as JSON (see
            :func:`fingerprint_hash`).
        columns: sequence of equal-length 1-D arrays, written column-wise
            (``np.column_stack(columns)``).
        col_header: optional human-readable column-name line, written before
            the fingerprint lines (e.g. ``"T[K] rate[1/s]"``). May be several
            lines separated by ``"\\n"`` -- ``np.savetxt`` prefixes each with
            ``"# "`` -- which is how the QED pressure tables keep their
            three-line physics provenance block above the fingerprint.
        fmt: optional ``np.savetxt`` format string. ``None`` (default) uses
            numpy's own default (``"%.18e"``), which is what every cache
            written from scratch here uses. It is overridden only by the QED
            pressure tables, whose shipped files predate this fingerprinting
            and were written with ``"%.6E"``: passing that format lets them
            gain a fingerprint header with their data rows *byte-identical*,
            so adding the header is provably a header-only change (see
            :func:`primat.qed_pressure.save_qed_tables`).
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
        without crashing read-only installs.

    Example:
        >>> write_cache_with_fingerprint(
        ...     "nTOp_frwrd.txt",
        ...     {"format_version": 1, "sampling_nTOp_per_decade": 80},
        ...     [T_all, frwrd], col_header="T[K] rate[1/s]")
    """
    fp_hash = fingerprint_hash(fingerprint)
    fp_json = _canonical_json(fingerprint)
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
        # dirname is "" for a bare filename (a cache written into the current
        # working directory); os.makedirs("") raises FileNotFoundError, which
        # the except below would have turned into a spurious "could not write
        # cache" for a perfectly writable target -- so only create a parent
        # directory when there is one to create.
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        savetxt_kw = {} if fmt is None else {"fmt": fmt}
        np.savetxt(tmp_path, np.column_stack(columns), header="\n".join(header_lines),
                   **savetxt_kw)
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
# Plasma+weak cache tree location. Both regenerable cache trees -- the
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
# Writable cache directories: inspection / cleanup (`primat --cache-info` /
# `--cache-clear`). Every new PRIMATConfig fingerprint run with
# spectral_distortions/incomplete_decoupling etc. drops another
# nTOp_<hash>.txt / nTOp_thermal_<hash>.txt file under the weak/ subdir, AND --
# since the electron-thermo cache became hash-named too -- another
# electron_thermo_<hash>.txt under plasma/. Both trees are therefore swept.
#
# Everything these helpers touch is regenerable on demand (a fresh run just
# recomputes and re-caches), so it is always safe to delete -- including the
# copies that ship with the package, which are a precomputed convenience rather
# than irreplaceable data (see clear_cache for what a recompute does and does
# not reproduce).
#
# The QED pressure tables (QED_pressure_correction_e{2,3}.txt) are deliberately
# NOT swept: they keep fixed filenames and are written only when
# recompute_qed_corrections asks for it, so unlike the two hash-named families
# they cannot proliferate. Nothing accumulates, so there is nothing to clean.
# ---------------------------------------------------------------------------

# Cache-file basename prefixes swept by list_cache_files/clear_cache, per
# subdirectory of the cache tree. Prefix + ".txt" suffix, matching how each
# family is named at its write site (weak_rates/api.py, plasma.py).
_CACHE_PREFIXES = {
    "weak":   ("nTOp_",),            # covers nTOp_thermal_<hash>.txt as well
    "plasma": ("electron_thermo_",),
}

def weak_cache_dir(cfg) -> str:    # write dir; used by the CLI cleanup helpers
    """Return the ``weak/`` WRITE directory (``cache_dir`` if set, else the
    package's ``cache_plasma_weak/weak``)."""
    return cache_write_dir(cfg, "weak")


def plasma_cache_dir(cfg) -> str:
    """Return the ``plasma/`` WRITE directory (``cache_dir`` if set, else the
    package's ``cache_plasma_weak/plasma``)."""
    return cache_write_dir(cfg, "plasma")


def list_cache_files(cfg, subdirs=None):
    """Return the sorted list of hash-named cache file paths on disk.

    Sweeps the two families of fingerprint-named caches -- ``weak/nTOp_*.txt``
    (which includes ``nTOp_thermal_*.txt``) and
    ``plasma/electron_thermo_*.txt`` -- so ``--cache-info``/``--cache-clear``
    see every file that can accumulate. The fixed-name QED pressure tables are
    excluded by construction: see the section comment above.

    Iterates over the overlay bases (the ``cache_dir`` redirect, if set, and
    the shipped package tree) so both are visible, deduplicated by
    ``<subdir>/<basename>`` with the redirect winning. Keying the dedup on the
    subdir as well as the name matters now that two subdirs are swept: two
    unrelated caches could otherwise collide on a shared basename.

    Args:
        cfg: PRIMATConfig instance.
        subdirs: optional iterable restricting the sweep, e.g. ``("weak",)``
            to reproduce the historical weak-only behaviour. ``None``
            (default) sweeps every subdir in :data:`_CACHE_PREFIXES`.

    Returns:
        list[str], sorted absolute paths.

    Example:
        >>> len(list_cache_files(cfg))                 # doctest: +SKIP
        59
        >>> len(list_cache_files(cfg, subdirs=("plasma",)))   # doctest: +SKIP
        1
    """
    wanted = _CACHE_PREFIXES if subdirs is None else {
        s: _CACHE_PREFIXES[s] for s in subdirs}
    seen = {}
    for base in _cache_bases(cfg):
        for subdir, prefixes in wanted.items():
            d = os.path.join(base, subdir)
            if not os.path.isdir(d):
                continue
            for name in os.listdir(d):
                if name.endswith(".txt") and name.startswith(prefixes):
                    seen.setdefault((subdir, name), os.path.join(d, name))
    return sorted(seen.values())


def clear_cache(cfg, subdirs=None) -> int:
    """Delete every hash-named cache file. Returns the count removed.

    Covers both families swept by :func:`list_cache_files`: the n<->p
    weak-rate tables under ``weak/`` and the e± thermodynamic tables under
    ``plasma/``.

    The cache is purely an optimisation -- every entry is reproducible from
    ``cfg`` by recomputing -- so clearing all of it is safe: the next run
    simply pays the one-time recompute cost again per configuration touched.
    That includes the copies shipped with the package, which are a
    precomputed convenience, not irreplaceable data:

    * a non-thermal ``nTOp_<hash>.txt`` recomputes deterministically (measured
      agreement with the shipped file: 9.4e-11 relative, i.e. below the
      ~1e-6 adaptive-step jitter the default ``numerical_precision=1e-7``
      already leaves in the observables);
    * an ``electron_thermo_<hash>.txt`` recomputes deterministically too, in
      ~0.7 s of quad calls;
    * a thermal ``nTOp_thermal_<hash>.txt`` is a vegas Monte-Carlo estimate, so
      the recompute agrees only to its own MC noise -- and costs minutes. Without
      vegas installed the recompute still succeeds, falling back to
      ``scipy.integrate.dblquad`` with a warning (see
      ``weak_rates.corrections``), so no install is left unable to regenerate.

    Cost, not correctness, is therefore the only thing a user forfeits here.

    Args:
        cfg: PRIMATConfig instance.
        subdirs: optional iterable restricting the sweep (see
            :func:`list_cache_files`).

    Returns:
        int, number of files actually removed.

    Example:
        >>> n = clear_cache(cfg)
        >>> print(f"removed {n} cache file(s)")
    """
    removed = 0
    for path in list_cache_files(cfg, subdirs=subdirs):
        try:
            os.remove(path)
            removed += 1
        except OSError:
            # A read-only install cannot be pruned; skip the file rather than
            # crash the --cache-clear command.
            pass
    return removed
