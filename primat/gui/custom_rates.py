# -*- coding: utf-8 -*-
"""
primat.gui.custom_rates
==========================

Helpers backing the GUI's "Customise Reactions" panel: parsing a user-uploaded
rate table, resampling it the same way :func:`primat.network_data.load_network`
does, and packing/unpacking a customisation as an in-memory zip so it can be
exported/re-imported without ever touching disk (the GUI may be running on a
read-only deployment).

The customisation itself is the small JSON-serialisable structure threaded
through ``params_items``/``PRIMAT(custom_network=...)`` (see
``primat.network_data.UpdateNuclearRates``):

    {"removed": [name, ...], "replaced": {name: raw_table_text, ...}}

``raw_table_text`` is the verbatim text of the uploaded file (2 or 3
whitespace-separated columns: T9 [GK], rate, optional uncertainty) -- *not*
pre-resampled -- so :func:`primat.network_data.load_network`'s
``_resample_rate_table`` remains the single interpolation path, applied
exactly once, both at solve time and after a zip round trip. Exports are
verbatim-on-the-original-grid for the same reason (see
:func:`verbatim_table_text`): resampling at export time and again at load
would drift from the GUI's own run by ~1e-6.
"""
import io
import json
import math
import os
import re
import zipfile

import numpy as np
import streamlit as st

from primat.network_data import (
    reaction_stoichiometry, reaction_display_name,
    load_reaction_names, _load_decay_table,
)


def sanitize_filename(title):
    """Turn a free-text network title into a safe zip/filename stem.

    Every user-visible network title (the "Create custom network" dialog's
    free-text "Network title" field, and hence the ``networks/<stem>.txt``
    entry of an exported zip, the download's own filename, and the
    ``network=`` value written into a reproduction bundle's ``.py``/``.ini``)
    goes through here, so a title containing a path separator or a space
    cannot produce a zip whose internal path no longer matches the
    ``network=`` name that is supposed to select it.

    Lives here rather than in ``params_form`` so that ``panels`` can use it
    without an import cycle (``params_form`` imports ``panels``).

    Parameters
    ----------
    title : str
        Free-text network title, possibly empty.

    Returns
    -------
    str
        ``title`` with every character outside ``[A-Za-z0-9_.-]`` collapsed to
        a single underscore and stripped from both ends; ``"custom"`` when
        nothing usable remains.

    Example
    -------
    >>> sanitize_filename("my net/v2")
    'my_net_v2'
    """
    cleaned = re.sub(r'[^A-Za-z0-9_.-]+', '_', (title or "").strip())
    return cleaned.strip("_") or "custom"


def decode_upload_text(raw):
    """Decode an uploaded file's bytes as UTF-8 text, or raise ``ValueError``.

    ``st.file_uploader`` accepts any file the user picks (the rate-table
    uploaders deliberately do not restrict ``type=``, since a rate table may
    legitimately be named ``.txt``/``.dat``/``.csv``/anything). A binary or
    non-UTF-8 file would therefore reach a bare ``bytes.decode()`` and raise
    ``UnicodeDecodeError``, which is *not* a ``ValueError`` and so escaped the
    upload handlers' own error handling as a raw traceback in the GUI.
    Funnelling every upload through here turns that into the same clean,
    catchable ``ValueError`` every other malformed-upload case already raises.

    Parameters
    ----------
    raw : bytes or str
        Raw upload contents (``st.file_uploader``'s ``getvalue()``), or text
        that is already decoded (returned unchanged).

    Returns
    -------
    str

    Raises
    ------
    ValueError
        If ``raw`` is not valid UTF-8 text.

    Example
    -------
    >>> decode_upload_text(b"1.0  2.0\\n")
    '1.0  2.0\\n'
    """
    if not isinstance(raw, bytes):
        return raw
    try:
        return raw.decode()
    except UnicodeDecodeError as exc:
        raise ValueError(
            "this file is not UTF-8 text -- a rate table must be a plain-text "
            "file with 2 or 3 numeric columns, not a binary file "
            f"(decoding failed at byte {exc.start})."
        ) from exc


def validate_new_reaction(name, data_dir=None):
    """Validate a brand-new reaction name and return its readable equation.

    Backs the GUI's "Add a new reaction" pop-up: a user types a reaction name
    in the ``a_b__c_d`` syntax (reactants and products separated by ``__``,
    nuclides within a side by ``_``; ``g`` denotes a photon, ``Bm``/``Bp`` an
    emitted electron/positron, and ``d``/``t``/``a`` alias H2/H3/He4).  The
    name need not exist in the shipped catalog: its stoichiometry is derived
    from the name itself by :func:`primat.network_data.reaction_stoichiometry`,
    which also checks baryon-number and electric-charge conservation and that
    every nuclide token is known.

    Parameters
    ----------
    name : str
        Candidate reaction name, e.g. ``"He3_d__He4_p"``.
    data_dir : str, optional
        Catalog root supplying ``nuclides.csv``/``detailed_balance.csv``,
        forwarded to :func:`primat.network_data.reaction_stoichiometry`.
        Defaults to the shipped tree; the GUI passes its own resolved root so
        a ``data_dir`` override validates against *its* nuclide table.

    Returns
    -------
    str
        A human-readable equation such as ``"He3 + d -> He4 + p"`` (reactants
        and products joined with ``+`` and separated by ``->``), suitable for
        confirming back to the user what was parsed.

    Raises
    ------
    ValueError
        If the name is empty, cannot be tokenised, has no ``__``/``TO``
        separator, references an unknown nuclide, or does not conserve A/Z.

    Example
    -------
    >>> validate_new_reaction("He3_d__He4_p")
    'He3 + d -> He4 + p'
    """
    name = (name or "").strip()
    if not name:
        raise ValueError("enter a reaction name (e.g. 'He3_d__He4_p').")
    if "__" not in name:
        raise ValueError(
            "name must use the 'a_b__c_d' syntax: reactants and products "
            "separated by a double underscore '__'.")
    try:
        react, prod = reaction_stoichiometry(name, data_dir)
    except (ValueError, KeyError) as exc:
        raise ValueError(str(exc)) from exc

    def _side(counts):
        return " + ".join(s for s, c in counts.items() for _ in range(int(c)))

    return f"{_side(react)} -> {_side(prod)}"


# Shown to the user whenever an uploaded rate table fails to parse (see
# parse_rate_upload's ValueError, caught at every upload site in
# params_form.py) -- the opening lines of a real shipped table
# (B12_t__C15_g), illustrating the expected layout: optional leading
# '#'-comments, then whitespace-separated columns T9 [GK], rate[, error].
RATE_TABLE_FORMAT_EXAMPLE = """\
# B12 + t > C15 + g   [B12_t__C15_g]   ref=TALYS2, Koning et al. 2023
# detailed balance: alpha=1.1104e+11 beta=1.5 gamma=-214.052  Q=18.4456
# T9                 rate                error
1.000000e-03   1.000000e-35   1.000000e+02
1.018629e-03   1.014102e-35   1.000000e+02
1.037605e-03   1.025274e-35   1.000000e+02
"""


def show_rate_format_help():
    """Explain the expected rate-table layout, with a real shipped example.

    Called wherever an upload fails :func:`parse_rate_upload`'s validation,
    so the user immediately sees what is expected instead of just the bare
    parse error.
    """
    st.info(
        "Expected format: optional leading lines starting with `#` "
        "(comments/provenance, ignored), followed by 2 or 3 "
        "whitespace-separated numeric columns -- `T9 [GK]`, `rate`, and an "
        "optional `error` (uncertainty factor). Example (first lines of a "
        "shipped table):"
    )
    st.code(RATE_TABLE_FORMAT_EXAMPLE, language=None)


def _validate_rate_domain(T9, rate, err):
    """Reject a parsed rate table the log-log resampler cannot use.

    :func:`primat.network_data._resample_rate_table` interpolates in
    log10-log10, so it needs finite entries and a strictly increasing,
    strictly positive ``T9`` column; scipy's cubic additionally refuses
    duplicate abscissae. Without this check the bad table reaches the solver,
    which either fails with a message that never mentions the upload or --
    for an unsorted ``T9`` or a negative rate -- extrapolates nonsense and
    reports it as an ordinary result (pinned in
    ``tests/test_gui_robustness.py``).

    Parameters
    ----------
    T9, rate, err : np.ndarray
        The three columns as parsed by :func:`parse_rate_upload`.

    Raises
    ------
    ValueError
        Naming the offending row (1-based, as a text editor counts) and value.
    """
    for label, col in (("T9", T9), ("rate", rate), ("error", err)):
        bad = np.flatnonzero(~np.isfinite(col))
        if bad.size:
            raise ValueError(
                f"{label} column has a non-finite value ({col[bad[0]]}) at "
                f"data row {bad[0] + 1}; every entry must be a finite number."
            )
    bad = np.flatnonzero(T9 <= 0.0)
    if bad.size:
        raise ValueError(
            f"T9 must be strictly positive (rates are interpolated in "
            f"log-log), but data row {bad[0] + 1} has T9 = {T9[bad[0]]:g}."
        )
    step = np.diff(T9)
    bad = np.flatnonzero(step <= 0.0)
    if bad.size:
        i = int(bad[0])
        how = "repeats" if step[i] == 0.0 else "goes backwards"
        raise ValueError(
            f"the T9 column must increase strictly down the file, but it "
            f"{how} at data row {i + 2} ({T9[i]:g} -> {T9[i + 1]:g}). "
            "Sort the table by ascending T9."
        )
    bad = np.flatnonzero(rate < 0.0)
    if bad.size:
        raise ValueError(
            f"the rate column must not be negative, but data row "
            f"{bad[0] + 1} has rate = {rate[bad[0]]:g}."
        )


def parse_rate_upload(fh, cfg=None, warn=True):
    """Parse an uploaded rate-table file into raw ``(T9, rate, err, header)``.

    Parameters
    ----------
    fh : file-like, bytes or str
        2- or 3-column whitespace-separated text (as produced by
        ``st.file_uploader``, or a plain file object): ``T9 [GK]``, ``rate``,
        and an optional third uncertainty column.  Leading ``#``-prefixed
        lines are the uploader's own header/provenance comment, preserved
        verbatim (see ``header``) rather than discarded.
    cfg : PRIMATConfig, optional
        Supplies the master grid's span (``rate_grid_T9_min``/
        ``rate_grid_T9_max``) for the coverage warning below.  Omit it (the
        default) to fall back on the shipped defaults of 0.001-10 GK -- the
        warning is then only advisory, so a caller with no ``cfg`` to hand
        still gets a sensible message.
    warn : bool
        Whether to emit the coverage ``st.warning``.  Callers that are merely
        *re-parsing* an already-accepted table (e.g. :func:`export_zip`
        recovering a stored upload's arrays, which runs on every Streamlit
        rerun) pass ``False``: the warning belongs to the moment the user
        actually uploads the file, not to every later render that happens to
        touch it again.

    Returns
    -------
    (T9, rate, err, header) : tuple of (np.ndarray, np.ndarray, np.ndarray, list[str])
        ``err`` is an all-zero array of the same length when the upload has
        only 2 columns.  ``header`` is the list of the upload's own leading
        ``#``-prefixed lines (possibly empty), preserved so a re-exported zip
        carries the original provenance rather than a generic "custom upload"
        label (see :func:`verbatim_table_text`).

    Raises
    ------
    ValueError
        If the file is not UTF-8 text (see :func:`decode_upload_text`), does
        not parse as 2 or 3 numeric columns, or fails
        :func:`_validate_rate_domain` (non-finite entry, non-positive or
        non-increasing ``T9``, negative rate).

    Notes
    -----
    Emits an ``st.warning`` (not an error -- the table is still usable) if its
    T9 range does not cover the master grid's span: outside the upload's own
    range, ``_resample_rate_table`` extrapolates by continuing the table's end
    slope in log-log, which drifts from the truth the further out it goes.
    (That function raises its own ``UserWarning`` too; this one surfaces the
    problem in the GUI at upload time, before a run is launched.)

    Example
    -------
    >>> T9, rate, err, header = parse_rate_upload("1.0  2.0\\n2.0  3.0\\n")
    """
    if hasattr(fh, "read"):
        text = fh.read()
    else:
        text = fh
    text = decode_upload_text(text)
    header = [line for line in text.splitlines() if line.startswith("#")]
    data = np.loadtxt(io.StringIO(text), unpack=True)
    if data.ndim != 2 or data.shape[0] not in (2, 3):
        raise ValueError(
            f"expected 2 or 3 columns (T9, rate[, err]), got shape {data.shape}"
        )
    T9, rate = data[0], data[1]
    err = data[2] if data.shape[0] == 3 else np.zeros_like(rate)
    _validate_rate_domain(T9, rate, err)
    # Read the grid span off cfg rather than hard-coding it, so the warning
    # stays truthful when rate_grid_T9_min/max are overridden.
    T9_min = getattr(cfg, "rate_grid_T9_min", 1.0e-3)
    T9_max = getattr(cfg, "rate_grid_T9_max", 10.0)
    if warn and (T9.min() > T9_min or T9.max() < T9_max):
        st.warning(
            f"Uploaded table spans T9 = [{T9.min():.3g}, {T9.max():.3g}] GK, "
            f"narrower than the standard grid [{T9_min:.3g}, {T9_max:.3g}] GK "
            "-- values outside this range are extrapolated."
        )
    return T9, rate, err, header


def stamp_upload(name, raw_text):
    """Prepend a provenance header to a freshly uploaded rate-table text.

    Called right where a user's uploaded file is first accepted as a custom
    rate table (the "New rate table for <name>" and "Add a new rate"
    uploaders in ``params_form``), so that *every* later view of this text --
    the "Show rate table" preview popup (:func:`primat.gui.params_form
    ._current_table_text`), the per-reaction "Source" column
    (:func:`primat.network_data._reaction_source_from_lines`, which reads
    this very header back), and a re-exported zip -- carries an unambiguous
    "this is a primat-loaded custom table for reaction X" label, even when
    the uploaded file itself had no ``#`` header at all.

    Parameters
    ----------
    name : str
        Bare reaction name this table is for.
    raw_text : str
        The verbatim uploaded file contents (already validated by
        :func:`parse_rate_upload`).

    Returns
    -------
    str
        ``raw_text`` with two new leading lines: a one-line provenance
        comment naming the reaction (in the same human-readable
        ``"react1 + react2 > prod1 + prod2   [name]"`` form as the shipped
        tables' own headers, see :func:`reaction_display_name`), then a full
        line of ``#`` as a visual separator from whatever header the upload
        itself carried.
    """
    lines = [
        f"# {reaction_display_name(name)}   [{name}]   (custom rate)",
        "#" * 70,
    ]
    return "\n".join(lines) + "\n" + raw_text


def _strip_own_stamp(name, header):
    """Drop :func:`stamp_upload`'s own two-line preamble from ``header``.

    ``header`` (as returned by :func:`parse_rate_upload`) is read straight
    off an already-:func:`stamp_upload`-ed table -- e.g. when
    :func:`export_zip` re-parses a stored "kept" table to recover its arrays
    -- so it starts with ``stamp_upload``'s own bookkeeping lines
    (``"# {react} > {prod}   [{name}]   (custom rate)"`` + a ``"#"*70``
    fence). Passing those straight through to :func:`verbatim_table_text` as
    ``source_header`` would duplicate that bookkeeping underneath the new
    header it writes itself; only a genuine header carried by the *original*
    upload (if any), following that preamble, is worth preserving.

    Parameters
    ----------
    name : str
        Bare reaction name (must match the preamble's own ``name``).
    header : sequence[str]
        Leading ``#``-lines as returned by :func:`parse_rate_upload`.

    Returns
    -------
    list[str]
    """
    header = list(header)
    own_stamp = f"# {reaction_display_name(name)}   [{name}]   (custom rate)"
    if len(header) >= 2 and header[0] == own_stamp and header[1] == "#" * 70:
        return header[2:]
    return header


def verbatim_table_text(T9, rate, err, name="custom", source_header=()):
    """Return table text with the upload on its ORIGINAL grid, full precision.

    Writes the parsed ``(T9, rate, err)`` arrays verbatim -- same grid points
    the user uploaded, at full float64 precision (``%.17e``). It is what an
    exported network zip and a reproduction bundle carry for a
    replaced/added reaction, so that ``load_network``'s
    ``_resample_rate_table`` runs *once* on exactly the data the GUI's own
    live run resampled -- reproducing the run bit-for-bit.

    Pre-resampling onto the master grid at export time instead would break
    that bit-for-bit match two ways: the exported values are rounded
    (``%.6e``), and resampling a coarse upload onto the wider master grid
    *extrapolates*, so the overlay would then resample that
    extrapolated+rounded table a second time -- diverging from the GUI's
    single raw resample by ~1e-6. This is why the export is deliberately
    verbatim-on-the-original-grid rather than "the effective on-grid table".

    Parameters
    ----------
    T9, rate, err : np.ndarray
        Raw uploaded arrays, as returned by :func:`parse_rate_upload`.
    name : str
        Reaction name, for the header's ``ref=`` field.
    source_header : sequence[str]
        The uploader's own ``#``-prefixed header lines, preserved verbatim.

    Returns
    -------
    str
        Table text: ``#`` header line(s) then ``T9 rate err`` rows on the
        upload's own grid, full precision.
    """
    lines = [
        f"# {reaction_display_name(name)}   [{name}]   "
        "(custom rate, original grid, verbatim)",
        "#" * 70,
    ]
    lines.extend(source_header)
    for t9, r, e in zip(np.asarray(T9), np.asarray(rate), np.asarray(err)):
        lines.append(f"{t9:.17e}   {r:.17e}   {e:.17e}")
    return "\n".join(lines)


def decay_override_table_text(name, rate_s):
    """Synthetic constant-rate table text for a user-overridden decay rate.

    Decay rates are T9-independent (see ``_load_decay_table``); routing an
    override through ``load_network``'s existing ``custom_tables`` mechanism
    (checked *before* the decays.txt branch, so an override always wins, see
    ``load_network``'s rate-loading loop) needs at least 4 points for the
    cubic log-log resampling in :func:`primat.network_data._resample_rate_table`,
    so this repeats the same rate across a handful of grid points rather than
    using decays.txt's single-row format. Used both to feed the solver and
    (rarely, if a decay override could not be expressed via the dedicated
    ``tables/decays.txt`` zip entry -- see :func:`export_zip`) as a fallback
    table representation.
    """
    grid = (1.0e-3, 1.0e-2, 1.0e-1, 1.0, 10.0)
    lines = [f"# {name}: decay rate overridden in primat (was log(2)/halflife)"]
    lines += [f"{t9:.6e}   {rate_s:.6e}   {1.0:.6e}" for t9 in grid]
    return "\n".join(lines)


def _base_network_filenames(cfg):
    """``{bare_reaction_name: table_filename}`` for ``cfg.network``'s own list.

    A network file may pin a *specific* rate table per reaction with the
    ``name, filename.txt`` syntax -- ``small_parthenope.txt`` does so for all
    12 of its entries (``n_p__d_g, n_p__d_g_parthenope3.0.txt``, ...), and
    ``large.txt`` spells out ``*_primat.txt`` explicitly. :func:`export_zip`
    needs that mapping to copy the table the run *actually used* for a
    reaction the user never customised; assuming ``<name>_primat.txt`` instead
    silently substituted primat's own rates for Parthenope's, exporting a
    ``small_parthenope`` network that re-imported to different abundances
    (D/H 2.4999622e-05 -> 2.4358771e-05, Li7/H 4.812238e-10 -> 5.557655e-10).

    Reactions listed by bare name (no comma) are absent from the returned map:
    they have no pinned filename, so the caller's ``<name>_primat.txt``
    default is the right answer for them.

    Parameters
    ----------
    cfg : PRIMATConfig
        Its ``network`` names the list to read. A network with no file on disk
        ('small', or a name driven entirely by ``custom_network``) simply
        yields an empty map.

    Returns
    -------
    dict[str, str]

    Example
    -------
    >>> _base_network_filenames(PRIMATConfig({"network": "small_parthenope"}))
    {'n_p__d_g': 'n_p__d_g_parthenope3.0.txt', ...}
    """
    try:
        entries = load_reaction_names(cfg, cfg.network)
    except (ValueError, KeyError, OSError):
        # An unreadable/absent network list is not worth failing an export
        # over -- fall back to the per-reaction default for every reaction.
        return {}
    pinned = {}
    for entry in entries:
        parts = re.split(r'[, ]+', entry, maxsplit=1)
        if len(parts) > 1 and parts[1].strip():
            pinned[parts[0].strip()] = parts[1].strip()
    return pinned


def _safe_basename(candidate, fallback):
    """``candidate`` as a bare filename safe to write into a zip, else ``fallback``.

    :func:`export_zip` writes this straight into a ``tables/<name>/<file>``
    archive path. ``candidate`` may have come from an imported zip, so a value
    like ``../..`` would make primat author an archive that escapes its own
    directory on extraction. :func:`sanitize_filename` already collapses path
    separators, but leaves a pure ``..`` intact -- hence the explicit check.

    Parameters
    ----------
    candidate : str or None
        Preferred basename (the filename agreed at upload time, or an imported
        zip's own).
    fallback : str
        Used when ``candidate`` is empty or does not survive sanitisation.

    Returns
    -------
    str
    """
    if not candidate:
        return fallback
    cleaned = sanitize_filename(candidate)
    if cleaned in (os.curdir, os.pardir):
        return fallback
    return cleaned


def _shipped_table_dir(cfg, name):
    return os.path.join(cfg._resolved_data_dir, "nuclear", "tables", name)


def _match_shipped_file(cfg, name, raw_text):
    """If ``raw_text`` is byte-identical to an on-disk ``tables/<name>/*.txt``
    file, return that file's basename; else ``None``.

    Distinguishes "the user picked an existing alternate shipped table from
    the dropdown" (e.g. a ``*_parthenope3.0.txt`` sibling) -- which keeps its
    real name and content unaltered -- from "the user actually uploaded new
    content", which gets the ``_newnetwork`` treatment in :func:`export_zip`.
    """
    folder = _shipped_table_dir(cfg, name)
    try:
        candidates = os.listdir(folder)
    except OSError:
        return None
    for fname in candidates:
        if not fname.endswith(".txt"):
            continue
        try:
            with open(os.path.join(folder, fname)) as f:
                if f.read() == raw_text:
                    return fname
        except OSError:
            continue
    return None




# Fixed timestamp for every zip entry (1980-01-01, the earliest a DOS-format
# zip can express). Without it zipfile stamps each entry with the wall clock,
# so two exports of the same network are byte-different and a user cannot
# checksum a download against a re-export.
_ZIP_TIMESTAMP = (1980, 1, 1, 0, 0, 0)


def export_zip(cfg, custom_network, kept_names, network_filename="custom"):
    """Pack a customisation into an in-memory zip mirroring the repo layout.

    Parameters
    ----------
    cfg : PRIMATConfig
        Supplies the master-grid parameters (for resampling replaced tables).
    custom_network : dict
        ``{"removed": [...], "replaced": {name: raw_text, ...},
        "added": {name: raw_text, ...}}``.
    kept_names : sequence[str]
        The full ordered list of reaction names actually in the network after
        removal *and* additions (i.e. ``cfg.network``'s list minus
        ``custom_network["removed"]``, plus ``custom_network["added"]``).  In
        the GUI this is read off the solved network's reaction list, so added
        reactions are already included.
    network_filename : str
        Basename (without ``.txt``) for the network file under ``networks/``,
        e.g. the user-chosen custom-network title (sanitised). Defaults to
        ``"custom"`` for callers that don't have a user-chosen title (e.g. the
        post-run Reactions-tab export of a legacy "Customise Reactions"
        session).

    Returns
    -------
    bytes
        Zip file contents with ``networks/<network_filename>.txt`` (one
        reaction name per line; every per-reaction-table reaction is written
        as ``name, <filename>`` -- the filename is *always* explicit, never
        implied, even for an unmodified shipped default) and one
        ``tables/<name>/<filename>`` per such reaction:

        * An unmodified reaction (still using a shipped table, default or an
          alternate like ``*_parthenope3.0.txt``) keeps its real, unaltered
          filename and content -- so picking an existing alternate from the
          dropdown is never confused with a genuine edit. Which table that is
          comes from ``cfg.network``'s own ``name, filename`` pairing (see
          :func:`_base_network_filenames`), *not* from assuming
          ``<name>_primat.txt``.
        * A genuinely new/uploaded/edited table is written by
          :func:`verbatim_table_text` (the upload on its own original grid, at
          full precision, under the basename agreed at upload time), so that
          re-importing it resamples exactly once and reproduces the run
          bit-for-bit.
        * A decay reaction (Bm/Bp, rate from the shared ``decays.txt``, not a
          per-reaction file) has no table file at all; if its rate has been
          overridden the network-file line is instead ``name, <rate_s>`` --
          a bare number where every other reaction has a filename. There is
          no separate ``decays.txt`` entry in the zip: ``rate_s`` *is* the
          override, right there in the line that names the reaction.

        The zip is fully self-contained: *every* kept reaction's table is
        included, not just replaced/added ones, so it reproduces the exact
        network even on an install whose shipped tables might differ.
    """
    # Replaced (override a kept reaction) and added (brand-new) reactions are
    # both backed by an uploaded table, so they are written to the zip the same
    # way -- merge them into one map of custom tables.
    custom_tables = {**custom_network.get("replaced", {}),
                     **custom_network.get("added", {})}
    decay_table = _load_decay_table(os.path.join(cfg._resolved_data_dir, "nuclear", "tables"))
    # Which on-disk table each *uncustomised* reaction actually uses, per
    # cfg.network's own list -- e.g. "*_parthenope3.0.txt" for
    # small_parthenope. Assuming "<name>_primat.txt" here used to silently
    # export the wrong rates for any such network.
    pinned_filenames = _base_network_filenames(cfg)

    # Decay-reaction overrides are pulled out of custom_tables here: they get
    # their rate written inline in the network file, not a per-reaction
    # table file.
    decay_overrides = {}
    for name in list(custom_tables):
        if name in decay_table:
            raw_text = custom_tables.pop(name)
            try:
                T9, rate, err, header = parse_rate_upload(
                    raw_text, cfg=cfg, warn=False)
                decay_overrides[name] = float(np.asarray(rate).reshape(-1)[0])
            except (ValueError, IndexError):
                pass

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        # Entries are written through _write_entry, which stamps a fixed
        # timestamp: zipfile defaults each entry to the wall clock, so two
        # exports of the same network minutes apart differ in bytes and cannot
        # be checksummed against each other -- which is the whole point of the
        # verbatim, original-grid export (see this module's docstring on
        # bit-for-bit round trips).
        def zf_writestr(name, data):
            info = zipfile.ZipInfo(name, date_time=_ZIP_TIMESTAMP)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, data)

        lines = []
        for name in kept_names:
            if name in decay_table:
                if name in decay_overrides:
                    lines.append(f"{name}, {decay_overrides[name]:.6e}")
                else:
                    lines.append(name)
                continue
            if name in custom_tables:
                raw_text = custom_tables[name]
                shipped_name = _match_shipped_file(cfg, name, raw_text)
                if shipped_name is not None:
                    # Just an existing shipped table picked from the
                    # dropdown -- not an edit. The shipped default already
                    # carries the "_primat" suffix on disk (see
                    # convert_ac2024_rates.py), so it reads unambiguously as
                    # "primat's own rate"; an already-distinctly-named
                    # alternate (e.g. "*_parthenope3.0.txt") keeps its name.
                    lines.append(f"{name}, {shipped_name}")
                    zf_writestr(f"tables/{name}/{shipped_name}", raw_text)
                    continue
                try:
                    T9, rate, err, header = parse_rate_upload(
                        raw_text, cfg=cfg, warn=False)
                    # Write the upload on its ORIGINAL grid at full precision --
                    # NOT pre-resampled onto the master grid -- so load_network
                    # resamples it exactly once, identically to the GUI's own
                    # live run, reproducing the run bit-for-bit. Pre-resampling
                    # here (effective_table_text) rounds + extrapolates and then
                    # gets resampled again at load, drifting ~1e-6.
                    table_text = verbatim_table_text(
                        T9, rate, err, name=name,
                        source_header=_strip_own_stamp(name, header))
                except ValueError:
                    table_text = raw_text
                # Prefer the basename already agreed on at upload time
                # ("<name>_custom_<uploaded filename>", see the "New rate
                # table for <name>" uploader in params_form.py) so a
                # downloaded zip's filename matches what the dialog itself
                # showed throughout editing. Only legacy callers with no
                # "filenames" entry (e.g. the post-run Reactions-tab export
                # of an old-style "Customise Reactions" session) fall back
                # to a generic name suffixed with this network's own title.
                fname = _safe_basename(
                    custom_network.get("filenames", {}).get(name),
                    f"{name}_{network_filename}.txt")
                lines.append(f"{name}, {fname}")
                zf_writestr(f"tables/{name}/{fname}", table_text)
                continue
            # Unmodified shipped reaction: copy its on-disk table verbatim (no
            # resampling needed, it is already a valid rate file) under its own
            # real name, so the zip does not depend on the importing install's
            # own shipped tables/<name>/ folder. The filename comes from the
            # base network's own pairing when it pins one (small_parthenope's
            # "*_parthenope3.0.txt"), falling back to primat's default table
            # for a network that lists the reaction by bare name.
            fname = pinned_filenames.get(name, f"{name}_primat.txt")
            path = os.path.join(cfg._resolved_data_dir, "nuclear", "tables", name, fname)
            try:
                with open(path) as f:
                    table_text = f.read()
            except OSError:
                lines.append(name)
                continue
            lines.append(f"{name}, {fname}")
            zf_writestr(f"tables/{name}/{fname}", table_text)
        zf_writestr(f"networks/{network_filename}.txt", "\n".join(lines) + "\n")
    return buf.getvalue()


@st.cache_data(show_spinner=False, max_entries=8)
def _export_zip_cached(custom_network_json, kept_names, network_filename,
                       network, data_dir, _cfg):
    """Memoised :func:`export_zip`. See :func:`export_zip_cached`.

    ``_cfg`` is underscore-prefixed so Streamlit excludes it from the cache
    key (a ``PRIMATConfig`` is not hashable); the two ``cfg`` attributes that
    actually change the output -- ``network`` (which table each unmodified
    reaction uses, see :func:`_base_network_filenames`) and the resolved data
    root -- are passed separately as ordinary hashable key components.
    """
    return export_zip(_cfg, json.loads(custom_network_json), list(kept_names),
                      network_filename=network_filename)


def export_zip_cached(cfg, custom_network, kept_names, network_filename="custom"):
    """:func:`export_zip`, memoised on its inputs.

    Every "Download network (zip)" button builds its bytes *eagerly* --
    Streamlit's ``st.download_button`` takes the data up front, it cannot pull
    them lazily at click time -- and the buttons live inside panels/dialogs
    that re-render on every widget interaction. For the large network that is
    428 rate-table files read off disk and deflated into a 4.3 MB archive,
    measured at ~0.47 s, paid again on every single reaction toggle inside the
    "Create custom network" dialog. Memoising collapses that to once per
    distinct (customisation, kept set, title, network, data root).

    Parameters are exactly :func:`export_zip`'s; ``custom_network`` must be
    JSON-serialisable (it always is -- the GUI already round-trips it through
    ``json.dumps`` to reach ``PRIMAT(custom_network=...)``).

    Returns
    -------
    bytes
        The zip contents, identical to what :func:`export_zip` returns.

    Example
    -------
    >>> export_zip_cached(cfg, {"removed": [], "replaced": {}, "added": {}},
    ...                   ["n_p__d_g"], network_filename="small")
    """
    return _export_zip_cached(
        json.dumps(custom_network, sort_keys=True), tuple(kept_names),
        network_filename, cfg.network, cfg._resolved_data_dir, _cfg=cfg)


# Budget for an uploaded zip, enforced by _check_zip_budget before anything is
# decompressed. A zip's compression ratio is unbounded, so without these a
# small upload can allocate arbitrarily much -- fatal on the ~1 GB public demo
# at primat.streamlit.app. Each is roughly 4x the corresponding figure for the
# largest archive primat itself produces (the full large network), so no
# legitimate export comes close; test_gui_robustness.py pins both directions.
_ZIP_MAX_TOTAL_BYTES = 64_000_000
_ZIP_MAX_ENTRY_BYTES = 4_000_000
_ZIP_MAX_ENTRIES = 2000


def _check_zip_budget(zf):
    """Reject an uploaded zip that would decompress to more than the budget.

    Reads only the central directory (``infolist``), so an over-budget archive
    is refused *before* a single byte is decompressed. The declared sizes can
    lie, which is why :func:`_read_zip_text` caps each read as well.

    Raises
    ------
    ValueError
        Quoting the offending size against its limit, in MB.
    """
    infos = zf.infolist()
    if len(infos) > _ZIP_MAX_ENTRIES:
        raise ValueError(
            f"the archive has {len(infos)} entries, more than the "
            f"{_ZIP_MAX_ENTRIES} a network zip may contain."
        )
    total = sum(info.file_size for info in infos)
    if total > _ZIP_MAX_TOTAL_BYTES:
        raise ValueError(
            f"the archive expands to {total / 1e6:.1f} MB, more than the "
            f"{_ZIP_MAX_TOTAL_BYTES / 1e6:.0f} MB a network zip may contain."
        )
    for info in infos:
        if info.file_size > _ZIP_MAX_ENTRY_BYTES:
            raise ValueError(
                f"'{info.filename}' expands to {info.file_size / 1e6:.1f} MB, "
                f"more than the {_ZIP_MAX_ENTRY_BYTES / 1e6:.0f} MB a single "
                "rate table may occupy."
            )


def _read_zip_text(zf, name):
    """Read one zip member as UTF-8 text, capped at ``_ZIP_MAX_ENTRY_BYTES``.

    The cap backs up :func:`_check_zip_budget`, which can only trust the
    central directory's *declared* sizes. Decoding and corruption errors are
    re-raised naming the member, so the GUI shows "tables/x/y.txt is not UTF-8
    text" rather than a bare codec message about a byte offset in an unnamed
    file, or a ``BadZipFile`` that is not a ``ValueError`` at all.

    Raises
    ------
    ValueError
        If the member overruns the cap, is corrupt, or is not UTF-8.
    """
    try:
        with zf.open(name) as f:
            raw = f.read(_ZIP_MAX_ENTRY_BYTES + 1)
    except zipfile.BadZipFile as exc:
        raise ValueError(
            f"'{name}' is corrupt and could not be read from the archive "
            f"({exc}). Re-download or re-export the zip."
        ) from exc
    if len(raw) > _ZIP_MAX_ENTRY_BYTES:
        raise ValueError(
            f"'{name}' is larger than the {_ZIP_MAX_ENTRY_BYTES / 1e6:.0f} MB "
            "a single entry may occupy (its declared size understated it)."
        )
    try:
        return raw.decode()
    except UnicodeDecodeError as exc:
        raise ValueError(
            f"'{name}' is not UTF-8 text -- a network zip must contain only "
            f"plain-text rate tables (decoding failed at byte {exc.start})."
        ) from exc


def _check_member_component(component, filename):
    """Reject a zip member whose path component escapes its own directory.

    ``import_zip`` splits ``tables/<name>/<file>`` on ``/``, so a member named
    ``tables/../evil.txt`` yields the component ``..``. primat never extracts
    an upload to disk, but :func:`export_zip` copies that component back into
    the path it writes, so an unchecked one would leave primat handing the
    user an archive that escapes its own directory on any naive extractor.

    Raises
    ------
    ValueError
    """
    if component in ("", os.curdir, os.pardir) or "\\" in component:
        raise ValueError(
            f"'{filename}' has an unsafe path component {component!r}; a "
            "network zip must use plain 'tables/<reaction>/<file>' entries."
        )


def import_zip(fh):
    """Rebuild a ``custom_network`` dict from a zip produced by :func:`export_zip`.

    Parameters
    ----------
    fh : file-like
        The uploaded zip file.

    Returns
    -------
    dict
        ``{"kept": [name, ...], "replaced": {name: raw_text, ...},
        "filenames": {name: basename, ...}, "decay_overrides":
        {name: rate_s, ...}, "title": str}``.
        Removal is implicit: the single file under ``networks/`` only lists
        the reactions that were *kept*, so any reaction of the selected
        network absent from ``kept`` is treated as removed by the caller.
        Brand-new (added) reactions also appear in ``kept`` with their table
        in ``replaced``; the caller tells them apart from replacements by
        checking which ``kept`` names do *not* belong to the selected
        network.  ``replaced`` carries one entry per kept reaction that has a
        per-reaction table file in the zip (shipped-default, alternate, or
        genuinely new -- see :func:`export_zip`); decay reactions have no
        such file and instead contribute to ``decay_overrides`` if their
        network-file line carries a bare number (the overridden rate)
        instead of a filename. ``filenames`` carries the zip's own basename
        for each such reaction (e.g. ``"B8_d__Be7_He3_primat.txt"``), so the
        Reactions tab's "File" column shows *something* meaningful after a
        round trip rather than ``None`` -- this is purely the in-zip
        filename, not a real on-disk path. ``title`` is the network file's
        basename (without ``.txt``), recovered without needing a separate
        metadata file.

    Raises
    ------
    ValueError
        If the upload is not a zip, is over the size budget (see
        ``_ZIP_MAX_TOTAL_BYTES``), has no single ``networks/`` file, contains
        an unsafe member path or a non-UTF-8/corrupt member, or names a rate
        table it does not carry. Every caller renders this as an ``st.error``.
    """
    replaced = {}
    filenames = {}
    decay_overrides = {}
    try:
        zf = zipfile.ZipFile(fh)
    except zipfile.BadZipFile:
        raise ValueError(
            "the uploaded file is not a valid zip archive (expected one "
            "produced by the 'Download network details' button)."
        ) from None
    with zf:
        _check_zip_budget(zf)
        net_files = [info.filename for info in zf.infolist()
                    if info.filename.startswith("networks/")
                    and info.filename.endswith(".txt")]
        if len(net_files) != 1:
            raise ValueError(
                f"expected exactly one file under 'networks/', found {len(net_files)}."
            )
        net_filename = net_files[0]
        title = os.path.basename(net_filename)[: -len(".txt")]
        net_text = _read_zip_text(zf, net_filename)
        kept_names = []
        declared = {}
        for line in net_text.splitlines():
            line = line.strip()
            if not line:
                continue
            parts = re.split(r'[, ]+', line, maxsplit=1)
            bare = parts[0].strip()
            kept_names.append(bare)
            if len(parts) > 1:
                # A decay reaction's overridden rate is written directly as
                # a bare number in the spot every other reaction uses for a
                # filename (see export_zip) -- no separate decays.txt entry.
                try:
                    decay_overrides[bare] = float(parts[1].strip())
                except ValueError:
                    declared[bare] = parts[1].strip()  # a filename, not a rate
        for info in zf.infolist():
            if info.filename.startswith("tables/") and info.filename.count("/") == 2:
                # "tables/<name>/<filename>" -- any per-reaction table file,
                # default-named, alternate-shipped, or genuinely new.
                bare, fname = info.filename.split("/")[1:3]
                _check_member_component(bare, info.filename)
                _check_member_component(fname, info.filename)
                replaced[bare] = _read_zip_text(zf, info.filename)
                filenames[bare] = fname
        # A line naming a table file the archive does not actually carry is a
        # corrupt or hand-edited zip. Left unreported, the import succeeds and
        # the run falls back to *this install's* shipped table -- exactly the
        # substitution the verbatim, original-grid export exists to prevent.
        missing = sorted(set(declared) - set(replaced))
        if missing:
            shown = ", ".join(f"{n} ({declared[n]})" for n in missing[:3])
            more = f" (+{len(missing) - 3} more)" if len(missing) > 3 else ""
            raise ValueError(
                f"the network file names {len(missing)} rate table(s) that are "
                f"not in the archive: {shown}{more}. The zip is incomplete -- "
                "re-export it from the 'Download network details' button."
            )
    return {"kept": kept_names, "replaced": replaced, "filenames": filenames,
            "decay_overrides": decay_overrides, "title": title}


def kept_to_custom_network(cfg, kept, replaced, decay_overrides=None, filenames=None):
    """Build the ``{"removed", "replaced", "added"}`` dict from an imported zip.

    Shared by both the sidebar's "Import custom network" dialog
    (``primat.gui.params_form``) and the post-run Reactions tab's own
    importer (``primat.gui.panels``) -- lives here, not in ``params_form``,
    so ``panels`` can call it without a circular import (``params_form``
    already imports from ``panels``).

    ``removed`` is computed against the *full, unfiltered* large-network
    reaction list -- not some amax-restricted view -- so that every
    catalog reaction absent from ``kept`` is actually excluded from the
    solved network. (An earlier version derived an "implied amax" from the
    heaviest category among ``kept`` and only marked reactions *within* that
    band as removed; every reaction above it was then neither removed nor
    kept, so ``UpdateNuclearRates`` silently treated "not removed" as "keep"
    and the imported network solved with hundreds of unwanted extra
    reactions -- this is why that derivation is gone.)

    Parameters
    ----------
    cfg : PRIMATConfig
        Used only to resolve ``data/nuclear/networks/large.txt``.
    kept : sequence[str]
        Reaction names kept in the imported network.
    replaced : dict[str, str]
        ``{name: raw_table_text}`` for every reaction the zip carried a table
        for (the exported zip's format includes one for *every* kept
        reaction, not just genuinely customised ones -- see ``export_zip``).
    decay_overrides : dict[str, float], optional
        ``{name: rate_s}`` parsed from the zip's ``tables/decays.txt`` (see
        :func:`import_zip`). Only entries that actually differ from the
        shipped ``decays.txt`` rate are turned into a synthetic
        ``replaced`` table entry (:func:`decay_override_table_text`) --
        an unmodified decay reaction needs no override at all.
    filenames : dict[str, str], optional
        ``{name: basename}`` from :func:`import_zip`, the in-zip filename
        for each reaction in ``replaced`` -- threaded through into the
        returned dict's own ``"filenames"`` key purely so the Reactions
        tab's "File" column has *something* to show after a round trip
        (``UpdateNuclearRates`` reads it via ``custom_network["filenames"]``;
        see its docstring) instead of ``None``.

    Returns
    -------
    dict
        ``{"removed": [...], "replaced": {...}, "added": {...},
        "filenames": {...}}``, the shape ``UpdateNuclearRates`` expects.
    """
    entries = load_reaction_names(cfg, "large")
    bare_names = {re.split(r'[, ]+', e, maxsplit=1)[0].strip() for e in entries}
    kept_set = set(kept)
    removed = sorted(bare_names - kept_set)
    added = {n: replaced[n] for n in kept_set - bare_names if n in replaced}
    true_replaced = {n: t for n, t in replaced.items() if n not in added}
    if decay_overrides:
        shipped = _load_decay_table(os.path.join(cfg._resolved_data_dir, "nuclear", "tables"))
        for name, rate_s in decay_overrides.items():
            shipped_entry = shipped.get(name)
            if shipped_entry is None or not math.isclose(
                rate_s, shipped_entry[0], rel_tol=1e-9, abs_tol=0.0
            ):
                true_replaced[name] = decay_override_table_text(name, rate_s)
    true_filenames = {n: f for n, f in (filenames or {}).items() if n in true_replaced}
    return {"removed": removed, "replaced": true_replaced, "added": added,
            "filenames": true_filenames}
