# -*- coding: utf-8 -*-
"""Adversarial GUI input handling: bad uploads must not become bad physics.

Goal: every hostile or careless upload the GUI accepts -- a malformed rate
table, a zip bomb, an archive whose manifest disagrees with its contents --
surfaces as a clean ``ValueError`` the dialogs already render as an
``st.error``, never as a traceback and never as a plausible-looking wrong
number. Pass 20's headline was that an unsorted T9 column ran to completion
and reported YP = 0.00000711 against a true 0.24699907.

These exercise ``primat.gui.custom_rates`` directly rather than through
``AppTest``: the checks are pure functions of the upload, so a unit test pins
them in a second instead of a minute.
"""
import io
import zipfile

import numpy as np
import pytest

pytest.importorskip("streamlit")

from primat.gui import custom_rates as cr


VALID_TABLE = "1e-3 1e-3\n1e-2 2e-3\n1e-1 3e-3\n1.0 4e-3\n10.0 5e-3\n"


def _zip(entries, compression=zipfile.ZIP_DEFLATED):
    """In-memory zip of ``{member_name: bytes_or_text}``."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression) as zf:
        for name, data in entries.items():
            zf.writestr(name, data)
    buf.seek(0)
    return buf


# ---------------------------------------------------------------------------
# Rate-table domain validation
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,text,expected", [
    ("unsorted T9", "1.0 4e-3\n1e-3 1e-3\n10.0 5e-3\n", "goes backwards"),
    ("duplicate T9", "1e-3 1e-3\n1e-3 2e-3\n1e-1 3e-3\n", "repeats"),
    ("zero T9", "0.0 1e-3\n1e-2 2e-3\n1e-1 3e-3\n", "strictly positive"),
    ("negative T9", "-1.0 1e-3\n1e-2 2e-3\n1e-1 3e-3\n", "strictly positive"),
    ("negative rate", "1e-3 1e-3\n1e-2 -2e-3\n1e-1 3e-3\n", "must not be negative"),
    ("NaN rate", "1e-3 1e-3\n1e-2 nan\n1e-1 3e-3\n", "non-finite"),
    ("inf T9", "1e-3 1e-3\n inf 2e-3\n1e-1 3e-3\n", "non-finite"),
])
def test_degenerate_rate_table_is_rejected_at_upload(label, text, expected):
    """A table the log-log resampler cannot use is refused, naming the row.

    Every one of these used to be accepted. Three then aborted the solve with
    ``cpr_ode_bdf: step size underflowed below machine precision``, which
    never mentions the upload; the unsorted and negative-rate cases ran to
    completion and reported YP ~ 7e-06 against a true 0.24699907.
    """
    with pytest.raises(ValueError, match=expected):
        cr.parse_rate_upload(text)


def test_valid_rate_tables_still_parse():
    """The validation rejects only degenerate input, not legitimate tables.

    A zero rate is legitimate (shipped tables floor at 1e-35, but 0 resamples
    finitely) and a 2-column table is the documented minimum.
    """
    T9, rate, err, header = cr.parse_rate_upload(VALID_TABLE)
    assert T9.tolist() == [1e-3, 1e-2, 1e-1, 1.0, 10.0]
    assert np.all(err == 0.0)

    T9, rate, _err, _h = cr.parse_rate_upload("1e-3 0.0\n1e-2 2e-3\n1e-1 3e-3\n")
    assert rate[0] == 0.0


def test_rate_table_error_message_points_at_the_offending_row():
    """The message names the 1-based data row, as a text editor counts it."""
    with pytest.raises(ValueError, match=r"data row 3"):
        cr.parse_rate_upload("1e-3 1e-3\n1e-2 2e-3\n1e-1 -3e-3\n")


# ---------------------------------------------------------------------------
# Zip budget
# ---------------------------------------------------------------------------

def test_zip_bomb_is_refused_before_decompression():
    """A highly compressible archive is refused on its declared sizes alone.

    The measured attack was a 0.10 MB upload expanding to +200 MB of resident
    memory in 0.06 s (685x) -- fatal on the ~1 GB public demo. The refusal
    must come from the central directory, so nothing is decompressed.
    """
    payload = b"1.0 2.0\n" * (8 * 1024 * 1024)   # 64 MB of plausible table text
    bomb = _zip({"networks/b.txt": "n_p__d_g, t.txt\n",
                 "tables/n_p__d_g/t.txt": payload})
    assert len(bomb.getvalue()) < 1_000_000, "test's own bomb should be small"
    with pytest.raises(ValueError, match="expands to"):
        cr.import_zip(bomb)


def test_zip_with_too_many_entries_is_refused():
    """An archive with more members than any real network is refused."""
    many = {"networks/x.txt": "n_p__d_g\n"}
    many.update({f"tables/r{i}/t.txt": "1e-3 1e-3\n" for i in range(cr._ZIP_MAX_ENTRIES + 1)})
    with pytest.raises(ValueError, match="entries"):
        cr.import_zip(_zip(many))


def test_understated_entry_size_surfaces_as_a_clean_error():
    """A central directory that lies about a member's size is refused cleanly.

    ``zipfile`` stops the read at the declared size and then fails the CRC, so
    the lie cannot get past the budget -- but it raises ``BadZipFile``, which
    is not a ``ValueError`` and so bypassed every ``except ValueError`` in the
    dialogs. It must arrive as the same clean message as any other bad upload.
    """
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("networks/x.txt", "n_p__d_g, t.txt\n")
        zf.writestr("tables/n_p__d_g/t.txt", b"1.0 2.0\n" * 700_000)
        for info in zf.infolist():          # lie about the big one
            if info.filename.endswith("t.txt") and info.file_size > 1000:
                info.file_size = 10
    buf.seek(0)
    with pytest.raises(ValueError, match="corrupt"):
        cr.import_zip(buf)


def test_truncated_archive_surfaces_as_a_clean_error():
    """A half-downloaded zip is refused with a message, not a BadZipFile."""
    good = _zip({"networks/x.txt": "n_p__d_g, t.txt\n",
                 "tables/n_p__d_g/t.txt": VALID_TABLE}).getvalue()
    with pytest.raises(ValueError):
        cr.import_zip(io.BytesIO(good[: len(good) // 2]))


# ---------------------------------------------------------------------------
# Manifest vs contents
# ---------------------------------------------------------------------------

def test_manifest_naming_an_absent_table_is_refused():
    """A zip listing a rate table it does not carry must not import silently.

    Left unreported, the run substituted *this install's* shipped table --
    the exact substitution the verbatim, original-grid export exists to
    prevent -- and reported YP = 0.00000000 with no warning.
    """
    with pytest.raises(ValueError, match="not in the archive"):
        cr.import_zip(_zip({"networks/x.txt": "n_p__d_g, missing.txt\n"}))


def test_bare_reaction_lines_need_no_table():
    """A line with no filename is legitimate and must still import.

    Decay reactions have no per-reaction table, and ``export_zip`` falls back
    to a bare name when a shipped table cannot be read -- neither is a
    manifest/contents disagreement.
    """
    result = cr.import_zip(_zip({"networks/x.txt": "n_p__d_g\nBm_H3\n"}))
    assert result["kept"] == ["n_p__d_g", "Bm_H3"]
    assert result["replaced"] == {}


def test_decay_override_line_is_not_mistaken_for_a_missing_table():
    """``name, <number>`` is an overridden decay rate, not a filename."""
    result = cr.import_zip(_zip({"networks/x.txt": "Bm_H3, 1.784e-09\n"}))
    assert result["decay_overrides"]["Bm_H3"] == pytest.approx(1.784e-09)


# ---------------------------------------------------------------------------
# Zip slip
# ---------------------------------------------------------------------------

def test_traversal_member_is_refused():
    """``tables/../evil.txt`` splits to the component '..' and must be refused."""
    with pytest.raises(ValueError, match="unsafe path component"):
        cr.import_zip(_zip({"networks/x.txt": "n_p__d_g, u.txt\n",
                            "tables/../evil.txt": VALID_TABLE,
                            "tables/n_p__d_g/u.txt": VALID_TABLE}))


def test_export_never_writes_an_escaping_member(tmp_path):
    """A hostile 'filenames' entry cannot make primat author a zip-slip archive.

    ``export_zip`` copies that basename straight into the archive path, which
    produced ``tables/n_p__d_g/../../../../evil.txt`` -- escaping on any naive
    extractor.
    """
    from primat.config import PRIMATConfig
    cfg = PRIMATConfig({})
    custom = {"removed": [], "added": {}, "replaced": {"n_p__d_g": VALID_TABLE},
              "filenames": {"n_p__d_g": "../../../../evil.txt"}}
    data = cr.export_zip(cfg, custom, ["n_p__d_g"], network_filename="pwn")
    names = zipfile.ZipFile(io.BytesIO(data)).namelist()
    assert not any(".." in name.split("/") for name in names), names
    # And extracting it stays inside the destination.
    zipfile.ZipFile(io.BytesIO(data)).extractall(tmp_path)
    for path in tmp_path.rglob("*"):
        assert tmp_path in path.parents or path.parent == tmp_path


@pytest.mark.parametrize("candidate,expected", [
    ("../../evil.txt", ".._.._evil.txt"),
    ("..", "fallback.txt"),
    (".", "fallback.txt"),
    ("", "fallback.txt"),
    (None, "fallback.txt"),
    ("normal_primat.txt", "normal_primat.txt"),
])
def test_safe_basename(candidate, expected):
    """Path separators collapse; a bare '.'/'..' falls back entirely."""
    assert cr._safe_basename(candidate, "fallback.txt") == expected


# ---------------------------------------------------------------------------
# Non-UTF-8 members
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("member", ["networks/x.txt", "tables/n_p__d_g/t.txt"])
def test_non_utf8_zip_member_names_the_file(member):
    """A binary member is refused with a message naming it, not a codec dump.

    The raw ``UnicodeDecodeError`` reached the dialog as "'utf-8' codec can't
    decode byte 0xff in position 0", which identifies neither the file nor
    the fix.
    """
    entries = {"networks/x.txt": "n_p__d_g, t.txt\n",
               "tables/n_p__d_g/t.txt": VALID_TABLE}
    entries[member] = b"\xff\xfe\x00 binary"
    with pytest.raises(ValueError, match=f"'{member}' is not UTF-8"):
        cr.import_zip(_zip(entries))


# ---------------------------------------------------------------------------
# Non-regression: the real shipped networks still round-trip
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("network", ["small", "small_parthenope", "large"])
def test_shipped_network_zip_still_round_trips(network):
    """The hardening refuses only hostile archives, not primat's own exports.

    The full ``large`` export -- 4.33 MB compressed, 16.86 MB uncompressed,
    391 entries -- is the widest legitimate archive and must stay comfortably
    inside every budget.
    """
    import re
    from primat.config import PRIMATConfig
    from primat.network_data import load_reaction_names

    cfg = PRIMATConfig({"network": network})
    kept = [re.split(r"[, ]+", e, maxsplit=1)[0].strip()
            for e in load_reaction_names(cfg, network)]
    data = cr.export_zip(cfg, {"removed": [], "replaced": {}, "added": {}},
                         kept, network_filename=network)
    result = cr.import_zip(io.BytesIO(data))
    assert result["kept"] == kept
    assert result["title"] == network
