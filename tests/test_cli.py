"""
Tests for the ``primat`` console-script CLI.

``primat.cli.main()`` is invoked in-process (no subprocess) with an
explicit ``argv`` list, which is exactly what the ``primat`` console
script does at startup -- except where the test is about what reaches the
process's own stdout/stderr, which only a subprocess can observe.  Each
invocation runs one full small-network solve (~1.2 s), so these tests are
marked ``slow``/``solve`` like the other single-solve tests in the "solve"
tier.
"""
import json
import os
import re
import subprocess
import sys

import pytest

from primat.cli import main
from primat.credits import cli_credits_text
from tests.reference_values import (
    DOH_ABS_TOL,
    DOH_REFERENCE,
    NEFF_ABS_TOL,
    NEFF_REFERENCE,
    YPBBN_ABS_TOL,
    YPBBN_REFERENCE,
)

pytestmark = [pytest.mark.slow, pytest.mark.solve]


def test_cli_default_summary(capsys):
    """No flags: default (small-network) run, human-readable summary.

    Parses the printed values rather than matching a literal string, and
    compares against the centralised constants and tolerances in
    tests/reference_values.py, so a routine default-parameter tweak (e.g.
    commit e00f062's rate_grid_npts/sampling_temperature_per_decade bump)
    does not require refreshing a hard-coded pin here.
    """
    rc = main([])
    assert rc == 0
    out = capsys.readouterr().out
    neff = float(re.search(r"Neff\s*=\s*([\d.]+)", out).group(1))
    yp   = float(re.search(r"YP \(BBN\)\s*=\s*([\d.]+)", out).group(1))
    doh  = float(re.search(r"D/H\s*=\s*([\d.eE+-]+)", out).group(1))
    assert neff == pytest.approx(NEFF_REFERENCE, abs=NEFF_ABS_TOL)
    assert yp   == pytest.approx(YPBBN_REFERENCE, abs=YPBBN_ABS_TOL)
    assert doh  == pytest.approx(DOH_REFERENCE, abs=DOH_ABS_TOL)
    assert "Li6/Li7" not in out


def test_cli_json_matches_default_summary(capsys):
    """--json prints the full results dict, parseable and consistent."""
    rc = main(["--json"])
    assert rc == 0
    results = json.loads(capsys.readouterr().out)
    assert results["Neff"]   == pytest.approx(NEFF_REFERENCE, abs=NEFF_ABS_TOL)
    assert results["YPBBN"]  == pytest.approx(YPBBN_REFERENCE, abs=YPBBN_ABS_TOL)
    assert results["DoH"]    == pytest.approx(DOH_REFERENCE, abs=DOH_ABS_TOL)
    assert "Li6oLi7" not in results


def test_cli_omegabh2_override_changes_doh(capsys):
    """--Omegabh2 is forwarded to PRIMATConfig and changes the result."""
    rc = main(["--Omegabh2", "0.024", "--json"])
    assert rc == 0
    results = json.loads(capsys.readouterr().out)
    # A higher baryon density measurably increases D/H away from the
    # Omegabh2=0.02242 reference value above.
    assert results["DoH"] != pytest.approx(2.4349347363779478e-05, rel=1e-6)


def test_cli_network_accepts_any_network_file(capsys):
    """--network accepts any name with a rates/nuclear/networks/<name>.txt
    file, not just 'small'/'small_parthenope'/'large'.

    'small_parthenope' (12-reaction network using Parthenope 3.0 rate
    tables) uses different reaction rates from 'small', so YPBBN differs
    from the default-network reference value above while remaining a
    physically reasonable abundance.
    """
    rc = main(["--network", "small_parthenope", "--json"])
    assert rc == 0
    results = json.loads(capsys.readouterr().out)
    assert 0.24 < results["YPBBN"] < 0.25
    assert results["YPBBN"] != pytest.approx(0.24699534223598402, rel=1e-6)


def test_cli_network_rejects_unknown_name(capsys):
    """An unknown --network name is a *user* error, not a crash.

    PRIMATConfig's ValueError must reach the user as the C CLI reports it --
    one "error: ..." line on stderr and exit code 2 (primat-c/src/cli.c) --
    rather than as a Python traceback, which reads as an internal failure.
    """
    assert main(["--network", "no_such_network"]) == 2
    assert "error: network must be" in capsys.readouterr().err


def test_cli_config_error_traceback_escape_hatch(monkeypatch):
    """PRIMAT_TRACEBACK=1 restores the raw exception for debugging primat."""
    monkeypatch.setenv("PRIMAT_TRACEBACK", "1")
    with pytest.raises(ValueError, match="network must be"):
        main(["--network", "no_such_network"])


def test_cli_network_error_mentions_data_tree(capsys):
    """The missing-network error should point at data/nuclear/networks."""
    assert main(["--network", "no_such_network"]) == 2
    assert re.search(r"data/nuclear/networks", capsys.readouterr().err)


def test_cli_network_error_lists_overlay_candidates(tmp_path, capsys):
    """A custom overlay should be named explicitly in the missing-network error."""
    overlay = tmp_path / "overlay"
    overlay.mkdir()
    expected = overlay / "networks" / "custom.txt"

    assert main(["--set", f"user_nuclear_dir={overlay}", "--network", "custom"]) == 2
    # The error quotes each searched path with repr(), which doubles
    # backslashes on Windows ('C:\\a' -> "'C:\\\\a'"); collapse them so the
    # substring check is separator-portable (no-op on POSIX forward slashes).
    message = capsys.readouterr().err.replace("\\\\", "\\")
    assert str(expected) in message


def test_cli_rejects_non_positive_mc():
    """--mc 0 / --mc -5 must be argparse errors, not a "value +/- 0" report.

    A negative count used to underflow the C sampler's size_t buffer size and
    abort the process with "out of memory (1.8e19 bytes)", while the Python
    backend silently reported every observable with a zero sigma.
    """
    for bad in ("0", "-5"):
        with pytest.raises(SystemExit) as excinfo:
            main(["--mc", bad])
        assert excinfo.value.code == 2


def test_cli_json_with_time_evolution_is_serialisable(capsys, tmp_path):
    """--json alongside --output_time_evolution must still print valid JSON.

    Both backends attach a non-serialisable EvolutionResult under
    results["evolution"]; the JSON writer drops it (the series has its own TSV)
    instead of dying with "Object of type EvolutionResult is not JSON
    serializable".
    """
    assert main(["--json", "--output_time_evolution",
                 "--output_file", str(tmp_path / "evo.tsv")]) == 0
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert "evolution" not in payload
    assert payload["DoH"] > 0
    # The omission is announced rather than silent.
    assert "time-evolution" in captured.err


@pytest.mark.parametrize("backend", ["c", "python"])
def test_cli_json_stdout_is_pure_json_in_a_subprocess(backend, tmp_path):
    """`primat --json … | jq` must work: nothing but JSON reaches stdout.

    GOAL: pin the redirection contract that the in-process test above cannot
    see. ``capsys`` captures Python-level writes only, so the C backend's
    ``fprintf`` to the underlying file descriptor is invisible to it -- the
    "[output] Time-evolution data …" progress line went to stdout on both
    backends and corrupted the JSON document while that test stayed green.
    """
    if backend == "c":
        from primat.backend import HAS_C_BACKEND
        if not HAS_C_BACKEND:
            pytest.skip("primat._primat_c C extension is not built")
    proc = subprocess.run(
        [sys.executable, "-m", "primat.cli", "--json", "--backend", backend,
         "--output_time_evolution", "--output_file", str(tmp_path / "evo.tsv")],
        capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode == 0, proc.stderr
    payload = json.loads(proc.stdout)          # fails if anything else got in
    assert payload["DoH"] > 0
    assert "[output]" in proc.stderr           # announced, on the right stream
    assert "[output]" not in proc.stdout


def test_cli_set_expands_tilde_in_path_values(monkeypatch, tmp_path, capsys):
    """Quoted ``~`` paths passed through ``--set`` should resolve to HOME.

    The CLI forwards ``--set`` values as raw strings, so path parameters must
    normalize home-directory prefixes inside the config layer rather than
    relying on shell expansion.
    """
    monkeypatch.setenv("HOME", str(tmp_path))
    # os.path.expanduser("~") consults HOME on POSIX but USERPROFILE (then
    # HOMEDRIVE+HOMEPATH) on Windows, so pin all of them at tmp_path to make
    # the expansion deterministic cross-platform.
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    monkeypatch.setenv("HOMEDRIVE", os.path.splitdrive(str(tmp_path))[0])
    monkeypatch.setenv("HOMEPATH", os.path.splitdrive(str(tmp_path))[1])
    (tmp_path / "custom").mkdir()

    rc = main(["--set", "user_nuclear_dir=~/custom", "--json"])
    assert rc == 0
    # The [init] overlay note quotes the resolved path with repr(), which
    # doubles backslashes on Windows; collapse them so the substring check
    # matches str(...) (single backslashes). No-op on POSIX forward slashes.
    err = capsys.readouterr().err.replace("\\\\", "\\")
    assert "nuclear networks and rate tables" in err
    assert str((tmp_path / "custom").resolve()) in err


def test_cli_help_shows_named_output_path_flags(capsys):
    """``primat --help`` documents the four output-path flags as basic options.

    These paths are user-facing CLI knobs, so they must appear in the printed
    help instead of being buried under the hidden ``--set KEY=VALUE`` escape
    hatch.
    """
    with pytest.raises(SystemExit) as excinfo:
        main(["--help"])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert "--output_file FILE" in out
    assert "--output_final_file FILE" in out
    assert "--output_background_file FILE" in out
    assert "--output_mc_file_prefix PREFIX" in out
    assert not re.search(r"(?m)^\s+--set\b", out)


def test_cli_list_params_covers_every_default_params_key(capsys):
    """--list-params must print every DEFAULT_PARAMS key with its
    default value, without running a solve. Output is grouped by
    primat.config.PARAM_GROUPS (a "# ---- <group> ----" header line per
    group, plus a blank separator line), so the key lines are a subset of
    the total, not every printed line."""
    from primat.config import DEFAULT_PARAMS, PARAM_GROUPS

    rc = main(["--list-params"])
    assert rc == 0
    out = capsys.readouterr().out
    for key in DEFAULT_PARAMS:
        assert re.search(rf"(?m)^{re.escape(key)}\s*=", out)
    for group in PARAM_GROUPS:
        assert f"# ---- {group} ----" in out
    # A few keys with a known one-line trailing comment in config.py should
    # carry it through verbatim, so the descriptions are not silently dropped.
    assert "Omega_b h^2" in out
    assert "Coulomb + T=0 resummed radiative corrections" in out


def test_cli_help_named_flag_defaults_match_default_params(capsys):
    """Named --help flags that quote a "(PRIMATConfig default: ...)" value
    (Omegabh2, DeltaNeff, network, numerical_precision, munuOverTnu, and
    every BooleanOptionalAction flag) build that string from DEFAULT_PARAMS
    at parser-construction time rather than a hand-typed literal, so it
    cannot silently drift from the real default if config.py changes."""
    from primat.config import DEFAULT_PARAMS

    with pytest.raises(SystemExit):
        main(["--help"])
    # argparse word-wraps --help text to the terminal width, so a literal
    # value can land on its own line; collapse all whitespace runs to a
    # single space before substring-matching.
    out = re.sub(r"\s+", " ", capsys.readouterr().out)
    for key in ("Omegabh2", "DeltaNeff", "numerical_precision", "munuOverTnu",
                "QED_corrections", "nuclear_qed_corrections",
                "radiative_corrections", "finite_mass_corrections",
                "thermal_corrections", "spectral_distortions",
                "output_time_evolution", "output_final_result",
                "output_background_evolution", "output_mc_samples",
                "output_mc_covariance", "output_mc_correlation",
                "show_progress"):
        assert f"PRIMATConfig default: {DEFAULT_PARAMS[key]}" in out, key
    assert f"PRIMATConfig default: {DEFAULT_PARAMS['network']!r}" in out


def test_cli_help_mentions_list_params_next_to_set(capsys):
    """--help's epilog should point power users at --list-params, the
    documented way to discover every --set-able parameter."""
    with pytest.raises(SystemExit) as excinfo:
        main(["--help"])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert "--list-params" in out
    assert not re.search(r"(?m)^\s+--set\b", out)


def test_cli_version_reports_backend_availability(capsys):
    """--version must also state whether the C backend is available, so
    users can tell why --backend c would (not) work without a separate run."""
    from primat.backend import HAS_C_BACKEND

    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert "C backend:" in out
    assert ("available" if HAS_C_BACKEND else "unavailable") in out


def test_cli_credits_prints_short_text(capsys):
    """--credits prints the attribution text without install/run guidance."""
    rc = main(["--credits"])
    assert rc == 0
    out = capsys.readouterr().out
    assert out.rstrip("\n") == cli_credits_text()
    assert "pip install primat" not in out


def test_cli_mc_output_announces_path(capsys, tmp_path):
    """The MC TSV writer must also emit a visible [output] line. With the
    prefix scheme the samples file is ``<prefix>_samples.tsv``."""
    prefix = tmp_path / "mc"
    rc = main(["--mc", "1", "--output_mc_samples",
               "--output_mc_file_prefix", str(prefix)])
    assert rc == 0
    out = capsys.readouterr().err          # [output] lines go to stderr
    assert "[output] MC samples (1 sample) written to" in out
    assert str((tmp_path / "mc_samples.tsv").resolve()) in out


def test_cli_mc_covariance_correlation_files(capsys, tmp_path):
    """--output_mc_covariance/--output_mc_correlation write the two matrix
    files under the shared prefix, each announced with an [output] line."""
    prefix = tmp_path / "mc"
    rc = main(["--mc", "3", "--output_mc_covariance", "--output_mc_correlation",
               "--output_mc_file_prefix", str(prefix)])
    assert rc == 0
    captured = capsys.readouterr()
    out = captured.err                     # [output] lines go to stderr
    assert "[output] MC covariance matrix written to" in out
    assert "[output] MC correlation matrix written to" in out
    assert (tmp_path / "mc_covariance.tsv").exists()
    assert (tmp_path / "mc_correlation.tsv").exists()
    # Header wording (author spec: one '#' line naming the file).
    cov_head = (tmp_path / "mc_covariance.tsv").read_text().splitlines()[0]
    assert cov_head.startswith("# Covariance matrix of the N=3 primat MC samples")
    assert "ddof=1" in cov_head
    # The CLI also prints the 4x4 correlation/covariance matrices of the main
    # products after the value +/- sigma block.
    assert "Correlation matrix (YPBBN, DoH, He3oHe4, Li7oH):" in captured.out
    assert "Covariance matrix (YPBBN, DoH, He3oHe4, Li7oH):" in captured.out


def test_cli_mc_output_file_without_enable_flag_does_not_write(capsys, tmp_path):
    """The filename option alone should not force MC output."""
    prefix = tmp_path / "mc"
    rc = main(["--mc", "1", "--output_mc_file_prefix", str(prefix)])
    assert rc == 0
    out = capsys.readouterr().err
    assert "[output] MC samples" not in out
    assert not (tmp_path / "mc_samples.tsv").exists()


def test_cli_mc_summary_includes_all_displayed_sigmas(capsys):
    """The human-readable MC summary should print sigma for every displayed
    ratio, not only the first few observables.

    This exercises a network that actually produces Li6/Li7 and CNO so the
    optional lines are present in the output.
    """
    rc = main(["--network", "large", "--mc", "1"])
    assert rc == 0
    out = capsys.readouterr().out
    assert re.search(r"He3/He4\s*=\s*[\d.eE+-]+\s+\+/-\s+[\d.eE+-]+", out)
    assert re.search(r"Li6/Li7\s*=\s*[\d.eE+-]+\s+\+/-\s+[\d.eE+-]+", out)
    assert re.search(r"CNO \(mass\)\s*=\s*[\d.eE+-]+\s+\+/-\s+[\d.eE+-]+", out)


def test_cli_unknown_rate_variation_key_is_fatal_under_strict(capsys):
    """A mistyped p_<reaction> exits 2 under strict_params, like any config error.

    Reaction names are long and underscore-heavy, so a typo is the archetypal
    silent no-op: the run otherwise proceeds with that rate unvaried. Both CLIs
    report it with the same wording and the same exit code -- the C side's
    check lives in primat-c/src/cli.c, since runs arriving through
    primat/backend.py are already validated by PRIMATConfig.
    """
    assert main(["--set", "strict_params=True", "--set", "p_not_a_reaction=1.0"]) == 2
    assert "does not match any reaction" in capsys.readouterr().err
