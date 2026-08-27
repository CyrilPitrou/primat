"""Dump what each CLI tells the user, for every condition worth comparing.

Usage:  python tests/message_probes.py [--timeout S] > probes.txt

A helper, not a test. ``test_verbose_parity.py`` pins the two backends' verbose
narration line for line; this covers everything else a user is told -- rejected
settings, warnings, ``--help``, the JSON and plain-text reports -- which is too
open-ended to assert but goes stale silently. Run it before and after a change
that touches messages and diff the two dumps.

Each probe is a command line both CLIs accept. Exit status, stdout and stderr
are recorded separately, so stream discipline (errors on stderr, results on
stdout) and wording parity can both be read off one file. Probes that would run
a full solve are cut short by the timeout, by which point the message under test
has already been printed.

Needs the standalone C CLI built (``make`` in ``primat-c/``).
"""
import argparse, subprocess, sys, os

PY_CLI = [sys.executable, "-m", "primat.cli"]
C_CLI = ["primat-c/build/primat-c"]

PROBES = [
    ("unknown parameter key",        ["--set", "nosuchparam=1"]),
    ("wrong type for a float",       ["--set", "numerical_precision=abc"]),
    ("wrong type for a bool",        ["--set", "QED_corrections=maybe"]),
    ("unknown network name",         ["--network", "nosuchnet"]),
    ("amax = 0",                     ["--amax", "0"]),
    ("amax negative",                ["--amax", "-3"]),
    ("fEDE out of range",            ["--set", "fEDE=1.5"]),
    ("wnEDE out of range",           ["--set", "fEDE=0.1", "--set", "wnEDE=0.1"]),
    ("DeltaNeff below -3",           ["--set", "DeltaNeff=-5"]),
    ("T9 grid bounds reversed",      ["--set", "rate_grid_T9_min=20"]),
    ("T_end above T_start",          ["--set", "T_end_MeV=100"]),
    ("mn - mp below me",             ["--set", "mn=938.5"]),
    ("data_dir missing",             ["--data_dir", "/nonexistent-primat"]),
    ("user_nuclear_dir missing",     ["--user_nuclear_dir", "/nonexistent-primat"]),
    ("nevo_file missing",            ["--set", "nevo_file=nope.csv"]),
    ("nevo_file_prefix missing",     ["--set", "nevo_file_prefix=NOPE"]),
    ("mc sample count zero",         ["--mc", "0"]),
    ("mc sample count not a number", ["--mc", "abc"]),
    ("external_scale_factor alone",  ["--set", "external_scale_factor=True",
                                      "--set", "incomplete_decoupling=False"]),
    ("unknown reaction variation",   ["--set", "p_no_such__reaction=1"]),
    ("unknown flag",                 ["--bogus"]),
    ("cache info",                   ["--cache-info"]),
    ("list reactions",               ["--list-reactions"]),
    ("verbose json run",             ["--json", "--verbose"]),
    ("json run",                     ["--json"]),
    ("coarse rate grid warning",     ["--set", "rate_grid_npts=2"]),
    ("loose tolerance warning",      ["--set", "numerical_precision=1e-3"]),
    ("decay era without large",      ["--set", "decay_era=True"]),
    ("chemical potential warning",   ["--set", "munuOverTnu=0.1"]),
    ("amax drops helium",            ["--amax", "3"]),
]


def run(cmd, timeout):
    if not os.path.exists(cmd[0]) and os.sep in cmd[0]:
        return "not built", "", f"{cmd[0]} does not exist -- run `make` in primat-c/"
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        return "timeout", (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""), \
               (e.stderr or b"").decode() if isinstance(e.stderr, bytes) else (e.stderr or "")


def show(tag, rc, out, err, head):
    print(f"[{tag}] exit={rc}")
    for stream, text in (("stdout", out), ("stderr", err)):
        lines = [l for l in text.splitlines() if l.strip()]
        if not lines:
            print(f"[{tag}] {stream}: <empty>")
            continue
        for line in lines[:head]:
            print(f"[{tag}] {stream}: {line[:400]}")
        if len(lines) > head:
            print(f"[{tag}] {stream}: ... ({len(lines) - head} more lines)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timeout", type=float, default=25.0)
    ap.add_argument("--head", type=int, default=8)
    ap.add_argument("--only", default=None)
    args = ap.parse_args()

    if not os.path.exists(C_CLI[0]):
        print(f"# {C_CLI[0]} is not built; run `make` in primat-c/ to compare "
              "both CLIs. Showing the Python side only.\n")

    for label, probe in PROBES:
        if args.only and args.only not in label:
            continue
        print("=" * 78)
        print(f"### {label}:  {' '.join(probe)}")
        for tag, cli in (("py", PY_CLI + ["--backend", "python"]),
                         ("c", C_CLI)):
            rc, out, err = run(cli + probe, args.timeout)
            show(tag, rc, out, err, args.head)
        print()


if __name__ == "__main__":
    main()
