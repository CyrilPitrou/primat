#!/usr/bin/env python3
"""Report the source coverage each fuzz target reached.

Replays a target's evolved corpus through the -fprofile-instr-generate build
(``make fuzz-coverage``) and prints llvm-cov's line and region coverage for
the source files that target is aimed at -- what was exercised, rather than
how long the run took.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
COV = ROOT / "build" / "fuzzcov"
DATA_DIR = ROOT.parent / "primat" / "data"

# The files each target is aimed at, so the report is about the parser under
# test rather than the whole binary it happens to be linked into.
SCOPE = {
    "fuzz_rate_table": ["src/table_io.c"],
    "fuzz_nevo": ["src/table_io.c", "src/neutrino_history.c"],
    "fuzz_nevo_spectral": ["src/table_io.c", "src/neutrino_history.c"],
    "fuzz_ini": ["src/ini.c", "src/config.c"],
    "fuzz_network_list": ["src/network_data.c"],
    "fuzz_decays": ["src/network_data.c"],
    "fuzz_detailed_balance": ["src/network_data.c"],
    "fuzz_reactions_large": ["src/network_data.c"],
    "fuzz_nuclides_csv": ["src/config.c"],
    "fuzz_cache": ["src/cache.c", "src/table_io.c"],
    "fuzz_cli": ["src/cli.c", "src/config.c"],
}

# Functions that are the actual entry point under test, reported separately:
# a whole-file percentage understates a target that saturates its own parser.
FUNCS = {
    "fuzz_rate_table": ["cpr_table_read", "split_fields", "cpr_validate_rate_table"],
    "fuzz_nevo": ["cpr_neutrino_history_init", "build_nevo_table"],
    "fuzz_nevo_spectral": ["cpr_neutrino_history_init", "build_nevo_table"],
    "fuzz_ini": ["cpr_ini_load", "cpr_parse_literal", "cpr_config_set_by_name"],
    "fuzz_network_list": ["cpr_load_network_list"],
    "fuzz_decays": ["cpr_load_decays"],
    "fuzz_detailed_balance": ["cpr_load_detailed_balance", "csv_split"],
    "fuzz_reactions_large": ["cpr_load_reactions_large", "csv_split"],
    "fuzz_nuclides_csv": ["load_nuclides"],
    "fuzz_cache": ["cpr_cache_read_fingerprint_hash", "cpr_table_read"],
    "fuzz_cli": ["cpr_cli_main", "apply_param"],
}


def xcrun(tool: str, *args: str) -> str:
    return subprocess.run(["xcrun", tool, *args], capture_output=True,
                          text=True, check=True).stdout


def main() -> int:
    rows = []
    for target, scope in SCOPE.items():
        binary = COV / target
        corpus = HERE / "corpus" / target.replace("fuzz_", "")
        if not binary.exists():
            print(f"{target}: not built (make fuzz-coverage)", file=sys.stderr)
            continue
        raw = COV / f"{target}.profraw"
        env = dict(os.environ)
        env["CPRIMAT_FUZZ_DATA_DIR"] = str(DATA_DIR)
        env["LLVM_PROFILE_FILE"] = str(raw)
        subprocess.run([str(binary), str(corpus), "-runs=0"], env=env, cwd=ROOT,
                       capture_output=True)
        data = COV / f"{target}.profdata"
        xcrun("llvm-profdata", "merge", "-sparse", str(raw), "-o", str(data))
        report = json.loads(xcrun("llvm-cov", "export", str(binary),
                                  f"-instr-profile={data}", "-summary-only"))
        files = {f["filename"]: f["summary"] for f in report["data"][0]["files"]}
        for src in scope:
            key = next((k for k in files if k.endswith(src)), None)
            if key is None:
                continue
            s = files[key]
            rows.append((target, src, s["lines"]["percent"], s["regions"]["percent"],
                         s["lines"]["covered"], s["lines"]["count"]))
        # Per-function figures for the entry points themselves.
        fn_report = json.loads(xcrun("llvm-cov", "export", str(binary),
                                     f"-instr-profile={data}"))
        by_name = {}
        for fn in fn_report["data"][0]["functions"]:
            # llvm-cov region tuple: [l0, c0, l1, c1, count, file_id, ...]
            covered = sum(1 for r in fn["regions"] if r[4] > 0)
            # A file-static function is exported as "path/to/file.c:name".
            short = fn["name"].rsplit(":", 1)[-1].lstrip("_")
            by_name[short] = (covered, len(fn["regions"]))
        for name in FUNCS.get(target, []):
            if name in by_name:
                c, t = by_name[name]
                rows.append((target, f"  {name}()", 100.0 * c / t if t else 0.0,
                             100.0 * c / t if t else 0.0, c, t))

    print(f"{'target':22} {'scope':34} {'lines%':>7} {'regions%':>9}  covered/total")
    for target, src, lp, rp, cov, tot in rows:
        print(f"{target:22} {src:34} {lp:7.1f} {rp:9.1f}  {cov}/{tot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
