#!/usr/bin/env python3
"""Seed each fuzz corpus from the shipped data of the matching kind.

Every corpus also gets the degenerate shapes passes 20 and 21 collected by
hand -- empty, comments-only, unsorted, negative, NaN, duplicate, truncated --
so the fuzzer starts from the boundary rather than having to rediscover it.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent.parent / "primat" / "data"
EXAMPLES = HERE.parent.parent / "examples"


def head(path: Path, n: int) -> bytes:
    return b"".join(path.open("rb").readlines()[:n])


DEGENERATE = {
    "empty": b"",
    "comments_only": b"# only a comment\n# and another\n",
    "blank_lines": b"\n\n   \n\t\n",
    "unsorted": b"2.0 1.0 0.1\n1.0 2.0 0.1\n3.0 0.5 0.1\n",
    "negative": b"1.0 -1.0 0.1\n2.0 -2.0 0.1\n",
    "nonfinite": b"1.0 nan 0.1\n2.0 inf 0.1\n3.0 -inf 0.1\n",
    "duplicate": b"1.0 1.0 0.1\n1.0 1.0 0.1\n",
    "single_row": b"1.0 1.0 0.1\n",
    "truncated": b"1.0 1.0 0.1\n2.0 2.0",
    "ragged": b"1.0 1.0 0.1\n2.0 2.0\n3.0 3.0 0.3 0.4\n",
    "overflow": b"1e999 1e999 1e999\n-1e999 2 3\n",
    "no_newline_eof": b"1.0 1.0 0.1",
    "crlf": b"1.0 1.0 0.1\r\n2.0 2.0 0.2\r\n",
    "long_line": b"1.0 " + b"9" * 9000 + b" 0.1\n",
    "commas": b"1.0,1.0,0.1\n2.0,2.0,0.2\n",
    "nul": b"1.0 1.0\x000.1\n",
}


def seeds() -> dict[str, dict[str, bytes]]:
    rate = DATA / "nuclear/tables/n_p__d_g/n_p__d_g_primat.txt"
    out: dict[str, dict[str, bytes]] = {}

    out["rate_table"] = {"shipped": head(rate, 40)}
    out["nevo"] = {"shipped": head(DATA / "NEVO/NEVOPRIMAT_col_1_7.csv", 40)}
    out["nevo_spectral"] = {
        "shipped": head(DATA / "NEVO/NEVOPRIMAT.csv", 6) + b"%%%"
                   + head(DATA / "NEVO/NEVOGrid.csv", 200),
        "width_mismatch": b"1 1 1 1 1 1 2 3\n2 1 1 1 1 1 2 3\n%%%1\n2\n3\n4\n5\n",
    }
    out["ini"] = {
        "shipped": (EXAMPLES / "run_basic.ini").read_bytes()[:4000],
        "small": (EXAMPLES / "run_small.ini").read_bytes(),
        "shapes": b"network = small\namax=8\nverbose true\nOmegabh2 = 0.022\n"
                  b"; comment\n# comment\nbad_key = 1\nnetwork =\n",
    }
    out["network_list"] = {
        "shipped": (DATA / "nuclear/networks/small.txt").read_bytes(),
        "bare": b"n_decay\nn_p__d_g, n_p__d_g_primat.txt\n",
        "dup": b"n_p__d_g, a.txt\nn_p__d_g, b.txt\n",
    }
    out["decays"] = {"shipped": head(DATA / "nuclear/tables/decays.txt", 40)}
    out["detailed_balance"] = {"shipped": head(DATA / "csv/detailed_balance.csv", 40)}
    out["reactions_large"] = {"shipped": head(DATA / "csv/reactions_large.csv", 40)}
    out["nuclides_csv"] = {"shipped": head(DATA / "csv/nuclides.csv", 40)}

    weak = sorted((DATA / "cache_plasma_weak/weak").glob("nTOp_*.txt"))
    out["cache"] = {"shipped": head(weak[0], 40)} if weak else {}
    out["cache"]["header_only"] = (
        b"# T9 forward backward\n# fingerprint_hash: 0123456789abcdef\n"
        b"# fingerprint: {\"a\": 1}\n")

    out["cli"] = {
        "flags": b"\0".join([b"--network", b"small", b"--amax", b"8",
                             b"--verbose", b"--json"]),
        "set": b"\0".join([b"--set", b"Omegabh2=0.022", b"--set", b"network=small"]),
        "consts": b"\0".join([b"--tau_n", b"879.4", b"--no-QED_corrections"]),
        "typo": b"\0".join([b"--netwrok", b"small"]),
        "empty_val": b"\0".join([b"--set", b"network="]),
    }
    return out


def main() -> int:
    per_target = seeds()
    for target, items in per_target.items():
        d = HERE / "corpus" / target
        d.mkdir(parents=True, exist_ok=True)
        payloads = dict(items)
        if target != "cli":
            payloads.update(DEGENERATE)
        else:
            payloads.update({"empty": b"", "nul": b"\0\0\0",
                             "dashes": b"--\0-\0---"})
        for name, blob in payloads.items():
            (d / f"seed_{name}").write_bytes(blob)
        print(f"{target}: {len(payloads)} seeds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
