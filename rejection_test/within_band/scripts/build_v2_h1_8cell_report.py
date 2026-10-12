#!/usr/bin/env python3
"""Build V2_H1_REPORT.md from per-cell JSONs in within_band/Results/raw/.

Pass --live to write a partial-progress report (V2_H1_LIVE_REPORT.md).
"""
from __future__ import annotations
import argparse
import glob
import json
import math
import sys
from datetime import datetime
from pathlib import Path

EXP = Path(__file__).resolve().parent.parent
RAW = EXP / "Results" / "raw"
OUT_DIR = EXP / "Results" / "tables"

PAIRS = [
    ("eta=0.0, chi=0.0 (H0)", "0", "0"),
    ("eta=0.1, chi=0.5",      "0.1", "0.5"),
    ("eta=0.3, chi=0.5",      "0.3", "0.5"),
    ("eta=0.3, chi=1.0",      "0.3", "1"),
    ("eta=0.5, chi=1.0",      "0.5", "1"),
    ("eta=1.0, chi=1.0",      "1", "1"),
]
REGIMES = [
    ("n=2000, p=800, beta=2",   2000, 800,  "2"),
    ("n=2000, p=800, beta=0.9", 2000, 800,  "0.9"),
    ("n=800, p=2000, beta=2",   800,  2000, "2"),
    ("n=800, p=2000, beta=0.9", 800,  2000, "0.9"),
]


def find_cell_json(eta: str, chi: str, n: int, p: int, beta: str) -> Path | None:
    if eta == "0" and chi == "0":
        pat = str(RAW / "**" / f"withinband_ce_n{n}_p{p}_beta{beta}_s010_eta0_chi0_Cdelta*_null_*_w*.json")
    else:
        pat = str(RAW / "**" / f"withinband_ce_n{n}_p{p}_beta{beta}_s010_eta{eta}_chi{chi}_Cdelta*_h1_*.json")
    matches = sorted(glob.glob(pat, recursive=True))
    return Path(matches[-1]) if matches else None


def load_cell(eta: str, chi: str, n: int, p: int, beta: str) -> dict | None:
    p_path = find_cell_json(eta, chi, n, p, beta)
    if p_path is None:
        return None
    with p_path.open() as f:
        return json.load(f)


def fmt(x, nd=3) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "nan"
    return f"{x:.{nd}f}"


def build_table(stat_key: str, nd: int, header_label: str, missing_marker: str = "n/a") -> list[str]:
    lines = [f"### {header_label}", ""]
    head = "| (eta, chi) |" + "|".join(f" {label} " for label, *_ in REGIMES) + "|"
    sep = "|------------|" + "|".join("--------" for _ in REGIMES) + "|"
    lines.append(head)
    lines.append(sep)
    for label, eta, chi in PAIRS:
        row = [f"| {label} "]
        for r_label, n, p, beta in REGIMES:
            cell = load_cell(eta, chi, n, p, beta)
            if cell is None:
                row.append(f"| {missing_marker} ")
                continue
            s = cell.get("summary", {})
            row.append(f"| {fmt(s.get(stat_key, float('nan')), nd)} ")
        row.append("|")
        lines.append("".join(row))
    lines.append("")
    return lines


def cell_status_table() -> list[str]:
    lines = ["### Cell status", ""]
    head = "| (eta, chi) |" + "|".join(f" {label} " for label, *_ in REGIMES) + "|"
    sep = "|------------|" + "|".join("--------" for _ in REGIMES) + "|"
    lines += [head, sep]
    total_done = 0
    total = 0
    for label, eta, chi in PAIRS:
        row = [f"| {label} "]
        for r_label, n, p, beta in REGIMES:
            total += 1
            cell = load_cell(eta, chi, n, p, beta)
            row.append("| done " if cell is not None else "| running ")
            if cell is not None:
                total_done += 1
        row.append("|")
        lines.append("".join(row))
    lines.append("")
    lines.append(f"_Cells finished so far: {total_done}/{total}._")
    lines.append("")
    return lines


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--live", action="store_true")
    args = parser.parse_args()
    if not RAW.is_dir():
        print(f"ERROR: missing {RAW}", file=sys.stderr)
        return 2

    title = "v2 H1 report - live" if args.live else "v2 H1 report"
    missing_marker = "running" if args.live else "n/a"

    lines: list[str] = []
    lines.append(f"# {title}")
    lines.append("")
    lines.append(f"_Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S %Z')}_")
    lines.append("")
    lines.append(
        "Calibration: theta0=1.0, rho=1.0, delta_n_constant=0.05, C_delta=0.5 "
        "(the only combo that kept CE rejection rate in [0.05, 0.15] under H0 "
        "across all four regimes during the v2 null delta-sweep). The H0 "
        "(eta=0, chi=0) row sources from that calibration cell."
    )
    lines.append("")
    if args.live:
        lines += cell_status_table()

    lines += build_table("ce_rejection_rate",                                    3, "CE rejection rate (POWER under H1)", missing_marker)
    lines += build_table("direct_rejection_rate",                                3, "Direct test rejection rate (baseline POWER)", missing_marker)
    lines += build_table("mean_number_pumped_total",                             1, "Mean number of pumps per epoch", missing_marker)
    lines += build_table("mean_T_tilde",                                         3, "Mean T_tilde (CE statistic)", missing_marker)
    lines += build_table("mean_T_hat",                                           3, "Mean T_hat (direct statistic)", missing_marker)
    lines += build_table("mean_q_hat",                                           3, "Mean q_hat (direct-bootstrap critical value)", missing_marker)
    lines += build_table("mean_number_pumped_weakened_edges",                    1, "Mean pumps on weakened edges (H1 signal)", missing_marker)
    lines += build_table("mean_number_pumped_unweakened_edges",                  1, "Mean pumps on unweakened edges (false pumps)", missing_marker)

    lines.append("## Power-difference reading")
    lines.append("")
    lines.append(
        "Compare 'CE rejection rate' to 'Direct test rejection rate' cell-by-cell. "
        "If CE > direct, the consistency-empowerment delivers extra power on this "
        "alternative. The 'mean pumps on weakened edges' table is the cleanest "
        "signal of CE engagement - it is the count of (i,j) pairs where the local "
        "inverse magnitude fell below tau * delta_n on a truly weakened edge, "
        "which is exactly where the pump should help under H1."
    )
    lines.append("")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / ("V2_H1_LIVE_REPORT.md" if args.live else "V2_H1_REPORT.md")
    out.write_text("\n".join(lines))
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
