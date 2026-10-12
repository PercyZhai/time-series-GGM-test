#!/usr/bin/env python3
"""Build V2_NULL_DELTA_SWEEP_REPORT.md from per-cell JSONs of the v2 null delta-sweep.

Pass --live to write a partial-progress report (V2_NULL_DELTA_SWEEP_LIVE_REPORT.md)
that tolerates missing cells.
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

# Each row is (label, list-of-subdirs-to-search).  Multiple subdirs allow data
# for a single (delta, C) combo to be sourced from per-regime probe runs.
DELTAS = [
    ("d=0.01 C=1.00", ["delta0.01_C1.0"]),
    ("d=0.05 C=0.25", ["delta0.05_C0.25"]),
    ("d=0.05 C=0.50", ["delta0.05_C0.5", "delta0.05_C0.5_n800p2000b2", "delta0.05_C0.5_n2000p800b2"]),
    ("d=0.05 C=1.00", ["delta0.05", "delta0.05_C1.0_n2000p800b2"]),
    ("d=0.10 C=1.00", ["delta0.1", "delta0.1_C1.0_n2000p800b2"]),
    ("d=0.15 C=1.00", ["delta0.15"]),
    ("d=0.20 C=1.00", ["delta0.2"]),
    ("d=0.50 C=1.00", ["delta0.5"]),
    ("d=0.80 C=1.00", ["delta0.8"]),
    ("d=1.20 C=1.00", ["delta1.2"]),
]
REGIMES = [
    ("n=2000, p=800, beta=2",   2000, 800,  "2"),
    ("n=2000, p=800, beta=0.9", 2000, 800,  "0.9"),
    ("n=800, p=2000, beta=2",   800,  2000, "2"),
    ("n=800, p=2000, beta=0.9", 800,  2000, "0.9"),
]
ALPHA = 0.10
PASS_LO, PASS_HI = 0.05, 0.15

GAUSSIAN_TAILS = {
    "0.2": 0.16,
    "0.5": 0.38,
    "0.8": 0.58,
    "1.2": 0.77,
}


def find_cell_json(d_subdirs, n: int, p: int, beta: str) -> Path | None:
    if isinstance(d_subdirs, str):
        d_subdirs = [d_subdirs]
    for d in d_subdirs:
        pat = str(RAW / "**" / d / f"withinband_ce_n{n}_p{p}_beta{beta}_s010_eta0_chi0_Cdelta*_null_*_w*.json")
        matches = sorted(glob.glob(pat, recursive=True))
        if matches:
            return Path(matches[-1])
    return None


def load_cell(d_subdirs, n: int, p: int, beta: str) -> dict | None:
    p_path = find_cell_json(d_subdirs, n, p, beta)
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
    head = "| delta_n |" + "|".join(f" {label} " for label, *_ in REGIMES) + "|"
    sep = "|---------|" + "|".join("--------" for _ in REGIMES) + "|"
    lines.append(head)
    lines.append(sep)
    for d_label, d_subdir in DELTAS:
        row = [f"| {d_label} "]
        for label, n, p, beta in REGIMES:
            cell = load_cell(d_subdir, n, p, beta)
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
    head = "| delta_n |" + "|".join(f" {label} " for label, *_ in REGIMES) + "|"
    sep = "|---------|" + "|".join("--------" for _ in REGIMES) + "|"
    lines += [head, sep]
    total_done = 0
    total = 0
    for d_label, d_subdir in DELTAS:
        row = [f"| {d_label} "]
        for label, n, p, beta in REGIMES:
            total += 1
            cell = load_cell(d_subdir, n, p, beta)
            row.append("| done " if cell is not None else "| running ")
            if cell is not None:
                total_done += 1
        row.append("|")
        lines.append("".join(row))
    lines.append("")
    lines.append(f"_Cells finished so far: {total_done}/{total}._")
    lines.append("")
    return lines


def passing_deltas() -> list[str]:
    out = []
    for d_label, d_subdir in DELTAS:
        ok = True
        for label, n, p, beta in REGIMES:
            cell = load_cell(d_subdir, n, p, beta)
            if cell is None:
                ok = False
                break
            s = cell.get("summary", {})
            ce = s.get("ce_rejection_rate", float("nan"))
            if not (PASS_LO <= ce <= PASS_HI):
                ok = False
                break
        if ok:
            out.append(d_label)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--live", action="store_true")
    args = parser.parse_args()
    if not RAW.is_dir():
        print(f"ERROR: missing {RAW}", file=sys.stderr)
        return 2

    winners = passing_deltas()
    largest = winners[-1] if winners else None
    smallest = winners[0] if winners else None
    complete = sum(
        1 for d_label, d_subdir in DELTAS for label, n, p, beta in REGIMES
        if load_cell(d_subdir, n, p, beta) is not None
    )
    all_done = complete == len(DELTAS) * len(REGIMES)
    title = "v2 null delta-sweep — live report" if args.live else "v2 null delta-sweep report"
    missing_marker = "running" if args.live else "n/a"

    lines: list[str] = []
    lines.append(f"# {title}")
    lines.append("")
    lines.append(f"_Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S %Z')}_")
    lines.append("")
    if args.live and not all_done:
        if smallest is not None:
            lines.append(
                f"**Provisional verdict (partial data):** delta_n in [{PASS_LO:.2f}, {PASS_HI:.2f}] "
                f"across all completed regimes: {', '.join(winners)} (smallest = {smallest}, largest = {largest})."
            )
        else:
            lines.append("**Verdict:** pending; no delta_n has cleared all completed regimes yet.")
    else:
        if smallest is not None:
            lines.append(
                f"**Verdict:** smallest delta_n in [{PASS_LO:.2f}, {PASS_HI:.2f}] across **all four regimes**: "
                f"**delta_n = {smallest}** (passing: {', '.join(winners)}; largest passing = {largest})."
            )
        else:
            lines.append(
                f"**Verdict:** no delta_n in {{0.2, 0.5, 0.8, 1.2}} clears alpha = {ALPHA:.2f} "
                f"(target band [{PASS_LO:.2f}, {PASS_HI:.2f}]) at all four regimes."
            )
    lines.append("")
    lines.append(
        "DGP: within-band weakening CE, v2. theta0 = 1.0, rho = 1.0, eta = 0, chi = 0, "
        "epochs = 200, alpha = 0.10, C_n = sqrt(log p), C_delta = 1, delta0 = 0.10, "
        "b0 = 4, s0 = 10. Block lengths 160 / 238 / 86 / 86. delta_n constant per row."
    )
    lines.append("")

    if args.live:
        lines += cell_status_table()

    lines += build_table("ce_rejection_rate",                                    3, "CE rejection rate (target ~ alpha = 0.10)", missing_marker)
    lines += build_table("mean_number_pumped_total",                             1, "Mean number of pumps per epoch (CE Delta_hat firings)", missing_marker)
    lines += build_table("mean_T_tilde",                                         3, "Mean T_tilde (CE statistic)", missing_marker)
    lines += build_table("mean_T_hat",                                           3, "Mean T_hat (direct statistic)", missing_marker)
    lines += build_table("mean_q_hat",                                           3, "Mean q_hat (direct-bootstrap critical value)", missing_marker)
    lines += build_table("mean_min_ratio_abs_omega_over_tau_true_nonzero",       3, "Mean min |omega_hat_0| / tau_hat over true non-zero coords", missing_marker)
    lines += build_table("mean_median_ratio_abs_omega_over_tau_true_nonzero",    3, "Mean median |omega_hat_0| / tau_hat over true non-zero coords", missing_marker)

    lines.append("## Size / power tradeoff for delta_n")
    lines.append("")
    lines.append(
        "Smaller delta_n lowers the H0 pump count (and so the H0 false-positive "
        "rate), but it also shrinks the chance that a truly deleted edge will "
        "fire the pump under H1. For a true zero, |omega_hat|/tau is "
        "approximately |Z| with Z standard normal, so "
        "P(pump fires) = P(|Z| <= delta_n). Reference values: "
        f"P(|Z| <= 0.2) = {GAUSSIAN_TAILS['0.2']:.2f}, "
        f"P(|Z| <= 0.5) = {GAUSSIAN_TAILS['0.5']:.2f}, "
        f"P(|Z| <= 0.8) = {GAUSSIAN_TAILS['0.8']:.2f}, "
        f"P(|Z| <= 1.2) = {GAUSSIAN_TAILS['1.2']:.2f}. "
        "Pick the largest delta that still keeps CE rejection rate near alpha "
        "under H0 -- that maximises the H1 pump probability while preserving "
        "the type-I error."
    )
    lines.append("")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / ("V2_NULL_DELTA_SWEEP_LIVE_REPORT.md" if args.live else "V2_NULL_DELTA_SWEEP_REPORT.md")
    out.write_text("\n".join(lines))
    print(f"Wrote {out}")
    if smallest is not None:
        print(f"Passing deltas: {winners}; largest = {largest}; smallest = {smallest}")
    else:
        print("No delta_n passes the null sanity yet.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
