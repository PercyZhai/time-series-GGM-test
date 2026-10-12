#!/usr/bin/env python3
"""Rebuild the published within-band rejection-rate table from per-cell JSON."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REGIMES = ((2000, 800, 2.0, 160), (2000, 800, 0.9, 238), (800, 2000, 2.0, 86), (800, 2000, 0.9, 86))
PAIRS = ((0.0, 0.0), (0.1, 0.5), (0.3, 1.0), (0.5, 1.0))


def build(raw_dir: Path) -> str:
    cells = {}
    for path in raw_dir.rglob("withinband_ce_*.json"):
        with path.open() as handle:
            result = json.load(handle)
        params = result["params"]
        if (float(result["C_delta"]), float(result["delta_n"])) != (0.5, 0.05):
            continue
        if int(params["seed"]) != 101 or int(params["epochs"]) != 200:
            continue
        regime = (int(params["n"]), int(params["p"]), float(params["beta"]), int(params["block_length"]))
        pair = (float(result["eta"]), float(result["chi"]))
        if regime not in REGIMES or pair not in PAIRS:
            continue
        key = (regime, pair)
        if key in cells:
            raise ValueError(f"Duplicate cell {key}: {path}")
        summary = result["summary"]
        cells[key] = (float(summary["ce_rejection_rate"]), float(summary["direct_rejection_rate"]))
    expected = {(regime, pair) for regime in REGIMES for pair in PAIRS}
    if cells.keys() != expected:
        raise ValueError(f"Missing cells: {expected - cells.keys()}; extra cells: {cells.keys() - expected}")
    lines = ["# Within-band empirical rejection rates", "", "Each cell is CE / direct rejection rate over 200 epochs.", ""]
    lines.append("| (eta, chi) | " + " | ".join(f"n={n}, p={p}, beta={beta:g}" for n, p, beta, _ in REGIMES) + " |")
    lines.append("|---|" + "---:|" * len(REGIMES))
    for pair in PAIRS:
        row = [f"| ({pair[0]:g}, {pair[1]:g})"]
        row += [f"{cells[(regime, pair)][0]:.3f} / {cells[(regime, pair)][1]:.3f}" for regime in REGIMES]
        lines.append(" | ".join(row) + " |")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=ROOT / "Results" / "raw")
    parser.add_argument("--output", type=Path, default=ROOT / "Results" / "tables" / "withinband_paper_table.md")
    args = parser.parse_args()
    report = build(args.raw_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
