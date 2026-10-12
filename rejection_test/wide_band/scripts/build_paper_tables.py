#!/usr/bin/env python3
"""Rebuild the four published kappa tables from per-cell JSON and epoch CSV files."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
KAPPAS = (0.03, 0.06, 0.09, 0.12, 0.15, 0.18)
TABLES = (
    ("table1_n2000_p800_beta09_l238_m80_200_400_800_1200", 2000, 800, 0.9, 238, (80, 200, 400, 800, 1200)),
    ("table2_n2000_p800_beta2_l160_m80_200_400_800_1200", 2000, 800, 2.0, 160, (80, 200, 400, 800, 1200)),
    ("table3_n800_p2000_beta2_l86_m200_500_1000_2000_3000", 800, 2000, 2.0, 86, (200, 500, 1000, 2000, 3000)),
    ("table4_n800_p2000_beta09_l86_m200_500_1000_2000_3000", 800, 2000, 0.9, 86, (200, 500, 1000, 2000, 3000)),
)


def cell_result(path: Path) -> tuple[tuple[int, float], str]:
    with path.open() as handle:
        result = json.load(handle)
    params = result["params"]
    epoch_path = path.parent / Path(result["epoch_results_csv"]).name
    with epoch_path.open(newline="") as handle:
        rejects = [row["reject"].lower() == "true" for row in csv.DictReader(handle)]
    epochs = int(params["epochs"])
    if len(rejects) != epochs or epochs < 2:
        raise ValueError(f"Expected {epochs} epochs in {epoch_path}, found {len(rejects)}")
    rate = sum(rejects) / epochs
    if not math.isclose(rate, float(result["summary"]["rejection_rate"]), abs_tol=1e-12):
        raise ValueError(f"JSON and epoch CSV disagree: {path}")
    # Sample standard deviation of the reject indicators divided by sqrt(epochs).
    mc_sd = math.sqrt(sum((float(x) - rate) ** 2 for x in rejects) / (epochs - 1)) / math.sqrt(epochs)
    return (int(params["m_extra"]), float(params["kappa"])), f"{rate:.3f} ({mc_sd:.3f})"


def build(raw_dir: Path) -> str:
    lines = [
        "# Kappa rejection-rate tables", "",
        "Each cell is rejection rate (Monte Carlo standard error), computed as the sample SD of epoch-level reject indicators divided by sqrt(epochs).", "",
    ]
    for table_no, (folder, n, p, beta, block_length, m_values) in enumerate(TABLES, 1):
        source = raw_dir / folder
        if not source.is_dir():
            raise FileNotFoundError(source)
        cells: dict[tuple[int, float], str] = {}
        files = list(source.glob("gpu_tgg_*.json"))
        for path in files:
            with path.open() as handle:
                params = json.load(handle)["params"]
            expected = (int(params["n"]), int(params["p"]), float(params["beta"]), int(params["block_length"]))
            if expected != (n, p, beta, block_length):
                raise ValueError(f"Unexpected parameters in {path}: {expected}")
            key, value = cell_result(path)
            if key in cells:
                raise ValueError(f"Duplicate cell {key} in {source}")
            cells[key] = value
        expected_cells = {(m, k) for m in m_values for k in KAPPAS}
        if cells.keys() != expected_cells:
            raise ValueError(f"{source}: missing={expected_cells - cells.keys()}, extra={cells.keys() - expected_cells}")
        lines += [
            f"## Table {table_no}. n={n}, p={p}, beta={beta:g}, block length={block_length}", "",
            "| m / kappa | " + " | ".join(f"{k:.2f}" for k in KAPPAS) + " |",
            "|---|" + "---:|" * len(KAPPAS),
        ]
        for m in m_values:
            lines.append(f"| {m} | " + " | ".join(cells[(m, k)] for k in KAPPAS) + " |")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=ROOT / "Results" / "raw")
    parser.add_argument("--output", type=Path, default=ROOT / "Results" / "tables" / "kappa_paper_tables.md")
    args = parser.parse_args()
    report = build(args.raw_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
