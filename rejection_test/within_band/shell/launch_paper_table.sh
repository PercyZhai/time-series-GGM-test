#!/usr/bin/env bash
set -euo pipefail

EXP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="$(command -v "${TGG_PYTHON:-python3}")"
OUT_DIR="${EXP}/Results/raw/paper_table"
mkdir -p "${OUT_DIR}"

# Four reported regimes and four reported (eta, chi) pairs: 16 cells.
for regime in 2000:800:2:160 2000:800:0.9:238 800:2000:2:86 800:2000:0.9:86; do
  for pair in 0:0 0.1:0.5 0.3:1.0 0.5:1.0; do
    IFS=: read -r eta chi <<< "${pair}"
    "${PYTHON}" "${EXP}/scripts/run_withinband_weakening_ce_v2_gpu.py" \
      --regime_grid "${regime}" --s0 10 --theta0 1.0 --rho 1.0 --delta0 0.10 \
      --epochs 200 --alpha 0.10 --seed 101 --batch_epochs 5 \
      --eta_grid "${eta}" --chi_grid "${chi}" \
      --C_delta 0.5 --edge_seed 202 --gpu_chunk_cols 4 \
      --delta_n_constant 0.05 \
      --run_label "paper_${regime//:/_}_eta${eta}_chi${chi}" \
      --out_dir "${OUT_DIR}"
  done
done

"${PYTHON}" "${EXP}/scripts/build_paper_table.py"
