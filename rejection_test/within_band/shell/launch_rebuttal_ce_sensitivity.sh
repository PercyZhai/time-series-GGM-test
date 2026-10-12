#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXP="${ROOT}/within_band"
PYTHON_SCRIPT="${EXP}/scripts/run_withinband_weakening_ce_v2_gpu.py"
OUT_BASE="${EXP}/Results/raw/withinband_ce_v2_rebuttal_sensitivity_seed303"
LOG_DIR="${OUT_BASE}/logs"
SCRIPT_REMOTE="${EXP}/shell/launch_rebuttal_ce_sensitivity.sh"

REGIME_GRID="2000:800:2:160,2000:800:0.9:238,800:2000:2:86,800:2000:0.9:86"
DELTAS=(0.01 0.01 0.01 0.05 0.05 0.05 0.10 0.10 0.10)
CDELTAS=(0.25 0.50 1.00 0.25 0.50 1.00 0.25 0.50 1.00)
N_WORKERS=4
SEED=303

mkdir -p "${OUT_BASE}" "${LOG_DIR}"

run_cell() {
  local delta="$1"
  local cdelta="$2"
  local state="$3"
  local eta="$4"
  local chi="$5"
  local tag="delta${delta}_Cdelta${cdelta}_${state}"
  local out_dir="${OUT_BASE}/${tag}"

  mkdir -p "${out_dir}"
  python3 "${PYTHON_SCRIPT}" \
    --regime_grid "${REGIME_GRID}" \
    --s0 10 \
    --theta0 1.0 \
    --rho 1.0 \
    --delta0 0.10 \
    --epochs 200 \
    --alpha 0.10 \
    --seed "${SEED}" \
    --batch_epochs 5 \
    --eta_grid "${eta}" \
    --chi_grid "${chi}" \
    --C_delta "${cdelta}" \
    --edge_seed 202 \
    --gpu_chunk_cols 4 \
    --combo_offset 0 \
    --combo_stride 1 \
    --delta_n_constant "${delta}" \
    --run_label "rebuttal_seed${SEED}_${tag}" \
    --out_dir "${out_dir}"
}

run_worker() {
  local worker_idx="$1"
  local pair_idx

  if [ -f "${ROOT}/.venv_gpu/bin/activate" ]; then source "${ROOT}/.venv_gpu/bin/activate"; fi
  cd "${EXP}"

  for pair_idx in "${!DELTAS[@]}"; do
    if (( pair_idx % N_WORKERS != worker_idx )); then
      continue
    fi

    local delta="${DELTAS[$pair_idx]}"
    local cdelta="${CDELTAS[$pair_idx]}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] worker=${worker_idx} pair=${pair_idx} delta=${delta} C_delta=${cdelta} H0"
    run_cell "${delta}" "${cdelta}" "H0" "0" "0"
    echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] worker=${worker_idx} pair=${pair_idx} delta=${delta} C_delta=${cdelta} H1"
    run_cell "${delta}" "${cdelta}" "H1_eta0.3_chi1.0" "0.3" "1.0"
  done

  echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] worker=${worker_idx} complete"
}

launch_workers() {
  local worker_idx

  for worker_idx in $(seq 0 $((N_WORKERS - 1))); do
    local session="ce_rebuttal_sensitivity_w${worker_idx}"
    local log_file="${LOG_DIR}/${session}.log"

    if tmux has-session -t "${session}" 2>/dev/null; then
      echo "${session} already exists; leaving it unchanged."
      continue
    fi

    tmux new-session -d -s "${session}" \
      "bash -lc '${SCRIPT_REMOTE} worker ${worker_idx}' > '${log_file}' 2>&1"
    echo "Started ${session}; log=${log_file}"
  done
}

case "${1:-launch}" in
  worker)
    run_worker "$2"
    ;;
  launch)
    launch_workers
    ;;
  *)
    echo "Usage: $0 [launch|worker INDEX]" >&2
    exit 2
    ;;
esac
