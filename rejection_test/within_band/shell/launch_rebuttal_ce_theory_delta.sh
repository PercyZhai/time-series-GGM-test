#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXP="${ROOT}/within_band"
PYTHON_SCRIPT="${EXP}/scripts/run_withinband_weakening_ce_v2_gpu.py"
OUT_BASE="${EXP}/Results/raw/withinband_ce_v2_rebuttal_theory_delta_seed303"
LOG_DIR="${OUT_BASE}/logs"
SCRIPT_REMOTE="${EXP}/shell/launch_rebuttal_ce_theory_delta.sh"
WAIT_SESSIONS=(
  ce_rebuttal_sensitivity_w0
  ce_rebuttal_sensitivity_w1
  ce_rebuttal_sensitivity_w2
  ce_rebuttal_sensitivity_w3
)
REGIMES=(
  "2000:800:2:160"
  "2000:800:0.9:238"
  "800:2000:2:86"
  "800:2000:0.9:86"
)
DELTA_VALUES=(
  "0.055139468476009"
  "0.055139468476009"
  "0.051709232164742"
  "0.051709232164742"
)
C_DELTA="0.50"
SEED=303

mkdir -p "${OUT_BASE}" "${LOG_DIR}"

run_state() {
  local worker_idx="$1"
  local state="$2"
  local eta="$3"
  local chi="$4"
  local regime="${REGIMES[$worker_idx]}"
  local delta="${DELTA_VALUES[$worker_idx]}"
  local tag="c0.02sqrtlogn_Cdelta${C_DELTA}_regime${worker_idx}_${state}"
  local out_dir="${OUT_BASE}/${tag}"

  mkdir -p "${out_dir}"
  python3 "${PYTHON_SCRIPT}" \
    --regime_grid "${regime}" \
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
    --C_delta "${C_DELTA}" \
    --edge_seed 202 \
    --gpu_chunk_cols 4 \
    --combo_offset 0 \
    --combo_stride 1 \
    --delta_n_constant "${delta}" \
    --run_label "rebuttal_theory_delta_seed${SEED}_${tag}" \
    --out_dir "${out_dir}"
}

run_worker() {
  local worker_idx="$1"

  if [ -f "${ROOT}/.venv_gpu/bin/activate" ]; then source "${ROOT}/.venv_gpu/bin/activate"; fi
  cd "${EXP}"

  echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] worker=${worker_idx} H0"
  run_state "${worker_idx}" "H0" "0" "0"
  echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] worker=${worker_idx} H1"
  run_state "${worker_idx}" "H1_eta0.3_chi1.0" "0.3" "1.0"
  echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] worker=${worker_idx} complete"
}

launch_workers() {
  local worker_idx

  for worker_idx in 0 1 2 3; do
    local session="ce_rebuttal_theory_delta_w${worker_idx}"
    local log_file="${LOG_DIR}/${session}.log"

    if tmux has-session -t "${session}" 2>/dev/null; then
      continue
    fi

    tmux new-session -d -s "${session}" \
      "bash -lc '${SCRIPT_REMOTE} worker ${worker_idx}' > '${log_file}' 2>&1"
  done
}

wait_and_launch() {
  local session

  while true; do
    local any_running=0

    for session in "${WAIT_SESSIONS[@]}"; do
      if tmux has-session -t "${session}" 2>/dev/null; then
        any_running=1
      fi
    done

    if (( any_running == 0 )); then
      break
    fi

    sleep 60
  done

  launch_workers
}

case "${1:-queue}" in
  worker)
    run_worker "$2"
    ;;
  launch)
    launch_workers
    ;;
  queue)
    wait_and_launch
    ;;
  *)
    echo "Usage: $0 [queue|launch|worker INDEX]" >&2
    exit 2
    ;;
esac
