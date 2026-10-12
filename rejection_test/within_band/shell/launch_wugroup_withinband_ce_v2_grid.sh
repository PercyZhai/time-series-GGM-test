#!/usr/bin/env bash
set -euo pipefail

# v2 full eta x chi grid (5 x 5 = 25), all four regimes, 4 workers.
# DGP fix vs v1: --rho 1.0 (was 0.5).
# Calibration set after the v2 null delta-sweep:
#   --theta0 1.0
#   --C_delta 0.5
#   --delta_n_constant 0.05
# These are the only (delta, C_delta) combo that kept CE rejection rate
# in [0.05, 0.15] under H0 across all four regimes (see V2_NULL_DELTA_SWEEP_REPORT.md).
#

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXP="${ROOT}/within_band"
RESULTS="${EXP}/Results"
LOG_DIR="${RESULTS}/logs"
OUT_DIR="${RESULTS}/raw"
SCRIPT_REMOTE="${EXP}/shell/launch_wugroup_withinband_ce_v2_grid.sh"
PYTHON_SCRIPT="${EXP}/scripts/run_withinband_weakening_ce_v2_gpu.py"
QUEUE_SESSION="withinband_ce_v2_full_wugroup_queue"
QUEUE_LOG="${LOG_DIR}/${QUEUE_SESSION}.log"

mkdir -p "${LOG_DIR}" "${OUT_DIR}"

timestamp() { date "+%Y-%m-%d %H:%M:%S %Z"; }
log() { echo "[$(timestamp)] $*" | tee -a "${QUEUE_LOG}"; }

start_controller_if_needed() {
  if [ "${IN_TMUX_CONTROLLER:-0}" = "1" ]; then
    return
  fi
  if tmux has-session -t "${QUEUE_SESSION}" 2>/dev/null; then
    echo "Controller session ${QUEUE_SESSION} already exists; leaving it alone."
    echo "Log: ${QUEUE_LOG}"
    exit 0
  fi
  tmux new-session -d -s "${QUEUE_SESSION}" \
    "bash -lc 'IN_TMUX_CONTROLLER=1 ${SCRIPT_REMOTE}' > ${QUEUE_LOG} 2>&1"
  echo "Started controller ${QUEUE_SESSION}; log=${QUEUE_LOG}; out_dir=${OUT_DIR}"
  exit 0
}

start_workers() {
  log "Starting four within-band CE v2 full-grid worker tmux sessions."
  for off in 0 1 2 3; do
    local session="withinband_ce_v2_full_w${off}"
    local log_file="${LOG_DIR}/${session}.log"
    if tmux has-session -t "${session}" 2>/dev/null; then
      log "Worker session ${session} already exists; leaving it alone."
      continue
    fi

    tmux new-session -d -s "${session}" \
      "bash -lc 'cd ${EXP} && if [ -f ${ROOT}/.venv_gpu/bin/activate ]; then source ${ROOT}/.venv_gpu/bin/activate; fi && \
python ${PYTHON_SCRIPT} \
  --regime_grid 2000:800:2:160,2000:800:0.9:238,800:2000:2:86,800:2000:0.9:86 \
  --s0 10 --theta0 1.0 --rho 1.0 --delta0 0.10 \
  --epochs 200 --alpha 0.10 --seed 101 --batch_epochs 5 \
  --eta_grid 0.02,0.05,0.10,0.20,0.30 \
  --chi_grid 0.20,0.40,0.60,0.80,1.00 \
  --C_delta 0.5 --edge_seed 202 --gpu_chunk_cols 4 \
  --delta_n_constant 0.05 \
  --combo_offset ${off} --combo_stride 4 \
  --run_label v2_full_w${off} \
  --out_dir ${OUT_DIR}' > ${log_file} 2>&1"
    log "Started ${session}; log=${log_file}"
  done
}

wait_for_workers() {
  log "Waiting for workers to finish."
  while true; do
    local alive=()
    for session in withinband_ce_v2_full_w0 withinband_ce_v2_full_w1 withinband_ce_v2_full_w2 withinband_ce_v2_full_w3; do
      if tmux has-session -t "${session}" 2>/dev/null; then
        alive+=("${session}")
      fi
    done
    if [ "${#alive[@]}" -eq 0 ]; then
      log "All workers finished."
      break
    fi
    log "Still running: ${alive[*]}"
    sleep 300
  done
}

main() {
  start_controller_if_needed
  log "Controller started for within-band weakening CE v2 full grid."
  log "Output directory: ${OUT_DIR}"
  start_workers
  wait_for_workers
  log "Within-band weakening CE v2 full grid complete."
}

main "$@"
