#!/usr/bin/env bash
set -euo pipefail

# v2 null delta_n-sweep: 4 workers, one delta_n_constant each in
# {0.50, 1.00, 1.50, 2.00}, theta0 = 1.0, rho = 1.0, eta = 0, chi = 0,
# all four regimes per worker.  Self-detaching: forks a controller tmux
# session and returns immediately.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXP="${ROOT}/within_band"
RESULTS="${EXP}/Results"
LOG_DIR="${RESULTS}/logs"
OUT_BASE="${RESULTS}/raw"
SCRIPT_REMOTE="${EXP}/shell/launch_wugroup_withinband_ce_v2_null_delta_sweep.sh"
PYTHON_SCRIPT="${EXP}/scripts/run_withinband_weakening_ce_v2_gpu.py"
QUEUE_SESSION="withinband_ce_v2_delta_sweep_queue"
QUEUE_LOG="${LOG_DIR}/withinband_ce_v2_delta_sweep_queue.log"

DELTAS=(0.2 0.5 0.8 1.2)
DELTA_TAGS=(0.2 0.5 0.8 1.2)
THETA0="1.0"
RHO="1.0"
REGIME_GRID="2000:800:2:160,2000:800:0.9:238,800:2000:2:86,800:2000:0.9:86"

mkdir -p "${LOG_DIR}" "${OUT_BASE}"

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
  echo "Started controller ${QUEUE_SESSION}; log=${QUEUE_LOG}; out_base=${OUT_BASE}"
  exit 0
}

start_workers() {
  log "Starting four within-band CE v2 null delta-sweep worker tmux sessions."
  log "delta_n_constants: ${DELTAS[*]}; theta0=${THETA0}; rho=${RHO}; regimes=${REGIME_GRID}; epochs=200; eta=0; chi=0."
  for idx in 0 1 2 3; do
    local delta="${DELTAS[$idx]}"
    local tag="${DELTA_TAGS[$idx]}"
    local session="withinband_ce_v2_dsweep_w${idx}"
    local log_file="${LOG_DIR}/${session}.log"
    local out_dir="${OUT_BASE}/delta${tag}"
    mkdir -p "${out_dir}"

    if tmux has-session -t "${session}" 2>/dev/null; then
      log "Worker session ${session} already exists; leaving it alone."
      continue
    fi

    tmux new-session -d -s "${session}" \
      "bash -lc 'cd ${EXP} && if [ -f ${ROOT}/.venv_gpu/bin/activate ]; then source ${ROOT}/.venv_gpu/bin/activate; fi && \
echo \"[\$(date +%Y-%m-%d_%H:%M:%S_%Z)] starting v2 delta=${delta} theta=${THETA0} rho=${RHO}, all 4 regimes\"; \
python ${PYTHON_SCRIPT} \
  --regime_grid ${REGIME_GRID} \
  --s0 10 --theta0 ${THETA0} --rho ${RHO} --delta0 0.10 \
  --epochs 200 --alpha 0.10 --seed 101 --batch_epochs 5 \
  --eta_grid 0 --chi_grid 0 \
  --C_delta 1.0 --edge_seed 202 --gpu_chunk_cols 4 \
  --combo_offset 0 --combo_stride 1 \
  --delta_n_constant ${delta} \
  --run_label null_delta${tag}_w${idx} \
  --out_dir ${out_dir}; \
echo \"[\$(date +%Y-%m-%d_%H:%M:%S_%Z)] finished v2 delta=${delta}\"' > ${log_file} 2>&1"
    log "Started ${session} (delta=${delta}); log=${log_file}; out_dir=${out_dir}"
  done
}

wait_for_workers() {
  log "Waiting for workers to finish."
  while true; do
    local alive=()
    for idx in 0 1 2 3; do
      local session="withinband_ce_v2_dsweep_w${idx}"
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

summarize_results() {
  log "Summarizing per-delta results."
  python3 - "${OUT_BASE}" <<'PYEOF' | tee -a "${QUEUE_LOG}"
import json, glob, sys
from pathlib import Path

base = Path(sys.argv[1])
delta_dirs = ["delta0.2", "delta0.5", "delta0.8", "delta1.2"]
regimes = [
    ("n=2000, p=800,  beta=2  ", 2000, 800,  "2"),
    ("n=2000, p=800,  beta=0.9", 2000, 800,  "0.9"),
    ("n=800,  p=2000, beta=2  ", 800,  2000, "2"),
    ("n=800,  p=2000, beta=0.9", 800,  2000, "0.9"),
]

def find_cell_json(d_dir, n, p, beta):
    matches = sorted(glob.glob(str(d_dir / f"withinband_ce_n{n}_p{p}_beta{beta}_s010_eta0_chi0_Cdelta*_null_delta*_w*.json")))
    return matches[-1] if matches else None

for d in delta_dirs:
    d_dir = base / d
    print(f"\n=== {d} ===")
    print(f"{'regime':28s} {'CE_rej':>8s} {'pumps':>8s} {'T_tilde':>10s} {'T_hat':>10s} {'q_hat':>10s}")
    for label, n, p, beta in regimes:
        path = find_cell_json(d_dir, n, p, beta)
        if path is None:
            print(f"{label:28s} {'MISSING':>8s}")
            continue
        with open(path) as f:
            obj = json.load(f)
        s = obj.get("summary", {})
        ce = s.get("ce_rejection_rate", float("nan"))
        pumps = s.get("mean_number_pumped_total", float("nan"))
        tt = s.get("mean_T_tilde", float("nan"))
        th = s.get("mean_T_hat", float("nan"))
        qh = s.get("mean_q_hat", float("nan"))
        print(f"{label:28s} {ce:8.3f} {pumps:8.1f} {tt:10.3f} {th:10.3f} {qh:10.3f}")
PYEOF
}

main() {
  start_controller_if_needed
  log "Controller started for within-band weakening CE v2 null delta-sweep."
  log "Output base: ${OUT_BASE}"
  start_workers
  wait_for_workers
  summarize_results
  log "Within-band weakening CE v2 null delta-sweep complete."
}

main "$@"
