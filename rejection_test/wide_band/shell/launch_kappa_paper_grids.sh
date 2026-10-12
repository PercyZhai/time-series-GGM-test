#!/usr/bin/env bash
set -euo pipefail

# Canonical launcher for the kappa grids reported in the paper tables.
# It intentionally includes only the m values shown in the paper tables.

ROOT="${TGG_WIDEBAND_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RESULTS="${ROOT}/Results/raw"
LOG_DIR="${ROOT}/Results/logs"
EDGE_DIR="${LOG_DIR}/paper_edge_files"
JOB_CSV="${LOG_DIR}/kappa_paper_grids_jobs.csv"
QUEUE_SESSION="${TGG_QUEUE_SESSION:-tgg_kappa_paper_grids_queue}"
QUEUE_LOG="${LOG_DIR}/${QUEUE_SESSION}.log"
WORKERS="${WORKERS:-4}"
SCRIPT_PATH="${ROOT}/shell/launch_kappa_paper_grids.sh"
TGG_PYTHON="$(command -v "${TGG_PYTHON:-python3}")"

mkdir -p "${RESULTS}" "${LOG_DIR}" "${EDGE_DIR}"

timestamp() {
  date "+%Y-%m-%d %H:%M:%S %Z"
}

log() {
  echo "[$(timestamp)] $*" | tee -a "${QUEUE_LOG}"
}

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
    "bash -lc 'IN_TMUX_CONTROLLER=1 TGG_WIDEBAND_ROOT=${ROOT@Q} TGG_PYTHON=${TGG_PYTHON@Q} WORKERS=${WORKERS@Q} ${SCRIPT_PATH@Q}' > ${QUEUE_LOG@Q} 2>&1"
  echo "Started controller ${QUEUE_SESSION}; log=${QUEUE_LOG}"
  exit 0
}

write_edge_files() {
  log "Preparing nested edge-prefix files outside Results/raw."
  cd "${ROOT}/scripts"
  "${TGG_PYTHON}" - <<'PY'
from pathlib import Path
import run_rejection_test_gpu as g

edge_dir = Path("../Results/logs/paper_edge_files")
edge_dir.mkdir(parents=True, exist_ok=True)

specs = [
    # These m_max values reproduce the historical nested-edge masters used
    # for the reported cells. Do not collapse them: the sampler's prefix can
    # depend on m_max through its chunking.
    (800, 4, 2000, [80, 200, 400]),
    (800, 4, 1600, [800, 1200]),
    (2000, 4, 5000, [200, 500, 1000]),
    (2000, 4, 4000, [2000, 3000]),
    (2000, 4, 3000, [200, 500, 1000, 2000, 3000]),
]

for p, b, m_max, prefixes in specs:
    edges = g.sample_ordered_nested_edges(
        p=p,
        b=b,
        m_max=m_max,
        edge_seed=202,
        sign_seed=303,
        edge_sampling="uniform_nested",
        degree_cap=None,
    )
    for m in prefixes:
        out = edge_dir / f"nested_edges_p{p}_b{b}_m{m}_from_m{m_max}_edge202_sign303.csv"
        g.write_ordered_edges_csv(out, edges[:m])
        print(f"wrote {out}")
PY
}

append_jobs_for_grid() {
  local n="$1"
  local p="$2"
  local beta="$3"
  local block_len="$4"
  local out_dir="$5"
  local m_values="$6"
  local edge_m_max="$7"
  local label_prefix="$8"

  IFS=',' read -r -a m_array <<< "${m_values}"
  for m_extra in "${m_array[@]}"; do
    local edge_file="${EDGE_DIR}/nested_edges_p${p}_b4_m${m_extra}_from_m${edge_m_max}_edge202_sign303.csv"
    for kappa in 0.03 0.06 0.09 0.12 0.15 0.18; do
      JOB_ID=$((JOB_ID + 1))
      echo "${JOB_ID},${n},${p},${beta},${block_len},${m_extra},${kappa},${edge_file},${out_dir},${label_prefix}_m${m_extra}_k${kappa}" >> "${JOB_CSV}"
    done
  done
}

write_job_csv() {
  log "Writing canonical paper-grid job CSV."
  echo "job_id,n,p,beta,block_len,m_extra,kappa,edge_file,out_dir,label" > "${JOB_CSV}"
  JOB_ID=0

  append_jobs_for_grid 2000 800 0.9 238 \
    "${RESULTS}/table1_n2000_p800_beta09_l238_m80_200_400_800_1200" \
    "80,200,400" 2000 "n2000_p800_beta09_l238"
  append_jobs_for_grid 2000 800 0.9 238 \
    "${RESULTS}/table1_n2000_p800_beta09_l238_m80_200_400_800_1200" \
    "800,1200" 1600 "n2000_p800_beta09_l238"

  append_jobs_for_grid 2000 800 2 160 \
    "${RESULTS}/table2_n2000_p800_beta2_l160_m80_200_400_800_1200" \
    "80,200,400" 2000 "n2000_p800_beta2_l160"
  append_jobs_for_grid 2000 800 2 160 \
    "${RESULTS}/table2_n2000_p800_beta2_l160_m80_200_400_800_1200" \
    "800,1200" 1600 "n2000_p800_beta2_l160"

  append_jobs_for_grid 800 2000 2 86 \
    "${RESULTS}/table3_n800_p2000_beta2_l86_m200_500_1000_2000_3000" \
    "200,500,1000" 5000 "n800_p2000_beta2_l86"
  append_jobs_for_grid 800 2000 2 86 \
    "${RESULTS}/table3_n800_p2000_beta2_l86_m200_500_1000_2000_3000" \
    "2000,3000" 4000 "n800_p2000_beta2_l86"

  append_jobs_for_grid 800 2000 0.9 86 \
    "${RESULTS}/table4_n800_p2000_beta09_l86_m200_500_1000_2000_3000" \
    "200,500,1000,2000,3000" 3000 "n800_p2000_beta09_l86"

  log "Wrote ${JOB_CSV} with ${JOB_ID} jobs."
}

start_workers() {
  log "Starting ${WORKERS} worker tmux sessions."
  for ((off = 0; off < WORKERS; off++)); do
    local session="tgg_kappa_paper_w${off}"
    local log_file="${LOG_DIR}/${session}.log"
    if tmux has-session -t "${session}" 2>/dev/null; then
      log "Worker session ${session} already exists; leaving it alone."
      continue
    fi
    tmux new-session -d -s "${session}" \
      "bash -lc 'IN_WORKER=1 TGG_WIDEBAND_ROOT=${ROOT@Q} TGG_PYTHON=${TGG_PYTHON@Q} WORKERS=${WORKERS@Q} WORKER_OFFSET=${off@Q} ${SCRIPT_PATH@Q} worker' > ${log_file@Q} 2>&1"
    log "Started ${session}; log=${log_file}"
  done
}

worker_main() {
  cd "${ROOT}"
  local off="${WORKER_OFFSET:?WORKER_OFFSET is required}"
  local workers="${WORKERS:?WORKERS is required}"

  tail -n +2 "${JOB_CSV}" | while IFS=, read -r job_id n p beta block_len m_extra kappa edge_file out_dir label; do
    local mod=$(( (job_id - 1) % workers ))
    if [ "${mod}" -ne "${off}" ]; then
      continue
    fi
    mkdir -p "${out_dir}"
    echo "[$(timestamp)] worker=${off} job=${job_id} label=${label}"
    "${TGG_PYTHON}" scripts/run_rejection_test_gpu.py \
      --n "${n}" --p "${p}" --beta "${beta}" --s0 10 --epochs 200 --alpha 0.10 \
      --seed 101 --batch_epochs 5 --block_length "${block_len}" --gpu_chunk_cols 4 \
      --diag_mode strict --signal_mode kappa \
      --m_extra_grid "${m_extra}" \
      --kappa_grid "${kappa}" \
      --edge_seed 202 --sign_seed 303 \
      --extra_edges_csv "${edge_file}" \
      --run_label "paper_w${off}_job${job_id}_${label}" \
      --out_dir "${out_dir}"
  done
}

wait_for_workers() {
  log "Waiting for workers to finish."
  while true; do
    local alive=()
    for ((off = 0; off < WORKERS; off++)); do
      local session="tgg_kappa_paper_w${off}"
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

controller_main() {
  start_controller_if_needed
  log "Canonical kappa paper-grid queue starting."
  log "Workers=${WORKERS}; root=${ROOT}"
  write_edge_files
  write_job_csv
  start_workers
  wait_for_workers
  log "Canonical kappa paper-grid queue complete."
}

if [ "${1:-}" = "worker" ]; then
  worker_main
else
  controller_main
fi
