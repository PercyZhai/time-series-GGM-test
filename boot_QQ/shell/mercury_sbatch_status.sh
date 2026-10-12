#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOG_DIR="${PROJECT_ROOT}/Results/logs"

manifest="${1:-}"
if [[ -z "${manifest}" ]]; then
  manifest="$(ls -1t "${LOG_DIR}"/sbatch_manifest_*.csv 2>/dev/null | head -n 1 || true)"
fi

if [[ -z "${manifest}" || ! -f "${manifest}" ]]; then
  echo "No sbatch manifest found." >&2
  exit 1
fi

job_ids="$(awk -F, 'NR > 1 { print $2 }' "${manifest}" | paste -sd, -)"

echo "manifest=${manifest}"
echo "job_count=$(awk 'END { print NR - 1 }' "${manifest}")"
echo

echo "squeue:"
if [[ -n "${job_ids}" ]]; then
  squeue -j "${job_ids}" -o "%.18i %.9T %.12M %.12l %.8C %.10m %.25j %.20R" || true
else
  echo "  no job ids"
fi
echo

echo "manifest tail:"
tail -n 12 "${manifest}"
echo

latest_job_log="$(ls -1t "${LOG_DIR}"/job_*.log 2>/dev/null | head -n 1 || true)"
if [[ -n "${latest_job_log}" ]]; then
  echo "latest_job_log=${latest_job_log}"
  grep -E '^(job_tag|start_iso|end_iso|real |user |sys |rscript_status|exit_status=)' "${latest_job_log}" || true
  grep -E 'Completed [0-9]+/[0-9]+' "${latest_job_log}" | tail -n 3 || true
  echo
  tail -n 20 "${latest_job_log}"
fi
