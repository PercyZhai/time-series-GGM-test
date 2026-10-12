#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00

set -euo pipefail

if [[ -n "${BOOT_QQ_ROOT:-}" ]]; then
  PROJECT_ROOT="${BOOT_QQ_ROOT}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi
SHELL_DIR="${PROJECT_ROOT}/shell"

cores="${SLURM_CPUS_PER_TASK:-16}"
fft_cores="8"
blas_threads="${BLAS_THREADS:-1}"
rscript_path="${RSCRIPT_PATH:-Rscript}"

b="4"
n_mc="200"
alpha="0.10"
seed="123"
ntrunc_mode="n2"
overwrite_arg="--overwrite"

args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --n|--p|--beta)
      args+=("$1" "$2")
      shift 2
      ;;
    --b)
      b="$2"
      shift 2
      ;;
    --n-mc|--mc)
      n_mc="$2"
      shift 2
      ;;
    --alpha)
      alpha="$2"
      shift 2
      ;;
    --seed)
      seed="$2"
      shift 2
      ;;
    --cores)
      cores="$2"
      shift 2
      ;;
    --fft-cores)
      fft_cores="$2"
      shift 2
      ;;
    --blas-threads)
      blas_threads="$2"
      shift 2
      ;;
    --ntrunc-mode)
      ntrunc_mode="$2"
      shift 2
      ;;
    --rscript)
      rscript_path="$2"
      shift 2
      ;;
    --no-overwrite)
      overwrite_arg=""
      shift
      ;;
    --overwrite)
      overwrite_arg="--overwrite"
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

n=""
p=""
beta=""
for ((idx = 0; idx < ${#args[@]}; idx += 2)); do
  case "${args[idx]}" in
    --n) n="${args[idx + 1]}" ;;
    --p) p="${args[idx + 1]}" ;;
    --beta) beta="${args[idx + 1]}" ;;
  esac
done

if [[ -z "${n}" || -z "${p}" || -z "${beta}" ]]; then
  echo "Missing required --n, --p, or --beta." >&2
  exit 2
fi

beta_tag="${beta//./p}"
mkdir -p "${PROJECT_ROOT}/Results/logs"
job_log="${PROJECT_ROOT}/Results/logs/job_n${n}_p${p}_b${b}_beta${beta_tag}_slurm${SLURM_JOB_ID:-manual}.log"

echo "slurm_job_id=${SLURM_JOB_ID:-manual}"
echo "slurm_job_name=${SLURM_JOB_NAME:-manual}"
echo "slurm_nodelist=${SLURM_NODELIST:-unknown}"
echo "slurm_cpus_per_task=${cores}"
echo "job_log=${job_log}"

bash "${SHELL_DIR}/run_one.sh" \
  --n "${n}" \
  --p "${p}" \
  --beta "${beta}" \
  --b "${b}" \
  --n-mc "${n_mc}" \
  --alpha "${alpha}" \
  --seed "${seed}" \
  --cores "${cores}" \
  --fft-cores "${fft_cores}" \
  --blas-threads "${blas_threads}" \
  --ntrunc-mode "${ntrunc_mode}" \
  --rscript "${rscript_path}" \
  --log "${job_log}" \
  ${overwrite_arg}
