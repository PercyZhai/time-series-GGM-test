#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  mercury_run_one.sh --n N --p P --beta BETA [options]

Options:
  --b B                 Graph bandwidth. Default: 4
  --n-mc N              Monte Carlo replications. Default: 200
  --b-gauss B           Gaussian target draws. Default: 2000
  --seed SEED           RNG seed. Default: 123
  --cores C             Monte Carlo/target worker processes. Default: 16
  --fft-cores C         Coordinate FFT worker processes. Default: 8
  --blas-threads T      BLAS/OpenMP threads per worker. Default: 1
  --n-trunc N           Optional truncation length override.
  --rscript PATH        Rscript path. Default: /apps/R-4.4.3/bin/Rscript
  --log PATH            Optional log path. Default: Results/logs/job_*.log
  --overwrite           Pass --overwrite to the R runner.
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

n=""
p=""
beta=""
b="4"
n_mc="200"
b_gauss="2000"
seed="123"
cores="16"
fft_cores="8"
blas_threads="1"
n_trunc=""
rscript_path="/apps/R-4.4.3/bin/Rscript"
run_log=""
overwrite="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --n)
      n="$2"
      shift 2
      ;;
    --p)
      p="$2"
      shift 2
      ;;
    --beta)
      beta="$2"
      shift 2
      ;;
    --b)
      b="$2"
      shift 2
      ;;
    --n-mc)
      n_mc="$2"
      shift 2
      ;;
    --b-gauss)
      b_gauss="$2"
      shift 2
      ;;
    --seed)
      seed="$2"
      shift 2
      ;;
    --threads|--cores)
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
    --n-trunc)
      n_trunc="$2"
      shift 2
      ;;
    --rscript)
      rscript_path="$2"
      shift 2
      ;;
    --log)
      run_log="$2"
      shift 2
      ;;
    --overwrite)
      overwrite="true"
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "${n}" || -z "${p}" || -z "${beta}" ]]; then
  echo "Missing required --n, --p, or --beta." >&2
  usage >&2
  exit 2
fi

if [[ ! -x "${rscript_path}" ]]; then
  echo "Rscript not executable: ${rscript_path}" >&2
  exit 127
fi

timestamp="$(date +%Y%m%d_%H%M%S)"
beta_tag="${beta//./p}"
mkdir -p "${PROJECT_ROOT}/Results/logs"

if [[ -z "${run_log}" ]]; then
  run_log="${PROJECT_ROOT}/Results/logs/job_n${n}_p${p}_b${b}_beta${beta_tag}_${timestamp}.log"
fi
mkdir -p "$(dirname "${run_log}")"

export OMP_NUM_THREADS="${blas_threads}"
export OPENBLAS_NUM_THREADS="${blas_threads}"
export MKL_NUM_THREADS="${blas_threads}"
export VECLIB_MAXIMUM_THREADS="${blas_threads}"
export NUMEXPR_NUM_THREADS="${blas_threads}"

cmd=(
  "${rscript_path}" "${PROJECT_ROOT}/scripts/run_gaussianity_sim.R"
  --preset custom
  --n-list "${n}"
  --p-list "${p}"
  --beta-list "${beta}"
  --b-list "${b}"
  --n-mc "${n_mc}"
  --b-gauss "${b_gauss}"
  --seed "${seed}"
  --cores "${cores}"
  --fft-cores "${fft_cores}"
)

if [[ -n "${n_trunc}" ]]; then
  cmd+=(--n-trunc "${n_trunc}")
fi

if [[ "${overwrite}" == "true" ]]; then
  cmd+=(--overwrite)
fi

finish() {
  local status=$?
  {
    echo "end_iso=$(date -Is)"
    echo "exit_status=${status}"
  } >> "${run_log}"
  exit "${status}"
}
trap finish EXIT

{
  echo "job_tag=n${n}_p${p}_b${b}_beta${beta}"
  echo "start_iso=$(date -Is)"
  echo "host=$(hostname)"
  echo "cwd=${PROJECT_ROOT}"
  echo "rscript=${rscript_path}"
  echo "cores=${cores}"
  echo "fft_cores=${fft_cores}"
  echo "blas_threads=${blas_threads}"
  echo "n=${n}"
  echo "p=${p}"
  echo "b=${b}"
  echo "beta=${beta}"
  echo "n_mc=${n_mc}"
  echo "b_gauss=${b_gauss}"
  echo "seed=${seed}"
  echo "n_trunc=${n_trunc:-default}"
  echo "overwrite=${overwrite}"
  free -h || true
  echo "command=${cmd[*]}"
  /usr/bin/time -p "${cmd[@]}"
  r_status=$?
  echo "rscript_status=${r_status}"
  exit "${r_status}"
} 2>&1 | tee -a "${run_log}"
