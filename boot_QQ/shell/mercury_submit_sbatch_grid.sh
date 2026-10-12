#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  mercury_submit_sbatch_grid.sh [options]

Submit the block-bootstrap QQ grid as separate Slurm jobs on Mercury.

Options:
  --cores C             CPUs per Slurm task and MC workers. Default: 16
  --fft-cores C         Coordinate FFT worker processes. Default: 8
  --blas-threads T      BLAS/OpenMP threads per worker. Default: 1
  --mem MEM             Slurm memory request for all jobs. Default: p-dependent
  --mem-p2000 MEM       Memory for p=2000 jobs. Default: 32G
  --mem-p4000 MEM       Memory for p=4000 jobs. Default: 48G
  --mem-p8000 MEM       Memory for p=8000 jobs. Default: 64G
  --time TIME           Slurm time limit. Default: 7-00:00:00
  --account ACCOUNT     Slurm account. Default: phd
  --qos QOS             Slurm QOS. Default: bronze
  --partition PARTITION Slurm partition. Default: standard
  --b B                 Graph bandwidth. Default: 4
  --n-mc N              Monte Carlo replications. Default: 200
  --alpha A             Bootstrap test level. Default: 0.10
  --seed SEED           RNG seed. Default: 123
  --ntrunc-mode MODE    Truncation mode. Default: n2
  --rscript PATH        Rscript path. Default: /apps/R-4.4.3/bin/Rscript
  --n-list VALUES       Space-delimited n list in quotes. Default: "200 400 800"
  --p-list VALUES       Space-delimited p list in quotes. Default: "2000 4000 8000"
  --beta-list VALUES    Space-delimited beta list in quotes. Default: "2 0.9 0.55"
  --no-overwrite        Do not pass overwrite to the R runner.
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cores="16"
fft_cores="8"
blas_threads="1"
mem=""
mem_p2000="32G"
mem_p4000="48G"
mem_p8000="64G"
time_limit="7-00:00:00"
account="phd"
qos="bronze"
partition="standard"
b="4"
n_mc="200"
alpha="0.10"
seed="123"
ntrunc_mode="n2"
rscript_path="/apps/R-4.4.3/bin/Rscript"
n_list="200 400 800"
p_list="2000 4000 8000"
beta_list="2 0.9 0.55"
overwrite="true"

while [[ $# -gt 0 ]]; do
  case "$1" in
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
    --mem)
      mem="$2"
      shift 2
      ;;
    --mem-p2000)
      mem_p2000="$2"
      shift 2
      ;;
    --mem-p4000)
      mem_p4000="$2"
      shift 2
      ;;
    --mem-p8000)
      mem_p8000="$2"
      shift 2
      ;;
    --time)
      time_limit="$2"
      shift 2
      ;;
    --account)
      account="$2"
      shift 2
      ;;
    --qos)
      qos="$2"
      shift 2
      ;;
    --partition)
      partition="$2"
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
    --ntrunc-mode)
      ntrunc_mode="$2"
      shift 2
      ;;
    --rscript)
      rscript_path="$2"
      shift 2
      ;;
    --n-list)
      n_list="$2"
      shift 2
      ;;
    --p-list)
      p_list="$2"
      shift 2
      ;;
    --beta-list)
      beta_list="$2"
      shift 2
      ;;
    --no-overwrite)
      overwrite="false"
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

timestamp="$(date +%Y%m%d_%H%M%S)"
log_dir="${PROJECT_ROOT}/Results/logs"
mkdir -p "${log_dir}"
manifest_csv="${log_dir}/sbatch_manifest_${timestamp}.csv"

echo "submit_ts,slurm_job_id,job_name,n,p,b,beta,n_mc,alpha,seed,cores,fft_cores,blas_threads,mem,time_limit,account,qos,partition,slurm_out" > "${manifest_csv}"

overwrite_arg="--overwrite"
if [[ "${overwrite}" != "true" ]]; then
  overwrite_arg="--no-overwrite"
fi

choose_mem() {
  local p_value="$1"
  if [[ -n "${mem}" ]]; then
    echo "${mem}"
  elif [[ "${p_value}" == "2000" ]]; then
    echo "${mem_p2000}"
  elif [[ "${p_value}" == "4000" ]]; then
    echo "${mem_p4000}"
  elif [[ "${p_value}" == "8000" ]]; then
    echo "${mem_p8000}"
  else
    echo "${mem_p8000}"
  fi
}

for n in ${n_list}; do
  for p in ${p_list}; do
    for beta in ${beta_list}; do
      beta_tag="${beta//./p}"
      job_name="bootqq_n${n}_p${p}_b${b}_be${beta_tag}"
      job_mem="$(choose_mem "${p}")"
      slurm_out="${log_dir}/slurm_${job_name}_%j.out"
      submit_output="$(
        sbatch \
          --job-name="${job_name}" \
          --cpus-per-task="${cores}" \
          --mem="${job_mem}" \
          --time="${time_limit}" \
          --account="${account}" \
          --qos="${qos}" \
          --partition="${partition}" \
          --chdir="${PROJECT_ROOT}" \
          --export="ALL,BOOT_QQ_ROOT=${PROJECT_ROOT},RSCRIPT_PATH=${rscript_path},BLAS_THREADS=${blas_threads}" \
          --output="${slurm_out}" \
          "${SCRIPT_DIR}/mercury_sbatch_job.sh" \
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
          ${overwrite_arg}
      )"
      slurm_job_id="$(awk '{print $4}' <<< "${submit_output}")"
      printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$(date -Is)" "${slurm_job_id}" "${job_name}" "${n}" "${p}" "${b}" "${beta}" \
        "${n_mc}" "${alpha}" "${seed}" "${cores}" "${fft_cores}" "${blas_threads}" "${job_mem}" "${time_limit}" \
        "${account}" "${qos}" "${partition}" "${slurm_out//%j/${slurm_job_id}}" >> "${manifest_csv}"
      echo "${submit_output} ${job_name}"
    done
  done
done

echo "manifest_csv=${manifest_csv}"
