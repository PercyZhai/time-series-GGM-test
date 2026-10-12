#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  submit_sbatch_grid.sh [options]

Submit the Gaussian approximation grid as separate Slurm jobs.

Options:
  --threads T           CPUs per Slurm task. Default: 16
  --mem MEM             Slurm memory request. Default: 56G
  --time TIME           Slurm time limit. Default: 7-00:00:00
  --account ACCOUNT     Optional Slurm account.
  --qos QOS             Optional Slurm QOS.
  --partition PARTITION Optional Slurm partition.
  --b B                 Graph bandwidth. Default: 4
  --n-mc N              Monte Carlo replications. Default: 200
  --b-gauss B           Gaussian target draws. Default: 2000
  --seed SEED           RNG seed. Default: 123
  --rscript PATH        Rscript command or path. Default: Rscript
  --n-list VALUES       Space-delimited n list in quotes. Default: "200 400 800"
  --p-list VALUES       Space-delimited p list in quotes. Default: "2000 4000 8000"
  --beta-list VALUES    Space-delimited beta list in quotes. Default: "2 0.9 0.55"
  --no-overwrite        Do not pass overwrite to the R runner.
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

threads="16"
mem="56G"
time_limit="7-00:00:00"
account="${SLURM_ACCOUNT:-}"
qos="${SLURM_QOS:-}"
partition="${SLURM_PARTITION:-}"
b="4"
n_mc="200"
b_gauss="2000"
seed="123"
rscript_path="${RSCRIPT_PATH:-Rscript}"
n_list="200 400 800"
p_list="2000 4000 8000"
beta_list="2 0.9 0.55"
overwrite="true"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads)
      threads="$2"
      shift 2
      ;;
    --mem)
      mem="$2"
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

echo "submit_ts,slurm_job_id,job_name,n,p,b,beta,n_mc,b_gauss,seed,threads,mem,time_limit,account,qos,partition,slurm_out" > "${manifest_csv}"

overwrite_arg="--overwrite"
if [[ "${overwrite}" != "true" ]]; then
  overwrite_arg="--no-overwrite"
fi

slurm_site_opts=("--chdir=${PROJECT_ROOT}")
[[ -n "${account}" ]] && slurm_site_opts+=("--account=${account}")
[[ -n "${qos}" ]] && slurm_site_opts+=("--qos=${qos}")
[[ -n "${partition}" ]] && slurm_site_opts+=("--partition=${partition}")

for n in ${n_list}; do
  for p in ${p_list}; do
    for beta in ${beta_list}; do
      beta_tag="${beta//./p}"
      job_name="gauss_n${n}_p${p}_b${b}_be${beta_tag}"
      slurm_out="${log_dir}/slurm_${job_name}_%j.out"
      submit_output="$(
        sbatch \
          --job-name="${job_name}" \
          --cpus-per-task="${threads}" \
          --mem="${mem}" \
          --time="${time_limit}" \
          "${slurm_site_opts[@]}" \
          --export="ALL,GAUSSIAN_SIM_ROOT=${PROJECT_ROOT},RSCRIPT_PATH=${rscript_path}" \
          --output="${slurm_out}" \
          "${SCRIPT_DIR}/sbatch_job.sh" \
          --n "${n}" \
          --p "${p}" \
          --beta "${beta}" \
          --b "${b}" \
          --n-mc "${n_mc}" \
          --b-gauss "${b_gauss}" \
          --seed "${seed}" \
          --rscript "${rscript_path}" \
          ${overwrite_arg}
      )"
      slurm_job_id="$(awk '{print $4}' <<< "${submit_output}")"
      printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${slurm_job_id}" "${job_name}" "${n}" "${p}" "${b}" "${beta}" \
        "${n_mc}" "${b_gauss}" "${seed}" "${threads}" "${mem}" "${time_limit}" \
        "${account}" "${qos}" "${partition}" "${slurm_out//%j/${slurm_job_id}}" >> "${manifest_csv}"
      echo "${submit_output} ${job_name}"
    done
  done
done

echo "manifest_csv=${manifest_csv}"
