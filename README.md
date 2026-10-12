# TGG experiments — code release

Code for *High-dimensional Gaussian Graphical Model Testing for Long-Memory Time Series* and its rebuttal experiments. This branch contains source and launch scripts only. Datasets, generated samples, figures, tables, and logs remain in the authors' local copy and are not distributed here.

| Directory | Experiment | Main entry point |
| --- | --- | --- |
| `Gaussian_sim/` | Gaussian approximation under the null | `scripts/run_gaussianity_sim.R` |
| `boot_QQ/` | Block bootstrap QQ calibration | `scripts/run_boot_qq.R` |
| `rejection_test/wide_band/` | Off-band edge rejection (κ grids) | `shell/launch_kappa_paper_grids.sh` |
| `rejection_test/within_band/` | Consistency-empowered within-band weakening and rebuttal sensitivity | `scripts/run_withinband_weakening_ce_v2_gpu.py` |
| `sherlock/` | Sherlock movie/recall two-sample analysis and long-memory diagnostics | `scripts/run_omega_compute.R` |

## Requirements

- R for the null experiments; `boot_QQ` also uses `matrixStats`. The Sherlock pipeline uses `Matrix`, `data.table`, `future`, and `future.apply`.
- Python with NumPy and a CuPy build compatible with the installed CUDA runtime for the GPU rejection experiments. Launchers use `tmux`; the null-grid launchers use Slurm on a cluster.

Each script writes into an experiment-local `Results/` directory, which Git ignores. For Sherlock, place the independently obtained input files in `sherlock/processed_fMRI_data/` before running. See each experiment's README for commands and the parameter grids. The rebuttal CE launchers are `rejection_test/within_band/shell/launch_rebuttal_ce_sensitivity.sh` and `launch_rebuttal_ce_theory_delta.sh`; they use seed 303 and write under `Results/raw/`.

The code reflects the local working versions of the Gaussian and bootstrap runners, including the changes made after the original TGG archive commit. Historical run snapshots under `mercury_pull_2026-09-28/` and legacy helpers are outside the paper/rebuttal code release.
