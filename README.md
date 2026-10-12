# Time-series GGM experiments

Code for *High-dimensional Gaussian Graphical Model Testing for Long-Memory Time Series*, including the CE tuning checks discussed in the rebuttal. Generated data, fMRI inputs, figures, tables, and logs are excluded. All outputs are written beneath ignored `Results/` directories.

| Directory | Paper result | Run |
| --- | --- | --- |
| `Gaussian_sim/` | Null Gaussian approximation, QQ plots and distance tables | `bash shell/submit_sbatch_grid.sh` or a direct R cell |
| `boot_QQ/` | Null block-bootstrap QQ plots and calibration | `bash shell/submit_sbatch_grid.sh` or a direct R cell |
| `rejection_test/wide_band/` | Four off-band κ rejection tables | `bash shell/launch_kappa_paper_grids.sh` |
| `rejection_test/within_band/` | Within-band CE rejection table and rebuttal sensitivity | `bash shell/launch_paper_table.sh` and rebuttal launchers |
| `sherlock/` | fMRI graph comparison and long-memory diagnostics | Four R scripts in `scripts/` |

## Environment

- The Gaussian and bootstrap experiments use R; the bootstrap runner can use `matrixStats`. The reported cluster runs used R 4.4.3. Each directory's README gives a direct command and the full parameter grid. The Slurm wrappers require a Slurm installation. Set `SLURM_ACCOUNT`, `SLURM_QOS`, and `SLURM_PARTITION` only if your cluster requires them; `RSCRIPT_PATH` can select a specific R binary.
- The two rejection experiments require Python, NumPy, a CUDA-compatible CuPy installation, and a CUDA GPU. Recorded per-cell results used NumPy 2.4.4 and CuPy 14.0.1 with Python 3.11 or 3.12. Activate an appropriate environment before launching; the wide-band script also accepts `TGG_PYTHON=/path/to/python3`. The wide-band launcher requires `tmux`.
- Sherlock requires R packages `Matrix`, `data.table`, `future`, and `future.apply`, plus separately obtained, processed 268-ROI input data. See [`sherlock/README.md`](sherlock/README.md) for the exact file contract and data source.

The simulation code creates its own synthetic data. For each table, the directory README gives the report builder that converts the ignored per-cell outputs to the published rows. The Sherlock analysis is **not reproducible from this repository alone**: the source-to-268-ROI preprocessing procedure and processed inputs are not included. Full paper reproduction also requires the specified compute environments and full-size runs, which this release has not verified on a fresh machine.
