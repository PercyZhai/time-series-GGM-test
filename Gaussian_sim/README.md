# Gaussian approximation

Section 4.1 null Gaussian approximation. The paper grid is `n={200,400,800}`, `p={2000,4000,8000}`, `beta={2,0.9,0.55}`, bandwidth `b=4`, 200 Monte Carlo replications, 2000 Gaussian draws, and seed 123.

Run one cell directly from this directory:

```bash
Rscript scripts/run_gaussianity_sim.R --preset custom --n-list 400 --p-list 4000 --beta-list 0.9 --b-list 4 --n-mc 200 --b-gauss 2000 --seed 123
```

On a Slurm cluster, `bash shell/submit_sbatch_grid.sh` submits the complete grid. Account, QOS, and partition can be supplied through the matching command options or `SLURM_ACCOUNT`, `SLURM_QOS`, and `SLURM_PARTITION`. Set `RSCRIPT_PATH` if `Rscript` is not on `PATH`. Generated per-cell QQ PDFs and summary CSV files are in ignored `Results/`. After all jobs finish, run `Rscript scripts/build_summary.R` to combine the per-cell summaries into one table; this avoids concurrent writes from separate Slurm jobs.
