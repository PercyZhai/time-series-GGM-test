# Gaussian approximation

The direct-statistic Gaussian approximation experiment in Section 4.1 of the paper. `scripts/run_gaussianity_sim.R` sources the other R modules in `scripts/`. The reported grid is `n={200,400,800}`, `p={2000,4000,8000}`, `beta={2,0.9,0.55}`, bandwidth `b=4`, 200 Monte Carlo replications, 2000 Gaussian draws, and seed 123.

Run a cell directly with:

```bash
Rscript scripts/run_gaussianity_sim.R --preset custom --n-list 400 --p-list 4000 --beta-list 0.9 --b-list 4 --n-mc 200 --b-gauss 2000 --seed 123
```

On Mercury, `bash shell/mercury_submit_sbatch_grid.sh` submits the paper grid. The Slurm scripts assume a cluster R path and account settings; adjust these to the target environment. All generated outputs go to ignored `Results/`.
