# Block bootstrap QQ calibration

The null block-bootstrap experiment in Section 4.1 of the paper. `scripts/run_boot_qq.R` sources the other R modules in `scripts/`. The reported grid is `n={200,400,800}`, `p={2000,4000,8000}`, `beta={2,0.9,0.55}`, `b=4`, 200 Monte Carlo replications, block length `floor(n^(2/3))`, and seed 123.

Run a cell directly with:

```bash
Rscript scripts/run_boot_qq.R --n 400 --p 4000 --beta 0.9 --b 4 --mc 200 --alpha 0.10 --seed 123 --ntrunc-mode n2
```

On Mercury, `bash shell/mercury_submit_sbatch_grid.sh` submits the paper grid. Adjust Slurm paths and account settings for other clusters. All generated outputs go to ignored `Results/`.
