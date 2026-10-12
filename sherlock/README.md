# Sherlock fMRI analysis

The movie-to-recall graph containment analysis and long-memory diagnostics in the paper. The data are not included. Obtain the Sherlock inputs independently and place the processed movie array and per-subject recall files under `processed_fMRI_data/`, matching the paths expected by the scripts. This directory is ignored by Git.

Run in order from this directory:

```bash
N_WORKERS=17 Rscript scripts/run_omega_compute.R
N_WORKERS=17 Rscript scripts/run_omega_bisect.R
N_WORKERS=17 Rscript scripts/run_omega_alpha_grid.R
N_WORKERS=8 Rscript scripts/run_long_memory_diagnostics.R
```

The first three stages share outputs in ignored `Results/`; the long-memory diagnostics can run separately. Required R packages: `Matrix`, `data.table`, `future`, `future.apply`.
