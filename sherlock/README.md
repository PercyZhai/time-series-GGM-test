# Sherlock fMRI analysis

Movie-versus-recall graph comparison and movie-phase long-memory diagnostics. Input data are excluded. The [Sherlock data author's page](https://jchenlab.johnshopkins.edu/) links to the public raw BIDS and preprocessed NIFTI datasets. Those downloads are **not** the exact 268-ROI files expected by these scripts; the ROI extraction, subject matching, and any filtering used to create this analysis input have not been supplied in this code release.

Place these 18 processed inputs in ignored `processed_fMRI_data/`:

- `moviefmrip268n17T1976.RData`: an R array with dimensions 17 subjects × 268 ROIs × 1976 movie time points, in any axis order.
- Seventeen headered text matrices matching `recallfmrip268T*subject*.txt`, one per subject, with 268 ROI columns (or rows); the scripts orient them as time × ROI.

With R packages `Matrix`, `data.table`, `future`, and `future.apply` installed, run from this directory:

```bash
N_WORKERS=17 Rscript scripts/run_omega_compute.R
N_WORKERS=17 Rscript scripts/run_omega_bisect.R
N_WORKERS=17 Rscript scripts/run_omega_alpha_grid.R
N_WORKERS=8 Rscript scripts/run_long_memory_diagnostics.R
```

All generated files are under ignored `Results/`. The first three stages are sequential; long-memory diagnostics can run independently. Reproducing the paper's Sherlock result requires the same processed ROI inputs or their exact preprocessing recipe.
