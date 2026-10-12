# Within-band consistency-empowered experiment

GPU/CuPy implementation for the within-band weakening experiment in Section 4.2 and the rebuttal's CE tuning sensitivity. The runner `scripts/run_withinband_weakening_ce_v2_gpu.py` imports DGP/bootstrap utilities from `../wide_band/scripts/run_rejection_test_gpu.py`. Requires NumPy and CUDA-compatible CuPy.

The paper grid uses `n,p,beta,block_length` values `(2000,800,2,160)`, `(2000,800,0.9,238)`, `(800,2000,2,86)`, `(800,2000,0.9,86)` with `theta0=rho=1`, `C_delta=0.5`, and constant `delta_n=0.05`. Run `bash shell/launch_wugroup_withinband_ce_v2_grid.sh` in a `tmux`-capable environment. The null calibration launcher is `shell/launch_wugroup_withinband_ce_v2_null_delta_sweep.sh`.

For the rebuttal's `(delta_n,C_delta)` sensitivity grid, run `bash shell/launch_rebuttal_ce_sensitivity.sh`; for the theory-scaled threshold cells, run `bash shell/launch_rebuttal_ce_theory_delta.sh`. Both use seed 303. The second script's default `queue` mode waits for sensitivity worker sessions to finish; use `launch` to start directly.

Launchers locate the repository from their own paths and optionally source `rejection_test/.venv_gpu/bin/activate`. Otherwise they use the current Python environment. All outputs go to ignored `Results/`; report builders in `scripts/` read those local outputs.
