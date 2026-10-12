# Within-band CE experiment

GPU/CuPy implementation for the four-regime within-band weakening table in Section 4.2 and the rebuttal's CE tuning checks. The runner imports simulation utilities from `../wide_band/scripts/run_rejection_test_gpu.py`.

From this directory, activate a CUDA-compatible Python environment, then run:

```bash
bash shell/launch_paper_table.sh
```

This runs the 16 published combinations of four `(n,p,beta,block_length)` regimes and `(eta,chi)={(0,0),(0.1,0.5),(0.3,1),(0.5,1)}`. It uses `theta0=rho=1`, `delta_n=0.05`, `C_delta=0.5`, 200 epochs, and seed 101. Results and the table are generated under ignored `Results/`. To rebuild the paper table from saved per-cell JSON, use `python3 scripts/build_paper_table.py`.

The rebuttal discusses sensitivity to `delta_n` and `C_delta`. `bash shell/launch_rebuttal_ce_sensitivity.sh` runs the nine parameter pairs on four regimes under H0 and one H1 alternative. `bash shell/launch_rebuttal_ce_theory_delta.sh launch` runs four theory-scaled threshold cells under each state; its default `queue` mode waits for sensitivity sessions to finish. These launchers use `tmux`, seed 303, and write under ignored `Results/raw/`. They will use `rejection_test/.venv_gpu` when present; otherwise Python 3 must be available on `PATH`.
