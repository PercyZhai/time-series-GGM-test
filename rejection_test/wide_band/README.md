# Off-band rejection experiment

GPU/CuPy implementation for the four κ tables in Section 4.2 of the paper. `scripts/run_rejection_test_gpu.py` contains the DGP, bootstrap, and rejection-rate calculation. `shell/launch_kappa_paper_grids.sh` runs the reported `(n,p,beta,m,κ)` grid and preserves the nested edge-list seeds used in the paper.

From this directory, run:

```bash
bash shell/launch_kappa_paper_grids.sh
```

The script uses `tmux` and writes data, summaries, and logs under ignored `Results/`. Set `WORKERS` to change concurrency and `TGG_WIDEBAND_ENV_CMD` to activate a Python environment. Requires NumPy and CUDA-compatible CuPy.
