# Off-band rejection experiment

GPU/CuPy implementation for the four κ tables in Section 4.2. From this directory, with a CUDA-compatible Python environment activated, run:

```bash
bash shell/launch_kappa_paper_grids.sh
python3 scripts/build_paper_tables.py
```

The launcher uses `tmux`; it starts a controller and returns immediately. Follow progress in `Results/logs/`. `WORKERS` sets the number of worker sessions (default 4). Set `TGG_PYTHON` to a Python executable if needed. The launcher generates 120 cells with 200 epochs each, uses the paper's nested edge-list seeds, and writes per-cell files to ignored `Results/raw/`. Run the table builder after the workers finish; it checks that all 120 cells and epoch CSV files exist before writing `Results/tables/kappa_paper_tables.md`.
