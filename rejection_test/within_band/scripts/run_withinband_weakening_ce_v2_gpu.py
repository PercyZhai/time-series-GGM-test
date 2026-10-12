#!/usr/bin/env python3

"""GPU runner for the within-band weakening CE experiment.

This experiment intentionally lives in its own folder.  It keeps the accepted
FFT DGP, GPU data generation, direct TGG bootstrap critical value, and block
lengths from the latest kappa runs, but changes the alternative to weaken a
nested fraction of null within-band precision edges.

The rejection statistic is the consistency-empowered statistic from Draft 0505:

    T_tilde = sqrt(n) max_{i,j} | e_i' Sigma_hat omega_tilde_0j - 1{i=j} |.

The direct statistic T_hat is computed only as a diagnostic.  Rejection uses
T_tilde > q_hat, where q_hat is still the direct-statistic bootstrap critical
value.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
import time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any

import cupy as cp
import numpy as np

EXPERIMENT_DIR = Path(__file__).resolve().parents[1]
REJECTION_TEST_DIR = EXPERIMENT_DIR.parent
WIDE_BAND_SCRIPTS_DIR = REJECTION_TEST_DIR / "wide_band" / "scripts"

# Reuse the accepted GPU/DGP/bootstrap utilities from the retained wide-band
# experiment.  Keeping this import explicit avoids a hidden current-working-dir
# dependency when the within-band runner is launched from tmux or SSH.
for import_dir in (WIDE_BAND_SCRIPTS_DIR, EXPERIMENT_DIR / "scripts"):
    if str(import_dir) not in sys.path:
        sys.path.insert(0, str(import_dir))

from run_rejection_test_gpu import (  # noqa: E402
    RunParams,
    empirical_quantile_type1,
    format_seconds,
    generate_observed_batch,
    prepare_fft_filter,
    rolling_block_max_abs,
    save_json,
    write_csv,
)


THETA0_DEFAULT = 0.25
RHO_DEFAULT = 0.5
DELTA0_DEFAULT = 0.10


def parse_float_grid(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def parse_regime_grid(text: str) -> list[tuple[int, int, float, int]]:
    out: list[tuple[int, int, float, int]] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        pieces = item.split(":")
        if len(pieces) != 4:
            raise ValueError(
                "--regime_grid entries must be n:p:beta:block_len, e.g. 2000:800:2:160"
            )
        n, p, beta, block_len = pieces
        out.append((int(n), int(p), float(beta), int(block_len)))
    return out


def regime_label(n: int, p: int, beta: float) -> str:
    beta_label = "beta09" if abs(beta - 0.9) < 1e-12 else f"beta{beta:g}"
    return f"n{n}_p{p}_{beta_label}"


def make_params(
    args: argparse.Namespace,
    n: int,
    p: int,
    beta: float,
    block_length: int,
) -> RunParams:
    b0 = (args.s0 - 1) // 2
    return RunParams(
        n=n,
        p=p,
        beta=beta,
        s0=args.s0,
        b=b0,
        m_extra=0,
        theta_extra=float("nan"),
        h=float("nan"),
        epochs=args.epochs,
        alpha=args.alpha,
        seed=args.seed,
        batch_epochs=args.batch_epochs,
        block_length=block_length,
        filter_length=args.fft_len if args.fft_len is not None else n * n,
        theta0=args.theta0,
        rho=args.rho,
        delta0=args.delta0,
        out_dir=args.out_dir,
        gpu_chunk_cols=args.gpu_chunk_cols,
        diag_mode="strict_diagonal_fixed_after_weakening",
        signal_mode="withinband_weakening_ce",
        kappa=float("nan"),
        edge_seed=args.edge_seed,
        sign_seed=0,
        edge_sampling="within_band_nested",
        degree_cap=None,
        kappa_tol=float("nan"),
        kappa_max_iter=0,
        kappa_damping=float("nan"),
    )


def baseline_precision_components(params: RunParams) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    omega = np.zeros((params.p, params.p), dtype=np.float64)
    row_abs = np.zeros(params.p, dtype=np.float64)
    for dist in range(1, params.b + 1):
        val = params.theta0 * (params.rho ** (dist - 1))
        idx = np.arange(params.p - dist)
        omega[idx, idx + dist] = val
        omega[idx + dist, idx] = val
        row_abs[idx] += abs(val)
        row_abs[idx + dist] += abs(val)
    diag = 1.0 + row_abs + params.delta0
    np.fill_diagonal(omega, diag)
    return omega, diag, row_abs


def ordered_withinband_edges(p: int, b0: int, seed: int) -> np.ndarray:
    rows = [(i, i + dist) for dist in range(1, b0 + 1) for i in range(p - dist)]
    edges = np.asarray(rows, dtype=np.int64)
    rng = np.random.default_rng(seed)
    order = rng.permutation(edges.shape[0])
    return edges[order, :]


def write_edges_csv(path: Path, edges: np.ndarray) -> None:
    rows = [{"rank": idx + 1, "i": int(i + 1), "j": int(j + 1)} for idx, (i, j) in enumerate(edges)]
    write_csv(path, rows, ["rank", "i", "j"])


def build_withinband_precision(
    params: RunParams,
    eta: float,
    chi: float,
    ordered_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    omega0, diag0, row_abs0 = baseline_precision_components(params)
    total_edges = ordered_edges.shape[0]
    m_weakened = int(math.floor(eta * total_edges + 1e-12))
    selected = ordered_edges[:m_weakened, :].copy()
    omega1 = omega0.copy()
    if selected.size:
        i = selected[:, 0]
        j = selected[:, 1]
        omega1[i, j] = (1.0 - chi) * omega0[i, j]
        omega1[j, i] = omega1[i, j]
    # Important: keep the original strict diagonal fixed.
    np.fill_diagonal(omega1, diag0)

    offdiag_abs = np.sum(np.abs(omega1), axis=1) - np.abs(np.diag(omega1))
    diagnostics = {
        "b0": params.b,
        "theta0": params.theta0,
        "rho": params.rho,
        "delta0": params.delta0,
        "eta": eta,
        "chi": chi,
        "total_withinband_edges": int(total_edges),
        "number_weakened_edges": int(m_weakened),
        "diagonal_fixed": True,
        "diagonal_min": float(np.min(diag0)),
        "diagonal_median": float(np.median(diag0)),
        "diagonal_max": float(np.max(diag0)),
        "max_abs_offdiag_row_sum_after_weakening": float(np.max(offdiag_abs)),
        "strict_diagonal_dominance_margin_after_weakening": float(np.min(diag0 - offdiag_abs)),
    }
    return omega1, selected, diagnostics


def e0_edges_ordered(p: int, b0: int) -> np.ndarray:
    return np.asarray([(i, i + dist) for dist in range(1, b0 + 1) for i in range(p - dist)], dtype=np.int64)


def make_edge_masks(
    p: int,
    b0: int,
    selected_edges: np.ndarray,
    chi: float,
) -> dict[str, Any]:
    e0_edges = e0_edges_ordered(p, b0)
    selected_mask_cpu = np.zeros((p, p), dtype=bool)
    if selected_edges.size:
        si = selected_edges[:, 0]
        sj = selected_edges[:, 1]
        selected_mask_cpu[si, sj] = True
        selected_mask_cpu[sj, si] = True

    e0_mask_cpu = np.zeros((p, p), dtype=bool)
    ei = e0_edges[:, 0]
    ej = e0_edges[:, 1]
    e0_mask_cpu[ei, ej] = True
    e0_mask_cpu[ej, ei] = True

    true_nonzero_cpu = e0_mask_cpu.copy()
    deleted_cpu = np.zeros((p, p), dtype=bool)
    if abs(chi - 1.0) < 1e-12 and selected_edges.size:
        deleted_cpu[selected_mask_cpu] = True
        true_nonzero_cpu[selected_mask_cpu] = False

    unselected_edges = e0_edges[~selected_mask_cpu[ei, ej], :]
    return {
        "e0_edges": e0_edges,
        "selected_edges": selected_edges,
        "unselected_edges": unselected_edges,
        "weakened_mask_gpu": cp.asarray(selected_mask_cpu),
        "true_nonzero_mask_gpu": cp.asarray(true_nonzero_cpu),
        "deleted_mask_gpu": cp.asarray(deleted_cpu),
    }


def local_block_lrv_tau(score_local: cp.ndarray, block_length: int, n: int) -> cp.ndarray:
    cs = cp.cumsum(score_local, axis=0)
    block_sums = cs[block_length - 1 :, :].copy()
    if block_length < score_local.shape[0]:
        block_sums[1:, :] -= cs[:-block_length, :]
    normalized = block_sums / math.sqrt(block_length)
    lrv = cp.mean(normalized * normalized, axis=0)
    return cp.sqrt(cp.maximum(lrv, 0.0) / n)


def direct_and_ce_test_one_epoch_gpu(
    x: cp.ndarray,
    params: RunParams,
    masks: dict[str, Any],
    c_delta: float,
    delta_n: float,
) -> dict[str, float | int | bool]:
    n, p, b0, block_length = params.n, params.p, params.b, params.block_length
    sigma_hat = (x.T @ x) / n
    n_blocks = n - block_length + 1
    boot_block_max = cp.zeros(n_blocks, dtype=cp.float64)
    max_abs_direct = cp.asarray(0.0, dtype=cp.float64)
    max_abs_ce = cp.asarray(0.0, dtype=cp.float64)
    invalid_local_solves = 0

    c_n = c_delta * math.sqrt(math.log(p))
    eps = cp.asarray(1e-300, dtype=cp.float64)

    pumped_matrix = cp.zeros((p, p), dtype=cp.bool_)
    ratio_matrix = cp.full((p, p), cp.nan, dtype=cp.float64)

    for j in range(p):
        s_left = max(0, j - b0)
        s_right = min(p - 1, j + b0)
        s_idx = cp.arange(s_left, s_right + 1)
        local_j = j - s_left
        s_size = int(s_idx.size)

        sigma_ss = sigma_hat[s_idx[:, None], s_idx]
        rhs = cp.concatenate((cp.eye(s_size, dtype=cp.float64), sigma_hat[s_idx, :]), axis=1)
        try:
            local_sol = cp.linalg.solve(sigma_ss, rhs)
        except Exception:
            invalid_local_solves += 1
            continue
        if not bool(cp.all(cp.isfinite(local_sol)).item()):
            invalid_local_solves += 1
            continue

        inv_ss = local_sol[:, :s_size]
        c_j = local_sol[:, s_size:]
        omega_hat_j = inv_ss[:, local_j]

        # Direct statistic and direct-statistic bootstrap critical value.
        residual_direct = c_j[local_j, :].copy()
        residual_direct[s_idx] = 0.0
        max_abs_direct = cp.maximum(max_abs_direct, cp.max(cp.abs(residual_direct)))

        x_s = x[:, s_idx]
        xomega = x_s @ omega_hat_j
        score_mat = x - (x_s @ c_j)
        score_mat *= xomega[:, None]
        score_mat[:, s_idx] = 0.0
        boot_block_max = cp.maximum(boot_block_max, rolling_block_max_abs(score_mat, block_length))

        # Consistency-empowered pump on the local inverse coordinates.
        v_scores = x_s @ inv_ss
        score_local = v_scores * xomega[:, None] - omega_hat_j[None, :]
        tau_hat = local_block_lrv_tau(score_local, block_length, n)
        ratio = cp.abs(omega_hat_j) / cp.maximum(tau_hat, eps)
        pump = ratio <= delta_n
        omega_tilde_j = omega_hat_j + c_n * pump.astype(cp.float64)

        ce_residual = sigma_hat[:, s_idx] @ omega_tilde_j
        ce_residual[j] -= 1.0
        max_abs_ce = cp.maximum(max_abs_ce, cp.max(cp.abs(ce_residual)))

        pumped_matrix[s_idx, j] = pump
        ratio_matrix[s_idx, j] = ratio

        del (
            sigma_ss,
            rhs,
            local_sol,
            inv_ss,
            c_j,
            omega_hat_j,
            residual_direct,
            x_s,
            xomega,
            score_mat,
            v_scores,
            score_local,
            tau_hat,
            ratio,
            pump,
            omega_tilde_j,
            ce_residual,
        )

    q_hat = empirical_quantile_type1(boot_block_max / math.sqrt(block_length), 1.0 - params.alpha)
    t_hat = math.sqrt(n) * float(max_abs_direct.item())
    t_tilde = math.sqrt(n) * float(max_abs_ce.item())

    pumped_total = int(cp.sum(pumped_matrix).item())
    pumped_diagonal = int(cp.sum(cp.diag(pumped_matrix)).item())
    pumped_offdiag = pumped_total - pumped_diagonal

    selected_edges = masks["selected_edges"]
    unselected_edges = masks["unselected_edges"]
    if selected_edges.size:
        si = cp.asarray(selected_edges[:, 0])
        sj = cp.asarray(selected_edges[:, 1])
        pumped_weakened_edges = int(cp.sum(cp.logical_or(pumped_matrix[si, sj], pumped_matrix[sj, si])).item())
    else:
        pumped_weakened_edges = 0
    if unselected_edges.size:
        ui = cp.asarray(unselected_edges[:, 0])
        uj = cp.asarray(unselected_edges[:, 1])
        pumped_unweakened_edges = int(cp.sum(cp.logical_or(pumped_matrix[ui, uj], pumped_matrix[uj, ui])).item())
    else:
        pumped_unweakened_edges = 0

    true_ratios = ratio_matrix[masks["true_nonzero_mask_gpu"]]
    true_ratios = true_ratios[cp.isfinite(true_ratios)]
    if int(true_ratios.size) > 0:
        min_true_ratio = float(cp.min(true_ratios).item())
        median_true_ratio = float(cp.median(true_ratios).item())
    else:
        min_true_ratio = float("nan")
        median_true_ratio = float("nan")

    deleted_ratios = ratio_matrix[masks["deleted_mask_gpu"]]
    deleted_ratios = deleted_ratios[cp.isfinite(deleted_ratios)]
    max_deleted_ratio = float(cp.max(deleted_ratios).item()) if int(deleted_ratios.size) > 0 else float("nan")

    return {
        "T_tilde": t_tilde,
        "T_hat": t_hat,
        "q_hat": q_hat,
        "reject_CE": bool(t_tilde > q_hat),
        "reject_direct": bool(t_hat > q_hat),
        "number_pumped_total": pumped_total,
        "number_pumped_offdiag": pumped_offdiag,
        "number_pumped_diagonal": pumped_diagonal,
        "number_pumped_weakened_edges": pumped_weakened_edges,
        "number_pumped_unweakened_edges": pumped_unweakened_edges,
        "min_ratio_abs_omega_over_tau_true_nonzero": min_true_ratio,
        "median_ratio_abs_omega_over_tau_true_nonzero": median_true_ratio,
        "max_ratio_abs_omega_over_tau_deleted": max_deleted_ratio,
        "invalid_local_solves": invalid_local_solves,
    }


def combo_stem(params: RunParams, eta: float, chi: float, run_id: str) -> str:
    return (
        f"withinband_ce_n{params.n}_p{params.p}_beta{params.beta:g}_s0{params.s0:02d}"
        f"_eta{eta:g}_chi{chi:g}_Cdelta{run_id}"
    )


def summarize_epoch_rows(epoch_rows: list[dict[str, Any]]) -> dict[str, Any]:
    ce_rejects = np.asarray([bool(row["reject_CE"]) for row in epoch_rows], dtype=bool)
    direct_rejects = np.asarray([bool(row["reject_direct"]) for row in epoch_rows], dtype=bool)
    numeric_fields = [
        "T_tilde",
        "T_hat",
        "q_hat",
        "number_pumped_total",
        "number_pumped_offdiag",
        "number_pumped_diagonal",
        "number_pumped_weakened_edges",
        "number_pumped_unweakened_edges",
        "min_ratio_abs_omega_over_tau_true_nonzero",
        "median_ratio_abs_omega_over_tau_true_nonzero",
        "max_ratio_abs_omega_over_tau_deleted",
        "invalid_local_solves",
    ]
    summary: dict[str, Any] = {
        "ce_rejection_rate": float(np.mean(ce_rejects)),
        "direct_rejection_rate": float(np.mean(direct_rejects)),
        "ce_reject_count": int(np.sum(ce_rejects)),
        "direct_reject_count": int(np.sum(direct_rejects)),
    }
    for field in numeric_fields:
        vals = np.asarray([float(row[field]) for row in epoch_rows], dtype=np.float64)
        summary[f"mean_{field}"] = float(np.nanmean(vals))
        summary[f"sd_{field}"] = float(np.nanstd(vals, ddof=1)) if vals.size > 1 else 0.0
    return summary


def run_one_combo(
    args: argparse.Namespace,
    n: int,
    p: int,
    beta: float,
    block_length: int,
    eta: float,
    chi: float,
    run_id: str,
    combo_index: int,
    total_combos: int,
) -> dict[str, Any]:
    params = make_params(args, n=n, p=p, beta=beta, block_length=block_length)
    out_dir = Path(params.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ordered_edges = ordered_withinband_edges(params.p, params.b, params.edge_seed)
    edge_file = out_dir / f"withinband_ordered_edges_p{params.p}_b{params.b}_seed{params.edge_seed}.csv"
    if not edge_file.exists():
        write_edges_csv(edge_file, ordered_edges)

    omega, selected_edges, precision_diagnostics = build_withinband_precision(params, eta, chi, ordered_edges)
    masks = make_edge_masks(params.p, params.b, selected_edges, chi)
    omega_gpu = cp.asarray(omega)
    cholesky_success = False
    r_upper = cp.linalg.cholesky(omega_gpu).T
    cholesky_success = True
    b_fft = prepare_fft_filter(params)
    cp.random.seed(params.seed + 100000 * combo_index)

    if args.delta_n_constant is not None and args.delta_n_constant > 0:
        delta_n_value = float(args.delta_n_constant)
        delta_n_mode = "constant"
    else:
        delta_n_value = math.sqrt(math.log(params.n))
        delta_n_mode = "sqrt_log_n"

    epoch_rows: list[dict[str, Any]] = []
    combo_start = time.perf_counter()
    epoch_start = 1
    while epoch_start <= params.epochs:
        epoch_end = min(params.epochs, epoch_start + params.batch_epochs - 1)
        batch_size = epoch_end - epoch_start + 1
        rows_needed = params.n * batch_size
        x_batch = generate_observed_batch(rows_needed, params, b_fft, r_upper)

        for offset in range(batch_size):
            global_epoch = epoch_start + offset
            row_left = offset * params.n
            row_right = (offset + 1) * params.n
            test_out = direct_and_ce_test_one_epoch_gpu(
                x_batch[row_left:row_right, :],
                params=params,
                masks=masks,
                c_delta=args.C_delta,
                delta_n=delta_n_value,
            )
            epoch_rows.append({"epoch": global_epoch, **test_out})
            if global_epoch == params.epochs or global_epoch % max(1, params.batch_epochs) == 0:
                elapsed = time.perf_counter() - combo_start
                print(
                    f"  combo {combo_index}/{total_combos} "
                    f"n={n}, p={p}, beta={beta:g}, eta={eta:g}, chi={chi:g}: "
                    f"epoch {global_epoch}/{params.epochs}, elapsed {format_seconds(elapsed)}",
                    flush=True,
                )

        del x_batch
        cp.get_default_memory_pool().free_all_blocks()
        epoch_start = epoch_end + 1

    cp.cuda.Stream.null.synchronize()
    runtime_sec = time.perf_counter() - combo_start
    summary = summarize_epoch_rows(epoch_rows)
    summary["runtime_sec"] = float(runtime_sec)
    summary["cholesky_success"] = cholesky_success

    stem = combo_stem(params, eta=eta, chi=chi, run_id=run_id)
    epoch_csv = out_dir / f"{stem}_epochs.csv"
    result_json = out_dir / f"{stem}.json"
    epoch_fields = [
        "epoch",
        "T_tilde",
        "T_hat",
        "q_hat",
        "reject_CE",
        "reject_direct",
        "number_pumped_total",
        "number_pumped_offdiag",
        "number_pumped_diagonal",
        "number_pumped_weakened_edges",
        "number_pumped_unweakened_edges",
        "min_ratio_abs_omega_over_tau_true_nonzero",
        "median_ratio_abs_omega_over_tau_true_nonzero",
        "max_ratio_abs_omega_over_tau_deleted",
        "invalid_local_solves",
    ]
    write_csv(epoch_csv, epoch_rows, epoch_fields)

    result_obj = {
        "run_id": run_id,
        "backend": "cupy",
        "design": "withinband_weakening_ce",
        "params": asdict(params),
        "eta": eta,
        "chi": chi,
        "C_delta": args.C_delta,
        "delta_n": delta_n_value,
        "delta_n_mode": delta_n_mode,
        "delta_n_default_sqrt_log_n": math.sqrt(math.log(params.n)),
        "C_n": args.C_delta * math.sqrt(math.log(params.p)),
        "precision_diagnostics": {**precision_diagnostics, "cholesky_success": cholesky_success},
        "ordered_edge_file": str(edge_file),
        "epoch_results_csv": str(epoch_csv),
        "summary": summary,
        "package_versions": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "cupy": cp.__version__,
            "cuda_device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        },
    }
    save_json(result_json, result_obj)

    del omega_gpu, r_upper, b_fft
    cp.get_default_memory_pool().free_all_blocks()

    return {
        "combo_index": combo_index,
        "n": n,
        "p": p,
        "beta": beta,
        "block_length": block_length,
        "eta": eta,
        "chi": chi,
        "C_delta": args.C_delta,
        "number_weakened_edges": precision_diagnostics["number_weakened_edges"],
        **summary,
        "status": "ok",
        "error_message": "",
        "result_json": str(result_json),
        "epoch_csv": str(epoch_csv),
    }


def run_grid(args: argparse.Namespace) -> dict[str, Any]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.run_label:
        run_id = f"{run_id}_{args.run_label}"

    regimes = parse_regime_grid(args.regime_grid)
    eta_grid = parse_float_grid(args.eta_grid)
    chi_grid = parse_float_grid(args.chi_grid)
    all_combos = [
        (n, p, beta, block_length, eta, chi)
        for (n, p, beta, block_length) in regimes
        for eta in eta_grid
        for chi in chi_grid
    ]
    if args.combo_stride < 1:
        raise ValueError("--combo_stride must be at least 1")
    if args.combo_offset < 0 or args.combo_offset >= args.combo_stride:
        raise ValueError("--combo_offset must be in [0, combo_stride)")
    combos = [
        (global_idx, *combo)
        for global_idx, combo in enumerate(all_combos, start=1)
        if ((global_idx - 1) % args.combo_stride) == args.combo_offset
    ]
    if not combos:
        raise ValueError("No grid cells assigned to this worker")

    device_name = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    print(f"Starting within-band weakening CE grid {run_id} on {device_name}")
    print(
        f"Base config: s0={args.s0}, theta0={args.theta0:g}, rho={args.rho:g}, "
        f"delta0={args.delta0:g}, epochs={args.epochs}, batch_epochs={args.batch_epochs}, "
        f"alpha={args.alpha:g}, C_delta={args.C_delta:g}"
    )
    print(f"Regimes: {regimes}")
    print(f"Grid: eta={eta_grid}, chi={chi_grid}")
    print(
        f"Worker assignment: offset={args.combo_offset}, stride={args.combo_stride}, "
        f"assigned_cells={len(combos)}/{len(all_combos)}"
    )

    rows: list[dict[str, Any]] = []
    grid_start = time.perf_counter()
    for worker_idx, (global_idx, n, p, beta, block_length, eta, chi) in enumerate(combos, start=1):
        try:
            print(
                f"\n[worker {worker_idx}/{len(combos)} | global {global_idx}/{len(all_combos)}] "
                f"Running n={n}, p={p}, beta={beta:g}, block={block_length}, eta={eta:g}, chi={chi:g}",
                flush=True,
            )
            row = run_one_combo(
                args=args,
                n=n,
                p=p,
                beta=beta,
                block_length=block_length,
                eta=eta,
                chi=chi,
                run_id=run_id,
                combo_index=global_idx,
                total_combos=len(all_combos),
            )
        except Exception as exc:
            row = {
                "combo_index": global_idx,
                "n": n,
                "p": p,
                "beta": beta,
                "block_length": block_length,
                "eta": eta,
                "chi": chi,
                "C_delta": args.C_delta,
                "number_weakened_edges": float("nan"),
                "ce_rejection_rate": float("nan"),
                "direct_rejection_rate": float("nan"),
                "ce_reject_count": 0,
                "direct_reject_count": 0,
                "mean_T_tilde": float("nan"),
                "mean_T_hat": float("nan"),
                "mean_q_hat": float("nan"),
                "mean_number_pumped_total": float("nan"),
                "mean_number_pumped_weakened_edges": float("nan"),
                "runtime_sec": float("nan"),
                "cholesky_success": False,
                "status": "error",
                "error_message": repr(exc),
                "result_json": "",
                "epoch_csv": "",
            }
            print(f"  ERROR: {exc!r}", flush=True)
        rows.append(row)
        elapsed = time.perf_counter() - grid_start
        eta_remaining = elapsed / worker_idx * (len(combos) - worker_idx)
        print(
            f"  Finished worker cell {worker_idx}/{len(combos)} "
            f"(global {global_idx}/{len(all_combos)}). Total elapsed {format_seconds(elapsed)}. "
            f"ETA {format_seconds(eta_remaining)}.",
            flush=True,
        )

    total_runtime_sec = time.perf_counter() - grid_start
    grid_csv = out_dir / f"withinband_grid_summary_{run_id}.csv"
    grid_json = out_dir / f"withinband_grid_summary_{run_id}.json"
    grid_txt = out_dir / f"withinband_grid_summary_{run_id}.txt"
    fieldnames = [
        "combo_index",
        "n",
        "p",
        "beta",
        "block_length",
        "eta",
        "chi",
        "C_delta",
        "number_weakened_edges",
        "ce_rejection_rate",
        "direct_rejection_rate",
        "ce_reject_count",
        "direct_reject_count",
        "mean_T_tilde",
        "sd_T_tilde",
        "mean_T_hat",
        "sd_T_hat",
        "mean_q_hat",
        "sd_q_hat",
        "mean_number_pumped_total",
        "sd_number_pumped_total",
        "mean_number_pumped_offdiag",
        "sd_number_pumped_offdiag",
        "mean_number_pumped_diagonal",
        "sd_number_pumped_diagonal",
        "mean_number_pumped_weakened_edges",
        "sd_number_pumped_weakened_edges",
        "mean_number_pumped_unweakened_edges",
        "sd_number_pumped_unweakened_edges",
        "mean_min_ratio_abs_omega_over_tau_true_nonzero",
        "sd_min_ratio_abs_omega_over_tau_true_nonzero",
        "mean_median_ratio_abs_omega_over_tau_true_nonzero",
        "sd_median_ratio_abs_omega_over_tau_true_nonzero",
        "mean_max_ratio_abs_omega_over_tau_deleted",
        "sd_max_ratio_abs_omega_over_tau_deleted",
        "mean_invalid_local_solves",
        "sd_invalid_local_solves",
        "runtime_sec",
        "cholesky_success",
        "status",
        "error_message",
        "result_json",
        "epoch_csv",
    ]
    # Ensure every row has every field, even for caught errors.
    for row in rows:
        for field in fieldnames:
            row.setdefault(field, "")
    write_csv(grid_csv, rows, fieldnames)
    grid_obj = {
        "run_id": run_id,
        "backend": "cupy",
        "design": "withinband_weakening_ce",
        "total_runtime_sec": total_runtime_sec,
        "total_runtime_hms": format_seconds(total_runtime_sec),
        "base_args": vars(args),
        "grid_summary_csv": str(grid_csv),
        "grid_summary": rows,
        "package_versions": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "cupy": cp.__version__,
            "cuda_device": device_name,
        },
    }
    save_json(grid_json, grid_obj)
    with grid_txt.open("w") as f:
        f.write(f"Within-band weakening CE GPU grid run {run_id}\n")
        f.write(f"backend=cupy\n")
        f.write(f"cuda_device={device_name}\n")
        f.write(f"total_runtime_sec={total_runtime_sec:.3f}\n")
        f.write(f"total_runtime_hms={format_seconds(total_runtime_sec)}\n")
        f.write(f"grid_summary_csv={grid_csv}\n")
        f.write(f"grid_summary_json={grid_json}\n")

    print("\nWithin-band weakening CE grid complete.")
    print(f"Total runtime: {format_seconds(total_runtime_sec)} ({total_runtime_sec:.3f} sec)")
    print(f"Saved grid CSV: {grid_csv}")
    print(f"Saved grid JSON: {grid_json}")
    print(f"Saved grid TXT: {grid_txt}")
    return grid_obj


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="GPU/CuPy within-band weakening CE experiment.")
    parser.add_argument(
        "--regime_grid",
        type=str,
        default="2000:800:2:160,2000:800:0.9:238,800:2000:2:86,800:2000:0.9:86",
    )
    parser.add_argument("--s0", type=int, default=10)
    parser.add_argument("--theta0", type=float, default=THETA0_DEFAULT)
    parser.add_argument("--rho", type=float, default=RHO_DEFAULT)
    parser.add_argument("--delta0", type=float, default=DELTA0_DEFAULT)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--alpha", type=float, default=0.10)
    parser.add_argument("--seed", type=int, default=101)
    parser.add_argument("--batch_epochs", type=int, default=5)
    parser.add_argument("--fft_len", type=int, default=None)
    parser.add_argument("--eta_grid", type=str, default="0.02,0.05,0.10,0.20,0.30")
    parser.add_argument("--chi_grid", type=str, default="0.20,0.40,0.60,0.80,1.00")
    parser.add_argument("--C_delta", type=float, default=1.0)
    parser.add_argument("--edge_seed", type=int, default=202)
    parser.add_argument("--gpu_chunk_cols", type=int, default=4)
    parser.add_argument(
        "--out_dir",
        type=str,
        default=str(EXPERIMENT_DIR / "Results" / "raw" / "withinband_ce_gpu"),
    )
    parser.add_argument("--combo_offset", type=int, default=0)
    parser.add_argument("--combo_stride", type=int, default=1)
    parser.add_argument("--run_label", type=str, default="")
    parser.add_argument(
        "--delta_n_constant",
        type=float,
        default=None,
        help="Use this constant in place of sqrt(log n) as the CE pump threshold delta_n. "
             "If unset (default), delta_n = sqrt(log n) per Theorem 4.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    run_grid(args)


if __name__ == "__main__":
    main()
