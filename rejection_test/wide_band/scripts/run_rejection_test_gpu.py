#!/usr/bin/env python3

"""GPU/CuPy runner for the TGG precision FFT rejection experiment.

This is intentionally separate from the R implementation.  It keeps the same
DGP and statistic, but uses CuPy for the FFT generation, dense Cholesky solve,
and the per-epoch TGG/bootstrap matrix operations.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import cupy as cp
import numpy as np


THETA0_DEFAULT = 0.25
RHO_DEFAULT = 0.5
DELTA0_DEFAULT = 0.10


@dataclass
class RunParams:
    n: int
    p: int
    beta: float
    s0: int
    b: int
    m_extra: int
    theta_extra: float
    h: float
    epochs: int
    alpha: float
    seed: int
    batch_epochs: int
    block_length: int
    filter_length: int
    theta0: float
    rho: float
    delta0: float
    out_dir: str
    gpu_chunk_cols: int
    diag_mode: str
    signal_mode: str
    kappa: float
    edge_seed: int
    sign_seed: int
    edge_sampling: str
    degree_cap: int | None
    kappa_tol: float
    kappa_max_iter: int
    kappa_damping: float


def parse_float_grid(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def parse_int_grid(text: str) -> list[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def format_seconds(seconds: float) -> str:
    seconds = int(round(seconds))
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60
    return f"{h:02d}:{m:02d}:{s:02d}"


def make_params(args: argparse.Namespace, m_extra: int, theta_extra: float, kappa: float = float("nan")) -> RunParams:
    b = (args.s0 - 1) // 2
    block_length = args.block_length
    if block_length is None:
        block_length = int(math.floor(args.n ** (2.0 / 3.0)))
    if block_length < 1 or block_length > args.n:
        raise ValueError("block_length must be between 1 and n")

    scale = math.sqrt(math.log(args.p) / args.n)
    h = theta_extra / scale if math.isfinite(theta_extra) else float("nan")
    return RunParams(
        n=args.n,
        p=args.p,
        beta=args.beta,
        s0=args.s0,
        b=b,
        m_extra=m_extra,
        theta_extra=theta_extra,
        h=h,
        epochs=args.epochs,
        alpha=args.alpha,
        seed=args.seed,
        batch_epochs=args.batch_epochs,
        block_length=block_length,
        filter_length=args.fft_len if args.fft_len is not None else args.n * args.n,
        theta0=args.theta0,
        rho=args.rho,
        delta0=args.delta0,
        out_dir=args.out_dir,
        gpu_chunk_cols=args.gpu_chunk_cols,
        diag_mode=args.diag_mode,
        signal_mode=args.signal_mode,
        kappa=kappa,
        edge_seed=args.edge_seed,
        sign_seed=args.sign_seed,
        edge_sampling=args.edge_sampling,
        degree_cap=args.degree_cap,
        kappa_tol=args.kappa_tol,
        kappa_max_iter=args.kappa_max_iter,
        kappa_damping=args.kappa_damping,
    )


def sample_extra_edges(p: int, b: int, m_extra: int, rng: np.random.Generator) -> np.ndarray:
    if m_extra == 0:
        return np.empty((0, 3), dtype=np.int64)

    total_pairs = p * (p - 1) // 2
    excluded_band = sum(p - dist for dist in range(1, b + 1))
    allowed_pairs = total_pairs - excluded_band
    if m_extra > allowed_pairs:
        raise ValueError(f"m_extra={m_extra} exceeds available off-band edges={allowed_pairs}")

    cum_counts = np.concatenate(([0], np.cumsum(np.arange(p - 1, 0, -1, dtype=np.int64))))
    seen: set[tuple[int, int]] = set()
    rows: list[tuple[int, int, int]] = []
    while len(rows) < m_extra:
        remaining = m_extra - len(rows)
        chunk_size = max(1000, 5 * remaining)
        pair_ids = rng.integers(0, total_pairs, size=chunk_size, dtype=np.int64)
        i_vals = np.searchsorted(cum_counts, pair_ids, side="right") - 1
        j_vals = i_vals + 1 + (pair_ids - cum_counts[i_vals])
        keep = np.flatnonzero((j_vals - i_vals) > b)

        for idx in keep:
            i = int(i_vals[idx])
            j = int(j_vals[idx])
            if (i, j) in seen:
                continue
            seen.add((i, j))
            rows.append((i, j, 0))
            if len(rows) == m_extra:
                break

    signs = rng.choice(np.array([-1, 1], dtype=np.int64), size=m_extra, replace=True)
    rows = [(i, j, int(signs[k])) for k, (i, j, _) in enumerate(rows)]
    return np.asarray(rows, dtype=np.int64)


def sample_edge_pairs_uniform(p: int, b: int, m_extra: int, rng: np.random.Generator) -> list[tuple[int, int]]:
    if m_extra == 0:
        return []

    total_pairs = p * (p - 1) // 2
    excluded_band = sum(p - dist for dist in range(1, b + 1))
    allowed_pairs = total_pairs - excluded_band
    if m_extra > allowed_pairs:
        raise ValueError(f"m_extra={m_extra} exceeds available off-band edges={allowed_pairs}")

    cum_counts = np.concatenate(([0], np.cumsum(np.arange(p - 1, 0, -1, dtype=np.int64))))
    seen: set[tuple[int, int]] = set()
    rows: list[tuple[int, int]] = []
    while len(rows) < m_extra:
        remaining = m_extra - len(rows)
        chunk_size = max(1000, 5 * remaining)
        pair_ids = rng.integers(0, total_pairs, size=chunk_size, dtype=np.int64)
        i_vals = np.searchsorted(cum_counts, pair_ids, side="right") - 1
        j_vals = i_vals + 1 + (pair_ids - cum_counts[i_vals])
        keep = np.flatnonzero((j_vals - i_vals) > b)
        for idx in keep:
            i = int(i_vals[idx])
            j = int(j_vals[idx])
            if (i, j) in seen:
                continue
            seen.add((i, j))
            rows.append((i, j))
            if len(rows) == m_extra:
                break
    return rows


def sample_edge_pairs_degree_capped(
    p: int,
    b: int,
    m_extra: int,
    degree_cap: int,
    rng: np.random.Generator,
) -> list[tuple[int, int]]:
    if m_extra == 0:
        return []
    if degree_cap <= 0:
        raise ValueError("degree_capped sampling requires positive degree_cap")
    if 2 * m_extra > p * degree_cap:
        raise ValueError("degree_cap is too small for requested m_extra")

    degree = np.zeros(p, dtype=np.int64)
    seen: set[tuple[int, int]] = set()
    rows: list[tuple[int, int]] = []
    attempts = 0
    max_attempts = max(1_000_000, 2000 * m_extra)
    while len(rows) < m_extra and attempts < max_attempts:
        attempts += 1
        i = int(rng.integers(0, p - 1))
        j = int(rng.integers(i + 1, p))
        if j - i <= b:
            continue
        if degree[i] >= degree_cap or degree[j] >= degree_cap:
            continue
        if (i, j) in seen:
            continue
        seen.add((i, j))
        degree[i] += 1
        degree[j] += 1
        rows.append((i, j))
    if len(rows) < m_extra:
        raise ValueError(f"Only sampled {len(rows)}/{m_extra} degree-capped edges")
    return rows


def sample_ordered_nested_edges(
    p: int,
    b: int,
    m_max: int,
    edge_seed: int,
    sign_seed: int,
    edge_sampling: str,
    degree_cap: int | None,
) -> np.ndarray:
    edge_rng = np.random.default_rng(edge_seed)
    sign_rng = np.random.default_rng(sign_seed)
    if edge_sampling == "uniform_nested":
        pairs = sample_edge_pairs_uniform(p, b, m_max, edge_rng)
    elif edge_sampling == "degree_capped":
        if degree_cap is None:
            raise ValueError("degree_capped sampling requires --degree_cap")
        pairs = sample_edge_pairs_degree_capped(p, b, m_max, degree_cap, edge_rng)
    else:
        raise ValueError(f"Unknown edge_sampling={edge_sampling!r}")
    signs = sign_rng.choice(np.array([-1, 1], dtype=np.int64), size=m_max, replace=True)
    rows = [(i, j, int(signs[k])) for k, (i, j) in enumerate(pairs)]
    return np.asarray(rows, dtype=np.int64)


def baseline_precision_components(params: RunParams) -> tuple[np.ndarray, np.ndarray]:
    omega = np.zeros((params.p, params.p), dtype=np.float64)
    row_abs = np.zeros(params.p, dtype=np.float64)
    for dist in range(1, params.b + 1):
        val = params.theta0 * (params.rho ** (dist - 1))
        idx = np.arange(params.p - dist)
        omega[idx, idx + dist] = val
        omega[idx + dist, idx] = val
        row_abs[idx] += abs(val)
        row_abs[idx + dist] += abs(val)
    d0 = 1.0 + row_abs + params.delta0
    return omega, d0


def extra_adjacency_times_vec(p: int, extra_edges: np.ndarray, v: np.ndarray) -> np.ndarray:
    out = np.zeros(p, dtype=np.float64)
    if extra_edges.shape[0] == 0:
        return out
    i = extra_edges[:, 0].astype(np.int64)
    j = extra_edges[:, 1].astype(np.int64)
    np.add.at(out, i, v[j])
    np.add.at(out, j, v[i])
    return out


def estimate_extra_adjacency_radius(p: int, extra_edges: np.ndarray, max_iter: int = 2000, tol: float = 1e-10) -> float:
    if extra_edges.shape[0] == 0:
        return 0.0
    v = np.full(p, 1.0 / math.sqrt(p), dtype=np.float64)
    lambda_old = 0.0
    for _ in range(max_iter):
        av = extra_adjacency_times_vec(p, extra_edges, v)
        norm = float(np.linalg.norm(av))
        if not math.isfinite(norm) or norm == 0.0:
            return 0.0
        v = av / norm
        lambda_new = float(v @ extra_adjacency_times_vec(p, extra_edges, v))
        if abs(lambda_new - lambda_old) / max(1.0, abs(lambda_old)) < tol:
            return lambda_new
        lambda_old = lambda_new
    return lambda_old


def solve_kappa_fixed_point(params: RunParams, d0: np.ndarray, extra_edges: np.ndarray) -> dict[str, object]:
    y = np.sqrt(d0)
    if extra_edges.shape[0] == 0 or params.kappa == 0.0:
        return {"y": y, "D": d0.copy(), "converged": True, "iterations": 0, "error": 0.0, "diverged": False}

    err = float("inf")
    converged = False
    diverged = False
    for iteration in range(1, params.kappa_max_iter + 1):
        s = extra_adjacency_times_vec(params.p, extra_edges, y)
        y_target = 0.5 * (params.kappa * s + np.sqrt((params.kappa * s) ** 2 + 4.0 * d0))
        y_new = (1.0 - params.kappa_damping) * y + params.kappa_damping * y_target
        err = float(np.max(np.abs(y_new - y) / np.maximum(1.0, np.abs(y))))
        if not np.all(np.isfinite(y_new)) or float(np.max(y_new)) > 1e12:
            diverged = True
            y = y_new
            break
        y = y_new
        if err < params.kappa_tol:
            converged = True
            break
    return {
        "y": y,
        "D": y * y,
        "converged": converged,
        "iterations": iteration,
        "error": err,
        "diverged": diverged,
    }


def edge_and_precision_diagnostics(params: RunParams, extra_edges: np.ndarray, diag_vals: np.ndarray, extra_abs: np.ndarray) -> dict[str, object]:
    degree = np.zeros(params.p, dtype=np.int64)
    if extra_edges.shape[0] > 0:
        np.add.at(degree, extra_edges[:, 0], 1)
        np.add.at(degree, extra_edges[:, 1], 1)
        actual_kappa = extra_abs / np.sqrt(diag_vals[extra_edges[:, 0]] * diag_vals[extra_edges[:, 1]])
    else:
        actual_kappa = np.asarray([], dtype=np.float64)
    rho_extra = estimate_extra_adjacency_radius(params.p, extra_edges)
    return {
        "m_extra": int(extra_edges.shape[0]),
        "mean_added_degree": float(np.mean(degree)),
        "max_added_degree": int(np.max(degree)),
        "degree_quantiles": [float(x) for x in np.quantile(degree, [0, 0.25, 0.5, 0.75, 0.95, 1.0])],
        "number_touched_nodes": int(np.sum(degree > 0)),
        "edge_sampling": params.edge_sampling,
        "degree_cap": params.degree_cap,
        "rho_extra_adj": float(rho_extra),
        "kappa_times_rho": float(params.kappa * rho_extra) if params.signal_mode == "kappa" else float("nan"),
        "actual_kappa_min": float(np.min(actual_kappa)) if actual_kappa.size else float("nan"),
        "actual_kappa_median": float(np.median(actual_kappa)) if actual_kappa.size else float("nan"),
        "actual_kappa_mean": float(np.mean(actual_kappa)) if actual_kappa.size else float("nan"),
        "actual_kappa_max": float(np.max(actual_kappa)) if actual_kappa.size else float("nan"),
        "actual_raw_x_min": float(np.min(extra_abs)) if extra_abs.size else float("nan"),
        "actual_raw_x_median": float(np.median(extra_abs)) if extra_abs.size else float("nan"),
        "actual_raw_x_mean": float(np.mean(extra_abs)) if extra_abs.size else float("nan"),
        "actual_raw_x_max": float(np.max(extra_abs)) if extra_abs.size else float("nan"),
        "diag_min": float(np.min(diag_vals)),
        "diag_median": float(np.median(diag_vals)),
        "diag_mean": float(np.mean(diag_vals)),
        "diag_max": float(np.max(diag_vals)),
    }


def apply_extra_diagonal_updates(
    row_abs: np.ndarray,
    extra_i: np.ndarray,
    extra_j: np.ndarray,
    abs_x: np.ndarray,
    diag_mode: str,
) -> None:
    if diag_mode == "strict":
        np.add.at(row_abs, extra_i, abs_x)
        np.add.at(row_abs, extra_j, abs_x)
        return

    if diag_mode == "r_legacy":
        # Match the current R code's duplicated-index subassignment behavior.
        rhs_i = row_abs[extra_i] + abs_x
        for idx, value in zip(extra_i, rhs_i):
            row_abs[idx] = value
        rhs_j = row_abs[extra_j] + abs_x
        for idx, value in zip(extra_j, rhs_j):
            row_abs[idx] = value
        return

    raise ValueError(f"Unknown diag_mode={diag_mode!r}")


def build_precision_matrix(params: RunParams, extra_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    omega, d0 = baseline_precision_components(params)
    diag_vals = d0.copy()
    extra_abs = np.zeros(extra_edges.shape[0], dtype=np.float64)
    fixed_point = None

    if params.signal_mode == "kappa":
        fixed_point = solve_kappa_fixed_point(params, d0, extra_edges)
        if not fixed_point["converged"]:
            raise RuntimeError(
                f"kappa fixed point did not converge: iterations={fixed_point['iterations']} "
                f"error={fixed_point['error']}"
            )
        y = fixed_point["y"]
        diag_vals = fixed_point["D"]
        if extra_edges.shape[0] > 0:
            extra_abs = params.kappa * y[extra_edges[:, 0]] * y[extra_edges[:, 1]]
    else:
        row_abs = d0 - 1.0 - params.delta0
        if extra_edges.shape[0] > 0:
            extra_abs = np.full(extra_edges.shape[0], abs(params.theta_extra), dtype=np.float64)
            apply_extra_diagonal_updates(
                row_abs=row_abs,
                extra_i=extra_edges[:, 0].astype(np.int64),
                extra_j=extra_edges[:, 1].astype(np.int64),
                abs_x=extra_abs,
                diag_mode=params.diag_mode,
            )
        diag_vals = 1.0 + row_abs + params.delta0

    for pos, (i, j, sign) in enumerate(extra_edges):
        val = extra_abs[pos] * sign
        omega[i, j] = val
        omega[j, i] = val

    np.fill_diagonal(omega, diag_vals)
    omega_diag = edge_and_precision_diagnostics(params, extra_edges, diag_vals, extra_abs)
    if fixed_point is not None:
        omega_diag.update(
            {
                "fixed_point_converged": bool(fixed_point["converged"]),
                "fixed_point_iterations": int(fixed_point["iterations"]),
                "fixed_point_error": float(fixed_point["error"]),
                "fixed_point_diverged": bool(fixed_point["diverged"]),
            }
        )
    return omega, extra_abs, omega_diag


def combo_label(params: RunParams) -> str:
    if params.signal_mode == "kappa":
        return f"m={params.m_extra}, kappa={params.kappa:g}"
    return f"m={params.m_extra}, theta={params.theta_extra:g}"


def combo_stem(params: RunParams, run_id: str) -> str:
    if params.signal_mode == "kappa":
        signal_part = f"kappa{params.kappa:g}"
    else:
        signal_part = f"theta{params.theta_extra:g}"
    return (
        f"gpu_tgg_n{params.n}_p{params.p}_beta{params.beta:g}_s0{params.s0:02d}"
        f"_m{params.m_extra}_{params.signal_mode}_{signal_part}_seed{params.seed}_{run_id}"
    )


def prepare_fft_filter(params: RunParams) -> cp.ndarray:
    n_fft = params.filter_length
    a_vec = cp.arange(1, n_fft + 1, dtype=cp.float64) ** (-params.beta)
    a_vec = a_vec / cp.sqrt(cp.sum(a_vec * a_vec))
    return cp.fft.fft(a_vec) / math.sqrt(n_fft)


def generate_temporal_batch(rows_needed: int, p: int, b_fft: cp.ndarray, chunk_cols: int) -> cp.ndarray:
    n_fft = int(b_fft.shape[0])
    y_batch = cp.empty((rows_needed, p), dtype=cp.float64)
    sqrt_n_fft = math.sqrt(n_fft)

    for start in range(0, p, chunk_cols):
        end = min(p, start + chunk_cols)
        cols = end - start
        eps = cp.random.standard_normal((n_fft, cols), dtype=cp.float64)
        f_eps = cp.fft.fft(eps, axis=0) / sqrt_n_fft
        y = cp.fft.ifft(b_fft[:, None] * f_eps, axis=0) * n_fft
        y_batch[:, start:end] = cp.real(y[:rows_needed, :])
        del eps, f_eps, y
    return y_batch


def generate_observed_batch(rows_needed: int, params: RunParams, b_fft: cp.ndarray, r_upper: cp.ndarray) -> cp.ndarray:
    y_batch = generate_temporal_batch(
        rows_needed=rows_needed,
        p=params.p,
        b_fft=b_fft,
        chunk_cols=params.gpu_chunk_cols,
    )
    x_t = cp.linalg.solve(r_upper, y_batch.T)
    x_batch = x_t.T
    del y_batch, x_t
    return x_batch


def empirical_quantile_type1(x: cp.ndarray, prob: float) -> float:
    flat = cp.sort(x.ravel())
    n = int(flat.size)
    idx = max(0, min(n - 1, int(math.ceil(prob * n)) - 1))
    return float(flat[idx].item())


def rolling_block_max_abs(score_mat: cp.ndarray, block_length: int) -> cp.ndarray:
    cs = cp.cumsum(score_mat, axis=0)
    block_sums = cs[block_length - 1 :, :].copy()
    if block_length < score_mat.shape[0]:
        block_sums[1:, :] -= cs[:-block_length, :]
    return cp.max(cp.abs(block_sums), axis=1)


def tgg_test_one_epoch_gpu(x: cp.ndarray, params: RunParams) -> dict[str, float | int | bool]:
    n, p, b, l = params.n, params.p, params.b, params.block_length
    sigma_hat = (x.T @ x) / n
    n_blocks = n - l + 1
    boot_block_max = cp.zeros(n_blocks, dtype=cp.float64)
    max_abs_residual = cp.asarray(0.0, dtype=cp.float64)
    invalid_local_solves = 0

    for i in range(p):
        s_left = max(0, i - b)
        s_right = min(p - 1, i + b)
        s_idx = cp.arange(s_left, s_right + 1)
        local_i = i - s_left

        sigma_ss = sigma_hat[s_idx[:, None], s_idx]
        rhs_left = cp.zeros((s_idx.size, 1), dtype=cp.float64)
        rhs_left[local_i, 0] = 1.0
        rhs = cp.concatenate((rhs_left, sigma_hat[s_idx, :]), axis=1)

        try:
            local_sol = cp.linalg.solve(sigma_ss, rhs)
        except Exception:
            invalid_local_solves += 1
            continue
        if not bool(cp.all(cp.isfinite(local_sol)).item()):
            invalid_local_solves += 1
            continue

        w_i = local_sol[:, 0]
        c_i = local_sol[:, 1:]

        residual_vec = c_i[local_i, :].copy()
        residual_vec[s_idx] = 0.0
        max_abs_residual = cp.maximum(max_abs_residual, cp.max(cp.abs(residual_vec)))

        x_s = x[:, s_idx]
        xw = x_s @ w_i
        score_mat = x - (x_s @ c_i)
        score_mat *= xw[:, None]
        score_mat[:, s_idx] = 0.0
        boot_block_max = cp.maximum(boot_block_max, rolling_block_max_abs(score_mat, l))
        del sigma_ss, rhs_left, rhs, local_sol, w_i, c_i, residual_vec, x_s, xw, score_mat

    q_hat = empirical_quantile_type1(boot_block_max / math.sqrt(l), 1.0 - params.alpha)
    t_hat = math.sqrt(n) * float(max_abs_residual.item())
    return {
        "T_hat": t_hat,
        "q_hat": q_hat,
        "reject": bool(t_hat > q_hat),
        "invalid_local_solves": invalid_local_solves,
    }


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_json(path: Path, obj: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(obj, f, indent=2)


def load_extra_edges_csv(path: str) -> np.ndarray:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
    if data.size == 0:
        return np.empty((0, 3), dtype=np.int64)
    data = np.atleast_1d(data)
    out = np.empty((data.size, 3), dtype=np.int64)
    out[:, 0] = np.asarray(data["i"], dtype=np.int64) - 1
    out[:, 1] = np.asarray(data["j"], dtype=np.int64) - 1
    out[:, 2] = np.asarray(data["sign"], dtype=np.int64)
    return out


def edge_file_for_combo(extra_edges_dir: str, m_extra: int, theta_extra: float) -> str:
    return str(Path(extra_edges_dir) / f"edges_m{m_extra}_theta{theta_extra:g}.csv")


def write_ordered_edges_csv(path: Path, edges: np.ndarray) -> None:
    rows = [
        {"rank": idx + 1, "i": int(i + 1), "j": int(j + 1), "sign": int(sign)}
        for idx, (i, j, sign) in enumerate(edges)
    ]
    write_csv(path, rows, ["rank", "i", "j", "sign"])


def run_one_combo(
    params: RunParams,
    run_id: str,
    combo_index: int,
    total_combos: int,
    extra_edges_override: np.ndarray | None = None,
) -> dict[str, object]:
    out_dir = Path(params.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np_rng = np.random.default_rng(params.seed)
    cp.random.seed(params.seed)

    if extra_edges_override is None:
        extra_edges = sample_extra_edges(params.p, params.b, params.m_extra, np_rng)
    else:
        extra_edges = extra_edges_override.copy()
        if extra_edges.shape[0] != params.m_extra:
            raise ValueError(
                f"extra_edges_csv has {extra_edges.shape[0]} edges but m_extra={params.m_extra}"
            )
    omega, extra_abs, precision_diagnostics = build_precision_matrix(params, extra_edges)
    omega_gpu = cp.asarray(omega)
    r_upper = cp.linalg.cholesky(omega_gpu).T
    b_fft = prepare_fft_filter(params)

    epoch_rows: list[dict[str, object]] = []
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
            test_out = tgg_test_one_epoch_gpu(x_batch[row_left:row_right, :], params)
            epoch_rows.append({"epoch": global_epoch, **test_out})
            if global_epoch == params.epochs or global_epoch % max(1, params.batch_epochs) == 0:
                elapsed = time.perf_counter() - combo_start
                print(
                    f"  combo {combo_index}/{total_combos} "
                    f"{combo_label(params)}: "
                    f"epoch {global_epoch}/{params.epochs}, elapsed {format_seconds(elapsed)}",
                    flush=True,
                )

        del x_batch
        cp.get_default_memory_pool().free_all_blocks()
        epoch_start = epoch_end + 1

    cp.cuda.Stream.null.synchronize()
    runtime_sec = time.perf_counter() - combo_start
    reject_values = [bool(row["reject"]) for row in epoch_rows]
    rejection_rate = float(np.mean(reject_values))
    mean_t_hat = float(np.mean([float(row["T_hat"]) for row in epoch_rows]))
    mean_q_hat = float(np.mean([float(row["q_hat"]) for row in epoch_rows]))
    invalid_total = int(np.sum([int(row["invalid_local_solves"]) for row in epoch_rows]))

    stem = combo_stem(params, run_id)
    epoch_csv = out_dir / f"{stem}_epochs.csv"
    result_json = out_dir / f"{stem}.json"
    write_csv(
        epoch_csv,
        epoch_rows,
        ["epoch", "T_hat", "q_hat", "reject", "invalid_local_solves"],
    )

    sampled_extra_edges = [
        {
            "i": int(i + 1),
            "j": int(j + 1),
            "sign": int(sign),
            "weight": float(extra_abs[pos] * sign),
            "abs_weight": float(extra_abs[pos]),
        }
        for pos, (i, j, sign) in enumerate(extra_edges)
    ]
    result_obj = {
        "run_id": run_id,
        "backend": "cupy",
        "params": asdict(params),
        "package_versions": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "cupy": cp.__version__,
            "cuda_device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        },
        "sampled_extra_edges": sampled_extra_edges,
        "edge_diagnostics": precision_diagnostics,
        "precision_diagnostics": precision_diagnostics,
        "epoch_results_csv": str(epoch_csv),
        "summary": {
            "rejection_rate": rejection_rate,
            "runtime_sec": runtime_sec,
            "mean_T_hat": mean_t_hat,
            "mean_q_hat": mean_q_hat,
            "total_invalid_local_solves": invalid_total,
        },
    }
    save_json(result_json, result_obj)

    return {
        "combo_index": combo_index,
        "m_extra": params.m_extra,
        "signal_mode": params.signal_mode,
        "theta_extra": params.theta_extra,
        "kappa": params.kappa,
        "h": params.h,
        "rejection_rate": rejection_rate,
        "runtime_sec": runtime_sec,
        "mean_T_hat": mean_t_hat,
        "mean_q_hat": mean_q_hat,
        "total_invalid_local_solves": invalid_total,
        "status": "ok",
        "error_message": "",
        "result_json": str(result_json),
        "epoch_csv": str(epoch_csv),
    }


def run_grid(args: argparse.Namespace) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.run_label:
        run_id = f"{run_id}_{args.run_label}"
    elif args.combo_stride > 1:
        run_id = f"{run_id}_offset{args.combo_offset}_stride{args.combo_stride}"
    m_grid = parse_int_grid(args.m_extra_grid)
    if args.signal_mode == "kappa":
        kappa_grid = parse_float_grid(args.kappa_grid)
        theta_grid: list[float] = []
        all_combos = [(m, float("nan"), kappa) for kappa in kappa_grid for m in m_grid]
    else:
        theta_grid = parse_float_grid(args.theta_extra_grid)
        kappa_grid = []
        all_combos = [(m, theta, float("nan")) for theta in theta_grid for m in m_grid]
    if args.combo_limit is not None:
        all_combos = all_combos[: args.combo_limit]
    if args.combo_stride < 1:
        raise ValueError("--combo_stride must be at least 1")
    if args.combo_offset < 0 or args.combo_offset >= args.combo_stride:
        raise ValueError("--combo_offset must be in [0, combo_stride)")
    combos = [
        (global_idx, m, theta, kappa)
        for global_idx, (m, theta, kappa) in enumerate(all_combos, start=1)
        if ((global_idx - 1) % args.combo_stride) == args.combo_offset
    ]
    if not combos:
        raise ValueError("No grid cells assigned to this worker; check --combo_offset/--combo_stride")
    extra_edges_override = None
    if args.extra_edges_csv:
        if len(combos) != 1:
            raise ValueError("--extra_edges_csv is only supported for a single grid cell")
        extra_edges_override = load_extra_edges_csv(args.extra_edges_csv)
    if args.extra_edges_csv and args.extra_edges_dir:
        raise ValueError("Use only one of --extra_edges_csv or --extra_edges_dir")
    if args.signal_mode == "kappa" and args.extra_edges_dir:
        raise ValueError("--extra_edges_dir is only supported for theta mode; kappa mode generates one nested edge list")

    ordered_nested_edges = None
    ordered_edges_csv = ""
    if args.signal_mode == "kappa" and extra_edges_override is None:
        b = (args.s0 - 1) // 2
        m_max = max((m for _, m, _, _ in combos), default=0)
        ordered_nested_edges = sample_ordered_nested_edges(
            p=args.p,
            b=b,
            m_max=m_max,
            edge_seed=args.edge_seed,
            sign_seed=args.sign_seed,
            edge_sampling=args.edge_sampling,
            degree_cap=args.degree_cap,
        )
        ordered_edges_path = out_dir / f"nested_edges_p{args.p}_b{b}_m{m_max}_edge{args.edge_seed}_sign{args.sign_seed}_{run_id}.csv"
        write_ordered_edges_csv(ordered_edges_path, ordered_nested_edges)
        ordered_edges_csv = str(ordered_edges_path)

    device_name = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    print(f"Starting GPU TGG grid {run_id} on {device_name}")
    print(
        f"Base config: n={args.n}, p={args.p}, beta={args.beta:g}, s0={args.s0}, "
        f"epochs={args.epochs}, batch_epochs={args.batch_epochs}, alpha={args.alpha:g}"
    )
    print(
        f"Worker assignment: offset={args.combo_offset}, stride={args.combo_stride}, "
        f"assigned_cells={len(combos)}/{len(all_combos)}"
    )
    if args.signal_mode == "kappa":
        print(f"Grid: m={m_grid}, kappa={kappa_grid}")
        if ordered_edges_csv:
            print(f"Nested edge list: {ordered_edges_csv}")
    else:
        print(f"Grid: m={m_grid}, theta={theta_grid}")

    grid_start = time.perf_counter()
    rows: list[dict[str, object]] = []
    for worker_idx, (global_idx, m_extra, theta_extra, kappa) in enumerate(combos, start=1):
        params = make_params(args, m_extra=m_extra, theta_extra=theta_extra, kappa=kappa)
        try:
            print(
                f"\n[worker {worker_idx}/{len(combos)} | global {global_idx}/{len(all_combos)}] "
                f"Running {combo_label(params)}",
                flush=True,
            )
            combo_edges_override = extra_edges_override
            if args.extra_edges_dir:
                combo_edges_override = load_extra_edges_csv(
                    edge_file_for_combo(args.extra_edges_dir, m_extra, theta_extra)
                )
            elif ordered_nested_edges is not None:
                combo_edges_override = ordered_nested_edges[:m_extra, :]
            row = run_one_combo(
                params,
                run_id,
                global_idx,
                len(all_combos),
                extra_edges_override=combo_edges_override,
            )
        except Exception as exc:
            row = {
                "combo_index": global_idx,
                "m_extra": m_extra,
                "signal_mode": args.signal_mode,
                "theta_extra": theta_extra,
                "kappa": kappa,
                "h": float("nan"),
                "rejection_rate": float("nan"),
                "runtime_sec": float("nan"),
                "mean_T_hat": float("nan"),
                "mean_q_hat": float("nan"),
                "total_invalid_local_solves": float("nan"),
                "status": "error",
                "error_message": repr(exc),
                "result_json": "",
                "epoch_csv": "",
            }
            print(f"  ERROR: {exc!r}", flush=True)
        rows.append(row)
        elapsed = time.perf_counter() - grid_start
        eta = elapsed / worker_idx * (len(combos) - worker_idx)
        print(
            f"  Finished worker cell {worker_idx}/{len(combos)} "
            f"(global {global_idx}/{len(all_combos)}). Total elapsed {format_seconds(elapsed)}. "
            f"ETA {format_seconds(eta)}.",
            flush=True,
        )

    total_runtime_sec = time.perf_counter() - grid_start
    grid_csv = out_dir / f"gpu_grid_summary_{run_id}.csv"
    grid_json = out_dir / f"gpu_grid_summary_{run_id}.json"
    grid_txt = out_dir / f"gpu_grid_summary_{run_id}.txt"
    fieldnames = [
        "combo_index",
        "m_extra",
        "signal_mode",
        "theta_extra",
        "kappa",
        "h",
        "rejection_rate",
        "runtime_sec",
        "mean_T_hat",
        "mean_q_hat",
        "total_invalid_local_solves",
        "status",
        "error_message",
        "result_json",
        "epoch_csv",
    ]
    write_csv(grid_csv, rows, fieldnames)
    grid_obj = {
        "run_id": run_id,
        "backend": "cupy",
        "total_runtime_sec": total_runtime_sec,
        "total_runtime_hms": format_seconds(total_runtime_sec),
        "base_args": vars(args),
        "ordered_edges_csv": ordered_edges_csv,
        "package_versions": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "cupy": cp.__version__,
            "cuda_device": device_name,
        },
        "grid_summary_csv": str(grid_csv),
        "grid_summary": rows,
    }
    save_json(grid_json, grid_obj)
    with grid_txt.open("w") as f:
        f.write(f"GPU TGG grid run {run_id}\n")
        f.write(f"backend=cupy\n")
        f.write(f"cuda_device={device_name}\n")
        f.write(f"total_runtime_sec={total_runtime_sec:.3f}\n")
        f.write(f"total_runtime_hms={format_seconds(total_runtime_sec)}\n")
        f.write(f"grid_summary_csv={grid_csv}\n")
        f.write(f"grid_summary_json={grid_json}\n")

    print("\nGPU grid complete.")
    print(f"Total runtime: {format_seconds(total_runtime_sec)} ({total_runtime_sec:.3f} sec)")
    print(f"Saved grid CSV: {grid_csv}")
    print(f"Saved grid JSON: {grid_json}")
    print(f"Saved grid TXT: {grid_txt}")
    return grid_obj


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="GPU/CuPy TGG precision FFT simulation grid.")
    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--p", type=int, default=200)
    parser.add_argument("--beta", type=float, default=0.9)
    parser.add_argument("--s0", type=int, default=10)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--alpha", type=float, default=0.10)
    parser.add_argument("--seed", type=int, default=233)
    parser.add_argument("--batch_epochs", type=int, default=10)
    parser.add_argument("--block_length", type=int, default=None)
    parser.add_argument("--fft_len", type=int, default=None)
    parser.add_argument("--theta0", type=float, default=THETA0_DEFAULT)
    parser.add_argument("--rho", type=float, default=RHO_DEFAULT)
    parser.add_argument("--delta0", type=float, default=DELTA0_DEFAULT)
    parser.add_argument("--m_extra_grid", type=str, default="10,25,50,100,250,500")
    parser.add_argument("--theta_extra_grid", type=str, default="0.1,0.2,0.3,0.4")
    parser.add_argument("--kappa_grid", type=str, default="0.03,0.06,0.09,0.12,0.15")
    parser.add_argument("--gpu_chunk_cols", type=int, default=8)
    parser.add_argument("--diag_mode", choices=["strict", "r_legacy"], default="strict")
    parser.add_argument("--signal_mode", choices=["theta", "kappa"], default="theta")
    parser.add_argument("--edge_seed", type=int, default=202)
    parser.add_argument("--sign_seed", type=int, default=303)
    parser.add_argument("--edge_sampling", choices=["uniform_nested", "degree_capped"], default="uniform_nested")
    parser.add_argument("--degree_cap", type=int, default=None)
    parser.add_argument("--kappa_tol", type=float, default=1e-10)
    parser.add_argument("--kappa_max_iter", type=int, default=10000)
    parser.add_argument("--kappa_damping", type=float, default=0.5)
    project_root = Path(__file__).resolve().parent.parent
    parser.add_argument("--out_dir", type=str, default=str(project_root / "Results" / "raw" / "gpu_raw"))
    parser.add_argument("--combo_limit", type=int, default=None)
    parser.add_argument("--combo_offset", type=int, default=0)
    parser.add_argument("--combo_stride", type=int, default=1)
    parser.add_argument("--run_label", type=str, default="")
    parser.add_argument("--extra_edges_csv", type=str, default=None)
    parser.add_argument("--extra_edges_dir", type=str, default=None)
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    run_grid(args)


if __name__ == "__main__":
    main()
