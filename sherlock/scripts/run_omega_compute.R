#!/usr/bin/env Rscript
# run_omega_compute.R
#
# Stage 1 of the |Omega| > nu thresholding pipeline.
#
# For each subject:
#   - Option A ROI handling (drop zero-recall-variance ROIs from both phases;
#     z-standardize movie phase on the active set).
#   - Compute Sigma_movie = crossprod(X_movie) / n_movie and
#     Omega_movie = solve(Sigma_movie). With p < n_movie no ridge is needed
#     for the Sherlock data; a ridge fallback is included for safety.
#   - Build a 19-point + boundary nu grid by inverting target densities
#     {0.80, 0.70, ..., 0.001} from the empirical |Omega| upper-triangle
#     distribution.
#   - For each candidate, compute degree quantile statistics and run the
#     in-sample movie self-diagnostic (block-bootstrap containment test on
#     the movie data itself).
#
# Outputs (under Results/):
#   omega_matrices/subject_XX_omega.rds   one per subject; contains
#                                         Omega, active_rois, dropped_rois,
#                                         ridge_used, p_alive, n_movie
#   candidate_path.csv                    every (subject, nu) pair with degree
#                                         stats and movie self-diag p-values;
#                                         used by run_omega_bisect.R and
#                                         run_omega_alpha_grid.R to seed the
#                                         bisection bracket.
#
# Run from the sherlock folder on percyland:
#   N_WORKERS=17 Rscript scripts/run_omega_compute.R

Sys.setenv(
  OMP_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  MKL_NUM_THREADS        = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS    = "1"
)

suppressPackageStartupMessages({
  library(Matrix)
  library(parallel)
  library(data.table)
  library(future)
  library(future.apply)
})

script_path <- {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
  } else {
    NA_character_
  }
}
PROJECT_ROOT <- if (is.na(script_path)) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

cfg <- list(
  alpha_test     = 0.05,
  target_density_grid = c(
    0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.25, 0.20,
    0.15, 0.10, 0.075, 0.05, 0.035, 0.025, 0.015, 0.01,
    0.005, 0.0025, 0.001
  ),
  min_offsupport = 500L,
  ridge_test     = 1e-6,
  ridge_invert   = 1e-4,
  block_length   = NULL,
  data_dir       = file.path(PROJECT_ROOT, "processed_fMRI_data"),
  out_dir        = file.path(PROJECT_ROOT, "Results")
)
set.seed(1)

n_workers <- {
  v <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", unset = NA)))
  if (is.na(v) || v < 1L) v <- min(8L, max(1L, parallel::detectCores() - 1L), 17L)
  max(1L, min(v, 17L))
}

dir.create(file.path(cfg$out_dir, "raw", "omega_matrices"),
           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "logs"),   showWarnings = FALSE, recursive = TRUE)

main_log <- file.path(cfg$out_dir, "logs", "compute.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = main_log, append = TRUE)
}
log_msg("=== run_omega_compute.R ===  n_workers=", n_workers,
        "  density_grid_len=", length(cfg$target_density_grid))

# ---------------------------------------------------------------------------
# Loaders / helpers
# ---------------------------------------------------------------------------
orient_matrix <- function(X, p = 268L) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (ncol(X) == p && nrow(X) != p) return(X)
  if (nrow(X) == p && ncol(X) != p) return(t(X))
  stop("ambiguous orientation; got ", nrow(X), " x ", ncol(X))
}
standardize_ts <- function(X, eps = 1e-12) {
  mu <- colMeans(X)
  Xc <- sweep(X, 2L, mu, FUN = "-")
  sds <- sqrt(colSums(Xc * Xc) / (nrow(Xc) - 1L))
  if (any(sds < eps)) stop("zero variance in standardize_ts")
  sweep(Xc, 2L, sds, FUN = "/")
}
parse_subject_id <- function(fname) {
  m <- regmatches(fname, regexpr("subject[0-9]+", fname))
  if (length(m) == 0) return(NA_integer_)
  as.integer(sub("subject", "", m))
}
load_movie_data <- function(path, p = 268L, n = 17L, T_expected = 1976L) {
  e <- new.env(); load(path, envir = e)
  for (nm in ls(e)) {
    obj <- get(nm, envir = e)
    if (is.array(obj) && length(dim(obj)) == 3L &&
        all(c(p, T_expected, n) %in% dim(obj))) {
      d <- dim(obj); pos_n <- which(d == n)[1]; pos_p <- which(d == p)[1]
      out <- vector("list", n)
      perm <- function(s_idx) {
        idx <- list(quote(expr=), quote(expr=), quote(expr=))
        idx[[pos_n]] <- s_idx
        do.call(`[`, c(list(obj), idx, list(drop = TRUE)))
      }
      for (s in seq_len(n)) {
        slab <- perm(s)
        out[[s]] <- if (nrow(slab) == p) t(slab) else slab
      }
      log_msg("Movie '", nm, "' detected_dim=", paste(d, collapse=" x "),
              " axes: subj@", pos_n, " ROI@", pos_p)
      return(out)
    }
  }
  stop("Could not locate movie array in ", path)
}
load_recall_data <- function(folder, p = 268L) {
  files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
  files <- files[grepl("recall", basename(files), ignore.case = TRUE)]
  if (length(files) != 17L) stop("Expected 17 recall files; found ", length(files))
  ids <- vapply(files, parse_subject_id, integer(1L))
  ord <- order(ids); files <- files[ord]; ids <- ids[ord]
  out <- vector("list", length(files))
  meta <- data.table(subject_id = ids, recall_file = basename(files),
                     n_recall = NA_integer_)
  for (k in seq_along(files)) {
    raw <- as.matrix(data.table::fread(files[k], header = TRUE))
    storage.mode(raw) <- "double"
    Xt <- orient_matrix(raw, p = p)
    out[[k]] <- Xt
    meta$n_recall[k] <- nrow(Xt)
  }
  list(data = out, meta = meta)
}

# Block-bootstrap graph containment test (deterministic single-pass).
graph_test_block_bootstrap <- function(X, E0,
                                       alpha_test  = 0.05,
                                       block_length = NULL,
                                       ridge_test  = 1e-6) {
  n <- nrow(X); p <- ncol(X)
  if (is.null(block_length)) {
    bl <- floor(n^(2/3))
    if (bl < 5L) bl <- max(5L, floor(sqrt(n)))
    if (bl >= n) bl <- floor(n / 2L)
    block_length <- as.integer(bl)
  }
  Sigma_hat <- crossprod(X) / n
  E0L <- (as.matrix(E0) != 0); diag(E0L) <- FALSE
  n_blocks <- n - block_length + 1L
  block_stats <- numeric(n_blocks); inv_sqrt_l <- 1 / sqrt(block_length)
  T_hat <- 0
  for (i in seq_len(p)) {
    nbrs <- which(E0L[, i])
    S_i <- sort(unique(c(i, nbrs)))
    off_i <- setdiff(seq_len(p), S_i)
    if (length(off_i) == 0L) next
    A_i <- Sigma_hat[S_i, S_i, drop = FALSE]
    A_i_inv <- tryCatch(chol2inv(chol(A_i)), error = function(e)
                        chol2inv(chol(A_i + ridge_test * diag(nrow(A_i)))))
    pos_i <- match(i, S_i)
    w_Si <- A_i_inv[, pos_i, drop = TRUE]
    w_full <- numeric(p); w_full[S_i] <- w_Si
    r <- as.numeric(Sigma_hat[off_i, , drop = FALSE] %*% w_full)
    T_hat <- max(T_hat, max(abs(r)))
    C_i <- A_i_inv %*% Sigma_hat[S_i, off_i, drop = FALSE]
    a_i <- as.numeric(X[, S_i, drop = FALSE] %*% w_Si)
    Uproj <- X[, off_i, drop = FALSE] - X[, S_i, drop = FALSE] %*% C_i
    Sc <- Uproj * a_i
    cs <- rbind(0, apply(Sc, 2L, cumsum))
    bs <- cs[(block_length + 1L):(n + 1L), , drop = FALSE] -
          cs[1L:(n - block_length + 1L), , drop = FALSE]
    block_stats <- pmax(block_stats, inv_sqrt_l * apply(abs(bs), 1L, max))
  }
  T_hat <- sqrt(n) * T_hat
  q  <- as.numeric(quantile(block_stats, probs = 1 - alpha_test, type = 1L))
  pv <- (1 + sum(block_stats >= T_hat)) / (length(block_stats) + 1L)
  list(T_hat = T_hat, q = q, pvalue = pv,
       reject = T_hat > q, block_length = block_length)
}

deg_stats <- function(A) {
  diag(A) <- FALSE
  d <- rowSums(A)
  e <- sum(A) %/% 2L
  qs <- as.numeric(quantile(d, c(0.10, 0.25, 0.50, 0.75, 0.90, 0.99), type = 7))
  list(edge_count = e, density = e / max(1L, choose(nrow(A), 2L)),
       mean_degree = mean(d),
       degree_q10 = qs[1], degree_q25 = qs[2], degree_q50 = qs[3],
       degree_q75 = qs[4], degree_q90 = qs[5], degree_q99 = qs[6],
       max_degree = max(d), s0 = max(d) + 1L,
       offsupport_pair_count = (sum(!A) - nrow(A)) %/% 2L)
}

# ---------------------------------------------------------------------------
# Per-subject worker: Option A -> Omega -> candidate path with self-diag
# ---------------------------------------------------------------------------
process_subject <- function(input) {
  sid <- input$subject_id
  Xm_raw <- input$Xm; Xr_raw <- input$Xr
  log_path <- file.path(cfg$out_dir, "logs", sprintf("compute_subject_%02d.log", sid))
  wmsg <- function(...) {
    cat("[", format(Sys.time(), "%H:%M:%S"), " pid=", Sys.getpid(), "] ",
        paste0(..., collapse = ""), "\n",
        sep = "", file = log_path, append = TRUE)
  }

  rec_var <- apply(Xr_raw, 2L, var)
  dead <- which(rec_var == 0)
  active <- setdiff(seq_len(ncol(Xr_raw)), dead)
  Xm <- standardize_ts(Xm_raw[, active, drop = FALSE])
  pA <- ncol(Xm); n_movie <- nrow(Xm)
  wmsg("=== subject ", sid, " === pA=", pA, "  n_movie=", n_movie,
       "  dead=", length(dead))

  Sigma_movie <- crossprod(Xm) / n_movie
  Omega <- tryCatch(chol2inv(chol(Sigma_movie)), error = function(e) NULL)
  ridge_used <- 0
  if (is.null(Omega) || any(diag(Omega) <= 0)) {
    Omega <- chol2inv(chol(Sigma_movie + cfg$ridge_invert * diag(pA)))
    ridge_used <- cfg$ridge_invert
    wmsg("ridge ", ridge_used, " used in inversion")
  }

  saveRDS(list(subject_id = sid, p_alive = pA, n_movie = n_movie,
               Omega = Omega, active_rois = active, dropped_rois = dead,
               ridge_used = ridge_used),
          file.path(cfg$out_dir, "raw", "omega_matrices",
                    sprintf("subject_%02d_omega.rds", sid)))

  ut <- which(upper.tri(Omega))
  absOm_ut <- abs(Omega)[ut]
  tg <- cfg$target_density_grid
  nus <- as.numeric(quantile(absOm_ut, probs = 1 - tg,
                             type = 1L, names = FALSE))
  nus <- sort(unique(c(nus, 0)), decreasing = TRUE)  # sparsest first

  rows <- vector("list", length(nus))
  for (k in seq_along(nus)) {
    nu <- nus[k]
    A <- abs(Omega) > nu; diag(A) <- FALSE
    ds <- deg_stats(A)
    if (ds$offsupport_pair_count == 0L) {
      mt <- list(T_hat = NA_real_, q = NA_real_, pvalue = NA_real_)
    } else {
      mt <- tryCatch(
        graph_test_block_bootstrap(Xm, A,
                                   alpha_test  = cfg$alpha_test,
                                   block_length = cfg$block_length,
                                   ridge_test  = cfg$ridge_test),
        error = function(e) {
          wmsg("self-diag error nu=", signif(nu, 3), ": ", conditionMessage(e))
          list(T_hat = NA_real_, q = NA_real_, pvalue = NA_real_)
        })
    }
    rows[[k]] <- data.table(
      subject_id = sid, p_alive = pA, cand_idx = k, nu = nu,
      edge_count = ds$edge_count, density = ds$density,
      mean_degree = ds$mean_degree,
      degree_q10 = ds$degree_q10, degree_q25 = ds$degree_q25,
      degree_q50 = ds$degree_q50, degree_q75 = ds$degree_q75,
      degree_q90 = ds$degree_q90, degree_q99 = ds$degree_q99,
      max_degree = ds$max_degree, s0_estim = ds$s0,
      offsupport_pair_count = ds$offsupport_pair_count,
      movie_T = mt$T_hat, movie_q = mt$q, movie_p_raw = mt$pvalue,
      movie_self_reject_005 = if (is.na(mt$pvalue)) NA else mt$pvalue < 0.05,
      movie_self_reject_010 = if (is.na(mt$pvalue)) NA else mt$pvalue < 0.10)
  }
  wmsg("done subject ", sid, "; ", length(nus), " candidates")
  rbindlist(rows)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
log_msg("Loading data ...")
movie_path <- file.path(cfg$data_dir, "moviefmrip268n17T1976.RData")
movie_list <- load_movie_data(movie_path)
recall <- load_recall_data(cfg$data_dir)
recall_data <- recall$data; recall_meta <- recall$meta

subject_inputs <- lapply(seq_len(17L), function(i) {
  list(subject_id = recall_meta$subject_id[i],
       Xm = movie_list[[recall_meta$subject_id[i]]],
       Xr = recall_data[[i]])
})

log_msg("Setting up future plan: multisession with ", n_workers, " workers")
plan(multisession, workers = n_workers)
on.exit(plan(sequential), add = TRUE)

log_msg("Launching parallel subject loop (chunk.size = 1) ...")
t0 <- Sys.time()
all_rows <- future.apply::future_lapply(
  subject_inputs, process_subject,
  future.seed = TRUE, future.chunk.size = 1L
)
cand_dt <- rbindlist(all_rows)
setorder(cand_dt, subject_id, -nu)
fwrite(cand_dt, file.path(cfg$out_dir, "tables", "candidate_path.csv"))
log_msg("Subject loop done in ",
        sprintf("%.1f", as.numeric(Sys.time() - t0, units = "secs")), "s")
log_msg("Wrote ", file.path(cfg$out_dir, "tables", "candidate_path.csv"))
log_msg("Wrote 17 RDS files in ", file.path(cfg$out_dir, "raw", "omega_matrices/"))
log_msg("Done.")
