#!/usr/bin/env Rscript
# run_omega_alpha_grid.R
#
# For each subject and each alpha in {0.05, 0.025, 0.01, 0.005, 0.001}:
#   bisect on nu to find the LARGEST nu such that the |Omega| > nu movie
#   graph passes the in-sample movie self-diagnostic at level alpha
#   (i.e. movie_p_raw >= alpha) with offsupport_pair_count >= 500.
#   Then run the block-bootstrap recall containment test and Holm-correct
#   across the 17 subjects WITHIN each alpha.
#
# Lower alpha -> more lenient acceptance -> sparser selected graphs.
#
# Brackets chain: alpha=0.05 boundary -> lower bound for alpha=0.025; etc.
#
# Run from the sherlock folder:
#   N_WORKERS=17 Rscript scripts/run_omega_alpha_grid.R

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
  alphas         = c(0.05, 0.025, 0.01, 0.005, 0.001),
  min_offsupport = 500L,
  ridge_test     = 1e-6,
  bisect_tol     = 1e-5,
  bisect_maxit   = 30L,
  block_length   = NULL,
  data_dir       = file.path(PROJECT_ROOT, "processed_fMRI_data"),
  base_dir       = file.path(PROJECT_ROOT, "Results"),
  out_dir        = file.path(PROJECT_ROOT, "Results")
)
set.seed(1)

n_workers <- {
  v <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", unset = NA)))
  if (is.na(v) || v < 1L) v <- min(8L, max(1L, parallel::detectCores() - 1L), 17L)
  max(1L, min(v, 17L))
}
dir.create(file.path(cfg$out_dir, "raw", "selected_graphs_alpha_grid"),
           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "logs"),   showWarnings = FALSE, recursive = TRUE)
main_log <- file.path(cfg$out_dir, "logs", "main_alpha_grid.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = main_log, append = TRUE)
}
log_msg("=== run_omega_alpha_grid.R === alphas=", paste(cfg$alphas, collapse=","))

# ---------------------------------------------------------------------------
# Helpers (lifted from prior scripts)
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
      d <- dim(obj); pos_n <- which(d == n)[1]
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
      return(out)
    }
  }
  stop("Could not locate movie array in ", path)
}
load_recall_data <- function(folder, p = 268L) {
  files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
  files <- files[grepl("recall", basename(files), ignore.case = TRUE)]
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
graph_test_block_bootstrap <- function(X, E0,
                                       alpha_test = 0.05,
                                       block_length = NULL,
                                       ridge_test = 1e-6) {
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
log_msg("Loading data + previous candidates ...")
movie_path <- file.path(cfg$data_dir, "moviefmrip268n17T1976.RData")
movie_list <- load_movie_data(movie_path)
recall <- load_recall_data(cfg$data_dir)
recall_data <- recall$data; recall_meta <- recall$meta
cand_all <- fread(file.path(cfg$base_dir, "tables", "candidate_path.csv"))

# Initial bracket from the 19-pt grid (sparsest passing at alpha=0.05, next-rejecting)
bracket_for <- function(sid) {
  sub <- cand_all[subject_id == sid &
                  offsupport_pair_count >= cfg$min_offsupport]
  pass <- sub[movie_self_reject_005 == FALSE]
  rej  <- sub[movie_self_reject_005 == TRUE]
  if (nrow(pass) == 0L) stop("subj ", sid, ": no passing candidate")
  nu_lo <- max(pass$nu)
  rej_above <- rej[nu > nu_lo]
  nu_hi <- if (nrow(rej_above) == 0L) nu_lo * 1.5 else min(rej_above$nu)
  list(nu_lo = nu_lo, nu_hi = nu_hi)
}

subject_inputs <- lapply(seq_len(17L), function(i) {
  sid <- recall_meta$subject_id[i]
  list(subject_id = sid,
       Xm = movie_list[[sid]],
       Xr = recall_data[[i]],
       sel_path = file.path(cfg$base_dir, "raw", "omega_matrices",
                            sprintf("subject_%02d_omega.rds", sid)),
       bracket = bracket_for(sid))
})

# ---------------------------------------------------------------------------
# Per-subject worker: bisect for each alpha sequentially, chaining brackets
# ---------------------------------------------------------------------------
process_subject <- function(input) {
  sid <- input$subject_id
  prev <- readRDS(input$sel_path)
  Omega <- prev$Omega
  active <- prev$active_rois
  Xm <- standardize_ts(input$Xm[, active, drop = FALSE])
  Xr <- standardize_ts(input$Xr[, active, drop = FALSE])
  pA <- ncol(Omega)
  log_path <- file.path(cfg$out_dir, "logs",
                        sprintf("subject_%02d_alpha_grid.log", sid))
  wmsg <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%H:%M:%S"),
                  " pid=", Sys.getpid(), "] ", paste0(...))
    cat(msg, "\n", sep = "", file = log_path, append = TRUE)
  }

  test_at <- function(nu) {
    A <- abs(Omega) > nu; diag(A) <- FALSE
    ds <- deg_stats(A)
    if (ds$offsupport_pair_count < cfg$min_offsupport)
      return(list(p_raw = NA_real_, A = A, ds = ds,
                  T_hat = NA_real_, q = NA_real_,
                  note = "offsupport_below_min"))
    rt <- graph_test_block_bootstrap(Xm, A, alpha_test = 0.05,
                                     block_length = cfg$block_length,
                                     ridge_test  = cfg$ridge_test)
    list(p_raw = rt$pvalue, A = A, ds = ds,
         T_hat = rt$T_hat, q = rt$q, note = "")
  }

  bisect_for_alpha <- function(alpha, nu_lo_init, nu_hi_init) {
    nu_lo <- nu_lo_init; nu_hi <- nu_hi_init
    res_lo <- test_at(nu_lo); res_hi <- test_at(nu_hi)

    # Walk down nu_lo to find passing if needed
    walks_down <- 0L
    while (!is.na(res_lo$p_raw) && res_lo$p_raw < alpha &&
           res_lo$ds$offsupport_pair_count >= cfg$min_offsupport &&
           walks_down < 8L) {
      nu_lo <- nu_lo * 0.85
      res_lo <- test_at(nu_lo)
      walks_down <- walks_down + 1L
    }
    # Walk up nu_hi to find rejecting if needed
    walks_up <- 0L
    while (!is.na(res_hi$p_raw) && res_hi$p_raw >= alpha && walks_up < 8L) {
      nu_hi <- nu_hi * 1.2
      res_hi <- test_at(nu_hi)
      walks_up <- walks_up + 1L
    }
    if (is.na(res_lo$p_raw) || res_lo$p_raw < alpha) {
      return(list(success = FALSE, sel = res_lo,
                  nu = nu_lo, alpha = alpha, iters = 0L,
                  note = "no_passing_graph_found"))
    }
    if (is.na(res_hi$p_raw) || res_hi$p_raw >= alpha) {
      return(list(success = TRUE, sel = res_hi,
                  nu = nu_hi, alpha = alpha, iters = 0L,
                  note = "even_widest_nu_passes"))
    }

    best_pass <- res_lo
    iter <- 0L
    while ((nu_hi - nu_lo) > cfg$bisect_tol && iter < cfg$bisect_maxit) {
      iter <- iter + 1L
      nu_mid <- 0.5 * (nu_lo + nu_hi)
      r <- test_at(nu_mid)
      if (!is.na(r$p_raw) && r$p_raw >= alpha &&
          r$ds$offsupport_pair_count >= cfg$min_offsupport) {
        nu_lo <- nu_mid; best_pass <- r
      } else {
        nu_hi <- nu_mid
      }
    }
    list(success = TRUE, sel = best_pass, nu = nu_lo,
         alpha = alpha, iters = iter, note = "")
  }

  rows <- vector("list", length(cfg$alphas))
  prev_nu_lo <- input$bracket$nu_lo  # alpha=0.05 grid lower bound
  for (a_idx in seq_along(cfg$alphas)) {
    alpha <- cfg$alphas[a_idx]
    wmsg("--- alpha=", alpha, " ---")
    # For lower alpha (more lenient), boundary moves to LARGER nu.
    # Lower bound = previous alpha's selected nu (same or larger graph passes here).
    # Upper bound = original grid hi (rejected at alpha=0.05, so rejects at any larger alpha;
    # for lower alpha may pass — bisect_for_alpha will widen if needed).
    nu_lo_init <- prev_nu_lo
    nu_hi_init <- input$bracket$nu_hi
    res <- bisect_for_alpha(alpha, nu_lo_init, nu_hi_init)
    wmsg("  nu_sel=", signif(res$nu, 5),
         " p_movie=", signif(res$sel$p_raw, 4),
         " density=", signif(res$sel$ds$density, 4),
         " iters=", res$iters,
         " note=", res$note)
    if (!res$success) {
      rows[[a_idx]] <- data.table(
        subject_id = sid, p_alive = pA, alpha = alpha,
        selected_nu = NA_real_, density = NA_real_, edge_count = NA_integer_,
        mean_degree = NA_real_,
        degree_q10 = NA_real_, degree_q25 = NA_real_, degree_q50 = NA_real_,
        degree_q75 = NA_real_, degree_q90 = NA_real_, degree_q99 = NA_real_,
        max_degree = NA_integer_, s0_estim = NA_integer_,
        offsupport_pair_count = NA_integer_,
        movie_T = NA_real_, movie_q = NA_real_, movie_p_raw = NA_real_,
        movie_self_reject_alpha = NA,
        recall_T = NA_real_, recall_q = NA_real_,
        recall_p_raw = NA_real_, recall_reject_raw = NA,
        bisect_iters = res$iters, notes = res$note)
      next
    }
    # Recall test at this selection
    A_sel <- res$sel$A; nu_sel <- res$nu
    rt <- tryCatch(
      graph_test_block_bootstrap(Xr, A_sel, alpha_test = 0.05,
                                 block_length = cfg$block_length,
                                 ridge_test  = cfg$ridge_test),
      error = function(e) NULL)
    if (is.null(rt)) {
      rT <- NA_real_; rq <- NA_real_; rp <- NA_real_; rr <- NA
    } else {
      rT <- rt$T_hat; rq <- rt$q; rp <- rt$pvalue
      rr <- rp < alpha
    }
    ds <- res$sel$ds
    saveRDS(list(A = A_sel, Omega = Omega, nu = nu_sel,
                 alpha = alpha, density = ds$density,
                 max_degree = ds$max_degree,
                 movie_p_raw = res$sel$p_raw,
                 active_rois = active, dropped_rois = prev$dropped_rois,
                 subject_id = sid, bisect_iters = res$iters),
            file.path(cfg$out_dir, "raw", "selected_graphs_alpha_grid",
                      sprintf("subject_%02d_alpha_%s_graph.rds",
                              sid, gsub("\\.", "p", as.character(alpha)))))
    rows[[a_idx]] <- data.table(
      subject_id = sid, p_alive = pA, alpha = alpha,
      selected_nu = nu_sel,
      edge_count = ds$edge_count, density = ds$density,
      mean_degree = ds$mean_degree,
      degree_q10 = ds$degree_q10, degree_q25 = ds$degree_q25,
      degree_q50 = ds$degree_q50, degree_q75 = ds$degree_q75,
      degree_q90 = ds$degree_q90, degree_q99 = ds$degree_q99,
      max_degree = ds$max_degree, s0_estim = ds$s0,
      offsupport_pair_count = ds$offsupport_pair_count,
      movie_T = res$sel$T_hat, movie_q = res$sel$q,
      movie_p_raw = res$sel$p_raw,
      movie_self_reject_alpha = res$sel$p_raw < alpha,
      recall_T = rT, recall_q = rq,
      recall_p_raw = rp, recall_reject_raw = rr,
      bisect_iters = res$iters, notes = res$note)
    # Chain: this nu becomes the lower bound for the next (lower) alpha
    prev_nu_lo <- nu_sel
  }
  rbindlist(rows, fill = TRUE)
}

log_msg("Setting up future plan: multisession with ", n_workers, " workers")
plan(multisession, workers = n_workers)
on.exit(plan(sequential), add = TRUE)
t0 <- Sys.time()
all_rows <- future.apply::future_lapply(
  subject_inputs, process_subject,
  future.seed = TRUE, future.chunk.size = 1L
)
dt <- rbindlist(all_rows, fill = TRUE)
log_msg("Subject loop done in ",
        sprintf("%.1f", as.numeric(Sys.time() - t0, units = "secs")), "s")

# Holm + Benjamini-Hochberg correction within each alpha
dt[, recall_p_holm    := NA_real_]
dt[, recall_reject_holm := NA]
dt[, recall_p_bh      := NA_real_]
dt[, recall_reject_bh   := NA]
for (a in cfg$alphas) {
  rows_idx <- which(dt$alpha == a & !is.na(dt$recall_p_raw))
  if (length(rows_idx) > 0L) {
    pv <- dt$recall_p_raw[rows_idx]
    ph <- p.adjust(pv, method = "holm", n = 17L)
    pb <- p.adjust(pv, method = "BH",   n = 17L)
    dt$recall_p_holm[rows_idx]      <- ph
    dt$recall_reject_holm[rows_idx] <- ph < a
    dt$recall_p_bh[rows_idx]        <- pb
    dt$recall_reject_bh[rows_idx]   <- pb < a
  }
}

setorder(dt, alpha, subject_id)
fwrite(dt, file.path(cfg$out_dir, "tables", "alpha_grid_subject_results.csv"))

# Per-alpha summary
agg <- dt[, .(
  n_subj                  = .N,
  n_self_pass             = sum(!is.na(selected_nu)),
  n_recall_reject_raw     = sum(recall_reject_raw   == TRUE, na.rm = TRUE),
  n_recall_reject_holm    = sum(recall_reject_holm  == TRUE, na.rm = TRUE),
  n_recall_reject_bh      = sum(recall_reject_bh    == TRUE, na.rm = TRUE),
  median_density          = median(density,    na.rm = TRUE),
  min_density             = min(density,       na.rm = TRUE),
  max_density             = max(density,       na.rm = TRUE),
  median_edge_count       = as.integer(median(edge_count, na.rm = TRUE)),
  median_max_degree       = as.integer(median(max_degree, na.rm = TRUE)),
  median_q50              = median(degree_q50, na.rm = TRUE),
  median_q90              = median(degree_q90, na.rm = TRUE),
  median_q99              = median(degree_q99, na.rm = TRUE),
  median_movie_p_raw      = median(movie_p_raw, na.rm = TRUE)
), by = alpha][order(-alpha)]
fwrite(agg, file.path(cfg$out_dir, "tables", "alpha_grid_summary.csv"))

# ---------------------------------------------------------------------------
sink(file.path(cfg$out_dir, "tables", "summary_alpha_grid.txt"))
cat("============================================================\n")
cat("  |Omega|>nu thresholding, BISECTED nu per subject across an alpha grid\n")
cat("  alphas = ", paste(cfg$alphas, collapse = ", "), "\n", sep = "")
cat("============================================================\n\n")

cat("--- Per-alpha summary (cohort) ---\n")
print(agg)
cat("\n")

cat("--- Headline ---\n")
for (k in seq_len(nrow(agg))) {
  with(agg[k], {
    cat(sprintf(
      "alpha=%6.4f : %2d/17 self-pass | recall %2d/17 raw, %2d/17 Holm, %2d/17 BH | median_density=%.4f (range %.4f-%.4f) | median_max_deg=%3d\n",
      alpha, n_self_pass, n_recall_reject_raw, n_recall_reject_holm, n_recall_reject_bh,
      median_density, min_density, max_density, median_max_degree))
  })
}
cat("\n")

cat("--- Per-subject density across alpha levels ---\n")
wide_dens <- dcast(dt, subject_id ~ alpha,
                   value.var = "density", fun.aggregate = mean)
setnames(wide_dens, old = setdiff(names(wide_dens), "subject_id"),
         new = paste0("dens_a", setdiff(names(wide_dens), "subject_id")))
print(wide_dens)
cat("\n")

cat("--- Per-subject max_degree across alpha levels ---\n")
wide_max <- dcast(dt, subject_id ~ alpha,
                  value.var = "max_degree", fun.aggregate = function(x) as.integer(round(mean(x))))
setnames(wide_max, old = setdiff(names(wide_max), "subject_id"),
         new = paste0("max_a", setdiff(names(wide_max), "subject_id")))
print(wide_max)
cat("\n")

cat("--- Per-subject recall_p_holm across alpha levels ---\n")
wide_holm <- dcast(dt, subject_id ~ alpha,
                   value.var = "recall_p_holm", fun.aggregate = mean)
setnames(wide_holm, old = setdiff(names(wide_holm), "subject_id"),
         new = paste0("holm_a", setdiff(names(wide_holm), "subject_id")))
print(wide_holm)
cat("\n")

cat("============================================================\n")
cat("Subject-level CSV   : ", normalizePath(file.path(cfg$out_dir, "tables", "alpha_grid_subject_results.csv"), winslash="/"), "\n", sep="")
cat("Per-alpha summary   : ", normalizePath(file.path(cfg$out_dir, "tables", "alpha_grid_summary.csv"),         winslash="/"), "\n", sep="")
sink()
log_msg("Done.")
