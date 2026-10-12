#!/usr/bin/env Rscript
# run_omega_bisect.R
#
# Refine the |Omega| > nu primary by bisecting on nu per subject:
# find the largest nu such that the |Omega|>nu movie graph passes the
# in-sample movie self-diagnostic at alpha=0.05 with offsupport_pair_count
# >= 500. Then run the block-bootstrap recall containment test against
# the bisected graph, Holm-correct across 17 subjects.
#
# Reuses Omega matrices saved in omega_matrices/subject_XX_omega.rds and the
# initial bracket from candidate_path.csv (both produced by run_omega_compute.R).
#
# Run from the sherlock folder:
#   N_WORKERS=17 Rscript scripts/run_omega_bisect.R

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
  min_offsupport = 500L,
  ridge_test     = 1e-6,
  bisect_tol     = 1e-5,    # absolute tolerance on nu
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

dir.create(file.path(cfg$out_dir, "raw", "selected_graphs_omega_bisect"),
           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "logs"),   showWarnings = FALSE, recursive = TRUE)
main_log <- file.path(cfg$out_dir, "logs", "main_bisect.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = main_log, append = TRUE)
}
log_msg("=== run_omega_bisect.R ===  tol=", cfg$bisect_tol, "  maxit=", cfg$bisect_maxit)

# ---------------------------------------------------------------------------
# Helpers
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
      d <- dim(obj)
      pos_n <- which(d == n)[1]; pos_p <- which(d == p)[1]
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
# Initial bracket lookup from candidate path
# ---------------------------------------------------------------------------
log_msg("Loading data + previous candidates ...")
movie_path <- file.path(cfg$data_dir, "moviefmrip268n17T1976.RData")
movie_list <- load_movie_data(movie_path)
recall <- load_recall_data(cfg$data_dir)
recall_data <- recall$data; recall_meta <- recall$meta
cand_all <- fread(file.path(cfg$base_dir, "tables", "candidate_path.csv"))

bracket_for <- function(sid) {
  sub <- cand_all[subject_id == sid &
                  offsupport_pair_count >= cfg$min_offsupport]
  pass <- sub[movie_self_reject_005 == FALSE]
  rej  <- sub[movie_self_reject_005 == TRUE]
  if (nrow(pass) == 0L) stop("subj ", sid, ": no passing candidate in bracket")
  nu_lo <- max(pass$nu)                             # largest passing nu
  rej_above <- rej[nu > nu_lo]
  if (nrow(rej_above) == 0L) {
    nu_hi <- nu_lo * 1.5                            # no rejecting above; widen
  } else {
    nu_hi <- min(rej_above$nu)                       # smallest rejecting nu above
  }
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
# Per-subject bisection worker
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
                        sprintf("subject_%02d_bisect.log", sid))
  wmsg <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%H:%M:%S"),
                  " pid=", Sys.getpid(), "] ", paste0(...))
    cat(msg, "\n", sep = "", file = log_path, append = TRUE)
  }

  test_at <- function(nu) {
    A <- abs(Omega) > nu; diag(A) <- FALSE
    ds <- deg_stats(A)
    if (ds$offsupport_pair_count < cfg$min_offsupport)
      return(list(pass = FALSE, p_raw = NA_real_, A = A, ds = ds,
                  T_hat = NA_real_, q = NA_real_,
                  note = "offsupport_below_min"))
    rt <- graph_test_block_bootstrap(Xm, A, alpha_test = cfg$alpha_test,
                                     block_length = cfg$block_length,
                                     ridge_test  = cfg$ridge_test)
    list(pass = (rt$pvalue >= cfg$alpha_test),
         p_raw = rt$pvalue, A = A, ds = ds,
         T_hat = rt$T_hat, q = rt$q, note = "")
  }

  nu_lo <- input$bracket$nu_lo
  nu_hi <- input$bracket$nu_hi
  wmsg("=== subject ", sid, " === pA=", pA,
       " bracket=[", signif(nu_lo, 5), ", ", signif(nu_hi, 5), "]")

  # Verify endpoints
  res_lo <- test_at(nu_lo)
  res_hi <- test_at(nu_hi)
  wmsg(" check lo nu=", signif(nu_lo, 5), " p=", signif(res_lo$p_raw, 4),
       " pass=", res_lo$pass)
  wmsg(" check hi nu=", signif(nu_hi, 5), " p=", signif(res_hi$p_raw, 4),
       " pass=", res_hi$pass)
  if (!res_lo$pass) {
    wmsg(" WARNING: nu_lo does not pass; widening bracket downward")
    # Walk down to find a passing nu_lo (cap to a few tries)
    for (try in seq_len(6L)) {
      nu_lo <- nu_lo * 0.85
      res_lo <- test_at(nu_lo)
      if (res_lo$pass) break
    }
  }
  if (res_hi$pass) {
    wmsg(" WARNING: nu_hi already passes; widening bracket upward")
    for (try in seq_len(6L)) {
      nu_hi <- nu_hi * 1.2
      res_hi <- test_at(nu_hi)
      if (!res_hi$pass) break
    }
  }

  best_pass <- res_lo  # largest known passing
  iters <- data.table(iter = 0L, nu = nu_lo, p_raw = res_lo$p_raw,
                      pass = res_lo$pass, density = res_lo$ds$density,
                      offsupport = res_lo$ds$offsupport_pair_count)
  iters <- rbind(iters, data.table(iter = 0L, nu = nu_hi, p_raw = res_hi$p_raw,
                                   pass = res_hi$pass, density = res_hi$ds$density,
                                   offsupport = res_hi$ds$offsupport_pair_count))

  iter <- 0L
  while ((nu_hi - nu_lo) > cfg$bisect_tol && iter < cfg$bisect_maxit) {
    iter <- iter + 1L
    nu_mid <- 0.5 * (nu_lo + nu_hi)
    r <- test_at(nu_mid)
    iters <- rbind(iters,
                   data.table(iter = iter, nu = nu_mid, p_raw = r$p_raw,
                              pass = r$pass, density = r$ds$density,
                              offsupport = r$ds$offsupport_pair_count))
    if (r$pass) {
      nu_lo <- nu_mid
      best_pass <- r
    } else {
      nu_hi <- nu_mid
    }
    wmsg(sprintf(" iter=%2d nu_mid=%.6f p=%.4f pass=%s | bracket=[%.6f,%.6f]",
                 iter, nu_mid, r$p_raw, r$pass, nu_lo, nu_hi))
  }
  wmsg(" CONVERGED iter=", iter, " final nu_lo=", signif(nu_lo, 6),
       " density=", signif(best_pass$ds$density, 4),
       " p_movie=", signif(best_pass$p_raw, 4))

  # Recall test against the bisected (largest passing) graph
  A_sel <- best_pass$A; nu_sel <- nu_lo
  rt <- tryCatch(
    graph_test_block_bootstrap(Xr, A_sel,
                               alpha_test  = cfg$alpha_test,
                               block_length = cfg$block_length,
                               ridge_test  = cfg$ridge_test),
    error = function(e) {
      wmsg("recall test error: ", conditionMessage(e)); NULL
    })
  if (is.null(rt)) {
    rT <- NA_real_; rq <- NA_real_; rp <- NA_real_; rr <- NA
  } else {
    rT <- rt$T_hat; rq <- rt$q; rp <- rt$pvalue; rr <- rt$reject
  }

  ds <- best_pass$ds
  saveRDS(list(A = A_sel, Omega = Omega, nu = nu_sel,
               density = ds$density, edge_count = ds$edge_count,
               max_degree = ds$max_degree, s0_estim = ds$s0,
               offsupport_pair_count = ds$offsupport_pair_count,
               movie_p_raw = best_pass$p_raw,
               active_rois = active, dropped_rois = prev$dropped_rois,
               subject_id = sid, bisect_iters = iter,
               bracket_initial = c(input$bracket$nu_lo, input$bracket$nu_hi)),
          file.path(cfg$out_dir, "raw", "selected_graphs_omega_bisect",
                    sprintf("subject_%02d_thresholded_omega_bisect_graph.rds", sid)))

  fwrite(iters, file.path(cfg$out_dir, "logs",
                          sprintf("subject_%02d_bisect_path.csv", sid)))

  data.table(
    subject_id = sid, p_alive = pA,
    selected_nu = nu_sel,
    bracket_lo_initial = input$bracket$nu_lo,
    bracket_hi_initial = input$bracket$nu_hi,
    bisect_iters = iter,
    edge_count = ds$edge_count, density = ds$density,
    mean_degree = ds$mean_degree,
    degree_q10 = ds$degree_q10, degree_q25 = ds$degree_q25,
    degree_q50 = ds$degree_q50, degree_q75 = ds$degree_q75,
    degree_q90 = ds$degree_q90, degree_q99 = ds$degree_q99,
    max_degree = ds$max_degree, s0_estim = ds$s0,
    offsupport_pair_count = ds$offsupport_pair_count,
    movie_T = best_pass$T_hat, movie_q = best_pass$q,
    movie_p_raw = best_pass$p_raw,
    movie_self_reject_005 = best_pass$p_raw < cfg$alpha_test,
    movie_self_reject_010 = best_pass$p_raw < 0.10,
    recall_T = rT, recall_q = rq,
    recall_p_raw = rp, recall_reject_raw_005 = rr)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
log_msg("Setting up future plan: multisession with ", n_workers, " workers")
plan(multisession, workers = n_workers)
on.exit(plan(sequential), add = TRUE)
t0 <- Sys.time()
rows <- future.apply::future_lapply(
  subject_inputs, process_subject,
  future.seed = TRUE, future.chunk.size = 1L
)
dt <- rbindlist(rows); setorder(dt, subject_id)
log_msg("Subject loop done in ",
        sprintf("%.1f", as.numeric(Sys.time() - t0, units = "secs")), "s")

dt[, recall_p_holm := p.adjust(recall_p_raw, method = "holm")]
dt[, recall_reject_holm_005 := recall_p_holm < cfg$alpha_test]
dt[, recall_p_bh   := p.adjust(recall_p_raw, method = "BH")]
dt[, recall_reject_bh_005   := recall_p_bh   < cfg$alpha_test]

fwrite(dt, file.path(cfg$out_dir, "tables", "selected_subject_results_omega_bisect.csv"))
fwrite(dt[, .(subject_id, p_alive, selected_nu, bisect_iters,
              density, edge_count, max_degree, s0_estim,
              movie_p_raw, movie_self_reject_005, movie_self_reject_010,
              recall_T, recall_q, recall_p_raw,
              recall_p_holm, recall_reject_holm_005,
              recall_p_bh, recall_reject_bh_005,
              recall_reject_raw_005)],
       file.path(cfg$out_dir, "tables", "recall_test_results_omega_bisect.csv"))

# ---------------------------------------------------------------------------
# Summary file
# ---------------------------------------------------------------------------
sink(file.path(cfg$out_dir, "tables", "summary_omega_bisect.txt"))
cat("============================================================\n")
cat("  |Omega| > nu thresholding, BISECTED nu per subject\n")
cat("============================================================\n\n")
cat("alpha_test     = ", cfg$alpha_test,     "\n", sep = "")
cat("min_offsupport = ", cfg$min_offsupport, "\n", sep = "")
cat("bisect_tol     = ", cfg$bisect_tol,     "\n", sep = "")
cat("bisect_maxit   = ", cfg$bisect_maxit,   "\n", sep = "")
cat("n_subjects     = 17\n\n")

n_self_005     <- sum(dt$movie_self_reject_005 == FALSE, na.rm = TRUE)
n_self_010     <- sum(dt$movie_self_reject_010 == FALSE, na.rm = TRUE)
n_recall_raw   <- sum(dt$recall_reject_raw_005 == TRUE, na.rm = TRUE)
n_recall_holm  <- sum(dt$recall_reject_holm_005 == TRUE, na.rm = TRUE)
n_recall_bh    <- sum(dt$recall_reject_bh_005   == TRUE, na.rm = TRUE)

cat("--- Headline counts ---\n")
cat(sprintf("Subjects with movie_self_reject_005 == FALSE       : %2d / 17\n", n_self_005))
cat(sprintf("Subjects with movie_self_reject_010 == FALSE       : %2d / 17\n", n_self_010))
cat(sprintf("Recall raw rejections (alpha=0.05)                 : %2d / 17\n", n_recall_raw))
cat(sprintf("Recall Holm rejections (alpha=0.05)                : %2d / 17\n", n_recall_holm))
cat(sprintf("Recall BH rejections (alpha=0.05)                  : %2d / 17\n", n_recall_bh))
cat("\n")

cat("--- Cohort medians (bisected selected graph) ---\n")
cat(sprintf("median selected nu         : %.5f\n", median(dt$selected_nu)))
cat(sprintf("median selected density    : %.5f\n", median(dt$density)))
cat(sprintf("median edge_count          : %d\n",   as.integer(median(dt$edge_count))))
cat(sprintf("median mean_degree         : %.1f\n", median(dt$mean_degree)))
cat(sprintf("median max_degree          : %d\n",   as.integer(median(dt$max_degree))))
cat(sprintf("median s0 (max+1)          : %d\n",   as.integer(median(dt$s0_estim))))
cat(sprintf("median q90                 : %.1f\n", median(dt$degree_q90)))
cat(sprintf("median q99                 : %.1f\n", median(dt$degree_q99)))
cat(sprintf("median ratio max/q50       : %.2f\n",
            median(dt$max_degree / pmax(dt$degree_q50, 1))))
cat(sprintf("median ratio q99/q50       : %.2f\n",
            median(dt$degree_q99 / pmax(dt$degree_q50, 1))))
cat(sprintf("median bisect iters        : %.0f  (range %d-%d)\n",
            median(dt$bisect_iters), min(dt$bisect_iters), max(dt$bisect_iters)))
cat("\n")

cat("--- Per-subject decisions (BISECTED) ---\n")
print(dt[, .(subject_id, selected_nu = signif(selected_nu, 5),
             density = signif(density, 4), edges = edge_count,
             mean_deg = signif(mean_degree, 4),
             q50 = degree_q50, q90 = degree_q90, q99 = degree_q99,
             max = max_degree, s0 = s0_estim,
             m_p = signif(movie_p_raw, 3),
             s5 = ifelse(movie_self_reject_005, "rej", "pass"),
             s10 = ifelse(movie_self_reject_010, "rej", "pass"),
             r_p_raw = signif(recall_p_raw, 3),
             r_p_holm = signif(recall_p_holm, 3),
             r_holm = recall_reject_holm_005)])
cat("\n")
cat("--- Bracket -> bisected nu (compare grid-based vs bisected) ---\n")
print(dt[, .(subject_id,
             bracket_lo_initial = signif(bracket_lo_initial, 5),
             bracket_hi_initial = signif(bracket_hi_initial, 5),
             selected_nu_bisect = signif(selected_nu, 5),
             bracket_width_initial = signif(bracket_hi_initial - bracket_lo_initial, 4),
             bracket_width_final  = signif(bracket_hi_initial - selected_nu, 4),
             gain_in_nu = signif(selected_nu - bracket_lo_initial, 4),
             bisect_iters)])
cat("\n")

cat("============================================================\n")
cat("Subject-level CSV : ",
    normalizePath(file.path(cfg$out_dir, "tables", "selected_subject_results_omega_bisect.csv"), winslash="/"), "\n", sep="")
cat("Recall test CSV   : ",
    normalizePath(file.path(cfg$out_dir, "tables", "recall_test_results_omega_bisect.csv"), winslash="/"), "\n", sep="")
sink()

log_msg("Done.")
