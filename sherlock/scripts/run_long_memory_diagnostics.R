#!/usr/bin/env Rscript
# run_long_memory_diagnostics.R
#
# Long-memory diagnostics for Sherlock movie-phase fMRI (analogous to
# IEEE long-memory paper Figs 10-11):
#   - DFA Hurst exponent for every (subject, ROI) movie time series
#   - GPH fractional memory parameter d (robustness)
#   - Distribution figures + ACF panels at percentile-selected series
#
# Outputs under Results/:
#   long_memory_estimates.csv     subject_id, roi_id, roi_label, n_timepoints,
#                                 dropped_zero_var, hurst_dfa, d_gph
#   figures/fig1_hurst_distribution.pdf
#   figures/fig2_acf_panels.pdf
#   figures/fig3_d_distribution.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
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
  data_path = file.path(PROJECT_ROOT, "processed_fMRI_data", "moviefmrip268n17T1976.RData"),
  out_dir   = file.path(PROJECT_ROOT, "Results"),
  detrend   = TRUE,
  standardize = TRUE,
  acf_max_lag = 100L
)
set.seed(1)

dir.create(file.path(cfg$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load movie array (auto-orient: subj, ROI, time axes)
# ---------------------------------------------------------------------------
e <- new.env()
load(cfg$data_path, envir = e)
arr <- NULL
for (nm in ls(e)) {
  obj <- get(nm, envir = e)
  if (is.array(obj) && length(dim(obj)) == 3L &&
      all(c(268L, 1976L, 17L) %in% dim(obj))) {
    arr <- obj; break
  }
}
stopifnot(!is.null(arr))
d <- dim(arr)
pos_subj <- which(d == 17L)[1]
pos_p    <- which(d == 268L)[1]
pos_T    <- setdiff(1:3, c(pos_subj, pos_p))[1]
cat(sprintf("[load] dims=%s  subj@%d  ROI@%d  time@%d\n",
            paste(d, collapse="x"), pos_subj, pos_p, pos_T))

get_subject_TxP <- function(s) {
  idx <- list(quote(expr=), quote(expr=), quote(expr=))
  idx[[pos_subj]] <- s
  slab <- do.call(`[`, c(list(arr), idx, list(drop = TRUE)))
  if (nrow(slab) == 268L) slab <- t(slab)  # ensure T x p
  storage.mode(slab) <- "double"
  slab
}

# ---------------------------------------------------------------------------
# Helpers: detrend, standardize, DFA, GPH
# ---------------------------------------------------------------------------
detrend_linear <- function(x) {
  t_ <- seq_along(x)
  beta <- .lm.fit(cbind(1, t_), x)$coefficients
  as.numeric(x - (beta[1] + beta[2] * t_))
}
zstd <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NULL)
  (x - mean(x)) / s
}

# DFA: Detrended Fluctuation Analysis (Peng et al. 1994 / Kantelhardt 2001)
# Returns alpha = slope of log F(s) vs log s; alpha ~ Hurst H for stationary.
dfa_hurst <- function(x, scales = NULL, min_scale = 8L) {
  n <- length(x)
  if (is.null(scales)) {
    s_max <- floor(n / 4)
    if (s_max < min_scale + 1) return(NA_real_)
    scales <- unique(round(2 ^ seq(log2(min_scale), log2(s_max),
                                   length.out = 20L)))
  }
  y <- cumsum(x - mean(x))
  Fvals <- numeric(length(scales))
  for (i in seq_along(scales)) {
    s <- scales[i]
    n_seg <- floor(n / s)
    if (n_seg < 1) { Fvals[i] <- NA_real_; next }
    rms_sq <- numeric(2L * n_seg)
    cnt <- 0L
    t_seg <- seq_len(s)
    Xmat <- cbind(1, t_seg)
    # forward
    for (k in seq_len(n_seg)) {
      idx <- ((k - 1L) * s + 1L):(k * s)
      seg <- y[idx]
      beta <- .lm.fit(Xmat, seg)$coefficients
      r <- seg - (beta[1] + beta[2] * t_seg)
      cnt <- cnt + 1L; rms_sq[cnt] <- mean(r * r)
    }
    # backward
    for (k in seq_len(n_seg)) {
      idx <- (n - k * s + 1L):(n - (k - 1L) * s)
      seg <- y[idx]
      beta <- .lm.fit(Xmat, seg)$coefficients
      r <- seg - (beta[1] + beta[2] * t_seg)
      cnt <- cnt + 1L; rms_sq[cnt] <- mean(r * r)
    }
    Fvals[i] <- sqrt(mean(rms_sq[seq_len(cnt)]))
  }
  ok <- is.finite(Fvals) & Fvals > 0
  if (sum(ok) < 4L) return(NA_real_)
  fit <- .lm.fit(cbind(1, log(scales[ok])), log(Fvals[ok]))
  as.numeric(fit$coefficients[2])
}

# GPH (Geweke-Porter-Hudak 1983): regress log periodogram on log(4 sin^2(lambda/2))
# d_gph = -slope; H_gph = d_gph + 0.5 for stationary fractional processes.
gph_d <- function(x, m = NULL) {
  n <- length(x)
  if (is.null(m)) m <- floor(n^0.5)
  if (m < 5L) return(NA_real_)
  z <- fft(x - mean(x))
  I <- (Mod(z)^2) / (2 * pi * n)
  j <- 1:m
  lambda <- 2 * pi * j / n
  Y <- log(I[j + 1L])
  X <- log(4 * sin(lambda / 2)^2)
  ok <- is.finite(Y) & is.finite(X)
  if (sum(ok) < 4L) return(NA_real_)
  fit <- .lm.fit(cbind(1, X[ok]), Y[ok])
  -as.numeric(fit$coefficients[2])
}

# ---------------------------------------------------------------------------
# Per-subject driver
# ---------------------------------------------------------------------------
process_subject <- function(s) {
  X <- get_subject_TxP(s)  # T x p
  n <- nrow(X); p <- ncol(X)
  rows <- vector("list", p)
  for (j in seq_len(p)) {
    x <- X[, j]
    v <- var(x)
    dropped <- !is.finite(v) || v == 0
    if (dropped) {
      rows[[j]] <- data.table(subject_id = s, roi_id = j, roi_label = NA_character_,
                              n_timepoints = n, dropped_zero_var = TRUE,
                              hurst_dfa = NA_real_, d_gph = NA_real_)
      next
    }
    if (cfg$detrend) x <- detrend_linear(x)
    if (cfg$standardize) {
      xs <- zstd(x)
      if (is.null(xs)) {
        rows[[j]] <- data.table(subject_id = s, roi_id = j, roi_label = NA_character_,
                                n_timepoints = n, dropped_zero_var = TRUE,
                                hurst_dfa = NA_real_, d_gph = NA_real_)
        next
      }
      x <- xs
    }
    H  <- tryCatch(dfa_hurst(x), error = function(e) NA_real_)
    dh <- tryCatch(gph_d(x),     error = function(e) NA_real_)
    rows[[j]] <- data.table(subject_id = s, roi_id = j, roi_label = NA_character_,
                            n_timepoints = n, dropped_zero_var = FALSE,
                            hurst_dfa = H, d_gph = dh)
  }
  rbindlist(rows)
}

cat("[run] estimating long-memory parameters across 17 x 268 series ...\n")
t0 <- Sys.time()
n_workers <- min(8L, parallel::detectCores() - 1L, 17L)
all_rows <- mclapply(seq_len(17L), process_subject, mc.cores = n_workers)
dt <- rbindlist(all_rows)
cat(sprintf("[run] done in %.1f s  (n series=%d, n dropped=%d)\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            nrow(dt), sum(dt$dropped_zero_var)))

setorder(dt, subject_id, roi_id)
fwrite(dt, file.path(cfg$out_dir, "tables", "long_memory_estimates.csv"))

# ---------------------------------------------------------------------------
# Summary stats
# ---------------------------------------------------------------------------
H_all <- dt$hurst_dfa[!is.na(dt$hurst_dfa)]
d_all <- dt$d_gph[!is.na(dt$d_gph)]
H_med <- median(H_all); H_iqr <- quantile(H_all, c(0.25, 0.75), names = FALSE)
H_frac_gt_05 <- mean(H_all > 0.5)
d_med <- median(d_all); d_iqr <- quantile(d_all, c(0.25, 0.75), names = FALSE)
d_frac_gt_0  <- mean(d_all > 0)
cat(sprintf("[summary] H median=%.3f  IQR=[%.3f, %.3f]  P(H>0.5)=%.3f  n=%d\n",
            H_med, H_iqr[1], H_iqr[2], H_frac_gt_05, length(H_all)))
cat(sprintf("[summary] d median=%.3f  IQR=[%.3f, %.3f]  P(d>0)=%.3f  n=%d\n",
            d_med, d_iqr[1], d_iqr[2], d_frac_gt_0, length(d_all)))

# ---------------------------------------------------------------------------
# Figure 1: Hurst distribution (pooled hist + per-subject boxplot)
# ---------------------------------------------------------------------------
fig1_caption <- sprintf(
  "Figure 1. Distribution of DFA Hurst exponents across all subject-ROI movie time series\n(n = %d valid series). Vertical dashed line: H = 0.5. Median H = %.3f, IQR [%.3f, %.3f],\nfraction with H > 0.5 = %.3f. Values H > 0.5 are consistent with long-range dependence.",
  length(H_all), H_med, H_iqr[1], H_iqr[2], H_frac_gt_05)

draw_fig1 <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.4, 1))
  par(mar = c(4.2, 4.2, 2.0, 0.6))
  hist(H_all, breaks = 60, col = "grey80", border = "grey40",
       xlab = "Hurst exponent (DFA)", ylab = "count",
       main = "Pooled histogram (all 17 x 268 series)")
  abline(v = 0.5, lty = 2, col = "red", lwd = 1.5)
  abline(v = H_med, lty = 1, col = "blue", lwd = 1.2)
  legend("topright", legend = c("H = 0.5", sprintf("median = %.3f", H_med)),
         lty = c(2, 1), col = c("red", "blue"), bty = "n", cex = 0.9)

  par(mar = c(4.2, 4.2, 2.0, 0.6))
  boxplot(hurst_dfa ~ subject_id, data = dt[!is.na(hurst_dfa)],
          xlab = "subject_id", ylab = "Hurst exponent (DFA)",
          col = "grey80", main = "Per-subject distribution",
          outpch = 20, outcex = 0.4)
  abline(h = 0.5, lty = 2, col = "red", lwd = 1.2)
}

pdf(file.path(cfg$out_dir, "figures", "fig1_hurst_distribution.pdf"),
    width = 11, height = 5)
draw_fig1()
dev.off()
cat(sprintf("[fig1] %s\n%s\n", "fig1_hurst_distribution", fig1_caption))

# ---------------------------------------------------------------------------
# Figure 2: ACF panels at percentiles of H distribution
# ---------------------------------------------------------------------------
percs <- c(0.50, 0.75, 0.90, 0.95)
H_targets <- quantile(H_all, percs, names = FALSE)
valid <- dt[!is.na(hurst_dfa)]
sel_rows <- vector("list", length(percs))
for (i in seq_along(percs)) {
  k <- which.min(abs(valid$hurst_dfa - H_targets[i]))
  sel_rows[[i]] <- valid[k, ]
}
sel_dt <- rbindlist(sel_rows)
sel_dt[, percentile := percs * 100]

# Compute ACFs for selected series
acf_max_lag <- min(cfg$acf_max_lag, floor(1976 / 2))
acf_list <- vector("list", nrow(sel_dt))
for (i in seq_len(nrow(sel_dt))) {
  s <- sel_dt$subject_id[i]; j <- sel_dt$roi_id[i]
  X <- get_subject_TxP(s)
  x <- X[, j]
  if (cfg$detrend) x <- detrend_linear(x)
  if (cfg$standardize) x <- zstd(x)
  a <- acf(x, lag.max = acf_max_lag, plot = FALSE)
  acf_list[[i]] <- list(lag = drop(a$lag), acf = drop(a$acf), n = length(x))
}
y_min <- min(sapply(acf_list, function(z) min(z$acf)))
y_max <- max(sapply(acf_list, function(z) max(z$acf)))
y_lim <- range(c(y_min - 0.05, y_max + 0.02, -0.2, 1.0))

fig2_caption <- sprintf(
  "Figure 2. Sample autocorrelation functions for movie time series at H-distribution percentiles\n(50th, 75th, 90th, 95th). Each panel shows ACF up to lag %d; dashed bands at +/- 1.96/sqrt(n).\nSlowly-decaying significant ACF coefficients are evidence consistent with long-range dependence.",
  acf_max_lag)

draw_fig2 <- function() {
  par(mfrow = c(2, 2), mar = c(4.0, 4.0, 2.4, 0.6), oma = c(0, 0, 0, 0))
  for (i in seq_len(nrow(sel_dt))) {
    z <- acf_list[[i]]
    band <- 1.96 / sqrt(z$n)
    plot(z$lag, z$acf, type = "h", col = "grey20", lwd = 1.2,
         xlim = c(0, acf_max_lag), ylim = y_lim,
         xlab = "lag", ylab = "ACF",
         main = sprintf("p%02d  subj %02d  ROI %3d  H=%.3f  n=%d",
                        sel_dt$percentile[i],
                        sel_dt$subject_id[i],
                        sel_dt$roi_id[i],
                        sel_dt$hurst_dfa[i],
                        z$n))
    abline(h = 0, col = "black", lwd = 0.6)
    abline(h = c(-band, band), col = "red", lty = 2, lwd = 1.0)
  }
}

pdf(file.path(cfg$out_dir, "figures", "fig2_acf_panels.pdf"),
    width = 10, height = 7.5)
draw_fig2()
dev.off()
cat(sprintf("[fig2] %s\n%s\n", "fig2_acf_panels", fig2_caption))
cat("[fig2 selected series]\n"); print(sel_dt[, .(percentile, subject_id, roi_id, hurst_dfa, n_timepoints)])

# ---------------------------------------------------------------------------
# Figure 3: GPH d distribution (robustness)
# ---------------------------------------------------------------------------
fig3_caption <- sprintf(
  "Figure 3. Distribution of GPH fractional memory parameter d across all subject-ROI movie\ntime series (n = %d). H ~ d + 0.5 for stationary fractional processes. Vertical dashed line at d = 0.\nMedian d = %.3f, IQR [%.3f, %.3f], fraction with d > 0 = %.3f.",
  length(d_all), d_med, d_iqr[1], d_iqr[2], d_frac_gt_0)

draw_fig3 <- function() {
  par(mar = c(4.2, 4.2, 2.0, 0.6))
  hist(d_all, breaks = 60, col = "grey80", border = "grey40",
       xlab = "GPH fractional memory parameter d",
       ylab = "count",
       main = "Distribution of d_GPH (movie time series)")
  abline(v = 0, lty = 2, col = "red", lwd = 1.5)
  abline(v = d_med, lty = 1, col = "blue", lwd = 1.2)
  legend("topright", legend = c("d = 0", sprintf("median = %.3f", d_med)),
         lty = c(2, 1), col = c("red", "blue"), bty = "n", cex = 0.9)
}

pdf(file.path(cfg$out_dir, "figures", "fig3_d_distribution.pdf"),
    width = 7, height = 5)
draw_fig3()
dev.off()
cat(sprintf("[fig3] %s\n%s\n", "fig3_d_distribution", fig3_caption))

# ---------------------------------------------------------------------------
# Captions written to file
# ---------------------------------------------------------------------------
cat(paste(fig1_caption, "\n\n", fig2_caption, "\n\n", fig3_caption, sep = ""),
    file = file.path(cfg$out_dir, "figures", "captions.txt"))

cat("\n[done]  CSV: ", file.path(cfg$out_dir, "tables", "long_memory_estimates.csv"),
    "\n        figs: ", file.path(cfg$out_dir, "figures"), "\n", sep = "")
