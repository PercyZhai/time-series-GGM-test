comparison_title <- function(
  statistic_name,
  target_name,
  n,
  p,
  b,
  beta,
  N_mc,
  B_gauss
) {
  sprintf(
    "%s vs %s\nn=%d, p=%d, b=%d, beta=%s, N_mc=%d, B_gauss=%d",
    statistic_name,
    target_name,
    n,
    p,
    b,
    format_beta(beta),
    N_mc,
    B_gauss
  )
}

comparison_caption <- function(metrics) {
  sprintf(
    "Kolmogorov distance = %.4f    Wasserstein-1 = %.4f",
    metrics$kolmogorov_distance,
    metrics$wasserstein1_distance
  )
}

open_plot_device <- function(output_path, width, height, res) {
  if (grepl("\\.pdf$", output_path, ignore.case = TRUE)) {
    pdf(file = output_path, width = width / res, height = height / res)
    return(invisible(NULL))
  }
  png_args <- list(filename = output_path, width = width, height = height, res = res)
  if (isTRUE(capabilities("cairo"))) png_args$type <- "cairo"
  do.call(png, png_args)
}

save_qq_plot <- function(
  statistic_samples,
  target_samples,
  statistic_name,
  target_name,
  metrics,
  n,
  p,
  b,
  beta,
  N_mc,
  B_gauss,
  output_path
) {
  q_count <- min(length(statistic_samples), length(target_samples))
  probs <- ((seq_len(q_count) - 0.5) / q_count)
  x_quantiles <- as.numeric(quantile(target_samples, probs = probs, names = FALSE, type = 8))
  y_quantiles <- as.numeric(quantile(statistic_samples, probs = probs, names = FALSE, type = 8))

  lims <- range(c(x_quantiles, y_quantiles))
  if (diff(lims) == 0) {
    padding <- if (lims[1] == 0) 1 else 0.05 * abs(lims[1])
    lims <- c(lims[1] - padding, lims[2] + padding)
  }

  open_plot_device(output_path, width = 1400, height = 1400, res = 180)
  par(mar = c(5.2, 5.2, 5.5, 2.2))
  plot(
    x_quantiles,
    y_quantiles,
    pch = 19,
    cex = 0.8,
    col = "#1f77b4",
    xlab = sprintf("%s quantiles", target_name),
    ylab = sprintf("%s quantiles", statistic_name),
    xlim = lims,
    ylim = lims,
    asp = 1,
    main = comparison_title(statistic_name, target_name, n, p, b, beta, N_mc, B_gauss),
    sub = comparison_caption(metrics)
  )
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "#d62728")
  dev.off()
}

save_cdf_plot <- function(
  statistic_samples,
  target_samples,
  statistic_name,
  target_name,
  metrics,
  n,
  p,
  b,
  beta,
  N_mc,
  B_gauss,
  output_path
) {
  statistic_sorted <- sort(as.numeric(statistic_samples))
  target_sorted <- sort(as.numeric(target_samples))
  grid <- sort(unique(c(statistic_sorted, target_sorted)))
  statistic_cdf <- findInterval(grid, statistic_sorted, rightmost.closed = TRUE) / length(statistic_sorted)
  target_cdf <- findInterval(grid, target_sorted, rightmost.closed = TRUE) / length(target_sorted)

  open_plot_device(output_path, width = 1400, height = 1000, res = 180)
  par(mar = c(5.2, 5.2, 5.5, 2.2))
  plot(
    grid,
    target_cdf,
    type = "s",
    lwd = 2,
    col = "#ff7f0e",
    ylim = c(0, 1),
    xlab = "Value",
    ylab = "Empirical CDF",
    main = comparison_title(statistic_name, target_name, n, p, b, beta, N_mc, B_gauss),
    sub = comparison_caption(metrics)
  )
  lines(grid, statistic_cdf, type = "s", lwd = 2, col = "#1f77b4")
  legend(
    "bottomright",
    legend = c(target_name, statistic_name),
    col = c("#ff7f0e", "#1f77b4"),
    lwd = 2,
    bty = "n"
  )
  dev.off()
}
