kolmogorov_distance <- function(sample_a, sample_b) {
  a_sorted <- sort(as.numeric(sample_a))
  b_sorted <- sort(as.numeric(sample_b))
  grid <- sort(unique(c(a_sorted, b_sorted)))
  if (length(grid) == 0L) {
    return(0)
  }
  cdf_a <- findInterval(grid, a_sorted, rightmost.closed = TRUE) / length(a_sorted)
  cdf_b <- findInterval(grid, b_sorted, rightmost.closed = TRUE) / length(b_sorted)
  as.numeric(max(abs(cdf_a - cdf_b)))
}

wasserstein1_distance <- function(sample_a, sample_b) {
  a_sorted <- sort(as.numeric(sample_a))
  b_sorted <- sort(as.numeric(sample_b))
  grid <- sort(unique(c(a_sorted, b_sorted)))
  if (length(grid) <= 1L) {
    return(0)
  }
  cdf_a <- findInterval(grid, a_sorted, rightmost.closed = TRUE) / length(a_sorted)
  cdf_b <- findInterval(grid, b_sorted, rightmost.closed = TRUE) / length(b_sorted)
  interval_widths <- diff(grid)
  as.numeric(sum(abs(cdf_a[-length(grid)] - cdf_b[-length(grid)]) * interval_widths))
}

sample_sd <- function(x) {
  x <- as.numeric(x)
  if (length(x) <= 1L) {
    return(0)
  }
  as.numeric(sd(x))
}

make_summary_row <- function(
  config,
  statistic_samples,
  bootstrap_samples,
  q_values,
  reject_values,
  fallback_count,
  seed
) {
  list(
    n = as.integer(config$n),
    p = as.integer(config$p),
    b = as.integer(config$b),
    beta = as.numeric(config$beta),
    N_mc = as.integer(config$N_mc),
    N_trunc = as.integer(config$N_trunc),
    alpha = as.numeric(config$alpha),
    block_length = as.integer(config$block_length),
    N_boot = as.integer(length(bootstrap_samples)),
    statistic = "T_hat",
    target = "Bootstrap",
    empirical_size = as.numeric(mean(reject_values)),
    kolmogorov_distance = kolmogorov_distance(statistic_samples, bootstrap_samples),
    wasserstein1_distance = wasserstein1_distance(statistic_samples, bootstrap_samples),
    mean_statistic = as.numeric(mean(statistic_samples)),
    mean_bootstrap = as.numeric(mean(bootstrap_samples)),
    mean_q_hat = as.numeric(mean(q_values)),
    sd_statistic = sample_sd(statistic_samples),
    sd_bootstrap = sample_sd(bootstrap_samples),
    sd_q_hat = sample_sd(q_values),
    fallback_count = as.integer(fallback_count),
    seed = as.integer(seed)
  )
}

make_quantile_table <- function(statistic_samples, bootstrap_samples) {
  probs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
  statistic_quantiles <- as.numeric(quantile(statistic_samples, probs = probs, names = FALSE, type = 8))
  bootstrap_quantiles <- as.numeric(quantile(bootstrap_samples, probs = probs, names = FALSE, type = 8))
  data.frame(
    probability = probs,
    T_hat_quantile = statistic_quantiles,
    bootstrap_quantile = bootstrap_quantiles,
    difference_T_minus_bootstrap = statistic_quantiles - bootstrap_quantiles,
    ratio_T_over_bootstrap = statistic_quantiles / bootstrap_quantiles,
    stringsAsFactors = FALSE
  )
}
