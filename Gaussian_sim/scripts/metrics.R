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

summarize_comparison <- function(
  n,
  p,
  b,
  beta,
  N_mc,
  B_gauss,
  N_trunc,
  L_n,
  L,
  statistic,
  target,
  statistic_samples,
  target_samples,
  fallback_count,
  seed
) {
  list(
    n = as.integer(n),
    p = as.integer(p),
    b = as.integer(b),
    beta = as.numeric(beta),
    N_mc = as.integer(N_mc),
    B_gauss = as.integer(B_gauss),
    N_trunc = as.integer(N_trunc),
    L_n = as.numeric(L_n),
    L = as.numeric(L),
    statistic = statistic,
    target = target,
    kolmogorov_distance = kolmogorov_distance(statistic_samples, target_samples),
    wasserstein1_distance = wasserstein1_distance(statistic_samples, target_samples),
    mean_statistic = as.numeric(mean(statistic_samples)),
    mean_target = as.numeric(mean(target_samples)),
    sd_statistic = sample_sd(statistic_samples),
    sd_target = sample_sd(target_samples),
    fallback_count = as.integer(fallback_count),
    seed = as.integer(seed)
  )
}
