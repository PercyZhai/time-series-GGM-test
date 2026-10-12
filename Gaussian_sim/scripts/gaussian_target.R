right_multiply_upper_banded <- function(left_matrix, r_upper, b) {
  p <- nrow(r_upper)
  out <- matrix(0, nrow = nrow(left_matrix), ncol = p)

  for (j in seq_len(p)) {
    rows <- max(1L, j - b):j
    out[, j] <- as.numeric(left_matrix[, rows, drop = FALSE] %*% r_upper[rows, j, drop = FALSE])
  }

  out
}

simulate_base_gaussian_suprema <- function(b_gauss, r_upper, graph, cores = 1L) {
  p <- nrow(r_upper)
  samples <- numeric(b_gauss)
  cores <- max(1L, min(as.integer(cores), b_gauss))
  pb <- txtProgressBar(min = 0, max = b_gauss, style = 3)

  evaluate_draw <- function(a_mat) {
    z_mat <- (a_mat + t(a_mat)) / sqrt(2)
    right_product <- right_multiply_upper_banded(
      left_matrix = z_mat,
      r_upper = r_upper,
      b = graph$bandwidth
    )
    transformed <- solve_upper_banded(
      rhs = right_product,
      r_upper = r_upper,
      b = graph$bandwidth
    )
    evaluate_off_support_max(transformed, graph)
  }

  draw_batches <- split(seq_len(b_gauss), ceiling(seq_len(b_gauss) / cores))
  completed <- 0L
  for (draw_idx in draw_batches) {
    # Parent-side RNG preserves the serial draw sequence for a fixed seed.
    gaussian_matrices <- lapply(draw_idx, function(i) {
      matrix(rnorm(p * p), nrow = p, ncol = p)
    })
    draw_values <- if (cores == 1L || length(draw_idx) == 1L) {
      lapply(gaussian_matrices, evaluate_draw)
    } else {
      parallel::mclapply(
        gaussian_matrices,
        evaluate_draw,
        mc.cores = min(cores, length(draw_idx)),
        mc.preschedule = TRUE,
        mc.set.seed = FALSE
      )
    }
    samples[draw_idx] <- unlist(draw_values, use.names = FALSE)
    completed <- completed + length(draw_idx)
    setTxtProgressBar(pb, completed)
  }

  close(pb)
  samples
}

simulate_gaussian_targets <- function(
  b_gauss,
  l_finite,
  l_asymptotic,
  r_upper,
  graph,
  use_asymptotic_g = FALSE,
  cores = 1L
) {
  base_samples <- simulate_base_gaussian_suprema(
    b_gauss = b_gauss,
    r_upper = r_upper,
    graph = graph,
    cores = cores
  )

  out <- list(
    G_finite_samples = sqrt(max(l_finite, 0)) * base_samples,
    G_asymptotic_samples = NULL
  )

  if (isTRUE(use_asymptotic_g)) {
    out$G_asymptotic_samples <- sqrt(max(l_asymptotic, 0)) * base_samples
  }

  out
}
