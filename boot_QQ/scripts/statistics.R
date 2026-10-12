matrix_stats_available <- requireNamespace("matrixStats", quietly = TRUE)

column_cumsums <- function(x) {
  if (matrix_stats_available) {
    return(matrixStats::colCumsums(x))
  }
  apply(x, 2L, cumsum)
}

row_max_abs <- function(x) {
  if (matrix_stats_available) {
    return(matrixStats::rowMaxs(abs(x)))
  }
  apply(abs(x), 1L, max)
}

build_graph_structure <- function(p, b) {
  supports <- vector("list", p)
  support_positions <- integer(p)

  for (i in seq_len(p)) {
    support <- max(1L, i - b):min(p, i + b)
    supports[[i]] <- support
    support_positions[i] <- match(i, support)
  }

  list(
    p = p,
    b = b,
    supports = supports,
    support_positions = support_positions
  )
}

off_support_for_index <- function(i, p, b) {
  left_end <- i - b - 1L
  right_start <- i + b + 1L

  if (left_end >= 1L && right_start <= p) {
    return(c(seq_len(left_end), seq.int(right_start, p)))
  }
  if (left_end >= 1L) {
    return(seq_len(left_end))
  }
  if (right_start <= p) {
    return(seq.int(right_start, p))
  }
  integer(0)
}

compute_plugin_estimator <- function(sigma_hat, graph, ridge = 1e-8) {
  p <- ncol(sigma_hat)
  local_inverses <- vector("list", p)
  local_solutions <- vector("list", p)
  fallback_count <- 0L

  for (i in seq_len(p)) {
    support <- graph$supports[[i]]
    block <- sigma_hat[support, support, drop = FALSE]
    block_inverse <- tryCatch(
      solve(block),
      error = function(...) NULL
    )
    if (is.null(block_inverse)) {
      block_inverse <- solve(block + ridge * diag(length(support)))
      fallback_count <- fallback_count + 1L
    }

    e_local <- rep(0, length(support))
    e_local[graph$support_positions[i]] <- 1
    solution <- as.numeric(block_inverse %*% e_local)
    local_solutions[[i]] <- solution
    local_inverses[[i]] <- block_inverse
  }

  list(
    local_inverses = local_inverses,
    local_solutions = local_solutions,
    fallback_count = fallback_count
  )
}

compute_t_hat <- function(sigma_hat, plugin, graph, n) {
  max_value <- 0
  for (i in seq_along(graph$supports)) {
    off_support <- off_support_for_index(i, graph$p, graph$b)
    if (length(off_support) == 0L) {
      next
    }
    support <- graph$supports[[i]]
    residual_i <- as.numeric(
      sigma_hat[off_support, support, drop = FALSE] %*% plugin$local_solutions[[i]]
    )
    max_value <- max(max_value, max(abs(residual_i)))
  }
  sqrt(n) * max_value
}

compute_bootstrap_stats <- function(x, sigma_hat, plugin, graph, block_length) {
  n <- nrow(x)
  n_blocks <- n - block_length + 1L
  if (n_blocks <= 0L) {
    stop("block_length must be at most n.", call. = FALSE)
  }

  block_max <- numeric(n_blocks)

  for (i in seq_along(graph$supports)) {
    support <- graph$supports[[i]]
    off_support <- off_support_for_index(i, graph$p, graph$b)
    if (length(off_support) == 0L) {
      next
    }

    x_support <- x[, support, drop = FALSE]
    q_i <- as.numeric(x_support %*% plugin$local_solutions[[i]])
    sigma_support_off <- sigma_hat[support, off_support, drop = FALSE]
    correction <- plugin$local_inverses[[i]] %*% sigma_support_off
    utx <- x[, off_support, drop = FALSE] - x_support %*% correction
    scores <- utx * q_i

    csum <- rbind(0, column_cumsums(scores))
    block_sums <- csum[(block_length + 1L):(n + 1L), , drop = FALSE] -
      csum[seq_len(n_blocks), , drop = FALSE]
    block_max <- pmax(block_max, row_max_abs(block_sums))
  }

  block_max / sqrt(block_length)
}

compute_replication_statistics <- function(x, graph, alpha, block_length, center = FALSE) {
  if (isTRUE(center)) {
    x <- scale(x, center = TRUE, scale = FALSE)
  }
  n <- nrow(x)
  sigma_hat <- crossprod(x) / n
  plugin <- compute_plugin_estimator(sigma_hat = sigma_hat, graph = graph)
  t_hat <- compute_t_hat(sigma_hat = sigma_hat, plugin = plugin, graph = graph, n = n)
  boot_stats <- compute_bootstrap_stats(
    x = x,
    sigma_hat = sigma_hat,
    plugin = plugin,
    graph = graph,
    block_length = block_length
  )
  q_hat <- as.numeric(quantile(boot_stats, probs = 1 - alpha, names = FALSE, type = 8))

  list(
    T_hat = as.numeric(t_hat),
    q_hat = q_hat,
    reject = as.integer(t_hat > q_hat),
    boot_stats = boot_stats,
    fallback_count = plugin$fallback_count
  )
}
