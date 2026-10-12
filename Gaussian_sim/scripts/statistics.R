pinv_base <- function(x, tol = sqrt(.Machine$double.eps)) {
  s <- svd(x)
  if (length(s$d) == 0L) {
    return(matrix(0, nrow = ncol(x), ncol = nrow(x)))
  }
  threshold <- tol * max(s$d)
  keep <- s$d > threshold
  if (!any(keep)) {
    return(matrix(0, nrow = ncol(x), ncol = nrow(x)))
  }
  s$v[, keep, drop = FALSE] %*%
    diag(1 / s$d[keep], nrow = sum(keep), ncol = sum(keep)) %*%
    t(s$u[, keep, drop = FALSE])
}

build_graph_structure <- function(sigma, omega, b) {
  p <- ncol(sigma)
  supports <- vector("list", p)
  off_supports <- vector("list", p)
  support_positions <- integer(p)
  u_transforms <- vector("list", p)
  off_support_mask <- matrix(FALSE, nrow = p, ncol = p)

  for (i in seq_len(p)) {
    left <- max(1L, i - b)
    right <- min(p, i + b)
    support <- left:right
    off_support <- setdiff(seq_len(p), support)
    sigma_ss <- sigma[support, support, drop = FALSE]
    inv_sigma_ss <- solve(sigma_ss)

    supports[[i]] <- support
    off_supports[[i]] <- off_support
    support_positions[i] <- match(i, support)
    u_transforms[[i]] <- sigma[, support, drop = FALSE] %*% inv_sigma_ss
    if (length(off_support) > 0L) {
      off_support_mask[off_support, i] <- TRUE
    }
  }

  list(
    bandwidth = b,
    supports = supports,
    off_supports = off_supports,
    support_positions = support_positions,
    u_transforms = u_transforms,
    w_true = omega,
    off_support_mask = off_support_mask
  )
}

compute_plugin_estimator <- function(sigma_hat, graph) {
  p <- ncol(sigma_hat)
  w_hat <- matrix(0, nrow = p, ncol = p)
  fallback_count <- 0L

  for (i in seq_len(p)) {
    support <- graph$supports[[i]]
    block <- sigma_hat[support, support, drop = FALSE]
    e_local <- rep(0, length(support))
    e_local[graph$support_positions[[i]]] <- 1

    solution <- tryCatch(
      solve(block, e_local),
      error = function(...) NULL
    )
    if (is.null(solution)) {
      solution <- pinv_base(block) %*% e_local
      fallback_count <- fallback_count + 1L
    }
    w_hat[support, i] <- as.numeric(solution)
  }

  list(w_hat = w_hat, fallback_count = fallback_count)
}

right_multiply_supported_columns <- function(left_matrix, column_matrix, supports) {
  p <- ncol(column_matrix)
  out <- matrix(0, nrow = nrow(left_matrix), ncol = p)

  for (i in seq_len(p)) {
    support <- supports[[i]]
    out[, i] <- as.numeric(left_matrix[, support, drop = FALSE] %*% column_matrix[support, i, drop = FALSE])
  }

  out
}

evaluate_off_support_max <- function(column_matrix, graph) {
  max_value <- 0
  p <- ncol(column_matrix)

  for (i in seq_len(p)) {
    off_support <- graph$off_supports[[i]]
    if (length(off_support) == 0L) {
      next
    }
    column_i <- column_matrix[, i]
    adjusted <- column_i - graph$u_transforms[[i]] %*% column_i[graph$supports[[i]]]
    current <- max(abs(adjusted[off_support]))
    if (current > max_value) {
      max_value <- current
    }
  }

  as.numeric(max_value)
}

compute_replication_statistics <- function(x, sigma, graph) {
  n <- nrow(x)
  p <- ncol(x)
  sigma_hat <- crossprod(x) / n

  plugin <- compute_plugin_estimator(sigma_hat = sigma_hat, graph = graph)
  residual <- right_multiply_supported_columns(
    left_matrix = sigma_hat,
    column_matrix = plugin$w_hat,
    supports = graph$supports
  )
  if (any(graph$off_support_mask)) {
    t_hat <- sqrt(n) * max(abs(residual[graph$off_support_mask]))
  } else {
    t_hat <- 0
  }

  delta <- sigma_hat - sigma
  linear_columns <- right_multiply_supported_columns(
    left_matrix = delta,
    column_matrix = graph$w_true,
    supports = graph$supports
  )
  t_l <- sqrt(n) * evaluate_off_support_max(linear_columns, graph)

  list(
    T_hat = as.numeric(t_hat),
    T_L = as.numeric(t_l),
    fallback_count = plugin$fallback_count
  )
}
