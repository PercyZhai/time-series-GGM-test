build_banded_precision <- function(p, b, theta0 = 0.25, rho = 0.5, delta0 = 0.10) {
  omega <- matrix(0, nrow = p, ncol = p)
  if (b > 0L) {
    for (distance in seq_len(b)) {
      value <- theta0 * rho^(distance - 1L)
      idx <- seq_len(p - distance)
      omega[cbind(idx, idx + distance)] <- value
      omega[cbind(idx + distance, idx)] <- value
    }
  }
  diag(omega) <- 1 + rowSums(abs(omega)) + delta0
  omega
}

build_banded_cholesky_upper <- function(omega, b) {
  p <- nrow(omega)
  r_upper <- matrix(0, nrow = p, ncol = p)

  for (j in seq_len(p)) {
    prev <- if (b > 0L && j > 1L) max(1L, j - b):(j - 1L) else integer(0)
    diag_adjustment <- if (length(prev) > 0L) sum(r_upper[prev, j]^2) else 0
    r_upper[j, j] <- sqrt(omega[j, j] - diag_adjustment)

    right <- if (b > 0L && j < p) (j + 1L):min(p, j + b) else integer(0)
    for (k in right) {
      shared_start <- if (j > 1L) max(1L, k - b, j - b) else 1L
      shared <- if (j > 1L && shared_start <= j - 1L) shared_start:(j - 1L) else integer(0)
      cross_term <- if (length(shared) > 0L) {
        sum(r_upper[shared, j] * r_upper[shared, k])
      } else {
        0
      }
      r_upper[j, k] <- (omega[j, k] - cross_term) / r_upper[j, j]
    }
  }

  r_upper
}

solve_upper_banded <- function(rhs, r_upper, b) {
  p <- nrow(r_upper)
  out <- rhs

  for (i in p:1L) {
    right <- if (b > 0L && i < p) (i + 1L):min(p, i + b) else integer(0)
    if (length(right) > 0L) {
      out[i, ] <- out[i, ] - as.numeric(r_upper[i, right, drop = FALSE] %*% out[right, , drop = FALSE])
    }
    out[i, ] <- out[i, ] / r_upper[i, i]
  }

  out
}

solve_lower_banded_identity <- function(r_upper, b) {
  p <- nrow(r_upper)
  out <- matrix(0, nrow = p, ncol = p)

  for (i in seq_len(p)) {
    out[i, i] <- 1
    left <- if (b > 0L && i > 1L) max(1L, i - b):(i - 1L) else integer(0)
    if (length(left) > 0L) {
      out[i, ] <- out[i, ] - as.numeric(t(r_upper[left, i, drop = FALSE]) %*% out[left, , drop = FALSE])
    }
    out[i, ] <- out[i, ] / r_upper[i, i]
  }

  out
}

invert_from_banded_cholesky <- function(r_upper, b) {
  lower_solution <- solve_lower_banded_identity(r_upper = r_upper, b = b)
  sigma <- solve_upper_banded(rhs = lower_solution, r_upper = r_upper, b = b)
  (sigma + t(sigma)) / 2
}

right_solve_upper_banded <- function(y, r_upper, b) {
  p <- ncol(y)
  out <- y

  for (j in seq_len(p)) {
    left <- if (b > 0L && j > 1L) max(1L, j - b):(j - 1L) else integer(0)
    if (length(left) > 0L) {
      out[, j] <- out[, j] - as.numeric(out[, left, drop = FALSE] %*% r_upper[left, j, drop = FALSE])
    }
    out[, j] <- out[, j] / r_upper[j, j]
  }

  out
}

build_model_matrices <- function(p, b, theta0 = 0.25, rho = 0.5, delta0 = 0.10) {
  omega <- build_banded_precision(p = p, b = b, theta0 = theta0, rho = rho, delta0 = delta0)
  r_upper <- build_banded_cholesky_upper(omega = omega, b = b)
  sigma <- invert_from_banded_cholesky(r_upper = r_upper, b = b)
  list(Omega = omega, Sigma = sigma, R = r_upper)
}

make_filter_coefficients <- function(beta, n_trunc) {
  lags <- 0:(n_trunc - 1L)
  a_tilde <- (lags + 1)^(-beta)
  a_tilde / sqrt(sum(a_tilde^2))
}

compute_gamma_sequence_fft <- function(a) {
  n_trunc <- length(a)
  fft_length <- 2^ceiling(log2(2 * n_trunc - 1L))
  a_padded <- c(a, rep(0, fft_length - n_trunc))
  f_a <- fft(a_padded)
  power <- Mod(f_a)^2
  autocorr <- Re(fft(power, inverse = TRUE)) / fft_length
  autocorr[1:n_trunc]
}

compute_l_factors <- function(gamma, n) {
  gamma_sq <- gamma^2
  l_asymptotic <- gamma_sq[1] + 2 * sum(gamma_sq[-1])

  max_lag <- min(n - 1L, length(gamma) - 1L)
  if (max_lag <= 0L) {
    return(list(L_n = gamma_sq[1], L = l_asymptotic))
  }

  k_vec <- seq_len(max_lag)
  weights <- 1 - k_vec / n
  l_finite <- gamma_sq[1] + 2 * sum(weights * gamma_sq[k_vec + 1L])
  list(L_n = l_finite, L = l_asymptotic)
}

generate_fft_sample_block <- function(
  n,
  n_mc,
  r_upper,
  b,
  a,
  p,
  verbose = TRUE,
  fft_cores = 1L
) {
  rows_needed <- n * n_mc
  n_trunc <- length(a)
  sequence_length <- rows_needed + n_trunc - 1L
  fft_length <- 2^ceiling(log2(sequence_length + n_trunc - 1L))
  fft_cores <- max(1L, min(as.integer(fft_cores), p))

  if (fft_cores > 1L && .Platform$OS.type != "unix") {
    stop("Parallel FFT generation requires a Unix-like platform.", call. = FALSE)
  }

  if (verbose) {
    message("Generating data using FFT-based scalar filtering...")
  }

  # This follows the E3B_quick.R pattern: pre-transform the filter once,
  # then generate each temporally correlated series by FFT multiplication.
  a_padded <- c(a, rep(0, fft_length - n_trunc))
  b_vec <- fft(a_padded)

  y_block <- matrix(0, nrow = rows_needed, ncol = p)

  filter_one_coordinate <- function(eps_j) {
    eps_padded <- c(eps_j, rep(0, fft_length - sequence_length))
    f_j <- fft(eps_padded)
    w_j_complex <- fft(b_vec * f_j, inverse = TRUE) / fft_length
    Re(w_j_complex)[n_trunc:(n_trunc + rows_needed - 1L)]
  }

  column_batches <- split(seq_len(p), ceiling(seq_len(p) / fft_cores))
  for (column_idx in column_batches) {
    innovations <- lapply(column_idx, function(j) rnorm(sequence_length))
    filtered_columns <- if (fft_cores == 1L || length(column_idx) == 1L) {
      lapply(innovations, filter_one_coordinate)
    } else {
      parallel::mclapply(
        innovations,
        filter_one_coordinate,
        mc.cores = min(fft_cores, length(column_idx)),
        mc.preschedule = TRUE,
        mc.set.seed = FALSE
      )
    }
    y_block[, column_idx] <- do.call(cbind, filtered_columns)

    finished <- max(column_idx)
    if (verbose && (finished %% 10L == 0L || finished == p)) {
      message(sprintf("Finished FFT series %d of %d", finished, p))
      invisible(gc())
    }
  }

  x_block <- right_solve_upper_banded(y = y_block, r_upper = r_upper, b = b)
  rm(y_block, b_vec, a_padded)
  invisible(gc())
  x_block
}
