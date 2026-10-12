build_banded_precision <- function(p, b, theta0 = 0.25, rho = 0.5, delta0 = 0.10) {
  omega <- matrix(0, nrow = p, ncol = p)
  if (b > 0L) {
    for (distance in seq_len(b)) {
      value <- theta0 * rho^(distance - 1L)
      rows <- seq_len(p - distance)
      cols <- rows + distance
      omega[cbind(rows, cols)] <- value
      omega[cbind(cols, rows)] <- value
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

build_model_matrices <- function(p, b, theta0 = 0.25, rho = 0.5, delta0 = 0.10) {
  omega <- build_banded_precision(p = p, b = b, theta0 = theta0, rho = rho, delta0 = delta0)
  r_upper <- build_banded_cholesky_upper(omega = omega, b = b)
  rm(omega)
  invisible(gc())
  list(R = r_upper, sparse = FALSE)
}

make_filter_coefficients <- function(beta, n_trunc) {
  lags <- 0:(n_trunc - 1L)
  a_tilde <- (lags + 1)^(-beta)
  a_tilde / sqrt(sum(a_tilde^2))
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

  a_padded <- c(a, rep(0, fft_length - n_trunc))
  b_vec <- fft(a_padded)
  y_block <- matrix(0, nrow = rows_needed, ncol = p)

  filter_one_coordinate <- function(eps_j) {
    eps_padded <- c(eps_j, rep(0, fft_length - sequence_length))
    f_j <- fft(eps_padded)
    filtered_j <- fft(b_vec * f_j, inverse = TRUE) / fft_length
    Re(filtered_j)[n_trunc:(n_trunc + rows_needed - 1L)]
  }

  column_batches <- split(
    seq_len(p),
    ceiling(seq_len(p) / fft_cores)
  )

  for (column_idx in column_batches) {
    # Generate innovations in the parent and in the original column order.
    # Workers perform no RNG calls, so the result is identical to the serial
    # generator for the same seed while the coordinate FFTs run in parallel.
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
    if (verbose && (finished %% 25L == 0L || finished == p)) {
      message(sprintf("Finished FFT series %d of %d", finished, p))
      invisible(gc())
    }
  }

  x_block <- right_solve_upper_banded(y = y_block, r_upper = r_upper, b = b)
  rm(y_block, b_vec, a_padded)
  invisible(gc())
  x_block
}
