#!/usr/bin/env Rscript

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
}

script_dir <- dirname(get_script_path())
project_root <- if (basename(script_dir) == "scripts") dirname(script_dir) else script_dir

source(file.path(project_root, "scripts", "io_utils.R"))
source(file.path(project_root, "scripts", "dgp.R"))
source(file.path(project_root, "scripts", "statistics.R"))
source(file.path(project_root, "scripts", "metrics.R"))
source(file.path(project_root, "scripts", "plots.R"))

normalize_flag <- function(flag) {
  switch(
    flag,
    "bandwidth" = "b",
    "mc" = "N_mc",
    "n-mc" = "N_mc",
    "block-length" = "block_length",
    "ell" = "block_length",
    "n-trunc" = "n_trunc",
    "ntrunc-mode" = "ntrunc_mode",
    "trunc-mult" = "trunc_mult",
    "theta0" = "theta0",
    "rho" = "rho",
    "delta0" = "delta0",
    "mc-cores" = "cores",
    "fft-cores" = "fft_cores",
    "progress-every" = "progress_every",
    gsub("-", "_", flag, fixed = TRUE)
  )
}

parse_args <- function() {
  tokens <- commandArgs(trailingOnly = TRUE)
  args <- list(
    quick = FALSE,
    n = NULL,
    p = NULL,
    b = NULL,
    beta = NULL,
    beta_list = NULL,
    N_mc = NULL,
    alpha = 0.10,
    block_length = NULL,
    seed = 20260424L,
    trunc_mult = 10,
    n_trunc = NULL,
    ntrunc_mode = "mult",
    theta0 = 0.25,
    rho = 0.5,
    delta0 = 0.10,
    center = FALSE,
    overwrite = FALSE,
    cores = 1L,
    fft_cores = 1L,
    progress_every = NULL
  )

  boolean_flags <- c("quick", "center", "overwrite")
  multi_flags <- c("beta-list")
  scalar_flags <- c(
    "n", "p", "b", "bandwidth", "beta", "mc", "n-mc", "alpha",
    "block-length", "ell", "seed", "trunc-mult", "n-trunc",
    "ntrunc-mode", "theta0", "rho", "delta0", "cores", "mc-cores", "fft-cores",
    "progress-every"
  )

  i <- 1L
  while (i <= length(tokens)) {
    token <- tokens[i]
    if (!startsWith(token, "--")) {
      stop(sprintf("Unexpected argument: %s", token), call. = FALSE)
    }

    flag <- sub("^--", "", token)
    key <- normalize_flag(flag)

    if (flag %in% boolean_flags) {
      args[[key]] <- TRUE
      i <- i + 1L
      next
    }

    if (!(flag %in% c(multi_flags, scalar_flags))) {
      stop(sprintf("Unknown flag: --%s", flag), call. = FALSE)
    }

    j <- i + 1L
    values <- character(0)
    while (j <= length(tokens) && !startsWith(tokens[j], "--")) {
      values <- c(values, tokens[j])
      j <- j + 1L
    }
    if (length(values) == 0L) {
      stop(sprintf("Flag --%s requires a value.", flag), call. = FALSE)
    }

    if (flag %in% multi_flags) {
      args[[key]] <- values
    } else {
      if (length(values) != 1L) {
        stop(sprintf("Flag --%s expects a single value.", flag), call. = FALSE)
      }
      args[[key]] <- values[[1]]
    }
    i <- j
  }

  if (!is.null(args$n)) args$n <- as.integer(args$n)
  if (!is.null(args$p)) args$p <- as.integer(args$p)
  if (!is.null(args$b)) args$b <- as.integer(args$b)
  if (!is.null(args$beta)) args$beta <- as.numeric(args$beta)
  if (!is.null(args$beta_list)) args$beta_list <- as.numeric(args$beta_list)
  if (!is.null(args$N_mc)) args$N_mc <- as.integer(args$N_mc)
  args$alpha <- as.numeric(args$alpha)
  if (!is.null(args$block_length)) args$block_length <- as.integer(args$block_length)
  args$seed <- as.integer(args$seed)
  args$trunc_mult <- as.numeric(args$trunc_mult)
  if (!is.null(args$n_trunc)) args$n_trunc <- as.integer(args$n_trunc)
  args$theta0 <- as.numeric(args$theta0)
  args$rho <- as.numeric(args$rho)
  args$delta0 <- as.numeric(args$delta0)
  args$cores <- as.integer(args$cores)
  args$fft_cores <- as.integer(args$fft_cores)
  if (!is.null(args$progress_every)) args$progress_every <- as.integer(args$progress_every)
  args
}

resolve_configs <- function(args) {
  n <- if (!is.null(args$n)) args$n else if (isTRUE(args$quick)) 120L else 300L
  p <- if (!is.null(args$p)) args$p else if (isTRUE(args$quick)) 40L else 100L
  b <- if (!is.null(args$b)) args$b else 2L
  N_mc <- if (!is.null(args$N_mc)) args$N_mc else if (isTRUE(args$quick)) 20L else 200L
  block_length <- if (!is.null(args$block_length)) args$block_length else floor(n^(2 / 3))
  beta_values <- if (!is.null(args$beta_list)) {
    args$beta_list
  } else if (!is.null(args$beta)) {
    args$beta
  } else {
    0.9
  }

  if (!args$ntrunc_mode %in% c("mult", "n2")) {
    stop("--ntrunc-mode must be one of mult or n2.", call. = FALSE)
  }
  if (!is.null(args$n_trunc)) {
    N_trunc <- args$n_trunc
  } else if (identical(args$ntrunc_mode, "n2")) {
    N_trunc <- n * n
  } else {
    N_trunc <- ceiling(args$trunc_mult * n)
  }

  if (n <= 1L || p <= 1L) {
    stop("n and p must both be greater than 1.", call. = FALSE)
  }
  if (b < 0L || b >= p) {
    stop("b must satisfy 0 <= b < p.", call. = FALSE)
  }
  if (N_mc <= 0L) {
    stop("N_mc must be positive.", call. = FALSE)
  }
  if (args$cores <= 0L) {
    stop("cores must be positive.", call. = FALSE)
  }
  if (args$fft_cores <= 0L) {
    stop("fft_cores must be positive.", call. = FALSE)
  }
  if (!(args$alpha > 0 && args$alpha < 1)) {
    stop("alpha must lie in (0, 1).", call. = FALSE)
  }
  if (block_length < 1L || block_length > n) {
    stop("block_length must be between 1 and n.", call. = FALSE)
  }
  if (N_trunc <= 0L) {
    stop("N_trunc must be positive.", call. = FALSE)
  }

  lapply(beta_values, function(beta) {
    list(
      n = as.integer(n),
      p = as.integer(p),
      b = as.integer(b),
      beta = as.numeric(beta),
      N_mc = as.integer(N_mc),
      alpha = as.numeric(args$alpha),
      block_length = as.integer(block_length),
      N_trunc = as.integer(N_trunc),
      theta0 = as.numeric(args$theta0),
      rho = as.numeric(args$rho),
      delta0 = as.numeric(args$delta0),
      center = isTRUE(args$center),
      trunc_mult = as.numeric(args$trunc_mult),
      ntrunc_mode = args$ntrunc_mode
    )
  })
}

make_config_seed <- function(base_seed, config) {
  beta_part <- as.integer(round(config$beta * 1e6))
  raw_seed <- as.numeric(base_seed) +
    1009 * as.numeric(config$n) +
    9176 * as.numeric(config$p) +
    1291 * as.numeric(config$b) +
    as.numeric(beta_part)
  positive_mod(raw_seed, .Machine$integer.max - 1L) + 1L
}

run_config <- function(config, args, results_paths) {
  base_name <- config_tag(config$n, config$p, config$b, config$beta)
  raw_path <- file.path(
    results_paths$raw,
    sprintf("boot_%s_seed%d.rds", base_name, args$seed)
  )
  if (file.exists(raw_path) && !isTRUE(args$overwrite)) {
    message(sprintf("Skipping existing configuration %s; use --overwrite to rerun.", base_name))
    return(NULL)
  }

  message(sprintf(
    "Running config n=%d, p=%d, b=%d, beta=%s, N_mc=%d, l=%d",
    config$n,
    config$p,
    config$b,
    format_beta(config$beta),
    config$N_mc,
    config$block_length
  ))

  config_seed <- make_config_seed(args$seed, config)
  set.seed(config_seed)

  matrices <- build_model_matrices(
    p = config$p,
    b = config$b,
    theta0 = config$theta0,
    rho = config$rho,
    delta0 = config$delta0
  )
  graph <- build_graph_structure(p = config$p, b = config$b)
  a <- make_filter_coefficients(beta = config$beta, n_trunc = config$N_trunc)

  n_blocks <- config$n - config$block_length + 1L
  t_hat_samples <- numeric(config$N_mc)
  q_hat_samples <- numeric(config$N_mc)
  reject_values <- integer(config$N_mc)
  fallback_counts <- integer(config$N_mc)
  pooled_bootstrap <- numeric(config$N_mc * n_blocks)

  progress_every <- args$progress_every
  if (is.null(progress_every)) {
    progress_every <- max(1L, min(10L, config$N_mc %/% 10L))
  }
  progress_every <- max(1L, as.integer(progress_every))

  n_cores <- min(as.integer(args$cores), config$N_mc)
  fft_cores <- min(as.integer(args$fft_cores), config$p)
  if (n_cores > 1L && .Platform$OS.type != "unix") {
    stop("Parallel replication mode requires a Unix-like platform.", call. = FALSE)
  }
  message(sprintf("Replication backend: %d core(s)", n_cores))
  message(sprintf("FFT generation backend: %d core(s)", fft_cores))
  message(sprintf("Statistics batch size: %d replication(s)", progress_every))

  run_replication <- function(rep_idx, x_rep) {
    rep_stats <- compute_replication_statistics(
      x = x_rep,
      graph = graph,
      alpha = config$alpha,
      block_length = config$block_length,
      center = config$center
    )
    list(
      rep_idx = rep_idx,
      T_hat = rep_stats$T_hat,
      q_hat = rep_stats$q_hat,
      reject = rep_stats$reject,
      fallback_count = rep_stats$fallback_count,
      boot_stats = rep_stats$boot_stats
    )
  }

  store_replication <- function(rep_result) {
    rep_idx <- rep_result$rep_idx
    t_hat_samples[rep_idx] <<- rep_result$T_hat
    q_hat_samples[rep_idx] <<- rep_result$q_hat
    reject_values[rep_idx] <<- rep_result$reject
    fallback_counts[rep_idx] <<- rep_result$fallback_count
    boot_start <- (rep_idx - 1L) * n_blocks + 1L
    boot_end <- rep_idx * n_blocks
    pooled_bootstrap[boot_start:boot_end] <<- rep_result$boot_stats
  }

  rep_pb <- txtProgressBar(min = 0, max = config$N_mc, style = 3)
  started <- proc.time()[["elapsed"]]
  rep_indices <- seq_len(config$N_mc)
  chunks <- split(rep_indices, ceiling(seq_along(rep_indices) / progress_every))
  completed <- 0L

  generation_started <- proc.time()[["elapsed"]]
  x_bank <- generate_fft_sample_block(
    n = config$n,
    n_mc = config$N_mc,
    r_upper = matrices$R,
    b = config$b,
    a = a,
    p = config$p,
    verbose = TRUE,
    fft_cores = fft_cores
  )
  message(sprintf(
    "Generated the full %d-replication data bank in %.2f seconds.",
    config$N_mc,
    proc.time()[["elapsed"]] - generation_started
  ))

  if (n_cores == 1L) {
    for (chunk in chunks) {
      for (rep_idx in chunk) {
        row_start <- (rep_idx - 1L) * config$n + 1L
        row_end <- rep_idx * config$n
        rep_result <- run_replication(rep_idx, x_bank[row_start:row_end, , drop = FALSE])
        store_replication(rep_result)
        completed <- completed + 1L
        setTxtProgressBar(rep_pb, completed)
        if (completed %% progress_every == 0L || completed == config$N_mc) {
          message(sprintf(
            "\nCompleted %d/%d replications: latest T_hat=%.4f, latest q=%.4f, running size=%.3f",
            completed,
            config$N_mc,
            rep_result$T_hat,
            rep_result$q_hat,
            mean(reject_values[seq_len(completed)])
          ))
        }
      }
      invisible(gc())
    }
  } else {
    for (chunk in chunks) {
      chunk_results <- parallel::mclapply(
        seq_along(chunk),
        function(pos) {
          rep_idx <- chunk[pos]
          row_start <- (rep_idx - 1L) * config$n + 1L
          row_end <- rep_idx * config$n
          run_replication(
            rep_idx = rep_idx,
            x_rep = x_bank[row_start:row_end, , drop = FALSE]
          )
        },
        mc.cores = min(n_cores, length(chunk)),
        mc.preschedule = FALSE
      )
      for (rep_result in chunk_results) {
      store_replication(rep_result)
      }
      completed <- completed + length(chunk_results)
      setTxtProgressBar(rep_pb, completed)
      latest <- chunk_results[[length(chunk_results)]]
      message(sprintf(
        "\nCompleted %d/%d replications: latest T_hat=%.4f, latest q=%.4f, running size=%.3f",
        completed,
        config$N_mc,
        latest$T_hat,
        latest$q_hat,
        mean(reject_values[seq_len(completed)])
      ))
      rm(chunk_results)
      invisible(gc())
    }
  }
  close(rep_pb)
  elapsed <- proc.time()[["elapsed"]] - started
  rm(x_bank)
  invisible(gc())

  fallback_count <- sum(fallback_counts)
  summary_row <- make_summary_row(
    config = config,
    statistic_samples = t_hat_samples,
    bootstrap_samples = pooled_bootstrap,
    q_values = q_hat_samples,
    reject_values = reject_values,
    fallback_count = fallback_count,
    seed = args$seed
  )
  config_summary_path <- file.path(
    results_paths$tables,
    sprintf("summary_T_hat_vs_Bootstrap_%s.csv", base_name)
  )
  write.csv(
    data.frame(summary_row, stringsAsFactors = FALSE, check.names = FALSE),
    config_summary_path,
    row.names = FALSE
  )
  metrics <- list(
    kolmogorov_distance = summary_row$kolmogorov_distance,
    wasserstein1_distance = summary_row$wasserstein1_distance
  )
  quantile_table <- make_quantile_table(t_hat_samples, pooled_bootstrap)

  raw_payload <- list(
    T_hat_samples = t_hat_samples,
    bootstrap_samples = pooled_bootstrap,
    q_hat_samples = q_hat_samples,
    reject_values = reject_values,
    fallback_counts = fallback_counts,
    n = config$n,
    p = config$p,
    b = config$b,
    beta = config$beta,
    N_mc = config$N_mc,
    N_trunc = config$N_trunc,
    alpha = config$alpha,
    block_length = config$block_length,
    N_boot = length(pooled_bootstrap),
    seed = args$seed,
    config_seed = config_seed,
    data_bank_size = config$N_mc,
    fft_cores = fft_cores,
    elapsed_seconds = elapsed
  )
  save_rds_file(raw_path, raw_payload, overwrite = args$overwrite)

  quantile_path <- file.path(
    results_paths$tables,
    sprintf("quantiles_T_hat_vs_Bootstrap_%s.csv", base_name)
  )
  write.csv(quantile_table, quantile_path, row.names = FALSE)

  qq_pdf <- file.path(results_paths$figures, sprintf("qq_T_hat_vs_Bootstrap_%s.pdf", base_name))

  save_qq_plot(
    statistic_samples = t_hat_samples,
    target_samples = pooled_bootstrap,
    statistic_name = "T_hat",
    target_name = "Bootstrap",
    metrics = metrics,
    n = config$n,
    p = config$p,
    b = config$b,
    beta = config$beta,
    N_mc = config$N_mc,
    N_boot = length(pooled_bootstrap),
    block_length = config$block_length,
    output_path = qq_pdf
  )

  list(
    base_name = base_name,
    raw_path = raw_path,
    config_summary_path = config_summary_path,
    quantile_path = quantile_path,
    qq_pdf = qq_pdf,
    summary_row = summary_row,
    elapsed_seconds = elapsed
  )
}

main <- function() {
  args <- parse_args()
  configs <- resolve_configs(args)
  results_paths <- ensure_results_tree(project_root)
  timestamp <- timestamp_string()
  started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

  config_log <- list(
    timestamp = timestamp,
    project_root = project_root,
    cli_args = args,
    resolved_configs = configs
  )
  write_json_file(
    file.path(results_paths$logs, sprintf("config_used_%s.json", timestamp)),
    config_log
  )

  summary_rows <- list()
  run_log_lines <- c(
    sprintf("Started: %s", started_at),
    sprintf("Seed: %d", args$seed),
    sprintf("Overwrite: %s", tolower(as.character(args$overwrite))),
    sprintf("Number of configurations: %d", length(configs)),
    ""
  )

  total_fallback_count <- 0L
  completed_configs <- 0L
  skipped_configs <- 0L

  for (config in configs) {
    result <- run_config(config, args, results_paths)
    if (is.null(result)) {
      skipped_configs <- skipped_configs + 1L
      next
    }

    summary_rows[[length(summary_rows) + 1L]] <- result$summary_row
    total_fallback_count <- total_fallback_count + result$summary_row$fallback_count
    completed_configs <- completed_configs + 1L
    run_log_lines <- c(
      run_log_lines,
      sprintf("Completed %s", result$base_name),
      sprintf("  Raw: %s", result$raw_path),
      sprintf("  Summary: %s", result$config_summary_path),
      sprintf("  Quantiles: %s", result$quantile_path),
      sprintf("  QQ PDF: %s", result$qq_pdf),
      sprintf(
        "  empirical_size=%.4f, KS=%.4f, W1=%.4f, fallback_count=%d, elapsed=%.1fs",
        result$summary_row$empirical_size,
        result$summary_row$kolmogorov_distance,
        result$summary_row$wasserstein1_distance,
        result$summary_row$fallback_count,
        result$elapsed_seconds
      ),
      ""
    )
  }

  summary_path <- file.path(results_paths$tables, "bootstrap_summary.csv")
  if (length(summary_rows) > 0L) {
    append_rows_csv(
      csv_path = summary_path,
      rows = summary_rows,
      key_columns = c("n", "p", "b", "beta", "N_mc", "N_trunc", "alpha", "block_length", "seed")
    )
  }

  ended_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  run_log_lines <- c(
    run_log_lines,
    sprintf("Finished: %s", ended_at),
    sprintf("Completed configurations: %d", completed_configs),
    sprintf("Skipped configurations: %d", skipped_configs),
    sprintf("Total fallback count: %d", total_fallback_count),
    sprintf("Summary CSV: %s", summary_path)
  )
  run_log_path <- file.path(results_paths$logs, sprintf("run_log_%s.txt", timestamp))
  write_text_file(run_log_path, paste(run_log_lines, collapse = "\n"))

  if (length(summary_rows) > 0L) {
    message("\nSummary metrics:")
    for (row in summary_rows) {
      message(sprintf(
        "n=%d, p=%d, b=%d, beta=%s: empirical_size=%.4f, KS=%.4f, W1=%.4f, fallback=%d",
        row$n,
        row$p,
        row$b,
        format_beta(row$beta),
        row$empirical_size,
        row$kolmogorov_distance,
        row$wasserstein1_distance,
        row$fallback_count
      ))
    }
  }
  message(sprintf("Run log: %s", run_log_path))
}

main()
