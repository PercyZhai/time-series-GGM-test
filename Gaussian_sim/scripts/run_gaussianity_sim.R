#!/usr/bin/env Rscript

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
}

script_dir <- dirname(get_script_path())
project_root <- if (basename(script_dir) == "scripts") dirname(script_dir) else script_dir

source(file.path(project_root, "scripts", "io_utils.R"))
source(file.path(project_root, "scripts", "dgp.R"))
source(file.path(project_root, "scripts", "statistics.R"))
source(file.path(project_root, "scripts", "gaussian_target.R"))
source(file.path(project_root, "scripts", "metrics.R"))
source(file.path(project_root, "scripts", "plots.R"))

parse_args <- function() {
  tokens <- commandArgs(trailingOnly = TRUE)
  args <- list(
    preset = "smoke",
    n_list = NULL,
    p_list = NULL,
    beta_list = NULL,
    b_list = NULL,
    n_mc = NULL,
    b_gauss = NULL,
    n_trunc = NULL,
    seed = 123L,
    theta0 = 0.25,
    rho = 0.5,
    delta0 = 0.10,
    cores = 1L,
    fft_cores = 1L,
    use_asymptotic_g = FALSE,
    overwrite = FALSE
  )

  multi_flags <- c("n-list", "p-list", "beta-list", "b-list")
  scalar_flags <- c("preset", "n-mc", "b-gauss", "n-trunc", "seed", "theta0", "rho", "delta0", "cores", "fft-cores")
  boolean_flags <- c("use-asymptotic-g", "overwrite")

  i <- 1L
  while (i <= length(tokens)) {
    token <- tokens[i]
    if (!startsWith(token, "--")) {
      stop(sprintf("Unexpected argument: %s", token), call. = FALSE)
    }

    flag <- sub("^--", "", token)
    key <- gsub("-", "_", flag, fixed = TRUE)

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
    if (length(values) == 0) {
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

  args$n_list <- as_integer_vector(args$n_list)
  args$p_list <- as_integer_vector(args$p_list)
  args$beta_list <- as_numeric_vector(args$beta_list)
  args$b_list <- as_integer_vector(args$b_list)
  args$n_mc <- as_nullable_integer(args$n_mc)
  args$b_gauss <- as_nullable_integer(args$b_gauss)
  args$n_trunc <- as_nullable_integer(args$n_trunc)
  args$seed <- as.integer(args$seed)
  args$theta0 <- as.numeric(args$theta0)
  args$rho <- as.numeric(args$rho)
  args$delta0 <- as.numeric(args$delta0)
  args$cores <- as.integer(args$cores)
  args$fft_cores <- as.integer(args$fft_cores)
  if (args$cores <= 0L || args$fft_cores <= 0L) {
    stop("cores and fft-cores must be positive.", call. = FALSE)
  }
  args
}

resolve_grid <- function(args) {
  smoke_defaults <- list(
    n_list = c(100L),
    p_list = c(20L),
    beta_list = c(2.0, 0.9, 0.55),
    b_list = c(2L),
    n_mc = 20L,
    b_gauss = 200L
  )
  paper_defaults <- list(
    n_list = c(200L, 500L, 1000L),
    p_list = c(20L, 50L, 100L),
    beta_list = c(2.0, 0.9, 0.55),
    b_list = c(2L, 3L),
    n_mc = 200L,
    b_gauss = 2000L
  )

  if (!(args$preset %in% c("smoke", "paper", "custom"))) {
    stop("`--preset` must be one of smoke, paper, or custom.", call. = FALSE)
  }

  if (identical(args$preset, "smoke")) {
    defaults <- smoke_defaults
  } else if (identical(args$preset, "paper")) {
    defaults <- paper_defaults
  } else {
    defaults <- list(
      n_list = args$n_list,
      p_list = args$p_list,
      beta_list = args$beta_list,
      b_list = args$b_list,
      n_mc = args$n_mc,
      b_gauss = args$b_gauss
    )
    missing_keys <- names(defaults)[vapply(defaults, is.null, logical(1))]
    if (length(missing_keys) > 0) {
      missing_flags <- paste0("--", gsub("_", "-", missing_keys))
      stop(
        sprintf("Custom preset requires: %s", paste(missing_flags, collapse = ", ")),
        call. = FALSE
      )
    }
  }

  n_list <- if (!is.null(args$n_list)) args$n_list else defaults$n_list
  p_list <- if (!is.null(args$p_list)) args$p_list else defaults$p_list
  beta_list <- if (!is.null(args$beta_list)) args$beta_list else defaults$beta_list
  b_list <- if (!is.null(args$b_list)) args$b_list else defaults$b_list
  n_mc <- if (!is.null(args$n_mc)) args$n_mc else defaults$n_mc
  b_gauss <- if (!is.null(args$b_gauss)) args$b_gauss else defaults$b_gauss

  configs <- list()
  idx <- 1L
  for (n in n_list) {
    if (!is.null(args$n_trunc)) {
      n_trunc <- args$n_trunc
    } else if (identical(args$preset, "smoke")) {
      n_trunc <- min(n * n, 20000L)
    } else {
      n_trunc <- n * n
    }
    for (p in p_list) {
      for (b in b_list) {
        for (beta in beta_list) {
          configs[[idx]] <- list(
            n = as.integer(n),
            p = as.integer(p),
            b = as.integer(b),
            beta = as.numeric(beta),
            N_mc = as.integer(n_mc),
            B_gauss = as.integer(b_gauss),
            N_trunc = as.integer(n_trunc),
            theta0 = as.numeric(args$theta0),
            rho = as.numeric(args$rho),
            delta0 = as.numeric(args$delta0)
          )
          idx <- idx + 1L
        }
      }
    }
  }
  configs
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

main <- function() {
  args <- parse_args()
  results_paths <- ensure_results_tree(project_root)
  timestamp <- timestamp_string()
  started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  configs <- resolve_grid(args)

  config_log <- list(
    timestamp = timestamp,
    project_root = project_root,
    results_paths = as.character_list(results_paths),
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
    sprintf("Preset: %s", args$preset),
    sprintf("Seed: %d", args$seed),
    sprintf("Overwrite: %s", tolower(as.character(args$overwrite))),
    sprintf("Use asymptotic G: %s", tolower(as.character(args$use_asymptotic_g))),
    sprintf("Number of configurations: %d", length(configs)),
    ""
  )

  total_fallback_count <- 0L
  completed_configs <- 0L
  skipped_configs <- 0L

  if (length(configs) > 0) {
    config_pb <- txtProgressBar(min = 0, max = length(configs), style = 3)
  }

  for (config_idx in seq_along(configs)) {
    config <- configs[[config_idx]]
    setTxtProgressBar(config_pb, config_idx - 1L)

    raw_path <- raw_sample_path(
      results_paths$raw,
      n = config$n,
      p = config$p,
      b = config$b,
      beta = config$beta,
      seed = args$seed
    )

    if (file.exists(raw_path) && !isTRUE(args$overwrite)) {
      skipped_configs <- skipped_configs + 1L
      run_log_lines <- c(
        run_log_lines,
        sprintf("Skipped existing configuration: %s", config_tag(config$n, config$p, config$b, config$beta))
      )
      next
    }

    message(sprintf(
      "Running config n=%d, p=%d, b=%d, beta=%s",
      config$n,
      config$p,
      config$b,
      format_beta(config$beta)
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
    graph <- build_graph_structure(
      sigma = matrices$Sigma,
      omega = matrices$Omega,
      b = config$b
    )
    a <- make_filter_coefficients(beta = config$beta, n_trunc = config$N_trunc)
    gamma <- compute_gamma_sequence_fft(a)
    l_factors <- compute_l_factors(gamma = gamma, n = config$n)

    x_block <- generate_fft_sample_block(
      n = config$n,
      n_mc = config$N_mc,
      r_upper = matrices$R,
      b = config$b,
      a = a,
      p = config$p,
      verbose = TRUE,
      fft_cores = args$fft_cores
    )

    t_hat_samples <- numeric(config$N_mc)
    t_l_samples <- numeric(config$N_mc)
    fallback_count <- 0L

    rep_pb <- txtProgressBar(min = 0, max = config$N_mc, style = 3)
    rep_indices <- seq_len(config$N_mc)
    rep_batches <- split(rep_indices, ceiling(rep_indices / args$cores))
    completed <- 0L
    for (rep_idx in rep_batches) {
      rep_results <- parallel::mclapply(
        rep_idx,
        function(idx) {
          row_start <- (idx - 1L) * config$n + 1L
          row_end <- idx * config$n
          compute_replication_statistics(
            x = x_block[row_start:row_end, , drop = FALSE],
            sigma = matrices$Sigma,
            graph = graph
          )
        },
        mc.cores = min(args$cores, length(rep_idx)),
        mc.preschedule = FALSE,
        mc.set.seed = FALSE
      )
      t_hat_samples[rep_idx] <- vapply(rep_results, `[[`, numeric(1), "T_hat")
      t_l_samples[rep_idx] <- vapply(rep_results, `[[`, numeric(1), "T_L")
      fallback_count <- fallback_count + sum(vapply(rep_results, `[[`, numeric(1), "fallback_count"))
      completed <- completed + length(rep_idx)
      setTxtProgressBar(rep_pb, completed)
    }
    close(rep_pb)
    rm(x_block)
    invisible(gc())

    set.seed(positive_mod(config_seed + 1000003L, .Machine$integer.max - 1L) + 1L)
    gaussian_samples <- simulate_gaussian_targets(
      b_gauss = config$B_gauss,
      l_finite = l_factors$L_n,
      l_asymptotic = l_factors$L,
      r_upper = matrices$R,
      graph = graph,
      use_asymptotic_g = args$use_asymptotic_g,
      cores = args$cores
    )

    raw_payload <- list(
      T_hat_samples = t_hat_samples,
      T_L_samples = t_l_samples,
      G_finite_samples = gaussian_samples$G_finite_samples,
      G_asymptotic_samples = gaussian_samples$G_asymptotic_samples,
      n = config$n,
      p = config$p,
      b = config$b,
      beta = config$beta,
      N_mc = config$N_mc,
      B_gauss = config$B_gauss,
      N_trunc = config$N_trunc,
      L_n = l_factors$L_n,
      L = l_factors$L,
      fallback_count = fallback_count
    )
    save_rds_file(raw_path, raw_payload, overwrite = args$overwrite)

    comparisons <- list(
      list(statistic = "T_hat", target = "G_finite", statistic_samples = t_hat_samples, target_samples = gaussian_samples$G_finite_samples),
      list(statistic = "T_L", target = "G_finite", statistic_samples = t_l_samples, target_samples = gaussian_samples$G_finite_samples)
    )
    if (!is.null(gaussian_samples$G_asymptotic_samples)) {
      comparisons <- c(
        comparisons,
        list(
          list(statistic = "T_hat", target = "G_asymptotic", statistic_samples = t_hat_samples, target_samples = gaussian_samples$G_asymptotic_samples),
          list(statistic = "T_L", target = "G_asymptotic", statistic_samples = t_l_samples, target_samples = gaussian_samples$G_asymptotic_samples)
        )
      )
    }

    base_name <- config_tag(config$n, config$p, config$b, config$beta)
    for (comparison in comparisons) {
      summary_row <- summarize_comparison(
        n = config$n,
        p = config$p,
        b = config$b,
        beta = config$beta,
        N_mc = config$N_mc,
        B_gauss = config$B_gauss,
        N_trunc = config$N_trunc,
        L_n = l_factors$L_n,
        L = l_factors$L,
        statistic = comparison$statistic,
        target = comparison$target,
        statistic_samples = comparison$statistic_samples,
        target_samples = comparison$target_samples,
        fallback_count = fallback_count,
        seed = args$seed
      )

      metrics <- list(
        kolmogorov_distance = summary_row$kolmogorov_distance,
        wasserstein1_distance = summary_row$wasserstein1_distance
      )

      qq_path <- file.path(
        results_paths$figures,
        sprintf("qq_%s_vs_%s_%s.pdf", comparison$statistic, comparison$target, base_name)
      )
      cdf_path <- file.path(
        results_paths$figures,
        sprintf("cdf_%s_vs_%s_%s.pdf", comparison$statistic, comparison$target, base_name)
      )

      save_qq_plot(
        statistic_samples = comparison$statistic_samples,
        target_samples = comparison$target_samples,
        statistic_name = comparison$statistic,
        target_name = comparison$target,
        metrics = metrics,
        n = config$n,
        p = config$p,
        b = config$b,
        beta = config$beta,
        N_mc = config$N_mc,
        B_gauss = config$B_gauss,
        output_path = qq_path
      )
      save_cdf_plot(
        statistic_samples = comparison$statistic_samples,
        target_samples = comparison$target_samples,
        statistic_name = comparison$statistic,
        target_name = comparison$target,
        metrics = metrics,
        n = config$n,
        p = config$p,
        b = config$b,
        beta = config$beta,
        N_mc = config$N_mc,
        B_gauss = config$B_gauss,
        output_path = cdf_path
      )
      summary_rows[[length(summary_rows) + 1L]] <- summary_row
    }

    total_fallback_count <- total_fallback_count + fallback_count
    completed_configs <- completed_configs + 1L
    run_log_lines <- c(
      run_log_lines,
      sprintf("Completed %s", base_name),
      sprintf("  Raw: %s", raw_path),
      sprintf("  L_n=%.8f, L=%.8f", l_factors$L_n, l_factors$L),
      sprintf("  fallback_count=%d", fallback_count)
    )
  }

  if (length(configs) > 0) {
    setTxtProgressBar(config_pb, length(configs))
    close(config_pb)
  }

  if (length(summary_rows) > 0) {
    summary_path <- file.path(results_paths$tables, "gaussianity_summary.csv")
    append_summary_rows(summary_path, summary_rows)
  }

  ended_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  run_log_lines <- c(
    run_log_lines,
    "",
    sprintf("Finished: %s", ended_at),
    sprintf("Completed configurations: %d", completed_configs),
    sprintf("Skipped configurations: %d", skipped_configs),
    sprintf("Total fallback count: %d", total_fallback_count),
    sprintf("Summary CSV: %s", file.path(results_paths$tables, "gaussianity_summary.csv"))
  )
  write_text_file(
    file.path(results_paths$logs, sprintf("run_log_%s.txt", timestamp)),
    paste(run_log_lines, collapse = "\n")
  )
}

main()
