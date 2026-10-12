ensure_results_tree <- function(project_root) {
  results_dir <- file.path(project_root, "Results")
  paths <- list(
    results = results_dir,
    figures = file.path(results_dir, "figures"),
    tables = file.path(results_dir, "tables"),
    raw = file.path(results_dir, "raw"),
    logs = file.path(results_dir, "logs")
  )
  for (path in paths) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  paths
}

timestamp_string <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

positive_mod <- function(x, modulus) {
  as.integer((((as.numeric(x) %% modulus) + modulus) %% modulus))
}

as_integer_vector <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  as.integer(x)
}

as_numeric_vector <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  as.numeric(x)
}

as_nullable_integer <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  as.integer(x)
}

format_beta <- function(beta) {
  out <- format(beta, scientific = FALSE, trim = TRUE, digits = 15)
  sub("\\.?0+$", "", out)
}

config_tag <- function(n, p, b, beta) {
  sprintf("n%d_p%d_b%d_beta%s", n, p, b, format_beta(beta))
}

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  x <- gsub("\b", "\\\\b", x)
  x <- gsub("\f", "\\\\f", x)
  x <- gsub("\n", "\\\\n", x)
  x <- gsub("\r", "\\\\r", x)
  gsub("\t", "\\\\t", x)
}

json_scalar <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) {
    return("null")
  }
  if (is.logical(x)) {
    return(if (isTRUE(x)) "true" else "false")
  }
  if (is.numeric(x)) {
    if (!is.finite(x)) {
      return("null")
    }
    return(format(x, scientific = FALSE, trim = TRUE, digits = 17))
  }
  sprintf("\"%s\"", json_escape(as.character(x)))
}

to_json_string <- function(x, indent = 0L) {
  pad <- paste(rep(" ", indent), collapse = "")
  child_pad <- paste(rep(" ", indent + 2L), collapse = "")

  if (is.null(x) || length(x) == 0L) {
    return("null")
  }

  if (is.atomic(x) && length(x) == 1L) {
    return(json_scalar(x))
  }

  if (is.atomic(x) && length(x) > 1L) {
    items <- vapply(
      as.list(x),
      function(item) paste0(child_pad, to_json_string(item, indent + 2L)),
      character(1)
    )
    return(sprintf("[\n%s\n%s]", paste(items, collapse = ",\n"), pad))
  }

  if (is.list(x)) {
    item_names <- names(x)
    if (!is.null(item_names) && all(nzchar(item_names))) {
      pieces <- vapply(
        seq_along(x),
        function(i) {
          sprintf(
            "%s\"%s\": %s",
            child_pad,
            json_escape(item_names[i]),
            to_json_string(x[[i]], indent + 2L)
          )
        },
        character(1)
      )
      return(sprintf("{\n%s\n%s}", paste(pieces, collapse = ",\n"), pad))
    }
    pieces <- vapply(
      x,
      function(item) paste0(child_pad, to_json_string(item, indent + 2L)),
      character(1)
    )
    return(sprintf("[\n%s\n%s]", paste(pieces, collapse = ",\n"), pad))
  }

  json_scalar(as.character(x))
}

write_json_file <- function(path, payload) {
  writeLines(to_json_string(payload), con = path, useBytes = TRUE)
}

write_text_file <- function(path, text) {
  writeLines(text, con = path, useBytes = TRUE)
}

save_rds_file <- function(path, object, overwrite = FALSE) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(sprintf("Refusing to overwrite existing file: %s", path), call. = FALSE)
  }
  saveRDS(object, path)
}

append_rows_csv <- function(csv_path, rows, key_columns = NULL) {
  new_df <- do.call(rbind, lapply(rows, function(row) {
    data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  }))

  if (file.exists(csv_path)) {
    old_df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
    combined <- rbind(old_df, new_df)
  } else {
    combined <- new_df
  }

  if (!is.null(key_columns)) {
    key <- do.call(paste, c(combined[, key_columns, drop = FALSE], sep = "|"))
    combined <- combined[!duplicated(key, fromLast = TRUE), , drop = FALSE]
  }

  write.csv(combined, csv_path, row.names = FALSE)
}
