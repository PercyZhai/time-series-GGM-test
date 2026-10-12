#!/usr/bin/env Rscript
# Combine completed per-cell summaries after all Slurm jobs have finished.
args <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1L])
root <- normalizePath(file.path(dirname(self), ".."), mustWork = TRUE)
tables <- file.path(root, "Results", "tables")
files <- list.files(tables, pattern = "^summary_T_hat_vs_Bootstrap_.*[.]csv$", full.names = TRUE)
if (!length(files)) stop("No per-cell bootstrap summaries found in ", tables)
rows <- do.call(rbind, lapply(files, function(path) read.csv(path, check.names = FALSE)))
rows <- rows[order(rows$n, rows$p, rows$b, rows$beta), , drop = FALSE]
out <- file.path(tables, "bootstrap_summary.csv")
write.csv(rows, out, row.names = FALSE)
cat("Wrote", nrow(rows), "rows to", out, "\n")
