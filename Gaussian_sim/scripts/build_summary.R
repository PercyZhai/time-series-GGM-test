#!/usr/bin/env Rscript
# Combine completed per-cell summaries after all Slurm jobs have finished.
args <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1L])
root <- normalizePath(file.path(dirname(self), ".."), mustWork = TRUE)
tables <- file.path(root, "Results", "tables")
files <- list.files(tables, pattern = "^summary_n[0-9]+_p[0-9]+_b[0-9]+_beta.*[.]csv$", full.names = TRUE)
if (!length(files)) stop("No per-cell Gaussian summaries found in ", tables)
rows <- do.call(rbind, lapply(files, function(path) read.csv(path, check.names = FALSE)))
rows <- rows[order(rows$n, rows$p, rows$b, rows$beta, rows$statistic, rows$target), , drop = FALSE]
out <- file.path(tables, "gaussianity_summary.csv")
write.csv(rows, out, row.names = FALSE)
cat("Wrote", nrow(rows), "rows to", out, "\n")
