#!/usr/bin/env Rscript
# Stage 6b — adaptive top-up: bring failing ridges to the 100k-draw cap
# (Austin 2026-07-04: cap 100k; ridges still unstable there are flagged
# "tail not stably calibrated"). Same production pipeline as stage 6.
# Seeds: 20260708 + cell_id (distinct stream from stage 6).
ROOT <- Sys.getenv("GRASS_SIM_ROOT", "/Users/austinsemmel/Desktop/PABAK_Investigation")
setwd(ROOT)
CAND_LIB <- file.path(ROOT, "grassr/simulation/v070_program/tier1/stage3_threshold_grid/lib")
suppressPackageStartupMessages(library(future.apply))
.libPaths(CAND_LIB)
suppressPackageStartupMessages(library(grassr, lib.loc = CAND_LIB))
if (requireNamespace("lme4", quietly = TRUE)) stop("lme4 visible; park it first.")
S6B <- file.path(ROOT, "grassr/simulation/v070_program/tier1/stage6c_highq")
PER_CELL <- file.path(S6B, "per_cell"); dir.create(PER_CELL, showWarnings = FALSE)
PROGRESS <- file.path(S6B, "PROGRESS.txt")
grid <- readRDS(file.path(S6B, "topup_grid.rds"))
mod <- as.integer(Sys.getenv("SPLIT_MOD", "1")); rem <- as.integer(Sys.getenv("SPLIT_REM", "0"))
grid <- grid[grid$cell_id %% mod == rem, , drop = FALSE]
done <- as.integer(sub("^cell_([0-9]+)\\.rds$", "\\1", list.files(PER_CELL)))
todo <- grid[!grid$cell_id %in% done, , drop = FALSE]
cat(sprintf("stage6c | split cells: %d | todo: %d | draws: %s\n",
    nrow(grid), nrow(todo), format(sum(todo$reps), big.mark = ",")))
gen_panel <- function(N, k, prev, q) {
  C <- rbinom(N, 1L, prev); Y <- matrix(0L, N, k)
  for (j in seq_len(k)) Y[, j] <- rbinom(N, 1L, ifelse(C == 1L, q, 1 - q)); Y
}
plan(multicore, workers = max(1L, min(as.integer(Sys.getenv("WORKERS", "10")),
                                      parallel::detectCores() - 2L)))
t0 <- Sys.time()
invisible(future_lapply(seq_len(nrow(todo)), function(i) {
  row <- todo[i, ]; set.seed(20260710L + row$cell_id)
  d <- vapply(seq_len(row$reps), function(r) {
    Y <- gen_panel(row$N, row$k, row$prev, row$q)
    suppressWarnings(suppressMessages(
      tryCatch(check_asymmetry(Y)$delta_hat, error = function(e) NA_real_)))
  }, numeric(1L))
  saveRDS(list(cell = row, delta = d, n_finite = sum(is.finite(d)),
               r_version = R.version.string),
          file.path(PER_CELL, sprintf("cell_%04d.rds", row$cell_id)))
  if (i %% 20L == 0L) cat(sprintf("[%s] %d/%d\n", format(Sys.time()), i, nrow(todo)),
                          file = PROGRESS, append = TRUE)
  NULL
}, future.seed = TRUE))
cat(sprintf("[%s] stage6c split %d/%d COMPLETE. Wall: %.1f min.\n",
    format(Sys.time()), rem, mod, as.numeric(Sys.time() - t0, units = "mins")),
    file = PROGRESS, append = TRUE)
