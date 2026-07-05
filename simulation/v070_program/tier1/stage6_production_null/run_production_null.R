#!/usr/bin/env Rscript
# Stage 6 — PRODUCTION-FAITHFUL null program for delta-hat (v0.7.0).
#
# Ratified 2026-07-04: delta-hat is reported as its percentile on the
# matched null ECDF; the shipped null must therefore be generated through
# the production pipeline itself. Each draw is check_asymmetry(Y)$delta_hat
# — production normalization, constant-rater dropping, clamp-guard, and
# the 3-coefficient delta family, exactly as a user's card computes it.
#
# lme4 is deliberately HIDDEN from .libPaths: compute_panel() then returns
# icc = NA (which delta never uses) and skips the per-rep glmer, cutting
# rep cost up to 78x. Requires the v0.7.0 NA-skip fix in check_asymmetry
# (candidate lib rebuilt 2026-07-04).
#
# Budget (Austin 2026-07-04, small-N weighted):
#   N in {15,20,30,50,75,100}: 10,000 reps per (prev,k,N,q)  [50k pooled]
#   N in {150,200,300}:         4,000                        [20k pooled]
#   N in {500,1000}:            2,000                        [10k pooled]
# Grid: prev {.05,.20,.50,.80,.95} x k {2,3,5,6,8,10,15,25} x 11 N x 5 q
# = 2,200 cells, ~15.2M production draws. Null only (A = 0); the power
# sweep re-runs in phase 2. Seeds: 20260707 + cell_id.

ROOT <- Sys.getenv("GRASS_SIM_ROOT",
                   "/Users/austinsemmel/Desktop/PABAK_Investigation")
setwd(ROOT)
CAND_LIB <- file.path(ROOT,
  "grassr/simulation/v070_program/tier1/stage3_threshold_grid/lib")
# Load the parallel machinery BEFORE hiding the user library (loaded
# namespaces persist across the .libPaths change); then restrict paths so
# lme4 is invisible and compute_panel()'s per-rep glmer is skipped.
suppressPackageStartupMessages(library(future.apply))
.libPaths(CAND_LIB)
suppressPackageStartupMessages(library(grassr, lib.loc = CAND_LIB))
if (requireNamespace("lme4", quietly = TRUE))
  stop("lme4 is still visible (base/site library); per-rep glmer would ",
       "run. Relocate lme4 or adjust .libPaths before launching stage 6.")

STAGE6   <- file.path(ROOT, "grassr/simulation/v070_program/tier1/stage6_production_null")
PER_CELL <- file.path(STAGE6, "per_cell")
PROGRESS <- file.path(STAGE6, "PROGRESS.txt")
dir.create(PER_CELL, recursive = TRUE, showWarnings = FALSE)

MASTER_SEED <- 20260707L
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

reps_for <- function(N) if (N <= 100L) 10000L else if (N <= 300L) 4000L else 2000L

grid <- expand.grid(prev = c(0.05, 0.20, 0.50, 0.80, 0.95),
                    k = c(2L, 3L, 5L, 6L, 8L, 10L, 15L, 25L),
                    N = c(15L, 20L, 30L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L),
                    q = c(0.65, 0.75, 0.85, 0.92, 0.97),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$reps <- vapply(grid$N, reps_for, integer(1L))
grid$cell_id <- seq_len(nrow(grid))
if (SMOKE) { grid <- grid[c(1L, nrow(grid) %/% 2L, nrow(grid)), ]; grid$reps <- 30L }

mod <- as.integer(Sys.getenv("SPLIT_MOD", "1"))
rem <- as.integer(Sys.getenv("SPLIT_REM", "0"))
grid <- grid[grid$cell_id %% mod == rem, , drop = FALSE]
done <- as.integer(sub("^cell_([0-9]+)\\.rds$", "\\1",
                       list.files(PER_CELL, pattern = "^cell_[0-9]+\\.rds$")))
todo <- grid[!grid$cell_id %in% done, , drop = FALSE]
cat(sprintf("stage6 | cells(split): %d | done: %d | todo: %d | draws(split): %s\n",
            nrow(grid), sum(grid$cell_id %in% done), nrow(todo),
            format(sum(todo$reps), big.mark = ",")))

gen_panel <- function(N, k, prev, q) {
  C <- rbinom(N, 1L, prev)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) Y[, j] <- rbinom(N, 1L, ifelse(C == 1L, q, 1 - q))
  Y
}

run_cell <- function(row) {
  set.seed(MASTER_SEED + row$cell_id)
  d <- vapply(seq_len(row$reps), function(r) {
    Y <- gen_panel(row$N, row$k, row$prev, row$q)
    suppressWarnings(suppressMessages(
      tryCatch(check_asymmetry(Y)$delta_hat, error = function(e) NA_real_)))
  }, numeric(1L))
  list(cell = row, delta = d,
       n_finite = sum(is.finite(d)), r_version = R.version.string)
}

n_workers <- max(1L, min(as.integer(Sys.getenv("WORKERS", "10")),
                         parallel::detectCores() - 2L))
plan(multicore, workers = n_workers)
cat(sprintf("[%s] stage6 start: %d todo, %d workers, split %d/%d, %s, lme4 hidden\n",
            format(Sys.time()), nrow(todo), n_workers, rem, mod, R.version.string),
    file = PROGRESS, append = file.exists(PROGRESS))

t0 <- Sys.time()
invisible(future_lapply(seq_len(nrow(todo)), function(i) {
  row <- todo[i, ]
  saveRDS(run_cell(row), file.path(PER_CELL, sprintf("cell_%04d.rds", row$cell_id)))
  if (i %% 50L == 0L)
    cat(sprintf("[%s] %d/%d (last k=%d N=%d q=%.2f prev=%.2f)\n",
                format(Sys.time()), i, nrow(todo), row$k, row$N, row$q, row$prev),
        file = PROGRESS, append = TRUE)
  NULL
}, future.seed = TRUE))
cat(sprintf("[%s] stage6 split %d/%d COMPLETE. Wall: %.1f min.\n",
            format(Sys.time()), rem, mod, as.numeric(Sys.time() - t0, units = "mins")),
    file = PROGRESS, append = TRUE)
cat(sprintf("Done. Wall: %.1f min.\n", as.numeric(Sys.time() - t0, units = "mins")))
