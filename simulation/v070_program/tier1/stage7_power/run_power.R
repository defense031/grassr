#!/usr/bin/env Rscript
# Stage 7 — production-path POWER sweep for the percentile-flag convention.
# A > 0 draws (half-pattern asymmetry, mean-norm A) through
# check_asymmetry(); nulls came from stage 6. TPR is later computed as
# P(delta-hat percentile on matched null >= 95th/99th) per (k,N,q,A).
# Grid: 2,200 design cells x A in {0.05,0.10,0.15,0.20,0.30}, 2,000 reps.
# Seeds: 20260709 + cell_id.
ROOT <- Sys.getenv("GRASS_SIM_ROOT", "/Users/austinsemmel/Desktop/PABAK_Investigation")
setwd(ROOT)
CAND_LIB <- file.path(ROOT, "grassr/simulation/v070_program/tier1/stage3_threshold_grid/lib")
suppressPackageStartupMessages(library(future.apply))
.libPaths(CAND_LIB)
suppressPackageStartupMessages(library(grassr, lib.loc = CAND_LIB))
if (requireNamespace("lme4", quietly = TRUE)) stop("lme4 visible; park it first.")
S7 <- file.path(ROOT, "grassr/simulation/v070_program/tier1/stage7_power")
PER_CELL <- file.path(S7, "per_cell"); dir.create(PER_CELL, showWarnings = FALSE)
PROGRESS <- file.path(S7, "PROGRESS.txt")
grid <- expand.grid(prev = c(0.05, 0.20, 0.50, 0.80, 0.95),
                    k = c(2L, 3L, 5L, 6L, 8L, 10L, 15L, 25L),
                    N = c(15L, 20L, 30L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L),
                    q = c(0.65, 0.75, 0.85, 0.92, 0.97),
                    A = c(0.05, 0.10, 0.15, 0.20, 0.30),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$cell_id <- seq_len(nrow(grid))
REPS <- 2000L
SMOKE <- identical(Sys.getenv("SMOKE"), "1")
if (SMOKE) { grid <- grid[c(1L, nrow(grid) %/% 2L, nrow(grid)), ] }
mod <- as.integer(Sys.getenv("SPLIT_MOD", "1")); rem <- as.integer(Sys.getenv("SPLIT_REM", "0"))
grid <- grid[grid$cell_id %% mod == rem, , drop = FALSE]
done <- as.integer(sub("^cell_([0-9]+)\\.rds$", "\\1", list.files(PER_CELL)))
todo <- grid[!grid$cell_id %in% done, , drop = FALSE]
reps <- if (SMOKE) 25L else REPS
cat(sprintf("stage7 | split cells: %d | todo: %d | draws: %s\n",
    nrow(grid), nrow(todo), format(sum(nrow(todo)) * reps, big.mark = ",")))
gen_panel <- function(N, k, prev, q, A) {
  s  <- sample(c(-1, 1), k, replace = TRUE)
  Se <- pmin(pmax(q + s * A / 2, 0.001), 0.999)
  Sp <- pmin(pmax(q - s * A / 2, 0.001), 0.999)
  C  <- rbinom(N, 1L, prev); Y <- matrix(0L, N, k)
  for (j in seq_len(k)) Y[, j] <- rbinom(N, 1L, ifelse(C == 1L, Se[j], 1 - Sp[j]))
  Y
}
plan(multicore, workers = max(1L, min(as.integer(Sys.getenv("WORKERS", "10")),
                                      parallel::detectCores() - 2L)))
t0 <- Sys.time()
invisible(future_lapply(seq_len(nrow(todo)), function(i) {
  row <- todo[i, ]; set.seed(20260709L + row$cell_id)
  d <- vapply(seq_len(reps), function(r) {
    Y <- gen_panel(row$N, row$k, row$prev, row$q, row$A)
    suppressWarnings(suppressMessages(
      tryCatch(check_asymmetry(Y)$delta_hat, error = function(e) NA_real_)))
  }, numeric(1L))
  saveRDS(list(cell = row, delta = d, n_finite = sum(is.finite(d))),
          file.path(PER_CELL, sprintf("cell_%05d.rds", row$cell_id)))
  if (i %% 100L == 0L) cat(sprintf("[%s] %d/%d\n", format(Sys.time()), i, nrow(todo)),
                           file = PROGRESS, append = TRUE)
  NULL
}, future.seed = TRUE))
cat(sprintf("[%s] stage7 split %d/%d COMPLETE. Wall: %.1f min.\n",
    format(Sys.time()), rem, mod, as.numeric(Sys.time() - t0, units = "mins")),
    file = PROGRESS, append = TRUE)
