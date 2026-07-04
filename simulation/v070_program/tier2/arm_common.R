# Shared harness for Tier 2 arms A-C (sourced by arm_*.R).
# Each arm script defines: ARM (dir name), GRID (data.frame with cell_id
# and whatever columns its gen_panel needs), gen_panel(row) -> Y matrix,
# then calls run_arm(). Mechanics mirror tier1 stage 3: candidate lib,
# per-cell RDS resume, SPLIT_MOD/SPLIT_REM, multicore fork, per-cell seeds.
#
# Per rep we record the three delta-family coefficients, their
# candidate-surface percentiles, and delta_hat — so every arm supports
# both deliverables (percentile drift vs its null anchor, and delta-hat
# size inflation) from one output format.

suppressPackageStartupMessages({ library(future.apply) })

ROOT <- Sys.getenv("GRASS_SIM_ROOT",
                   "/Users/austinsemmel/Desktop/PABAK_Investigation")
setwd(ROOT)
CAND_LIB <- file.path(ROOT,
  "grassr/simulation/v070_program/tier1/stage3_threshold_grid/lib")
stopifnot(dir.exists(file.path(CAND_LIB, "grassr")))
.libPaths(c(CAND_LIB, .libPaths()))
suppressPackageStartupMessages(library(grassr, lib.loc = CAND_LIB))

N_REPS      <- 1000L
MASTER_SEED <- 20260706L
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

mean_pabak <- function(Y) {
  ij <- utils::combn(ncol(Y), 2)
  mean(apply(ij, 2, function(p) 2 * mean(Y[, p[1]] == Y[, p[2]]) - 1))
}
mean_ac1 <- function(Y) {
  ij <- utils::combn(ncol(Y), 2)
  vals <- apply(ij, 2, function(p) {
    pa     <- mean(Y[, p[1]] == Y[, p[2]])
    pi_bar <- (mean(Y[, p[1]]) + mean(Y[, p[2]])) / 2
    pe     <- 2 * pi_bar * (1 - pi_bar)
    if (pe >= 0.999999) NA else (pa - pe) / (1 - pe)
  })
  mean(vals, na.rm = TRUE)
}
fleiss_kappa_fast <- function(Y) {
  N <- nrow(Y); k <- ncol(Y)
  n1 <- rowSums(Y); n0 <- k - n1
  Pbar <- mean((n1 * (n1 - 1) + n0 * (n0 - 1)) / (k * (k - 1)))
  p1 <- mean(Y); pe <- p1^2 + (1 - p1)^2
  if (pe >= 0.999999) return(NA_real_)
  (Pbar - pe) / (1 - pe)
}
invert_pct <- function(obs, metric, pi_hat, k, N) {
  if (is.na(obs)) return(NA_real_)
  res <- try(position_on_surface(obs_value = obs, metric = metric,
                                 pi_hat = pi_hat, k = k, N = N,
                                 method = "empirical"), silent = TRUE)
  if (inherits(res, "try-error") || is.null(res$percentile)) return(NA_real_)
  res$percentile * 100
}

run_arm <- function(ARM, GRID, gen_panel) {
  ARM_DIR  <- file.path(ROOT, "grassr/simulation/v070_program/tier2", ARM)
  PER_CELL <- file.path(ARM_DIR, "per_cell")
  PROGRESS <- file.path(ARM_DIR, "PROGRESS.txt")
  dir.create(PER_CELL, recursive = TRUE, showWarnings = FALSE)

  grid <- GRID
  if (SMOKE) grid <- grid[c(1L, nrow(grid) %/% 2L, nrow(grid)), ]
  mod <- as.integer(Sys.getenv("SPLIT_MOD", "1"))
  rem <- as.integer(Sys.getenv("SPLIT_REM", "0"))
  grid <- grid[grid$cell_id %% mod == rem, , drop = FALSE]
  done <- as.integer(sub("^cell_([0-9]+)\\.rds$", "\\1",
                         list.files(PER_CELL, pattern = "^cell_[0-9]+\\.rds$")))
  todo <- grid[!grid$cell_id %in% done, , drop = FALSE]
  reps <- if (SMOKE) 25L else N_REPS
  cat(sprintf("%s | cells(split): %d | done: %d | todo: %d | reps: %d\n",
              ARM, nrow(grid), sum(grid$cell_id %in% done), nrow(todo), reps))

  run_cell <- function(row) {
    set.seed(MASTER_SEED + row$cell_id)
    out <- vector("list", reps)
    for (r in seq_len(reps)) {
      Y  <- gen_panel(row)
      pp <- mean(Y)
      if (pp <= 0 || pp >= 1) next
      pi_hat <- max(min(pp, 0.999), 0.001)
      obs <- c(pabak        = mean_pabak(Y),
               fleiss_kappa = fleiss_kappa_fast(Y),
               mean_ac1     = mean_ac1(Y))
      pcts <- mapply(invert_pct, obs, names(obs),
                     MoreArgs = list(pi_hat = pi_hat, k = row$k, N = row$N))
      fin <- pcts[is.finite(pcts)]
      out[[r]] <- data.frame(rep = r, t(obs), t(pcts),
        delta = if (length(fin) >= 2L) diff(range(fin)) else NA_real_)
    }
    cbind(row, do.call(rbind, out), row.names = NULL)
  }

  n_workers <- max(1L, min(as.integer(Sys.getenv("WORKERS", "10")),
                           parallel::detectCores() - 2L))
  plan(multicore, workers = n_workers)
  cat(sprintf("[%s] %s start: %d todo, %d workers, split %d/%d, %s\n",
              format(Sys.time()), ARM, nrow(todo), n_workers, rem, mod,
              R.version.string),
      file = PROGRESS, append = file.exists(PROGRESS))
  t0 <- Sys.time()
  invisible(future_lapply(seq_len(nrow(todo)), function(i) {
    row <- todo[i, ]
    saveRDS(run_cell(row),
            file.path(PER_CELL, sprintf("cell_%04d.rds", row$cell_id)))
    if (i %% 20L == 0L)
      cat(sprintf("[%s] %d/%d\n", format(Sys.time()), i, nrow(todo)),
          file = PROGRESS, append = TRUE)
    NULL
  }, future.seed = TRUE))
  cat(sprintf("[%s] %s split %d/%d COMPLETE. Wall: %.1f min.\n",
              format(Sys.time()), ARM, rem, mod,
              as.numeric(Sys.time() - t0, units = "mins")),
      file = PROGRESS, append = TRUE)
  cat(sprintf("%s done. Wall: %.1f min.\n", ARM,
              as.numeric(Sys.time() - t0, units = "mins")))
}
