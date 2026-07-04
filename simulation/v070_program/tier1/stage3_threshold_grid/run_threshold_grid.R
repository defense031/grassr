#!/usr/bin/env Rscript
# Stage 3 — delta-hat threshold grid, q-conditioned (grassr v0.7.0 Tier 1).
#
# Full grid: prev {0.05,0.20,0.50,0.80,0.95} x k {2,3,5,6,8,10,15,25}
#            x N {15,20,30,50,75,100,150,200,300,500,1000}
#            x A {0,0.05,0.10,0.15,0.20,0.30} x q {0.65,0.75,0.85,0.92,0.97}
# = 13,200 cells x 1,000 reps. A = 0 rows give the null (size calibration);
# A > 0 rows give TPR. The q axis is the bias check: the shipped table was
# calibrated at q = 0.85 only.
#
# delta-hat family matches the SHIPPED 3-coefficient definition (PABAK,
# mean AC1, Fleiss kappa) — the pre-Pass-6 sweep included krippendorff_a.
# Percentile inversion runs against the STAGE-2 CANDIDATE surface via a
# grassr build installed into an isolated library (prepare_candidate_lib.R).
#
# Env: GRASS_SIM_ROOT, SPLIT_MOD/SPLIT_REM, WORKERS, SMOKE=1.
# Seeds: master_seed 20260705 + cell_id (per design doc).
# Resume: per-cell RDS in per_cell/.

suppressPackageStartupMessages({ library(future.apply) })

ROOT <- Sys.getenv("GRASS_SIM_ROOT",
                   "/Users/austinsemmel/Desktop/PABAK_Investigation")
setwd(ROOT)
STAGE3   <- file.path(ROOT, "grassr/simulation/v070_program/tier1/stage3_threshold_grid")
PER_CELL <- file.path(STAGE3, "per_cell")
PROGRESS <- file.path(STAGE3, "PROGRESS.txt")
CAND_LIB <- file.path(STAGE3, "lib")
dir.create(PER_CELL, recursive = TRUE, showWarnings = FALSE)

stopifnot(dir.exists(file.path(CAND_LIB, "grassr")))
.libPaths(c(CAND_LIB, .libPaths()))
suppressPackageStartupMessages(library(grassr, lib.loc = CAND_LIB))
cat(sprintf("grassr candidate lib loaded: %s (surface cells: %d)\n",
            as.character(packageVersion("grassr")),
            nrow(get("empirical_q_hat_surface",
                     envir = asNamespace("grassr"))$index)))

PI_GRID <- c(0.05, 0.20, 0.50, 0.80, 0.95)
K_GRID  <- c(2L, 3L, 5L, 6L, 8L, 10L, 15L, 25L)
N_GRID  <- c(15L, 20L, 30L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L)
A_GRID  <- c(0, 0.05, 0.10, 0.15, 0.20, 0.30)
Q_GRID  <- c(0.65, 0.75, 0.85, 0.92, 0.97)
N_REPS  <- 1000L
MASTER_SEED <- 20260705L

SMOKE <- identical(Sys.getenv("SMOKE"), "1")

gen_panel <- function(N, k, prev, q, A) {
  s  <- sample(c(-1, 1), k, replace = TRUE)
  Se <- pmin(pmax(q + s * A / 2, 0.001), 0.999)
  Sp <- pmin(pmax(q - s * A / 2, 0.001), 0.999)
  C  <- rbinom(N, 1L, prev)
  Y  <- matrix(0L, N, k)
  for (j in seq_len(k)) {
    Y[, j] <- rbinom(N, 1L, ifelse(C == 1L, Se[j], 1 - Sp[j]))
  }
  Y
}

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
  # Fleiss 1971 direct computation (no irr dependency; identical value)
  N <- nrow(Y); k <- ncol(Y)
  n1 <- rowSums(Y); n0 <- k - n1
  P_i <- (n1 * (n1 - 1) + n0 * (n0 - 1)) / (k * (k - 1))
  Pbar <- mean(P_i)
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

compute_delta <- function(Y, pi_hat, k, N) {
  obs <- c(pabak        = mean_pabak(Y),
           fleiss_kappa = fleiss_kappa_fast(Y),
           mean_ac1     = mean_ac1(Y))
  pcts <- mapply(invert_pct, obs, names(obs),
                 MoreArgs = list(pi_hat = pi_hat, k = k, N = N))
  pcts <- pcts[is.finite(pcts)]
  if (length(pcts) < 2L) return(NA_real_)
  diff(range(pcts))
}

run_cell <- function(row) {
  set.seed(MASTER_SEED + row$cell_id)
  reps <- if (SMOKE) 25L else N_REPS
  out <- numeric(reps)
  for (r in seq_len(reps)) {
    Y  <- gen_panel(row$N, row$k, row$prev, row$q, row$A)
    pp <- mean(Y)
    if (pp <= 0 || pp >= 1) { out[r] <- NA_real_; next }
    out[r] <- compute_delta(Y, max(min(pp, 0.999), 0.001), row$k, row$N)
  }
  data.frame(cell_id = row$cell_id, prev = row$prev, k = row$k, N = row$N,
             A = row$A, q = row$q, rep = seq_len(reps), delta = out)
}

grid <- expand.grid(prev = PI_GRID, k = K_GRID, N = N_GRID,
                    A = A_GRID, q = Q_GRID,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$cell_id <- seq_len(nrow(grid))
if (SMOKE) grid <- grid[c(1L, nrow(grid) %/% 2L, nrow(grid)), ]

mod <- as.integer(Sys.getenv("SPLIT_MOD", "1"))
rem <- as.integer(Sys.getenv("SPLIT_REM", "0"))
grid <- grid[grid$cell_id %% mod == rem, , drop = FALSE]

done_ids <- as.integer(sub("^cell_([0-9]+)\\.rds$", "\\1",
                           list.files(PER_CELL, pattern = "^cell_[0-9]+\\.rds$")))
todo <- grid[!grid$cell_id %in% done_ids, , drop = FALSE]
cat(sprintf("stage3 | cells(this split): %d | done: %d | todo: %d | reps/cell: %d\n",
            nrow(grid), sum(grid$cell_id %in% done_ids), nrow(todo),
            if (SMOKE) 25L else N_REPS))

n_workers <- max(1L, min(as.integer(Sys.getenv("WORKERS", "10")),
                         parallel::detectCores() - 2L))
plan(multicore, workers = n_workers)   # fork: workers inherit candidate lib
cat(sprintf("Workers: %d (multicore) | split %d/%d | %s\n",
            n_workers, rem, mod, R.version.string))
cat(sprintf("[%s] stage3 start: %d todo, %d workers, split %d/%d\n",
            format(Sys.time()), nrow(todo), n_workers, rem, mod),
    file = PROGRESS, append = file.exists(PROGRESS))

t0 <- Sys.time()
invisible(future_lapply(seq_len(nrow(todo)), function(i) {
  row <- todo[i, ]
  res <- run_cell(row)
  saveRDS(res, file.path(PER_CELL, sprintf("cell_%05d.rds", row$cell_id)))
  if (i %% 100L == 0L)
    cat(sprintf("[%s] %d/%d (last cell k=%d N=%d q=%.2f A=%.2f)\n",
                format(Sys.time()), i, nrow(todo), row$k, row$N, row$q, row$A),
        file = PROGRESS, append = TRUE)
  NULL
}, future.seed = TRUE))
t1 <- Sys.time()

cat(sprintf("[%s] stage3 split %d/%d COMPLETE. Wall: %.1f min.\n",
            format(Sys.time()), rem, mod, as.numeric(t1 - t0, units = "mins")),
    file = PROGRESS, append = TRUE)
cat(sprintf("Done. Wall: %.1f min.\n", as.numeric(t1 - t0, units = "mins")))
