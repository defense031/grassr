#!/usr/bin/env Rscript
# Stage 7B -- Option-B delta_hat POWER sweep (grassr v0.7.1).
# =====================================================================
# A > 0 alternative sweep through the NEW check_asymmetry() pipeline; TPR is
# computed later (analyze_power_deltaB.R) as the rate at which each draw's
# delta_hat lands >= 95th / >= 99th percentile of the MATCHED Option-B null
# ridge (stage6_null_deltaB). The 0.7.0 power program (stage7_power) stored
# only old-units delta, so it regenerates here; this time we also store the
# three per-coefficient implied q_hats (never again re-simulate for storage).
#
# GRID: SAME 11,000-cell (prev x k x N x q x A) grid as run_power.R, 2,000
# reps per cell = 22,000,000 draws. Same half-pattern (mean-norm A) gen_panel.
# DETERMINISM: seed = MASTER_SEED + cell_id set ONCE per cell; batching and
# WORKERS do not change rbinom call order -> bit-identical across machines and
# resume. lme4 irrelevant (no RNG in the panel/surface path).
#
# LAUNCH (per machine): CANDIDATE_LIB, SPLIT_MOD, SPLIT_REM, WORKERS.

suppressWarnings(suppressMessages({

get_script_path <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[[1]])))
  of <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(of)) return(normalizePath(of))
  stop("cannot locate script; run via Rscript or set GRASSR_SIM_ROOT")
}
S7       <- dirname(get_script_path())        # .../tier1/stage7_power_deltaB
PER_CELL <- file.path(S7, "per_cell")
PROGRESS <- file.path(S7, "PROGRESS.txt")
dir.create(PER_CELL, recursive = TRUE, showWarnings = FALSE)

CAND_LIB <- Sys.getenv("CANDIDATE_LIB", "")
if (!nzchar(CAND_LIB) || !dir.exists(CAND_LIB))
  stop("CANDIDATE_LIB must point to the built Option-B candidate grassr ",
       "library. Got: '", CAND_LIB, "'")
library(future.apply)
.libPaths(c(CAND_LIB, .libPaths()))
library(grassr, lib.loc = CAND_LIB)

if (requireNamespace("lme4", quietly = TRUE))
  message("[advisory] lme4 is visible: per-rep glmer will run (much slower). ",
          "It does not change delta_hat or determinism; park it for speed.")

{
  set.seed(1L)
  probe <- tryCatch(check_asymmetry(matrix(rbinom(30L * 4L, 1L, 0.7), 30L, 4L)),
                    error = function(e) e)
  if (inherits(probe, "error"))
    stop("candidate probe failed (", conditionMessage(probe),
         "). Densified surface present in CANDIDATE_LIB sysdata?")
  if (!("implied_q" %in% names(probe$panel)))
    stop("loaded grassr is NOT the Option-B build (no `implied_q`). ",
         "Rebuild CANDIDATE_LIB from the v0.7.1 branch.")
}

}))  # end setup

MASTER_SEED <- 20260707L                       # stage 7B regen master seed
REPS        <- 2000L
N_PER_BATCH <- as.integer(Sys.getenv("N_PER_BATCH", "2000"))
SMOKE       <- identical(Sys.getenv("SMOKE"), "1")

# ---- 11,000-cell grid (SAME as run_power.R). DO NOT reorder: fixes cell_id.
grid <- expand.grid(prev = c(0.05, 0.20, 0.50, 0.80, 0.95),
                    k = c(2L, 3L, 5L, 6L, 8L, 10L, 15L, 25L),
                    N = c(15L, 20L, 30L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L),
                    q = c(0.65, 0.75, 0.85, 0.92, 0.97),
                    A = c(0.05, 0.10, 0.15, 0.20, 0.30),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$cell_id <- seq_len(nrow(grid))
if (SMOKE) grid <- grid[c(1L, nrow(grid) %/% 2L, nrow(grid)), , drop = FALSE]
reps <- if (SMOKE) 25L else REPS

mod <- as.integer(Sys.getenv("SPLIT_MOD", "1"))
rem <- as.integer(Sys.getenv("SPLIT_REM", "0"))
grid <- grid[grid$cell_id %% mod == rem, , drop = FALSE]

cell_file <- function(id) file.path(PER_CELL, sprintf("cell_%05d.rds", id))
is_complete <- function(id) {
  f <- cell_file(id); if (!file.exists(f)) return(FALSE)
  ok <- tryCatch({ x <- readRDS(f)
    is.list(x) && !is.null(x$draws) && is.data.frame(x$draws) &&
      !is.null(x$target_draws) && nrow(x$draws) == x$target_draws &&
      # complete means complete AT THIS RUN'S TARGET (guards stale
      # SMOKE/partial-target files from being skipped as done)
      x$target_draws == reps
  }, error = function(e) FALSE)
  isTRUE(ok)
}
keep <- !vapply(grid$cell_id, is_complete, logical(1L))
todo <- grid[keep, , drop = FALSE]
cat(sprintf("stage7B | cells(split %d/%d): %d | done: %d | todo: %d | draws(todo): %s\n",
            rem, mod, nrow(grid), sum(!keep), nrow(todo),
            format(nrow(todo) * reps, big.mark = ",")))

# ---- alternative panel: half-pattern, mean-norm A (SAME as run_power.R) ---
gen_panel <- function(N, k, prev, q, A) {
  s  <- sample(c(-1, 1), k, replace = TRUE)
  Se <- pmin(pmax(q + s * A / 2, 0.001), 0.999)
  Sp <- pmin(pmax(q - s * A / 2, 0.001), 0.999)
  C  <- rbinom(N, 1L, prev); Y <- matrix(0L, N, k)
  for (j in seq_len(k)) Y[, j] <- rbinom(N, 1L, ifelse(C == 1L, Se[j], 1 - Sp[j]))
  Y
}

one_draw <- function(N, k, prev, q, A) {
  res <- tryCatch(suppressWarnings(suppressMessages(
    check_asymmetry(gen_panel(N, k, prev, q, A)))),
    error = function(e) NULL)
  if (is.null(res)) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  p <- res$panel
  q_of <- function(nm) { i <- which(p$coefficient == nm)
    if (length(i) == 1L) p$implied_q[i] else NA_real_ }
  c(as.numeric(res$delta_hat), q_of("pabak"), q_of("mean_ac1"), q_of("fleiss_kappa"))
}

run_cell <- function(row) {
  set.seed(MASTER_SEED + row$cell_id)   # ONCE, before any RNG for this cell
  delta <- numeric(reps); qp <- numeric(reps); qa <- numeric(reps); qf <- numeric(reps)
  pos <- 0L; left <- reps
  while (left > 0L) {
    nb <- min(N_PER_BATCH, left)
    b <- vapply(seq_len(nb), function(r) one_draw(row$N, row$k, row$prev, row$q, row$A),
                numeric(4L))
    ix <- (pos + 1L):(pos + nb)
    delta[ix] <- b[1L, ]; qp[ix] <- b[2L, ]; qa[ix] <- b[3L, ]; qf[ix] <- b[4L, ]
    pos <- pos + nb; left <- left - nb
  }
  list(cell = row[, c("prev", "k", "N", "q", "A", "cell_id")],
       draws = data.frame(delta = delta, q_pabak = qp, q_mean_ac1 = qa,
                          q_fleiss_kappa = qf),
       target_draws = reps,
       n_finite = sum(is.finite(delta)),
       master_seed = MASTER_SEED,
       r_version = R.version.string)
}

save_cell <- function(obj, id) {
  f <- cell_file(id); tmp <- paste0(f, ".tmp", Sys.getpid())
  saveRDS(obj, tmp); file.rename(tmp, f)
}

n_workers <- max(1L, min(as.integer(Sys.getenv("WORKERS", "8")),
                         parallel::detectCores() - 1L))
plan(multicore, workers = n_workers)
cat(sprintf("[%s] stage7B start: %d todo, %d workers, split %d/%d, seed %d, %s\n",
            format(Sys.time()), nrow(todo), n_workers, rem, mod, MASTER_SEED,
            R.version.string),
    file = PROGRESS, append = file.exists(PROGRESS))

t0 <- Sys.time()
if (nrow(todo) > 0L)
  invisible(future_lapply(seq_len(nrow(todo)), function(i) {
    row <- todo[i, ]
    save_cell(run_cell(row), row$cell_id)
    if (i %% 100L == 0L)
      cat(sprintf("[%s] %d/%d\n", format(Sys.time()), i, nrow(todo)),
          file = PROGRESS, append = TRUE)
    NULL
  }, future.seed = TRUE))

cat(sprintf("[%s] stage7B split %d/%d COMPLETE. Wall: %.1f min.\n",
            format(Sys.time()), rem, mod, as.numeric(Sys.time() - t0, units = "mins")),
    file = PROGRESS, append = TRUE)
writeLines(sprintf("stage7B split %d/%d complete %s", rem, mod, format(Sys.time())),
           file.path(S7, sprintf("split_%d_of_%d.DONE", rem, mod)))
cat(sprintf("Done. Wall: %.1f min.\n", as.numeric(Sys.time() - t0, units = "mins")))
