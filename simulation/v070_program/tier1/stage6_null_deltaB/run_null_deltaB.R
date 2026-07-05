#!/usr/bin/env Rscript
# Stage 6B -- Option-B delta_hat NULL regeneration (grassr v0.7.1).
# =====================================================================
# WHY THIS EXISTS.  delta_hat was redefined 2026-07-05 (Option B, ratified;
# design/v0.7.1_position_redesign.md): it is now the IMPLIED-QUALITY spread
# of the agreement family (PABAK / mean AC1 / Fleiss kappa) in quality
# percentage points, computed inside check_asymmetry(). The 0.7.0 null
# (stage6_production_null) stored ONLY old-units delta draws, so the whole
# null regenerates through the NEW pipeline. Each draw records
# check_asymmetry(Y)$delta_hat AND the three per-coefficient implied q_hats
# (never again shipping a null we cannot re-derive).
#
# GRID (Option B, ratified reading).  The shipped null is indexed by the 440
# (k, N, q) ridges. The 0.7.0 program enumerated prevalence as a separate
# grid dimension (2,200 prev x k x N x q cells) and POOLED over prevalence at
# extract, reaching 50k/ridge only after two top-up passes (6b, 6c). This
# regen does it cleanly in one program: 440 cells keyed by (k, N, q), each
# running 50,000 draws = 10,000 at each of the five calibration prevalences
# {.05,.20,.50,.80,.95} swept INSIDE the cell. Same k/N/q/prev vectors and
# same null gen_panel() as run_production_null.R; only the cell bookkeeping
# differs (prev marginalized into the ridge instead of into extract-time
# pooling). Total: 440 x 50,000 = 22,000,000 production draws.
#
# DETERMINISM.  seed = MASTER_SEED + cell_id is set ONCE per cell before any
# RNG. The sequence of gen_panel() draws within a cell (prev-major, then
# sequential) is INDEPENDENT of WORKERS and of N_PER_BATCH -- batching only
# chunks the vapply, never the rbinom call order -- so results are bit-
# identical across machines, worker counts, and resume/rebalance. lme4 does
# not enter: compute_panel()/position_on_surface() use no RNG, so whether
# lme4 is parked changes only speed and whether ICC appears on the panel
# (ICC never enters delta_hat).
#
# LAUNCH (per machine): supply CANDIDATE_LIB (Option-B candidate grassr),
# SPLIT_MOD, SPLIT_REM, WORKERS. See stage6_null_deltaB/README.md.

suppressWarnings(suppressMessages({

# ---- paths: derived from THIS script's location, never from cwd ----------
get_script_path <- function() {
  a <- commandArgs(FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[[1]])))
  of <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(of)) return(normalizePath(of))
  stop("cannot locate script; run via Rscript or set GRASSR_SIM_ROOT")
}
SCRIPT_DIR <- dirname(get_script_path())            # .../tier1/stage6_null_deltaB
STAGE6     <- SCRIPT_DIR
PER_CELL   <- file.path(STAGE6, "per_cell")
PROGRESS   <- file.path(STAGE6, "PROGRESS.txt")
dir.create(PER_CELL, recursive = TRUE, showWarnings = FALSE)

# ---- candidate library (Option-B grassr) ---------------------------------
CAND_LIB <- Sys.getenv("CANDIDATE_LIB", "")
if (!nzchar(CAND_LIB) || !dir.exists(CAND_LIB))
  stop("CANDIDATE_LIB must point to the built Option-B candidate grassr ",
       "library. Got: '", CAND_LIB, "'")
# Load parallel machinery BEFORE the path change (loaded namespaces persist).
library(future.apply)
.libPaths(c(CAND_LIB, .libPaths()))
library(grassr, lib.loc = CAND_LIB)

# lme4 advisory only (never required; does not affect delta_hat/determinism).
if (requireNamespace("lme4", quietly = TRUE))
  message("[advisory] lme4 is visible: per-rep glmer will run (much slower). ",
          "It does not change delta_hat or determinism; park it for speed ",
          "and to match the shipped-null build.")

# ---- Option-B guard: fail fast if the loaded grassr predates implied q ----
{
  set.seed(1L)  # isolated, deterministic probe; cells re-seed independently
  probe_Y <- matrix(rbinom(30L * 4L, 1L, 0.7), nrow = 30L, ncol = 4L)
  probe <- tryCatch(check_asymmetry(probe_Y), error = function(e) e)
  if (inherits(probe, "error"))
    stop("candidate probe failed (", conditionMessage(probe), "). Is the ",
         "densified q_hat surface present in CANDIDATE_LIB's sysdata?")
  if (!("implied_q" %in% names(probe$panel)))
    stop("loaded grassr is NOT the Option-B build (panel has no `implied_q`). ",
         "Rebuild CANDIDATE_LIB from the v0.7.1 branch before launching.")
}

}))  # end suppressWarnings/suppressMessages of setup

# ---- run parameters ------------------------------------------------------
MASTER_SEED    <- 20260706L                 # stage 6B regen master seed
DRAWS_PER_PREV <- 10000L                     # x 5 prev = 50,000 per ridge
PREVS          <- c(0.05, 0.20, 0.50, 0.80, 0.95)
N_PER_BATCH    <- as.integer(Sys.getenv("N_PER_BATCH", "2000"))
SMOKE          <- identical(Sys.getenv("SMOKE"), "1")

# ---- 440-cell (k, N, q) grid. DO NOT reorder: row order fixes cell_id -----
# fixes seed. Same k/N/q vectors as run_production_null.R.
grid <- expand.grid(k = c(2L, 3L, 5L, 6L, 8L, 10L, 15L, 25L),
                    N = c(15L, 20L, 30L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L),
                    q = c(0.65, 0.75, 0.85, 0.92, 0.97),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$cell_id <- seq_len(nrow(grid))
if (SMOKE) {
  grid <- grid[c(1L, nrow(grid) %/% 2L, nrow(grid)), , drop = FALSE]
  DRAWS_PER_PREV <- 40L
}

# ---- two-machine split ---------------------------------------------------
mod <- as.integer(Sys.getenv("SPLIT_MOD", "1"))
rem <- as.integer(Sys.getenv("SPLIT_REM", "0"))
grid <- grid[grid$cell_id %% mod == rem, , drop = FALSE]

cell_file <- function(id) file.path(PER_CELL, sprintf("cell_%04d.rds", id))

# ---- resume: skip cells whose RDS exists, is READABLE, and is COMPLETE ----
is_complete <- function(id) {
  f <- cell_file(id)
  if (!file.exists(f)) return(FALSE)
  ok <- tryCatch({
    x <- readRDS(f)
    is.list(x) && !is.null(x$draws) && is.data.frame(x$draws) &&
      !is.null(x$target_draws) && nrow(x$draws) == x$target_draws &&
      # guard against stale SMOKE/partial-target files: complete means
      # complete AT THIS RUN'S TARGET, not at whatever the file recorded
      x$target_draws == length(PREVS) * DRAWS_PER_PREV
  }, error = function(e) FALSE)   # corrupt / partial -> recompute
  isTRUE(ok)
}
keep <- !vapply(grid$cell_id, is_complete, logical(1L))
todo <- grid[keep, , drop = FALSE]
target_per_cell <- length(PREVS) * DRAWS_PER_PREV
cat(sprintf("stage6B | cells(split %d/%d): %d | done: %d | todo: %d | draws(todo): %s\n",
            rem, mod, nrow(grid), sum(!keep), nrow(todo),
            format(nrow(todo) * target_per_cell, big.mark = ",")))

# ---- draw machinery (null: A = 0; same gen_panel as production null) ------
gen_panel <- function(N, k, prev, q) {
  C <- rbinom(N, 1L, prev)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) Y[, j] <- rbinom(N, 1L, ifelse(C == 1L, q, 1 - q))
  Y
}

# one draw -> c(delta_hat, q_pabak, q_mean_ac1, q_fleiss_kappa)
# implied q_hats are $panel$implied_q (0-1); delta_hat is in quality pp.
one_draw <- function(N, k, prev, q) {
  res <- tryCatch(suppressWarnings(suppressMessages(
    check_asymmetry(gen_panel(N, k, prev, q)))),
    error = function(e) NULL)
  if (is.null(res)) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  p <- res$panel
  q_of <- function(nm) { i <- which(p$coefficient == nm)
    if (length(i) == 1L) p$implied_q[i] else NA_real_ }
  c(as.numeric(res$delta_hat), q_of("pabak"), q_of("mean_ac1"), q_of("fleiss_kappa"))
}

run_cell <- function(row) {
  set.seed(MASTER_SEED + row$cell_id)   # ONCE, before any RNG for this cell
  tgt <- length(PREVS) * DRAWS_PER_PREV
  delta <- numeric(tgt); qp <- numeric(tgt); qa <- numeric(tgt)
  qf <- numeric(tgt); pv <- numeric(tgt)
  pos <- 0L
  for (prev in PREVS) {              # prev swept INSIDE the cell (equal weight)
    left <- DRAWS_PER_PREV
    while (left > 0L) {              # batched only to bound intermediates
      nb <- min(N_PER_BATCH, left)
      b <- vapply(seq_len(nb), function(r) one_draw(row$N, row$k, prev, row$q),
                  numeric(4L))       # 4 x nb
      ix <- (pos + 1L):(pos + nb)
      delta[ix] <- b[1L, ]; qp[ix] <- b[2L, ]; qa[ix] <- b[3L, ]; qf[ix] <- b[4L, ]
      pv[ix] <- prev
      pos <- pos + nb; left <- left - nb
    }
  }
  list(cell = row[, c("k", "N", "q", "cell_id")],
       draws = data.frame(delta = delta, q_pabak = qp, q_mean_ac1 = qa,
                          q_fleiss_kappa = qf, prev = pv),
       target_draws = tgt,
       n_finite = sum(is.finite(delta)),
       master_seed = MASTER_SEED,
       r_version = R.version.string)
}

# atomic per-cell write (temp + rename) so a crash never leaves a partial RDS
save_cell <- function(obj, id) {
  f <- cell_file(id); tmp <- paste0(f, ".tmp", Sys.getpid())
  saveRDS(obj, tmp); file.rename(tmp, f)
}

# ---- parallel over cells (each cell fully self-seeded) -------------------
n_workers <- max(1L, min(as.integer(Sys.getenv("WORKERS", "8")),
                         parallel::detectCores() - 1L))
plan(multicore, workers = n_workers)
cat(sprintf("[%s] stage6B start: %d todo, %d workers, split %d/%d, seed %d, %s\n",
            format(Sys.time()), nrow(todo), n_workers, rem, mod, MASTER_SEED,
            R.version.string),
    file = PROGRESS, append = file.exists(PROGRESS))

t0 <- Sys.time()
if (nrow(todo) > 0L)
  invisible(future_lapply(seq_len(nrow(todo)), function(i) {
    row <- todo[i, ]
    save_cell(run_cell(row), row$cell_id)
    if (i %% 10L == 0L)
      cat(sprintf("[%s] %d/%d (k=%d N=%d q=%.2f)\n",
                  format(Sys.time()), i, nrow(todo), row$k, row$N, row$q),
          file = PROGRESS, append = TRUE)
    NULL
  }, future.seed = TRUE))

cat(sprintf("[%s] stage6B split %d/%d COMPLETE. Wall: %.1f min.\n",
            format(Sys.time()), rem, mod, as.numeric(Sys.time() - t0, units = "mins")),
    file = PROGRESS, append = TRUE)
writeLines(sprintf("stage6B split %d/%d complete %s", rem, mod, format(Sys.time())),
           file.path(STAGE6, sprintf("split_%d_of_%d.DONE", rem, mod)))
cat(sprintf("Done. Wall: %.1f min.\n", as.numeric(Sys.time() - t0, units = "mins")))
