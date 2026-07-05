#!/usr/bin/env Rscript
# Tier 2 delta-B re-run -- Arms A/B/C through the v0.7.1 (Option-B) pipeline.
# =====================================================================
# WHY THIS EXISTS.  The 0.7.0 Tier-2 arm outputs (tier2/arm_*/per_cell) stored
# raw coefficients but NOT each rep's realized pi-hat, which the AC1/Fleiss
# quality inversions need, so the "arithmetic recompute" of TIER2_DESIGN.md is
# replaced by a RE-RUN of the three arms through the Option-B check_asymmetry()
# pipeline (TIER2_DESIGN.md, "Arms A/B/C under v0.7.1", 2026-07-05). The DGP
# extensions are reproduced EXACTLY from the 0.7.0 arm scripts
# (tier2/arm_a_item_difficulty.R / arm_b_correlated_errors.R /
# arm_c_asymmetry_patterns.R); ONLY the pipeline and per-rep storage modernize.
# This time every rep stores its implied q_hats and pooled percentiles, so the
# deliverables never require a re-simulate for storage again.
#
# GRIDS (EXACT 0.7.0 cell counts, verified):
#   Arm A -- item difficulty : sd_d in {0,0.5,1.0,1.5}   -> 180 cells
#   Arm B -- correlated errs  : rho  in {0,0.1,0.25,0.5}  -> 180 cells
#   Arm C -- asymmetry pattern: single/half/graded @A=.20 -> 108 cells
#   Total 468 cells x 1,000 reps = 468,000 panels through the Option-B pipeline.
#
# DETERMINISM.  seed = MASTER_SEED + arm_offset + cell_id is set ONCE per cell
# before any RNG (arm_offset 0/1000/2000 keeps seeds globally unique across
# arms). The split (SPLIT_MOD/SPLIT_REM) is over a GLOBAL job index and never
# touches the seed, so cells are bit-identical across machines and resume.
#
# LAUNCH (per machine): supply CANDIDATE_LIB (Option-B candidate grassr),
# SPLIT_MOD, SPLIT_REM, WORKERS. See tier2_deltaB/README.md. Chained behind the
# stage-7B split by chain_tier2_deltaB.sh.

suppressWarnings(suppressMessages({

# ---- paths: derived from THIS script's location, never from cwd ----------
get_script_path <- function() {
  a <- commandArgs(FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[[1]])))
  of <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(of)) return(normalizePath(of))
  stop("cannot locate script; run via Rscript")
}
SCRIPT_DIR <- dirname(get_script_path())            # .../tier2_deltaB
TIER2DB    <- SCRIPT_DIR
PER_CELL   <- file.path(TIER2DB, "per_cell")
PROGRESS   <- file.path(TIER2DB, "PROGRESS.txt")
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

# lme4 advisory only (ICC never enters delta_hat or the stored columns).
if (requireNamespace("lme4", quietly = TRUE))
  message("[advisory] lme4 is visible: per-rep glmer will run (much slower). ",
          "It does not change delta_hat or the stored q_hats; park it for speed.")

# ---- Option-B guard: fail fast if the loaded grassr predates implied q ----
{
  set.seed(1L)  # isolated, deterministic probe; cells re-seed independently
  probe_Y <- matrix(rbinom(30L * 4L, 1L, 0.7), nrow = 30L, ncol = 4L)
  probe <- tryCatch(check_asymmetry(probe_Y), error = function(e) e)
  if (inherits(probe, "error"))
    stop("candidate probe failed (", conditionMessage(probe), "). Is the ",
         "densified q_hat surface present in CANDIDATE_LIB's sysdata?")
  if (!all(c("implied_q", "percentile_pp") %in% names(probe$panel)))
    stop("loaded grassr is NOT the Option-B build (panel lacks implied_q / ",
         "percentile_pp). Rebuild CANDIDATE_LIB from the v0.7.1 branch.")
}

}))  # end suppressWarnings/suppressMessages of setup

# ---- run parameters ------------------------------------------------------
MASTER_SEED <- 20260708L                      # tier2 delta-B master seed
SMOKE       <- identical(Sys.getenv("SMOKE"), "1")
N_REPS      <- if (SMOKE) 25L else 1000L
reps        <- N_REPS

# ---- arm grids (EXACT reproductions of the 0.7.0 arm scripts) -------------
build_grid_A <- function() {
  g <- expand.grid(sd_d = c(0, 0.5, 1.0, 1.5),
                   kN   = c("3_50", "5_200", "8_100", "10_500", "25_1000"),
                   q    = c(0.75, 0.85, 0.92),
                   prev = c(0.10, 0.30, 0.50),
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g$k <- as.integer(sub("_.*", "", g$kN))
  g$N <- as.integer(sub(".*_", "", g$kN))
  g$kN <- NULL
  g$cell_id <- seq_len(nrow(g))
  g
}
build_grid_B <- function() {
  g <- expand.grid(rho  = c(0, 0.1, 0.25, 0.5),
                   kN   = c("3_50", "5_200", "8_100", "10_500", "25_1000"),
                   q    = c(0.75, 0.85, 0.92),
                   prev = c(0.10, 0.30, 0.50),
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g$k <- as.integer(sub("_.*", "", g$kN))
  g$N <- as.integer(sub(".*_", "", g$kN))
  g$kN <- NULL
  g$cell_id <- seq_len(nrow(g))
  g
}
A_FIX <- 0.20
build_grid_C <- function() {
  g <- expand.grid(pattern = c("half", "single", "ramp"),
                   k    = c(3L, 5L, 8L, 10L, 15L, 25L),
                   N    = c(50L, 200L, 1000L),
                   prev = c(0.20, 0.50),
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g$q <- 0.85
  g$cell_id <- seq_len(nrow(g))
  g
}

# ---- DGP generators (copied VERBATIM from the 0.7.0 arm scripts) ----------
gen_panel_A <- function(row) {
  d  <- rnorm(row$N, 0, row$sd_d)
  Se <- plogis(qlogis(row$q) - d)   # per-item accuracy, all raters
  C  <- rbinom(row$N, 1L, row$prev)
  p  <- ifelse(C == 1L, Se, 1 - Se)  # Sp_i = Se_i under symmetric raters
  Y  <- matrix(0L, row$N, row$k)
  for (j in seq_len(row$k)) Y[, j] <- rbinom(row$N, 1L, ifelse(C == 1L, Se, 1 - Se))
  Y
}
gen_panel_B <- function(row) {
  C <- rbinom(row$N, 1L, row$prev)
  z <- rnorm(row$N)
  thr <- qnorm(1 - row$q)   # error prob = 1 - q on both margins
  Y <- matrix(0L, row$N, row$k)
  for (j in seq_len(row$k)) {
    err <- (sqrt(row$rho) * z + sqrt(1 - row$rho) * rnorm(row$N)) < thr
    Y[, j] <- ifelse(err, 1L - C, C)
  }
  Y
}
gen_panel_C <- function(row) {
  k <- row$k; q <- row$q
  off <- switch(row$pattern,
    half   = sample(c(-1, 1), k, replace = TRUE) * A_FIX / 2,
    single = c(k * A_FIX / 2, rep(0, k - 1L)),
    ramp   = { v <- seq(-1, 1, length.out = k); v * (A_FIX / 2) / mean(abs(v)) })
  Se <- pmin(pmax(q + off, 0.001), 0.999)
  Sp <- pmin(pmax(q - off, 0.001), 0.999)
  C  <- rbinom(row$N, 1L, row$prev)
  Y  <- matrix(0L, row$N, k)
  for (j in seq_len(k)) Y[, j] <- rbinom(row$N, 1L, ifelse(C == 1L, Se[j], 1 - Sp[j]))
  Y
}

GRIDS <- list(A = build_grid_A(), B = build_grid_B(), C = build_grid_C())
GENS  <- list(A = gen_panel_A,    B = gen_panel_B,    C = gen_panel_C)
OFFSET <- c(A = 0L, B = 1000L, C = 2000L)   # keeps per-arm seeds globally unique

if (SMOKE)
  GRIDS <- lapply(GRIDS, function(g)
    g[unique(c(1L, nrow(g) %/% 2L, nrow(g))), , drop = FALSE])

# ---- global manifest (arm, cell_id, seed) + two-machine split -------------
manifest <- do.call(rbind, lapply(names(GRIDS), function(arm)
  data.frame(arm = arm, cell_id = GRIDS[[arm]]$cell_id,
             seed = MASTER_SEED + OFFSET[[arm]] + GRIDS[[arm]]$cell_id,
             stringsAsFactors = FALSE)))
manifest$job_id <- seq_len(nrow(manifest))
mod <- as.integer(Sys.getenv("SPLIT_MOD", "1"))
rem <- as.integer(Sys.getenv("SPLIT_REM", "0"))
manifest <- manifest[manifest$job_id %% mod == rem, , drop = FALSE]

cell_file <- function(arm, id) file.path(PER_CELL, sprintf("arm%s_cell_%03d.rds", arm, id))

# ---- resume: skip cells whose RDS exists, is READABLE, and is COMPLETE ----
# complete means complete AT THIS RUN'S TARGET (guards stale SMOKE files).
is_complete <- function(arm, id) {
  f <- cell_file(arm, id)
  if (!file.exists(f)) return(FALSE)
  ok <- tryCatch({
    x <- readRDS(f)
    is.list(x) && !is.null(x$draws) && is.data.frame(x$draws) &&
      !is.null(x$target_draws) && nrow(x$draws) == x$target_draws &&
      x$target_draws == reps
  }, error = function(e) FALSE)
  isTRUE(ok)
}
keep <- !mapply(is_complete, manifest$arm, manifest$cell_id)
todo <- manifest[keep, , drop = FALSE]
cat(sprintf("tier2B | jobs(split %d/%d): %d | done: %d | todo: %d | reps: %d\n",
            rem, mod, nrow(manifest), sum(!keep), nrow(todo), reps))

# ---- one rep -> c(delta, obs x3, implied_q x3, percentile_pp x3) ----------
# Order: delta, pabak, fleiss_kappa, mean_ac1,
#        q_pabak, q_fleiss_kappa, q_mean_ac1,
#        pp_pabak, pp_fleiss_kappa, pp_mean_ac1
COEFS <- c("pabak", "fleiss_kappa", "mean_ac1")
NA10  <- rep(NA_real_, 10L)
extract_rep <- function(Y) {
  res <- tryCatch(suppressWarnings(suppressMessages(check_asymmetry(Y))),
                  error = function(e) NULL)
  if (is.null(res)) return(NA10)
  p <- res$panel
  val <- function(nm, col) { i <- which(p$coefficient == nm)
    if (length(i) == 1L) p[[col]][i] else NA_real_ }
  c(as.numeric(res$delta_hat),
    val("pabak", "observed"),      val("fleiss_kappa", "observed"),      val("mean_ac1", "observed"),
    val("pabak", "implied_q"),     val("fleiss_kappa", "implied_q"),     val("mean_ac1", "implied_q"),
    val("pabak", "percentile_pp"), val("fleiss_kappa", "percentile_pp"), val("mean_ac1", "percentile_pp"))
}

run_cell <- function(arm, cell_id, seed) {
  row <- GRIDS[[arm]][GRIDS[[arm]]$cell_id == cell_id, , drop = FALSE]
  gen <- GENS[[arm]]
  set.seed(seed)                    # ONCE, before any RNG for this cell
  M <- matrix(NA_real_, reps, 10L)
  for (r in seq_len(reps)) M[r, ] <- extract_rep(gen(row))
  draws <- data.frame(
    rep             = seq_len(reps),
    delta           = M[, 1L],
    pabak           = M[, 2L], fleiss_kappa   = M[, 3L], mean_ac1        = M[, 4L],
    q_pabak         = M[, 5L], q_fleiss_kappa = M[, 6L], q_mean_ac1      = M[, 7L],
    pp_pabak        = M[, 8L], pp_fleiss_kappa = M[, 9L], pp_mean_ac1    = M[, 10L])
  list(arm = arm, cell = row, draws = draws, target_draws = reps,
       n_finite = sum(is.finite(M[, 1L])), master_seed = MASTER_SEED,
       seed = seed, r_version = R.version.string)
}

# atomic per-cell write (temp + rename) so a crash never leaves a partial RDS
save_cell <- function(obj, arm, id) {
  f <- cell_file(arm, id); tmp <- paste0(f, ".tmp", Sys.getpid())
  saveRDS(obj, tmp); file.rename(tmp, f)
}

# ---- parallel over cells (each cell fully self-seeded) --------------------
n_workers <- max(1L, min(as.integer(Sys.getenv("WORKERS", "8")),
                         parallel::detectCores() - 1L))
plan(multicore, workers = n_workers)
cat(sprintf("[%s] tier2B start: %d todo, %d workers, split %d/%d, seed %d, %s\n",
            format(Sys.time()), nrow(todo), n_workers, rem, mod, MASTER_SEED,
            R.version.string),
    file = PROGRESS, append = file.exists(PROGRESS))

t0 <- Sys.time()
if (nrow(todo) > 0L)
  invisible(future_lapply(seq_len(nrow(todo)), function(i) {
    row <- todo[i, ]
    save_cell(run_cell(row$arm, row$cell_id, row$seed), row$arm, row$cell_id)
    if (i %% 20L == 0L)
      cat(sprintf("[%s] %d/%d (arm%s cell %d)\n",
                  format(Sys.time()), i, nrow(todo), row$arm, row$cell_id),
          file = PROGRESS, append = TRUE)
    NULL
  }, future.seed = TRUE))

cat(sprintf("[%s] tier2B split %d/%d COMPLETE. Wall: %.1f min.\n",
            format(Sys.time()), rem, mod, as.numeric(Sys.time() - t0, units = "mins")),
    file = PROGRESS, append = TRUE)
writeLines(sprintf("tier2B split %d/%d complete %s", rem, mod, format(Sys.time())),
           file.path(TIER2DB, sprintf("split_%d_of_%d.DONE", rem, mod)))
cat(sprintf("Done. Wall: %.1f min.\n", as.numeric(Sys.time() - t0, units = "mins")))
