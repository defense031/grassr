#!/usr/bin/env Rscript
# Stage 1 — surface N-densification + q-unification (grassr v0.7.0 Tier 1).
#
# Adds N in {15, 20, 30, 75, 100, 150, 300, 500} to the reference-surface
# grid and unifies the q-grid across k to union(per-k existing set,
# {0.92, 0.94, 0.97}). The grid is built by anti-join against the SHIPPED
# empirical_q_hat_surface index, so existing cells are never re-simulated
# and new cells inherit the exact (F_key, q, k) structure of production.
#
# Mirrors run_q094_k2_k25_augmentation.R (v0.5.2 augmentation pattern):
# per-scenario RDS with resume, future_lapply, oracle ICC (glmer disabled).
#
# Env:
#   GRASS_SIM_ROOT  bundle root containing paper2/code and this program
#                   (default: the PABAK_Investigation repo root)
#   SPLIT_MOD, SPLIT_REM  optional worker split: run scenarios where
#                   scenario_id %% SPLIT_MOD == SPLIT_REM
#   SMOKE=1         4 scenarios at 50 reps (pipeline check only)
#
# Seeds: master_seed 20260704 + scenario_id (per design doc).

suppressPackageStartupMessages({
  library(future.apply)
})

ROOT <- Sys.getenv("GRASS_SIM_ROOT",
                   "/Users/austinsemmel/Desktop/PABAK_Investigation")
setwd(ROOT)
PROG     <- file.path(ROOT, "grassr/simulation/v070_program/tier1/stage1_surface_ndense")
PER_SCEN <- file.path(PROG, "per_scenario")
PROGRESS <- file.path(PROG, "PROGRESS.txt")
dir.create(PER_SCEN, recursive = TRUE, showWarnings = FALSE)

source("paper2/code/01_dgp.R")
source("paper2/code/02_metrics_observed.R")
source("paper2/code/04_reference_closed_form.R")
source("paper2/code/06_grid.R")

SMOKE <- identical(Sys.getenv("SMOKE"), "1")

# ---- target grid: anti-join against the shipped surface index ----
e <- new.env()
load(file.path(ROOT, "grassr/R/sysdata.rda"), envir = e)
idx <- e$empirical_q_hat_surface$index   # F_key, scenario_id, F_family, q_true, k, N, M1

N_ALL <- c(15L, 20L, 30L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L)
Q_AUG <- c(0.92, 0.94, 0.97)

target <- do.call(rbind, lapply(sort(unique(idx$k)), function(kk) {
  q_target <- sort(union(unique(idx$q_true[idx$k == kk]), Q_AUG))
  expand.grid(k = kk, q = q_target, N = N_ALL,
              F_key = unique(idx$F_key),
              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}))
have <- paste(idx$F_key, idx$q_true, idx$k, idx$N)
target <- target[!(paste(target$F_key, target$q, target$k, target$N) %in% have), ]

# F_family from the registry key prefix; reps fixed at 2000
target$F_family <- ifelse(grepl("^LN_", target$F_key),
                          "logit_normal", "discrete_mixture")
target$reps <- 2000L
target <- target[order(target$k, target$N, target$q, target$F_key), ]
target$scenario_id <- seq_len(nrow(target))

if (SMOKE) {
  target <- target[c(1L, nrow(target) %/% 3L, 2L * nrow(target) %/% 3L, nrow(target)), ]
  target$reps <- 50L
}

# optional cross-machine split
mod <- as.integer(Sys.getenv("SPLIT_MOD", "1"))
rem <- as.integer(Sys.getenv("SPLIT_REM", "0"))
target <- target[target$scenario_id %% mod == rem, , drop = FALSE]

existing <- list.files(PER_SCEN, pattern = "^scenario_[0-9]+\\.rds$")
done_ids <- as.integer(sub("^scenario_([0-9]+)\\.rds$", "\\1", existing))
todo <- target[!target$scenario_id %in% done_ids, , drop = FALSE]

cat(sprintf("stage1 | target(this split): %d | done: %d | todo: %d\n",
            nrow(target), sum(target$scenario_id %in% done_ids), nrow(todo)))

run_scenario <- function(row, master_seed = 20260704L) {
  sid <- row$scenario_id
  set.seed(master_seed + sid)
  q <- row$q; k <- row$k; N <- row$N; R <- row$reps
  F_params <- F_params_for(row$F_key)

  metric_names <- c("fleiss_kappa", "mean_pabak", "mean_ac1",
                    "krippendorff_a", "logit_mixed_icc",
                    "logit_mixed_icc_oracle")
  M <- matrix(NA_real_, nrow = R, ncol = length(metric_names),
              dimnames = list(NULL, metric_names))
  for (r in seq_len(R)) {
    sim <- simulate_dataset(N = N, k = k, q = q,
                            F_family = row$F_family, F_params = F_params)
    Y <- sim$Y
    M[r, "fleiss_kappa"]   <- obs_fleiss_kappa(Y)
    M[r, "mean_pabak"]     <- obs_mean_pairwise_pabak(Y)
    M[r, "mean_ac1"]       <- obs_mean_pairwise_ac1(Y)
    M[r, "krippendorff_a"] <- obs_krippendorff_alpha(Y)
    M[r, "logit_mixed_icc_oracle"] <- obs_logit_mixed_icc_oracle(sim)
  }
  means <- colMeans(M, na.rm = TRUE)
  sds   <- apply(M, 2, function(x) sd(x, na.rm = TRUE))
  n_ok  <- colSums(!is.na(M))
  mcse  <- sds / sqrt(pmax(n_ok, 1))
  list(
    summary = data.frame(
      scenario_id = sid,
      F_family = row$F_family, F_key = row$F_key,
      q = q, k = k, N = N, reps = R,
      glmer_fitted = FALSE,
      mean_fleiss_kappa    = means["fleiss_kappa"],
      mean_mean_pabak      = means["mean_pabak"],
      mean_mean_ac1        = means["mean_ac1"],
      mean_krippendorff_a  = means["krippendorff_a"],
      mean_icc_glmm        = means["logit_mixed_icc"],
      mean_icc_oracle      = means["logit_mixed_icc_oracle"],
      mcse_fleiss_kappa    = mcse["fleiss_kappa"],
      mcse_mean_pabak      = mcse["mean_pabak"],
      mcse_mean_ac1        = mcse["mean_ac1"],
      mcse_krippendorff_a  = mcse["krippendorff_a"],
      mcse_icc_glmm        = mcse["logit_mixed_icc"],
      mcse_icc_oracle      = mcse["logit_mixed_icc_oracle"],
      n_ok_fleiss_kappa    = n_ok["fleiss_kappa"],
      n_ok_icc_glmm        = n_ok["logit_mixed_icc"],
      row.names = NULL, stringsAsFactors = FALSE),
    per_rep = M,
    r_version = R.version.string)
}

n_workers <- max(1L, min(as.integer(Sys.getenv("WORKERS", "10")),
                         parallel::detectCores() - 2L))
plan(multisession, workers = n_workers)
cat(sprintf("Workers: %d | split %d/%d | R %s\n",
            n_workers, rem, mod, R.version.string))

cat(sprintf("[%s] stage1 start: %d todo, %d workers, split %d/%d, %s\n",
            format(Sys.time()), nrow(todo), n_workers, rem, mod,
            R.version.string),
    file = PROGRESS, append = file.exists(PROGRESS))

t0 <- Sys.time()
invisible(future_lapply(seq_len(nrow(todo)), function(i) {
  row <- todo[i, ]
  res <- run_scenario(row)
  saveRDS(res, file.path(PER_SCEN,
                         sprintf("scenario_%05d.rds", row$scenario_id)))
  if (i %% 50L == 0L)
    cat(sprintf("[%s] %d/%d done (last sid=%05d k=%d N=%d q=%.2f)\n",
                format(Sys.time()), i, nrow(todo), row$scenario_id,
                row$k, row$N, row$q),
        file = PROGRESS, append = TRUE)
  NULL
}, future.seed = TRUE))
t1 <- Sys.time()

cat(sprintf("[%s] stage1 split %d/%d COMPLETE. Wall: %.1f min.\n",
            format(Sys.time()), rem, mod,
            as.numeric(t1 - t0, units = "mins")),
    file = PROGRESS, append = TRUE)
cat(sprintf("Done. Wall: %.1f min.\n", as.numeric(t1 - t0, units = "mins")))
