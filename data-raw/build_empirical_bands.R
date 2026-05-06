# build_empirical_bands.R — Precompute the empirical-band sysdata asset for
# grass::position_on_surface(method = "empirical").
#
# For every (metric, F_key, k, N, q_true) cell in the multirater_sim_v3 grid
# (6,240 scenarios x 5 metrics = 31,200 cells), invert each of the 2,000 per-
# rep observed metric values through the closed-form reference curve at that
# (F_key, k) to obtain a 2,000-length empirical q_hat_rep vector. Compress
# each vector to a quantile summary and persist into grass/R/sysdata.rda
# alongside the existing `reference_binary` object.
#
# Compression choice: OPTION (A) -- fine empirical quantile summary.
# We store quantiles at 13 probabilities
#   {0.01, 0.025, 0.05, 0.10, 0.25, 0.375, 0.5, 0.625, 0.75, 0.90, 0.95,
#    0.975, 0.99}
# per cell. This lets position_on_surface() reconstruct P(q_hat_rep <= t)
# for any user-specified band boundary t via piecewise-linear interpolation
# on the empirical CDF, while preserving tail behaviour that a 5-point
# summary (default in the spec) collapses. Size: 31,200 x 13 ~ 405,600 numeric
# cells ~ 3.2 MB uncompressed, ~1.5 MB on disk after xz.
#
# Run from PABAK_Investigation project root:
#   Rscript grass/data-raw/build_empirical_bands.R
#
# Inputs:
#   paper2/simulation_output/multirater_sim_v3/multirater_sim_v3_full.rds
#   paper2/code/04_reference_closed_form.R
#   paper2/code/06_grid.R
# Outputs:
#   grass/R/sysdata.rda  (updated -- existing objects preserved)

suppressMessages({
  source("paper2/code/04_reference_closed_form.R")
  source("paper2/code/06_grid.R")
})

FULL_PATH  <- "paper2/simulation_output/multirater_sim_v3/multirater_sim_v3_combined_full.rds"
# Swapped 2026-04-23 from multirater_sim_v3_full.rds to the combined
# (v3 + v3_k25) full.rds so the empirical bands cover k in
# {3,5,8,15,25} rather than just {3,5,8,15}.
SYSDATA    <- "grass/R/sysdata.rda"

Q_GRID <- seq(0.5, 1.0, length.out = 501L)
QUANTILE_PROBS <- c(0.01, 0.025, 0.05, 0.10, 0.25, 0.375, 0.5,
                    0.625, 0.75, 0.90, 0.95, 0.975, 0.99)

# Map from per_rep column names (sim output) to the exported metric names we
# use in grass::position_on_surface(). We drop the (fitted) logit_mixed_icc
# column in favour of the oracle version because that is what q_recovery.rds
# uses; fitted ICC has its own separate recovery table.
PER_REP_COL_FOR <- c(
  pabak          = "mean_pabak",
  fleiss_kappa   = "fleiss_kappa",
  mean_ac1       = "mean_ac1",
  krippendorff_a = "krippendorff_a",
  icc            = "logit_mixed_icc_oracle"
)

# Map from metric name to the CF-panel field name (paper2/code/04_... uses
# different naming again).
CF_FIELD_FOR <- c(
  pabak          = "mean_pabak",
  fleiss_kappa   = "fleiss_kappa",
  mean_ac1       = "mean_ac1",
  krippendorff_a = "krippendorff_a",
  icc            = "logit_mixed_icc_oracle"
)

# ---- Load existing sysdata so we can add rather than replace ---------------
existing_env <- new.env()
if (file.exists(SYSDATA)) {
  load(SYSDATA, envir = existing_env)
  cat("Existing sysdata objects:", paste(ls(existing_env), collapse = ", "),
      "\n")
} else {
  cat("No existing sysdata.rda -- creating fresh.\n")
}

# ---- Build closed-form lookup tables (501 x 5) per (F_key, k) --------------
cat("Loading full sim (~313 MB, 10-20s)...\n")
t0 <- Sys.time()
full <- readRDS(FULL_PATH)
cat(sprintf("  loaded in %.1f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))

s <- full$summaries
per_rep <- full$per_rep

uniq_Fk <- unique(s[, c("F_family", "F_key", "k")])
cat(sprintf("Building %d closed-form lookup tables...\n", nrow(uniq_Fk)))
t0 <- Sys.time()
lookups <- vector("list", nrow(uniq_Fk))
names(lookups) <- sprintf("%s|%s|%d", uniq_Fk$F_family, uniq_Fk$F_key,
                          uniq_Fk$k)
for (i in seq_len(nrow(uniq_Fk))) {
  row <- uniq_Fk[i, ]
  fp <- F_params_for(row$F_key)
  M <- matrix(NA_real_, nrow = length(Q_GRID), ncol = length(CF_FIELD_FOR),
              dimnames = list(NULL, names(CF_FIELD_FOR)))
  for (qi in seq_along(Q_GRID)) {
    cfp <- cf_panel(q = Q_GRID[qi], F_family = row$F_family,
                    F_params = fp, k = row$k)
    for (m in names(CF_FIELD_FOR)) M[qi, m] <- cfp[[CF_FIELD_FOR[[m]]]]
  }
  lookups[[i]] <- M
}
cat(sprintf("  done in %.1f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))

# ---- Invert each per-rep metric value to q_hat_rep, summarise --------------
#
# For a monotone-increasing reference curve on q in [0.5, 1] we invert via
# stats::approx. If the reference curve is flat or non-monotone (e.g. alpha
# at extreme prevalence), we fall back to clamping to the nearest q on the
# grid. Below-min obs -> q=0.5; above-max obs -> q=1.0.
invert_vec <- function(obs_vec, ref_col) {
  if (all(is.na(ref_col))) return(rep(NA_real_, length(obs_vec)))
  # Sort by ref_col for stable monotone interp; then undo nothing (we just
  # want q_hat for each obs, not ordered).
  keep <- is.finite(ref_col)
  if (sum(keep) < 2L) return(rep(NA_real_, length(obs_vec)))
  rc <- ref_col[keep]
  qg <- Q_GRID[keep]
  ord <- order(rc)
  rc_ord <- rc[ord]
  qg_ord <- qg[ord]
  rng <- range(rc_ord)
  q_hat <- rep(NA_real_, length(obs_vec))
  ok <- is.finite(obs_vec)
  # Below-min -> q at min (which is q=0.5 for the monotone-increasing case)
  below <- ok & obs_vec <= rng[1]
  above <- ok & obs_vec >= rng[2]
  mid   <- ok & !below & !above
  q_hat[below] <- qg_ord[1L]
  q_hat[above] <- qg_ord[length(qg_ord)]
  if (any(mid)) {
    q_hat[mid] <- stats::approx(x = rc_ord, y = qg_ord,
                                xout = obs_vec[mid],
                                rule = 2, ties = "ordered")$y
  }
  q_hat
}

cat("Inverting per-rep metrics (6,240 scenarios x 5 metrics x 2,000 reps)...\n")
t0 <- Sys.time()

n_scn <- nrow(s)
n_met <- length(PER_REP_COL_FOR)
n_q   <- length(QUANTILE_PROBS)

# Quantile tensor: [scenario_id, metric, quantile_probs].
# We keep scenario order aligned to s$scenario_id so lookup by index is O(1).
q_hat_quantiles <- array(
  NA_real_,
  dim = c(n_scn, n_met, n_q),
  dimnames = list(scenario_id = as.character(s$scenario_id),
                  metric      = names(PER_REP_COL_FOR),
                  prob        = as.character(QUANTILE_PROBS))
)

# Also capture mean + sd for diagnostic convenience (tiny: 31,200 x 2).
q_hat_mean <- matrix(NA_real_, nrow = n_scn, ncol = n_met,
                     dimnames = list(NULL, names(PER_REP_COL_FOR)))
q_hat_sd   <- q_hat_mean

for (i in seq_len(n_scn)) {
  row <- s[i, ]
  key <- sprintf("%s|%s|%d", row$F_family, row$F_key, row$k)
  M <- lookups[[key]]
  reps <- per_rep[[as.character(row$scenario_id)]]
  if (is.null(reps)) {
    # fallback: match by index position
    reps <- per_rep[[i]]
  }
  if (is.null(reps)) next
  for (mi in seq_along(PER_REP_COL_FOR)) {
    metric_name <- names(PER_REP_COL_FOR)[mi]
    col_name <- PER_REP_COL_FOR[[mi]]
    obs_vec <- reps[, col_name]
    q_hat_rep <- invert_vec(obs_vec, M[, metric_name])
    q_hat_rep <- q_hat_rep[is.finite(q_hat_rep)]
    if (length(q_hat_rep) >= 10L) {
      q_hat_quantiles[i, mi, ] <- stats::quantile(
        q_hat_rep, probs = QUANTILE_PROBS, names = FALSE, type = 7
      )
      q_hat_mean[i, mi] <- mean(q_hat_rep)
      q_hat_sd[i, mi]   <- stats::sd(q_hat_rep)
    }
  }
  if (i %% 500L == 0L) {
    cat(sprintf("  scenario %d / %d  (%.1f s elapsed)\n",
                i, n_scn, as.numeric(Sys.time() - t0, units = "secs")))
  }
}
cat(sprintf("  done in %.1f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))

# ---- Index table: one row per scenario for fast NN lookup ------------------
# Columns: scenario_id (int), F_key (chr), q_true, k, N, M1.
# We merge M1 from the F_params registry the same way 12_q_inversion does.
compute_M1 <- function(F_family, F_key) {
  fp <- F_params_for(F_key)
  if (F_family == "logit_normal") {
    e_p_logit_normal(mu = fp$mu, tau2 = fp$tau2)
  } else {
    sum(fp$w * fp$p)
  }
}
M1_by_Fkey <- vapply(
  seq_len(nrow(uniq_Fk)),
  function(i) compute_M1(uniq_Fk$F_family[i], uniq_Fk$F_key[i]),
  numeric(1L)
)
M1_map <- data.frame(F_family = uniq_Fk$F_family, F_key = uniq_Fk$F_key,
                     M1 = M1_by_Fkey, stringsAsFactors = FALSE)
# Dedupe (same F_key gets same M1 irrespective of k)
M1_map <- unique(M1_map[, c("F_key", "M1")])

empirical_band_index <- data.frame(
  scenario_id = s$scenario_id,
  F_family    = s$F_family,
  F_key       = s$F_key,
  q_true      = s$q,
  k           = s$k,
  N           = s$N,
  stringsAsFactors = FALSE
)
empirical_band_index <- merge(empirical_band_index,
                              M1_map, by = "F_key", all.x = TRUE, sort = FALSE)
# Restore scenario_id order so index matches the quantile array's first dim.
empirical_band_index <- empirical_band_index[
  order(empirical_band_index$scenario_id), ]
stopifnot(identical(as.integer(empirical_band_index$scenario_id),
                    as.integer(s$scenario_id[order(s$scenario_id)])))
rownames(empirical_band_index) <- NULL

# Reorder the quantile array to match the sorted scenario_id order as well,
# to guarantee index == scenario_id row-by-row.
ord <- order(as.integer(dimnames(q_hat_quantiles)$scenario_id))
q_hat_quantiles <- q_hat_quantiles[ord, , , drop = FALSE]
q_hat_mean      <- q_hat_mean[ord, , drop = FALSE]
q_hat_sd        <- q_hat_sd[ord, , drop = FALSE]

# ---- Bundle into a single named list for sysdata ---------------------------
empirical_q_hat_surface <- list(
  probs     = QUANTILE_PROBS,
  quantiles = q_hat_quantiles,  # [scenario_id, metric, prob]
  mean      = q_hat_mean,       # [scenario_id, metric]
  sd        = q_hat_sd,         # [scenario_id, metric]
  index     = empirical_band_index,  # scenario_id, F_family, F_key, q_true, k, N, M1
  metrics   = names(PER_REP_COL_FOR),
  source    = sprintf("%s (%s scenarios x 2,000 reps; k in {%s})",
                      FULL_PATH,
                      format(nrow(empirical_band_index), big.mark = ","),
                      paste(sort(unique(empirical_band_index$k)), collapse = ",")),
  built_on  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

# ---- Persist, PRESERVING existing sysdata objects --------------------------
new_env <- new.env()
# Copy all existing objects
for (n in ls(existing_env)) {
  assign(n, get(n, envir = existing_env), envir = new_env)
}
# Add ours (and overwrite if it exists)
assign("empirical_q_hat_surface", empirical_q_hat_surface, envir = new_env)

cat("sysdata objects (new):",
    paste(ls(new_env), collapse = ", "), "\n")

save(list = ls(new_env), file = SYSDATA, envir = new_env,
     compress = "xz", compression_level = 9)

info <- file.info(SYSDATA)
cat(sprintf("Wrote %s (%.2f KB)\n", SYSDATA, info$size / 1024))
cat("DONE.\n")
