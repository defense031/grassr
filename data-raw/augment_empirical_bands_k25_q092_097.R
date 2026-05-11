# augment_empirical_bands_k25_q092_097.R --- APPEND q in {0.92, 0.97}
# entries to the existing k=25 portion of empirical_q_hat_surface in
# grass/R/sysdata.rda.
#
# Background: the k=25 surface (post 2026-05-10 fix) has q_true grid
# {0.55, 0.60, ..., 0.95, 0.99}. Panels with q_hat in (0.92, 0.94) or
# (0.96, 0.98) snap to a q_true cohort whose extremely tight empirical
# distribution (k=25 N=1000 has q-hat sd ~0.003) does not bracket
# the panel's q_hat, clamping the percentile to 0 or 100. Adding
# q_true in {0.92, 0.97} closes the snap gap.
#
# This script is additive: it reads the current sysdata, finds k=25
# rows, appends 312 new k=25 q in {0.92, 0.97} rows, and writes back.
# It does not touch the k=2/3/5/8/15 portion of the surface.
#
# Run from PABAK_Investigation project root (after sim completes):
#   Rscript grass/data-raw/augment_empirical_bands_k25_q092_097.R
#
# Inputs:
#   paper1_2_merged/output/multirater_sim_k25_q092_097_aug/summaries.rds
#       (312 x 22)
#   paper1_2_merged/output/multirater_sim_k25_q092_097_aug/per_scenario/scenario_NNNN.rds
#   paper2/code/04_reference_closed_form.R  (cf_panel)
#   paper2/code/06_grid.R                   (F_params_for, e_p_logit_normal)
# Outputs:
#   grass/R/sysdata.rda  (updated --- existing entries preserved)

suppressMessages({
  source("paper2/code/04_reference_closed_form.R")
  source("paper2/code/06_grid.R")
})

AUG_DIR <- "paper1_2_merged/output/multirater_sim_k25_q092_097_aug"
SYSDATA <- "grass/R/sysdata.rda"

Q_GRID <- seq(0.5, 1.0, length.out = 501L)
QUANTILE_PROBS <- c(0.01, 0.025, 0.05, 0.10, 0.25, 0.375, 0.5,
                    0.625, 0.75, 0.90, 0.95, 0.975, 0.99)
PER_REP_COL_FOR <- c(
  pabak          = "mean_pabak",
  fleiss_kappa   = "fleiss_kappa",
  mean_ac1       = "mean_ac1",
  krippendorff_a = "krippendorff_a",
  icc            = "logit_mixed_icc_oracle"
)
CF_FIELD_FOR <- PER_REP_COL_FOR

# Offset: +220000 to keep distinct from base k=2 (+200000),
# k25 base (+100000), and k=2 q=0.97 augment (+210000).
AUG_SCENARIO_ID_OFFSET <- 220000L

# ---- Load existing sysdata ------------------------------------------------
existing_env <- new.env()
stopifnot(file.exists(SYSDATA))
load(SYSDATA, envir = existing_env)
cat("Existing sysdata objects:", paste(ls(existing_env), collapse = ", "), "\n")
existing_eqs <- existing_env$empirical_q_hat_surface
cat(sprintf("Existing surface: %d scenarios, k in {%s}\n",
            nrow(existing_eqs$index),
            paste(sort(unique(existing_eqs$index$k)), collapse = ",")))
existing_k25 <- existing_eqs$index[existing_eqs$index$k == 25L, ]
cat(sprintf("  k=25 q-grid: {%s}\n",
            paste(sort(unique(existing_k25$q_true)), collapse = ",")))
new_qs <- c(0.92, 0.97)
already <- new_qs[new_qs %in% existing_k25$q_true]
if (length(already) > 0L) {
  cat(sprintf("WARNING: q in {%s} already present at k=25; existing rows will be replaced.\n",
              paste(already, collapse = ",")))
}

# ---- Defensive precondition checks ---------------------------------------
stopifnot(file.exists(file.path(AUG_DIR, "summaries.rds")))

# ---- Load summaries + per_rep ---------------------------------------------
cat("Loading q in {0.92, 0.97} k=25 augmentation summaries...\n")
s <- readRDS(file.path(AUG_DIR, "summaries.rds"))
cat(sprintf("  %d scenarios, k=%d, q in {%s}, N in {%s}, F-keys=%d\n",
            nrow(s), unique(s$k),
            paste(sort(unique(s$q)), collapse = ","),
            paste(sort(unique(s$N)), collapse = ","),
            length(unique(s$F_key))))
stopifnot(unique(s$k) == 25L)
stopifnot(setequal(unique(s$q), c(0.92, 0.97)))
stopifnot(nrow(s) == 312L)

cat(sprintf("Loading %d per-scenario rds files...\n", nrow(s)))
per_rep <- vector("list", nrow(s))
for (i in seq_len(nrow(s))) {
  sid <- s$scenario_id[i]
  f <- file.path(AUG_DIR, "per_scenario",
                 sprintf("scenario_%04d.rds", sid))
  per_rep[[i]] <- readRDS(f)$per_rep
}
stopifnot(length(per_rep) == nrow(s))

s$scenario_id <- s$scenario_id + AUG_SCENARIO_ID_OFFSET
names(per_rep) <- as.character(s$scenario_id)

# Defensive: name and order must match (would catch a per_rep ordering bug).
stopifnot(setequal(names(per_rep), as.character(s$scenario_id)))
stopifnot(!any(duplicated(names(per_rep))))

# ---- Closed-form lookup tables --------------------------------------------
uniq_Fk <- unique(s[, c("F_family", "F_key", "k")])
cat(sprintf("Building %d closed-form lookup tables for k=25...\n", nrow(uniq_Fk)))
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

# ---- Inversion ------------------------------------------------------------
invert_vec <- function(obs_vec, ref_col) {
  if (all(is.na(ref_col))) return(rep(NA_real_, length(obs_vec)))
  keep <- is.finite(ref_col)
  if (sum(keep) < 2L) return(rep(NA_real_, length(obs_vec)))
  rc <- ref_col[keep]; qg <- Q_GRID[keep]
  ord <- order(rc); rc_ord <- rc[ord]; qg_ord <- qg[ord]
  rng <- range(rc_ord)
  q_hat <- rep(NA_real_, length(obs_vec))
  ok <- is.finite(obs_vec)
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

cat(sprintf("Inverting per-rep metrics (%d scenarios x %d metrics)...\n",
            nrow(s), length(PER_REP_COL_FOR)))
n_scn <- nrow(s); n_met <- length(PER_REP_COL_FOR); n_q <- length(QUANTILE_PROBS)
q_hat_quantiles <- array(NA_real_, dim = c(n_scn, n_met, n_q),
  dimnames = list(scenario_id = as.character(s$scenario_id),
                  metric      = names(PER_REP_COL_FOR),
                  prob        = as.character(QUANTILE_PROBS)))
q_hat_mean <- matrix(NA_real_, nrow = n_scn, ncol = n_met,
                     dimnames = list(NULL, names(PER_REP_COL_FOR)))
q_hat_sd   <- q_hat_mean

for (i in seq_len(n_scn)) {
  row <- s[i, ]
  key <- sprintf("%s|%s|%d", row$F_family, row$F_key, row$k)
  M <- lookups[[key]]
  reps <- per_rep[[as.character(row$scenario_id)]]
  if (is.null(reps)) next
  for (mi in seq_along(PER_REP_COL_FOR)) {
    metric_name <- names(PER_REP_COL_FOR)[mi]
    col_name <- PER_REP_COL_FOR[[mi]]
    obs_vec <- reps[, col_name]
    q_hat_rep <- invert_vec(obs_vec, M[, metric_name])
    q_hat_rep <- q_hat_rep[is.finite(q_hat_rep)]
    if (length(q_hat_rep) >= 10L) {
      q_hat_quantiles[i, mi, ] <- stats::quantile(
        q_hat_rep, probs = QUANTILE_PROBS, names = FALSE, type = 7)
      q_hat_mean[i, mi] <- mean(q_hat_rep)
      q_hat_sd[i, mi]   <- stats::sd(q_hat_rep)
    }
  }
}

compute_M1 <- function(F_family, F_key) {
  fp <- F_params_for(F_key)
  if (F_family == "logit_normal") {
    e_p_logit_normal(mu = fp$mu, tau2 = fp$tau2)
  } else {
    sum(fp$w * fp$p)
  }
}
M1_by_Fkey <- vapply(seq_len(nrow(uniq_Fk)),
                     function(i) compute_M1(uniq_Fk$F_family[i], uniq_Fk$F_key[i]),
                     numeric(1L))
M1_map <- unique(data.frame(F_key = uniq_Fk$F_key, M1 = M1_by_Fkey,
                            stringsAsFactors = FALSE))

new_index <- data.frame(
  scenario_id = s$scenario_id,
  F_family    = s$F_family,
  F_key       = s$F_key,
  q_true      = s$q,
  k           = s$k,
  N           = s$N,
  stringsAsFactors = FALSE
)
new_index <- merge(new_index, M1_map, by = "F_key", all.x = TRUE, sort = FALSE)
new_index <- new_index[order(new_index$scenario_id), ]
rownames(new_index) <- NULL

ord <- order(as.integer(dimnames(q_hat_quantiles)$scenario_id))
q_hat_quantiles <- q_hat_quantiles[ord, , , drop = FALSE]
q_hat_mean      <- q_hat_mean[ord, , drop = FALSE]
q_hat_sd        <- q_hat_sd[ord, , drop = FALSE]

# ---- Merge into existing surface ------------------------------------------
existing_idx <- existing_eqs$index
existing_quantiles <- existing_eqs$quantiles
existing_mean <- existing_eqs$mean
existing_sd <- existing_eqs$sd

keep_existing <- !(existing_idx$k == 25L & existing_idx$q_true %in% new_qs)
n_dropped <- sum(!keep_existing)
if (n_dropped > 0L) {
  cat(sprintf("Dropping %d existing k=25 q in {0.92, 0.97} rows\n", n_dropped))
}
existing_idx <- existing_idx[keep_existing, , drop = FALSE]
existing_quantiles <- existing_quantiles[keep_existing, , , drop = FALSE]
existing_mean <- existing_mean[keep_existing, , drop = FALSE]
existing_sd <- existing_sd[keep_existing, , drop = FALSE]

merged_index <- rbind(existing_idx, new_index)
merged_quantiles <- array(NA_real_,
  dim = c(nrow(existing_quantiles) + nrow(q_hat_quantiles),
          dim(existing_quantiles)[2], dim(existing_quantiles)[3]),
  dimnames = list(
    scenario_id = c(dimnames(existing_quantiles)$scenario_id,
                    dimnames(q_hat_quantiles)$scenario_id),
    metric = dimnames(existing_quantiles)$metric,
    prob   = dimnames(existing_quantiles)$prob))
merged_quantiles[seq_len(nrow(existing_quantiles)), , ] <- existing_quantiles
merged_quantiles[seq.int(nrow(existing_quantiles) + 1L,
                          dim(merged_quantiles)[1]), , ] <- q_hat_quantiles
merged_mean <- rbind(existing_mean, q_hat_mean)
merged_sd   <- rbind(existing_sd,   q_hat_sd)

ord_all <- order(merged_index$scenario_id)
merged_index <- merged_index[ord_all, , drop = FALSE]
rownames(merged_index) <- NULL
merged_quantiles <- merged_quantiles[ord_all, , , drop = FALSE]
merged_mean <- merged_mean[ord_all, , drop = FALSE]
merged_sd   <- merged_sd[ord_all, , drop = FALSE]
dimnames(merged_quantiles)$scenario_id <- as.character(merged_index$scenario_id)

empirical_q_hat_surface <- list(
  probs     = QUANTILE_PROBS,
  quantiles = merged_quantiles,
  mean      = merged_mean,
  sd        = merged_sd,
  index     = merged_index,
  metrics   = names(PER_REP_COL_FOR),
  source    = sprintf(
    "multirater_sim_v3_combined_full.rds + multirater_sim_k2_extension + multirater_sim_k2_q097_aug + multirater_sim_k25_q092_097_aug (%s scenarios; k in {%s})",
    format(nrow(merged_index), big.mark = ","),
    paste(sort(unique(merged_index$k)), collapse = ",")),
  built_on  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

cat(sprintf("Final empirical_q_hat_surface: %d scenarios, k in {%s}\n",
            nrow(merged_index),
            paste(sort(unique(merged_index$k)), collapse = ",")))
new_k25 <- merged_index[merged_index$k == 25L, ]
cat(sprintf("  k=25 q-grid: {%s} (%d scenarios)\n",
            paste(sort(unique(new_k25$q_true)), collapse = ","), nrow(new_k25)))

new_env <- new.env()
for (obj in ls(existing_env)) {
  assign(obj, get(obj, envir = existing_env), envir = new_env)
}
assign("empirical_q_hat_surface", empirical_q_hat_surface, envir = new_env)
save(list = ls(new_env), file = SYSDATA, envir = new_env, compress = "xz")
cat(sprintf("Wrote %s (%.1f MB)\n", SYSDATA,
            file.info(SYSDATA)$size / 1024 / 1024))
