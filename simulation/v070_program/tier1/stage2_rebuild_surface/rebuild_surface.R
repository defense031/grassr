#!/usr/bin/env Rscript
# Stage 2 — surface rebuild + verification (grassr v0.7.0 Tier 1).
#
# Merges the Stage-1 per-scenario outputs into the FULL-PRECISION surface
# archive (never the shipped 5-decimal sysdata), producing:
#   grassr/data-raw/sysdata_fullprecision_v0.7.0.rda   (new archive)
#   tier1/stage2_rebuild_surface/candidate_sysdata.rda (5-dec, shippable)
#   tier1/stage2_rebuild_surface/size_report.txt
# The package itself is NOT touched; Stage 5 does assembly.
#
# Verification: existing cells must be preserved exactly (append-only);
# the package test suite runs against the candidate in a scratch copy.
# Failures are CLASSIFIED, not auto-fatal: with a denser N/q grid,
# lookups that previously snapped (e.g. N=22 -> 50) now resolve to new
# cells, so tests pinning old-grid snap behavior may fail legitimately.
# DONE is written only on zero failures; otherwise NEEDS_REVIEW lists
# them for adjudication.
#
# Run from repo root after merging Sebastian's per_scenario into the
# laptop's (rsync); mirrors augment_empirical_bands_q094_k2_k25.R.

suppressMessages({
  source("paper2/code/04_reference_closed_form.R")
  source("paper2/code/06_grid.R")
})

ROOT      <- getwd()
STAGE1    <- "grassr/simulation/v070_program/tier1/stage1_surface_ndense/per_scenario"
STAGE2    <- "grassr/simulation/v070_program/tier1/stage2_rebuild_surface"
FULL_IN   <- "grassr/data-raw/sysdata_fullprecision_v0.6.0.rda"
FULL_OUT  <- "grassr/data-raw/sysdata_fullprecision_v0.7.0.rda"
CANDIDATE <- file.path(STAGE2, "candidate_sysdata.rda")
OFFSET    <- 300000L   # distinct from all prior scenario-id offsets

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

# ---- Collate stage-1 outputs ----------------------------------------------
files <- list.files(STAGE1, pattern = "^scenario_[0-9]+\\.rds$", full.names = TRUE)
cat(sprintf("Stage-1 per-scenario files: %d\n", length(files)))
stopifnot(length(files) == 34476L)

s_list <- vector("list", length(files))
per_rep <- vector("list", length(files))
for (i in seq_along(files)) {
  x <- readRDS(files[i])
  s_list[[i]] <- x$summary
  per_rep[[i]] <- x$per_rep
}
s <- do.call(rbind, s_list)
stopifnot(!any(duplicated(s$scenario_id)))
s$scenario_id <- s$scenario_id + OFFSET
names(per_rep) <- as.character(s$scenario_id)

# ---- Load full-precision base + anti-join assertion ------------------------
base_env <- new.env(); load(FULL_IN, envir = base_env)
eqs <- base_env$empirical_q_hat_surface
have <- paste(eqs$index$F_key, eqs$index$q_true, eqs$index$k, eqs$index$N)
new_keys <- paste(s$F_key, s$q, s$k, s$N)
stopifnot(!any(new_keys %in% have))     # append-only by construction
stopifnot(!any(duplicated(new_keys)))

# ---- Closed-form lookups + inversion (identical to augment mechanics) -----
uniq_Fk <- unique(s[, c("F_family", "F_key", "k")])
cat(sprintf("Building %d closed-form lookup tables...\n", nrow(uniq_Fk)))
lookups <- vector("list", nrow(uniq_Fk))
names(lookups) <- sprintf("%s|%s|%d", uniq_Fk$F_family, uniq_Fk$F_key, uniq_Fk$k)
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
    q_hat[mid] <- stats::approx(x = rc_ord, y = qg_ord, xout = obs_vec[mid],
                                rule = 2, ties = "ordered")$y
  }
  q_hat
}

n_scn <- nrow(s); n_met <- length(PER_REP_COL_FOR); n_q <- length(QUANTILE_PROBS)
cat(sprintf("Inverting %d scenarios x %d metrics...\n", n_scn, n_met))
q_hat_quantiles <- array(NA_real_, dim = c(n_scn, n_met, n_q),
  dimnames = list(scenario_id = as.character(s$scenario_id),
                  metric      = names(PER_REP_COL_FOR),
                  prob        = as.character(QUANTILE_PROBS)))
q_hat_mean <- matrix(NA_real_, nrow = n_scn, ncol = n_met,
                     dimnames = list(NULL, names(PER_REP_COL_FOR)))
q_hat_sd <- q_hat_mean

for (i in seq_len(n_scn)) {
  row <- s[i, ]
  M <- lookups[[sprintf("%s|%s|%d", row$F_family, row$F_key, row$k)]]
  reps <- per_rep[[as.character(row$scenario_id)]]
  if (is.null(reps)) next
  for (mi in seq_along(PER_REP_COL_FOR)) {
    obs_vec <- reps[, PER_REP_COL_FOR[[mi]]]
    q_hat_rep <- invert_vec(obs_vec, M[, names(PER_REP_COL_FOR)[mi]])
    q_hat_rep <- q_hat_rep[is.finite(q_hat_rep)]
    if (length(q_hat_rep) >= 10L) {
      q_hat_quantiles[i, mi, ] <- stats::quantile(
        q_hat_rep, probs = QUANTILE_PROBS, names = FALSE, type = 7)
      q_hat_mean[i, mi] <- mean(q_hat_rep)
      q_hat_sd[i, mi]   <- stats::sd(q_hat_rep)
    }
  }
  if (i %% 2000L == 0L) cat(sprintf("  %d/%d\n", i, n_scn))
}

compute_M1 <- function(F_family, F_key) {
  fp <- F_params_for(F_key)
  if (F_family == "logit_normal") e_p_logit_normal(mu = fp$mu, tau2 = fp$tau2)
  else sum(fp$w * fp$p)
}
M1_map <- unique(data.frame(
  F_key = uniq_Fk$F_key,
  M1 = vapply(seq_len(nrow(uniq_Fk)),
              function(i) compute_M1(uniq_Fk$F_family[i], uniq_Fk$F_key[i]),
              numeric(1L)),
  stringsAsFactors = FALSE))

new_index <- data.frame(scenario_id = s$scenario_id, F_family = s$F_family,
                        F_key = s$F_key, q_true = s$q, k = s$k, N = s$N,
                        stringsAsFactors = FALSE)
new_index <- merge(new_index, M1_map, by = "F_key", all.x = TRUE, sort = FALSE)
new_index <- new_index[order(new_index$scenario_id),
                       c("F_key", "scenario_id", "F_family", "q_true", "k", "N", "M1")]
rownames(new_index) <- NULL
ord <- order(as.integer(dimnames(q_hat_quantiles)$scenario_id))
q_hat_quantiles <- q_hat_quantiles[ord, , , drop = FALSE]
q_hat_mean <- q_hat_mean[ord, , drop = FALSE]
q_hat_sd   <- q_hat_sd[ord, , drop = FALSE]

# ---- Merge (append-only) ---------------------------------------------------
stopifnot(identical(names(eqs$index), names(new_index)))
merged_index <- rbind(eqs$index, new_index)
merged_quantiles <- array(NA_real_,
  dim = c(dim(eqs$quantiles)[1] + n_scn, dim(eqs$quantiles)[2], dim(eqs$quantiles)[3]),
  dimnames = list(
    scenario_id = c(dimnames(eqs$quantiles)$scenario_id,
                    dimnames(q_hat_quantiles)$scenario_id),
    metric = dimnames(eqs$quantiles)$metric,
    prob   = dimnames(eqs$quantiles)$prob))
merged_quantiles[seq_len(dim(eqs$quantiles)[1]), , ] <- eqs$quantiles
merged_quantiles[dim(eqs$quantiles)[1] + seq_len(n_scn), , ] <- q_hat_quantiles
merged_mean <- rbind(eqs$mean, q_hat_mean)
merged_sd   <- rbind(eqs$sd,   q_hat_sd)

# existing cells preserved exactly (bit-identical sub-arrays)
stopifnot(identical(merged_quantiles[seq_len(dim(eqs$quantiles)[1]), , ],
                    eqs$quantiles))

new_eqs <- eqs
new_eqs$index     <- merged_index
new_eqs$quantiles <- merged_quantiles
new_eqs$mean      <- merged_mean
new_eqs$sd        <- merged_sd
new_eqs$built_on  <- c(eqs$built_on,
                       sprintf("v0.7.0 stage1 N-densification %s", Sys.Date()))

out_env <- new.env()
for (nm in ls(base_env)) assign(nm, get(nm, base_env), out_env)
assign("empirical_q_hat_surface", new_eqs, out_env)
save(list = ls(out_env), envir = out_env, file = FULL_OUT, compress = "xz")
cat(sprintf("Full-precision archive: %s (%.2f MB, %d cells)\n",
            FULL_OUT, file.size(FULL_OUT)/1e6, nrow(merged_index)))

# ---- Candidate (5-decimal, matching the release rounding) ------------------
cs <- new_eqs
cs$quantiles <- round(cs$quantiles, 5)
cs$mean <- round(cs$mean, 5); cs$sd <- round(cs$sd, 5)
cand_env <- new.env()
for (nm in ls(out_env)) assign(nm, get(nm, out_env), cand_env)
assign("empirical_q_hat_surface", cs, cand_env)
f <- cand_env$fitted_icc_reference_curves
f$curves <- round(f$curves, 5)
assign("fitted_icc_reference_curves", f, cand_env)
d <- cand_env$delta_thresholds_lookup
attr(d, ".internal.selfref") <- NULL; class(d) <- "data.frame"
assign("delta_thresholds_lookup", d, cand_env)
save(list = ls(cand_env), envir = cand_env, file = CANDIDATE, compress = "xz")

rep_lines <- c(
  sprintf("candidate sysdata: %.2f MB (shipped 0.6.2: 2.83 MB)", file.size(CANDIDATE)/1e6),
  sprintf("surface cells: %d (was 10,140)", nrow(merged_index)),
  sprintf("N grid: %s", paste(sort(unique(merged_index$N)), collapse = ", ")),
  sprintf("q grid sizes by k: %s",
          paste(sapply(split(merged_index$q_true, merged_index$k),
                       function(x) length(unique(x))), collapse = ", ")))
writeLines(rep_lines, file.path(STAGE2, "size_report.txt"))
cat(paste(rep_lines, collapse = "\n"), "\n")
cat("Stage 2 build complete. Run verify_candidate.R next (writes DONE/NEEDS_REVIEW).\n")
