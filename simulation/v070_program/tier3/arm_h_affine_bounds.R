#!/usr/bin/env Rscript
# Tier 3 Arm H (analytic form, ratified 2026-07-05) — within-subject
# dependence bounds for the intra-rater axis, from shipped sysdata.
#
# At W = 2 occasions, ANY within-subject dependence (memory or a rater's
# own consistent rubric — unidentifiable from a single rater's matrix)
# maps to an effective carryover rho: with probability rho the second
# viewing repeats the first, so
#   match_rho = rho + (1 - rho) * match_0   =>
#   PABAK_rho = rho + (1 - rho) * PABAK_0     (PABAK = 2*match - 1)
# The bound on the flattering direction is therefore an affine transform
# of the calibrated k = 2 surface — no simulation needed (the 0.7.0 plan's
# 540k-panel arm collapses to this table; TIER3_DESIGN.md Arm H).
#
# For each (N, q) cell at k = 2 (balanced F key) and rho in {0.1, .25, .5}:
#   - PABAK_0   = cell median coefficient (dependence-free panel)
#   - PABAK_rho = affine transform
#   - inflation of the POOLED percentile (v0.7.1 convention)
#   - inflation of the implied quality q_hat and the consistency band
# Output: arm_h_affine_bounds.rds + plain-text table alongside.

CAND <- Sys.getenv("CANDIDATE_LIB",
  file.path(dirname(dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))),
    "candidate_lib_v071"))
.libPaths(c(CAND, .libPaths()))
suppressMessages(library(grassr, lib.loc = CAND))
`%||%` <- function(x, y) if (is.null(x)) y else x

OUT_DIR <- dirname(normalizePath(sub("^--file=",
  "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))

surf <- get("empirical_q_hat_surface", envir = asNamespace("grassr"))
idx  <- surf$index
mi   <- which(surf$metrics == "pabak")
FK   <- "LN_mu=+0.000_tau2=1.0000"          # balanced prevalence reference

rows <- which(idx$F_key == FK & idx$k == 2L)
cells <- idx[rows, c("q_true", "N", "scenario_id")]
arr_ids <- as.integer(dimnames(surf$quantiles)$scenario_id)

RHOS <- c(0.10, 0.25, 0.50)
res <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  sc  <- cells[i, ]
  r   <- match(as.integer(sc$scenario_id), arr_ids)
  qh_med <- surf$quantiles[r, mi, which.min(abs(surf$probs - 0.5))]
  if (!is.finite(qh_med)) return(NULL)
  # median q_hat -> median PABAK of the dependence-free cell
  pab0 <- (2 * qh_med - 1)^2
  base <- position_on_surface(obs_value = pab0, metric = "pabak",
                              pi_hat = 0.5, k = 2L, N = sc$N)
  do.call(rbind, lapply(RHOS, function(rho) {
    pabr <- rho + (1 - rho) * pab0
    pos  <- position_on_surface(obs_value = pabr, metric = "pabak",
                                pi_hat = 0.5, k = 2L, N = sc$N)
    data.frame(N = sc$N, q_true = sc$q_true, rho = rho,
               pabak_0 = pab0, pabak_rho = pabr,
               pct_0 = 100 * base$percentile,
               pct_rho = 100 * pos$percentile,
               pct_inflation = 100 * (pos$percentile - base$percentile),
               q_hat_0 = base$q_hat, q_hat_rho = pos$q_hat,
               q_inflation_pp = 100 * (pos$q_hat - base$q_hat),
               band_lo_rho = pos$band$lo %||% NA_real_,
               band_hi_rho = pos$band$hi %||% NA_real_)
  }))
}))

saveRDS(res, file.path(OUT_DIR, "arm_h_affine_bounds.rds"))

# summary: worst-case + typical inflation per rho
s <- do.call(rbind, lapply(split(res, res$rho), function(d) {
  data.frame(rho = d$rho[1],
             med_pct_inflation = median(d$pct_inflation, na.rm = TRUE),
             p95_pct_inflation = quantile(d$pct_inflation, 0.95, na.rm = TRUE),
             med_q_inflation_pp = median(d$q_inflation_pp, na.rm = TRUE),
             p95_q_inflation_pp = quantile(d$q_inflation_pp, 0.95, na.rm = TRUE))
}))
txt <- file.path(OUT_DIR, "arm_h_affine_bounds_summary.txt")
sink(txt)
cat("Arm H analytic dependence bounds (W = 2, k = 2 surface, balanced F)\n")
cat("PABAK_rho = rho + (1 - rho) * PABAK_0; pooled-percentile convention\n\n")
print(s, row.names = FALSE, digits = 3)
cat("\nFull grid:", nrow(res), "rows (N x q x rho) in arm_h_affine_bounds.rds\n")
sink()
cat(readLines(txt), sep = "\n")
