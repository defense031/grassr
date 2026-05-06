#!/usr/bin/env Rscript
# run_recovery_analysis.R — Post-hoc rater-quality recovery analysis.
#
# For each simulation scenario, invert the observed kappa_mean, PABAK_mean,
# and AC1_mean back to an estimated operating quality q under the Se = Sp
# diagonal analytical model (grass v0.1.2). Compare each q-estimate to the
# generating ground truth and report recovery error by prevalence x N.
#
# Also computes the per-scenario q-spread diagnostic:
#   q_spread = max(q_kappa, q_PABAK, q_AC1) - min(q_kappa, q_PABAK, q_AC1)
# which practitioners can compute on their own data as a single-number
# summary of metric discordance.
#
# Restricted to symmetric (Se = Sp) rater profiles where the generating q
# is unambiguous.
#
# Produces:
#   output/recovery_analysis.rds
#   output/figures/fig_recovery_heatmap.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------
# 1. Load unified sim output
# ------------------------------------------------------------------
uni <- readRDS("output/unified_sim/unified_results.rds")
sm  <- as.data.table(uni$scenario_means)

# Ground-truth q: each rater's Se and Sp. Restrict to scenarios where
# Se1 == Sp1 == Se2 == Sp2 so the generating q is unambiguous.
sym <- sm[Se1 == Sp1 & Se2 == Sp2 & Se1 == Se2]
sym[, q_true := Se1]

cat("Symmetric-scenario subset: ", nrow(sym), " of ", nrow(sm), " scenarios\n", sep = "")
cat("Unique q_true values in symmetric subset: ",
    paste(sort(unique(sym$q_true)), collapse = ", "), "\n")

# ------------------------------------------------------------------
# 2. Inversion helpers (analytical model: Se = Sp = q, cond. indep.)
# ------------------------------------------------------------------
# PABAK = 2 * (q^2 + (1-q)^2) - 1  => q = (1 + sqrt(PABAK)) / 2
# (competent-rater branch; PABAK must be >= 0)
invert_PABAK <- function(pabak) {
  ifelse(pabak < 0, NA_real_, (1 + sqrt(pmax(pabak, 0))) / 2)
}

# kappa and AC1 depend on prevalence too; closed forms exist but are messier.
# We solve numerically on [0.5, 1]. Under symmetry P0 = q^2 + (1-q)^2 does
# not depend on prevalence; so kappa_target(q) - kappa_obs has a unique root
# in q > 0.5 for each (p, kappa_obs).
expected_kappa <- function(q, p) {
  P0  <- q^2 + (1-q)^2
  pi_ <- p*q + (1-p)*(1-q)
  Pe  <- pi_^2 + (1-pi_)^2
  (P0 - Pe) / (1 - Pe)
}
expected_AC1 <- function(q, p) {
  P0  <- q^2 + (1-q)^2
  pi_ <- p*q + (1-p)*(1-q)
  Pe  <- 2*pi_*(1-pi_)
  (P0 - Pe) / (1 - Pe)
}
invert_on_diagonal <- function(obs, p, efun) {
  if (is.na(obs) || obs <= -1 || obs >= 1) return(NA_real_)
  f <- function(q) efun(q, p) - obs
  # bracket: at q = 0.5 metric = 0; at q -> 1 metric -> 1
  if (f(0.5) >= 0 && obs <= 0) return(0.5)
  tryCatch(uniroot(f, c(0.5, 0.9999), tol = 1e-6)$root,
           error = function(e) NA_real_)
}
invert_kappa <- function(kappa, p) invert_on_diagonal(kappa, p, expected_kappa)
invert_AC1   <- function(ac1, p)   invert_on_diagonal(ac1,   p, expected_AC1)

# ------------------------------------------------------------------
# 3. Invert each scenario
# ------------------------------------------------------------------
cat("Inverting metrics for ", nrow(sym), " scenarios...\n", sep = "")
sym[, q_PABAK := invert_PABAK(PABAK_mean)]
sym[, q_kappa := mapply(invert_kappa, kappa_mean, prevalence)]
sym[, q_AC1   := mapply(invert_AC1,   AC1_mean,   prevalence)]

# q-spread diagnostic: max - min across the three q-estimates
sym[, q_spread := pmax(q_PABAK, q_kappa, q_AC1, na.rm = FALSE) -
                  pmin(q_PABAK, q_kappa, q_AC1, na.rm = FALSE)]

# Recovery error per metric: q_hat - q_true
sym[, err_kappa := q_kappa - q_true]
sym[, err_PABAK := q_PABAK - q_true]
sym[, err_AC1   := q_AC1   - q_true]

# ------------------------------------------------------------------
# 4. Aggregate by (prevalence, N)
# ------------------------------------------------------------------
agg_by_prev_N <- sym[, .(
  n_scenarios    = .N,
  rmse_kappa     = sqrt(mean(err_kappa^2, na.rm = TRUE)),
  rmse_PABAK     = sqrt(mean(err_PABAK^2, na.rm = TRUE)),
  rmse_AC1       = sqrt(mean(err_AC1^2,   na.rm = TRUE)),
  mean_q_spread  = mean(q_spread, na.rm = TRUE),
  median_q_spread = median(q_spread, na.rm = TRUE)
), by = .(prevalence, N)]
setorder(agg_by_prev_N, prevalence, N)

# Also summarize overall by prevalence alone (across all N)
agg_by_prev <- sym[, .(
  n_scenarios    = .N,
  rmse_kappa     = sqrt(mean(err_kappa^2, na.rm = TRUE)),
  rmse_PABAK     = sqrt(mean(err_PABAK^2, na.rm = TRUE)),
  rmse_AC1       = sqrt(mean(err_AC1^2,   na.rm = TRUE)),
  mean_q_spread  = mean(q_spread, na.rm = TRUE)
), by = prevalence]
setorder(agg_by_prev, prevalence)

cat("\n=== Recovery RMSE by prevalence (symmetric subset, all N) ===\n")
print(agg_by_prev, digits = 3)

# ------------------------------------------------------------------
# 5. Figure: recovery RMSE heatmap (prevalence x N, faceted by metric)
# ------------------------------------------------------------------
long <- melt(agg_by_prev_N,
             id.vars = c("prevalence", "N"),
             measure.vars = c("rmse_kappa", "rmse_PABAK", "rmse_AC1"),
             variable.name = "metric",
             value.name = "rmse")
long[, metric := factor(metric,
                        levels = c("rmse_kappa", "rmse_PABAK", "rmse_AC1"),
                        labels = c("kappa", "PABAK", "AC1"))]

p <- ggplot(long, aes(x = factor(prevalence), y = factor(N), fill = rmse)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rmse)), size = 2.5) +
  scale_fill_gradient2(low = "#1A9641", mid = "#FFFFBF", high = "#D7191C",
                       midpoint = 0.10, limits = c(0, NA),
                       name = "RMSE of\nq-recovery") +
  facet_wrap(~ metric, nrow = 1) +
  labs(
    x = "Prevalence",
    y = "Sample size (N)",
    title = "Rater-quality recovery error by prevalence and N",
    subtitle = expression(paste(
      "RMSE of ", hat(q), " against the generating q across symmetric (Se = Sp) scenarios"))
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figures/fig_recovery_heatmap.png", p,
       width = 13, height = 4.5, dpi = 300)

# ------------------------------------------------------------------
# 6. Persist
# ------------------------------------------------------------------
saveRDS(list(
  scenarios      = sym,
  by_prev_N      = agg_by_prev_N,
  by_prev        = agg_by_prev
), "output/recovery_analysis.rds")

cat("\nSaved output/recovery_analysis.rds and output/figures/fig_recovery_heatmap.png\n")
