#!/usr/bin/env Rscript
# run_metric_recommendation.R — Which metric best recovers ground truth
# at each prevalence level and rater profile?
#
# Produces: output/figures/fig_metric_recommendation.png
#           output/metric_recommendation.rds

source("R/00_packages.R")
library(data.table)
library(ggplot2)

# ------------------------------------------------------------------
# 1. Load raw results and grid
# ------------------------------------------------------------------
cat("Loading raw results...\n")
results <- readRDS("output/sim_results/all_results.rds")
grid    <- readRDS("output/parameter_grid.rds")

# ------------------------------------------------------------------
# 2. Aggregate scenario-level means for kappa, PABAK, AC1
# ------------------------------------------------------------------
cat("Aggregating scenario-level means...\n")
scenario_means <- results[, .(
  kappa_mean = mean(kappa, na.rm = TRUE),
  PABAK_mean = mean(PABAK),
  AC1_mean   = mean(AC1, na.rm = TRUE)
), by = scenario_id]

# Merge with grid
scenario_means <- merge(scenario_means, grid, by = "scenario_id")

# ------------------------------------------------------------------
# 3. Define prevalence-conditioned reference values (grass v0.1.2)
# ------------------------------------------------------------------
# Analytical expected metric values on the Se = Sp = 0.85 diagonal under
# conditional independence. See grass::grass_reference_table(0.85) and
# Appendix C of the paper for the derivation.
thresholds <- data.table(
  prevalence    = c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                    0.48, 0.50, 0.52, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
                    0.90, 0.95, 0.99),
  kappa_thresh  = c(0.037, 0.154, 0.257, 0.329, 0.381, 0.419, 0.447, 0.466, 0.480, 0.487,
                    0.490, 0.490, 0.490, 0.487, 0.480, 0.466, 0.447, 0.419, 0.381, 0.329,
                    0.257, 0.154, 0.037),
  PABAK_thresh  = rep(0.490, 23),
  AC1_thresh    = c(0.653, 0.635, 0.612, 0.589, 0.566, 0.546, 0.527, 0.512, 0.500, 0.492,
                    0.490, 0.490, 0.490, 0.492, 0.500, 0.512, 0.527, 0.546, 0.566, 0.589,
                    0.612, 0.635, 0.653)
)

# Merge thresholds
scenario_means <- merge(scenario_means, thresholds, by = "prevalence")

# ------------------------------------------------------------------
# 4. Classify each scenario using each metric's threshold
# ------------------------------------------------------------------
# Ground truth: quality_band == "high" means min(Se, Sp) >= 0.90
scenario_means[, ground_truth := quality_band == "high"]

# Metric classification: does the metric exceed its prevalence-conditioned threshold?
scenario_means[, kappa_classifies_high := kappa_mean >= kappa_thresh]
scenario_means[, PABAK_classifies_high := PABAK_mean >= PABAK_thresh]
scenario_means[, AC1_classifies_high := AC1_mean >= AC1_thresh]

# Correct classification: metric agrees with ground truth
scenario_means[, kappa_correct := kappa_classifies_high == ground_truth]
scenario_means[, PABAK_correct := PABAK_classifies_high == ground_truth]
scenario_means[, AC1_correct := AC1_classifies_high == ground_truth]

# ------------------------------------------------------------------
# 5. Determine best metric per prevalence x profile
# ------------------------------------------------------------------
# Group profiles by type for cleaner visualization
profile_types <- data.table(
  profile_name = c(
    "symmetric_vlow", "symmetric_low", "symmetric_med", "symmetric_high", "symmetric_vhigh",
    "symmetric_high_alt", "high_se_dominant", "high_sp_dominant", "high_mixed", "high_between",
    "low_asym_se_high", "low_asym_sp_high", "low_between", "low_opposite", "low_moderate",
    "asym_within_se", "asym_within_sp", "asym_between_se", "asym_between_sp",
    "asym_between_both", "mixed_bias_1", "mixed_bias_2"
  ),
  profile_type = c(
    rep("Symmetric", 5),
    rep("High quality", 5),
    rep("Low quality", 5),
    rep("Asymmetric", 7)
  )
)

scenario_means <- merge(scenario_means, profile_types, by = "profile_name", all.x = TRUE)

# Accuracy by prevalence x profile, aggregated across sample sizes
accuracy_by_prev_profile <- scenario_means[, .(
  kappa_accuracy = mean(kappa_correct),
  PABAK_accuracy = mean(PABAK_correct),
  AC1_accuracy   = mean(AC1_correct),
  n_scenarios    = .N
), by = .(prevalence, profile_name, profile_type, quality_band)]

# Determine which metric is best for each prevalence x profile
accuracy_by_prev_profile[, best_metric := ifelse(
  kappa_accuracy >= PABAK_accuracy & kappa_accuracy >= AC1_accuracy, "kappa",
  ifelse(PABAK_accuracy >= AC1_accuracy, "PABAK", "AC1")
)]

# Handle ties: if all three are equal, mark as "All equal"
accuracy_by_prev_profile[, all_equal := (kappa_accuracy == PABAK_accuracy) &
                                         (PABAK_accuracy == AC1_accuracy)]
accuracy_by_prev_profile[all_equal == TRUE, best_metric := "All equal"]

# ------------------------------------------------------------------
# 6. Accuracy by prevalence x profile type (aggregated)
# ------------------------------------------------------------------
accuracy_by_prev_type <- scenario_means[, .(
  kappa_accuracy = mean(kappa_correct),
  PABAK_accuracy = mean(PABAK_correct),
  AC1_accuracy   = mean(AC1_correct),
  n_scenarios    = .N
), by = .(prevalence, profile_type)]

accuracy_by_prev_type[, best_metric := ifelse(
  kappa_accuracy >= PABAK_accuracy & kappa_accuracy >= AC1_accuracy, "kappa",
  ifelse(PABAK_accuracy >= AC1_accuracy, "PABAK", "AC1")
)]
accuracy_by_prev_type[, all_equal := (kappa_accuracy == PABAK_accuracy) &
                                      (PABAK_accuracy == AC1_accuracy)]
accuracy_by_prev_type[all_equal == TRUE, best_metric := "All equal"]

# ------------------------------------------------------------------
# 7. Heatmap: best metric by prevalence x profile
# ------------------------------------------------------------------

# Order profiles by quality band then name
accuracy_by_prev_profile[, profile_label := paste0(profile_name, " (", quality_band, ")")]
profile_order <- accuracy_by_prev_profile[, .(qb = quality_band[1]), by = profile_name]
profile_order <- profile_order[order(factor(qb, levels = c("high", "medium", "low")), profile_name)]
accuracy_by_prev_profile[, profile_label := factor(profile_label,
  levels = paste0(profile_order$profile_name, " (", profile_order$qb, ")"))]

p1 <- ggplot(accuracy_by_prev_profile, aes(x = factor(prevalence), y = profile_label, fill = best_metric)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("kappa" = "#E41A1C", "PABAK" = "#377EB8", "AC1" = "#4DAF4A", "All equal" = "#999999"),
    name = "Best metric"
  ) +
  labs(
    x = "Prevalence",
    y = "Rater profile (quality band)",
    title = "Which metric best recovers ground-truth rater quality?",
    subtitle = "Best = highest classification accuracy against ground truth, aggregated across sample sizes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  )

# ------------------------------------------------------------------
# 8. Heatmap: accuracy difference (best - worst) to show magnitude
# ------------------------------------------------------------------
accuracy_by_prev_profile[, max_accuracy := pmax(kappa_accuracy, PABAK_accuracy, AC1_accuracy)]
accuracy_by_prev_profile[, min_accuracy := pmin(kappa_accuracy, PABAK_accuracy, AC1_accuracy)]
accuracy_by_prev_profile[, accuracy_spread := max_accuracy - min_accuracy]

p2 <- ggplot(accuracy_by_prev_profile, aes(x = factor(prevalence), y = profile_label, fill = accuracy_spread)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "white", high = "#D7191C", name = "Accuracy\nspread") +
  labs(
    x = "Prevalence",
    y = "Rater profile (quality band)",
    title = "How much does metric choice matter?",
    subtitle = "Spread = classification accuracy of best metric minus worst metric"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  )

# ------------------------------------------------------------------
# 9. Summary table: accuracy by prevalence x metric (overall)
# ------------------------------------------------------------------
accuracy_overall <- scenario_means[, .(
  kappa_accuracy = mean(kappa_correct),
  PABAK_accuracy = mean(PABAK_correct),
  AC1_accuracy   = mean(AC1_correct)
), by = prevalence]
setorder(accuracy_overall, prevalence)

cat("\nClassification accuracy by prevalence (across all profiles and sample sizes):\n")
print(accuracy_overall, digits = 3)

# ------------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------------
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

saveRDS(list(
  by_profile    = accuracy_by_prev_profile,
  by_type       = accuracy_by_prev_type,
  overall       = accuracy_overall,
  scenario_data = scenario_means
), "output/metric_recommendation.rds")

ggsave("output/figures/fig_metric_recommendation_best.png", p1,
       width = 12, height = 10, dpi = 300)
ggsave("output/figures/fig_metric_recommendation_best.pdf", p1,
       width = 12, height = 10)

ggsave("output/figures/fig_metric_recommendation_spread.png", p2,
       width = 12, height = 10, dpi = 300)
ggsave("output/figures/fig_metric_recommendation_spread.pdf", p2,
       width = 12, height = 10)

cat("\nFigures saved to output/figures/fig_metric_recommendation.{png,pdf}\n")
cat("Data saved to output/metric_recommendation.rds\n")
