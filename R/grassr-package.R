#' grassr: Context-Conditioned Reporting for Binary Rater Reliability
#'
#' Generates a Report Card for rater reliability on binary outcomes from
#' an N x k subject-by-rater rating matrix. Each panel coefficient is
#' positioned on a simulation-calibrated reference surface conditioned on
#' the study's rater count, sample size, and prevalence; a
#' cross-coefficient discordance diagnostic (delta-hat) flags panels for
#' which no single coefficient is a stable summary. The package implements
#' the GRASS framework (Guide for Rater Agreement under Structural Skew).
#'
#' @section Entry points:
#' * [grass_report()] -- the Report Card: coefficient panel, surface
#'   percentiles, delta-hat flag, and divergent-path routing
#' * [position_on_surface()] -- position one observed coefficient on its
#'   calibrated reference surface
#' * [pairwise_agreement()] -- pairwise PABAK matrix on the k = 2 surface
#' * [latent_class_fit()] -- per-rater sensitivity/specificity via
#'   Dawid-Skene (Hui-Walter bounds at k = 2)
#' * [grass_compute()] -- raw metric panel from 2x2 or rater data
#' * [plot_surface()] -- plot an observed coefficient against its surface
#'
#' @section Roadmap:
#' See [grass_roadmap] for the framework taxonomy and the planned
#' ordinal / multirater / continuous families.
#'
#' @keywords internal
#' @name grassr-package
#' @aliases grassr
#' @importFrom stats approx plogis qnorm reshape setNames
"_PACKAGE"

# Silence R CMD check NOTEs for ggplot2 aesthetic columns.
utils::globalVariables(c(
  "prevalence", "threshold", "metric", "value", "quality",
  "lower", "upper", "verdict", "J", "x", "y", "label",
  "PI2", "BI2", "reference_binary"
))
