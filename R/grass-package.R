#' grass: Prevalence-Aware Interpretation of Binary Inter-Rater Agreement
#'
#' Computes Cohen's kappa, PABAK, and Gwet's AC1 for binary inter-rater
#' agreement, and interprets them against simulation-derived,
#' prevalence-conditioned thresholds from the GRASS framework.
#'
#' @section Entry points:
#' * [grass_report()] -- contextual profile: metrics + PI + BI + regime + reference
#' * [grass_spec_binary()] -- spec for the binary family (the one implemented today)
#' * [grass_compute()] -- raw metric panel from 2x2 or rater data
#' * [grass_reference()] -- look up reference curves at a given prevalence
#' * [grass_reference_table()] -- full internal reference-curve table
#' * [grass_prevalence()] -- estimate prevalence from rater marginals
#' * [grass_plot()] -- plot a grass object
#'
#' @section Roadmap:
#' See [grass_roadmap] for the framework taxonomy and the planned
#' ordinal / multirater / continuous families.
#'
#' @keywords internal
#' @name grass-package
#' @aliases grass
#' @importFrom stats approx plogis qnorm reshape setNames
"_PACKAGE"

# Silence R CMD check NOTEs for ggplot2 aesthetic columns.
utils::globalVariables(c(
  "prevalence", "threshold", "metric", "value", "quality",
  "lower", "upper", "verdict", "J", "x", "y", "label",
  "PI2", "BI2", "reference_binary"
))
