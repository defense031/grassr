#' @export
as.data.frame.grass_metrics <- function(x, row.names = NULL, optional = FALSE, ...) {
  v <- x$values
  df <- data.frame(
    N = x$n,
    positive_level = x$positive_level,
    t(as.matrix(v)),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = row.names
  )
  df
}

# --------------------------------------------------------------------------
# as.data.frame.grass_card -- v0.2.0 Target-2 Report Card tidy long-format
# --------------------------------------------------------------------------
# Returns x$panel augmented with an `is_primary` logical column. When the
# divergent flag fires and per_rater is populated, returns a list with two
# elements (panel, per_rater) so tidy consumers can pivot accordingly.
# Default returns the single panel data.frame for "aligned" / "caution"
# cards (which lack a per_rater table).

#' Coerce a grass_card to a tidy data.frame (panel; per_rater appended on
#' divergent)
#'
#' Returns the panel of observed coefficients with surface position metadata
#' as a data.frame. The primary coefficient is flagged via the `is_primary`
#' column. When the cross-coefficient flag is `divergent` and the per-rater
#' latent-class table is populated, returns a list `c(panel = ..., per_rater = ...)`
#' so tidy consumers can pivot the per-rater rows separately.
#'
#' @param x A `grass_card` object.
#' @param row.names,optional Standard `as.data.frame` arguments.
#' @param ... Unused.
#' @return A `data.frame` (the panel) when `x$per_rater` is `NULL` or empty.
#'   When the divergent flag fires, returns a list with `panel` and
#'   `per_rater` data frames.
#' @export
as.data.frame.grass_card <- function(x, row.names = NULL, optional = FALSE, ...) {
  panel <- x$panel
  panel$is_primary <- panel$coefficient == x$coefficient$primary
  if (!is.null(row.names)) rownames(panel) <- row.names

  if (!is.null(x$per_rater) && nrow(x$per_rater) > 0L) {
    return(list(panel = panel, per_rater = x$per_rater))
  }
  panel
}

#' @export
as.data.frame.grass_reference <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- x$reference
  if (!is.null(row.names)) rownames(df) <- row.names
  df
}

#' Coerce a grass_result to a tidy one-row data.frame
#'
#' @param x A `grass_result`.
#' @param row.names,optional Standard `as.data.frame` arguments.
#' @param compact If `TRUE`, drop `regime_note` from the output. The note is
#'   a long free-text string that wraps awkwardly when many rows are bound
#'   together. Default `FALSE` keeps it.
#' @param ... Unused.
#' @return A one-row data.frame with 20 columns (19 when `compact = TRUE`).
#'   Includes both `reference_level` (numeric Se = Sp band) and the
#'   legacy `reference_quality` character label; `reference_quality` is
#'   `NA` for the numeric-only bands (0.80, 0.90).
#' @export
as.data.frame.grass_result <- function(x, row.names = NULL, optional = FALSE,
                                       compact = FALSE, ...) {
  v <- x$metrics$values
  d <- x$distance
  ref_get <- function(metric, field = "reference") {
    if (is.null(d)) return(NA_real_)
    val <- d[[field]][d$metric == metric]
    if (length(val) == 0) NA_real_ else unname(val)
  }
  df <- data.frame(
    N                 = x$metrics$n,
    prevalence        = x$prevalence,
    prevalence_source = x$prevalence_source,
    regime            = x$regime,
    regime_note       = x$regime_note,
    kappa             = unname(v["kappa"]),
    kappa_ci_lower    = unname(v["kappa_wilson_lower"]),
    kappa_ci_upper    = unname(v["kappa_wilson_upper"]),
    kappa_ref         = ref_get("kappa"),
    kappa_distance    = ref_get("kappa", "distance"),
    PABAK             = unname(v["PABAK"]),
    PABAK_ref         = ref_get("PABAK"),
    PABAK_distance    = ref_get("PABAK", "distance"),
    AC1               = unname(v["AC1"]),
    AC1_ref           = ref_get("AC1"),
    AC1_distance      = ref_get("AC1", "distance"),
    prevalence_index  = unname(v["prevalence_index"]),
    bias_index        = unname(v["bias_index"]),
    reference_level   = if (!is.null(x$spec) && !is.null(x$spec$reference_level))
      x$spec$reference_level
    else if (!is.null(x$reference) && !is.null(x$reference$reference$reference_level))
      x$reference$reference$reference_level[1]
    else NA_real_,
    reference_quality = if (is.null(x$reference)) NA_character_ else x$reference$quality,
    stringsAsFactors = FALSE,
    row.names = row.names
  )
  if (isTRUE(compact)) df$regime_note <- NULL
  df
}
