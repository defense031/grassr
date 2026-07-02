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
