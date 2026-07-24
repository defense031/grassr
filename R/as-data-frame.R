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
#' Returns the panel of observed coefficients with surface-position metadata
#' as a data.frame. Each row carries the v0.7.1 panel columns
#' (`surface_percentile`, `band_lo` / `band_hi` / `band_open_low` /
#' `band_open_high`, `q_hat`, `se_q_hat`, `clamped`, `reference_used`,
#' `in_delta_hat`), an `is_primary` flag, and the panel-level delta
#' diagnostics recycled across rows (`delta_hat` implied-quality spread in
#' pp, `delta_percentile`, `delta_flag`, and the matched-null cell
#' `matched_null_k` / `matched_null_N` / `matched_null_q`). When the
#' cross-coefficient flag is `divergent` and the per-rater latent-class
#' table is populated, returns a list `c(panel = ..., per_rater = ...)` so
#' tidy consumers can pivot the per-rater rows separately.
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

  # Attach the panel-level delta diagnostics, recycled across coefficient
  # rows, so a single flat frame carries the flag and its basis alongside
  # each coefficient (v0.7.1: delta_hat is the implied-quality spread in pp;
  # the flag comes from delta_percentile on the matched-null cell).
  d  <- x$delta
  mn <- d$matched_null
  panel$delta_hat        <- as.numeric(d$delta_hat %||% NA_real_)
  panel$delta_percentile <- as.numeric(d$delta_percentile %||% NA_real_)
  panel$delta_flag       <- as.character(d$flag %||% NA_character_)
  panel$matched_null_k   <- if (!is.null(mn)) as.integer(mn$k) else NA_integer_
  panel$matched_null_N   <- if (!is.null(mn)) as.integer(mn$N) else NA_integer_
  panel$matched_null_q   <- if (!is.null(mn)) as.numeric(mn$q) else NA_real_

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
