# broom::tidy methods -- registered conditionally in zzz.R so `broom` is not
# a hard dependency.

tidy.grass_metrics <- function(x, ...) {
  v <- x$values
  data.frame(
    term     = names(v),
    estimate = unname(as.numeric(v)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

tidy.grass_reference <- function(x, ...) {
  x$reference
}

# tidy.grass_card -- one row per panel coefficient in the v0.7.1 shape.
# Exposes the pooled surface percentile, the consistency-band endpoints
# (band_lo / band_hi / band_open_low / band_open_high), implied quality
# q_hat / se_q_hat, and the panel-level delta diagnostics (implied-quality
# spread delta_hat in pp, delta_percentile, flag, and the matched-null
# cell) recycled across rows. The retired modal-band label / confidence
# qualifier columns are gone. Unlike as.data.frame.grass_card(), tidy()
# always returns a single flat data.frame (per-rater rows are not appended)
# so downstream broom consumers get a stable rectangular shape.
tidy.grass_card <- function(x, ...) {
  panel <- x$panel
  panel$is_primary <- panel$coefficient == x$coefficient$primary

  d  <- x$delta
  mn <- d$matched_null
  panel$delta_hat        <- as.numeric(d$delta_hat %||% NA_real_)
  panel$delta_percentile <- as.numeric(d$delta_percentile %||% NA_real_)
  panel$delta_flag       <- as.character(d$flag %||% NA_character_)
  panel$matched_null_k   <- if (!is.null(mn)) as.integer(mn$k) else NA_integer_
  panel$matched_null_N   <- if (!is.null(mn)) as.integer(mn$N) else NA_integer_
  panel$matched_null_q   <- if (!is.null(mn)) as.numeric(mn$q) else NA_real_

  rownames(panel) <- NULL
  panel
}
