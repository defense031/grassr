# delta_null.R — matched-null lookup and percentile positioning for the
# cross-coefficient diagnostic.
#
# The per-(k, N) threshold table is retired (0.7.0). delta_hat is reported
# as its percentile on the null distribution of delta_hat at the matched
# (k, N, q_hat) cell, calibrated by the production-pipeline null program.
# Flags are conventions on that percentile: caution at the 95th, divergent
# at the 99th. The percentile is computed by interpolation on a stored
# fine quantile grid (1% steps + 99.5%).
#
# v0.7.1 (Option B): delta_hat is the implied-quality spread in quality
# percentage points (design/v0.7.1_position_redesign.md). The stored null
# MUST be calibrated in the same units — the stage 6 regeneration under
# the Option-B pipeline replaces the 0.7.0 percentile-spread null.

# Resolve the nearest calibrated (k, N, q) cell. k and N snap by the same
# normalized (k, log10 N) metric the surface lookup uses; q snaps to the
# nearest calibrated quality level. Snap distances are returned so the
# card can disclose them.
lookup_delta_null <- function(k, N, q_hat) {
  obj <- tryCatch(
    get("delta_null_ecdf", envir = asNamespace("grassr"), inherits = FALSE),
    error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  idx <- obj$index
  qs <- sort(unique(idx$q))
  q_near <- qs[which.min(abs(qs - q_hat))]
  sub <- idx[idx$q == q_near, , drop = FALSE]
  dk <- (sub$k - k) / diff(range(idx$k))
  dN <- (log10(sub$N) - log10(N)) / diff(range(log10(idx$N)))
  i <- which.min(sqrt(dk^2 + dN^2))
  row <- sub[i, ]
  cell_i <- which(idx$k == row$k & idx$N == row$N & idx$q == row$q)
  list(
    k = row$k, N = row$N, q = row$q,
    n_draws = row$n_draws,
    unstable_tail = isTRUE(row$unstable_tail),
    snapped = (row$k != k || row$N != N || abs(row$q - q_hat) > 1e-9),
    probs = obj$probs,
    values = obj$values[cell_i, ],
    conventions = obj$flag_conventions
  )
}

# Percentile of an observed delta_hat on the matched null ECDF, in
# [1, 99.5+] percent. Values beyond the stored 99.5th report as > 99.5.
#
# Tied quantile runs (plateaus) are resolved to the TOP of the run: at
# small N the null often has a point mass at delta_hat = 0 that can reach
# well past the stored 1st percentile, and reporting a value inside the
# plateau at the run's lowest prob would understate P(D <= delta_hat)
# (0.7.0 shipped exactly that bug at the v[1] short-circuit: identical
# raters printed "1.0 percentile" on cells where P(D <= 0) was ~35%).
delta_null_percentile <- function(delta_hat, cell) {
  if (!is.finite(delta_hat) || is.null(cell)) return(NA_real_)
  v <- cell$values; p <- cell$probs
  if (delta_hat <= v[1L]) {
    # Highest stored prob whose quantile still equals (or exceeds by
    # numerical tie) the lowest stored value: the top of the zero/floor
    # plateau, i.e. P(D <= v[1]) as the grid resolves it.
    run <- which(v <= v[1L] + 1e-12)
    return(100 * p[max(run)])
  }
  if (delta_hat >= v[length(v)]) return(100 * p[length(p)])
  100 * stats::approx(x = v, y = p, xout = delta_hat,
                      rule = 2, ties = "ordered")$y
}

# Flag from the percentile conventions. A cell with an unstable tail
# still flags (the ECDF is stable); the instability is surfaced as a
# note, not a verdict change.
delta_flag_from_percentile <- function(pct, conventions = c(caution = 0.95,
                                                            divergent = 0.99)) {
  if (!is.finite(pct)) return("not_calibrated")
  if (pct >= 100 * conventions[["divergent"]]) return("divergent")
  if (pct >= 100 * conventions[["caution"]]) return("caution")
  "aligned"
}
