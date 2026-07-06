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
  # k = 2 is NOT calibrated and must never snap to k = 3: at k = 2 the
  # agreement family is PABAK + AC1, whose implied qualities coincide by
  # construction (delta_hat is identically zero — 2.75M Option-B null
  # draws produced no exception). delta_hat carries no information at
  # k = 2; check_asymmetry() reports not_applicable there.
  if (as.numeric(k) < 3) return(NULL)
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
# [~0.5, 99.5+] percent. Values beyond the stored 99.5th report as > 99.5.
#
# Tied quantile runs (point masses) use the MID-P convention:
#   percentile = 100 * (P(D < d) + 0.5 * P(D = d))
# the standard treatment for discrete nulls. Small-N Option-B nulls carry
# real point mass at delta_hat = 0 (41% of the grid at k=5/N=1000/q=0.97);
# reporting the BOTTOM of the tie run understates the position (the 0.7.0
# v[1] short-circuit bug: "1.0 percentile" where P(D <= 0) was ~35%) and
# reporting the TOP misfires the flag convention on heavily-plateaued
# cells (an observed 0 on a mostly-zero null would read as extreme).
delta_null_percentile <- function(delta_hat, cell) {
  if (!is.finite(delta_hat) || is.null(cell)) return(NA_real_)
  v <- cell$values; p <- cell$probs
  tol <- 1e-12
  if (delta_hat <= v[1L] + tol) {
    # Observation at (or below) the grid floor. Tie run = every stored
    # quantile equal to the floor value; P(D < d) ~ 0 below the stored
    # grid, P(D <= d) ~ prob at the top of the run.
    run_top <- max(which(v <= v[1L] + tol))
    return(100 * 0.5 * p[run_top])
  }
  if (delta_hat >= v[length(v)] - tol) {
    # At or beyond the stored ceiling. If the ceiling itself is a tie
    # run (mass at the max), mid-p between the run bottom and the cap;
    # strictly beyond the ceiling reports the cap.
    if (delta_hat > v[length(v)] + tol) return(100 * p[length(p)])
    run_bot <- min(which(v >= v[length(v)] - tol))
    lo <- if (run_bot == 1L) 0 else p[run_bot - 1L]
    return(100 * 0.5 * (lo + p[length(p)]))
  }
  # Interior: if d lands inside a tie run, mid-p across the run;
  # otherwise plain interpolation on the strictly-increasing segments.
  ties_at <- which(abs(v - delta_hat) <= tol)
  if (length(ties_at) > 1L) {
    lo <- if (min(ties_at) == 1L) 0 else p[min(ties_at) - 1L]
    return(100 * 0.5 * (lo + p[max(ties_at)]))
  }
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
