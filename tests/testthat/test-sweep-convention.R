# Tests for the v0.7.1 sweep / pooled-percentile / consistency-band
# convention (design/v0.7.1_position_redesign.md, ratified 2026-07-05).
# The monotonicity regression here is the test that would have caught the
# retired nearest-q_hat-cell sawtooth.

# ---------------------------------------------------------------------------
# (a) Monotonicity regression: the pooled percentile is non-decreasing in
#     the observed coefficient and spans a nontrivial range. Two designs --
#     one large (k=10, N=200), one tiny (k=3, N=20) -- because the sawtooth
#     was worst at large k / N but the reference must be informative at
#     small designs too.
# ---------------------------------------------------------------------------
test_that("pooled percentile is monotone non-decreasing in obs_value (sawtooth regression)", {
  obs_grid <- seq(0.30, 0.90, by = 0.02)
  for (design in list(c(k = 10, N = 200), c(k = 3, N = 20))) {
    pct <- vapply(obs_grid, function(o) {
      position_on_surface(obs_value = o, metric = "pabak",
                          pi_hat = 0.5,
                          k = design[["k"]], N = design[["N"]],
                          method = "empirical")$percentile
    }, numeric(1))
    expect_true(all(is.finite(pct)),
                info = sprintf("k=%d N=%d", design[["k"]], design[["N"]]))
    # Non-decreasing (allow tiny numerical slack).
    expect_true(all(diff(pct) >= -1e-9),
                info = sprintf("monotone failed at k=%d N=%d (min diff = %.3g)",
                               design[["k"]], design[["N"]], min(diff(pct))))
    # Nontrivial span -- the reference actually discriminates.
    expect_gt(max(pct) - min(pct), 0.2)
  }
})

# ---------------------------------------------------------------------------
# (b) Band sanity: a larger design (more raters, more subjects) gives a
#     tighter consistency band on quality; the band brackets q_hat when
#     both endpoints are finite and neither side is open; open flags are
#     logical.
# ---------------------------------------------------------------------------
test_that("consistency band is narrower at a larger design and brackets q_hat", {
  ov <- 0.60
  big   <- position_on_surface(ov, "pabak", pi_hat = 0.5, k = 8, N = 200,
                               method = "empirical")
  small <- position_on_surface(ov, "pabak", pi_hat = 0.5, k = 3, N = 20,
                               method = "empirical")

  expect_true(is.logical(big$band$open_low) && is.logical(big$band$open_high))
  expect_true(is.logical(small$band$open_low) && is.logical(small$band$open_high))

  big_width   <- big$band$hi - big$band$lo
  small_width <- small$band$hi - small$band$lo
  expect_true(is.finite(big_width) && is.finite(small_width))
  expect_lt(big_width, small_width)

  # Band brackets q_hat when two-sided and closed.
  if (is.finite(big$band$lo) && is.finite(big$band$hi) &&
      !isTRUE(big$band$open_low) && !isTRUE(big$band$open_high)) {
    expect_gte(big$q_hat, big$band$lo - 1e-6)
    expect_lte(big$q_hat, big$band$hi + 1e-6)
  }
})

# ---------------------------------------------------------------------------
# (c) Pooled percentile bounds + documented basis string.
# ---------------------------------------------------------------------------
test_that("pooled percentile is in [0, 1] with a documented basis", {
  r <- position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200,
                           method = "empirical")
  expect_gte(r$percentile, 0)
  expect_lte(r$percentile, 1)
  expect_true(r$percentile_basis %in%
              c("pooled-achievable-range",
                "pooled-achievable-range-delta-approx",
                "user-supplied-cohort"))
})

# ---------------------------------------------------------------------------
# (d) Deprecation: non-default bands warn "deprecated"; defaults are silent.
# ---------------------------------------------------------------------------
test_that("non-default bands warn deprecated; defaults are silent", {
  expect_warning(
    position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200,
                        bands = c(0.5, 0.6, 0.7, 0.8, 1.0)),
    "deprecated"
  )
  expect_no_warning(
    position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200)
  )
})

# ---------------------------------------------------------------------------
# (e) delta-B (Option B): delta_hat is the implied-quality spread in quality
#     pp. A symmetric panel (all raters same q) has a tiny spread (< 1 pp);
#     a strongly split-biased half-panel has a strictly larger spread. The
#     panel data.frame exposes the implied_q column.
# ---------------------------------------------------------------------------
test_that("check_asymmetry delta_hat separates symmetric from split-biased panels", {
  N <- 200L

  # Symmetric panel: every rater at Se = Sp = 0.90 against the same truth.
  set.seed(101L)
  truth_sym <- rbinom(N, 1, 0.4)
  sym <- matrix(0L, N, 5L)
  for (j in seq_len(5L)) {
    sym[, j] <- rbinom(N, 1, ifelse(truth_sym == 1L, 0.90, 1 - 0.90))
  }
  out_sym <- check_asymmetry(sym)
  expect_lt(out_sym$delta_hat, 1)                 # quality pp
  expect_true("implied_q" %in% names(out_sym$panel))

  # Split-biased half-panel: raters 1-2 symmetric (Se=Sp=0.95); raters 3-5
  # high-Se / low-Sp (Se=0.95, Sp=0.45). The cross-coefficient inversion
  # spreads the implied qualities.
  set.seed(202L)
  truth_asy <- rbinom(N, 1, 0.4)
  Se <- c(0.95, 0.95, 0.95, 0.95, 0.95)
  Sp <- c(0.95, 0.95, 0.45, 0.45, 0.45)
  asy <- matrix(0L, N, 5L)
  for (j in seq_len(5L)) {
    asy[, j] <- rbinom(N, 1, ifelse(truth_asy == 1L, Se[j], 1 - Sp[j]))
  }
  out_asy <- check_asymmetry(asy)
  expect_true("implied_q" %in% names(out_asy$panel))

  # The split-biased panel's implied-quality spread strictly exceeds the
  # symmetric panel's.
  expect_gt(out_asy$delta_hat, out_sym$delta_hat)
})

# ---------------------------------------------------------------------------
# (f) format_consistency_band renders two-sided, one-sided, and open cases.
# ---------------------------------------------------------------------------
test_that("format_consistency_band renders two-sided / one-sided / open bands", {
  fcb <- grassr:::format_consistency_band

  two_sided <- fcb(list(lo = 0.80, hi = 0.90, level = 0.95,
                        open_low = FALSE, open_high = FALSE))
  expect_match(two_sided, "0.80-0.90", fixed = TRUE)
  expect_match(two_sided, "quality")

  one_sided <- fcb(list(lo = NA_real_, hi = 0.90, level = 0.95,
                        open_low = FALSE, open_high = FALSE))
  expect_match(one_sided, "<= 0.90", fixed = TRUE)

  open_flagged <- fcb(list(lo = 0.80, hi = 0.90, level = 0.95,
                           open_low = TRUE, open_high = TRUE))
  expect_match(open_flagged, "<=", fixed = TRUE)  # open-low inequality marker
  expect_match(open_flagged, "+", fixed = TRUE)   # open-high marker

  # NULL band renders a graceful sentinel rather than erroring.
  expect_equal(fcb(NULL), "not derived")
})
