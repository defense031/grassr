# Tests for the v0.7.0 matched-null percentile convention.

test_that("lookup_delta_null resolves exact and snapped cells", {
  cell <- grassr:::lookup_delta_null(k = 5, N = 200, q_hat = 0.85)
  expect_equal(c(cell$k, cell$N), c(5, 200))
  expect_equal(cell$q, 0.85)
  expect_false(cell$snapped)
  expect_gte(cell$n_draws, 49000)
  snap <- grassr:::lookup_delta_null(k = 9, N = 220, q_hat = 0.80)
  expect_true(snap$snapped)
  expect_true(snap$k %in% c(8, 10))
})

test_that("delta_null_percentile is monotone and bounded", {
  # Probe values are taken RELATIVE to the cell's own quantile scale
  # (Option-B nulls are in quality pp, ~1000x smaller than the retired
  # percentile-pp units, so absolute probe constants would saturate).
  cell <- grassr:::lookup_delta_null(k = 5, N = 200, q_hat = 0.85)
  v_lo  <- cell$values[["0.05"]]
  v_mid <- cell$values[["0.5"]]
  v_hi  <- 2 * max(cell$values) + 1
  p_lo  <- grassr:::delta_null_percentile(v_lo, cell)
  p_mid <- grassr:::delta_null_percentile(v_mid, cell)
  p_hi  <- grassr:::delta_null_percentile(v_hi, cell)
  expect_lt(p_lo, p_mid); expect_lt(p_mid, p_hi)
  # mid-p at the median can drift a few points on tie runs (quantized
  # tight ridges); it must still sit near the middle.
  expect_gt(p_mid, 40); expect_lt(p_mid, 60)
  expect_lte(p_hi, 99.5 + 1e-9)
})

test_that("delta_null_percentile applies the MID-P convention at plateaus", {
  # Structural (unit-agnostic) test of the tie handling: point masses use
  # percentile = 100 * (P(D < d) + 0.5 * P(D = d)). Here v[1..3] all
  # equal 0, so P(D < 0) ~ 0 and P(D = 0) resolves to probs[3] = 0.10 ->
  # mid-p 5 (top-of-run 10 would misfire the flag convention on heavily
  # plateaued cells; bottom-of-run 1 understates the position).
  cell <- list(values = c(0, 0, 0, 0.5, 1),
               probs  = c(0.01, 0.05, 0.10, 0.5, 0.99))
  expect_equal(grassr:::delta_null_percentile(0, cell), 5)
  # A value at/below the grid floor resolves through the same tie run.
  expect_equal(grassr:::delta_null_percentile(-1, cell), 5)
  # A fully degenerate (all-zero) grid puts an observed zero mid-mass,
  # never in the extreme tail.
  cell0 <- list(values = rep(0, 5), probs = c(0.01, 0.05, 0.10, 0.5, 0.995))
  expect_lt(grassr:::delta_null_percentile(0, cell0), 60)
})

test_that("flag conventions map percentiles correctly", {
  f <- grassr:::delta_flag_from_percentile
  expect_equal(f(50), "aligned")
  expect_equal(f(96), "caution")
  expect_equal(f(99.2), "divergent")
  expect_equal(f(NA_real_), "not_calibrated")
})

test_that("check_asymmetry emits matched-null fields", {
  set.seed(11)
  Y <- matrix(rbinom(500 * 5, 1, 0.5), 500, 5)
  a <- check_asymmetry(Y)
  expect_true(is.finite(a$delta_percentile))
  expect_equal(a$thresholds_source, "matched_null_ecdf")
  expect_named(a$thresholds, c("caution", "divergent"))
})
