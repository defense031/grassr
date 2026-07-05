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
  cell <- grassr:::lookup_delta_null(k = 5, N = 200, q_hat = 0.85)
  p_lo <- grassr:::delta_null_percentile(0.5, cell)
  p_mid <- grassr:::delta_null_percentile(cell$values[["0.5"]], cell)
  p_hi <- grassr:::delta_null_percentile(60, cell)
  expect_lt(p_lo, p_mid); expect_lt(p_mid, p_hi)
  expect_equal(round(p_mid), 50)
  expect_lte(p_hi, 99.5 + 1e-9)
})

test_that("delta_null_percentile resolves a zero-plateau to the top of the tied run", {
  # Structural (unit-agnostic) test of the zero-plateau lower-tail fix:
  # at delta_hat <= v[1] the percentile is the prob at the TOP of the
  # tied run at v[1], not blindly p[1]. Here v[1..3] all equal 0, so
  # P(D <= 0) resolves to probs[3] = 0.10 -> 10 pct, not probs[1] = 1 pct.
  cell <- list(values = c(0, 0, 0, 0.5, 1),
               probs  = c(0.01, 0.05, 0.10, 0.5, 0.99))
  expect_equal(grassr:::delta_null_percentile(0, cell), 10)
  # A value strictly inside the plateau (still <= v[1]) resolves the same.
  expect_equal(grassr:::delta_null_percentile(-1, cell), 10)
})

test_that("flag conventions map percentiles correctly", {
  f <- grassr:::delta_flag_from_percentile
  expect_equal(f(50), "aligned")
  expect_equal(f(96), "caution")
  expect_equal(f(99.2), "divergent")
  expect_equal(f(NA_real_), "not_calibrated")
})

test_that("check_asymmetry emits matched-null fields and honors legacy arg", {
  set.seed(11)
  Y <- matrix(rbinom(500 * 5, 1, 0.5), 500, 5)
  a <- check_asymmetry(Y)
  expect_true(is.finite(a$delta_percentile))
  expect_equal(a$thresholds_source, "matched_null_ecdf")
  expect_named(a$thresholds, c("caution", "divergent"))
  expect_warning(b <- check_asymmetry(Y, delta_thresholds = c(9.25, 11.75)),
                 "deprecated")
  expect_equal(b$thresholds_source, "user_supplied_legacy")
})
