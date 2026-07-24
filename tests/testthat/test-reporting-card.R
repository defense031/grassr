# NOTE (v0.2.0): `check_asymmetry(se = ..., sp = ...)` was renamed to
# `check_rater_asymmetry()` and the new `check_asymmetry(ratings = Y)`
# computes the cross-coefficient panel spread (paper §3.2). The OLD
# per-rater calls below still resolve via a soft-deprecation route inside
# `check_asymmetry()`, so these tests continue to exercise the per-rater
# diagnostic via the same call sites. The new ratings-input path is
# exercised in test-check_asymmetry.R; the new grass_card flow lives in
# test-grass_report-card.R.
#
# 0.7.x: `classify()`, `emr_panel()`, and `grass_use_case_ladder()` were
# removed from the package outright (never released on CRAN; no
# deprecation shims), and their tests were removed with them. Only the
# per-rater soft-deprecation route below remains under test.

test_that("check_asymmetry tiers on delta_hat thresholds", {
  ok   <- check_asymmetry(se = c(0.86, 0.88, 0.84), sp = c(0.85, 0.87, 0.86))
  caut <- check_asymmetry(se = c(0.90, 0.88, 0.92), sp = c(0.82, 0.86, 0.85))
  unsf <- check_asymmetry(se = c(0.95, 0.93, 0.94), sp = c(0.78, 0.80, 0.79))

  expect_s3_class(ok, "grass_asymmetry")
  expect_equal(ok$regime,   "ok")
  expect_equal(caut$regime, "caution")
  expect_equal(unsf$regime, "unsafe")
  expect_equal(ok$tier,   1L)
  expect_equal(caut$tier, 2L)
  expect_equal(unsf$tier, 3L)
})

test_that("check_asymmetry 'max' is more conservative than 'mean'", {
  se <- c(0.92, 0.85, 0.85)
  sp <- c(0.80, 0.85, 0.85)
  mx <- check_asymmetry(se, sp, summary = "max")
  mn <- check_asymmetry(se, sp, summary = "mean")
  expect_equal(mx$delta_hat, 0.12, tolerance = 1e-8)
  expect_equal(mn$delta_hat, 0.04, tolerance = 1e-8)
  expect_equal(mx$regime, "unsafe")
  expect_equal(mn$regime, "ok")
})

test_that("check_asymmetry tiers at values safely inside each band", {
  r1 <- check_asymmetry(se = 0.90, sp = 0.86) # gap = 0.04 -> ok
  r2 <- check_asymmetry(se = 0.90, sp = 0.84) # gap = 0.06 -> caution
  r3 <- check_asymmetry(se = 0.93, sp = 0.82) # gap = 0.11 -> unsafe
  expect_equal(r1$regime, "ok")
  expect_equal(r2$regime, "caution")
  expect_equal(r3$regime, "unsafe")
})

test_that("check_asymmetry rejects bad input", {
  expect_error(check_asymmetry(se = 0.8, sp = c(0.8, 0.9)))
  expect_error(check_asymmetry(se = 1.1, sp = 0.5))
  expect_error(check_asymmetry(se = NA, sp = 0.5))
  expect_error(check_asymmetry(se = 0.8, sp = 0.8,
                               threshold_caution = 0.1,
                               threshold_unsafe  = 0.05))
})
