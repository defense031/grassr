# NOTE (v0.2.0): `check_asymmetry(se = ..., sp = ...)` was renamed to
# `check_rater_asymmetry()` and the new `check_asymmetry(ratings = Y)`
# computes the cross-coefficient panel spread (paper §3.2). The OLD
# per-rater calls below still resolve via a soft-deprecation route inside
# `check_asymmetry()`, so these tests continue to exercise the per-rater
# diagnostic via the same call sites. The new ratings-input path is
# exercised in test-check_asymmetry.R; the new grass_card flow lives in
# test-grass_report-card.R.
#
# Phase 5A (2026-05-03): `classify()`, `emr_panel()`, and
# `grass_use_case_ladder()` were unexported and marked retired in
# grass 0.2.0. The arithmetic tests below still cover the bodies so the
# code remains audit-traceable; calls are wrapped in `suppressWarnings()`
# to swallow the per-call `.Deprecated()` notice. See
# grass/design/v0.2.0_paper_alignment.md.

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

test_that("emr_panel matches pbinom closed form at odd k", {
  # P(panel wrong) = P(Bin(k, q) < ceiling(k/2)).
  suppressWarnings({
    expect_equal(emr_panel(q = 0.85, k = 5),
                 pbinom(2, size = 5, prob = 0.85))
    qs <- c(0.80, 0.85, 0.90)
    expect_equal(emr_panel(q = qs, k = 3),
                 pbinom(1, size = 3, prob = qs))
    # Monotone decreasing in q (higher q -> lower panel error)
    expect_true(all(diff(emr_panel(q = c(0.7, 0.8, 0.9), k = 5)) < 0))
    # Bounded in [0, 1]
    expect_true(emr_panel(q = 0, k = 5) == 1)
    expect_true(emr_panel(q = 1, k = 5) == 0)
  })
})

test_that("emr_panel rejects even k", {
  suppressWarnings({
    expect_error(emr_panel(q = 0.85, k = 4))
    expect_error(emr_panel(q = 0.85, k = 2))
    expect_error(emr_panel(q = 0.85, k = 5.5))
  })
})

test_that("classify partitions into five mutually exclusive tiers", {
  suppressWarnings({
    # Screening: T = 0.10, T/2 = 0.05, T_next = 0.05 (diagnostic)
    # U = 0.03  -> Exceeds (< 0.05 = T_next)
    exc <- classify(emr_upper = 0.03, emr_lower = 0.01, use_case = "screening")
    expect_equal(exc$tier, "Exceeds")

    # Diagnostic: T = 0.05, T/2 = 0.025, T_next = 0.02 (clinical)
    # U = 0.022 -> Met with distinction (< T/2 = 0.025, not < T_next = 0.02)
    mwd <- classify(emr_upper = 0.022, emr_lower = 0.010, use_case = "diagnostic")
    expect_equal(mwd$tier, "Met with distinction")

    # Screening with U = 0.07: T/2 = 0.05 <= 0.07 < T = 0.10 -> Meets
    mts <- classify(emr_upper = 0.07, emr_lower = 0.05, use_case = "screening")
    expect_equal(mts$tier, "Meets")

    # Screening with U = 0.12, L = 0.08 (straddles T = 0.10) -> Indeterminate
    ind <- classify(emr_upper = 0.12, emr_lower = 0.08, use_case = "screening")
    expect_equal(ind$tier, "Indeterminate")

    # Screening with L = 0.11 >= T = 0.10 -> Fails
    fl <- classify(emr_upper = 0.15, emr_lower = 0.11, use_case = "screening")
    expect_equal(fl$tier, "Fails")
  })
})

test_that("classify collapses Fails into Indeterminate without a lower bound", {
  suppressWarnings({
    r <- classify(emr_upper = 0.15, emr_lower = NA, use_case = "screening")
    expect_equal(r$tier, "Indeterminate")
  })
})

test_that("classify 'clinical' has no Exceeds (top of ladder)", {
  # Clinical: T = 0.02, T/2 = 0.01, T_next = NA
  # Any U < T/2 is Met with distinction, not Exceeds
  suppressWarnings({
    r <- classify(emr_upper = 0.005, emr_lower = 0.001, use_case = "clinical")
    expect_true(is.na(r$tolerance_next))
    expect_equal(r$tier, "Met with distinction")
  })
})

test_that("classify 'custom' tolerance works and requires T", {
  suppressWarnings({
    r <- classify(emr_upper = 0.03, use_case = "custom", tolerance = 0.08)
    expect_equal(r$tolerance, 0.08)
    expect_true(is.na(r$tolerance_next))
    expect_equal(r$tier, "Met with distinction") # 0.03 < 0.04 = T/2
    # Boundary: U = T/2 exactly lands in Meets (strict <)
    r2 <- classify(emr_upper = 0.04, use_case = "custom", tolerance = 0.08)
    expect_equal(r2$tier, "Meets")
    expect_error(classify(emr_upper = 0.04, use_case = "custom"))
  })
})

test_that("classify attaches Column A (regime) when supplied", {
  suppressWarnings({
    asym <- check_asymmetry(se = c(0.88, 0.86), sp = c(0.85, 0.87))
    r <- classify(emr_upper = 0.06, emr_lower = 0.04,
                  use_case = "screening", regime = asym)
    expect_identical(r$regime, asym)
    df <- as.data.frame(r)
    expect_equal(df$column_a_regime, "ok")
    expect_equal(df$column_a_tier, 1L)
    expect_error(classify(emr_upper = 0.06, use_case = "screening",
                          regime = list(bogus = TRUE)))
  })
})

test_that("grass_use_case_ladder returns the documented four rows", {
  suppressWarnings({
    ld <- grass_use_case_ladder()
    expect_equal(nrow(ld), 4L)
    expect_setequal(ld$use_case,
                    c("research", "screening", "diagnostic", "clinical"))
    expect_equal(ld$tolerance[ld$use_case == "screening"], 0.10)
    expect_equal(ld$tolerance[ld$use_case == "diagnostic"], 0.05)
  })
})
