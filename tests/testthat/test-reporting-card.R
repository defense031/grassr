# `check_asymmetry(ratings = Y)` computes the cross-coefficient panel
# spread (paper section 3.2). The pre-0.2.0 per-rater `se =`/`sp =` call
# form errors (0.8.0). The ratings-input path is exercised in
# test-check_asymmetry.R; the grass_card flow lives in
# test-grass_report-card.R.
#
# 0.7.x: `classify()`, `emr_panel()`, and `grass_use_case_ladder()` were
# removed from the package outright (never released on CRAN; no
# deprecation shims), and their tests were removed with them. Only the
# per-rater soft-deprecation route below remains under test.

test_that("legacy check_asymmetry(se, sp) errors with a pointer to latent_class_fit", {
  expect_error(check_asymmetry(se = c(0.95, 0.93, 0.94), sp = c(0.78, 0.80, 0.79)),
               "latent_class_fit")
})

test_that("check_asymmetry rejects bad input", {
  expect_error(check_asymmetry(se = 0.8, sp = c(0.8, 0.9)))
  expect_error(check_asymmetry(se = 1.1, sp = 0.5))
  expect_error(check_asymmetry(se = NA, sp = 0.5))
  expect_error(check_asymmetry(se = 0.8, sp = 0.8), "latent_class_fit")
})
