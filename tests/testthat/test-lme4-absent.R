# Guards the 0.7.0 regression where grass_report() hard-errored when the
# ICC coefficient could not be computed (lme4 absent / glmer failure).
# Under v0.7.1 a non-finite panel coefficient must be DROPPED with a note,
# never propagated into position_on_surface() as an NA (which errors).
#
# We simulate lme4 absence WITHOUT uninstalling it: obs_icc_glmer() is the
# package's ICC computer (it requireNamespace()s lme4 internally); mocking
# it to return NA reproduces the "ICC uncomputable" state deterministically
# and independently of whether lme4 is on the test machine.

test_that("grass_report drops the ICC row (with a note) when ICC is uncomputable", {
  skip_if_not_installed("testthat")  # local_mocked_bindings support
  set.seed(1L)
  Y <- matrix(rbinom(5 * 200, 1, 0.3), nrow = 200, ncol = 5)

  testthat::local_mocked_bindings(
    obs_icc_glmer = function(...) NA_real_,
    .package = "grassr"
  )

  card <- grass_report(Y, bootstrap_B = 0)

  # 1. It succeeds (does not hard-error) and returns a grass_card.
  expect_s3_class(card, "grass_card")

  # 2. The panel lacks the icc row.
  expect_false("icc" %in% card$panel$coefficient)
  expect_true(all(c("pabak", "mean_ac1", "fleiss_kappa") %in%
                  card$panel$coefficient))

  # 3. A note names the dropped coefficient.
  expect_true(any(grepl("dropped from the panel", card$notes, fixed = TRUE)))
  expect_true(any(grepl("icc", card$notes, fixed = TRUE)))
})

test_that("check_asymmetry also degrades gracefully when ICC is uncomputable", {
  # The 0.7.0 fix landed first in check_asymmetry(); pin it here too so the
  # ratings-input diagnostic never hard-errors on a missing ICC.
  set.seed(2L)
  Y <- matrix(rbinom(5 * 200, 1, 0.3), nrow = 200, ncol = 5)

  testthat::local_mocked_bindings(
    obs_icc_glmer = function(...) NA_real_,
    .package = "grassr"
  )

  out <- check_asymmetry(Y)
  expect_s3_class(out, "grass_asymmetry_panel")
  expect_false("icc" %in% out$panel$coefficient)
  expect_true(is.finite(out$delta_hat))
})
