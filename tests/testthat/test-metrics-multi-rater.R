# Tests for the internal multi-rater observed-metric ports in
# R/metrics-multi-rater.R. Cross-checks against `irr` and `irrCAC` skip
# when those packages are not installed.

# ---- Fixtures --------------------------------------------------------------

make_Y_basic <- function(seed = 1L, N = 100L, k = 5L, p = 0.3) {
  set.seed(seed)
  matrix(rbinom(N * k, 1, p), nrow = N, ncol = k)
}

# Fully crossed, balanced prevalence, q = 0.85 (truth->vote agreement). Build
# subject-specific true-class labels at p = 0.5, then flip each rater
# independently with prob 1 - q.
make_Y_q85 <- function(seed = 7L, N = 500L, k = 5L) {
  set.seed(seed)
  C  <- rbinom(N, 1, 0.5)
  q  <- 0.85
  Y  <- matrix(0L, nrow = N, ncol = k)
  for (j in seq_len(k)) {
    flips <- rbinom(N, 1, 1 - q)
    Y[, j] <- as.integer(ifelse(flips == 1, 1 - C, C))
  }
  Y
}

# ---- obs_fleiss_kappa ------------------------------------------------------

test_that("obs_fleiss_kappa returns a finite scalar", {
  Y <- make_Y_basic()
  v <- obs_fleiss_kappa(Y)
  expect_type(v, "double")
  expect_length(v, 1L)
  expect_true(is.finite(v))
})

test_that("obs_fleiss_kappa matches irr::kappam.fleiss when available", {
  skip_if_not_installed("irr")
  Y <- make_Y_basic()
  ours <- obs_fleiss_kappa(Y)
  ref  <- irr::kappam.fleiss(as.data.frame(Y))$value
  expect_lt(abs(ours - ref), 0.001)
})

# ---- obs_mean_pairwise_pabak ----------------------------------------------

test_that("obs_mean_pairwise_pabak returns a finite scalar", {
  Y <- make_Y_basic()
  v <- obs_mean_pairwise_pabak(Y)
  expect_type(v, "double")
  expect_length(v, 1L)
  expect_true(is.finite(v))
})

test_that("obs_mean_pairwise_pabak reduces to 2*P0 - 1 at k = 2", {
  set.seed(2)
  Y <- matrix(rbinom(200, 1, 0.4), nrow = 100, ncol = 2)
  P0 <- mean(Y[, 1] == Y[, 2])
  expect_equal(obs_mean_pairwise_pabak(Y), 2 * P0 - 1, tolerance = 1e-12)
})

# ---- obs_mean_pairwise_ac1 -------------------------------------------------

test_that("obs_mean_pairwise_ac1 returns a finite scalar", {
  Y <- make_Y_basic()
  v <- obs_mean_pairwise_ac1(Y)
  expect_type(v, "double")
  expect_length(v, 1L)
  expect_true(is.finite(v))
})

test_that("obs_mean_pairwise_ac1 matches irrCAC::gwet.ac1.raw when available", {
  skip_if_not_installed("irrCAC")
  Y <- make_Y_basic()
  ours <- obs_mean_pairwise_ac1(Y)
  ref  <- irrCAC::gwet.ac1.raw(Y)$est$coeff.val
  expect_lt(abs(ours - ref), 0.005)
})

# ---- obs_krippendorff_alpha ------------------------------------------------

test_that("obs_krippendorff_alpha returns a finite scalar on non-degenerate input", {
  Y <- make_Y_basic()
  v <- obs_krippendorff_alpha(Y)
  expect_type(v, "double")
  expect_length(v, 1L)
  expect_true(is.finite(v))
})

test_that("obs_krippendorff_alpha returns NA on all-0 / all-1 input", {
  Y <- matrix(0L, nrow = 50L, ncol = 4L)
  v <- obs_krippendorff_alpha(Y)
  expect_true(is.na(v))
  expect_match(attr(v, "note"), "undefined")
})

# ---- obs_icc_glmer ---------------------------------------------------------

test_that("obs_icc_glmer returns numeric in (0, 1) when lme4 available, else NA + note", {
  # Use a fixture with real subject heterogeneity (q = 0.85) so glmer fits
  # a strictly positive subject variance. On pure binary noise the variance
  # estimate hits the zero boundary, which is a legitimate fit result but
  # collides with the strict (0, 1) bound in the design spec.
  Y <- make_Y_q85()
  v <- obs_icc_glmer(Y)
  if (requireNamespace("lme4", quietly = TRUE)) {
    if (is.na(v)) {
      # glmer can fail to converge on small N / extreme prevalence — accept
      # the documented NA + note path here.
      expect_match(attr(v, "note"), "glmer")
    } else {
      expect_true(is.finite(v))
      expect_gt(v, 0)
      expect_lt(v, 1)
    }
  } else {
    expect_true(is.na(v))
    expect_equal(attr(v, "note"), "lme4 not available")
  }
})

# ---- Sanity: panel values in (-0.05, 1.0) at q = 0.85 ---------------------

test_that("all observed metrics fall in (-0.05, 1.0) at q = 0.85, balanced prevalence, 5x500", {
  Y <- make_Y_q85()
  vals <- c(
    obs_fleiss_kappa(Y),
    obs_mean_pairwise_pabak(Y),
    obs_mean_pairwise_ac1(Y),
    obs_krippendorff_alpha(Y)
  )
  for (v in vals) {
    expect_true(is.finite(v))
    expect_gt(v, -0.05)
    expect_lt(v, 1.0)
  }
  if (requireNamespace("lme4", quietly = TRUE)) {
    icc <- obs_icc_glmer(Y)
    if (!is.na(icc)) {
      expect_gt(icc, -0.05)
      expect_lt(icc, 1.0)
    }
  }
})

# ---- compute_observed dispatch --------------------------------------------

test_that("compute_observed dispatches to the right obs_* function", {
  Y <- make_Y_basic()
  expect_equal(compute_observed("pabak", Y),          obs_mean_pairwise_pabak(Y))
  expect_equal(compute_observed("fleiss_kappa", Y),   obs_fleiss_kappa(Y))
  expect_equal(compute_observed("mean_ac1", Y),       obs_mean_pairwise_ac1(Y))
  expect_equal(compute_observed("krippendorff_a", Y), obs_krippendorff_alpha(Y))
})

test_that("compute_observed errors on unknown metric", {
  Y <- make_Y_basic()
  expect_error(compute_observed("nope", Y))
})
