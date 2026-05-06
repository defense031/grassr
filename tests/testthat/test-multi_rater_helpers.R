# Tests for normalize_ratings(), compute_panel(), pick_primary_coefficient()
# in R/multi_rater_helpers.R.

# ---- Fixtures --------------------------------------------------------------

make_Y_k5 <- function(seed = 1L, N = 100L, k = 5L, p = 0.3) {
  set.seed(seed)
  matrix(rbinom(N * k, 1, p), nrow = N, ncol = k)
}

make_Y_k2 <- function(seed = 4L, N = 80L, p = 0.4) {
  set.seed(seed)
  matrix(rbinom(N * 2L, 1, p), nrow = N, ncol = 2L)
}

# ---- normalize_ratings -----------------------------------------------------

test_that("normalize_ratings accepts an N x k integer matrix and returns canonical form", {
  Y_in <- make_Y_k5()
  Y    <- normalize_ratings(Y_in)
  expect_true(is.matrix(Y))
  expect_equal(storage.mode(Y), "integer")
  expect_equal(dim(Y), c(100L, 5L))
  expect_true(all(Y %in% c(0L, 1L)))
})

test_that("normalize_ratings accepts a data.frame with k binary columns", {
  Y_in <- make_Y_k5()
  df   <- as.data.frame(Y_in)
  Y    <- normalize_ratings(df)
  expect_equal(dim(Y), c(100L, 5L))
  expect_true(all(Y %in% c(0L, 1L)))
  expect_equal(unname(Y), unname(Y_in))
})

test_that("normalize_ratings accepts a list of two equal-length 0/1 vectors (k = 2)", {
  set.seed(11)
  r1 <- rbinom(50, 1, 0.5)
  r2 <- rbinom(50, 1, 0.5)
  Y  <- normalize_ratings(list(r1, r2))
  expect_equal(dim(Y), c(50L, 2L))
  expect_equal(unname(Y[, 1]), as.integer(r1))
  expect_equal(unname(Y[, 2]), as.integer(r2))
})

test_that("normalize_ratings rejects non-binary input", {
  Y_bad <- matrix(c(0, 1, 2, 1), nrow = 2L, ncol = 2L)
  expect_error(normalize_ratings(Y_bad), "0/1")
})

test_that("normalize_ratings rejects NA values", {
  Y_bad <- make_Y_k5()
  Y_bad[1, 1] <- NA
  expect_error(normalize_ratings(Y_bad), "NA")
})

test_that("normalize_ratings rejects k < 2 (single-column matrix)", {
  Y_bad <- matrix(rbinom(20, 1, 0.5), nrow = 20, ncol = 1)
  expect_error(normalize_ratings(Y_bad), "k >= 2")
})

test_that("normalize_ratings rejects unsupported types", {
  expect_error(normalize_ratings("not a matrix"), "unsupported")
  expect_error(normalize_ratings(NULL), "NULL")
})

test_that("normalize_ratings rejects mismatched-length list of two vectors", {
  expect_error(normalize_ratings(list(c(0, 1, 0), c(1, 0))),
               "unequal")
})

# ---- compute_panel ---------------------------------------------------------

test_that("compute_panel at k = 2 returns pabak/ac1/kappa/krippendorff_a", {
  Y <- make_Y_k2()
  panel <- compute_panel(Y)
  expect_named(panel,
               c("pabak", "ac1", "kappa", "krippendorff_a"),
               ignore.order = TRUE)
  for (nm in names(panel)) {
    expect_true(is.finite(panel[[nm]]),
                info = paste("expected finite for", nm))
  }
})

test_that("compute_panel at k = 5 returns pabak/ac1/fleiss_kappa/krippendorff_a/icc", {
  Y <- make_Y_k5()
  panel <- compute_panel(Y)
  expect_named(panel,
               c("pabak", "ac1", "fleiss_kappa", "krippendorff_a", "icc"),
               ignore.order = TRUE)
  expect_true(is.finite(panel$pabak))
  expect_true(is.finite(panel$ac1))
  expect_true(is.finite(panel$fleiss_kappa))
  expect_true(is.finite(panel$krippendorff_a))
  # icc may be NA + note when lme4 unavailable; both branches are documented.
  if (requireNamespace("lme4", quietly = TRUE)) {
    if (!is.na(panel$icc)) expect_true(is.finite(panel$icc))
  } else {
    expect_true(is.na(panel$icc))
  }
})

test_that("compute_panel matches individual obs_* calls at k = 5", {
  Y <- make_Y_k5()
  panel <- compute_panel(Y)
  expect_equal(panel$pabak,          obs_mean_pairwise_pabak(Y))
  expect_equal(panel$ac1,            obs_mean_pairwise_ac1(Y))
  expect_equal(panel$fleiss_kappa,   obs_fleiss_kappa(Y))
  expect_equal(panel$krippendorff_a, obs_krippendorff_alpha(Y))
})

test_that("compute_panel at k = 2 matches the 2x2 metrics-core path", {
  Y <- make_Y_k2()
  panel <- compute_panel(Y)
  tab   <- build_table(Y[, 1], Y[, 2])
  base  <- compute_agreement_metrics(tab)
  expect_equal(panel$pabak, unname(base[["PABAK"]]))
  expect_equal(panel$ac1,   unname(base[["AC1"]]))
  expect_equal(panel$kappa, unname(base[["kappa"]]))
})

# ---- pick_primary_coefficient ---------------------------------------------

test_that("pick_primary_coefficient: k = 2, balanced prevalence -> pabak", {
  expect_equal(pick_primary_coefficient(2, 0.5, "inter"), "pabak")
})

test_that("pick_primary_coefficient: k = 2, low prevalence -> ac1", {
  expect_equal(pick_primary_coefficient(2, 0.05, "inter"), "ac1")
})

test_that("pick_primary_coefficient: k = 2, high prevalence -> ac1", {
  expect_equal(pick_primary_coefficient(2, 0.95, "inter"), "ac1")
})

test_that("pick_primary_coefficient: k = 5, balanced prevalence -> pabak", {
  expect_equal(pick_primary_coefficient(5, 0.5, "inter"), "pabak")
})

test_that("pick_primary_coefficient: k = 5, prevalence < 0.15 -> ac1", {
  expect_equal(pick_primary_coefficient(5, 0.10, "inter"), "ac1")
})

test_that("pick_primary_coefficient: axis = intra always returns icc", {
  expect_equal(pick_primary_coefficient(5, 0.5, "intra"), "icc")
  expect_equal(pick_primary_coefficient(2, 0.05, "intra"), "icc")
})

test_that("pick_primary_coefficient errors on invalid k or pi_hat", {
  expect_error(pick_primary_coefficient(1, 0.5, "inter"), "k")
  expect_error(pick_primary_coefficient(5, 1.5, "inter"), "pi_hat")
  expect_error(pick_primary_coefficient(5, -0.1, "inter"), "pi_hat")
})
