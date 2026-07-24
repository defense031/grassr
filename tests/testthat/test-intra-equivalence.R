# Intra-axis equivalence (v0.7.1, TIER3_DESIGN.md Arm F package half).
#
# Under the reference model -- a stable rater, occasions conditionally
# independent given the latent class -- an N x W intra-rater matrix has the
# same joint distribution as an inter-rater panel of k = W diagonal raters,
# by construction of the DGP (identical per-column Bernoulli(q | C_i) draws,
# conditional independence). The package therefore reads both axes against
# the same calibration; these tests pin the code-path half of that claim:
# the two axes must produce IDENTICAL panel numbers on the same matrix,
# differing only in labeling (primary coefficient, notes).

test_that("intra and inter axes position the same matrix identically", {
  set.seed(20260705)
  truth <- rbinom(120, 1, 0.45)
  Y <- sapply(1:3, function(j) {
    flip <- rbinom(120, 1, 0.12) == 1
    ifelse(flip, 1L - truth, truth)
  })

  card_inter <- suppressWarnings(grass_report(Y, axis = "inter",
                                              bootstrap_B = 0))
  card_intra <- suppressWarnings(grass_report(Y, axis = "intra",
                                              bootstrap_B = 0))

  # Same coefficients, same observed values, same pooled percentiles,
  # same consistency bands, same delta_hat: the intra surface IS the
  # inter diagonal at k = W under the reference model.
  expect_equal(card_intra$panel$coefficient, card_inter$panel$coefficient)
  expect_equal(card_intra$panel$observed_value,
               card_inter$panel$observed_value, tolerance = 1e-12)
  expect_equal(card_intra$panel$surface_percentile,
               card_inter$panel$surface_percentile, tolerance = 1e-12)
  expect_equal(card_intra$panel$band_lo, card_inter$panel$band_lo,
               tolerance = 1e-12)
  expect_equal(card_intra$panel$band_hi, card_inter$panel$band_hi,
               tolerance = 1e-12)
  expect_equal(card_intra$delta$delta_hat, card_inter$delta$delta_hat,
               tolerance = 1e-12)
})

test_that("intra axis note states exactness, not approximation", {
  set.seed(1)
  Y <- matrix(rbinom(200, 1, 0.5), nrow = 50, ncol = 4)
  card <- suppressWarnings(grass_report(Y, axis = "intra", bootstrap_B = 0))
  intra_notes <- grep("Intra-rater axis", card$notes, value = TRUE)
  expect_length(intra_notes, 1L)
  expect_match(intra_notes, "EXACT")
  expect_false(any(grepl("approximation|queued|provisional", intra_notes)))
})
