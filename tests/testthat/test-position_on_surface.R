# Tests for grass::position_on_surface() — the Target-2 reporting primitive.
# See grass/R/position_on_surface.R for the spec.

test_that("returns a grass_surface_position with the documented fields", {
  r <- position_on_surface(
    obs_value = 0.62,
    metric    = "pabak",
    pi_hat    = 0.42,
    k         = 5,
    N         = 50
  )
  expect_s3_class(r, "grass_surface_position")
  expected <- c("observed_value", "metric", "design", "q_hat", "se_q_hat",
                "percentile", "band_probabilities", "modal_band",
                "modal_band_label", "confidence", "sampling_method", "notes")
  expect_true(all(expected %in% names(r)))
  expect_equal(r$observed_value, 0.62)
  expect_equal(r$metric, "pabak")
  expect_equal(r$design, list(pi_hat = 0.42, k = 5L, N = 50L))
})

test_that("q_hat for PABAK matches the closed-form inverse 0.5 * (1 + sqrt(PABAK))", {
  # PABAK = (2q - 1)^2  =>  q = 0.5 * (1 + sqrt(PABAK))
  for (obs in c(0.20, 0.50, 0.80, 0.95)) {
    r <- position_on_surface(obs, "pabak", pi_hat = 0.5,
                             k = 5, N = 200)
    expect_equal(r$q_hat, 0.5 * (1 + sqrt(obs)),
                 tolerance = 2e-3,
                 info = sprintf("obs = %.2f", obs))
  }
})

test_that("all four agreement-family metrics invert to q_hat in [0.5, 1]", {
  for (m in c("pabak", "fleiss_kappa", "mean_ac1", "krippendorff_a")) {
    r <- position_on_surface(0.55, m, pi_hat = 0.3, k = 5, N = 100)
    expect_true(is.finite(r$q_hat),
                info = paste("metric =", m))
    expect_gte(r$q_hat, 0.5 - 1e-6)
    expect_lte(r$q_hat, 1.0 + 1e-6)
  }
})

test_that("band probabilities sum to 1 and the modal label matches", {
  r <- position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200)
  expect_equal(sum(r$band_probabilities), 1, tolerance = 1e-8)
  expect_equal(names(r$band_probabilities),
               c("Poor", "Moderate", "Strong", "Excellent"))
  expect_equal(r$modal_band_label,
               names(r$band_probabilities)[r$modal_band])
})

test_that("confidence qualifier thresholds on p_star are respected", {
  # Large N drives se_q_hat down -> decisive
  big <- position_on_surface(0.64, "pabak", pi_hat = 0.5, k = 15, N = 1000)
  expect_equal(big$confidence, "decisive")
  # Small N widens q_hat sampling distribution -> weak
  small <- position_on_surface(0.64, "pabak", pi_hat = 0.5, k = 3, N = 30)
  expect_true(small$confidence %in% c("weak", "moderate"))
  p_star <- max(small$band_probabilities)
  expect_lt(p_star, 0.90)
})

test_that("percentile lies in [0, 1] and is approximately 0.5 at the surface mean", {
  # At obs_value == E[metric](q_hat) the delta-method percentile should
  # sit near 0.5 by symmetry of the normal approximation.
  r <- position_on_surface(0.49, "pabak", pi_hat = 0.5, k = 5, N = 200)
  expect_gte(r$percentile, 0)
  expect_lte(r$percentile, 1)
  expect_equal(r$percentile, 0.5, tolerance = 0.1)
})

test_that("high q_hat lands in Excellent band with a decisive qualifier at large N", {
  r <- position_on_surface(0.92, "pabak", pi_hat = 0.5, k = 15, N = 1000)
  expect_equal(r$modal_band_label, "Excellent")
  expect_equal(r$confidence, "decisive")
})

test_that("low PABAK at balanced prevalence lands in a low-quality band", {
  # PABAK = 0.02 implies q_hat = 0.5 * (1 + sqrt(0.02)) ~ 0.57
  r <- position_on_surface(0.02, "pabak", pi_hat = 0.5, k = 5, N = 200)
  expect_true(r$modal_band_label %in% c("Poor", "Moderate"))
  expect_lt(r$q_hat, 0.625)
})

test_that("obs_value outside achievable range is clamped and flagged in notes", {
  # PABAK tops out at 1. Try an impossible value.
  r <- position_on_surface(1.05, "pabak", pi_hat = 0.5, k = 5, N = 200)
  expect_true(any(grepl("clamped", r$notes)))
  expect_equal(r$q_hat, 1.0, tolerance = 1e-3)
})

test_that("delta-method percentile agrees with empirical percentile from a large normal sample", {
  # Build a fake per_rep matrix consistent with the delta-method normal
  # approximation and check the two routes land at similar percentiles.
  # Pick an obs_value deliberately off the surface mean so qnorm is
  # well-defined. This test pins the caller-supplied per_rep path; the
  # bundled empirical q_hat surface is a separate code path covered by the
  # T3-specific tests below.
  set.seed(7L)
  r_delta <- position_on_surface(0.55, "pabak", pi_hat = 0.5,
                                 k = 5, N = 500, method = "delta")
  mu_m <- (2 * r_delta$q_hat - 1)^2
  # Use a sd informed by the delta-method implied se_q_hat and dE/dq.
  # For PABAK: dE/dq = 2*(2q-1)*2 = 4*(2q-1). sd_metric = se_q_hat * |dE/dq|.
  dEdq <- 4 * (2 * r_delta$q_hat - 1)
  implied_sd <- r_delta$se_q_hat * abs(dEdq)
  expect_true(is.finite(implied_sd) && implied_sd > 0)
  fake_reps <- stats::rnorm(10000, mean = mu_m, sd = implied_sd)
  r_emp <- position_on_surface(
    0.55, "pabak", pi_hat = 0.5, k = 5, N = 500,
    method = "empirical",
    surface_data = list(per_rep = fake_reps)
  )
  expect_equal(r_emp$sampling_method, "empirical")
  # Percentile from the caller-supplied per_rep should agree with the
  # delta-method pnorm percentile under matching mean / sd.
  expect_equal(r_emp$percentile, r_delta$percentile, tolerance = 0.05)
})

test_that("custom bands and band_labels are honored", {
  r <- position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200,
                           bands = c(0.5, 0.7, 0.85, 0.95, 1.0),
                           band_labels = c("A", "B", "C", "D"))
  expect_equal(names(r$band_probabilities), c("A", "B", "C", "D"))
  expect_equal(sum(r$band_probabilities), 1, tolerance = 1e-8)
  expect_true(r$modal_band_label %in% c("A", "B", "C", "D"))
})

test_that("input validation rejects bad inputs with informative stops", {
  expect_error(position_on_surface(0.5, "bogus_metric", 0.5, 5, 100),
               "metric")
  expect_error(position_on_surface("not-numeric", "pabak", 0.5, 5, 100),
               "obs_value")
  expect_error(position_on_surface(0.5, "pabak", 0, 5, 100),
               "pi_hat")
  expect_error(position_on_surface(0.5, "pabak", 1, 5, 100),
               "pi_hat")
  expect_error(position_on_surface(0.5, "pabak", 0.5, 1, 100),
               "k")
  expect_error(position_on_surface(0.5, "pabak", 0.5, 5, 0),
               "N")
  expect_error(position_on_surface(0.5, "pabak", 0.5, 5, 100,
                                   bands = c(0.5, 0.4, 0.75, 0.875, 1.0)),
               "increasing")
  expect_error(position_on_surface(0.5, "pabak", 0.5, 5, 100,
                                   bands = c(0.4, 0.6, 0.75, 0.875, 1.0)),
               "0\\.5")
  expect_error(position_on_surface(0.5, "pabak", 0.5, 5, 100,
                                   band_labels = c("a", "b", "c")),
               "band_labels")
})

test_that("empirical method for ICC uses bundled reference curves", {
  # Default ICC path uses the fitted reference (GLMM-gap corrected). The
  # bundled sysdata supplies the nearest (F_key, k, N) cell so practitioners
  # do not need to hand-build a reference curve.
  r <- position_on_surface(0.3, "icc", 0.5, 5, 100, method = "empirical")
  expect_s3_class(r, "grass_surface_position")
  expect_true(is.finite(r$q_hat))
  # Fitted-reference note should identify the nearest cell.
  expect_true(any(grepl("Fitted-ICC reference", r$notes)))

  # Oracle reference_type option still works.
  r_or <- position_on_surface(0.3, "icc", 0.5, 5, 100,
                              method = "empirical",
                              reference_type = "oracle")
  expect_s3_class(r_or, "grass_surface_position")
  expect_true(any(grepl("bundled F_key", r_or$notes)))
})

# ---- T3 empirical-band sysdata backfill tests ------------------------------

test_that("empirical method returns sensible bands for PABAK=0.62 at k=5,N=200,pi=0.5", {
  r <- position_on_surface(obs_value = 0.62, metric = "pabak",
                           pi_hat = 0.5, k = 5, N = 200,
                           method = "empirical")
  expect_equal(r$sampling_method, "empirical")
  expect_equal(sum(r$band_probabilities), 1, tolerance = 1e-6)
  expect_true(all(r$band_probabilities >= 0))
  # q_hat for PABAK=0.62 is 0.5*(1+sqrt(0.62)) ~ 0.894, so modal band at
  # k=5,N=200 should be Excellent or Strong (tight sampling distribution,
  # but q_hat is right at the 0.875 boundary so either is plausible).
  expect_true(r$modal_band_label %in% c("Strong", "Excellent"))
  expect_true(r$q_hat > 0.85 && r$q_hat < 0.94)
})

test_that("empirical and delta methods agree within +/- 0.15 on modal band probability", {
  r_e <- position_on_surface(0.62, "pabak", 0.5, 5, 200, method = "empirical")
  r_d <- position_on_surface(0.62, "pabak", 0.5, 5, 200, method = "delta")
  # Both land in the same broad region of the surface.
  expect_true(abs(max(r_e$band_probabilities) -
                  max(r_d$band_probabilities)) <= 0.15)
  # Modal band adjacent or equal.
  expect_true(abs(r_e$modal_band - r_d$modal_band) <= 1)
})

test_that("nearest-neighbor clamping for out-of-grid design is flagged in notes", {
  # k=4 is not in the sim grid {3, 5, 8, 15}; should clamp to 3 or 5 and flag.
  # N=1234 is not in {50, 200, 1000}; should clamp to 1000 and flag.
  r <- position_on_surface(0.62, "pabak", 0.5, k = 4, N = 1234,
                           method = "empirical")
  expect_true(any(grepl("clamped", r$notes)))
  expect_true(any(grepl("k=", r$notes) | grepl("N=", r$notes)))
})

test_that("ICC at k beyond fitted grid falls back to oracle reference with a prominent note", {
  # The v0.4 fitted-ICC bundle covers k in {3, 5, 8, 15, 25}. Beyond
  # max(k_grid) + 1 = 26 the GLMM-gap correction is not available; the
  # function should fall back to the oracle reference and surface a
  # user-visible note that names the k gap, rather than silently
  # clamping to the fitted reference at k=25.
  r <- position_on_surface(0.30, metric = "icc",
                           pi_hat = 0.5, k = 50, N = 100,
                           method = "empirical")
  expect_s3_class(r, "grass_surface_position")
  fallback_note <- grep("Fitted-ICC reference unavailable", r$notes,
                       value = TRUE)
  expect_length(fallback_note, 1L)
  expect_match(fallback_note, "k=50")
  # And no fitted-ICC reference cell note (which would indicate silent clamp).
  expect_false(any(grepl("Fitted-ICC reference \\(GLMM-gap corrected\\)",
                         r$notes)))
})

test_that("caller-supplied surface_data$per_rep still overrides bundled lookup", {
  # Force per_rep with a tight distribution around 0.90 so percentile of 0.62
  # collapses to 0 (well below the rep cloud).
  set.seed(11L)
  fake_reps <- stats::rnorm(5000, mean = 0.90, sd = 0.02)
  r <- position_on_surface(0.62, "pabak", 0.5, 5, 200,
                           method = "empirical",
                           surface_data = list(per_rep = fake_reps))
  expect_equal(r$sampling_method, "empirical")
  expect_lt(r$percentile, 0.01)
})

test_that("q_grid_per_rep hook is preserved (not consumed) for empirical callers", {
  q_grid <- seq(0.5, 1.0, length.out = 11L)
  fake_qgrid_reps <- matrix(stats::runif(11 * 100), nrow = 11, ncol = 100)
  r <- position_on_surface(0.62, "pabak", 0.5, 5, 200,
                           method = "empirical",
                           surface_data = list(
                             q_grid = q_grid,
                             q_grid_per_rep = fake_qgrid_reps
                           ))
  expect_true(any(grepl("q_grid_per_rep", r$notes)))
})

test_that("icc glmer path: ratings matrix pins down F_key via (mu, tau2)", {
  skip_if_not_installed("lme4")
  # Simulate a rating matrix with known mu and tau2. Draw N subject-specific
  # propensities from N(mu, tau2) on the logit scale; each rater's binary
  # rating is drawn Bernoulli(plogis(u_i)). With symmetric raters at q=1
  # (noise-free), obs_ICC is fully determined by tau2.
  set.seed(42L)
  N <- 80L
  k <- 5L
  mu_true   <- -0.5
  tau2_true <- 1.0
  u <- rnorm(N, mean = mu_true, sd = sqrt(tau2_true))
  p <- plogis(u)
  ratings <- matrix(rbinom(N * k, size = 1, prob = rep(p, each = k)),
                    nrow = k, ncol = N)
  # Empirical pi_hat from the rating matrix
  pi_hat <- mean(ratings)
  # Fit glmer via the helper directly to check it converges
  fit <- grass:::fit_tau2_from_ratings(ratings)
  expect_true(is.finite(fit$mu))
  expect_true(is.finite(fit$tau2) && fit$tau2 > 0)
  # mu_hat should be near mu_true (within ~0.3 given finite N)
  expect_lt(abs(fit$mu - mu_true), 0.5)
  # Now call position_on_surface with ratings; note should mention glmer
  r <- position_on_surface(
    obs_value = 0.55, metric = "icc",
    pi_hat = pi_hat, k = k, N = N,
    ratings = ratings
  )
  expect_s3_class(r, "grass_surface_position")
  expect_true(any(grepl("glmer", r$notes)))
  # Without ratings, the fallback note should mention nearest-M1 instead
  r_fallback <- position_on_surface(
    obs_value = 0.55, metric = "icc",
    pi_hat = pi_hat, k = k, N = N
  )
  expect_true(any(grepl("nearest-M1", r_fallback$notes)))
})

test_that("icc glmer fallback on tiny N or non-numeric ratings surfaces a note", {
  skip_if_not_installed("lme4")
  # Tiny N: helper should return NA estimates with explanatory note.
  tiny <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
  fit_tiny <- grass:::fit_tau2_from_ratings(tiny)
  expect_true(!is.finite(fit_tiny$mu))
  expect_true(grepl("N>=10", fit_tiny$note))
  # Bad input: character rating matrix.
  fit_bad <- grass:::fit_tau2_from_ratings(list(foo = 1))
  expect_true(grepl("k x N integer matrix", fit_bad$note))
})

test_that("icc resolves via bundled sysdata; caller reference_curve wins when supplied", {
  # Default path: bundled icc_reference_curves handles the inversion.
  r <- position_on_surface(0.3, "icc", 0.5, 5, 100)
  expect_s3_class(r, "grass_surface_position")
  expect_true(is.finite(r$q_hat))
  # Extreme pi_hat far from any sim F_key's M1 produces either a coarse-lookup
  # note or identifies the fitted-reference F_key cell in notes.
  r_far <- position_on_surface(0.3, "icc", pi_hat = 0.01, k = 5, N = 100)
  expect_true(any(grepl("coarse", r_far$notes) |
                  grepl("Fitted-ICC reference", r_far$notes)))
  # Caller-supplied reference curve of the right length still wins.
  q_grid <- seq(0.5, 1.0, length.out = 501)
  fake_icc <- (q_grid - 0.5) / 0.5 * 0.6
  r_custom <- position_on_surface(0.3, "icc", 0.5, 5, 200,
                                  surface_data = list(reference_curve = fake_icc))
  expect_s3_class(r_custom, "grass_surface_position")
  expect_true(is.finite(r_custom$q_hat))
  expect_true(any(grepl("supplied by caller", r_custom$notes)))
})

test_that("as.data.frame() flattens a surface_position into a one-row frame", {
  r <- position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200)
  df <- as.data.frame(r)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1L)
  expect_true(all(c("metric", "observed_value", "pi_hat", "k", "N",
                    "q_hat", "se_q_hat", "percentile", "modal_band",
                    "modal_band_label", "confidence", "sampling_method",
                    "p_band_1", "p_band_2", "p_band_3", "p_band_4") %in%
                  names(df)))
  expect_equal(sum(df[, paste0("p_band_", 1:4)]), 1, tolerance = 1e-8)
})

# ---- v0.2.0 ratings-primary path ------------------------------------------
# The ratings-primary entry: `position_on_surface(ratings = Y, metric = m)`
# auto-derives obs_value, pi_hat, k, N from an N x k binary matrix. The
# scalar-input path is preserved exactly; round-trip tests below pin the
# floating-point equivalence.

test_that("ratings-primary path round-trips against scalar path for non-ICC metrics", {
  set.seed(1L)
  Y <- matrix(rbinom(500, 1, 0.30), nrow = 100, ncol = 5)
  for (m in c("pabak", "fleiss_kappa", "mean_ac1", "krippendorff_a")) {
    a <- position_on_surface(ratings = Y, metric = m)
    b <- position_on_surface(obs_value = compute_observed(m, Y),
                             pi_hat   = mean(Y),
                             k        = ncol(Y),
                             N        = nrow(Y),
                             metric   = m)
    expect_equal(a$percentile, b$percentile, tolerance = 1e-10,
                 info = paste("metric =", m))
    expect_equal(a$q_hat, b$q_hat, tolerance = 1e-10,
                 info = paste("metric =", m))
  }
})

test_that("ratings-primary path round-trips against scalar path for ICC", {
  skip_if_not_installed("lme4")
  set.seed(3L)
  Y <- matrix(rbinom(500, 1, 0.30), nrow = 100, ncol = 5)
  obs_icc <- compute_observed("icc", Y)
  skip_if(!is.finite(as.numeric(obs_icc)),
          "obs_icc unavailable in this environment")
  a <- position_on_surface(ratings = Y, metric = "icc")
  # Scalar-path equivalent: pass the same N x k Y as `ratings` for the
  # glmer fit (transposed to k x N to match the legacy convention).
  b <- position_on_surface(obs_value = as.numeric(obs_icc),
                           pi_hat = mean(Y), k = ncol(Y), N = nrow(Y),
                           metric = "icc",
                           ratings = t(Y))
  # GLMM fits hit the same data; expected agreement at near-machine precision.
  expect_equal(a$percentile, b$percentile, tolerance = 1e-6)
  expect_equal(a$q_hat, b$q_hat, tolerance = 1e-6)
})

test_that("ratings-primary path works for k = 2", {
  set.seed(2L)
  Y2 <- matrix(rbinom(400, 1, 0.4), nrow = 200, ncol = 2)
  r <- position_on_surface(ratings = Y2, metric = "pabak")
  expect_s3_class(r, "grass_surface_position")
  expect_true(is.finite(r$percentile))
  expect_equal(r$design$k, 2L)
  expect_equal(r$design$N, 200L)
})

test_that("missing both ratings and complete scalar set errors helpfully", {
  expect_error(position_on_surface(metric = "pabak"),
               "Either supply")
})

test_that("non-binary ratings input fails with an informative error", {
  # Includes the value 2 — must be rejected, not silently coerced.
  expect_error(
    position_on_surface(ratings = matrix(c(0, 1, 2, 1), nrow = 2),
                        metric = "pabak"),
    regexp = NULL
  )
})
