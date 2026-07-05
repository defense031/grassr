# Tests for grassr::position_on_surface() — the Target-2 reporting primitive.
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
  # v0.7.1 sweep convention: band_probabilities / modal_band /
  # modal_band_label / confidence are retired; the object now carries the
  # pooled percentile, its basis, the consistency band, and the sweep.
  expected <- c("observed_value", "metric", "design", "q_hat", "se_q_hat",
                "percentile", "percentile_basis", "band", "sweep",
                "sampling_method", "reference_used", "notes")
  expect_true(all(expected %in% names(r)))
  # The retired fields must be gone.
  retired <- c("band_probabilities", "modal_band", "modal_band_label",
               "confidence")
  expect_false(any(retired %in% names(r)))
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

test_that("the sweep is a valid p(q) profile and the band is well-formed", {
  r <- position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200)
  expect_s3_class(r$sweep, "data.frame")
  expect_true(all(c("q", "p") %in% names(r$sweep)))
  expect_true(all(r$sweep$q >= 0.5 & r$sweep$q <= 1))
  expect_true(all(r$sweep$p >= 0 & r$sweep$p <= 1))
  # Band is a list with the documented endpoints / flags.
  expect_true(all(c("lo", "hi", "level", "open_low", "open_high") %in%
                  names(r$band)))
  expect_true(is.logical(r$band$open_low) && is.logical(r$band$open_high))
})

test_that("pooled percentile lies in [0, 1] with a documented basis", {
  r <- position_on_surface(0.49, "pabak", pi_hat = 0.5, k = 5, N = 200)
  expect_gte(r$percentile, 0)
  expect_lte(r$percentile, 1)
  expect_true(r$percentile_basis %in%
              c("pooled-achievable-range",
                "pooled-achievable-range-delta-approx",
                "user-supplied-cohort"))
})

test_that("higher PABAK ranks higher on the pooled percentile (monotone)", {
  # The pooled percentile is monotone in obs_value by construction.
  lo <- position_on_surface(0.20, "pabak", pi_hat = 0.5, k = 15, N = 1000)
  hi <- position_on_surface(0.92, "pabak", pi_hat = 0.5, k = 15, N = 1000)
  expect_lt(lo$percentile, hi$percentile)
  # High q_hat sits high in the achievable range; its band is on quality.
  expect_gt(hi$q_hat, 0.90)
  expect_true(is.finite(hi$band$lo) || isTRUE(hi$band$open_low) ||
              is.na(hi$band$lo))
})

test_that("low PABAK at balanced prevalence implies a low q_hat", {
  # PABAK = 0.02 implies q_hat = 0.5 * (1 + sqrt(0.02)) ~ 0.57
  r <- position_on_surface(0.02, "pabak", pi_hat = 0.5, k = 5, N = 200)
  expect_lt(r$q_hat, 0.625)
  # And a low pooled percentile relative to a high-PABAK observation.
  r_hi <- position_on_surface(0.80, "pabak", pi_hat = 0.5, k = 5, N = 200)
  expect_lt(r$percentile, r_hi$percentile)
})

test_that("obs_value outside achievable range is clamped and flagged in notes", {
  # PABAK tops out at 1. Try an impossible value.
  r <- position_on_surface(1.05, "pabak", pi_hat = 0.5, k = 5, N = 200)
  expect_true(any(grepl("clamped", r$notes)))
  expect_equal(r$q_hat, 1.0, tolerance = 1e-3)
})

test_that("caller-supplied per_rep yields a plain cohort percentile", {
  # The legacy per_rep hook is honored for reproducibility audits: it
  # gives a plain cohort percentile (mean(reps <= obs_value)) with basis
  # "user-supplied-cohort" and no sweep/band derived.
  set.seed(7L)
  fake_reps <- stats::rnorm(10000, mean = 0.55, sd = 0.05)
  r_emp <- position_on_surface(
    0.55, "pabak", pi_hat = 0.5, k = 5, N = 500,
    method = "empirical",
    surface_data = list(per_rep = fake_reps)
  )
  expect_equal(r_emp$sampling_method, "empirical")
  expect_equal(r_emp$percentile_basis, "user-supplied-cohort")
  expect_equal(r_emp$percentile, mean(fake_reps <= 0.55), tolerance = 1e-10)
  expect_null(r_emp$sweep)
  expect_null(r_emp$band)
})

test_that("deprecated bands / band_labels warn on non-default and are silent on default", {
  # Non-default bands / band_labels draw a deprecation warning and are
  # ignored (v0.7.1 sweep redesign retired the stipulated q-band partition).
  expect_warning(
    position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200,
                        bands = c(0.5, 0.7, 0.85, 0.95, 1.0)),
    "deprecated"
  )
  expect_warning(
    position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200,
                        band_labels = c("A", "B", "C", "D")),
    "deprecated"
  )
  # Defaults are silent.
  expect_no_warning(
    position_on_surface(0.62, "pabak", pi_hat = 0.5, k = 5, N = 200)
  )
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

test_that("empirical method returns a sensible sweep/band for PABAK=0.62 at k=5,N=200,pi=0.5", {
  r <- position_on_surface(obs_value = 0.62, metric = "pabak",
                           pi_hat = 0.5, k = 5, N = 200,
                           method = "empirical")
  expect_equal(r$sampling_method, "empirical")
  expect_s3_class(r$sweep, "data.frame")
  expect_true(all(r$sweep$p >= 0 & r$sweep$p <= 1))
  expect_gte(r$percentile, 0)
  expect_lte(r$percentile, 1)
  # q_hat for PABAK=0.62 is 0.5*(1+sqrt(0.62)) ~ 0.894; the consistency
  # band on quality should surround it (when finite / not open).
  expect_true(r$q_hat > 0.85 && r$q_hat < 0.94)
  if (is.finite(r$band$lo) && is.finite(r$band$hi) &&
      !isTRUE(r$band$open_low) && !isTRUE(r$band$open_high)) {
    expect_gte(r$q_hat, r$band$lo - 1e-6)
    expect_lte(r$q_hat, r$band$hi + 1e-6)
  }
})

test_that("empirical and delta methods give comparable pooled percentiles", {
  r_e <- position_on_surface(0.62, "pabak", 0.5, 5, 200, method = "empirical")
  r_d <- position_on_surface(0.62, "pabak", 0.5, 5, 200, method = "delta")
  # Both pooled percentiles are finite and land in the same broad region.
  expect_true(is.finite(r_e$percentile) && is.finite(r_d$percentile))
  expect_lt(abs(r_e$percentile - r_d$percentile), 0.25)
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
  fit <- grassr:::fit_tau2_from_ratings(ratings)
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
  fit_tiny <- grassr:::fit_tau2_from_ratings(tiny)
  expect_true(!is.finite(fit_tiny$mu))
  expect_true(grepl("N>=10", fit_tiny$note))
  # Bad input: character rating matrix.
  fit_bad <- grassr:::fit_tau2_from_ratings(list(foo = 1))
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
  # v0.7.1 one-row frame: pooled percentile + consistency-band endpoints,
  # not the retired band-probability / modal-band / confidence columns.
  expect_true(all(c("metric", "observed_value", "pi_hat", "k", "N",
                    "q_hat", "se_q_hat", "percentile", "percentile_basis",
                    "band_lo", "band_hi", "band_open_low", "band_open_high",
                    "band_level", "sampling_method") %in% names(df)))
  expect_false(any(c("modal_band", "modal_band_label", "confidence",
                     "p_band_1") %in% names(df)))
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
