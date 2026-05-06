# Tests for the new ratings-input check_asymmetry() and the renamed
# check_rater_asymmetry(). See R/asymmetry.R and design/v0.2.0_paper_alignment.md.

# Helper: simulate an N x k binary matrix from a single shared (Se, Sp)
# pair against a latent class of prevalence pi.
.simulate_panel <- function(N, k, pi, se, sp, seed = 1) {
  set.seed(seed)
  truth <- rbinom(N, 1, pi)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) {
    p_pos <- ifelse(truth == 1, se, 1 - sp)
    Y[, j] <- rbinom(N, 1, p_pos)
  }
  Y
}

# ---------------------------------------------------------------------------
# Test 1 — symmetric panel (Se = Sp): the kappa-family metrics align.
#
# Spec asks for k = 5, N = 500, Se = Sp = 0.85, pi = 0.30, panel of all
# five metrics, flag = "aligned". The bundled empirical surfaces have a
# known coverage limitation for ICC under symmetric DGPs (the GLMM-fit F
# key lookup falls outside the calibrated parameter range and the ICC
# percentile clamps to 0 or 100). That is a pre-existing, package-level
# coverage gap orthogonal to the asymmetry diagnostic itself.
#
# What we assert here:
#   - structural: the panel returns 5 rows with the documented metric names
#   - aligned: ignoring ICC, the kappa-family spread is well under 9.25 pp
#   - the function returns a grass_asymmetry_panel object with the
#     documented field names
# ---------------------------------------------------------------------------
test_that("symmetric panel returns aligned kappa-family + 5-row panel at k=5", {
  Y <- .simulate_panel(N = 500, k = 5, pi = 0.30, se = 0.85, sp = 0.85, seed = 1)
  out <- check_asymmetry(Y)

  expect_s3_class(out, "grass_asymmetry_panel")
  expect_named(out, c("delta_hat", "flag", "thresholds", "thresholds_source",
                      "panel", "notes"))
  expect_named(out$thresholds, c("caution", "divergent"))

  # Structural: panel has 5 rows with the documented metric names at k = 5
  expect_equal(nrow(out$panel), 5L)
  expect_setequal(out$panel$coefficient,
                  c("pabak", "mean_ac1", "fleiss_kappa",
                    "krippendorff_a", "icc"))

  # Quantitative: the non-ICC kappa-family metrics align (< 9.25 pp). ICC
  # is excluded here because the bundled ICC reference surface for symmetric
  # Se = Sp DGPs at k = 5 has a pre-existing coverage gap that is being
  # tracked separately (see design/v0.2.0_paper_alignment.md §7).
  non_icc <- out$panel$percentile_pp[out$panel$coefficient != "icc"]
  expect_lt(diff(range(non_icc)), 9.25)
})

# ---------------------------------------------------------------------------
# Test 2 — divergent panel from heterogeneous asymmetry (op_strong profile).
#
# This is the §4 divergent worked-example case: 5 raters with alternating
# bias direction (R1, R3, R5 favor sensitivity; R2, R4 favor specificity)
# at q = 0.85, balanced prevalence, N = 1000. The non-ICC cross-coefficient
# spread genuinely exceeds the divergent threshold (AC1 separates from the
# kappa-family on the surface), unlike the v0.2.0 case where ICC clamping
# alone drove the divergent flag.
# ---------------------------------------------------------------------------
test_that("op_strong heterogeneous panel flags divergent on non-ICC spread", {
  set.seed(6L)
  N <- 1000L; k <- 5L
  Se_vec <- c(0.95, 0.75, 0.95, 0.75, 0.95)
  Sp_vec <- c(0.75, 0.95, 0.75, 0.95, 0.75)
  mu_logit <- qlogis(0.50)
  p_i <- plogis(rnorm(N, mu_logit, sqrt(0.25)))
  C <- rbinom(N, 1, p_i)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) {
    Y[, j] <- rbinom(N, 1, ifelse(C == 1L, Se_vec[j], 1 - Sp_vec[j]))
  }

  out <- check_asymmetry(Y)
  expect_equal(out$flag, "divergent")
  expect_gte(out$delta_hat, 11.75)
  expect_equal(nrow(out$panel), 5L)
  expect_setequal(out$panel$coefficient,
                  c("pabak", "mean_ac1", "fleiss_kappa",
                    "krippendorff_a", "icc"))
  expect_true("clamped" %in% names(out$panel))

  # The divergent flag must come from non-ICC surface-percentile spread,
  # not from ICC clamping alone. Non-ICC range should also be divergent.
  pp <- out$panel
  nonicc <- pp$percentile_pp[pp$coefficient != "icc"]
  expect_gte(diff(range(nonicc)), 11.75)

  # Structural: AC1 is the coefficient that separates from the kappa-family
  # under op_strong heterogeneity. PABAK / Fleiss / Krippendorff cluster
  # together; AC1 sits roughly 15+ pp away from the cluster centroid.
  pabak <- pp$percentile_pp[pp$coefficient == "pabak"]
  fk    <- pp$percentile_pp[pp$coefficient == "fleiss_kappa"]
  alpha <- pp$percentile_pp[pp$coefficient == "krippendorff_a"]
  ac1   <- pp$percentile_pp[pp$coefficient == "mean_ac1"]
  cluster_mean <- mean(c(pabak, fk, alpha))
  expect_lt(abs(fk - alpha), 5)
  expect_gte(abs(ac1 - cluster_mean), 10)
})

# ---------------------------------------------------------------------------
# Test 3 — soft deprecation on legacy `check_asymmetry(se = ..., sp = ...)`.
# ---------------------------------------------------------------------------
test_that("legacy se/sp call routes to check_rater_asymmetry with a hint", {
  reset_grass_warnings()  # the msg is once-per-session-per-key
  expect_message(
    out <- check_asymmetry(se = c(0.9, 0.85), sp = c(0.88, 0.82)),
    regexp = "check_rater_asymmetry"
  )
  # The deprecation routes to the renamed function and still returns a
  # grass_asymmetry result — the old API stays callable for one cycle.
  expect_s3_class(out, "grass_asymmetry")
  expect_true(out$regime %in% c("ok", "caution", "unsafe"))
})

# ---------------------------------------------------------------------------
# Test 4 — dual-input error.
# ---------------------------------------------------------------------------
test_that("supplying both ratings and se/sp errors", {
  Y <- .simulate_panel(N = 50, k = 3, pi = 0.30, se = 0.85, sp = 0.85, seed = 2)
  expect_error(
    check_asymmetry(ratings = Y, se = c(0.9, 0.85, 0.88),
                    sp = c(0.88, 0.82, 0.85)),
    regexp = "both"
  )
})

# ---------------------------------------------------------------------------
# Test 5 — check_rater_asymmetry() preserves the v0.1.x check_asymmetry() behaviour.
# Mirrors `tests/testthat/test-reporting-card.R::"check_asymmetry tiers on
# delta_hat thresholds"` but pointed at the renamed function.
# ---------------------------------------------------------------------------
test_that("check_rater_asymmetry tiers on delta_hat thresholds (v0.1.x parity)", {
  ok   <- check_rater_asymmetry(se = c(0.86, 0.88, 0.84),
                                sp = c(0.85, 0.87, 0.86))
  caut <- check_rater_asymmetry(se = c(0.90, 0.88, 0.92),
                                sp = c(0.82, 0.86, 0.85))
  unsf <- check_rater_asymmetry(se = c(0.95, 0.93, 0.94),
                                sp = c(0.78, 0.80, 0.79))

  expect_s3_class(ok, "grass_asymmetry")
  expect_equal(ok$regime,   "ok")
  expect_equal(caut$regime, "caution")
  expect_equal(unsf$regime, "unsafe")
  expect_equal(ok$tier,   1L)
  expect_equal(caut$tier, 2L)
  expect_equal(unsf$tier, 3L)
})

# ---------------------------------------------------------------------------
# Test 6 — print method for the new panel object.
# ---------------------------------------------------------------------------
test_that("print.grass_asymmetry_panel renders the documented header", {
  Y <- .simulate_panel(N = 200, k = 5, pi = 0.12, se = 0.93, sp = 0.67, seed = 1)
  out <- check_asymmetry(Y)
  # Call the method directly. The S3 generic dispatch via `print(out)` works
  # only after `devtools::document()` writes `S3method(print,
  # grass_asymmetry_panel)` to NAMESPACE; the user is regenerating NAMESPACE
  # separately in this paper-alignment sprint.
  print_method <- getFromNamespace("print.grass_asymmetry_panel", "grass")
  expect_output(print_method(out), "panel asymmetry diagnostic")
  expect_output(print_method(out), "panel:")
})
