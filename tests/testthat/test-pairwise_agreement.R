# Tests for v0.2.3 pairwise_agreement().
#
# Covers input validation, output structure, the k = 2 omit-pool branch,
# pooled-reference recovery on the §4.2 op_strong worked example, and the
# print method.

# Shared fixture mirrors test-grass_report-card.R's .simulate_op_strong_panel
# so the two test files cite the same DGP.
.pairwise_op_strong_panel <- function(seed = 17L, N = 1000L) {
  set.seed(seed)
  k <- 5L
  Se <- c(0.95, 0.75, 0.95, 0.75, 0.95)
  Sp <- c(0.75, 0.95, 0.75, 0.95, 0.75)
  mu_logit <- qlogis(0.50)
  p_i <- plogis(rnorm(N, mu_logit, sqrt(0.25)))
  C <- rbinom(N, 1, p_i)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) {
    Y[, j] <- rbinom(N, 1, ifelse(C == 1L, Se[j], 1 - Sp[j]))
  }
  Y
}

# ---- Test 1 ---------------------------------------------------------------
# Input validation.

test_that("pairwise_agreement errors when k < 2", {
  Y <- matrix(rbinom(50, 1, 0.5), nrow = 50, ncol = 1)
  expect_error(pairwise_agreement(Y),
               regexp = "at least 2 rater columns",
               fixed = FALSE)
})

# ---- Test 2 ---------------------------------------------------------------
# Returned object structure: class, fields, sample list keys.

test_that("pairwise_agreement returns a grass_pairwise S3 with documented fields", {
  Y <- .pairwise_op_strong_panel(seed = 17L, N = 200L)
  pw <- pairwise_agreement(Y)

  expect_s3_class(pw, "grass_pairwise")
  expect_type(pw, "list")

  fields <- c("pabak_matrix", "percentile_matrix", "marginal_matrix",
              "band_matrix", "qualifier_matrix", "pooled_per_rater",
              "sample", "notes", "call")
  for (f in fields) expect_true(f %in% names(pw),
                                info = paste("missing field:", f))

  expect_true(all(c("k", "N", "pi_hat", "tau2_hat", "axis") %in%
                  names(pw$sample)))
  expect_equal(pw$sample$k, 5L)
  expect_equal(pw$sample$N, 200L)
  expect_true(is.finite(pw$sample$pi_hat))
  expect_true(is.finite(pw$sample$tau2_hat))
})

# ---- Test 3 ---------------------------------------------------------------
# Matrix shape + symmetry.

test_that("pabak_matrix and percentile_matrix are k x k and symmetric", {
  Y <- .pairwise_op_strong_panel(seed = 17L, N = 200L)
  pw <- pairwise_agreement(Y)
  k <- ncol(Y)

  expect_equal(dim(pw$pabak_matrix), c(k, k))
  expect_equal(dim(pw$percentile_matrix), c(k, k))
  expect_equal(dim(pw$marginal_matrix), c(k, k))

  # PABAK matrix: diagonal = 1, symmetric off-diagonal.
  expect_equal(unname(diag(pw$pabak_matrix)), rep(1, k))
  expect_equal(pw$pabak_matrix, t(pw$pabak_matrix))

  # Percentile matrix: diagonal NA, symmetric.
  expect_true(all(is.na(diag(pw$percentile_matrix))))
  expect_equal(pw$percentile_matrix[lower.tri(pw$percentile_matrix)],
               t(pw$percentile_matrix)[lower.tri(pw$percentile_matrix)])

  # Off-diagonal percentiles in [0, 100].
  off <- pw$percentile_matrix[!is.na(pw$percentile_matrix)]
  expect_true(all(off >= 0 & off <= 100))
})

# ---- Test 4 ---------------------------------------------------------------
# Pooled per-rater table at k >= 3.

test_that("pooled_per_rater has k rows of finite Se_tilde/Sp_tilde at k >= 3", {
  Y <- .pairwise_op_strong_panel(seed = 17L, N = 200L)
  pw <- pairwise_agreement(Y)

  pp <- pw$pooled_per_rater
  expect_s3_class(pp, "data.frame")
  expect_equal(nrow(pp), 5L)
  expect_true(all(c("rater", "se_tilde", "sp_tilde",
                    "n_pool_pos", "n_pool_neg",
                    "n_pool_excluded") %in% names(pp)))
  expect_true(all(is.finite(pp$se_tilde)))
  expect_true(all(is.finite(pp$sp_tilde)))
})

# ---- Test 5 ---------------------------------------------------------------
# Pooled-reference recovers per-rater bias direction on op_strong DGP.
# §3.3 / §4.2: Se-favoring raters (R1, R3, R5) should have Se_tilde > Sp_tilde;
# Sp-favoring raters (R2, R4) should have Sp_tilde > Se_tilde. This is the
# load-bearing methodological claim of the divergent recovery branch.

test_that("pooled-reference recovers op_strong per-rater bias direction", {
  Y <- .pairwise_op_strong_panel(seed = 17L, N = 1000L)
  pw <- pairwise_agreement(Y)
  pp <- pw$pooled_per_rater

  for (j in c(1L, 3L, 5L)) {
    expect_gt(pp$se_tilde[j], pp$sp_tilde[j],
              label = sprintf("R%d should be Se-favoring under pooled ref", j))
  }
  for (j in c(2L, 4L)) {
    expect_gt(pp$sp_tilde[j], pp$se_tilde[j],
              label = sprintf("R%d should be Sp-favoring under pooled ref", j))
  }
})

# ---- Test 6 ---------------------------------------------------------------
# k = 2 omits pooled per-rater table and emits a note explaining why.

test_that("pairwise_agreement omits pooled per-rater table at k = 2", {
  set.seed(11)
  Y <- matrix(rbinom(400, 1, 0.4), nrow = 200, ncol = 2)
  pw <- pairwise_agreement(Y)

  expect_equal(pw$sample$k, 2L)
  expect_true(all(is.na(pw$pooled_per_rater$se_tilde)))
  expect_true(all(is.na(pw$pooled_per_rater$sp_tilde)))
  expect_true(any(grepl("k < 3", pw$notes)))
})

# ---- Test 7 ---------------------------------------------------------------
# Print method runs and surfaces the expected headers.

test_that("print.grass_pairwise renders header and pooled per-rater table", {
  Y <- .pairwise_op_strong_panel(seed = 17L, N = 200L)
  pw <- pairwise_agreement(Y)

  out <- capture.output(print(pw))
  expect_true(any(grepl("GRASS Pairwise Reliability", out, fixed = TRUE)))
  expect_true(any(grepl("Pairwise PABAK", out, fixed = TRUE)))
  expect_true(any(grepl("Per-rater behavior", out, fixed = TRUE)))
})
