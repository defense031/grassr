# NOTE (v0.2.0): grass_report() was breaking-changed for the Target-2
# Report Card. Tests calling the OLD `grass_report(data, format = "matrix",
# spec, reference)` API or consuming the `grass_result$spec` /
# `grass_result$reference` slots fail under v0.2.0 and are skipped. See
# grass/design/v0.2.0_paper_alignment.md.

source(test_path("fixtures", "published-tables.R"))

# ---- spec-dispatch architecture ---------------------------------------

test_that("grass_spec_binary constructs with valid reference_level", {
  s <- grass_spec_binary()
  expect_s3_class(s, "grass_spec_binary")
  expect_s3_class(s, "grass_spec")
  expect_equal(s$family, "binary")
  expect_equal(s$reference_level, 0.85)
  for (lvl in c(0.70, 0.80, 0.85, 0.90)) {
    expect_equal(grass_spec_binary(reference_level = lvl)$reference_level, lvl)
  }
  expect_null(grass_spec_binary(reference_level = NULL)$reference_level)
})

test_that("grass_spec_binary rejects invalid reference_level", {
  expect_error(grass_spec_binary(reference_level = 0.75), "must be one of")
  expect_error(grass_spec_binary(reference_level = "high"), "must be one of")
  expect_error(grass_spec_binary(reference_level = c(0.85, 0.90)), "must be one of")
})

test_that("stub spec constructors return placeholder specs", {
  for (ctor in list(grass_spec_ordinal, grass_spec_multirater, grass_spec_continuous)) {
    s <- ctor()
    expect_s3_class(s, "grass_spec")
  }
})

test_that("stub specs error at dispatch with a ?grass_roadmap pointer", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("grass_report returns spec slot on result", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("reference_level = NULL skips the reference attachment", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

# ---- reference_level bands --------------------------------------------

test_that("Each reference_level matches the analytical diagonal closed-form", {
  # At p=0.5, Se=Sp=q gives P_agree = q^2 + (1-q)^2 and Pe_cohen = 0.5.
  # kappa = 2*(q^2 + (1-q)^2 - 0.5) = 2q^2 - 2q + ... simplifies to
  # (q-0.5)^2 * 4 basically; cross-check numerically.
  expected <- function(p, q) {
    m <- q*p + (1-q)*(1-p)
    P_agree <- p*(q^2 + (1-q)^2) + (1-p)*(q^2 + (1-q)^2)
    Pe <- m^2 + (1-m)^2
    (P_agree - Pe) / (1 - Pe)
  }
  for (q in c(0.70, 0.80, 0.85, 0.90)) {
    r <- grass_reference(0.5, reference_level = q)
    k <- r$reference$reference[r$reference$metric == "kappa"]
    expect_equal(k, expected(0.5, q), tolerance = 1e-6,
                 info = sprintf("kappa at q=%.2f", q))
  }
})

# ---- .parallel on grass_report_by -------------------------------------

test_that(".parallel matches sequential output", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that(".parallel errors cleanly when future.apply is unavailable", {
  # Simulate missing future.apply via namespace unloading check.
  has_future <- requireNamespace("future.apply", quietly = TRUE)
  has_progressr <- requireNamespace("progressr", quietly = TRUE)
  skip_if(has_future && has_progressr,
          "test requires at least one of future.apply / progressr to be absent")
  # When at least one of them is absent, .parallel = TRUE must error.
  set.seed(12)
  df <- data.frame(
    rA = rbinom(20, 1, 0.3),
    rB = rbinom(20, 1, 0.3),
    cohort = rep(c("s1","s2"), each = 10)
  )
  expect_error(grass_report_by(df, cohort, .parallel = TRUE), "install.packages")
})

# ---- grass_methods() --------------------------------------------------

test_that("grass_methods returns a non-trivial paragraph", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("grass_methods honours all three formats", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("grass_methods describes reference_level = NULL gracefully", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

# ---- legacy reference = ... still works with deprecation --------------

test_that("legacy reference = 'high' maps to reference_level = 0.85", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})
