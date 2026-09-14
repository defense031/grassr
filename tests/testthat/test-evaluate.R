# NOTE (v0.2.0): grass_report() was breaking-changed for the Target-2
# Report Card. The OLD signature `grass_report(data, format = "matrix",
# spec, reference, prevalence)` and its `grass_result` regime / distance /
# reference-level outputs are retired. Tests for that pre-Target-2 API
# are skipped under v0.2.0; the new flow is exercised in
# `test-grass_report-card.R`. See grass/design/v0.2.0_paper_alignment.md.
#
# This file retains the substrate-level checks that still apply:
# and skips the OLD-framework regime / reference / distance assertions.

source(test_path("fixtures", "published-tables.R"))

test_that("grass_report returns a grass_result with expected slots", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("grass_report has no verdict column anywhere", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("Regime classification matches the PI / BI structure", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("reference = 'none' drops the reference curve and distance", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("User-supplied prevalence overrides the marginal estimate", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("Distance column is a signed numeric, not a category", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})
