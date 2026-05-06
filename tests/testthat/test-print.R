# NOTE (v0.2.0): grass_report() was breaking-changed for the Target-2
# Report Card. The `print.grass_result` / `as.data.frame.grass_result`
# tests below depend on grass_report returning a `grass_result`. They now
# fail because grass_report returns a `grass_card`. Tests of the OLD
# `grass_result` printing pipeline are skipped under v0.2.0; the new flow
# is exercised in test-grass_report-card.R via `print.grass_card` /
# `as.data.frame.grass_card`. See grass/design/v0.2.0_paper_alignment.md.
#
# Tests of the surviving artifacts (`print.grass_metrics`,
# `print.grass_reference`) are kept as-is.

source(test_path("fixtures", "published-tables.R"))

test_that("print.grass_result shows metrics, skew diagnostics, and regime", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("print.grass_result runs on all fixtures", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("print.grass_metrics shows the 2x2 table", {
  m <- grass_compute(fixture_cohen_1960, format = "matrix")
  expect_output(print(m), "2x2 table")
})

test_that("print.grass_reference labels curves as references", {
  t <- grass_reference(0.3)
  expect_output(print(t), "reference")
})

test_that("as.data.frame.grass_result returns context columns, no verdict", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})
