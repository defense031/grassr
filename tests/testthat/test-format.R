# NOTE (v0.2.0): grass_report() was breaking-changed for the Target-2
# Report Card. The OLD `grass_format_report()` consumes a `grass_result`
# (returned by the OLD grass_report API). Those calls now fail because
# grass_report returns a `grass_card`. Tests of the OLD `grass_format_report`
# pipeline are skipped under v0.2.0; the new flow is exercised in
# test-grass_report-card.R via `print.grass_card` / `format.grass_card`.
# See grass/design/v0.2.0_paper_alignment.md.

source(test_path("fixtures", "published-tables.R"))

test_that("grass_format_report returns a single character string", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("grass_format_report includes all expected components", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("ascii = TRUE emits 'kappa' instead of Unicode", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("ascii = FALSE emits Unicode kappa", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("grass_format_report errors on non-grass_result input", {
  expect_error(grass_format_report(list(a = 1)), "grass_result")
})
