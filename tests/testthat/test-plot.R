# NOTE (v0.2.0): grass_report() was breaking-changed for the Target-2
# Report Card. plot() on a `grass_result` is retired; `plot.grass_card` is
# implemented separately by Phase 4B. The OLD-API plot test below is
# skipped under v0.2.0. See grass/design/v0.2.0_paper_alignment.md.

source(test_path("fixtures", "published-tables.R"))

test_that("plot methods return ggplot objects when ggplot2 is available", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("plot.grass_metrics errors with a redirect message", {
  skip_if_not_installed("ggplot2")
  m <- grass_compute(fixture_cohen_1960, format = "matrix")
  expect_error(plot(m), "grass_report")
})

test_that("theme_grass returns a ggplot2 theme", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(theme_grass(), "theme")
})
