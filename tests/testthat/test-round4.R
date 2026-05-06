# NOTE (v0.2.0): grass_report() was breaking-changed for the Target-2
# Report Card. Tests calling the OLD `grass_report(data, format = "matrix")`
# API or consuming the OLD `grass_result` (e.g., `grass_format_report`,
# `tidy.grass_result`, `grass_report_by`) now fail and are skipped under
# v0.2.0. The new flow is exercised in test-grass_report-card.R. See
# grass/design/v0.2.0_paper_alignment.md.

source(test_path("fixtures", "published-tables.R"))

# ---- response = alias ------------------------------------------------

test_that("`response =` is an alias for `rater_cols =` in wide format", {
  set.seed(1)
  df <- data.frame(
    xray = 1:20,
    r1   = rbinom(20, 1, 0.3),
    r2   = rbinom(20, 1, 0.3)
  )
  a <- suppressWarnings(
    grass_compute(df, format = "wide", id_col = "xray",
                  rater_cols = c("r1", "r2"))
  )
  b <- suppressWarnings(
    grass_compute(df, format = "wide", id_col = "xray",
                  response = c("r1", "r2"))
  )
  expect_equal(a$values, b$values)
})

test_that("`response =` and `rater_cols =` conflict errors", {
  set.seed(2)
  df <- data.frame(
    r1 = rbinom(20, 1, 0.3),
    r2 = rbinom(20, 1, 0.3),
    other = rbinom(20, 1, 0.3)
  )
  expect_error(
    grass_compute(df, format = "wide",
                  rater_cols = c("r1", "r2"),
                  response = c("r1", "other")),
    "disagree"
  )
})

test_that("`response =` works in long format as alias for `rating =`", {
  long_df <- data.frame(
    subject = rep(1:10, 2),
    rater   = rep(c("r1", "r2"), each = 10),
    score   = rbinom(20, 1, 0.3)
  )
  a <- suppressWarnings(grass_compute(long_df, format = "long", rating = "score"))
  b <- suppressWarnings(grass_compute(long_df, format = "long", response = "score"))
  expect_equal(a$values, b$values)
})

# ---- ci_width on grass_format_report ---------------------------------

test_that("ci_width = TRUE appends width and descriptor", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("ci_width = FALSE (default) does not append width", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

# ---- tidy.grass_result long form -------------------------------------

test_that("tidy.grass_result returns 9 rows with expected schema", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

# ---- grass_report_by -------------------------------------------------

test_that("grass_report_by returns one row per group with .cohort column", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

test_that("grass_report_by accepts a string for `group`", {
  skip("v0.2.0: old framework retired; see grass/design/v0.2.0_paper_alignment.md")
})

# ---- id_col hint in wide-format error --------------------------------

test_that("wide-format error softly suggests id_col when one column is non-binary", {
  df <- data.frame(
    subject_id = sprintf("S-%03d", 1:10),
    a = c(1,0,1,0,1,1,0,1,0,1),
    b = c(0,1,1,0,1,0,1,1,0,0)
  )
  msg <- tryCatch(grass_compute(df, format = "wide"),
                  error = function(e) conditionMessage(e))
  expect_true(grepl("subject_id", msg))
  expect_true(grepl("looks like an identifier", msg))
})

test_that("wide-format error hint is silent when ambiguous (no clear id column)", {
  df <- data.frame(
    a = c(1,0,1,0,1),
    b = c(0,1,1,0,0),
    c = c(1,1,0,0,1)
  )
  msg <- tryCatch(grass_compute(df, format = "wide"),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("looks like an identifier", msg))
})

# ---- reset_grass_warnings --------------------------------------------

test_that("reset_grass_warnings clears the once-per-session cache", {
  env <- grass:::.grass_env
  env$msg_seen[["test_key"]] <- TRUE
  expect_true(isTRUE(env$msg_seen[["test_key"]]))
  reset_grass_warnings()
  expect_null(env$msg_seen[["test_key"]])
})
