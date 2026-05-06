#' Report agreement metrics across groups in a tidy frame
#'
#' Splits `data` by `group`, runs [grass_report()] on each subset, and binds
#' the per-group results into a single tidy data.frame using the compact
#' `as.data.frame(..., compact = TRUE)` output (drops the multi-line
#' `regime_note` that wraps poorly in bound frames).
#'
#' @section Deprecated:
#' `grass_report_by()` is deprecated in grass 0.2.0. The new
#' `grass_report(ratings = Y)` headline API takes a single rating matrix and
#' returns a `grass_card`; for cohort-by-cohort analysis, loop over
#' subsets and bind the resulting cards via `do.call(rbind,
#' lapply(parts, function(Y) as.data.frame(grass_report(ratings = Y))))`.
#'
#' @param data A data.frame.
#' @param group Grouping column. Accepts either a bare symbol (e.g.,
#'   `cohort`) or a character string (e.g., `"cohort"`).
#' @param ... Passed through to [grass_report()] (e.g., `spec`, `format`,
#'   `positive`, `id_col`, `rater_cols`, `response`, `prevalence`).
#' @param .cohort_col Name of the column in the returned data.frame that
#'   carries the group value. Default `".cohort"`.
#' @param .parallel If `TRUE`, run per-group reports in parallel via
#'   `future.apply::future_lapply()` using the user's active
#'   `future::plan()`. When enabled, a `progressr::progressor()` reports
#'   progress — wrap the call in `progressr::with_progress({ ... })` to
#'   render a progress bar. Both `future.apply` and `progressr` are
#'   optional dependencies; `grass_report_by()` errors cleanly if either
#'   is missing and names the install command. Default `FALSE` (sequential).
#'
#' @return A data.frame with one row per group. Column `.cohort` carries the
#'   group value as character; remaining columns match
#'   `as.data.frame(result, compact = TRUE)`.
#' @export
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   xray_id = sprintf("CXR-%03d", 1:60),
#'   rater_A = sample(c("abnormal","normal"), 60, replace = TRUE),
#'   rater_B = sample(c("abnormal","normal"), 60, replace = TRUE),
#'   cohort  = rep(c("site1","site2","site3"), each = 20)
#' )
#' grass_report_by(df, cohort, id_col = "xray_id", positive = "abnormal")
#' \dontrun{
#' future::plan(future::multisession, workers = 2)
#' progressr::with_progress({
#'   grass_report_by(df, cohort, id_col = "xray_id",
#'                   positive = "abnormal", .parallel = TRUE)
#' })
#' }
grass_report_by <- function(data, group, ..., .cohort_col = ".cohort",
                            .parallel = FALSE) {
  msg_once(
    "deprecate_grass_report_by",
    paste0(
      "`grass_report_by()` is deprecated in grass 0.2.0. ",
      "The new `grass_report(ratings = Y)` API takes a single rating matrix; ",
      "for cohort splits, loop over subsets and bind the resulting cards. ",
      "See `vignette('reporting-card')` and `?grass_report`."
    )
  )
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  # Accept `group` as bare symbol OR string. substitute() captures the
  # caller's expression without evaluating it — that's the NSE trick.
  group_expr <- substitute(group)
  group_name <- if (is.character(group_expr)) {
    group_expr
  } else if (is.name(group_expr) || is.symbol(group_expr)) {
    as.character(group_expr)
  } else {
    stop("`group` must be a column name (bare or quoted), got: ",
         deparse(group_expr), ".", call. = FALSE)
  }
  if (length(group_name) != 1 || !group_name %in% names(data)) {
    stop("`group` must name a single column in `data`. Got ",
         shQuote(group_name), "; available: ",
         paste(names(data), collapse = ", "), ".", call. = FALSE)
  }

  # Remove the grouping column from the per-subset data so it is not treated
  # as a rater column. Users can still pass `id_col =` via `...` if their
  # subject identifier is a separate column.
  g <- as.character(data[[group_name]])
  data_no_group <- data[, setdiff(names(data), group_name), drop = FALSE]

  parts <- split(data_no_group, g)
  keys <- names(parts)
  # Materialize `...` into a concrete list so it survives the worker
  # boundary under future.apply (multisession serializes the closure but
  # not the calling frame's dots).
  extra_args <- list(...)
  cohort_col <- .cohort_col

  one_row <- function(key) {
    sub <- parts[[key]]
    result <- do.call(grass_report, c(list(sub), extra_args))
    row <- as.data.frame(result, compact = TRUE)
    cbind(setNames(data.frame(key, stringsAsFactors = FALSE), cohort_col),
          row, stringsAsFactors = FALSE)
  }

  if (isTRUE(.parallel)) {
    missing_pkgs <- c("future.apply", "progressr")[
      !vapply(c("future.apply", "progressr"),
              requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs)) {
      stop(".parallel = TRUE needs: ", paste(missing_pkgs, collapse = ", "),
           ". Install with install.packages(c(",
           paste(shQuote(missing_pkgs), collapse = ", "), ")).",
           call. = FALSE)
    }
    p <- progressr::progressor(steps = length(keys))
    rows <- future.apply::future_lapply(keys, function(key) {
      out <- one_row(key)
      p(sprintf("cohort %s", key))
      out
    }, future.seed = TRUE)
  } else {
    rows <- lapply(keys, one_row)
  }
  do.call(rbind, rows)
}
