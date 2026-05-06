#' Compute the full inter-rater agreement panel
#'
#' Computes Cohen's kappa, PABAK, Gwet's AC1, positive and negative agreement,
#' Byrt prevalence and bias indices, and two confidence intervals for kappa
#' (Wald, and a logit-transformed Wilson-type interval).
#'
#' @section Deprecated:
#' `grass_compute()` is deprecated in grass 0.2.0. The new headline API is
#' `grass_report(ratings = Y)` returning a `grass_card` object. The full
#' coefficient panel is available via `summary()` or `as.data.frame()`. See
#' `vignette("reporting-card")` and [grass_report()].
#'
#' @param data Input data. Form depends on `format`:
#'   * `"matrix"` — a 2x2 integer count matrix with `R1` rows, `R2` columns.
#'   * `"wide"` — a data.frame with exactly two rater columns (or pass
#'     `rater_cols = c("r1", "r2")` if the data contains extra columns).
#'   * `"long"` — a data.frame with `subject`, `rater`, and `rating` columns
#'     (names configurable via `subject =`, `rater =`, `rating =`). Must
#'     contain exactly two distinct raters.
#'   * `"paired"` — a list of two equal-length rating vectors, a 2-column
#'     matrix, or a 2-column data.frame.
#' @param format One of `"matrix"`, `"wide"`, `"long"`, `"paired"`.
#' @param positive Optional character string naming the level to treat as the
#'   positive (=1) class. Defaults to `"yes"` if present among levels, then
#'   `"1"` / `"true"` / `"positive"` / `"case"`, then first-encountered.
#' @param ... Format-specific arguments: `rater_cols` for `"wide"`;
#'   `subject`, `rater`, `rating` for `"long"`.
#'
#' @return A `grass_metrics` S3 object.
#' @export
#'
#' @examples
#' # 2x2 counts (Cohen 1960 Table 1)
#' tab <- matrix(c(88, 10, 14, 88), nrow = 2,
#'               dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
#' grass_compute(tab, format = "matrix")
#'
#' # Paired vectors
#' r1 <- c(1, 1, 0, 0, 1, 0, 1)
#' r2 <- c(1, 0, 0, 0, 1, 1, 1)
#' grass_compute(list(r1, r2), format = "paired")
grass_compute <- function(data, format = c("wide", "matrix", "long", "paired"),
                          positive = NULL, ...) {
  msg_once(
    "deprecate_grass_compute",
    paste0(
      "`grass_compute()` is deprecated in grass 0.2.0. ",
      "The new headline API is `grass_report(ratings = Y)`. ",
      "See `vignette('reporting-card')` and `?grass_report`."
    )
  )
  format <- match.arg(format)
  call <- match.call()
  norm <- normalize_input(data, format = format, positive = positive, ...)

  if (!is.null(norm$table)) {
    tab <- norm$table
    n_ratings <- sum(tab)
  } else {
    tab <- build_table(norm$r1, norm$r2)
    n_ratings <- length(norm$r1)
  }

  values <- compute_agreement_metrics(tab)

  # Two distinct small-N signals run at different severities:
  #   - Here (N < 10): the metrics themselves are unreliable. Warning.
  #   - In print.grass_result (N < 30): the reference deltas get noisy even
  #     if the metrics are defensible. Informational note only.
  if (n_ratings < 10) {
    warning("Small sample (N = ", n_ratings, "). Inference based on GRASS ",
            "simulations is unreliable below N = 10.",
            call. = FALSE)
  }

  new_grass_metrics(
    values = values,
    n = n_ratings,
    table = tab,
    positive_level = norm$positive_level,
    n_dropped = norm$n_dropped %||% 0L,
    call = call
  )
}
