#' Look up prevalence-conditioned reference values for agreement metrics
#'
#' Returns the prevalence-conditioned reference values for Cohen's kappa,
#' PABAK, and Gwet's AC1 at a given prevalence and a chosen rater-quality
#' band. The reference is the analytical Youden-J-optimal metric value on
#' the Se = Sp diagonal at the chosen `reference_level`, under conditional
#' independence of the two raters given the true class.
#'
#' In the paper these are called the *prevalence-conditioned thresholds*.
#' The package uses the name `reference` because its API reports the
#' values alongside prevalence index and bias index rather than applying
#' them as cutoffs.
#'
#' @section Deprecated:
#' `grass_reference()` is deprecated in grass 0.2.0. The percentile machinery
#' in [position_on_surface()] and the new [grass_report()] supersede the
#' fixed reference-curve lookup. See `vignette("reporting-card")` for the
#' v0.2.0 surface-position workflow.
#'
#' @param prevalence Numeric in `[0, 1]`. Values outside `[0.01, 0.99]` are
#'   clamped with a warning.
#' @param reference_level One of `0.70`, `0.80`, `0.85`, `0.90`. Default
#'   `0.85`.
#' @param quality Deprecated alias. `"high"` → `reference_level = 0.85`,
#'   `"medium"` → `0.70`. Emits a one-time soft deprecation warning.
#'
#' @return A `grass_reference` S3 object containing a 3-row data.frame
#'   (one row per metric) with columns `metric`, `reference`, `J`,
#'   `reference_level`, `quality`, `prevalence`. The `quality` column
#'   retains the legacy `"high"`/`"medium"` label for 0.85 and 0.70 bands
#'   and is `NA_character_` for the new 0.80 / 0.90 bands.
#' @export
#'
#' @examples
#' grass_reference(0.5)
#' grass_reference(0.1, reference_level = 0.70)
grass_reference <- function(prevalence,
                            reference_level = 0.85,
                            quality = NULL) {
  msg_once(
    "deprecate_grass_reference",
    paste0(
      "`grass_reference()` is deprecated in grass 0.2.0. ",
      "The new headline API is `grass_report(ratings = Y)`; ",
      "for percentile lookups use `position_on_surface(ratings = Y, metric = ...)`. ",
      "See `vignette('reporting-card')` and `?grass_report`."
    )
  )
  call <- match.call()
  if (!is.null(quality)) {
    msg_once("deprecate_quality",
             "grass_reference(): `quality =` is deprecated. Use `reference_level =` (numeric band).")
    reference_level <- switch(quality,
                              high   = 0.85,
                              medium = 0.70,
                              stop("`quality` must be \"high\" or \"medium\".",
                                   call. = FALSE))
  }
  level <- validate_reference_level(reference_level)
  p <- validate_prevalence(prevalence)
  reference_for_binary(p, level, call = call)
}

# Internal: build a grass_reference for the binary family at one
# prevalence and one level.
reference_for_binary <- function(p, level, call = NULL) {
  tbl <- reference_lookup_binary(p, level)
  new_grass_reference(prevalence = p,
                      quality = level_to_quality_label(level),
                      reference = tbl,
                      call = call)
}

# Internal: valid reference levels and their legacy string labels.
.valid_reference_levels <- c(0.70, 0.80, 0.85, 0.90)

validate_reference_level <- function(x) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) ||
      !any(abs(x - .valid_reference_levels) < 1e-8)) {
    stop("`reference_level` must be one of ",
         paste(format(.valid_reference_levels, nsmall = 2), collapse = ", "),
         ". Got: ", deparse(x), ".", call. = FALSE)
  }
  .valid_reference_levels[which.min(abs(x - .valid_reference_levels))]
}

# Map a numeric level to the legacy string label ("high"/"medium") or NA
# for the new bands. Preserves backward compatibility in printed output
# and in the $reference$quality slot.
level_to_quality_label <- function(level) {
  if (abs(level - 0.85) < 1e-8) return("high")
  if (abs(level - 0.70) < 1e-8) return("medium")
  NA_character_
}

# Internal: interpolate the reference_binary long table at a single
# prevalence and one quality level. Returns a plain data.frame with
# columns metric, reference, J, reference_level, quality, prevalence.
reference_lookup_binary <- function(p, level) {
  rb <- reference_binary  # internal dataset in R/sysdata.rda
  sub <- rb[abs(rb$reference_level - level) < 1e-8, , drop = FALSE]
  metrics <- c("kappa", "PABAK", "AC1")
  out <- data.frame(
    metric = metrics,
    reference = vapply(metrics, function(m) {
      rows <- sub[sub$metric == m, , drop = FALSE]
      as.numeric(approx(rows$prevalence, rows$reference, xout = p, rule = 2)$y)
    }, numeric(1)),
    J = vapply(metrics, function(m) {
      rows <- sub[sub$metric == m, , drop = FALSE]
      as.numeric(approx(rows$prevalence, rows$J, xout = p, rule = 2)$y)
    }, numeric(1)),
    reference_level = level,
    quality = level_to_quality_label(level),
    prevalence = p,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out
}

#' Return the full internal reference table
#'
#' For inspection, custom interpolation, or reproducibility. In the paper
#' this is the prevalence-conditioned threshold table.
#'
#' @param reference_level Optional numeric band (`0.70`, `0.80`, `0.85`,
#'   or `0.90`). If supplied, only that band is returned. Default `NULL`
#'   returns all four bands in long form.
#'
#' @return A long-form data.frame with columns `prevalence`,
#'   `reference_level`, `metric`, `reference`, and `J`.
#' @export
#'
#' @examples
#' head(grass_reference_table())
#' head(grass_reference_table(reference_level = 0.85))
grass_reference_table <- function(reference_level = NULL) {
  rb <- as.data.frame(reference_binary)
  if (is.null(reference_level)) return(rb)
  level <- validate_reference_level(reference_level)
  rb[abs(rb$reference_level - level) < 1e-8, , drop = FALSE]
}
