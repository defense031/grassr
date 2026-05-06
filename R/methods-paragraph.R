#' Manuscript-ready methods paragraph for a grass result
#'
#' Generates a GRRAS-compliant (Guidelines for Reporting Reliability and
#' Agreement Studies) methods paragraph pre-filled with the study's
#' numbers. The paragraph describes the inter-rater agreement analysis
#' precisely enough to drop into a manuscript methods section.
#'
#' @section Extension contract:
#'
#' `grass_methods()` is an S3 generic dispatching on
#' `result$spec$family`. Each metric family implements its own method
#' (`grass_methods.grass_spec_binary`, and — in future releases —
#' `grass_methods.grass_spec_ordinal`, `.grass_spec_multirater`,
#' `.grass_spec_continuous`). Every method must return a single
#' character string and honour the `format` argument.
#'
#' @section Deprecated:
#' `grass_methods()` is deprecated in grass 0.2.0. The Target-2 reporting
#' workflow does not generate a regime/PI-BI methods paragraph; the
#' published methods text now references `position_on_surface()` and
#' `check_asymmetry()` directly. See `vignette("reporting-card")` and
#' [grass_report()].
#'
#' @param result A `grass_result` object (from [grass_report()]).
#' @param format One of `"markdown"` (default), `"latex"`, or
#'   `"plain"`. Controls emphasis markers and a small amount of
#'   mathematical notation.
#' @param digits Number of decimal places for reported numerics.
#'   Default 2 for main numbers (agreement coefficients), 3 for
#'   derived indices (PI, BI).
#' @param ... Unused. Reserved for family-specific extensions.
#'
#' @return A single character string. Newlines are literal `\n`; the
#'   output can be pasted into an `.Rmd` or `.tex` document directly.
#' @export
#'
#' @examples
#' tab <- matrix(c(88, 10, 14, 88), nrow = 2,
#'               dimnames = list(R1 = c("0","1"), R2 = c("0","1")))
#' r <- grass_report(tab, format = "matrix")
#' cat(grass_methods(r))
grass_methods <- function(result, format = c("markdown", "latex", "plain"),
                          digits = 2, ...) {
  msg_once(
    "deprecate_grass_methods",
    paste0(
      "`grass_methods()` is deprecated in grass 0.2.0. ",
      "The new reporting workflow renders its own methods text via ",
      "`print(grass_report(ratings = Y))`. ",
      "See `vignette('reporting-card')` and `?grass_report`."
    )
  )
  if (!inherits(result, "grass_result")) {
    stop("`result` must be a grass_result, as returned by grass_report().",
         call. = FALSE)
  }
  format <- match.arg(format)
  spec <- result$spec
  if (is.null(spec)) {
    stop("`result$spec` is NULL; this grass_result was produced by a ",
         "pre-0.1.2 version. Re-run grass_report() to use grass_methods().",
         call. = FALSE)
  }
  UseMethod("grass_methods", spec)
}

#' @export
grass_methods.grass_spec_binary <- function(result, format = "markdown",
                                            digits = 2, ...) {
  v <- result$metrics$values
  n <- result$metrics$n
  prev <- result$prevalence
  prev_source <- result$prevalence_source
  regime <- result$regime
  pi_val <- unname(v["prevalence_index"])
  bi_val <- unname(v["bias_index"])
  kappa <- unname(v["kappa"])
  kappa_lo <- unname(v["kappa_wilson_lower"])
  kappa_hi <- unname(v["kappa_wilson_upper"])
  PABAK <- unname(v["PABAK"])
  AC1 <- unname(v["AC1"])
  lvl <- result$spec$reference_level

  emph <- function(s) switch(format,
                             markdown = paste0("*", s, "*"),
                             latex    = paste0("\\emph{", s, "}"),
                             plain    = s)
  kappa_sym <- switch(format,
                      markdown = "\u03ba",
                      latex    = "$\\kappa$",
                      plain    = "kappa")
  # GRRAS citation, format-aware. Used in the opening sentence so a user
  # who pastes grass_methods() output into a manuscript lands with the
  # reporting-guideline citation already in place.
  grras_cite <- switch(format,
                       markdown = "[Kottner et al. 2011](https://doi.org/10.1016/j.jclinepi.2010.03.002)",
                       latex    = "\\citep{kottner2011grras}",
                       plain    = "(Kottner et al. 2011)")

  ref_sentence <- if (is.null(lvl)) {
    "No prevalence-conditioned reference curve was attached for this analysis."
  } else {
    sprintf(
      paste0("Observed agreement was compared against the prevalence-conditioned ",
             "reference curve for raters operating at Se = Sp = %.2f (the analytical ",
             "expected value of each coefficient under conditional independence of ",
             "the two raters given the true class). Youden's J at this reference ",
             "operating point is %.2f. Because Cohen's %s, PABAK, and Gwet's AC1 are ",
             "known to disagree at skewed prevalences under Landis-Koch fixed cutoffs, ",
             "we report the signed distance from this analytical reference rather ",
             "than a single static interpretation band (Semmel 202X)."),
      lvl, 2 * lvl - 1, kappa_sym)
  }

  regime_sentence <- sprintf(
    "The (PI, BI) regime was classified as %s (PI = %.*f, BI = %.*f).",
    emph(regime), digits + 1, pi_val, digits + 1, bi_val)

  prev_sentence <- if (identical(prev_source, "user_supplied")) {
    sprintf("Disease prevalence was user-specified as %.*f.", digits, prev)
  } else {
    sprintf(
      "Disease prevalence was estimated from the two raters' marginal positive rates as %.*f.",
      digits, prev)
  }

  agreement_sentence <- sprintf(
    paste0("Observed agreement was: Cohen's %s = %.*f (95%% CI %.*f, %.*f), ",
           "PABAK = %.*f, Gwet's AC1 = %.*f."),
    kappa_sym, digits, kappa, digits, kappa_lo, digits, kappa_hi,
    digits, PABAK, digits, AC1)

  opening <- sprintf(
    paste0("Inter-rater agreement between two raters on a binary outcome (N = %d ",
           "rating pairs) was analysed under the GRASS framework (Guide for Rater ",
           "Agreement under Structural Skew), with reporting aligned to the ",
           "Guidelines for Reporting Reliability and Agreement Studies %s. ",
           "We computed Cohen's %s, PABAK (prevalence- and bias-adjusted %s), ",
           "and Gwet's AC1 with Wilson-logit 95%% confidence intervals on %s."),
    n, grras_cite, kappa_sym, kappa_sym, kappa_sym)

  paste(opening, prev_sentence, agreement_sentence, regime_sentence,
        ref_sentence, sep = " ")
}

# Stub dispatch for the future families so calling grass_methods() on a
# result built from a stub spec (if somehow constructed) gives a clear
# message rather than a cryptic no-applicable-method error.
#' @export
grass_methods.grass_spec_ordinal    <- function(result, ...) stop_family_unimplemented("ordinal")
#' @export
grass_methods.grass_spec_multirater <- function(result, ...) stop_family_unimplemented("multirater")
#' @export
grass_methods.grass_spec_continuous <- function(result, ...) stop_family_unimplemented("continuous")
