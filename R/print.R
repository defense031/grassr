#' @export
print.grass_metrics <- function(x, digits = 4, ...) {
  v <- x$values
  cat("grass metrics (N = ", x$n, "; positive level = ", shQuote(x$positive_level), ")\n",
      sep = "")
  if (isTRUE(x$n_dropped > 0)) {
    cat("  (", x$n_dropped, " rating pairs dropped for NA)\n", sep = "")
  }
  cat("  2x2 table\n")
  print(x$table)
  cat("\n")
  cat("  Observed agreement P0  : ", fmt_num(v["P0"], digits), "\n", sep = "")
  cat("  Expected agreement Pe  : ", fmt_num(v["Pe"], digits), "\n", sep = "")
  cat("  Cohen's kappa          : ", fmt_num(v["kappa"], digits),
      "   Wald 95% CI ", fmt_ci(v["kappa_wald_lower"], v["kappa_wald_upper"], digits), "\n",
      sep = "")
  cat("  PABAK                  : ", fmt_num(v["PABAK"], digits), "\n", sep = "")
  cat("  Gwet's AC1             : ", fmt_num(v["AC1"], digits), "\n", sep = "")
  cat("  Positive agreement     : ", fmt_num(v["pos_agreement"], digits), "\n", sep = "")
  cat("  Negative agreement     : ", fmt_num(v["neg_agreement"], digits), "\n", sep = "")
  cat("  Prevalence index (PI)  : ", fmt_num(v["prevalence_index"], digits), "\n", sep = "")
  cat("  Bias index (BI)        : ", fmt_num(v["bias_index"], digits), "\n", sep = "")
  invisible(x)
}

#' @export
print.grass_reference <- function(x, digits = 4, ...) {
  p <- x$prevalence
  ref <- x$reference
  lvl <- ref$reference_level[1]
  cat("grass reference curve at prevalence = ", fmt_num(p, digits),
      " (Se = Sp = ", format(lvl, nsmall = 2), ")\n", sep = "")
  for (i in seq_len(nrow(ref))) {
    cat(sprintf("  %-8s  reference = %s    (Youden's J = %s)\n",
                ref$metric[i],
                fmt_num(ref$reference[i], digits),
                fmt_num(ref$J[i], digits)))
  }
  invisible(x)
}

# --------------------------------------------------------------------------
# print.grass_card -- v0.2.0 Target-2 Report Card
# --------------------------------------------------------------------------
# Renders the four-field card matching paper Figure 2 right panel
# character-for-character. When `flag` is "aligned" or "caution", shows only
# the primary coefficient line + delta + band/qualifier. When "divergent",
# shows ALL panel coefficients with their percentiles + delta(divergent) +
# band=suppressed + per-rater table (if non-NULL).

#' @export
print.grass_card <- function(x, digits = 2, ...) {
  lines <- format(x, digits = digits, ...)
  cat(paste(lines, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
