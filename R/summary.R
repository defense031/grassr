#' @export
summary.grass_metrics <- function(object, ...) {
  structure(object, class = c("summary.grass_metrics", class(object)))
}

#' @export
print.summary.grass_metrics <- function(x, digits = 4, ...) {
  v <- x$values
  cat("grass metrics -- summary\n")
  cat("  N = ", x$n, "   positive level = ", shQuote(x$positive_level), "\n", sep = "")
  cat("  2x2 table\n")
  print(x$table)
  cat("\n")
  cat("  Agreement\n")
  cat("    P0 (observed)      : ", fmt_num(v["P0"], digits), "\n", sep = "")
  cat("    Pe (expected)      : ", fmt_num(v["Pe"], digits), "\n", sep = "")
  cat("    Cohen's kappa      : ", fmt_num(v["kappa"], digits), "\n", sep = "")
  cat("      Wald 95% CI      : ",
      fmt_ci(v["kappa_wald_lower"], v["kappa_wald_upper"], digits), "\n", sep = "")
  cat("      Wilson-logit 95% : ",
      fmt_ci(v["kappa_wilson_lower"], v["kappa_wilson_upper"], digits), "\n", sep = "")
  cat("    PABAK              : ", fmt_num(v["PABAK"], digits), "\n", sep = "")
  cat("    Gwet's AC1         : ", fmt_num(v["AC1"], digits), "\n", sep = "")
  cat("    Positive agreement : ", fmt_num(v["pos_agreement"], digits), "\n", sep = "")
  cat("    Negative agreement : ", fmt_num(v["neg_agreement"], digits), "\n", sep = "")
  cat("\n  Skew diagnostics\n")
  cat("    Prevalence index   : ", fmt_num(v["prevalence_index"], digits), "\n", sep = "")
  cat("    Bias index         : ", fmt_num(v["bias_index"], digits), "\n", sep = "")
  invisible(x)
}

#' @export
summary.grass_result <- function(object, ...) {
  structure(object, class = c("summary.grass_result", class(object)))
}

#' @export
print.summary.grass_result <- function(x, digits = 4, ...) {
  print.grass_result(x, digits = digits)
  invisible(x)
}

# --------------------------------------------------------------------------
# summary.grass_card -- v0.2.0 Target-2 Report Card detailed summary
# --------------------------------------------------------------------------
# Returns a summary.grass_card list with the full panel data frame, the
# per-rater table (when non-NULL), all notes, sample info, and the full
# delta. print.summary.grass_card formats this as a multi-section block.

#' @export
summary.grass_card <- function(object, ...) {
  out <- list(
    sample      = object$sample,
    coefficient = object$coefficient,
    delta       = object$delta,
    panel       = object$panel,
    per_rater   = object$per_rater,
    notes       = object$notes,
    grass_version = object$grass_version,
    timestamp   = object$timestamp
  )
  class(out) <- c("summary.grass_card", "list")
  out
}

#' @export
print.summary.grass_card <- function(x, digits = 3, ...) {
  s <- x$sample
  cat("GRASS Report Card -- summary\n\n")
  cat(sprintf("  sample       : k = %d raters, N = %d, pi_hat = %.*f, axis = %s\n",
              s$k, s$N, digits, s$pi_hat, s$axis))
  if (!is.null(s$tau2_hat) && is.finite(s$tau2_hat)) {
    cat(sprintf("  tau2_hat     : %.*f\n", digits, s$tau2_hat))
  }

  cat("\n  primary coefficient\n")
  co <- x$coefficient
  cat(sprintf("    %-12s : %s\n", "name", co$primary))
  cat(sprintf("    %-12s : %.*f\n", "observed",
              digits, co$observed_value))
  cat(sprintf("    %-12s : %.*f pp\n", "percentile",
              max(digits - 1L, 1L), co$surface_percentile))
  band_str <- if (is.null(co$band) || is.na(co$band)) "NA" else co$band
  qual_str <- if (is.null(co$qualifier) || is.na(co$qualifier)) "NA" else co$qualifier
  cat(sprintf("    %-12s : %s\n", "band", band_str))
  cat(sprintf("    %-12s : %s\n", "qualifier", qual_str))

  cat("\n  delta (cross-coefficient asymmetry)\n")
  cat(sprintf("    %-12s : %.*f pp\n", "delta_hat",
              max(digits - 1L, 1L), x$delta$delta_hat))
  cat(sprintf("    %-12s : %s\n", "flag", x$delta$flag))
  cat(sprintf("    %-12s : caution = %.*f, divergent = %.*f\n",
              "thresholds",
              max(digits - 1L, 1L), unname(x$delta$thresholds[["caution"]]),
              max(digits - 1L, 1L), unname(x$delta$thresholds[["divergent"]])))

  cat("\n  panel (full table)\n")
  pn <- x$panel
  for (i in seq_len(nrow(pn))) {
    cat(sprintf("    %-15s observed = %.*f  pct = %5.*f pp  band = %-12s qualifier = %s\n",
                pn$coefficient[i],
                digits, pn$observed_value[i],
                max(digits - 1L, 1L), pn$surface_percentile[i],
                if (is.na(pn$band[i])) "NA" else pn$band[i],
                if (is.na(pn$qualifier[i])) "NA" else pn$qualifier[i]))
  }

  if (!is.null(x$per_rater) && nrow(x$per_rater) > 0L) {
    cat("\n  per-rater (latent-class fit)\n")
    pr <- x$per_rater
    for (i in seq_len(nrow(pr))) {
      if (isTRUE(pr$bound_only[i])) {
        cat(sprintf("    %-4s  Se in [%.*f, %.*f]   Sp in [%.*f, %.*f]   (Hui-Walter bounds)\n",
                    pr$rater[i],
                    digits, pr$se_lower[i],
                    digits, pr$se_upper[i],
                    digits, pr$sp_lower[i],
                    digits, pr$sp_upper[i]))
      } else {
        cat(sprintf("    %-4s  Se = %.*f  (%.*f, %.*f)   Sp = %.*f  (%.*f, %.*f)\n",
                    pr$rater[i],
                    digits, pr$se_hat[i],
                    digits, pr$se_lower[i], digits, pr$se_upper[i],
                    digits, pr$sp_hat[i],
                    digits, pr$sp_lower[i], digits, pr$sp_upper[i]))
      }
    }
  }

  if (length(x$notes) > 0L) {
    cat("\n  notes\n")
    for (n in x$notes) cat("    - ", n, "\n", sep = "")
  }

  cat(sprintf("\n  grass version : %s\n", as.character(x$grass_version)))
  cat(sprintf("  timestamp     : %s\n", format(x$timestamp)))
  invisible(x)
}

#' @export
summary.grass_reference <- function(object, ...) {
  structure(object, class = c("summary.grass_reference", class(object)))
}

#' @export
print.summary.grass_reference <- function(x, digits = 4, ...) {
  print.grass_reference(x, digits = digits)
  invisible(x)
}
