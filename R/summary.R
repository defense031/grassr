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

# --------------------------------------------------------------------------
# summary.grass_card -- v0.2.0 Target-2 Report Card detailed summary
# --------------------------------------------------------------------------
# Returns a summary.grass_card list with the full panel data frame, the
# per-rater table (when non-NULL), all notes, sample info, and the full
# delta. print.summary.grass_card formats this as a multi-section block.

#' Summarize a GRASS Report Card
#'
#' The full detail behind the card: the primary coefficient with its
#' estimate `q_hat`, the `delta_hat` diagnostic and its matched-null
#' percentile, the whole panel table, the per-rater latent-class table on
#' a divergent card, and every note the card carries, including the
#' provenance notes the card print leaves out.
#'
#' @param object A `grass_card` object from [grass_report()].
#' @param x A `summary.grass_card` object.
#' @param digits Decimals for printed values. Default 3.
#' @param ... Ignored.
#' @return For `summary()`, a `summary.grass_card` list; its print method
#'   returns it invisibly.
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

#' @rdname summary.grass_card
#' @export
print.summary.grass_card <- function(x, digits = 3, ...) {
  s <- x$sample
  cat("GRASS Report Card -- summary\n\n")
  cat(sprintf("  sample       : k = %d raters, N = %d, pi_hat = %.*f, axis = %s\n",
              s$k, s$N, digits, s$pi_hat, s$axis))

  cat("\n  primary coefficient\n")
  co <- x$coefficient
  cat(sprintf("    %-12s : %s\n", "name", co$primary))
  cat(sprintf("    %-12s : %.*f\n", "observed",
              digits, co$observed_value))
  cat(sprintf("    %-12s : %.*f (of the achievable range in this study context)\n",
              "percentile", max(digits - 1L, 1L), co$surface_percentile))
  cat(sprintf("    %-12s : %.*f\n", "q_hat", digits, co$q_hat %||% NA_real_))
  # v0.7.1: co$band is a rendered consistency-band string (or "suppressed"
  # when the flag is divergent); the retired modal-band label and
  # confidence qualifier are gone.
  band_str <- if (is.null(co$band) ||
                  (length(co$band) == 1L && is.na(co$band))) "NA"
              else co$band
  cat(sprintf("    %-12s : %s\n", "band", band_str))

  cat("\n  delta (cross-coefficient asymmetry)\n")
  cat(sprintf("    %-16s : %.*f\n", "implied-q spread (pp)",
              max(digits - 1L, 1L), x$delta$delta_hat))
  cat(sprintf("    %-16s : %s\n", "flag", x$delta$flag))
  mn <- x$delta$matched_null
  dp <- x$delta$delta_percentile
  if (!is.null(mn) && is.finite(dp %||% NA_real_)) {
    cat(sprintf("    %-16s : %.1f (percentile on the matched null, k=%d, N=%d, q=%.2f)\n",
                "delta percentile", dp, mn$k, mn$N, mn$q))
  } else if (!is.null(mn)) {
    cat(sprintf("    %-16s : matched null (k=%d, N=%d, q=%.2f)\n",
                "matched null", mn$k, mn$N, mn$q))
  } else {
    cat(sprintf("    %-16s : unavailable (flag not calibrated)\n",
                "matched null"))
  }

  cat("\n  panel (full table)\n")
  pn <- x$panel
  for (i in seq_len(nrow(pn))) {
    band_str <- if (identical(x$delta$flag, "divergent")) {
      "band suppressed (divergent)"
    } else {
      .fmt_band_range(pn$band_lo[i], pn$band_hi[i],
                      pn$band_open_low[i], pn$band_open_high[i])
    }
    if (!nzchar(band_str)) band_str <- "NA"
    ref_used <- if ("reference_used" %in% names(pn)) pn$reference_used[i] else NA_character_
    cat(sprintf("    %-15s observed = %.*f  q_hat = %.*f  pct = %5.*f  %s  ref = %s\n",
                pn$coefficient[i],
                digits, pn$observed_value[i],
                digits, pn$q_hat[i],
                max(digits - 1L, 1L), pn$surface_percentile[i],
                band_str,
                if (is.na(ref_used)) "NA" else ref_used))
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
    for (n in x$notes) cat(.wrap_note_lines(n), sep = "\n")
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
