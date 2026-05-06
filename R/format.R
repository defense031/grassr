# --------------------------------------------------------------------------
# format.grass_card -- v0.2.0 Target-2 Report Card text rendering
# --------------------------------------------------------------------------
# Returns a character vector (one element per line) rendering the four-field
# Report Card. Used by `print.grass_card()` (which cats the lines) and by
# downstream embedders (Markdown, knitr) that want the lines as data.

# Pretty-print a percentile in pp with an ordinal suffix ("48th", "30th",
# "92nd"). Centralised so the panel and the headline coefficient share the
# same look.
.fmt_pp <- function(p) {
  if (!is.finite(p)) return("    NA")
  pn <- round(p)
  suf <- if (pn %% 100 %in% c(11, 12, 13)) {
    "th"
  } else if (pn %% 10 == 1L) {
    "st"
  } else if (pn %% 10 == 2L) {
    "nd"
  } else if (pn %% 10 == 3L) {
    "rd"
  } else {
    "th"
  }
  paste0(sprintf("%2d", pn), suf, " percentile")
}

# Render a one-line summary of the active (caution, divergent) threshold
# pair plus a short source tag. Slot it in after the "delta = X pp (flag)"
# line. Source tags map: calibrated_at_k_N -> "calibrated", snapped_from_nearest
# -> "calibrated, k=k', N=N' (snapped)", default_fallback -> "modal-design
# default", user_supplied -> "user-supplied".
.fmt_thresholds_line <- function(delta_list) {
  thr <- delta_list$thresholds
  if (is.null(thr) || length(thr) != 2L) return(character(0L))
  src <- delta_list$thresholds_source %||% "user_supplied"
  src_str <- switch(
    src,
    calibrated_at_k_N             = "calibrated at this (k, N)",
    snapped_to_nearest_calibrated = "calibrated at nearest grid cell (snapped)",
    default_fallback              = "modal-design default (calibration table missing)",
    user_supplied                 = "user-supplied",
    src
  )
  sprintf("  thresholds  = (%.2f, %.2f) [%s]",
          unname(thr[1]), unname(thr[2]), src_str)
}

# Render the bootstrap-CI and tier-probability lines when bootstrap_delta_B
# was > 0 in the grass_report() call. Returns 0 or 2 character lines that
# slot in below the "delta = X pp (flag)" line in either branch.
.fmt_bootstrap_tier_lines <- function(delta_list, digits) {
  if (is.null(delta_list$bootstrap_B) || isTRUE(delta_list$bootstrap_B == 0L)) {
    return(character(0L))
  }
  ci  <- delta_list$delta_hat_ci
  tp  <- delta_list$tier_probabilities
  if (is.null(ci) || is.null(tp)) return(character(0L))
  ci_str <- sprintf("[%.*f, %.*f]",
                     digits, ci[1], digits, ci[2])
  tp_str <- sprintf("aligned %.2f | caution %.2f | divergent %.2f",
                    tp["aligned"], tp["caution"], tp["divergent"])
  c(sprintf("                95%% bootstrap CI: %s (B = %d)",
            ci_str, as.integer(delta_list$bootstrap_B)),
    sprintf("                tier prob:       %s", tp_str))
}

# Pretty label-mapping for the panel coefficient names.
.coef_label <- function(name) {
  switch(name,
    pabak          = "PABAK",
    mean_ac1       = "AC1",
    ac1            = "AC1",
    fleiss_kappa   = "Fleiss kappa",
    krippendorff_a = "alpha",
    icc            = "ICC",
    kappa          = "Cohen's kappa",
    name
  )
}

#' Format a grass_card as a character vector
#'
#' Returns the four-field Report Card rendered as a character vector, one
#' line per element. Used by [print.grass_card()] and by downstream embedders
#' (Markdown, knitr) that want the lines as data.
#'
#' @param x A `grass_card` object (returned by [grass_report()]).
#' @param digits Numeric. Decimals for observed coefficient values. Default 2.
#' @param ... Ignored.
#' @return A character vector (one element per line).
#' @export
format.grass_card <- function(x, digits = 2, ...) {
  num_fmt <- paste0("%.", digits, "f")

  ax <- x$sample$axis
  k  <- x$sample$k
  N  <- x$sample$N
  pi_hat_str <- formatC(x$sample$pi_hat, format = "f", digits = digits)
  tau2_str <- if (!is.null(x$sample$tau2_hat) &&
                  is.finite(x$sample$tau2_hat)) {
    formatC(x$sample$tau2_hat, format = "f", digits = max(digits + 1L, 3L))
  } else NULL

  flag <- x$delta$flag
  delta_pp <- round(x$delta$delta_hat, digits = max(digits - 1L, 0L))

  lines <- character(0L)
  lines <- c(lines, "GRASS Report Card", "")
  sample_line <- sprintf("  sample      = %d raters, N = %d, pi_hat = %s",
                         k, N, pi_hat_str)
  if (!is.null(tau2_str)) {
    sample_line <- paste0(sample_line, ", tau2_hat = ", tau2_str)
  }
  lines <- c(lines, sample_line)

  if (identical(flag, "divergent")) {
    # Divergent: show ALL panel coefficients.
    pn <- x$panel
    primary <- x$coefficient$primary
    name_w <- max(nchar(vapply(pn$coefficient, .coef_label, character(1))))
    name_w <- max(name_w, nchar("Cohen's kappa"))
    for (i in seq_len(nrow(pn))) {
      label <- .coef_label(pn$coefficient[i])
      val_str <- formatC(pn$observed_value[i], format = "f", digits = digits)
      pct_str <- .fmt_pp(pn$surface_percentile[i])
      is_clamped <- isTRUE(pn$clamped[i])
      ref_used <- if ("reference_used" %in% names(pn)) pn$reference_used[i] else NA_character_
      ref_marker <- if (isTRUE(ref_used == "oracle-icc-fallback")) {
        "  [oracle-fallback]"
      } else if (isTRUE(ref_used == "oracle-icc-explicit")) {
        "  [oracle-icc]"
      } else {
        ""
      }
      marker <- if (identical(pn$coefficient[i], primary)) {
        paste0("  <- primary", ref_marker)
      } else if (is_clamped) {
        paste0("  [clamped - excluded from delta]", ref_marker)
      } else {
        ref_marker
      }
      lines <- c(lines,
                 sprintf("  %-*s = %s  ->  %s%s",
                         name_w, label, val_str, pct_str, marker))
    }
    lines <- c(lines,
               sprintf("  band        = suppressed"))
    lines <- c(lines,
               sprintf("  delta       = %s pp (divergent)",
                       formatC(delta_pp, format = "g", digits = max(digits, 2L))))
    lines <- c(lines, .fmt_thresholds_line(x$delta))
    lines <- c(lines, .fmt_bootstrap_tier_lines(x$delta, digits))
    # ---- Pairwise (paper §3.3 recommended primary deliverable) ----
    if (!is.null(x$pairwise)) {
      pw <- x$pairwise
      kk <- nrow(pw$pabak_matrix)
      lines <- c(lines, "",
                 "  pairwise PABAK / surface percentile (lower / upper):")
      rn <- rownames(pw$pabak_matrix)
      hdr_row <- paste0("    ", paste0(format(rn, width = 6), collapse = " "))
      lines <- c(lines, hdr_row)
      for (ii in seq_len(kk)) {
        cells <- character(kk)
        for (jj in seq_len(kk)) {
          if (ii == jj) {
            cells[jj] <- "  --  "
          } else if (ii > jj) {
            cells[jj] <- formatC(pw$pabak_matrix[ii, jj], format = "f",
                                  digits = digits, width = 6)
          } else {
            pct <- pw$percentile_matrix[ii, jj]
            cells[jj] <- if (is.finite(pct)) sprintf(" %3.0f%% ", pct) else "  NA  "
          }
        }
        lines <- c(lines,
          paste0("  ", format(rn[ii], width = 4), paste0(cells, collapse = " ")))
      }
      if (nrow(pw$pooled_per_rater) > 0L &&
          any(is.finite(pw$pooled_per_rater$se_tilde))) {
        lines <- c(lines, "",
                   "  per-rater vs panel-majority of OTHER raters (pooled-reference):")
        pp <- pw$pooled_per_rater
        for (ii in seq_len(nrow(pp))) {
          se_str <- if (is.finite(pp$se_tilde[ii]))
            formatC(pp$se_tilde[ii], format = "f", digits = digits) else "  NA"
          sp_str <- if (is.finite(pp$sp_tilde[ii]))
            formatC(pp$sp_tilde[ii], format = "f", digits = digits) else "  NA"
          lines <- c(lines,
            sprintf("    %-4s Se_tilde = %s  Sp_tilde = %s  (n_pos = %d, n_neg = %d, excl = %d)",
                    pp$rater[ii], se_str, sp_str,
                    pp$n_pool_pos[ii], pp$n_pool_neg[ii],
                    pp$n_pool_excluded[ii]))
        }
      }
    }
    if (!is.null(x$per_rater) && nrow(x$per_rater) > 0L) {
      lines <- c(lines, "")
      bound_only <- isTRUE(any(x$per_rater$bound_only))
      hdr <- if (bound_only) {
        "  per-rater (Hui-Walter bounds, k = 2; alongside pairwise):"
      } else {
        "  per-rater (latent-class fit; alongside pairwise):"
      }
      lines <- c(lines, hdr)
      pr <- x$per_rater
      for (i in seq_len(nrow(pr))) {
        if (isTRUE(pr$bound_only[i])) {
          lines <- c(lines,
                     sprintf("    %-4s Se in [%s, %s]   Sp in [%s, %s]",
                             pr$rater[i],
                             formatC(pr$se_lower[i], format = "f", digits = digits),
                             formatC(pr$se_upper[i], format = "f", digits = digits),
                             formatC(pr$sp_lower[i], format = "f", digits = digits),
                             formatC(pr$sp_upper[i], format = "f", digits = digits)))
        } else {
          lines <- c(lines,
                     sprintf("    %-4s Se = %s  (%s, %s)   Sp = %s  (%s, %s)",
                             pr$rater[i],
                             formatC(pr$se_hat[i], format = "f", digits = digits),
                             formatC(pr$se_lower[i], format = "f", digits = digits),
                             formatC(pr$se_upper[i], format = "f", digits = digits),
                             formatC(pr$sp_hat[i], format = "f", digits = digits),
                             formatC(pr$sp_lower[i], format = "f", digits = digits),
                             formatC(pr$sp_upper[i], format = "f", digits = digits)))
        }
      }
    }
  } else {
    # Aligned or caution: show only the primary coefficient + band + delta.
    coef <- x$coefficient
    label <- .coef_label(coef$primary)
    val_str <- formatC(coef$observed_value, format = "f", digits = digits)
    pct_str <- .fmt_pp(coef$surface_percentile)
    name_w <- max(nchar(label), 6L)
    primary_ref <- if (!is.null(x$panel) && "reference_used" %in% names(x$panel)) {
      i <- which(x$panel$coefficient == coef$primary)
      if (length(i) == 1L) x$panel$reference_used[i] else NA_character_
    } else NA_character_
    primary_marker <- if (isTRUE(primary_ref == "oracle-icc-fallback")) {
      "  [oracle-fallback]"
    } else if (isTRUE(primary_ref == "oracle-icc-explicit")) {
      "  [oracle-icc]"
    } else ""
    lines <- c(lines,
               sprintf("  %-*s = %s  ->  %s%s",
                       name_w, label, val_str, pct_str, primary_marker))
    band_str <- if (is.null(coef$band) || is.na(coef$band)) "NA" else coef$band
    qual_str <- if (is.null(coef$qualifier) || is.na(coef$qualifier)) "" else coef$qualifier
    if (nzchar(qual_str)) {
      lines <- c(lines,
                 sprintf("  band        = %s (%s)",
                         band_str, qual_str))
    } else {
      lines <- c(lines,
                 sprintf("  band        = %s",
                         band_str))
    }
    lines <- c(lines,
               sprintf("  delta       = %s pp (%s)",
                       formatC(delta_pp, format = "g", digits = max(digits, 2L)),
                       flag))
    lines <- c(lines, .fmt_thresholds_line(x$delta))
    lines <- c(lines, .fmt_bootstrap_tier_lines(x$delta, digits))
  }

  if (length(x$notes) > 0L) {
    lines <- c(lines, "")
    lines <- c(lines, "  Notes:")
    for (n in x$notes) {
      lines <- c(lines, paste0("    - ", n))
    }
  }

  lines <- c(lines, "")
  lines <- c(lines, "  See `summary(...)` for full panel and CI details.")
  lines <- c(lines, "  See `plot(...)` for a surface-position visualization.")
  lines
}

#' One-line paper-ready summary of a grass report
#'
#' Returns a single character string suitable for dropping into a
#' manuscript line or a results table cell. Contains the three agreement
#' metrics (with a confidence interval on kappa), the skew diagnostics
#' (PI, BI), sample size, and regime.
#'
#' @section Deprecated:
#' `grass_format_report()` is deprecated in grass 0.2.0. The Target-2
#' Report Card is rendered by [print.grass_card()] / [format.grass_card()]
#' on the object returned by `grass_report(ratings = Y)`. See
#' `vignette("reporting-card")` and [grass_report()].
#'
#' @param x A `grass_result` object.
#' @param digits Number of decimal places for rounding. Default 2.
#' @param ascii If `TRUE` (the default), emit `kappa` instead of the Unicode
#'   `\u03ba`. Safer for Slack, Markdown, and non-UTF-8 locales. Pass
#'   `ascii = FALSE` for a Unicode-rendering context (e.g., RStudio console).
#' @param ci_width If `TRUE`, append the kappa CI width and a one-word
#'   descriptor (`tight` / `moderate` / `wide`) based on Wilson-logit
#'   half-width cutpoints at 0.10 and 0.20. These cutpoints are calibrated
#'   for kappa CIs; Se/Sp CIs at the same N typically run roughly 2-2.5x
#'   tighter. Default `FALSE`.
#'
#' @return A single character string.
#' @export
#'
#' @examples
#' tab <- matrix(c(88, 10, 14, 88), nrow = 2,
#'               dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
#' grass_format_report(grass_report(tab, format = "matrix"))
#' grass_format_report(grass_report(tab, format = "matrix"), ci_width = TRUE)
grass_format_report <- function(x, digits = 2, ascii = TRUE, ci_width = FALSE) {
  msg_once(
    "deprecate_grass_format_report",
    paste0(
      "`grass_format_report()` is deprecated in grass 0.2.0. ",
      "The Target-2 Report Card is rendered by `print(grass_report(ratings = Y))`. ",
      "See `vignette('reporting-card')` and `?grass_report`."
    )
  )
  if (!inherits(x, "grass_result")) {
    stop("`x` must be a grass_result, as returned by grass_report().",
         call. = FALSE)
  }
  v <- x$metrics$values
  kappa_name <- if (isTRUE(ascii)) "kappa" else "\u03ba"
  fmt <- paste0("%.", digits, "f")
  ci  <- sprintf(paste0("[", fmt, ", ", fmt, "]"),
                 unname(v["kappa_wilson_lower"]),
                 unname(v["kappa_wilson_upper"]))
  base <- sprintf(
    paste0("%s = ", fmt, " ", "%s, PABAK = ", fmt, ", AC1 = ", fmt,
           ", PI = ", fmt, ", BI = ", fmt, ", N = %d, prevalence = ", fmt,
           ", %s regime"),
    kappa_name,
    unname(v["kappa"]),
    ci,
    unname(v["PABAK"]),
    unname(v["AC1"]),
    unname(v["prevalence_index"]),
    unname(v["bias_index"]),
    x$metrics$n,
    x$prevalence,
    x$regime
  )
  if (isTRUE(ci_width)) {
    lo <- unname(v["kappa_wilson_lower"])
    hi <- unname(v["kappa_wilson_upper"])
    width <- hi - lo
    # Round to display precision before comparing. (0.50 - 0.40)/2 computes
    # to 0.0499999... in binary float; a naive `half < 0.10` flips the
    # descriptor at exact cutpoints. Rounding to 4 decimals matches the
    # precision the user sees in the printed interval and removes the
    # boundary-flip surprise.
    half  <- round(width / 2, 4)
    # Cutpoints calibrated for kappa Wilson-logit CIs, which run roughly
    # 2-2.5x wider than component Se/Sp CIs at the same N. The 0.10 / 0.20
    # half-width thresholds match empirical kappa CI spread rather than the
    # Se/Sp conventions (0.05 / 0.10) from clinical diagnostic accuracy.
    descriptor <- if (is.na(half))       "unknown"
                  else if (half <  0.10) "tight"
                  else if (half <  0.20) "moderate"
                  else                   "wide"
    base <- paste0(base,
                   sprintf(paste0("; kappa CI width = ", fmt, ", %s"),
                           width, descriptor))
  }
  base
}
