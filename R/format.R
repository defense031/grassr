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
  # v0.7.0: the flag is delta_hat's percentile on the matched null; the
  # implied pp cuts (95th/99th of that null) are shown as context.
  mn <- delta_list$matched_null
  pct <- delta_list$delta_percentile
  if (!is.null(mn) && is.finite(pct %||% NA_real_)) {
    sprintf("  matched null = (k=%d, N=%d, q=%.2f): delta_hat at the %.1f percentile%s%s",
            mn$k, mn$N, mn$q, pct,
            if (isTRUE(mn$snapped)) " [design snapped]" else "",
            if (isTRUE(mn$unstable_tail)) " [tail not stably invertible]" else "")
  } else if (identical(delta_list$thresholds_source, "user_supplied_legacy")) {
    thr <- delta_list$thresholds
    sprintf("  thresholds  = (%.2f, %.2f) [user-supplied, legacy pp cuts]",
            thr[["caution"]], thr[["divergent"]])
  } else {
    "  matched null = unavailable (flag not calibrated)"
  }
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
#' line per element. Used by `print.grass_card()` and by downstream embedders
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
  flag <- x$delta$flag
  delta_pp <- round(x$delta$delta_hat, digits = max(digits - 1L, 0L))

  lines <- character(0L)
  lines <- c(lines, "GRASS Report Card", "")
  sample_line <- sprintf("  sample      = %d raters, N = %d, pi_hat = %s",
                         k, N, pi_hat_str)
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
      is_icc <- identical(pn$coefficient[i], "icc")
      ref_used <- if ("reference_used" %in% names(pn)) pn$reference_used[i] else NA_character_
      ref_marker <- if (isTRUE(ref_used == "oracle-icc-fallback")) {
        "  [oracle-fallback]"
      } else if (isTRUE(ref_used == "oracle-icc-explicit")) {
        "  [oracle-icc]"
      } else if (is_icc) {
        "  [distribution-sensitive]"
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
               sprintf("  panel-agg.  = suppressed (divergent)"))
    lines <- c(lines,
               sprintf("  delta       = %s pp (divergent)",
                       formatC(delta_pp, format = "g", digits = max(digits, 2L))))
    lines <- c(lines, .fmt_thresholds_line(x$delta))
    lines <- c(lines, .fmt_bootstrap_tier_lines(x$delta, digits))
    # ---- Pairwise (paper Sec.3.3 recommended primary deliverable) ----
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
    # Aligned or caution: show all panel coefficients with their
    # surface percentiles, with the qualifier on the primary row.
    # The percentile is the categorical score per the v0.5 manifesto;
    # no band line.
    primary <- x$coefficient$primary
    qual_str <- if (is.null(x$coefficient$qualifier) ||
                    is.na(x$coefficient$qualifier)) "" else x$coefficient$qualifier
    pn <- x$panel
    name_w <- max(nchar(vapply(pn$coefficient, .coef_label, character(1))))
    name_w <- max(name_w, 6L)
    for (i in seq_len(nrow(pn))) {
      label <- .coef_label(pn$coefficient[i])
      val_str <- formatC(pn$observed_value[i], format = "f", digits = digits)
      pct_str <- .fmt_pp(pn$surface_percentile[i])
      is_primary <- identical(pn$coefficient[i], primary)
      is_icc <- identical(pn$coefficient[i], "icc")
      qual_suffix <- if (is_primary && nzchar(qual_str)) {
        sprintf("  (%s)", qual_str)
      } else ""
      ref_used <- if ("reference_used" %in% names(pn)) pn$reference_used[i] else NA_character_
      ref_marker <- if (isTRUE(ref_used == "oracle-icc-fallback")) {
        "  [oracle-fallback]"
      } else if (isTRUE(ref_used == "oracle-icc-explicit")) {
        "  [oracle-icc]"
      } else if (is_icc) {
        "  [distribution-sensitive]"
      } else ""
      primary_tag <- if (is_primary) "  <- primary" else ""
      lines <- c(lines,
                 sprintf("  %-*s = %s  ->  %s%s%s%s",
                         name_w, label, val_str, pct_str,
                         qual_suffix, ref_marker, primary_tag))
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
