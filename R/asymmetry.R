# Asymmetry diagnostics for the GRASS framework.
#
# One exported function lives here:
#
#   * `check_rater_asymmetry()` -- REMOVED in 0.8.0. It shipped in CRAN
#     0.7.4 with a stipulated three-tier (0.05 / 0.10) scheme the framework
#     retired. Legacy `check_asymmetry(se =, sp =)` calls now error with a
#     pointer to latent_class_fit().
#
#   * `check_asymmetry(ratings, ...)` -- the new ratings-input diagnostic
#     defined by Sec.4.3 of the v0.2.0 paper-alignment design doc. Computes
#     the cross-coefficient IMPLIED-QUALITY spread `delta_hat` (in quality
#     percentage points): each agreement-family coefficient inverts to its
#     own implied panel quality `q_hat` on the shared (q, pi_+) reference,
#     and `delta_hat` is their max-min spread. The flag is `delta_hat`'s
#     percentile on the matched (k, N, q_hat) null ECDF (>= 95th caution,
#     >= 99th divergent); the per-(k, N) size-alpha threshold table is
#     retired (0.7.0/0.7.1). Class: `grass_asymmetry_panel`.
#
# v0.5.0 ICC scope decision (2026-05-05). delta_hat is computed over the
# AGREEMENT FAMILY ONLY -- PABAK, mean AC1, Fleiss kappa -- because each has
# a closed-form reference depending on (q, pi_+) only and
# is therefore DGP-robust at the panel level. ICC is reported alongside on
# the panel rows but does NOT enter delta_hat: ICC's reference surface
# depends on the full subject-prevalence distribution F, and a panel whose
# true F does not match the bundled logit-normal reference incurs surface-
# percentile drift on the order of 20 pp at small designs (sanity probe,
# 2026-05-05). The agreement family is mutually distribution-robust, so
# cross-family spread is a clean SPLIT-BIAS detector: it fires when raters
# tilt in different directions across (Se, Sp) so the coefficients
# disagree on surface position. Shared/uniform bias (every rater tilts
# the same direction) produces small delta_hat with low minimum surface
# percentile and is detected by the percentile, not by delta_hat.
#
# v0.6.0 alpha scope decision (2026-07-02, paper Pass 6). Krippendorff's
# alpha is removed from the Report Card panel and from delta_hat: in the
# binary fully-crossed case alpha coincides with Fleiss kappa
# asymptotically (paper App A.1; empirically median |Fleiss - alpha| =
# 0.00024 across 10,140 calibration cells), so it added a redundant row
# and could tip borderline delta_hat readings on small-sample noise
# alone. obs_krippendorff_alpha() and
# position_on_surface(metric = "krippendorff_a") remain available for
# manual use. delta_hat now spans the 3-coefficient agreement family
# (PABAK, mean AC1, Fleiss kappa) as an implied-quality spread; its
# matched null is calibrated for that 3-coefficient set.

# Coefficients that enter delta_hat. Anything else (currently just `icc`) is
# reported on the panel with `in_delta_hat = FALSE`.
.DELTA_AGREEMENT_COEFS <- c("pabak", "mean_ac1", "fleiss_kappa")
#
# If `check_asymmetry()` is called with the pre-0.2.0 signature
# (`se = ...`, `sp = ...`) it errors with a pointer to latent_class_fit().


# ---------------------------------------------------------------------------
# (2) New: check_asymmetry(ratings, ...)
# ---------------------------------------------------------------------------

#' Cross-coefficient panel asymmetry diagnostic
#'
#' `check_asymmetry()` takes an N x k binary ratings matrix and returns a
#' scalar `delta_hat` (in quality percentage points, "pp"): the max-min
#' spread of the *implied panel qualities* across the three agreement
#' coefficients (PABAK, mean AC1, Fleiss kappa). Each coefficient inverts
#' to its own `q_hat` on the shared (q, pi_+) reference; if the
#' calibration DGP held exactly, all three would imply the same quality,
#' so the spread measures cross-coefficient model discordance in
#' interpretable units of quality. ICC never enters `delta_hat`, since
#' its reference is distribution-sensitive in ways the agreement family
#' is not.
#'
#' `delta_hat` is a *split-bias* detector. It fires when raters tilt in
#' different directions across (Se, Sp) -- e.g., one rater high-Se / low-Sp,
#' another high-Sp / low-Se -- because the three coefficients respond to
#' heterogeneous per-rater behavior differently and end up implying
#' different panel qualities. The framework's other failure mode,
#' *shared/uniform bias* (every rater tilts the same direction, e.g., a
#' panel trained on one protocol all favoring specificity over
#' sensitivity), produces uniform degradation across coefficients: small
#' `delta_hat`, low implied quality. Shared bias is detected by the
#' panel's implied quality (and its consistency band), not by `delta_hat`.
#' A divergent flag therefore identifies a specific kind of disagreement
#' -- cross-coefficient inversion from heterogeneous rater behavior -- and
#' routes the user to the per-rater pairwise PABAK matrix and
#' pooled-reference (Se_tilde, Sp_tilde) diagnostic.
#'
#' Each coefficient is positioned on its reference surface via
#' [position_on_surface()], which reports its implied `q_hat`. The panel
#' diagnostic is `delta_hat = max(q_hat) - min(q_hat)` (in pp of quality),
#' computed over agreement-family coefficients whose observed value sits
#' within the achievable range of their reference surface (see
#' *Surface-envelope clamp* below).
#'
#' @section Surface-envelope clamp (v0.2.1+):
#' If an observed coefficient value falls outside the achievable range of
#' its reference surface at the study's design `(pi_hat, k, N)`, the
#' inversion to `q_hat` clamps to the boundary. Including such clamped
#' implied qualities in the max-min `delta_hat` would inflate the panel
#' spread purely because of the clamp, not because the panel disagrees on
#' quality. Since v0.2.1 the function therefore *excludes* clamped
#' coefficients from `delta_hat` whenever at least two unclamped
#' agreement-family coefficients remain. The affected coefficients are
#' still shown in the returned `panel` data.frame with `clamped = TRUE`,
#' and a note in `$notes` names which coefficients were excluded. This
#' matters most often for ICC (which never enters `delta_hat` anyway) and
#' at designs beyond the bundled reference range, where a coefficient
#' clamps to the achievable boundary. If fewer than two unclamped
#' agreement-family coefficients remain, `delta_hat` falls back to the raw
#' spread including clamped values and the note records the fallback.
#'
#' @section Flag from the matched null (v0.7.0/0.7.1):
#' The per-(k, N) size-alpha threshold table is retired. The flag is
#' `delta_hat`'s percentile on the null distribution of `delta_hat` at the
#' matched (k, N, q_hat) cell of the bundled `delta_null_ecdf`, with the
#' cut convention >= 95th caution, >= 99th divergent. The panel's `q_hat`
#' is resolved first (median of the agreement-family implied qualities),
#' then the matched null cell is looked up; the reported
#' `thresholds` carry the implied pp cuts (95th/99th of that null) as
#' context, and `thresholds_source` records how the flag was resolved. The
#' three flags are:
#'
#' - **`aligned`** (below the 95th percentile of the matched null): the
#'   panel agrees on the implied quality. Any single coefficient is a
#'   stable summary; the primary coefficient (Table 2) carries the
#'   headline.
#' - **`caution`** (>= 95th, < 99th): the panel is mildly inconsistent.
#'   Report the primary coefficient with a caution flag and the
#'   `delta_hat` value.
#' - **`divergent`** (>= 99th): no single coefficient is a stable summary.
#'   Use [latent_class_fit()] to recover per-rater `(Se, Sp)` and report
#'   those instead.
#'
#' The new `check_asymmetry(ratings, ...)` signature replaces the
#' v0.1.x `check_asymmetry(se, sp, ...)` per-rater signature. For
#' Calling `check_asymmetry()` with the pre-0.2.0 `se = ...` and
#' `sp = ...` arguments is an error; per-rater sensitivity and
#' specificity come from [latent_class_fit()].
#'
#' @param ratings User input: an `N x k` binary matrix, an `N` x `k`
#'   data.frame whose columns are 0/1 / logical / 2-level factor, or a
#'   list of two equal-length 0/1 vectors (k = 2 paired form); the Quick
#'   start section of `vignette("grassr")` lists the accepted shapes.
#' @param axis `"inter"` (default) or `"intra"`. Selects the surface
#'   family.
#' @param occasion Reserved for `axis = "intra"` (a vector / factor
#'   identifying viewing occasion); ignored when `axis = "inter"`.
#' @param fit_icc If `FALSE`, skip the `lme4::glmer` fit behind `icc` and drop
#'   ICC from the panel. `icc` never enters `delta_hat`, and the fit draws no
#'   random numbers, so a caller that reports only `delta_hat` and the implied
#'   qualities gets identical results at a fraction of the cost. The Monte
#'   Carlo null loop sets this; interactive users should not.
#' @param ... Forwarded to [position_on_surface()] (e.g. `reference_type`).
#'
#' @return An S3 object of class `grass_asymmetry_panel` with fields:
#' - `delta_hat`: scalar implied-quality spread, in pp of quality
#' - `delta_percentile`: `delta_hat`'s percentile on the matched
#'   (k, N, q_hat) null ECDF (`NA` if the null is uncalibrated at the
#'   design)
#' - `flag`: one of `"aligned"`, `"caution"`, `"divergent"`
#' - `matched_null`: list describing the matched null cell
#'   (`k`, `N`, `q`, `prev` — the bridged true-prevalence estimate the
#'   lookup conditioned on, `prev_apparent` — the panel's raw positive
#'   rate, `prev_bridged`, `q_hat_panel`, `n_draws`,
#'   `snapped`, `interpolated`, `unstable_tail`), or `NULL` if
#'   uncalibrated
#' - `thresholds`: named numeric vector of the implied (caution, divergent)
#'   pp cuts (95th/99th of the matched null)
#' - `thresholds_source`: one of `"matched_null_ecdf"`,
#'   `"not_applicable_k2"`, `"not_calibrated"`
#' - `panel`: data.frame with `coefficient`, `observed`, `implied_q`,
#'   `percentile_pp` (pooled percentile), `clamped`, `in_delta_hat`
#' - `notes`: character vector of unique caveats from the underlying
#'   surface positioning calls (e.g. nearest-neighbor gaps, ICC
#'   unavailability, matched-null provenance)
#'
#' @seealso [latent_class_fit()] for the divergent-branch recovery of per-rater
#'   `(Se, Sp)`; [position_on_surface()] for the underlying surface
#'   positioning.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' # Build a 5x200 symmetric panel -- should print as 'aligned'.
#' Y <- matrix(rbinom(5 * 200, 1, 0.30), nrow = 200, ncol = 5)
#' check_asymmetry(Y)
check_asymmetry <- function(ratings,
                            axis = c("inter", "intra"),
                            occasion = NULL,
                            fit_icc = TRUE,
                            ...) {
  # ---- Soft-deprecation dispatch for OLD `check_asymmetry(se, sp, ...)` ---
  # Detection: the OLD signature was `check_asymmetry(se, sp, ...)` where
  # `se` and `sp` are numeric vectors of length k (per-rater Se/Sp). The
  # NEW signature takes `ratings` as an N x k matrix / data.frame / list.
  # So a call routes to the OLD path if any of:
  #   (a) the user supplied `se =` and / or `sp =` as named args, or
  #   (b) the first positional arg (now bound to `ratings`) is a numeric
  #       vector that is NOT a matrix / data.frame / list -- i.e. the OLD
  #       per-rater Se vector got positionally bound to `ratings`.
  cl <- match.call()
  dots <- list(...)
  has_se <- "se" %in% names(dots)
  has_sp <- "sp" %in% names(dots)
  ratings_supplied <- !missing(ratings)

  ratings_is_per_rater_vec <- ratings_supplied &&
    is.numeric(ratings) && !is.matrix(ratings) &&
    !is.data.frame(ratings) && !is.list(ratings)

  is_legacy_call <- has_se || has_sp || ratings_is_per_rater_vec

  if (is_legacy_call) {
    # Block dual-input: ratings supplied as a real ratings matrix AND
    # per-rater se/sp also supplied.
    if (ratings_supplied && !ratings_is_per_rater_vec && (has_se || has_sp)) {
      stop("check_asymmetry(): you supplied both `ratings = ...` and ",
           "per-rater `se =` / `sp =`. Pick one. For panel-level percentile ",
           "spread use `check_asymmetry(ratings = Y)`; for the per-rater ",
           "Se/Sp table use `latent_class_fit()`.",
           call. = FALSE)
    }
    stop("check_asymmetry(se = ..., sp = ...) is no longer supported. ",
         "check_asymmetry() takes a rating matrix, `check_asymmetry(ratings = Y)`. ",
         "Per-rater sensitivity and specificity come from `latent_class_fit()`.",
         call. = FALSE)
  }

  # ---- Normal path: ratings-input panel diagnostic -----------------------
  axis <- match.arg(axis)
  Y <- normalize_ratings(ratings)
  panel_obs <- compute_panel(Y, axis = axis, occasion = occasion,
                             fit_icc = fit_icc)

  # The flag comes from delta_hat's percentile on the matched (k, N, q_hat)
  # null ECDF (>= 95th caution, >= 99th divergent), resolved after the
  # panel's q_hat is known.
  thresholds_note <- ""

  # The panel-name -> surface-metric mapping. compute_panel() returns
  # `ac1` (panel-friendly short name) but position_on_surface() expects
  # `mean_ac1`; Cohen's `kappa` (k=2 only) has no surface and is dropped.
  surface_metric_for <- c(
    pabak          = "pabak",
    ac1            = "mean_ac1",
    fleiss_kappa   = "fleiss_kappa",
    icc            = "icc"
  )
  panel_keep <- intersect(names(panel_obs), names(surface_metric_for))
  panel_obs <- panel_obs[panel_keep]
  # Re-key so the printed panel uses the surface metric names (paper labels).
  names(panel_obs) <- unname(surface_metric_for[panel_keep])
  # Drop non-finite panel entries before positioning. compute_panel()
  # returns icc = NA when lme4 (Suggests) is unavailable or glmer fails;
  # positioning an NA is an error, and a missing coefficient must degrade
  # to "absent from the panel", never to a hard failure (v0.7.0 fix).
  panel_obs <- panel_obs[vapply(panel_obs, function(v)
    is.numeric(v) && is.finite(v), logical(1L))]

  positions <- lapply(names(panel_obs), function(m) {
    position_on_surface(ratings = Y, metric = m, ...)
  })
  names(positions) <- names(panel_obs)

  # v0.7.1 (Option B, ratified 2026-07-05): delta_hat is the spread of the
  # IMPLIED QUALITIES, in quality percentage points. Each agreement-family
  # coefficient inverts to its own q_hat on the shared (q, pi_+) reference;
  # if the calibration DGP held exactly, all family members would imply the
  # same quality, so the spread measures cross-coefficient model
  # discordance in interpretable units. (The previous definition -- spread
  # of surface percentiles -- ran through the retired nearest-cell
  # percentile machinery, whose sawtooth inflated delta_hat with
  # quantization noise; see design/v0.7.1_position_redesign.md.)
  qhat_pp <- vapply(
    positions,
    function(p) p$q_hat * 100,
    numeric(1)
  )

  # delta_hat spans the AGREEMENT FAMILY only (PABAK, mean AC1, Fleiss
  # kappa). ICC is reported on the panel but excluded from delta_hat: see
  # the v0.5.0 scope note at the top of this file. Two further exclusions
  # inside the agreement family:
  #   1. Reference-envelope clamps: a coefficient whose obs_value falls
  #      outside the achievable range of its reference curve gets q_hat
  #      clamped to the boundary. Including those would inflate delta_hat
  #      purely from the clamp. Excluded when >= 2 unclamped remain.
  #   2. (Implicit) ICC, regardless of clamp status: never in delta_hat.
  clamped <- vapply(
    positions,
    function(p) any(grepl("q_hat clamped", p$notes, fixed = TRUE)),
    logical(1L)
  )

  # Set in_delta_hat = TRUE only for agreement-family coefficients that are
  # not currently clamped. ICC is always FALSE.
  in_delta_set <- names(positions) %in% .DELTA_AGREEMENT_COEFS & !clamped
  names(in_delta_set) <- names(positions)

  clamp_note <- character(0L)
  if (sum(in_delta_set) >= 2L) {
    delta_hat <- max(qhat_pp[in_delta_set]) - min(qhat_pp[in_delta_set])
    # Note clamped-but-otherwise-eligible agreement coefficients
    excluded_for_clamp <- names(positions)[
      names(positions) %in% .DELTA_AGREEMENT_COEFS & clamped
    ]
    if (length(excluded_for_clamp)) {
      clamp_note <- sprintf(
        "delta_hat (implied-quality spread) over %d agreement-family coefficients; excluded due to reference-envelope clamp: %s. ICC is shown but not included in delta_hat; its reference depends on how positive probability is spread across subjects, while the other three depend only on prevalence and quality.",
        sum(in_delta_set),
        paste(excluded_for_clamp, collapse = ", ")
      )
    } else if ("icc" %in% names(positions)) {
      # No clamps; just remind that ICC isn't in delta_hat.
      clamp_note <- "delta_hat is the implied-quality spread over the agreement family (PABAK, mean AC1, Fleiss kappa). ICC is shown but not included in delta_hat; its reference depends on how positive probability is spread across subjects, while the other three depend only on prevalence and quality."
    }
  } else {
    # Fall back to all available agreement-family implied qualities, even
    # if clamped, so the spread is still defined.
    agreement_idx <- names(positions) %in% .DELTA_AGREEMENT_COEFS
    if (sum(agreement_idx) >= 2L) {
      delta_hat <- max(qhat_pp[agreement_idx]) - min(qhat_pp[agreement_idx])
      clamp_note <- sprintf(
        "Fewer than 2 unclamped agreement-family coefficients; delta_hat uses raw implied-quality spread over %d agreement coefficients including clamped: %s.",
        sum(agreement_idx),
        paste(names(positions)[agreement_idx], collapse = ", ")
      )
    } else {
      delta_hat <- NA_real_
      clamp_note <- "delta_hat undefined: fewer than 2 agreement-family coefficients available."
    }
  }

  # ---- Matched-null resolution + flag ------------------------------------
  # k = 2: delta_hat is structurally uninformative (the two-coefficient
  # family implies identical quality by construction; the Option-B null
  # is a point mass at zero on all 2.75M draws). No lookup, no flag —
  # report not_applicable and route asymmetry assessment to the k = 2
  # identifiable-bounds / pairwise path. This restores the paper's
  # original k = 2 position (v0.7.1; see design/v0.7.1_position_redesign.md).
  k2_degenerate <- ncol(Y) == 2L
  qh <- vapply(positions[names(positions) %in% .DELTA_AGREEMENT_COEFS],
               function(p) p$q_hat, numeric(1L))
  q_hat_panel <- stats::median(qh[is.finite(qh)])
  # The null grid is indexed by TRUE prevalence; mean(Y) is the APPARENT
  # positive rate, which sits closer to 0.5 than truth whenever raters
  # err. Under the reference model the two are related by
  #   pi_apparent = pi (2q - 1) + (1 - q),
  # so the lookup conditions on the closed-form inversion at the panel's
  # estimated quality. Querying at apparent prevalence reads a too-narrow
  # null and inflates realized flag size ~2x at skewed prevalence
  # (G2 Tier B, 2026-07-24); the bridged query holds it at nominal. At
  # low estimated quality the inversion divisor 2q - 1 degenerates, so
  # the raw rate is kept there (the null is wide at low q regardless).
  pi_apparent <- mean(Y)
  pi_bridged <- if (is.finite(q_hat_panel) && (2 * q_hat_panel - 1) > 0.10)
    min(max((pi_apparent - (1 - q_hat_panel)) / (2 * q_hat_panel - 1), 0), 1)
  else pi_apparent
  null_cell <- if (!k2_degenerate && is.finite(q_hat_panel))
    lookup_delta_null(k = ncol(Y), N = nrow(Y), q_hat = q_hat_panel,
                      pi_hat = pi_bridged)
  else NULL

  delta_percentile <- NA_real_
  implied_cuts <- c(caution = NA_real_, divergent = NA_real_)
  matched_null <- NULL
  if (!is.null(null_cell)) {
    delta_percentile <- delta_null_percentile(delta_hat, null_cell)
    i95 <- which(abs(null_cell$probs - 0.95) < 1e-9)
    i99 <- which(abs(null_cell$probs - 0.99) < 1e-9)
    implied_cuts <- c(caution = unname(null_cell$values[i95]),
                      divergent = unname(null_cell$values[i99]))
    matched_null <- list(k = null_cell$k, N = null_cell$N, q = null_cell$q,
                         prev = null_cell$prev,
                         prev_apparent = pi_apparent,
                         prev_bridged = !identical(pi_bridged, pi_apparent),
                         q_hat_panel = unname(q_hat_panel),
                         n_draws = null_cell$n_draws,
                         snapped = null_cell$snapped,
                         interpolated = null_cell$interpolated,
                         unstable_tail = null_cell$unstable_tail)
    thresholds_note <- sprintf(
      "flag from delta_hat's percentile on the matched null (k=%d, N=%d, q=%.2f%s; %s draws)%s%s%s.",
      null_cell$k, null_cell$N, null_cell$q,
      sprintf(", prev=%.2f", null_cell$prev),
      format(null_cell$n_draws, big.mark = ","),
      if (null_cell$snapped) "; design snapped to the calibrated grid" else "",
      if (null_cell$interpolated)
        "; null interpolated between calibrated grid nodes" else "",
      if (null_cell$unstable_tail)
        "; this cell's extreme tail is flagged as not stably invertible (percentile reading unaffected)" else "")
  } else if (k2_degenerate) {
    thresholds_note <- paste0(
      "delta_hat is not applicable at k = 2: the two-coefficient agreement ",
      "family (PABAK, AC1) implies identical panel quality by construction, ",
      "so cross-coefficient discordance cannot be observed. Use the k = 2 ",
      "identifiable bounds and pairwise path for asymmetry assessment.")
  } else {
    thresholds_note <- "delta_null_ecdf unavailable; flag not calibrated."
  }

  flag <- if (k2_degenerate) {
    "not_applicable"
  } else {
    delta_flag_from_percentile(delta_percentile,
      if (!is.null(null_cell)) null_cell$conventions
      else c(caution = 0.95, divergent = 0.99))
  }
  thresholds_source <- if (k2_degenerate) "not_applicable_k2"
                       else if (!is.null(null_cell)) "matched_null_ecdf"
                       else "not_calibrated"
  report_cuts <- implied_cuts

  combined_notes <- unique(c(
    unlist(lapply(positions, `[[`, "notes")),
    clamp_note,
    if (nzchar(thresholds_note)) thresholds_note else character(0L)
  ))

  out <- list(
    delta_hat  = unname(delta_hat),
    delta_percentile = unname(delta_percentile),
    flag       = flag,
    matched_null = matched_null,
    thresholds = report_cuts,
    thresholds_source = thresholds_source,
    panel      = data.frame(
      coefficient   = names(panel_obs),
      observed      = unlist(panel_obs, use.names = FALSE),
      implied_q     = unname(qhat_pp) / 100,
      percentile_pp = vapply(positions, function(p)
        100 * (p$percentile %||% NA_real_), numeric(1L)),
      clamped       = unname(clamped),
      in_delta_hat  = unname(names(panel_obs) %in% .DELTA_AGREEMENT_COEFS),
      stringsAsFactors = FALSE
    ),
    notes      = combined_notes
  )
  class(out) <- c("grass_asymmetry_panel", "list")
  out
}

#' @export
print.grass_asymmetry_panel <- function(x, digits = 1, ...) {
  cat("GRASS panel asymmetry diagnostic\n\n", sep = "")
  cat(sprintf("  delta_hat = %.*f pp  (spread of the implied panel qualities)\n",
              digits, x$delta_hat))
  if (is.finite(x$delta_percentile %||% NA_real_)) {
    cat(sprintf("  flag      = %s  (%.1f percentile of matched null: k=%d, N=%d, q=%.2f)\n",
                x$flag, x$delta_percentile,
                x$matched_null$k, x$matched_null$N, x$matched_null$q))
  } else {
    cat(sprintf("  flag      = %s\n", x$flag))
  }
  cat("\n  panel:\n")
  cat("    coefficient        observed   implied q   pooled pctile   in delta_hat\n")
  pn <- x$panel
  for (i in seq_len(nrow(pn))) {
    in_d <- if (isTRUE(pn$in_delta_hat[i])) "yes"
            else if (identical(pn$coefficient[i], "icc"))
              "no [distribution-sensitive]"
            else "no"
    cat(sprintf("    %-18s %.2f       %.3f       %.*f            %s\n",
                pn$coefficient[i],
                pn$observed[i],
                pn$implied_q[i],
                digits,
                pn$percentile_pp[i],
                in_d))
  }
  if (x$flag == "divergent") {
    cat("\n  Note: at the divergent flag, no single coefficient is a stable summary.\n")
    cat("        Use latent_class_fit() to recover per-rater (Se, Sp).\n")
  } else if (x$flag == "caution") {
    cat("\n  Note: at the caution flag, panel coefficients are mildly inconsistent.\n")
    cat("        Inspect the panel before relying on a single headline.\n")
  }
  if (length(x$notes) > 0L) {
    cat("\n  Surface caveats:\n")
    for (n in x$notes) {
      cat(.wrap_note_lines(n), sep = "\n")
    }
  }
  invisible(x)
}
