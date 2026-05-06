# Asymmetry diagnostics for the GRASS framework.
#
# Two functions live here:
#
#   * `check_rater_asymmetry()` — the old per-rater Se/Sp diagnostic.
#     Renamed (formerly `check_asymmetry()`). Computes mean / max
#     `|Se_j - Sp_j|` and tiers by 0.05 / 0.10. Class: `grass_asymmetry`
#     (unchanged so the existing print method, as.data.frame method,
#     and reporting-card test fixtures keep working).
#
#   * `check_asymmetry(ratings, ...)` — the new ratings-input diagnostic
#     defined by §4.3 of the v0.2.0 paper-alignment design doc. Computes
#     the cross-coefficient percentile spread `delta_hat` (in pp) and
#     tiers by NP-motivated size-alpha thresholds (9.25, 11.75) per paper §3.2 / App G.
#     Class: `grass_asymmetry_panel`.
#
# v0.5.0 ICC scope decision (2026-05-05). delta_hat is computed over the
# AGREEMENT FAMILY ONLY -- PABAK, mean AC1, Fleiss kappa, Krippendorff alpha
# -- because each has a closed-form reference depending on (q, pi_+) only and
# is therefore DGP-robust at the panel level. ICC is reported alongside on
# the panel rows but does NOT enter delta_hat: ICC's reference surface
# depends on the full F-shape variance structure, and a panel whose true F
# does not match the bundled logit-normal reference incurs surface-percentile
# drift on the order of 20 pp at small designs (sanity probe, 2026-05-05).
# The agreement family is mutually F-shape-robust, so cross-family spread is
# a clean asymmetry detector.

# Coefficients that enter delta_hat. Anything else (currently just `icc`) is
# reported on the panel with `in_delta_hat = FALSE`.
.DELTA_AGREEMENT_COEFS <- c("pabak", "mean_ac1", "fleiss_kappa",
                            "krippendorff_a")
#
# Soft-deprecation: if `check_asymmetry()` is called with the OLD signature
# (`se = ...`, `sp = ...`), we route the call to `check_rater_asymmetry()`
# and emit a one-time `msg_once()` hint. Supplying both `ratings = ...` and
# `se =`/`sp =` is an error.

# ---------------------------------------------------------------------------
# (1) Renamed: check_rater_asymmetry()
# ---------------------------------------------------------------------------

#' Per-rater Se/Sp asymmetry diagnostic
#'
#' `check_rater_asymmetry()` is the rater-level companion to
#' [check_asymmetry()]. It takes per-rater sensitivity and specificity
#' estimates, computes the per-rater gap `|Se_j - Sp_j|`, reduces them to a
#' scalar `delta_hat = mean_norm |Se_j - Sp_j|` (or `max`), and assigns the
#' three-tier model-safety regime that governs whether `q_hat` (the GRASS
#' operating-quality projection onto the Se = Sp diagonal) is trustworthy
#' as the primary summary.
#'
#' This function is appropriate when the user already has per-rater
#' Se/Sp estimates — typically from a Hui-Walter / Dawid-Skene latent-class
#' fit (see [latent_class_fit()]), a simulation with known truth, or a
#' reference-standard comparison. It is the function called inside the
#' divergent branch of the new [check_asymmetry()] / `grass_report()`
#' workflow when per-rater (Se, Sp) become available.
#'
#' For ratings-matrix input (when only the N x k binary ratings are
#' available, no per-rater Se/Sp), use [check_asymmetry()] instead — it
#' computes the cross-coefficient percentile spread without requiring
#' identifiability of per-rater (Se, Sp).
#'
#' Tier thresholds follow the three-tier architecture:
#'
#' - **Tier 1 — `ok`** (`delta_hat < 0.05`): diagonal default. Report
#'   `q_hat +/- SE` plus EDR and EMR_panel without per-rater disclosure.
#' - **Tier 2 — `caution`** (`0.05 <= delta_hat < 0.10`): diagonal +
#'   diagnostic. Report `q_hat +/- SE` with a caution flag and the
#'   `delta_hat` value.
#' - **Tier 3 — `unsafe`** (`delta_hat >= 0.10`): full latent-class.
#'   `q_hat` is withheld as primary; report per-rater `(Se_j, Sp_j)` via a
#'   Hui-Walter / Dawid-Skene fit. PABAK's prevalence-flatness was
#'   conditional on `Se = Sp`; at Tier 3 that condition is visibly broken.
#'
#' `se` and `sp` are not identifiable from a single 2-rater binary table
#' without external ground truth. Supply them from: (a) a simulation with
#' known truth, (b) a Hui-Walter / Dawid-Skene latent-class fit under
#' `k >= 3` and prevalence heterogeneity, or (c) a reference-standard
#' comparison. Passing raw `table()` off-diagonals as `se` / `sp` is
#' **wrong** and the function cannot detect the error.
#'
#' @param se Numeric vector, one per rater. Per-rater sensitivity estimates
#'   in `[0, 1]`.
#' @param sp Numeric vector, one per rater. Per-rater specificity estimates
#'   in `[0, 1]`. Same length as `se`.
#' @param rater Optional character vector of rater labels, same length as
#'   `se`. Defaults to `"R1"`, `"R2"`, ... .
#' @param threshold_caution Boundary between Tier 1 (`ok`) and Tier 2
#'   (`caution`). Default `0.05`.
#' @param threshold_unsafe Boundary between Tier 2 (`caution`) and Tier 3
#'   (`unsafe`). Default `0.10`.
#' @param summary How to reduce per-rater gaps to the scalar `delta_hat`.
#'   `"max"` (default, conservative tripwire — any single rater above the
#'   threshold triggers escalation) or `"mean"` (panel-average asymmetry).
#'   The framework uses `"max"`; `"mean"` is provided for sensitivity
#'   checks.
#'
#' @return An S3 object of class `grass_asymmetry` with fields:
#' - `per_rater`: data.frame of `rater`, `se`, `sp`, and `gap = |se - sp|`
#' - `delta_hat`: the scalar `delta_hat` summary
#' - `summary`: which summary statistic was used (`"max"` or `"mean"`)
#' - `tier`: integer tier (`1`, `2`, or `3`)
#' - `regime`: character regime label (`"ok"`, `"caution"`, or `"unsafe"`)
#' - `thresholds`: named list of the caution and unsafe cutoffs
#'
#' @seealso [check_asymmetry()] for the ratings-matrix companion (panel-level
#'   percentile spread). [latent_class_fit()] for the recommended source of
#'   `(se, sp)` estimates.
#' @export
#'
#' @examples
#' # Tier 1: symmetric raters — safe to use q_hat as primary
#' check_rater_asymmetry(se = c(0.86, 0.88, 0.84), sp = c(0.85, 0.87, 0.86))
#'
#' # Tier 2: one rater pressing the Se = Sp assumption
#' check_rater_asymmetry(se = c(0.90, 0.88, 0.92), sp = c(0.82, 0.86, 0.85))
#'
#' # Tier 3: within-rater Se-favoring regime — requires latent-class fit
#' check_rater_asymmetry(se = c(0.95, 0.93, 0.94), sp = c(0.78, 0.80, 0.79))
check_rater_asymmetry <- function(se, sp,
                                  rater = NULL,
                                  threshold_caution = 0.05,
                                  threshold_unsafe = 0.10,
                                  summary = c("max", "mean")) {
  summary <- match.arg(summary)

  if (!is.numeric(se) || !is.numeric(sp)) {
    stop("`se` and `sp` must be numeric vectors.", call. = FALSE)
  }
  if (length(se) != length(sp)) {
    stop("`se` and `sp` must have the same length (one entry per rater).",
         call. = FALSE)
  }
  if (length(se) < 1L) {
    stop("Need at least one rater's (se, sp) estimate.", call. = FALSE)
  }
  if (any(!is.finite(se)) || any(!is.finite(sp))) {
    stop("`se` and `sp` must be finite.", call. = FALSE)
  }
  if (any(se < 0 | se > 1) || any(sp < 0 | sp > 1)) {
    stop("`se` and `sp` values must lie in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(threshold_caution) || !is.numeric(threshold_unsafe) ||
      length(threshold_caution) != 1L || length(threshold_unsafe) != 1L ||
      threshold_caution <= 0 || threshold_unsafe <= threshold_caution) {
    stop("Thresholds must satisfy 0 < threshold_caution < threshold_unsafe.",
         call. = FALSE)
  }

  if (is.null(rater)) {
    rater <- paste0("R", seq_along(se))
  } else if (length(rater) != length(se)) {
    stop("`rater` must be the same length as `se` and `sp`.", call. = FALSE)
  }

  gap <- abs(se - sp)
  delta_hat <- if (summary == "max") max(gap) else mean(gap)

  tier <- if (delta_hat < threshold_caution) {
    1L
  } else if (delta_hat < threshold_unsafe) {
    2L
  } else {
    3L
  }
  regime <- c("ok", "caution", "unsafe")[tier]

  out <- list(
    per_rater = data.frame(
      rater = as.character(rater),
      se = as.numeric(se),
      sp = as.numeric(sp),
      gap = as.numeric(gap),
      stringsAsFactors = FALSE
    ),
    delta_hat = as.numeric(delta_hat),
    summary = summary,
    tier = tier,
    regime = regime,
    thresholds = list(caution = as.numeric(threshold_caution),
                      unsafe = as.numeric(threshold_unsafe))
  )
  class(out) <- c("grass_asymmetry", "list")
  out
}

#' @export
print.grass_asymmetry <- function(x, digits = 3, ...) {
  cat("grass rater asymmetry diagnostic (per-rater |Se - Sp|)\n",
      "  raters              : ", nrow(x$per_rater), "\n",
      "  per-rater gaps |Se-Sp|\n", sep = "")
  pr <- x$per_rater
  for (i in seq_len(nrow(pr))) {
    cat(sprintf("    %-8s  Se = %.*f   Sp = %.*f   gap = %.*f\n",
                pr$rater[i],
                digits, pr$se[i],
                digits, pr$sp[i],
                digits, pr$gap[i]))
  }
  cat(sprintf("  delta_hat (%s)    : %.*f\n",
              x$summary, digits, x$delta_hat))
  cat(sprintf("  thresholds           : caution = %.*f,  unsafe = %.*f\n",
              digits, x$thresholds$caution,
              digits, x$thresholds$unsafe))
  cat("  tier                 : ", x$tier,
      "  (", x$regime, ")\n", sep = "")
  tier_msg <- switch(
    x$regime,
    ok      = "Diagonal default: report q_hat plus EDR / EMR_panel.",
    caution = "Diagonal + diagnostic: report q_hat with caution flag and delta_hat.",
    unsafe  = "Full latent-class: withhold q_hat as primary; report per-rater (Se_j, Sp_j)."
  )
  cat("  reporting guidance   : ", tier_msg, "\n", sep = "")
  invisible(x)
}

#' Coerce a `grass_asymmetry` result to a one-row data.frame
#'
#' @param x A `grass_asymmetry` object.
#' @param row.names,optional,... Standard arguments; ignored except
#'   `row.names`.
#' @return A one-row data.frame with `n_raters`, `delta_hat`, `summary`,
#'   `tier`, and `regime` columns.
#' @export
as.data.frame.grass_asymmetry <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    n_raters = nrow(x$per_rater),
    delta_hat = x$delta_hat,
    summary = x$summary,
    tier = x$tier,
    regime = x$regime,
    threshold_caution = x$thresholds$caution,
    threshold_unsafe = x$thresholds$unsafe,
    row.names = row.names,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# (2) New: check_asymmetry(ratings, ...)
# ---------------------------------------------------------------------------

#' Cross-coefficient panel asymmetry diagnostic
#'
#' `check_asymmetry()` takes an N x k binary ratings matrix and returns a
#' scalar `delta_hat` (in percentile points, "pp") summarising the spread
#' across the panel of agreement coefficients (PABAK, AC1, Fleiss kappa,
#' Krippendorff alpha, ICC) when each is positioned on its calibrated
#' reference surface. A wide spread signals that the panel is reporting
#' inconsistent stories about the same data — a structural symptom of
#' shared rater bias or other Se != Sp asymmetry — and routes the user to
#' the per-rater latent-class diagnostic.
#'
#' Each coefficient is positioned on its reference surface via
#' [position_on_surface()], yielding a percentile in `[0, 100]` pp. The
#' panel diagnostic is `delta_hat = max(percentile) - min(percentile)`,
#' computed over coefficients whose observed value sits within the
#' achievable range of their reference surface (see *Surface-envelope
#' clamp* below).
#'
#' @section Surface-envelope clamp (v0.2.1+):
#' If an observed coefficient value falls outside the achievable range of
#' its reference surface at the study's design `(pi_hat, k, N)`, the
#' inversion to `q_hat` clamps to the boundary and the percentile lands
#' at exactly `0` or `100`. Including such clamped percentiles in the
#' max-min `delta_hat` would inflate the panel-spread purely because of
#' the clamp, not because the panel disagrees on percentile. Since
#' v0.2.1 the function therefore *excludes* clamped coefficients from
#' `delta_hat` whenever at least two unclamped coefficients remain. The
#' affected coefficients are still shown in the returned `panel`
#' data.frame with `clamped = TRUE` and surface percentile `0` or `100`,
#' and a note in `$notes` names which coefficients were excluded. This
#' matters most often for ICC at `N > 200` (the bundled fitted-ICC
#' reference is calibrated through `N = 200`; beyond that the function
#' falls back to the oracle reference, which over-predicts the
#' practitioner's GLMM-fit ICC and therefore clamps). The behavior
#' guarantees that a `divergent` flag always reflects genuine
#' cross-coefficient surface-percentile disagreement rather than an
#' envelope artifact; the trade-off is that a coefficient with no
#' calibrated reference at the user's design contributes nothing to
#' the diagnostic. If fewer than two unclamped coefficients remain,
#' `delta_hat` falls back to the raw spread including clamped values
#' and the note records the fallback.
#'
#' Tier thresholds default to `c(9.25, 11.75)` pp, the Neyman-Pearson-motivated
#' size-alpha cutoffs at `alpha = 0.05` and `alpha = 0.01` from §3.2 of the
#' GRASS paper (see also App G operating characteristics). The construction is
#' NP-optimal within the test class of threshold rules on `delta_hat`; it is
#' not a likelihood-ratio test on the joint distribution of the panel
#' coefficients (`delta_hat` is a paired-margin difference, not a likelihood
#' ratio). The three flags are:
#'
#' - **`aligned`** (`delta_hat < 9.25` pp): the panel agrees on the
#'   surface position. Any single coefficient is a stable summary; the
#'   primary coefficient (Table 2) carries the headline.
#' - **`caution`** (`9.25 <= delta_hat < 11.75` pp): the panel is mildly
#'   inconsistent. Report the primary coefficient with a caution flag and
#'   the `delta_hat` value.
#' - **`divergent`** (`delta_hat >= 11.75` pp): no single coefficient is a
#'   stable summary. Use [latent_class_fit()] to recover per-rater
#'   `(Se, Sp)` and report those instead.
#'
#' The new `check_asymmetry(ratings, ...)` signature replaces the
#' v0.1.x `check_asymmetry(se, sp, ...)` per-rater signature. For
#' backward compatibility, calling `check_asymmetry()` with `se = ...`
#' and `sp = ...` named arguments emits a one-time deprecation hint and
#' routes the call to [check_rater_asymmetry()]. Supplying both
#' `ratings = ...` and per-rater `se =` / `sp =` is an error.
#'
#' @param ratings User input: an `N x k` binary matrix, an `N` x `k`
#'   data.frame whose columns are 0/1 / logical / 2-level factor, or a
#'   list of two equal-length 0/1 vectors (k = 2 paired form). See
#'   `?normalize_ratings` for accepted shapes.
#' @param axis `"inter"` (default) or `"intra"`. Selects the surface
#'   family.
#' @param occasion Reserved for `axis = "intra"` (a vector / factor
#'   identifying viewing occasion); ignored when `axis = "inter"`.
#' @param delta_thresholds Length-2 numeric vector `c(caution, divergent)`
#'   in percentile points. Default `c(9.25, 11.75)` from paper §3.2 (NP
#'   optimum at alpha = 0.05 / 0.01).
#' @param ... Forwarded to [position_on_surface()] (e.g. `bands`,
#'   `band_labels`, `reference_type`).
#'
#' @return An S3 object of class `grass_asymmetry_panel` with fields:
#' - `delta_hat`: scalar percentile-spread in pp
#' - `flag`: one of `"aligned"`, `"caution"`, `"divergent"`
#' - `thresholds`: named numeric vector of the (caution, divergent)
#'   cutoffs, in pp
#' - `panel`: data.frame with `coefficient`, `observed`, `percentile_pp`
#' - `notes`: character vector of unique caveats from the underlying
#'   surface positioning calls (e.g. nearest-neighbor gaps, ICC
#'   unavailability)
#'
#' @seealso [check_rater_asymmetry()] for the per-rater Se/Sp companion;
#'   [latent_class_fit()] for the divergent-branch recovery of per-rater
#'   `(Se, Sp)`; [position_on_surface()] for the underlying surface
#'   positioning.
#'
#' @references
#' Neyman, J. and Pearson, E. S. (1933). On the problem of the most
#' efficient tests of statistical hypotheses. *Philosophical Transactions
#' of the Royal Society of London A*, 231, 289-337.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' # Build a 5x200 symmetric panel — should print as 'aligned'.
#' Y <- matrix(rbinom(5 * 200, 1, 0.30), nrow = 200, ncol = 5)
#' check_asymmetry(Y)
#' }
check_asymmetry <- function(ratings,
                            axis = c("inter", "intra"),
                            occasion = NULL,
                            delta_thresholds = NULL,
                            ...) {
  # ---- Soft-deprecation dispatch for OLD `check_asymmetry(se, sp, ...)` ---
  # Detection: the OLD signature was `check_asymmetry(se, sp, ...)` where
  # `se` and `sp` are numeric vectors of length k (per-rater Se/Sp). The
  # NEW signature takes `ratings` as an N x k matrix / data.frame / list.
  # So a call routes to the OLD path if any of:
  #   (a) the user supplied `se =` and / or `sp =` as named args, or
  #   (b) the first positional arg (now bound to `ratings`) is a numeric
  #       vector that is NOT a matrix / data.frame / list — i.e. the OLD
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
           "Se/Sp diagnostic use `check_rater_asymmetry(se, sp, ...)`.",
           call. = FALSE)
    }
    msg_once(
      "deprecate_check_asymmetry_se_sp",
      paste0("grass: `check_asymmetry(se = ..., sp = ...)` is deprecated. ",
             "Use `check_rater_asymmetry(se, sp, ...)` for the per-rater ",
             "diagnostic, or `check_asymmetry(ratings = Y)` for the new ",
             "panel-level percentile-spread diagnostic. Routing this call ",
             "to `check_rater_asymmetry()`.")
    )

    # Forward to check_rater_asymmetry. Map both named-arg forms and the
    # positional form. In the positional form: `ratings` actually held the
    # per-rater Se vector and `axis` held the per-rater Sp vector (because
    # those were the OLD positional slots).
    if (has_se || has_sp) {
      # Named se / sp form: forward dots as-is.
      old_args <- dots
      # The OLD signature also had `summary` and `rater` and threshold args
      # which may live in dots. Pass them through.
      return(do.call(check_rater_asymmetry, old_args))
    } else {
      # Positional form: ratings = se vector, axis = sp vector.
      # `axis` may still be the default `c("inter","intra")` if user only
      # supplied one positional arg. In that case error out — a single
      # numeric vector is ambiguous.
      if (missing(axis) || identical(axis, c("inter", "intra"))) {
        stop("check_asymmetry(): single numeric vector input is ambiguous. ",
             "Did you mean `check_rater_asymmetry(se, sp, ...)` (per-rater ",
             "diagnostic) or `check_asymmetry(ratings = Y)` (ratings ",
             "matrix)?", call. = FALSE)
      }
      old_args <- c(list(se = ratings, sp = axis), dots)
      return(do.call(check_rater_asymmetry, old_args))
    }
  }

  # ---- Normal path: ratings-input panel diagnostic -----------------------
  axis <- match.arg(axis)
  Y <- normalize_ratings(ratings)
  panel_obs <- compute_panel(Y, axis = axis, occasion = occasion)

  # Resolve delta_thresholds. NULL (default) -> per-(k, N) lookup against
  # the calibration grid; non-NULL -> honor user-supplied pair after the
  # usual validation. Lookup result includes a source tag so callers can
  # surface the choice transparently.
  thresholds_source <- "user_supplied"
  thresholds_note   <- ""
  if (is.null(delta_thresholds)) {
    lk <- lookup_delta_thresholds(k = ncol(Y), N = nrow(Y))
    delta_thresholds  <- lk$thresholds
    thresholds_source <- lk$source
    thresholds_note   <- lk$note
  }
  if (length(delta_thresholds) != 2L ||
      !is.numeric(delta_thresholds) ||
      !all(is.finite(delta_thresholds)) ||
      delta_thresholds[1L] <= 0 ||
      delta_thresholds[2L] <= delta_thresholds[1L]) {
    stop("`delta_thresholds` must be a length-2 numeric vector with ",
         "0 < caution < divergent (in pp).", call. = FALSE)
  }

  # The panel-name -> surface-metric mapping. compute_panel() returns
  # `ac1` (panel-friendly short name) but position_on_surface() expects
  # `mean_ac1`; Cohen's `kappa` (k=2 only) has no surface and is dropped.
  surface_metric_for <- c(
    pabak          = "pabak",
    ac1            = "mean_ac1",
    fleiss_kappa   = "fleiss_kappa",
    krippendorff_a = "krippendorff_a",
    icc            = "icc"
  )
  panel_keep <- intersect(names(panel_obs), names(surface_metric_for))
  panel_obs <- panel_obs[panel_keep]
  # Re-key so the printed panel uses the surface metric names (paper labels).
  names(panel_obs) <- unname(surface_metric_for[panel_keep])

  positions <- lapply(names(panel_obs), function(m) {
    position_on_surface(ratings = Y, metric = m, ...)
  })
  names(positions) <- names(panel_obs)

  pcts_pp <- vapply(
    positions,
    function(p) p$percentile * 100,
    numeric(1)
  )

  # delta_hat measures cross-coefficient surface-percentile spread over the
  # AGREEMENT FAMILY only (PABAK, mean AC1, Fleiss kappa, Krippendorff alpha).
  # ICC is reported on the panel but excluded from delta_hat: see the v0.5.0
  # scope note at the top of this file. Two further exclusions inside the
  # agreement family:
  #   1. Surface-envelope clamps: a coefficient whose obs_value falls outside
  #      the achievable range of its reference surface gets q_hat clamped
  #      (percentile lands at 0 or 100). Including those would inflate
  #      delta_hat purely from the clamp. Excluded when >= 2 unclamped remain.
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
    delta_hat <- max(pcts_pp[in_delta_set]) - min(pcts_pp[in_delta_set])
    # Note clamped-but-otherwise-eligible agreement coefficients
    excluded_for_clamp <- names(positions)[
      names(positions) %in% .DELTA_AGREEMENT_COEFS & clamped
    ]
    if (length(excluded_for_clamp)) {
      clamp_note <- sprintf(
        "delta_hat over %d agreement-family coefficients; excluded due to surface-envelope clamp: %s. (ICC is always reported but never enters delta_hat per the v0.5.0 scope note; its surface percentile depends on the full F-shape and is reported separately.)",
        sum(in_delta_set),
        paste(excluded_for_clamp, collapse = ", ")
      )
    } else if ("icc" %in% names(positions)) {
      # No clamps; just remind that ICC isn't in delta_hat.
      clamp_note <- "delta_hat is computed over the agreement family (PABAK, AC1, Fleiss, alpha). ICC is reported on the panel but does not enter delta_hat (v0.5.0 scope: ICC's reference surface is F-shape-dependent and does not share the (q, pi_+) sufficient statistic the agreement family does)."
    }
  } else {
    # Fall back to all available agreement-family percentiles, even if
    # clamped, so the spread is still defined.
    agreement_idx <- names(positions) %in% .DELTA_AGREEMENT_COEFS
    if (sum(agreement_idx) >= 2L) {
      delta_hat <- max(pcts_pp[agreement_idx]) - min(pcts_pp[agreement_idx])
      clamp_note <- sprintf(
        "Fewer than 2 unclamped agreement-family coefficients; delta_hat uses raw spread over %d agreement coefficients including clamped: %s.",
        sum(agreement_idx),
        paste(names(positions)[agreement_idx], collapse = ", ")
      )
    } else {
      delta_hat <- NA_real_
      clamp_note <- "delta_hat undefined: fewer than 2 agreement-family coefficients available."
    }
  }

  flag <- if (delta_hat < delta_thresholds[1L]) {
    "aligned"
  } else if (delta_hat < delta_thresholds[2L]) {
    "caution"
  } else {
    "divergent"
  }

  combined_notes <- unique(c(
    unlist(lapply(positions, `[[`, "notes")),
    clamp_note,
    if (nzchar(thresholds_note)) thresholds_note else character(0L)
  ))

  out <- list(
    delta_hat  = unname(delta_hat),
    flag       = flag,
    thresholds = setNames(as.numeric(delta_thresholds),
                          c("caution", "divergent")),
    thresholds_source = thresholds_source,
    panel      = data.frame(
      coefficient   = names(panel_obs),
      observed      = unlist(panel_obs, use.names = FALSE),
      percentile_pp = unname(pcts_pp),
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
  cat(sprintf("  delta_hat = %.*f pp\n",
              digits, x$delta_hat))
  cat(sprintf("  flag      = %s  (%s %.*f)\n",
              x$flag,
              if (x$flag == "aligned") "<" else ">=",
              digits,
              if (x$flag == "aligned") x$thresholds[["caution"]]
              else if (x$flag == "caution") x$thresholds[["caution"]]
              else x$thresholds[["divergent"]]))
  cat("\n  panel:\n")
  cat("    coefficient        observed   percentile (pp)   in delta_hat\n")
  pn <- x$panel
  for (i in seq_len(nrow(pn))) {
    in_d <- if (isTRUE(pn$in_delta_hat[i])) "yes"
            else if (identical(pn$coefficient[i], "icc"))
              "no [F-shape sensitive]"
            else "no"
    cat(sprintf("    %-18s %.2f       %.*f             %s\n",
                pn$coefficient[i],
                pn$observed[i],
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
      cat("    - ", n, "\n", sep = "")
    }
  }
  invisible(x)
}
