# Asymmetry diagnostics for the GRASS framework.
#
# Two functions live here:
#
#   * `check_rater_asymmetry()` -- the old per-rater Se/Sp diagnostic.
#     Renamed (formerly `check_asymmetry()`). Computes mean / max
#     `|Se_j - Sp_j|` and tiers by 0.05 / 0.10. Class: `grass_asymmetry`
#     (unchanged so the existing print method, as.data.frame method,
#     and reporting-card test fixtures keep working).
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
#' Se/Sp estimates -- typically from a Hui-Walter / Dawid-Skene latent-class
#' fit (see [latent_class_fit()]), a simulation with known truth, or a
#' reference-standard comparison. It is the function called inside the
#' divergent branch of the new [check_asymmetry()] / `grass_report()`
#' workflow when per-rater (Se, Sp) become available.
#'
#' For ratings-matrix input (when only the N x k binary ratings are
#' available, no per-rater Se/Sp), use [check_asymmetry()] instead -- it
#' computes the cross-coefficient percentile spread without requiring
#' identifiability of per-rater (Se, Sp).
#'
#' Tier thresholds follow the three-tier architecture:
#'
#' - **Tier 1 -- `ok`** (`delta_hat < 0.05`): diagonal default. Report
#'   `q_hat +/- SE` plus EDR and EMR_panel without per-rater disclosure.
#' - **Tier 2 -- `caution`** (`0.05 <= delta_hat < 0.10`): diagonal +
#'   diagnostic. Report `q_hat +/- SE` with a caution flag and the
#'   `delta_hat` value.
#' - **Tier 3 -- `unsafe`** (`delta_hat >= 0.10`): full latent-class.
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
#'   `"max"` (default, conservative tripwire -- any single rater above the
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
#' # Tier 1: symmetric raters -- safe to use q_hat as primary
#' check_rater_asymmetry(se = c(0.86, 0.88, 0.84), sp = c(0.85, 0.87, 0.86))
#'
#' # Tier 2: one rater pressing the Se = Sp assumption
#' check_rater_asymmetry(se = c(0.90, 0.88, 0.92), sp = c(0.82, 0.86, 0.85))
#'
#' # Tier 3: within-rater Se-favoring regime -- requires latent-class fit
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
#' scalar `delta_hat` (in quality percentage points, "pp"): the max-min
#' spread of the *implied panel qualities* across the three agreement
#' coefficients (PABAK, mean AC1, Fleiss kappa). Each coefficient inverts
#' to its own `q_hat` on the shared (q, pi_+) reference; if the
#' calibration DGP held exactly, all three would imply the same quality,
#' so the spread measures cross-coefficient model discordance in
#' interpretable units of quality. (Option B, ratified 2026-07-05; the
#' previous definition -- spread of surface percentiles across four
#' coefficients including Krippendorff alpha -- ran through the retired
#' nearest-cell percentile machinery, whose sawtooth inflated `delta_hat`
#' with quantization noise. Alpha left the panel at 0.6.0; ICC never
#' enters `delta_hat`, its reference being distribution-sensitive in ways
#' the agreement family is not.)
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
#' @param delta_thresholds Deprecated (grassr 0.7.0). The per-(k, N)
#'   size-alpha threshold table is retired; the flag is now `delta_hat`'s
#'   percentile on the matched null distribution. Supplying a length-2
#'   `c(caution, divergent)` pair (in pp) is still honored for a single
#'   call with the legacy pp-cut semantics, but raises a deprecation
#'   warning. Default `NULL` (use the matched-null percentile convention).
#' @param ... Forwarded to [position_on_surface()] (e.g. `reference_type`).
#'
#' @return An S3 object of class `grass_asymmetry_panel` with fields:
#' - `delta_hat`: scalar implied-quality spread, in pp of quality
#' - `delta_percentile`: `delta_hat`'s percentile on the matched
#'   (k, N, q_hat) null ECDF (`NA` if the null is uncalibrated at the
#'   design)
#' - `flag`: one of `"aligned"`, `"caution"`, `"divergent"`
#' - `matched_null`: list describing the matched null cell
#'   (`k`, `N`, `q`, `q_hat_panel`, `n_draws`, `snapped`,
#'   `unstable_tail`), or `NULL` if uncalibrated
#' - `thresholds`: named numeric vector of the implied (caution, divergent)
#'   pp cuts (95th/99th of the matched null), or the legacy user-supplied
#'   cuts
#' - `thresholds_source`: one of `"matched_null_ecdf"`,
#'   `"user_supplied_legacy"`, `"not_calibrated"`
#' - `panel`: data.frame with `coefficient`, `observed`, `implied_q`,
#'   `percentile_pp` (pooled percentile), `clamped`, `in_delta_hat`
#' - `notes`: character vector of unique caveats from the underlying
#'   surface positioning calls (e.g. nearest-neighbor gaps, ICC
#'   unavailability, matched-null provenance)
#'
#' @seealso [check_rater_asymmetry()] for the per-rater Se/Sp companion;
#'   [latent_class_fit()] for the divergent-branch recovery of per-rater
#'   `(Se, Sp)`; [position_on_surface()] for the underlying surface
#'   positioning.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' # Build a 5x200 symmetric panel -- should print as 'aligned'.
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
      # supplied one positional arg. In that case error out -- a single
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

  # v0.7.0: the per-(k, N) threshold table is retired. The flag comes from
  # delta_hat's percentile on the matched (k, N, q_hat) null ECDF
  # (>= 95th caution, >= 99th divergent), resolved after the panel's
  # q_hat is known. A user-supplied `delta_thresholds` pair is honored
  # with the legacy pp-cut semantics, with a deprecation warning.
  legacy_thresholds <- NULL
  thresholds_note   <- ""
  if (!is.null(delta_thresholds)) {
    if (length(delta_thresholds) != 2L ||
        !is.numeric(delta_thresholds) ||
        !all(is.finite(delta_thresholds)) ||
        delta_thresholds[1L] <= 0 ||
        delta_thresholds[2L] <= delta_thresholds[1L]) {
      stop("`delta_thresholds` must be a length-2 numeric vector with ",
           "0 < caution < divergent (in pp).", call. = FALSE)
    }
    warning("`delta_thresholds` is deprecated as of grassr 0.7.0: the flag ",
            "is now delta_hat's percentile on the matched null ",
            "distribution. The supplied pp cuts are honored for this call.",
            call. = FALSE)
    legacy_thresholds <- delta_thresholds
  }

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
        "delta_hat (implied-quality spread) over %d agreement-family coefficients; excluded due to reference-envelope clamp: %s. (ICC is always reported but never enters delta_hat per the v0.5.0 scope note; its reference depends on the full subject-prevalence distribution F and is reported separately.)",
        sum(in_delta_set),
        paste(excluded_for_clamp, collapse = ", ")
      )
    } else if ("icc" %in% names(positions)) {
      # No clamps; just remind that ICC isn't in delta_hat.
      clamp_note <- "delta_hat is the implied-quality spread over the agreement family (PABAK, mean AC1, Fleiss kappa). ICC is reported on the panel but does not enter delta_hat (v0.5.0 scope: ICC's reference depends on the full subject-prevalence distribution F and does not share the (q, pi_+) sufficient statistic the agreement family does)."
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
  null_cell <- if (!k2_degenerate && is.finite(q_hat_panel))
    lookup_delta_null(k = ncol(Y), N = nrow(Y), q_hat = q_hat_panel)
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
                         q_hat_panel = unname(q_hat_panel),
                         n_draws = null_cell$n_draws,
                         snapped = null_cell$snapped,
                         unstable_tail = null_cell$unstable_tail)
    thresholds_note <- sprintf(
      "flag from delta_hat's percentile on the matched null (k=%d, N=%d, q=%.2f; %s draws)%s%s.",
      null_cell$k, null_cell$N, null_cell$q,
      format(null_cell$n_draws, big.mark = ","),
      if (null_cell$snapped) "; design snapped to nearest calibrated cell" else "",
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

  flag <- if (!is.null(legacy_thresholds)) {
    if (delta_hat < legacy_thresholds[1L]) "aligned"
    else if (delta_hat < legacy_thresholds[2L]) "caution"
    else "divergent"
  } else if (k2_degenerate) {
    "not_applicable"
  } else {
    delta_flag_from_percentile(delta_percentile,
      if (!is.null(null_cell)) null_cell$conventions
      else c(caution = 0.95, divergent = 0.99))
  }
  thresholds_source <- if (!is.null(legacy_thresholds)) "user_supplied_legacy"
                       else if (k2_degenerate) "not_applicable_k2"
                       else if (!is.null(null_cell)) "matched_null_ecdf"
                       else "not_calibrated"
  report_cuts <- if (!is.null(legacy_thresholds))
    setNames(as.numeric(legacy_thresholds), c("caution", "divergent"))
  else implied_cuts

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
      cat("    - ", n, "\n", sep = "")
    }
  }
  invisible(x)
}
