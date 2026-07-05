# Target-2 surface-position reporting convention (v0.7.1 sweep redesign).
# Takes an observed agreement coefficient together with the design context
# (pi_hat, k, N) and positions it against the DGP-calibrated reference
# surface at that design. Returns three read-outs of one sweep object: the
# pooled percentile (position within the design's achievable agreement
# range), the 95% test-inversion consistency band on panel quality, and
# the full p(q) sweep profile. See design/v0.7.1_position_redesign.md for
# the ratified spec this function operationalises.

#' Position an observed agreement coefficient on its DGP-calibrated surface
#'
#' `position_on_surface()` is the Target-2 reporting primitive for the
#' merged GRASS binary-rater-reliability paper. Given an observed coefficient
#' value and the study design `(pi_hat, k, N)`, it inverts the coefficient
#' to an implied panel quality `q_hat` (rater operating quality on the
#' Se = Sp diagonal under the clustered latent-class DGP) and evaluates the
#' observed value against EVERY calibrated quality level at the matched
#' design -- a sweep -- from which it derives three read-outs: the pooled
#' percentile, the consistency band on quality, and the p(q) sweep profile.
#'
#' The function implements the v0.7.1 sweep convention (ratified
#' 2026-07-05): the practitioner cites the observed coefficient, its pooled
#' percentile (position within the design's achievable range), and the
#' consistency band on panel quality. `q_hat` is promoted to the card via
#' the consistency band; it also carries the surface parameterization and
#' delta-method SE. The stipulated four-band adjective
#' (Poor/Moderate/Strong/Excellent) and the modal-band confidence
#' qualifier (decisive/moderate/weak) are retired; the `bands` /
#' `band_labels` arguments are deprecated (warned, ignored).
#'
#' @section The three read-outs of the sweep:
#' At the matched design `(F/pi_hat, k, N)` the observed coefficient is
#' evaluated against every calibrated quality level q with no cell
#' selection and no q snapping. Write
#' `p(q) = P(coefficient <= obs_value | panel of quality q, this design)`.
#'
#' - **`percentile`** -- the pooled percentile: a trapezoid-weighted
#'   average of `p(q)` over the calibrated quality axis (trapezoid because
#'   the grid is non-uniform). It reads as the observed coefficient's
#'   position within the design's full achievable agreement range and is
#'   monotone in `obs_value` by construction. Returned in `[0, 1]`;
#'   callers such as `grass_report()` print it on the `[0, 100]` scale.
#'   This replaces the retired nearest-q_hat-cell percentile, whose cohort
#'   selection by a statistic derived from the coefficient made it a
#'   non-monotone sawtooth (panel review 2026-07-05).
#' - **`band`** -- the 95% test-inversion consistency band on quality: the
#'   quality levels q with `0.025 <= p(q) <= 0.975`, endpoints
#'   interpolated where `p(q)` crosses 0.975 (lower) and 0.025 (upper).
#'   Open-ended at the grid boundary is reported with a boundary flag.
#' - **`sweep`** -- the full `data.frame(q, p)` profile, the object the
#'   sweep-ridgeline graphic renders.
#'
#' @section Sweep construction:
#' Two methods are implemented.
#'
#' **Empirical method** (`method = "empirical"`, default): at the matched
#' `(F_key, k, N)` cell of the bundled `empirical_q_hat_surface`, ranks
#' `q_hat` within each calibrated quality cell's `q_hat_rep` distribution
#' (monotone-equivalent to ranking `obs_value` within the cell's
#' coefficient distribution). The full per-rep data is not bundled
#' (~300 MB); the package ships a precomputed multi-point empirical-quantile
#' summary per cell. The whole quality axis is consulted -- there is
#' deliberately no q selection. When the design `(pi_hat, k, N)` falls
#' outside the simulated grid, nearest-neighbour clamping is applied and
#' flagged in `notes`.
#'
#' **Delta method** (`method = "delta"`): the summary-stats-only fallback
#' (also used for ICC, whose reference curve carries the F-shape
#' conditioning). At each swept q it approximates the sampling distribution
#' as `Normal(E[metric](q), sd_metric(q))` and evaluates
#' `p(q) = pnorm(obs_value; mean, sd)` on the calibrated q axis. A
#' caller-supplied `surface_data$per_rep` single-cohort vector is still
#' honored for reproducibility audits, yielding a plain cohort percentile
#' with no sweep or band.
#'
#' @section Internal reference-surface arithmetic:
#' Under the clustered latent-class DGP with symmetric raters (`Se = Sp = q`),
#' the large-N closed forms for PABAK, Fleiss kappa, AC1, and Krippendorff's
#' alpha depend on `q` and the marginal positive rate `pi_+` only. This
#' function uses `pi_hat` as a plug-in for `pi_+` and inverts the observed
#' value on a 501-point q-grid on `[0.5, 1]` (matching
#' `paper2/code/12_q_inversion.R` resolution). ICC requires the full
#' subject-prevalence distribution F; in the absence of `surface_data` containing an ICC lookup, ICC
#' requests fall through to a warning-noted delta-method approximation using
#' the caller-supplied `q_hat_override` / `se_q_hat_override` if present, or
#' stop with a clear message.
#'
#' @section Ratings-primary path:
#' From v0.2.0, the preferred entry point is to hand the rating matrix
#' directly: `position_on_surface(ratings = Y, metric = "pabak")`. When
#' `ratings` is supplied, the function auto-derives `obs_value` (via
#' `compute_observed(metric, Y)`), `pi_hat` (`mean(Y)`), `k` (`ncol(Y)`),
#' and `N` (`nrow(Y)`); any of those four arguments still supplied by the
#' caller wins. This collapses the audit-style scalar-input path used in
#' v0.1.x to a single matrix argument while keeping the scalar path callable
#' for reproducibility checks. `ratings` accepts an `N x k` integer matrix in
#' `{0, 1}`, a data.frame with `k` rater columns, or a length-2 list of
#' equal-length 0/1 vectors (k = 2). Round-trip equality with the scalar
#' path is a tested invariant.
#'
#' @param obs_value Numeric scalar. The observed agreement coefficient. Optional
#'   when `ratings` is supplied (auto-derived via
#'   `compute_observed(metric, Y)`).
#' @param metric Character scalar. One of `"pabak"`, `"fleiss_kappa"`,
#'   `"mean_ac1"`, `"krippendorff_a"`, `"icc"`.
#' @param pi_hat Numeric scalar in `(0, 1)`. The panel-identified marginal
#'   positive rate. Optional when `ratings` is supplied (auto-derived via
#'   `mean(Y)`); otherwise estimate from the rating matrix via
#'   `grass_prevalence()` or directly from rater marginals.
#' @param k Integer >= 2. Number of raters. Optional when `ratings` is
#'   supplied (auto-derived as `ncol(Y)`).
#' @param N Integer >= 1. Number of subjects. Optional when `ratings` is
#'   supplied (auto-derived as `nrow(Y)`).
#' @param method One of `"empirical"` (default; uses the bundled sim-derived
#'   empirical q_hat sampling distribution) or `"delta"` (closed-form normal
#'   approximation from delta-method SE).
#' @param bands Deprecated (grassr 0.7.1). The stipulated four-band
#'   partition of `q` is retired; non-default values draw a warning and
#'   are ignored.
#' @param band_labels Deprecated (grassr 0.7.1). See `bands`.
#' @param reference_type For `metric = "icc"` only. One of `"fitted"`
#'   (default; GLMM-gap-corrected reference matching what practitioners
#'   compute via `glmer`) or `"oracle"` (closed-form
#'   `sigma^2_subject / (sigma^2_subject + pi^2/3)` with `sigma^2_subject`
#'   known from F). Use `"oracle"` only if `obs_value` was computed via
#'   oracle variance decomposition (non-standard for applied work).
#'   For N beyond the fitted-reference sim range (currently N > 200), the
#'   function auto-falls-back to oracle with an explanatory note.
#' @param ratings Optional. From v0.2.0 this is the **primary input for all
#'   metrics** (not just ICC): supplying an `N x k` rating matrix auto-derives
#'   `obs_value`, `pi_hat`, `k`, and `N`. Accepts an `N x k` integer matrix
#'   of 0/1 values (rows = subjects, cols = raters), a data.frame with
#'   `k` rater columns, or a length-2 list of equal-length 0/1 vectors
#'   (`k = 2`). For `metric = "icc"`, supplying `ratings` (also accepted as
#'   a long data.frame with columns `subject` and `rating`) additionally
#'   enables a `glmer` fit for `(mu, tau2)` that pins down the correct
#'   `F_key` for ICC inversion; without `ratings`, ICC falls back to a
#'   nearest-M1 `F_key` lookup with a prominent caveat note (tau2 is
#'   unidentified from `pi_hat` alone). The `glmer` path requires `lme4`
#'   (Suggests).
#' @param surface_data Optional. A list with one or more of the following
#'   components, used when `method = "empirical"`:
#'   - `per_rep`: a vector of per-rep metric values at the caller's own
#'     `(q, pi_hat, k, N)` cell -- an empirical sampling distribution at
#'     that design. Honored for reproducibility audits; yields a plain
#'     cohort percentile with no sweep or consistency band.
#'   - `q_grid_per_rep`, `q_grid`: legacy per-q-grid empirical inputs.
#'     Retained for backward compatibility but no longer consumed
#'     internally (the sweep convention consults the bundled per-cell
#'     quantile surface); supplying them draws a note.
#'   When `surface_data` is `NULL` and `method = "empirical"`, the function
#'   uses the bundled empirical q_hat surface, falling back to the
#'   delta-method sweep when that surface is unavailable at the design.
#' @param ... Reserved for future extension.
#'
#' @return A list of class `grass_surface_position` with fields:
#' - `observed_value` -- echo of `obs_value`
#' - `metric` -- echo of `metric`
#' - `design` -- `list(pi_hat, k, N)`
#' - `q_hat` -- implied panel quality (coefficient inverted on the
#'   reference curve); the point estimate the consistency band surrounds
#' - `se_q_hat` -- delta-method SE of `q_hat`
#' - `percentile` -- POOLED percentile in `[0, 1]`: the observed
#'   coefficient's position within the design's full achievable range
#'   (trapezoid-weighted mixture over every calibrated quality level).
#'   Monotone in `obs_value` by construction.
#' - `percentile_basis` -- provenance string for `percentile`
#' - `band` -- 95% test-inversion consistency band on quality:
#'   `list(lo, hi, level, open_low, open_high, note)`. The quality levels
#'   whose sampling distributions are consistent with the observed value
#'   at this design.
#' - `sweep` -- `data.frame(q, p)`: the full profile
#'   `p(q) = P(coefficient <= obs_value | quality q, this design)`
#' - `sampling_method` -- which method was used
#' - `reference_used` -- which reference produced the curve
#' - `notes` -- character vector of caveats (e.g. nearest-neighbor gaps)
#'
#' @seealso [check_asymmetry()] for the companion Column A tier (rater
#'   asymmetry model-safety).
#' @export
#'
#' @examples
#' # Ratings-primary path: just hand it the matrix.
#' set.seed(1)
#' Y <- matrix(rbinom(1000, 1, 0.3), nrow = 200, ncol = 5)
#' position_on_surface(ratings = Y, metric = "pabak")
#'
#' # Equivalent scalar-input path (audit):
#' position_on_surface(
#'   obs_value = 2 * mean(Y[, 1] == Y[, 2]) - 1,  # PABAK on first pair
#'   metric = "pabak", pi_hat = mean(Y), k = ncol(Y), N = nrow(Y)
#' )
#'
#' # Scalar path -- the three read-outs of the sweep convention.
#' r <- position_on_surface(
#'   obs_value = 0.62,
#'   metric    = "pabak",
#'   pi_hat    = 0.42,
#'   k         = 5,
#'   N         = 50
#' )
#' r$percentile   # pooled percentile of the achievable range
#' r$band         # consistency band on panel quality
#' head(r$sweep)  # the full p(q) profile
#'
#' # Fleiss kappa at imbalanced prevalence.
#' position_on_surface(
#'   obs_value = 0.18,
#'   metric    = "fleiss_kappa",
#'   pi_hat    = 0.08,
#'   k         = 3,
#'   N         = 200
#' )
position_on_surface <- function(obs_value = NULL,
                                metric,
                                pi_hat = NULL,
                                k = NULL,
                                N = NULL,
                                method = c("empirical", "delta"),
                                bands = c(0.5, 0.625, 0.75, 0.875, 1.0),
                                band_labels = c("Poor", "Moderate",
                                                "Strong", "Excellent"),
                                surface_data = NULL,
                                ratings = NULL,
                                reference_type = c("fitted", "oracle"),
                                ...) {
  method <- match.arg(method)
  reference_type <- match.arg(reference_type)

  # ---- Ratings-primary path (v0.2.0) -------------------------------------
  # When `ratings` is supplied AND any of `obs_value`/`pi_hat`/`k`/`N` is
  # NULL, auto-derive the missing scalars from the N x k rating matrix.
  # When the caller already supplied all four scalars, leave `ratings`
  # untouched: this preserves the legacy k x N ICC pass-through (callers
  # who hand a k x N matrix alongside scalars get exactly v0.1.x behavior).
  needs_autoderive <- is.null(obs_value) || is.null(pi_hat) ||
    is.null(k) || is.null(N)
  if (!is.null(ratings) && needs_autoderive) {
    # Detect whether `ratings` is a long data.frame (subject/rating columns)
    # bound for the ICC `glmer` fit only -- that form can't be normalized to
    # an N x k matrix here. In that case we don't auto-derive; the caller
    # must still provide the missing scalars.
    is_long_df <- is.data.frame(ratings) &&
      all(c("subject", "rating") %in% names(ratings)) &&
      !all(vapply(ratings, function(col) all(col %in% c(0, 1, NA, TRUE, FALSE)),
                  logical(1)))
    if (!is_long_df) {
      # Validate `metric` early so `compute_observed()` doesn't error first.
      if (!is.character(metric) || length(metric) != 1L) {
        stop("`metric` must be a single string.", call. = FALSE)
      }
      Y <- normalize_ratings(ratings)
      if (is.null(obs_value)) obs_value <- compute_observed(metric, Y)
      if (is.null(pi_hat))    pi_hat    <- mean(Y)
      if (is.null(k))         k         <- ncol(Y)
      if (is.null(N))         N         <- nrow(Y)
      # The downstream ICC `glmer` machinery (fit_tau2_from_ratings) expects
      # `ratings` as a `k x N` matrix (rows = raters). The ratings-primary
      # convention is `N x k`. Transpose so existing ICC code keeps working.
      if (metric == "icc") {
        ratings <- t(Y)
      }
    }
  } else if (is.null(ratings) && needs_autoderive) {
    stop("Either supply `ratings = <N x k matrix>` or all of ",
         "`obs_value`, `pi_hat`, `k`, `N` as scalars.", call. = FALSE)
  }

  # ---- Input validation --------------------------------------------------
  allowed_metrics <- c("pabak", "fleiss_kappa", "mean_ac1",
                       "krippendorff_a", "icc")
  if (!is.character(metric) || length(metric) != 1L ||
      !metric %in% allowed_metrics) {
    stop("`metric` must be one of: ",
         paste(shQuote(allowed_metrics), collapse = ", "), ".",
         call. = FALSE)
  }
  if (!is.numeric(obs_value) || length(obs_value) != 1L ||
      !is.finite(obs_value)) {
    stop("`obs_value` must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(pi_hat) || length(pi_hat) != 1L ||
      !is.finite(pi_hat) || pi_hat <= 0 || pi_hat >= 1) {
    stop("`pi_hat` must be a numeric scalar in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 2 || k != as.integer(k)) {
    stop("`k` must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) ||
      N < 1 || N != as.integer(N)) {
    stop("`N` must be an integer >= 1.", call. = FALSE)
  }
  k <- as.integer(k)
  N <- as.integer(N)

  # v0.7.1: the stipulated four-band partition of q and the modal-band
  # qualifier are retired (they reproduced the labeled-band convention the
  # framework replaces). The arguments remain in the signature so 0.6.x
  # call sites keep working; non-default values draw a deprecation warning
  # and are ignored.
  default_bands  <- c(0.5, 0.625, 0.75, 0.875, 1.0)
  default_labels <- c("Poor", "Moderate", "Strong", "Excellent")
  if (!isTRUE(all.equal(bands, default_bands)) ||
      !identical(band_labels, default_labels)) {
    warning("`bands` / `band_labels` are deprecated as of grassr 0.7.1 and ",
            "ignored: the stipulated q-band partition and modal-band ",
            "qualifier are retired. The report now carries a consistency ",
            "band on q (test inversion) and a pooled percentile.",
            call. = FALSE)
  }

  # `method = "empirical"` is serviced by the bundled empirical q_hat
  # surface when `surface_data` is NULL. The ICC branch still requires a
  # caller-supplied reference curve because its closed form depends on
  # the full subject-prevalence distribution F (see the ICC block below).

  notes <- character(0L)
  # reference_used tracks which path actually produced the reference
  # curve (transparency-over-silence rule). For non-ICC metrics the
  # closed form is exact and reference_used = "closed-form". For ICC
  # the resolution order is fitted -> oracle -> caller-supplied.
  reference_used <- "closed-form"

  # ---- Closed-form reference curve on the q-grid -------------------------
  # For PABAK, Fleiss kappa, AC1, Krippendorff alpha, E[metric] depends on
  # (q, pi_+) only under the symmetric DGP (see
  # paper2/code/04_reference_closed_form.R). We use pi_hat as the plug-in
  # for pi_+ and build the 501-point q-grid lookup inline.
  q_grid <- seq(0.5, 1.0, length.out = 501L)

  if (metric == "icc") {
    # ICC depends on the full subject-prevalence distribution F (not pi_hat alone). Resolution order:
    #   1. Caller-supplied surface_data$reference_curve (wins if present)
    #   2. If reference_type = "fitted" (default) and N is in the fitted-
    #      reference sim range: use fitted_icc_reference_curves. glmer-fitted
    #      tau2 from `ratings` pins the (mu, tau2) F_key; without `ratings`,
    #      nearest-M1 is used and flagged.
    #   3. If reference_type = "oracle" or N exceeds fitted sim range:
    #      use icc_reference_curves (oracle closed form).
    #   4. Error with guidance.
    if (!is.null(surface_data) && !is.null(surface_data$reference_curve)) {
      ref_curve <- surface_data$reference_curve
      if (!is.numeric(ref_curve) || length(ref_curve) != length(q_grid)) {
        stop("`surface_data$reference_curve` must be numeric with length ",
             length(q_grid), " (one value per q-grid point on [0.5, 1]).",
             call. = FALSE)
      }
      reference_used <- "user-supplied"
      notes <- c(notes,
                 "ICC reference curve supplied by caller (overrides bundle).")
    } else if (reference_type == "fitted") {
      icc_lookup <- lookup_fitted_icc_reference_curve(
        pi_hat = pi_hat, k = k, N = N,
        q_grid = q_grid, ratings = ratings
      )
      if (is.null(icc_lookup)) {
        # Fitted unavailable (k or N out of range, or sysdata missing): fall
        # back to oracle with a note that names which dimension was the gap.
        reference_used <- "oracle-icc-fallback"
        icc_lookup <- lookup_icc_reference_curve(pi_hat = pi_hat,
                                                 q_grid = q_grid,
                                                 ratings = ratings)
        if (is.null(icc_lookup)) {
          stop("`metric = \"icc\"` cannot be resolved: neither ",
               "`fitted_icc_reference_curves` nor `icc_reference_curves` ",
               "sysdata is available and no `surface_data$reference_curve` ",
               "was supplied.", call. = FALSE)
        }
        gap_dim <- character(0L)
        bundle <- tryCatch(
          get("fitted_icc_reference_curves",
              envir = asNamespace("grassr"), inherits = FALSE),
          error = function(e) NULL
        )
        if (!is.null(bundle)) {
          if (as.numeric(k) > max(bundle$k_grid) + 1) {
            gap_dim <- c(gap_dim, sprintf("k=%s (fitted-ICC k_grid maxes at %d)",
                                          as.character(k), max(bundle$k_grid)))
          }
          if (as.numeric(N) > max(bundle$N_grid) + 1) {
            gap_dim <- c(gap_dim, sprintf("N=%s (fitted-ICC N_grid maxes at %d)",
                                          as.character(N), max(bundle$N_grid)))
          }
        }
        gap_msg <- if (length(gap_dim))
          paste0("Fitted-ICC reference unavailable at ",
                 paste(gap_dim, collapse = ", "),
                 "; using oracle ICC reference (GLMM-gap not corrected). ",
                 "Treat the surface position as an approximation.")
        else
          "Fitted ICC reference unavailable at this (k, N); falling back to oracle reference."
        notes <- c(notes, gap_msg)
      }
      ref_curve <- icc_lookup$reference_curve
      notes <- c(notes, icc_lookup$notes)
      # If fitted lookup succeeded (icc_lookup wasn't NULL on first attempt),
      # reference_used stays "closed-form" placeholder; promote to "fitted-icc"
      # whenever the fitted path produced the curve.
      if (identical(reference_used, "closed-form")) {
        reference_used <- "fitted-icc"
      }
    } else {
      # reference_type == "oracle"
      reference_used <- "oracle-icc-explicit"
      icc_lookup <- lookup_icc_reference_curve(pi_hat = pi_hat,
                                               q_grid = q_grid,
                                               ratings = ratings)
      if (is.null(icc_lookup)) {
        stop("`metric = \"icc\"` cannot be resolved with reference_type = ",
             "\"oracle\": bundled `icc_reference_curves` sysdata is unavailable.",
             call. = FALSE)
      }
      ref_curve <- icc_lookup$reference_curve
      notes <- c(notes, icc_lookup$notes)
    }
  } else {
    ref_curve <- closed_form_reference_curve(metric = metric,
                                             pi_plus = pi_hat,
                                             q_grid = q_grid)
  }

  # ---- Invert obs_value to q_hat + delta-method SE -----------------------
  inv <- invert_metric_to_q(obs_value = obs_value,
                            ref_curve = ref_curve,
                            q_grid = q_grid)
  q_hat <- inv$q_hat
  dEdq  <- inv$dEdq
  if (!is.null(inv$note)) notes <- c(notes, inv$note)

  # SE of the observed mean approximated from metric variance at q_hat
  # together with the sample-size context (k, N). Without a calibrated
  # Monte-Carlo SD we approximate by the Bernoulli-agreement variance:
  #   Var(observed) ~ [p_a(q_hat) * (1 - p_a(q_hat))] / effective_n
  # where p_a = 1 - 2q(1-q) under the symmetric DGP. Effective n uses
  # (k choose 2) * N, the number of within-subject rater pairs summed
  # over subjects. This is a rough Column-B fallback; empirical / sim-
  # derived SEs via `surface_data$sd_metric` override when supplied.
  sd_metric <- surface_data$sd_metric %||% approx_metric_sd(metric = metric,
                                                            q_hat  = q_hat,
                                                            pi_hat = pi_hat,
                                                            k = k, N = N)
  if (!is.finite(dEdq) || abs(dEdq) < 1e-10) {
    se_q_hat <- NA_real_
    notes <- c(notes,
               "Delta-method SE undefined: dE/dq near zero at q_hat.")
  } else {
    se_q_hat <- sd_metric / abs(dEdq)
  }

  # ---- Sweep positioning (v0.7.1 convention) ------------------------------
  # The observed coefficient is evaluated against EVERY calibrated quality
  # level q at the matched design (F/pi_hat, k, N):
  #   p(q) = P(coefficient <= obs_value | panel of quality q, this design)
  # Three read-outs of one object (design/v0.7.1_position_redesign.md):
  #   sweep      -- the p(q) profile across the calibrated q axis
  #   band       -- 95% test-inversion consistency band on q: the quality
  #                 levels with 0.025 <= p(q) <= 0.975
  #   percentile -- pooled percentile: trapezoid-weighted average of p(q),
  #                 i.e. the observed coefficient's position within the
  #                 design's full achievable range. Monotone in obs_value
  #                 by construction.
  # The nearest-q_hat-cell percentile is retired: selecting the reference
  # cohort by a statistic derived from the coefficient itself made the
  # percentile a non-monotone sawtooth in the coefficient (panel review
  # 2026-07-05). No q cell is selected here; the whole axis is consulted.
  sampling_method_used <- "delta"
  percentile <- NA_real_
  percentile_basis <- NA_character_
  sweep <- NULL
  band <- NULL

  if (!is.null(surface_data$q_grid_per_rep) &&
      !is.null(surface_data$q_grid)) {
    notes <- c(notes,
               "`surface_data$q_grid_per_rep` supplied but not consumed internally; ",
               "using bundled sweep / delta-method fallback.")
  }

  if (!is.null(surface_data$per_rep)) {
    # Legacy caller-supplied single-cohort hook: an empirical sampling
    # distribution at the caller's own design. Honored for reproducibility
    # audits; yields a plain cohort percentile with no sweep or band.
    pr <- as.numeric(surface_data$per_rep)
    pr <- pr[is.finite(pr)]
    if (length(pr) < 2L) {
      stop("`surface_data$per_rep` must contain >= 2 finite values.",
           call. = FALSE)
    }
    percentile <- mean(pr <= obs_value)
    percentile_basis <- "user-supplied-cohort"
    sampling_method_used <- "empirical"
    notes <- c(notes,
               "Percentile from caller-supplied per_rep cohort; sweep/band not derived.")
  } else {
    sweep_lookup <- NULL
    if (method == "empirical" && metric != "icc") {
      sweep_lookup <- lookup_empirical_q_sweep(metric = metric,
                                               pi_hat = pi_hat,
                                               k = k, N = N)
      if (is.null(sweep_lookup)) {
        notes <- c(notes,
                   "Bundled empirical q_hat surface unavailable; falling back to delta-method sweep.")
      }
    }

    if (!is.null(sweep_lookup)) {
      # Empirical sweep in q_hat space: rank q_hat within each calibrated
      # cell's q_hat_rep distribution (monotone-equivalent to ranking
      # obs_value within the cell's coefficient distribution).
      p_vals <- vapply(seq_along(sweep_lookup$q_true), function(j) {
        empirical_cdf_at(q_hat,
                         quantiles = sweep_lookup$quantiles[j, ],
                         probs     = sweep_lookup$probs)
      }, numeric(1L))
      sweep <- data.frame(q = sweep_lookup$q_true, p = p_vals)
      sampling_method_used <- "empirical"
      percentile_basis <- "pooled-achievable-range"
      if (length(sweep_lookup$clamp_notes)) {
        notes <- c(notes, sweep_lookup$clamp_notes)
      }
    } else {
      # Delta-method sweep in coefficient space: at each q, approximate the
      # sampling distribution as Normal(E[metric](q), sd_metric(q)). Serves
      # the summary-stats-only path and the ICC branch (whose reference
      # curve carries the F-shape conditioning).
      q_sweep <- c(seq(0.55, 0.90, by = 0.05), 0.92, 0.94, 0.95, 0.97, 0.99)
      p_vals <- vapply(q_sweep, function(qq) {
        mu_q <- approx_at(ref_curve, q_grid, qq)
        sd_q <- approx_metric_sd(metric = metric, q_hat = qq,
                                 pi_hat = pi_hat, k = k, N = N)
        if (!is.finite(mu_q) || !is.finite(sd_q) || sd_q <= 0) {
          return(NA_real_)
        }
        stats::pnorm(obs_value, mean = mu_q, sd = sd_q)
      }, numeric(1L))
      keep <- is.finite(p_vals)
      if (sum(keep) >= 3L) {
        sweep <- data.frame(q = q_sweep[keep], p = p_vals[keep])
        sampling_method_used <- "delta"
        percentile_basis <- "pooled-achievable-range-delta-approx"
      }
    }

    if (!is.null(sweep)) {
      percentile <- pooled_percentile_from_sweep(sweep$q, sweep$p)
      band <- consistency_band_from_sweep(sweep$q, sweep$p, level = 0.95)
      if (!is.null(band$note)) notes <- c(notes, band$note)
    } else {
      notes <- c(notes,
                 "Sweep unavailable (no bundled surface and delta approximation undefined); percentile NA.")
    }
  }

  out <- list(
    observed_value     = as.numeric(obs_value),
    metric             = metric,
    design             = list(pi_hat = as.numeric(pi_hat),
                              k = k, N = N),
    q_hat              = as.numeric(q_hat),
    se_q_hat           = as.numeric(se_q_hat),
    percentile         = as.numeric(percentile),
    percentile_basis   = percentile_basis,
    band               = band,
    sweep              = sweep,
    sampling_method    = sampling_method_used,
    reference_used     = reference_used,
    notes              = notes
  )
  class(out) <- c("grass_surface_position", "list")
  out
}

# ---- Internal: closed-form reference curve on q-grid ----------------------
# E[metric](q | pi_+) for the four agreement-family metrics, built by
# iterating the 501-point q-grid. pi_+ is the plug-in marginal positive rate
# from the observed rating matrix.
closed_form_reference_curve <- function(metric, pi_plus, q_grid) {
  # pi_+ implies an M1 via pi_+ = (1-q) + (2q-1)*M1 -> M1 = (pi_+ - (1-q)) /
  # (2q - 1). But pi_+ is nearly invariant across our q-grid (the design
  # holds F fixed), so we treat pi_+ as an exogenous sufficient statistic
  # here and evaluate E[metric] as if pi_+ were independent of q. This is
  # the same simplification used in the paper's large-N closed-form surface
  # construction: the surface is parameterised by (pi_+, k, N), not by mu /
  # tau2 directly.
  switch(metric,
    pabak = {
      # E[PABAK] = (2q - 1)^2, independent of pi_+.
      (2 * q_grid - 1)^2
    },
    fleiss_kappa = {
      P_bar <- 1 - 2 * q_grid * (1 - q_grid)
      P_e   <- pi_plus^2 + (1 - pi_plus)^2
      (P_bar - P_e) / (1 - P_e)
    },
    mean_ac1 = {
      p_a <- 1 - 2 * q_grid * (1 - q_grid)
      p_e <- 2 * pi_plus * (1 - pi_plus)
      (p_a - p_e) / (1 - p_e)
    },
    krippendorff_a = {
      1 - q_grid * (1 - q_grid) / (pi_plus * (1 - pi_plus))
    },
    stop("Unhandled metric in closed_form_reference_curve: ", metric,
         call. = FALSE)
  )
}

# ---- Internal: invert obs_value -> q_hat on the q-grid --------------------
# Mirrors paper2/code/12_q_inversion.R::invert_one but defensively handles
# non-monotone reference curves (e.g. alpha outside achievable range at
# extreme pi_plus).
invert_metric_to_q <- function(obs_value, ref_curve, q_grid) {
  note <- NULL
  # Strip non-finite rows (alpha can produce -Inf at pi_plus -> 0 or 1).
  keep <- is.finite(ref_curve)
  if (sum(keep) < 2L) {
    return(list(q_hat = NA_real_, dEdq = NA_real_,
                note = "Reference curve has insufficient finite values."))
  }
  rc <- ref_curve[keep]
  qg <- q_grid[keep]
  rng <- range(rc)
  if (obs_value <= rng[1]) {
    q_hat <- qg[which.min(rc)]
    note <- sprintf("obs_value %.4f below achievable minimum (%.4f); q_hat clamped.",
                    obs_value, rng[1])
    dEdq <- NA_real_
  } else if (obs_value >= rng[2]) {
    q_hat <- qg[which.max(rc)]
    note <- sprintf("obs_value %.4f above achievable maximum (%.4f); q_hat clamped.",
                    obs_value, rng[2])
    dEdq <- NA_real_
  } else {
    # Monotone inversion via linear interp. `approx` returns NA on ties;
    # we ensure strict monotonicity by sorting on rc when needed.
    ord <- order(rc)
    q_hat <- stats::approx(x = rc[ord], y = qg[ord],
                           xout = obs_value, rule = 2, ties = "ordered")$y
    # Numerical dE/dq at q_hat via central difference on the grid.
    idx <- findInterval(q_hat, qg, all.inside = TRUE)
    idx <- max(2L, min(length(qg) - 1L, idx))
    dq <- qg[idx + 1L] - qg[idx - 1L]
    dE <- rc[idx + 1L] - rc[idx - 1L]
    dEdq <- if (dq > 0) dE / dq else NA_real_
  }
  list(q_hat = q_hat, dEdq = dEdq, note = note)
}

# ---- Internal: evaluate a grid-based reference at a scalar q --------------
approx_at <- function(ref_curve, q_grid, q) {
  keep <- is.finite(ref_curve)
  stats::approx(q_grid[keep], ref_curve[keep], xout = q, rule = 2)$y
}

# ---- Internal: approximate metric SD at (q_hat, pi_hat, k, N) -------------
# Rough expression for the SD of the observed mean-pairwise metric under
# the symmetric DGP. p_a = 1 - 2q(1-q) is the pairwise agreement
# probability. Pairwise comparisons within the same subject are correlated
# (they share C_i), so the effective independent unit is the subject, not
# the rater-pair-within-subject. We use n_eff = N * k / 2 as a compromise
# between the k-pair overdispersion and the N-subject bound; this
# empirically matches simulated SE(PABAK) at (k=5, N=50, balanced) to
# within 10 per cent. Users with sim-derived sd_metric should pass it via
# surface_data$sd_metric for an exact value.
approx_metric_sd <- function(metric, q_hat, pi_hat, k, N) {
  if (!is.finite(q_hat)) return(NA_real_)
  p_a <- 1 - 2 * q_hat * (1 - q_hat)
  # Effective N: calibrated empirically against `multirater_sim_v3`'s
  # per-rep PABAK SD across (k=5, N in {50, 200, 1000}) cells. The naive
  # `n_pairs = choose(k, 2) * N` overstates independence because pairwise
  # comparisons within the same subject share C_i; a shrinkage factor of
  # roughly 0.33 reconciles the formula with the simulated SDs to within
  # ~10 per cent over k in {3, 5, 8, 15}. See
  # paper2/simulation_output/multirater_sim_v3/q_recovery.rds for the
  # calibration reference; callers with sim-derived SDs should pass
  # `surface_data$sd_metric` to bypass this approximation.
  n_eff <- max(choose(k, 2) * N * 0.33, 1)
  var_pa <- p_a * (1 - p_a) / n_eff
  if (!is.finite(var_pa) || var_pa < 0) return(NA_real_)
  # Local slope dMetric/dp_a under symmetric DGP.
  slope <- switch(metric,
    pabak = 2,
    fleiss_kappa = {
      P_e <- pi_hat^2 + (1 - pi_hat)^2
      1 / (1 - P_e)
    },
    mean_ac1 = {
      p_e <- 2 * pi_hat * (1 - pi_hat)
      1 / (1 - p_e)
    },
    krippendorff_a = {
      # alpha = 1 - q(1-q) / [pi_+(1-pi_+)] and p_a = 1 - 2q(1-q) so
      # alpha = 1 - (1 - p_a)/2 / [pi_+(1-pi_+)]. d alpha / d p_a = 1 /
      # (2 * pi_+(1-pi_+)).
      1 / (2 * pi_hat * (1 - pi_hat))
    },
    icc = 1   # placeholder: caller should override via surface_data$sd_metric
  )
  abs(slope) * sqrt(var_pa)
}

# ---- Internal: empirical q_hat lookup from bundled sysdata ----------------
# Resolves the nearest (F_key, k, N, q_true) scenario cell in the bundled
# `empirical_q_hat_surface` package dataset, for the requested metric. Returns
# NULL if the bundled data is unavailable (e.g. sysdata not loaded) or the
# metric is unsupported. The returned list carries the 13-point empirical
# quantile summary at that cell plus any clamping notes describing how far
# off-grid the query was.
lookup_empirical_q_sweep <- function(metric, pi_hat, k, N) {
  # Fetch package-internal data. `::` does not work for non-exported
  # datasets, so use get(); the sysdata is loaded into the package
  # namespace by R's normal data-loading mechanism.
  #
  # v0.7.1: returns ALL calibrated q cells at the matched (F_key, k, N)
  # design -- the sweep convention consults the whole quality axis, so
  # there is deliberately no q selection (and no q clamp note) here.
  surf <- tryCatch(
    get("empirical_q_hat_surface",
        envir = asNamespace("grassr"),
        inherits = FALSE),
    error = function(e) NULL
  )
  if (is.null(surf)) return(NULL)
  if (!metric %in% surf$metrics) return(NULL)

  idx <- surf$index

  k_grid <- sort(unique(idx$k))
  k_near <- k_grid[which.min(abs(k_grid - as.numeric(k)))]

  n_grid <- sort(unique(idx$N))
  n_near <- n_grid[which.min(abs(n_grid - as.numeric(N)))]

  # Nearest F_key via M1 (mean subject prevalence) ~ pi_hat. Restrict to
  # the (k_near, n_near) subgrid so that q_true grids are equal.
  sub <- idx[idx$k == k_near & idx$N == n_near, ]
  if (nrow(sub) == 0L) return(NULL)
  F_keys <- unique(sub[, c("F_key", "M1")])
  fi <- which.min(abs(F_keys$M1 - as.numeric(pi_hat)))
  m1_near <- F_keys$M1[fi]
  F_key_near <- F_keys$F_key[fi]

  sub2 <- sub[sub$F_key == F_key_near, ]
  sub2 <- sub2[order(sub2$q_true), ]
  arr_ids <- as.integer(dimnames(surf$quantiles)$scenario_id)
  rows_in_arr <- match(as.integer(sub2$scenario_id), arr_ids)
  keep <- !is.na(rows_in_arr)
  sub2 <- sub2[keep, , drop = FALSE]
  rows_in_arr <- rows_in_arr[keep]
  if (nrow(sub2) < 2L) return(NULL)

  qmat <- surf$quantiles[rows_in_arr, metric, , drop = FALSE]
  qmat <- matrix(qmat, nrow = nrow(sub2))
  finite_rows <- apply(qmat, 1L, function(v) sum(is.finite(v)) >= 2L)
  if (sum(finite_rows) < 2L) return(NULL)

  clamp_notes <- character(0L)
  if (as.numeric(k) != k_near) {
    clamp_notes <- c(clamp_notes,
                     sprintf("k=%s clamped to nearest sim-grid k=%d.",
                             as.character(k), k_near))
  }
  if (as.numeric(N) != n_near) {
    clamp_notes <- c(clamp_notes,
                     sprintf("N=%s clamped to nearest sim-grid N=%d.",
                             as.character(N), n_near))
  }
  if (abs(as.numeric(pi_hat) - m1_near) > 0.05) {
    clamp_notes <- c(clamp_notes,
                     sprintf("pi_hat=%.3f clamped to nearest sim F_key with M1=%.3f.",
                             as.numeric(pi_hat), m1_near))
  }

  list(
    q_true      = as.numeric(sub2$q_true[finite_rows]),
    quantiles   = qmat[finite_rows, , drop = FALSE],
    probs       = as.numeric(surf$probs),
    F_key       = F_key_near,
    k_nearest   = k_near,
    N_nearest   = n_near,
    M1_nearest  = m1_near,
    clamp_notes = clamp_notes
  )
}

# ---- Internal: ICC reference curve lookup from bundled sysdata ------------
# Resolves the nearest sim F_key and returns the 501-point E[ICC](q) curve
# at that F_key. The ICC closed form depends on the full subject-prevalence
# distribution F (mu, tau2 for logit-normal), not pi_hat alone.
#
# Selection path:
#   - If `ratings` is supplied AND lme4 is available: fit
#     glmer(rating ~ 1 + (1|subject), family=binomial) to estimate (mu, tau2)
#     from the practitioner's data, then pick the F_key at nearest (mu, tau2)
#     in (mu, log(tau2)) distance.
#   - Else: fall back to nearest-M1 F_key, with a prominent caveat note
#     warning that tau2 is unidentified and the chosen reference curve may
#     not span obs_value.
#
# Returns NULL if the bundled sysdata is unavailable. Returns a list with
# `reference_curve` and `notes` otherwise.
lookup_icc_reference_curve <- function(pi_hat, q_grid, ratings = NULL) {
  bundle <- tryCatch(
    get("icc_reference_curves",
        envir = asNamespace("grassr"),
        inherits = FALSE),
    error = function(e) NULL
  )
  if (is.null(bundle)) return(NULL)

  # Sanity: the bundled q-grid must match the caller's 501-point grid.
  if (length(bundle$q_grid) != length(q_grid) ||
      any(abs(bundle$q_grid - q_grid) > .Machine$double.eps^0.5)) {
    # Interpolate onto caller's q-grid. Should never happen under the
    # current design but keep the fallback.
    needs_interp <- TRUE
  } else {
    needs_interp <- FALSE
  }

  idx <- bundle$index
  # Parse (mu, tau2) from F_key strings for logit-normal entries. Format is
  # "LN_mu=%+.3f_tau2=%.4f" (see paper2/code/06_grid.R); discrete-mixture
  # entries get NA.
  parsed_mu   <- suppressWarnings(as.numeric(
    sub("^LN_mu=([+-]?[0-9.]+)_tau2=.*$", "\\1", idx$F_key)))
  parsed_tau2 <- suppressWarnings(as.numeric(
    sub("^LN_mu=[+-]?[0-9.]+_tau2=([0-9.]+)$", "\\1", idx$F_key)))
  idx$mu   <- parsed_mu
  idx$tau2 <- parsed_tau2

  fit_notes <- character(0L)
  fit <- NULL
  if (!is.null(ratings)) {
    fit <- fit_tau2_from_ratings(ratings)
    if (!is.null(fit$note)) fit_notes <- c(fit_notes, fit$note)
  }

  if (!is.null(fit) && is.finite(fit$mu) && is.finite(fit$tau2) &&
      fit$tau2 > 0) {
    # glmer path: pick F_key at nearest (mu, log(tau2)).
    valid <- is.finite(idx$mu) & is.finite(idx$tau2) & idx$tau2 > 0
    dist <- rep(Inf, nrow(idx))
    dist[valid] <- (idx$mu[valid] - fit$mu)^2 +
                   (log(idx$tau2[valid]) - log(fit$tau2))^2
    cand_idx <- which.min(dist)
    fit_notes <- c(fit_notes,
                   sprintf("ICC F_key picked via glmer: mu_hat=%.3f, tau2_hat=%.3f -> nearest F_key tau2=%.4f, mu=%.3f.",
                           fit$mu, fit$tau2,
                           idx$tau2[cand_idx], idx$mu[cand_idx]))
  } else {
    # Fallback: nearest-M1 only. tau2 is unidentified from pi_hat alone.
    cand_idx <- which.min(abs(idx$M1 - as.numeric(pi_hat)))
    fit_notes <- c(fit_notes,
                   "ICC F_key picked via nearest-M1 only; tau2 not estimated. ",
                   "Pass `ratings = <rating matrix>` for a glmer-fitted F_key, ",
                   "or `surface_data$reference_curve` to override directly.")
  }
  F_key_near <- idx$F_key[cand_idx]
  M1_near    <- idx$M1[cand_idx]
  F_family   <- idx$F_family[cand_idx]

  ref_row <- which(rownames(bundle$curves) == F_key_near)
  if (length(ref_row) == 0L) return(NULL)
  ref_curve <- as.numeric(bundle$curves[ref_row, ])

  if (needs_interp) {
    keep <- is.finite(ref_curve)
    if (sum(keep) < 2L) return(NULL)
    ref_curve <- stats::approx(bundle$q_grid[keep], ref_curve[keep],
                               xout = q_grid, rule = 2)$y
  }

  notes <- character(0L)
  notes <- c(notes, fit_notes)
  notes <- c(notes,
             sprintf("ICC reference curve from bundled F_key=%s (family=%s, M1=%.3f).",
                     F_key_near, F_family, M1_near))
  if (abs(as.numeric(pi_hat) - M1_near) > 0.05) {
    notes <- c(notes,
               sprintf("pi_hat=%.3f more than 0.05 from nearest sim F_key M1=%.3f; treat ICC position as coarse.",
                       as.numeric(pi_hat), M1_near))
  }

  list(
    reference_curve = ref_curve,
    F_key_nearest   = F_key_near,
    F_family        = F_family,
    M1_nearest      = M1_near,
    notes           = notes
  )
}

# ---- Internal: fitted-ICC reference lookup (GLMM-gap corrected) ----------
# Resolves the nearest sim (F_key, k, N) cell and returns the 501-point
# fitted E[ICC](q) reference curve built as oracle_ref + bias_emp correction
# from paper2/simulation_output/multirater_sim_v3/q_recovery_fitted_icc.rds.
#
# This is the correct reference for practitioners whose obs_ICC was computed
# via glmer (the standard workflow). The oracle reference routinely over-
# shoots the practitioner's scale due to the GLMM gap (framework_notes.md
# Sec.0.4.iii); the fitted reference pulls the curve up to match glmer output.
#
# Returns NULL when N is outside the fitted sim range ({50, 200}) or the
# bundled sysdata is missing; caller falls back to oracle with a note.
lookup_fitted_icc_reference_curve <- function(pi_hat, k, N, q_grid,
                                              ratings = NULL) {
  bundle <- tryCatch(
    get("fitted_icc_reference_curves",
        envir = asNamespace("grassr"),
        inherits = FALSE),
    error = function(e) NULL
  )
  if (is.null(bundle)) return(NULL)

  # If N or k exceeds the fitted sim range, fitted correction isn't available;
  # caller falls back to oracle. Small-N / in-range N and k get a nearest match.
  # The fitted-ICC bundle currently covers k in {3, 5, 8, 15}, N in {50, 200};
  # extrapolating the GLMM-gap correction much beyond k=15 would over-correct
  # because the gap shrinks with k.
  N_grid <- bundle$N_grid
  if (as.numeric(N) > max(N_grid) + 1) {
    return(NULL)  # signals caller to use oracle instead
  }
  N_near <- N_grid[which.min(abs(N_grid - as.numeric(N)))]

  k_grid <- bundle$k_grid
  if (as.numeric(k) > max(k_grid) + 1) {
    return(NULL)  # signals caller to use oracle instead
  }
  k_near <- k_grid[which.min(abs(k_grid - as.numeric(k)))]

  idx <- bundle$index
  idx$mu   <- suppressWarnings(as.numeric(
    sub("^LN_mu=([+-]?[0-9.]+)_tau2=.*$", "\\1", idx$F_key)))
  idx$tau2 <- suppressWarnings(as.numeric(
    sub("^LN_mu=[+-]?[0-9.]+_tau2=([0-9.]+)$", "\\1", idx$F_key)))

  fit_notes <- character(0L)
  fit <- NULL
  if (!is.null(ratings)) {
    fit <- fit_tau2_from_ratings(ratings)
    if (!is.null(fit$note)) fit_notes <- c(fit_notes, fit$note)
  }

  if (!is.null(fit) && is.finite(fit$mu) && is.finite(fit$tau2) &&
      fit$tau2 > 0) {
    valid <- is.finite(idx$mu) & is.finite(idx$tau2) & idx$tau2 > 0
    dist <- rep(Inf, nrow(idx))
    dist[valid] <- (idx$mu[valid] - fit$mu)^2 +
                   (log(idx$tau2[valid]) - log(fit$tau2))^2
    cand_idx <- which.min(dist)
    fit_notes <- c(fit_notes,
                   sprintf("Fitted-ICC F_key picked via glmer: mu_hat=%.3f, tau2_hat=%.3f -> F_key tau2=%.4f, mu=%.3f.",
                           fit$mu, fit$tau2,
                           idx$tau2[cand_idx], idx$mu[cand_idx]))
  } else {
    cand_idx <- which.min(abs(idx$M1 - as.numeric(pi_hat)))
    fit_notes <- c(fit_notes,
                   "Fitted-ICC F_key picked via nearest-M1 only; tau2 not estimated. ",
                   "Pass `ratings = <rating matrix>` for a glmer-fitted F_key.")
  }

  F_key_near <- idx$F_key[cand_idx]
  M1_near    <- idx$M1[cand_idx]
  F_family   <- idx$F_family[cand_idx]

  ref_curve <- bundle$curves[
    F_key_near,
    as.character(k_near),
    as.character(N_near),
  ]
  ref_curve <- as.numeric(ref_curve)

  # Interpolate onto caller q-grid if different from bundled.
  if (length(bundle$q_grid) != length(q_grid) ||
      any(abs(bundle$q_grid - q_grid) > .Machine$double.eps^0.5)) {
    keep <- is.finite(ref_curve)
    if (sum(keep) < 2L) return(NULL)
    ref_curve <- stats::approx(bundle$q_grid[keep], ref_curve[keep],
                               xout = q_grid, rule = 2)$y
  }

  if (all(!is.finite(ref_curve))) return(NULL)

  notes <- character(0L)
  notes <- c(notes, fit_notes)
  notes <- c(notes,
             sprintf("Fitted-ICC reference (GLMM-gap corrected) at F_key=%s, k=%d, N=%d (family=%s, M1=%.3f).",
                     F_key_near, k_near, N_near, F_family, M1_near))
  if (as.numeric(k) != k_near) {
    notes <- c(notes,
               sprintf("k=%s clamped to nearest sim-grid k=%d.",
                       as.character(k), k_near))
  }
  if (as.numeric(N) != N_near) {
    notes <- c(notes,
               sprintf("N=%s clamped to nearest sim-grid N=%d.",
                       as.character(N), N_near))
  }
  if (abs(as.numeric(pi_hat) - M1_near) > 0.05 && is.null(ratings)) {
    notes <- c(notes,
               sprintf("pi_hat=%.3f more than 0.05 from nearest sim F_key M1=%.3f; treat ICC position as coarse.",
                       as.numeric(pi_hat), M1_near))
  }

  list(
    reference_curve = ref_curve,
    F_key_nearest   = F_key_near,
    F_family        = F_family,
    M1_nearest      = M1_near,
    k_nearest       = k_near,
    N_nearest       = N_near,
    notes           = notes
  )
}

# ---- Internal: estimate (mu, tau2) from a rating matrix via glmer ---------
# Fits the one-way random-subject-intercept logistic mixed model
#   logit(P(rating = 1 | subject_i)) = mu + u_i,   u_i ~ N(0, tau2)
# on the practitioner's rating data, and returns the MLE (mu_hat, tau2_hat).
# Accepts a k x N integer matrix (rows = raters, cols = subjects) with 0/1
# entries, or a long data.frame with columns `subject` and `rating`. Graceful
# fallback on missing lme4, tiny N, or non-convergence: returns a list with
# NA estimates and a note explaining why.
fit_tau2_from_ratings <- function(ratings) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(list(
      mu   = NA_real_, tau2 = NA_real_,
      note = "Package `lme4` is not installed; tau2 cannot be estimated. Install lme4, or supply `surface_data$reference_curve` directly."
    ))
  }

  if (is.matrix(ratings)) {
    k <- nrow(ratings)
    n_subj <- ncol(ratings)
    if (k < 2L || n_subj < 10L) {
      return(list(
        mu = NA_real_, tau2 = NA_real_,
        note = sprintf("Rating matrix is %d x %d; need k>=2 and N>=10 for glmer. Falling back to nearest-M1 F_key.",
                       k, n_subj)
      ))
    }
    long <- data.frame(
      subject = rep(seq_len(n_subj), each = k),
      rating  = as.integer(as.vector(ratings))
    )
  } else if (is.data.frame(ratings) &&
             all(c("subject", "rating") %in% names(ratings))) {
    long <- data.frame(
      subject = as.integer(as.factor(ratings$subject)),
      rating  = as.integer(ratings$rating)
    )
  } else {
    return(list(
      mu = NA_real_, tau2 = NA_real_,
      note = "`ratings` must be a k x N integer matrix (rows=raters, cols=subjects) or a data.frame with `subject` and `rating` columns."
    ))
  }

  if (!all(long$rating %in% c(0L, 1L))) {
    return(list(
      mu = NA_real_, tau2 = NA_real_,
      note = "`ratings` must contain only 0/1 values after coercion."
    ))
  }

  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      lme4::glmer(rating ~ 1 + (1 | subject), data = long,
                  family = stats::binomial(link = "logit"),
                  control = lme4::glmerControl(optimizer = "bobyqa"))
    )),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      mu = NA_real_, tau2 = NA_real_,
      note = sprintf("glmer fit failed (%s); falling back to nearest-M1 F_key.",
                     conditionMessage(fit))
    ))
  }

  mu_hat   <- as.numeric(lme4::fixef(fit)[1])
  vc       <- lme4::VarCorr(fit)
  tau2_hat <- tryCatch(as.numeric(vc$subject[1, 1]), error = function(e) NA_real_)

  if (!is.finite(mu_hat) || !is.finite(tau2_hat) || tau2_hat <= 0) {
    return(list(
      mu = NA_real_, tau2 = NA_real_,
      note = "glmer returned non-finite or degenerate estimates (mu_hat or tau2_hat); falling back to nearest-M1 F_key."
    ))
  }

  list(mu = mu_hat, tau2 = tau2_hat, note = NULL)
}

# ---- Internal: empirical CDF at point x from a quantile summary -----------
# Given quantile values `q_at` at probabilities `probs`, compute P(X <= x)
# via piecewise-linear interpolation on the empirical CDF. Extrapolates
# flat outside the observed quantile range.
empirical_cdf_at <- function(x, quantiles, probs) {
  if (!is.finite(x)) return(NA_real_)
  # Anchor the CDF at the tails to keep the linear interp bounded in [0, 1].
  # We extend by extrapolating the first/last slope a small amount. For our
  # purposes clamping to [min(probs), max(probs)] (=~ 0.01 / 0.99) beyond
  # the observed range is acceptable because position_on_surface consumers
  # always supply q_hat values that live inside the tight (SD ~0.01-0.05)
  # sampling distributions; extreme extrapolation is a clamp/flag case
  # already handled upstream.
  ok <- is.finite(quantiles)
  if (sum(ok) < 2L) return(NA_real_)
  qx <- quantiles[ok]
  pr <- probs[ok]
  ord <- order(qx)
  qx <- qx[ord]
  pr <- pr[ord]
  # Strict monotonisation: if duplicates, jitter by a tiny epsilon.
  dups <- duplicated(qx)
  if (any(dups)) {
    eps <- .Machine$double.eps^0.5
    qx[dups] <- qx[dups] + cumsum(dups)[dups] * eps
  }
  # Below the empirical-quantile envelope: the stored tail is the 1st
  # percentile (probs[1] ~ 0.01), not an absolute lower bound. Treat
  # queries strictly below the lowest stored quantile as CDF ~ 0; above
  # the highest as CDF ~ 1. This avoids leaking ~probs[1] mass into the
  # leftmost band when the empirical distribution lies well inside the
  # band partition.
  if (x <= qx[1]) return(0)
  if (x >= qx[length(qx)]) return(1)
  stats::approx(x = qx, y = pr, xout = x, rule = 2, ties = "ordered")$y
}

# ---- Internal: pooled percentile from the q sweep --------------------------
# Trapezoid-weighted average of p(q) over the calibrated q axis: the CDF of
# the pooled mixture (equal density per unit q, so the non-uniform grid
# spacing does not distort the pool) evaluated at the observed value.
# Monotone in the observed coefficient by construction.
pooled_percentile_from_sweep <- function(q, p) {
  ok <- is.finite(p)
  q <- q[ok]; p <- p[ok]
  n <- length(q)
  if (n == 0L) return(NA_real_)
  if (n == 1L) return(p)
  w <- numeric(n)
  w[1L] <- (q[2L] - q[1L]) / 2
  w[n]  <- (q[n] - q[n - 1L]) / 2
  if (n > 2L) w[2:(n - 1L)] <- (q[3:n] - q[1:(n - 2L)]) / 2
  sum(w * p) / sum(w)
}

# ---- Internal: consistency band on q via test inversion --------------------
# Quality level q is CONSISTENT with the observed coefficient when
# alpha/2 <= p(q) <= 1 - alpha/2 (two-sided test inversion at `level`).
# p(q) decreases in q (a fixed observation ranks lower within higher-quality
# cohorts); band endpoints are interpolated linearly in q at the crossings.
# Monte-Carlo non-monotonicity is handled by taking the OUTERMOST crossings
# (conservative wide band). Open ends at the calibrated grid boundary are
# flagged rather than silently truncated.
consistency_band_from_sweep <- function(q, p, level = 0.95) {
  alpha <- 1 - level
  lo_p <- 1 - alpha / 2   # p above this: observation below quality q
  hi_p <- alpha / 2       # p below this: observation above quality q
  ok <- is.finite(p)
  q <- q[ok]; p <- p[ok]
  n <- length(q)
  if (n < 2L) {
    return(list(lo = NA_real_, hi = NA_real_, level = level,
                open_low = NA, open_high = NA,
                note = "Consistency band undefined: sweep too short."))
  }
  cons <- which(p >= hi_p & p <= lo_p)

  interp_cross <- function(j1, j2, target) {
    # Linear interpolation of the q at which p crosses `target` between
    # adjacent sweep points j1 < j2.
    if (abs(p[j2] - p[j1]) < 1e-12) return((q[j1] + q[j2]) / 2)
    q[j1] + (target - p[j1]) * (q[j2] - q[j1]) / (p[j2] - p[j1])
  }

  if (length(cons) == 0L) {
    if (all(p > lo_p)) {
      # Observation below even the lowest calibrated quality.
      return(list(lo = NA_real_, hi = q[1L], level = level,
                  open_low = TRUE, open_high = FALSE,
                  note = sprintf(
                    "Observed value below the sampling range of the lowest calibrated quality (q = %.2f); consistency band open below the calibrated grid.",
                    q[1L])))
    }
    if (all(p < hi_p)) {
      # Observation above even the highest calibrated quality.
      return(list(lo = q[n], hi = NA_real_, level = level,
                  open_low = FALSE, open_high = TRUE,
                  note = sprintf(
                    "Observed value above the sampling range of the highest calibrated quality (q = %.2f); consistency band open above the calibrated grid.",
                    q[n])))
    }
    # p jumps across the whole consistent range between two grid points
    # (very tight sampling distributions relative to grid spacing):
    # interpolate both crossings inside that gap.
    j <- which(p[-n] > lo_p & p[-1L] < hi_p)
    if (length(j)) {
      j <- j[1L]
      lo <- interp_cross(j, j + 1L, lo_p)
      hi <- interp_cross(j, j + 1L, hi_p)
      return(list(lo = min(lo, hi), hi = max(lo, hi), level = level,
                  open_low = FALSE, open_high = FALSE,
                  note = "Consistency band narrower than the calibrated q-grid spacing; endpoints interpolated within one grid gap."))
    }
    return(list(lo = NA_real_, hi = NA_real_, level = level,
                open_low = NA, open_high = NA,
                note = "Consistency band undefined: sweep profile irregular."))
  }

  j0 <- min(cons)
  j1 <- max(cons)
  if (j0 == 1L) {
    lo <- q[1L]; open_low <- TRUE
  } else {
    lo <- interp_cross(j0 - 1L, j0, lo_p); open_low <- FALSE
  }
  if (j1 == n) {
    hi <- q[n]; open_high <- TRUE
  } else {
    hi <- interp_cross(j1, j1 + 1L, hi_p); open_high <- FALSE
  }
  note <- NULL
  if (open_low && open_high) {
    note <- "Consistency band spans the entire calibrated quality grid: this design does not constrain panel quality."
  } else if (open_low) {
    note <- "Consistency band open at the low edge of the calibrated grid."
  } else if (open_high) {
    note <- "Consistency band open at the high edge of the calibrated grid."
  }
  list(lo = lo, hi = hi, level = level,
       open_low = open_low, open_high = open_high, note = note)
}

# ---- %||% coalescing helper (local, does not export) ---------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

#' @export
print.grass_surface_position <- function(x, digits = 3, ...) {
  cat("grass surface-position report (sweep convention, v0.7.1)\n",
      "  metric               : ", x$metric, "\n",
      "  observed value       : ", formatC(x$observed_value, digits = digits,
                                           format = "f"), "\n",
      "  design (pi_hat,k,N)  : (",
      formatC(x$design$pi_hat, digits = digits, format = "f"), ", ",
      x$design$k, ", ", x$design$N, ")\n",
      "  implied quality q_hat: ",
      formatC(x$q_hat, digits = digits, format = "f"), " +/- ",
      formatC(x$se_q_hat, digits = digits, format = "f"), "\n",
      "  pooled percentile    : ",
      if (is.finite(x$percentile))
        sprintf("%.1f (of the design's achievable range)",
                100 * x$percentile)
      else "NA", "\n", sep = "")
  if (!is.null(x$band)) {
    cat("  consistency band     : ", format_consistency_band(x$band),
        "\n", sep = "")
  }
  cat("  sampling method      : ", x$sampling_method, "\n", sep = "")
  if (length(x$notes)) {
    cat("  notes                :\n")
    for (n in x$notes) cat("    - ", n, "\n", sep = "")
  }
  invisible(x)
}

# ---- Internal: render a consistency band as one human-readable string -----
format_consistency_band <- function(band, digits = 2) {
  if (is.null(band)) return("not derived")
  lvl <- sprintf("%d%%", round(100 * (band$level %||% 0.95)))
  if (is.na(band$lo) && is.na(band$hi)) return(paste0(lvl, " band undefined"))
  if (is.na(band$lo)) {
    return(sprintf("quality <= %.*f (below calibrated grid; %s)",
                   digits, band$hi, lvl))
  }
  if (is.na(band$hi)) {
    return(sprintf("quality >= %.*f (above calibrated grid; %s)",
                   digits, band$lo, lvl))
  }
  lo_mark <- if (isTRUE(band$open_low))  "<=" else ""
  hi_mark <- if (isTRUE(band$open_high)) "+"  else ""
  sprintf("consistent with panel quality %s%.*f-%.*f%s (%s)",
          lo_mark, digits, band$lo, digits, band$hi, hi_mark, lvl)
}

#' Coerce a grass_surface_position to a one-row data.frame
#'
#' @param x A `grass_surface_position` object.
#' @param row.names,optional,... Standard arguments; ignored except
#'   `row.names`.
#' @return A one-row data.frame summarising the position: pooled
#'   percentile, consistency-band endpoints (`band_lo`, `band_hi`,
#'   `band_open_low`, `band_open_high`), implied quality, and method.
#' @export
as.data.frame.grass_surface_position <- function(x, row.names = NULL,
                                                 optional = FALSE, ...) {
  b <- x$band
  data.frame(
    metric             = x$metric,
    observed_value     = x$observed_value,
    pi_hat             = x$design$pi_hat,
    k                  = x$design$k,
    N                  = x$design$N,
    q_hat              = x$q_hat,
    se_q_hat           = x$se_q_hat,
    percentile         = x$percentile,
    percentile_basis   = x$percentile_basis %||% NA_character_,
    band_lo            = if (!is.null(b)) b$lo else NA_real_,
    band_hi            = if (!is.null(b)) b$hi else NA_real_,
    band_open_low      = if (!is.null(b)) isTRUE(b$open_low) else NA,
    band_open_high     = if (!is.null(b)) isTRUE(b$open_high) else NA,
    band_level         = if (!is.null(b)) (b$level %||% NA_real_) else NA_real_,
    sampling_method    = x$sampling_method,
    stringsAsFactors   = FALSE,
    row.names          = row.names
  )
}
