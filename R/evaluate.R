#' Generate a GRASS Report Card from a rating matrix
#'
#' `grass_report()` is the headline entry point for the v0.2.0 Target-2
#' framework. It takes an `N x k` binary rating matrix and returns a
#' four-field Report Card: the sample summary `(k, N, pi_hat)`, the primary
#' coefficient and its surface position, the cross-coefficient asymmetry
#' diagnostic `delta_hat` and flag, and (when `flag == "divergent"`) the
#' per-rater latent-class fit. The full panel of coefficients, surface
#' percentiles, band probabilities, and reference-surface artifacts ride
#' along on the same object for `summary()`, `as.data.frame()`, and
#' `plot()` access.
#'
#' The body, in order:
#'
#' 1. Normalize `ratings` to a canonical `N x k` integer matrix `Y` and
#'    derive `pi_hat = mean(Y)`, `k = ncol(Y)`, `N = nrow(Y)`. Validate
#'    `k >= 2`; warn at `N < 10`; note at `N < 30`.
#' 2. Compute the panel of observed coefficients
#'    ([compute_panel()] internal): at `k = 2`, PABAK / AC1 / Cohen's
#'    kappa / Krippendorff alpha; at `k >= 3`, PABAK / AC1 / Fleiss kappa /
#'    Krippendorff alpha / ICC.
#' 3. For each panel coefficient, position the observed value on its
#'    DGP-calibrated reference surface via [position_on_surface()].
#' 4. Pick the primary coefficient via Table 2 (`metric = "auto"`) or accept
#'    the user's override.
#' 5. Compute the cross-coefficient percentile spread `delta_hat` (in pp)
#'    via [check_asymmetry()] and tier into `aligned` /  `caution` /
#'    `divergent`.
#' 6. If `flag == "divergent"`: run a [latent_class_fit()] (Dawid-Skene EM
#'    at `k >= 3`; Hui-Walter bounds at `k = 2`) and attach the per-rater
#'    `(Se_j, Sp_j)` table.
#' 7. Assemble the `grass_card` S3 object.
#'
#' @param ratings User input: an `N x k` binary matrix, an `N x k` data.frame
#'   whose columns are 0/1 / logical / 2-level factor, or a list of two
#'   equal-length 0/1 vectors (`k = 2` paired form). See
#'   `?normalize_ratings` for accepted shapes.
#' @param axis One of `"inter"` (default) or `"intra"`. Selects the surface
#'   family. The intra-axis path uses `occasion` to identify viewings.
#' @param metric One of `"auto"` (default; calls [pick_primary_coefficient()]
#'   per Table 2), `"pabak"`, `"ac1"`, `"fleiss_kappa"`, `"krippendorff_a"`,
#'   `"icc"`. Selects which coefficient is the headline in the printed Report
#'   Card; the full panel is always populated.
#' @param occasion Reserved for `axis = "intra"`; ignored when `axis = "inter"`.
#' @param bands Numeric length-5 partition on `q in [0.5, 1]`. Default
#'   `c(0.5, 0.625, 0.75, 0.875, 1.0)`.
#' @param band_labels Character length-4 labels for the bands. Default
#'   `c("Poor", "Moderate", "Strong", "Excellent")`.
#' @param delta_thresholds Length-2 numeric vector `c(caution, divergent)`
#'   in percentile points. Default `c(9.25, 11.75)` (paper §3.2, NP-motivated size-alpha).
#' @param bootstrap_B Integer; bootstrap replicates for the divergent-branch
#'   latent-class CIs. Default `1000L`. Set lower for fast tests.
#' @param verbose Logical; emit progress messages on long calls. Default
#'   `FALSE`.
#' @param ... Reserved for future extension.
#'
#' @return An object of class `c("grass_card", "list")` with fields
#'   `sample`, `coefficient`, `delta`, `panel`, `per_rater`, `surface`,
#'   `call`, `grass_version`, `timestamp`, `inputs`, `notes`. See the
#'   v0.2.0 paper-alignment design doc §3.1 for the full structure.
#'
#'   **Percentile units.** `card$coefficient$surface_percentile` and
#'   `card$panel$surface_percentile` are reported on the 0-100 scale
#'   (e.g., `46.3` means the 46th percentile). The underlying
#'   `position_on_surface()` returns `percentile` on the 0-1 fraction
#'   scale; `grass_report()` multiplies by 100 to match the paper's
#'   prose convention. The `print` and `format` methods use ordinal
#'   notation ("46th percentile").
#' @export
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rbinom(1000, 1, 0.3), nrow = 200, ncol = 5)
#' card <- grass_report(ratings = Y)
#' card                       # print
#' summary(card)              # full panel + per-rater
#' as.data.frame(card)        # tidy long-format
grass_report <- function(ratings,
                         axis = c("inter", "intra"),
                         metric = "auto",
                         occasion = NULL,
                         bands = c(0.5, 0.625, 0.75, 0.875, 1.0),
                         band_labels = c("Poor", "Moderate",
                                         "Strong", "Excellent"),
                         delta_thresholds = NULL,
                         bootstrap_B = 1000L,
                         bootstrap_delta_B = 0L,
                         verbose = FALSE,
                         ...) {
  axis <- match.arg(axis)
  call <- match.call()

  if (isTRUE(verbose)) message("grass_report: normalizing ratings.")
  Y        <- normalize_ratings(ratings)
  pi_hat   <- mean(Y)
  tau2_hat <- compute_tau2_hat(Y)
  k        <- ncol(Y)
  N        <- nrow(Y)

  notes <- character(0L)
  if (k < 2L) {
    stop("`ratings` must have at least 2 rater columns (k >= 2). Got k = ",
         k, ".", call. = FALSE)
  }

  # Resolve delta_thresholds. NULL (default) -> per-(k, N) lookup against
  # the calibration grid bundled in sysdata.rda; non-NULL -> honor user
  # supplied pair. Source is recorded on the card so it can be surfaced
  # in print() and audit output.
  thresholds_source <- "user_supplied"
  thresholds_note   <- ""
  if (is.null(delta_thresholds)) {
    lk <- lookup_delta_thresholds(k = k, N = N)
    delta_thresholds  <- lk$thresholds
    thresholds_source <- lk$source
    thresholds_note   <- lk$note
    if (nzchar(thresholds_note)) {
      notes <- c(notes, thresholds_note)
    }
  }
  if (N < 10L) {
    warning("N = ", N, " is below 10; inference unreliable. ",
            "Reference percentiles and bootstrap CIs will be very noisy.",
            call. = FALSE)
  } else if (N < 30L) {
    notes <- c(notes,
               sprintf("N = %d is small (< 30); reference percentiles are sensitive to sampling noise.",
                       N))
  }
  if (identical(axis, "intra")) {
    notes <- c(notes,
               paste0("Intra-rater axis: surface percentile uses the bundled ",
                      "inter-rater diagonal calibration as an approximation. ",
                      "The dedicated intra-axis calibration cube (1,500-rep ",
                      "intra-rater extension from paper2/code/31_intra_dgp.R) ",
                      "is queued for a follow-on data release. Treat this ",
                      "Report Card's intra-axis percentile as approximate."))
  }

  # ---- Compute the observed panel ---------------------------------------
  if (isTRUE(verbose)) message("grass_report: computing observed-metric panel.")
  panel_obs <- compute_panel(Y, axis = axis, occasion = occasion)

  # Map compute_panel() short names to the surface-metric names accepted by
  # position_on_surface(). compute_panel() returns `ac1` (panel-friendly
  # short name), but position_on_surface() expects `mean_ac1`. Cohen's
  # `kappa` (k = 2 only) has no calibrated surface and is dropped from the
  # panel-positioning machinery.
  surface_metric_for <- c(
    pabak          = "pabak",
    ac1            = "mean_ac1",
    fleiss_kappa   = "fleiss_kappa",
    krippendorff_a = "krippendorff_a",
    icc            = "icc"
  )
  panel_keep <- intersect(names(panel_obs), names(surface_metric_for))
  panel_obs  <- panel_obs[panel_keep]
  names(panel_obs) <- unname(surface_metric_for[panel_keep])

  # ---- Position each panel coefficient on its reference surface ---------
  if (isTRUE(verbose)) message("grass_report: positioning panel on reference surfaces.")
  positions <- list()
  for (m in names(panel_obs)) {
    positions[[m]] <- position_on_surface(
      ratings     = Y,
      metric      = m,
      bands       = bands,
      band_labels = band_labels,
      ...
    )
    notes <- c(notes, positions[[m]]$notes %||% character(0L))
  }

  # ---- Pick the primary coefficient -------------------------------------
  if (identical(metric, "auto")) {
    primary_short <- pick_primary_coefficient(k, pi_hat, axis)
    primary <- unname(surface_metric_for[primary_short])
    if (is.na(primary)) primary <- primary_short
  } else {
    # Allow the user to pass either the panel short name or surface name.
    primary <- if (metric %in% names(surface_metric_for)) {
      unname(surface_metric_for[metric])
    } else {
      metric
    }
  }
  if (!primary %in% names(panel_obs)) {
    stop("Selected primary coefficient `", primary, "` is not in the ",
         "computed panel (", paste(shQuote(names(panel_obs)),
                                   collapse = ", "), "). ",
         "Either pass `metric = \"auto\"` or pick from the panel.",
         call. = FALSE)
  }

  # ---- Cross-coefficient delta_hat --------------------------------------
  if (isTRUE(verbose)) message("grass_report: computing cross-coefficient delta_hat.")
  asym <- check_asymmetry(
    ratings          = Y,
    axis             = axis,
    occasion         = occasion,
    delta_thresholds = delta_thresholds
  )
  delta_hat_pp <- asym$delta_hat
  flag         <- asym$flag

  # ---- Subject-bootstrap CI on delta_hat + tier-membership probabilities --
  # When bootstrap_delta_B > 0, resample subjects with replacement and
  # recompute delta_hat per resample. The bootstrap distribution gives
  # a 95% CI on delta_hat and an empirical tier-membership probability
  # vector P(aligned, caution, divergent). At small N the tier
  # probabilities spread across multiple tiers; at large N they
  # concentrate on the point-estimate tier. This is a frequentist
  # subject-bootstrap, not a Bayesian posterior.
  delta_boot_ci      <- NULL
  tier_probabilities <- NULL
  if (is.numeric(bootstrap_delta_B) && bootstrap_delta_B >= 50L) {
    if (isTRUE(verbose)) {
      message("grass_report: bootstrapping delta_hat (B = ",
              bootstrap_delta_B, ").")
    }
    delta_boot <- numeric(bootstrap_delta_B)
    for (b in seq_len(bootstrap_delta_B)) {
      idx  <- sample.int(N, N, replace = TRUE)
      Y_b  <- Y[idx, , drop = FALSE]
      asym_b <- tryCatch(
        check_asymmetry(
          ratings          = Y_b,
          axis             = axis,
          occasion         = occasion,
          delta_thresholds = delta_thresholds
        ),
        error = function(e) list(delta_hat = NA_real_))
      delta_boot[b] <- asym_b$delta_hat
    }
    delta_boot <- delta_boot[is.finite(delta_boot)]
    if (length(delta_boot) >= 50L) {
      delta_boot_ci <- as.numeric(stats::quantile(delta_boot,
                                                   probs = c(0.025, 0.975),
                                                   na.rm = TRUE))
      tier_probabilities <- c(
        aligned   = mean(delta_boot < delta_thresholds[1]),
        caution   = mean(delta_boot >= delta_thresholds[1] &
                         delta_boot <  delta_thresholds[2]),
        divergent = mean(delta_boot >= delta_thresholds[2])
      )
    } else {
      notes <- c(notes,
        sprintf("delta_hat bootstrap returned %d valid replicates (< 50); CI / tier probabilities suppressed.",
                length(delta_boot)))
    }
  }

  # ---- Divergent-branch per-rater latent-class fit ----------------------
  per_rater <- NULL
  pairwise  <- NULL
  if (identical(flag, "divergent")) {
    if (isTRUE(verbose)) {
      message("grass_report: divergent flag - running latent-class fit ",
              "(B = ", bootstrap_B, ").")
    }
    fit_method <- if (k == 2L) "hui_walter" else "dawid_skene_em"
    lc <- latent_class_fit(
      ratings = Y,
      B       = bootstrap_B,
      method  = fit_method
    )
    per_rater <- lc$per_rater

    # Recommended primary deliverable under divergent (paper §3.3): the
    # pairwise PABAK matrix on the k = 2 surface plus per-rater pooled-
    # reference Se/Sp from the panel-majority of the other k - 1 raters.
    # The latent-class fit above is kept alongside as the alternative
    # when external orientation is available.
    if (isTRUE(verbose)) {
      message("grass_report: divergent flag - computing pairwise reliability.")
    }
    pairwise <- tryCatch(
      pairwise_agreement(ratings = Y, axis = axis),
      error = function(e) {
        notes <<- c(notes,
          sprintf("pairwise_agreement() failed: %s", conditionMessage(e)))
        NULL
      })
  } else if (identical(flag, "caution")) {
    notes <- c(notes,
      "Caution: cross-coefficient spread is approaching the divergent threshold. Consider running pairwise_agreement(ratings = Y) for a pairwise PABAK matrix and per-rater pooled-reference Se/Sp; see paper §3.3.")
  }

  # ---- Assemble the panel data.frame ------------------------------------
  pct_pp <- vapply(positions,
                   function(p) p$percentile * 100,
                   numeric(1))
  bp_modal <- vapply(
    positions,
    function(p) {
      bp <- p$band_probabilities %||% NA_real_
      if (length(bp) == 0L || all(!is.finite(bp))) NA_real_ else max(bp)
    },
    numeric(1)
  )
  band_vec <- vapply(
    positions,
    function(p) p$modal_band_label %||% NA_character_,
    character(1)
  )
  qual_vec <- vapply(
    positions,
    function(p) p$confidence %||% NA_character_,
    character(1)
  )
  q_hat_vec <- vapply(
    positions,
    function(p) p$q_hat %||% NA_real_,
    numeric(1)
  )
  se_q_hat_vec <- vapply(
    positions,
    function(p) p$se_q_hat %||% NA_real_,
    numeric(1)
  )
  # Surface-envelope clamp flag: TRUE when the observed coefficient fell
  # outside the achievable range of its reference surface, in which case
  # invert_metric_to_q() clamped q_hat to the boundary and percentile
  # landed at 0 or 100. Such coefficients are excluded from delta_hat
  # (see asymmetry.R) and are surfaced in print() output for transparency.
  clamped_vec <- vapply(
    positions,
    function(p) any(grepl("q_hat clamped", p$notes %||% character(0L), fixed = TRUE)),
    logical(1)
  )
  # reference_used: which path produced each coefficient's reference
  # surface (closed-form / fitted-icc / oracle-icc-fallback /
  # oracle-icc-explicit / user-supplied). Surfaced on the card so a
  # practitioner can audit whether their ICC percentile came from the
  # fitted reference (default, GLMM-gap-corrected) or fell back to the
  # oracle reference.
  reference_used_vec <- vapply(
    positions,
    function(p) p$reference_used %||% "closed-form",
    character(1)
  )

  # in_delta_hat = TRUE for agreement-family coefficients (PABAK, AC1, Fleiss,
  # alpha) which contribute to delta_hat under the v0.5.0 scope decision; ICC
  # is always FALSE because its reference surface depends on the full F-shape
  # rather than (q, pi_+) and would inflate delta_hat under F-shape
  # misspecification. See asymmetry.R `.DELTA_AGREEMENT_COEFS` for the source
  # of truth on the coefficient set.
  in_delta_set <- names(panel_obs) %in% c("pabak", "mean_ac1",
                                          "fleiss_kappa", "krippendorff_a")
  panel_df <- data.frame(
    coefficient            = names(panel_obs),
    observed_value         = unlist(panel_obs, use.names = FALSE),
    surface_percentile     = unname(pct_pp),
    band                   = unname(band_vec),
    qualifier              = unname(qual_vec),
    band_probability_modal = unname(bp_modal),
    q_hat                  = unname(q_hat_vec),
    se_q_hat               = unname(se_q_hat_vec),
    clamped                = unname(clamped_vec),
    reference_used         = unname(reference_used_vec),
    in_delta_hat           = unname(in_delta_set),
    stringsAsFactors       = FALSE
  )

  # ---- Reference-surface artifacts (for plot() / audit) -----------------
  # position_on_surface() does not currently surface its 501-point q_grid or
  # the reference curve in its return value, so capture what's analogous
  # (NULL placeholders) per the design doc §3.1 NOTE.
  surface_payload <- list(
    q_grid           = positions[[primary]]$q_grid %||% NULL,
    reference_curves = lapply(positions, function(p) p$ref_curve %||% NULL),
    reference_type   = vapply(
      positions,
      function(p) p$reference_type %||% "closed-form",
      character(1)
    )
  )

  # ---- Assemble the grass_card ------------------------------------------
  primary_pos <- positions[[primary]]
  out <- list(
    sample = list(
      k        = as.integer(k),
      N        = as.integer(N),
      pi_hat   = as.numeric(pi_hat),
      tau2_hat = as.numeric(tau2_hat),
      axis     = axis
    ),
    coefficient = list(
      primary             = primary,
      observed_value      = unname(panel_obs[[primary]]),
      surface_percentile  = unname(primary_pos$percentile) * 100,
      band                = if (identical(flag, "divergent")) "suppressed"
                            else (primary_pos$modal_band_label %||% NA_character_),
      qualifier           = if (identical(flag, "divergent")) NA_character_
                            else (primary_pos$confidence %||% NA_character_)
    ),
    delta = list(
      delta_hat          = as.numeric(delta_hat_pp),
      flag               = flag,
      thresholds         = setNames(as.numeric(delta_thresholds),
                                    c("caution", "divergent")),
      thresholds_source  = thresholds_source,
      delta_hat_ci       = delta_boot_ci,
      tier_probabilities = tier_probabilities,
      bootstrap_B        = as.integer(if (is.numeric(bootstrap_delta_B))
                                        bootstrap_delta_B else 0L)
    ),
    panel     = panel_df,
    per_rater = per_rater,
    pairwise  = pairwise,
    surface   = surface_payload,
    call      = call,
    grass_version = utils::packageVersion("grass"),
    timestamp = Sys.time(),
    inputs = list(
      ratings_dim      = c(N = as.integer(N), k = as.integer(k)),
      axis             = axis,
      delta_thresholds = as.numeric(delta_thresholds),
      bands            = as.numeric(bands),
      band_labels      = as.character(band_labels)
    ),
    notes = unique(notes)
  )
  class(out) <- c("grass_card", "list")
  out
}

#' Estimate prevalence of the positive class from rater data
#'
#' Averages the marginal positive rates of the two raters.
#'
#' @inheritParams grass_compute
#'
#' @return A single numeric in `[0, 1]`.
#' @export
#'
#' @examples
#' tab <- matrix(c(88, 10, 14, 88), nrow = 2,
#'               dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
#' grass_prevalence(tab, format = "matrix")
grass_prevalence <- function(data, format = c("wide", "matrix", "long", "paired"),
                             positive = NULL, ...) {
  format <- match.arg(format)
  norm <- normalize_input(data, format = format, positive = positive, ...)
  tab <- if (!is.null(norm$table)) norm$table else build_table(norm$r1, norm$r2)
  estimate_prevalence(tab)
}

# Internal: average marginal positive rate from a 2x2 table.
estimate_prevalence <- function(tab) {
  N <- sum(tab)
  if (N == 0) return(NA_real_)
  p1 <- (tab[2, 1] + tab[2, 2]) / N  # R1 says 1
  p2 <- (tab[1, 2] + tab[2, 2]) / N  # R2 says 1
  (p1 + p2) / 2
}

# Internal: classify the (PI, BI) regime and attach the structural
# implication — what the mathematics forces on the metrics, stated as fact.
#   - balanced              : both indices small
#   - prevalence-dominated  : PI meaningfully exceeds BI
#   - bias-dominated        : BI meaningfully exceeds PI
#   - mixed                 : both indices non-trivial with neither dominant
#
# Called by the S3 method classify_regime.grass_spec_binary() in
# dispatch.R. Named `classify_regime_binary` (not `classify_regime`) so
# the generic of the same name can live in dispatch.R without shadowing.
classify_regime_binary <- function(pi, bi, small = 0.10, gap = 0.10) {
  pi <- unname(pi); bi <- unname(bi)
  if (is.na(pi) || is.na(bi)) {
    return(list(regime = NA_character_, note = NA_character_))
  }
  if (max(pi, bi) < small) {
    return(list(
      regime = "balanced",
      note = paste0(
        "Prevalence and bias indices are both small, so kappa, PABAK, and ",
        "AC1 are algebraically close. Substantial disagreement among the ",
        "three metrics in this regime points to model or data issues rather ",
        "than rater behaviour."
      )
    ))
  }
  if (pi - bi > gap) {
    return(list(
      regime = "prevalence-dominated",
      note = paste0(
        "Expected agreement under marginal independence is high, so kappa ",
        "is algebraically suppressed even when raters are accurate on the ",
        "minority class. PABAK and AC1 use different denominators and are ",
        "less sensitive to this effect; a gap between kappa and PABAK/AC1 ",
        "is a structural property of the prevalence."
      )
    ))
  }
  if (bi - pi > gap) {
    return(list(
      regime = "bias-dominated",
      note = paste0(
        "The two raters endorse the positive class at different marginal ",
        "rates. AC1 conditions on this disagreement and will trail kappa ",
        "and PABAK; the appropriate response is rater calibration, not ",
        "more training on individual cases."
      )
    ))
  }
  list(
    regime = "mixed",
    note = paste0(
      "Prevalence and bias effects are both present. No single metric ",
      "cleanly captures agreement here; report kappa, PABAK, and AC1 ",
      "alongside PI and BI so the reader sees the full structure."
    )
  )
}

# Internal: signed distance from observed metric to the reference curve at the
# evaluation prevalence. Positive = above the reference, negative = below.
distance_to_reference <- function(metrics, reference) {
  v <- metrics$values
  ref <- reference$reference
  observed <- c(kappa = unname(v["kappa"]),
                PABAK = unname(v["PABAK"]),
                AC1   = unname(v["AC1"]))
  data.frame(
    metric    = ref$metric,
    observed  = observed[ref$metric],
    reference = ref$reference,
    distance  = observed[ref$metric] - ref$reference,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
