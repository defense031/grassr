# Multi-rater helpers: input normalization, the observed coefficient panel,
# and Table-2 primary-coefficient pick. All internal (Phase 1A wiring).

#' Normalize ratings to a canonical N x k integer matrix
#'
#' Accepts:
#'   * `N x k` integer / numeric matrix with values in {0, 1}
#'   * `N x k` data.frame whose columns are all binary (0/1, logical, or
#'     factor with two levels)
#'   * a list of two equal-length 0/1 vectors (k = 2 paired form)
#'
#' Always returns an `N x k` integer matrix with values in `{0L, 1L}`.
#' Rejects non-binary input, NA values, k < 2, N < 1, or 2x2 count tables
#' (the latter requires the legacy `normalize_input(format = "matrix")`
#' path; we cannot recover the underlying N x 2 ratings without ordering).
#'
#' @param ratings User input.
#' @return `N x k` integer matrix.
#' @keywords internal
#' @noRd
normalize_ratings <- function(ratings) {
  if (is.null(ratings)) {
    stop("`ratings` is NULL. Supply an N x k binary matrix, a data.frame ",
         "with k rater columns, or a list of two equal-length 0/1 vectors.",
         call. = FALSE)
  }

  # Reject 2x2 count tables — ambiguous reconstruction.
  if (is.matrix(ratings) && all(dim(ratings) == c(2, 2)) &&
      !is.null(rownames(ratings)) &&
      identical(rownames(ratings), c("0", "1")) &&
      identical(colnames(ratings), c("0", "1"))) {
    stop("`ratings` looks like a 2x2 count table. The multi-rater path ",
         "needs the underlying N x k subject-by-rater matrix, not counts. ",
         "Use the legacy 2x2 entry point if you only have cell counts.",
         call. = FALSE)
  }

  # Branch 1: list of two vectors (k = 2 paired).
  if (is.list(ratings) && !is.data.frame(ratings) && length(ratings) == 2L) {
    r1 <- ratings[[1]]
    r2 <- ratings[[2]]
    if (length(r1) != length(r2)) {
      stop("`ratings` is a length-2 list but the two vectors have unequal ",
           "lengths (", length(r1), " vs ", length(r2), ").", call. = FALSE)
    }
    Y <- cbind(r1, r2)
    return(.coerce_to_binary_matrix(Y))
  }

  # Branch 2: data.frame with k rater columns.
  if (is.data.frame(ratings)) {
    cols <- lapply(ratings, function(x) {
      if (is.logical(x)) return(as.integer(x))
      if (is.factor(x)) {
        lv <- levels(droplevels(x))
        if (length(lv) > 2) {
          stop("data.frame column has more than 2 levels: ",
               paste(shQuote(lv), collapse = ", "),
               ". Recode to binary first.", call. = FALSE)
        }
        # Two-level factor: positive = level "1" if present, else second level.
        pos <- if ("1" %in% lv) "1" else if ("TRUE" %in% lv) "TRUE" else lv[length(lv)]
        return(as.integer(as.character(x) == pos))
      }
      if (is.numeric(x)) return(as.integer(x))
      stop("Unsupported data.frame column type: ",
           paste(class(x), collapse = "/"),
           ". Use logical, integer 0/1, or 2-level factor.",
           call. = FALSE)
    })
    Y <- do.call(cbind, cols)
    return(.coerce_to_binary_matrix(Y))
  }

  # Branch 3: matrix.
  if (is.matrix(ratings)) {
    return(.coerce_to_binary_matrix(ratings))
  }

  stop("`ratings` has unsupported type: ",
       paste(class(ratings), collapse = "/"),
       ". Use an N x k matrix, a data.frame with k rater columns, or a ",
       "list of two equal-length 0/1 vectors.", call. = FALSE)
}

# Internal: coerce a numeric / logical / integer matrix to N x k integer
# {0L, 1L} after validation.
.coerce_to_binary_matrix <- function(Y) {
  if (!is.matrix(Y)) {
    stop("Internal: expected a matrix here.", call. = FALSE)
  }
  if (anyNA(Y)) {
    stop("`ratings` contains NA values; drop or impute before calling.",
         call. = FALSE)
  }
  if (is.logical(Y)) Y <- Y * 1L
  if (!is.numeric(Y)) {
    stop("`ratings` must be numeric / logical after coercion (got ",
         paste(class(Y), collapse = "/"), ").", call. = FALSE)
  }
  vals <- unique(as.vector(Y))
  if (!all(vals %in% c(0, 1))) {
    stop("`ratings` must contain only 0/1 values. Saw: ",
         paste(sort(vals), collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(Y) < 1L) {
    stop("`ratings` must have at least one subject (N >= 1).", call. = FALSE)
  }
  if (ncol(Y) < 2L) {
    stop("`ratings` must have at least 2 rater columns (k >= 2). Got ",
         ncol(Y), ".", call. = FALSE)
  }
  storage.mode(Y) <- "integer"
  Y
}

#' Estimate subject-level prevalence variance from the rating matrix
#'
#' Method-of-moments on per-subject observed positive rates:
#'
#'   sigma2_p_hat = max(0, Var_i[p_hat_i] - pi_hat (1 - pi_hat) / k) * k / (k - 1)
#'
#' where p_hat_i = (sum_j Y_ij) / k is subject i's observed positive-call
#' rate. The subtracted term is the per-subject Bernoulli sampling
#' variance; what remains is the across-subject prevalence variance on
#' the prevalence scale (interpretable directly: sqrt(tau2_hat) is the
#' standard deviation of true positive rates across subjects, in the same
#' units as pi_hat).
#'
#' Returns 0 when k < 2, when N < 2, or when sample variance falls below
#' the binomial floor (degenerate at small k where per-subject rates take
#' too few distinct values for the method-of-moments correction to
#' separate prevalence variance from sampling noise).
#'
#' Note: F_key parameterization in the bundled reference surfaces uses
#' logit-scale tau2 (logit-normal F shape). The ICC reference-curve
#' selection path in `position_on_surface()` calls glmer separately to
#' fit (mu, tau2) on the logit scale; this prevalence-scale estimate is
#' the practitioner-visible summary surfaced in the Report Card sample
#' field, not the F_key selector.
#'
#' @param Y N x k binary matrix.
#' @return Single non-negative numeric: estimated prevalence-scale Var[p_i].
#' @keywords internal
#' @noRd
compute_tau2_hat <- function(Y) {
  k <- ncol(Y)
  N <- nrow(Y)
  if (k < 2L || N < 2L) return(0)
  pi_hat <- mean(Y)
  p_i <- rowMeans(Y)
  sample_var <- stats::var(p_i)
  binomial_floor <- pi_hat * (1 - pi_hat) / k
  max(0, sample_var - binomial_floor) * k / (k - 1)
}

#' Compute the observed coefficient panel
#'
#' At `k = 2`: returns `pabak`, `ac1`, `kappa` (Cohen's), `krippendorff_a`.
#' Cohen's kappa, PABAK, and AC1 at `k = 2` are derived from the existing
#' `compute_agreement_metrics()` 2x2 path so the package never has two
#' diverging implementations of the same coefficient.
#'
#' At `k >= 3`: returns `pabak`, `ac1`, `fleiss_kappa`, `krippendorff_a`,
#' `icc`. (`icc` may be `NA_real_` if `lme4` is not installed or `glmer`
#' fails; an `attr(., "note")` flags why.)
#'
#' @param ratings N x k binary matrix / data.frame / k=2 list.
#' @param axis "inter" (default) or "intra". Phase 1A treats both the same;
#'   the intra-axis cube-reshape lands in Phase 4+.
#' @param occasion Reserved for axis = "intra"; ignored in Phase 1A.
#' @return Named list of observed metric values.
#' @keywords internal
#' @noRd
compute_panel <- function(ratings, axis = "inter", occasion = NULL) {
  Y <- normalize_ratings(ratings)
  axis <- match.arg(axis, c("inter", "intra"))
  k <- ncol(Y)

  if (k == 2) {
    tab <- build_table(Y[, 1], Y[, 2])
    base <- compute_agreement_metrics(tab)
    out <- list(
      pabak          = unname(base[["PABAK"]]),
      ac1            = unname(base[["AC1"]]),
      kappa          = unname(base[["kappa"]]),
      krippendorff_a = obs_krippendorff_alpha(Y)
    )
  } else {
    out <- list(
      pabak          = obs_mean_pairwise_pabak(Y),
      ac1            = obs_mean_pairwise_ac1(Y),
      fleiss_kappa   = obs_fleiss_kappa(Y),
      krippendorff_a = obs_krippendorff_alpha(Y),
      icc            = obs_icc_glmer(Y)
    )
  }
  out
}

#' Pick the primary coefficient (Table 2 of the paper)
#'
#' For `axis = "intra"` the primary is always ICC (Koo & Li 2016).
#'
#' For `axis = "inter"`:
#'   * `k == 2`: PABAK at balanced prevalence; AC1 at extreme prevalence
#'     (`pi_hat < 0.20` or `> 0.80`).
#'   * `k >= 3`: PABAK at balanced; AC1 at extreme prevalence
#'     (`pi_hat < 0.15` or `> 0.85`), aligned with App F crossover evidence.
#'
#' Hard-coded thresholds are intentional. Override for advanced use by
#' calling `pick_primary_coefficient()` directly with a tweaked branch is
#' not supported in v0.2.0.
#'
#' @param k Integer number of raters.
#' @param pi_hat Observed positive rate in [0, 1].
#' @param axis "inter" (default) or "intra".
#' @return Single string in
#'   `c("pabak", "ac1", "fleiss_kappa", "krippendorff_a", "icc")`.
#' @keywords internal
#' @noRd
pick_primary_coefficient <- function(k, pi_hat, axis = "inter") {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2) {
    stop("`k` must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(pi_hat) || length(pi_hat) != 1L ||
      !is.finite(pi_hat) || pi_hat < 0 || pi_hat > 1) {
    stop("`pi_hat` must be a single numeric in [0, 1].", call. = FALSE)
  }
  axis <- match.arg(axis, c("inter", "intra"))

  if (axis == "intra") return("icc")

  if (k == 2) {
    return(if (pi_hat < 0.20 || pi_hat > 0.80) "ac1" else "pabak")
  }
  # k >= 3
  if (pi_hat < 0.15 || pi_hat > 0.85) "ac1" else "pabak"
}
