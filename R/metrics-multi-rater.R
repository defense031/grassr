# Multi-rater observed metrics on an N x k binary rating matrix.
#
# All `obs_*` functions take Y (N x k integer matrix, values in {0L, 1L}) and
# return a single numeric. Implementations are hand-rolled (transparent
# formulas) and ported verbatim from `paper2/code/02_metrics_observed.R`,
# where they are validated against `irr` / `irrCAC`.
#
# Notation:
#   N    = nrow(Y)               number of subjects
#   k    = ncol(Y)                raters per subject
#   r_i  = rowSums(Y)             positives among k raters for subject i
#   p_+  = mean(Y)                overall positive rate
#
# All entry points here are internal (not exported). Phase 1A wiring; Phase
# 2+ adds the public `position_on_surface(ratings = Y, metric = ...)` and
# `grass_report(ratings = Y)` paths that consume them.

.check_Y_multi <- function(Y) {
  if (!is.matrix(Y)) {
    stop("`Y` must be an integer matrix (use normalize_ratings()).", call. = FALSE)
  }
  if (!all(Y %in% c(0L, 1L))) {
    stop("`Y` must contain only 0/1 values.", call. = FALSE)
  }
  if (ncol(Y) < 2) {
    stop("`Y` must have at least 2 rater columns (k >= 2).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Fleiss' kappa (multi-rater nominal, binary case)
#'
#' `r_i = rowSums(Y)`. Per-subject agreement
#' `P_i = (r_i^2 + (k - r_i)^2 - k) / (k * (k - 1))`,
#' marginal `p_1 = mean(Y)`, expected agreement
#' `P_e = p_0^2 + p_1^2`, kappa = (P_bar - P_e) / (1 - P_e).
#'
#' @param Y N x k integer matrix in \{0L, 1L\}.
#' @return Single numeric.
#' @keywords internal
#' @noRd
obs_fleiss_kappa <- function(Y) {
  .check_Y_multi(Y)
  k <- ncol(Y)
  r <- rowSums(Y)
  P_i   <- (r^2 + (k - r)^2 - k) / (k * (k - 1))
  P_bar <- mean(P_i)
  p1    <- mean(Y)
  P_e   <- p1^2 + (1 - p1)^2
  (P_bar - P_e) / (1 - P_e)
}

#' Mean pairwise PABAK
#'
#' For each rater pair (j, j'): `PABAK_{j,j'} = 2 * mean(Y[, j] == Y[, j']) - 1`.
#' Mean over all `choose(k, 2)` pairs. Vectorized via per-subject agreeing-pair
#' count: subject i contributes
#' `(r_i * (r_i - 1) + (k - r_i) * (k - r_i - 1)) / 2` agreeing pairs out of
#' `choose(k, 2)`.
#'
#' @param Y N x k integer matrix in \{0L, 1L\}.
#' @return Single numeric.
#' @keywords internal
#' @noRd
obs_mean_pairwise_pabak <- function(Y) {
  .check_Y_multi(Y)
  k <- ncol(Y)
  r <- rowSums(Y)
  agree_pairs <- (r * (r - 1) + (k - r) * (k - r - 1)) / 2
  total_pairs <- k * (k - 1) / 2
  p_a <- mean(agree_pairs) / total_pairs
  2 * p_a - 1
}

#' Mean pairwise AC1 (Gwet)
#'
#' Per pair: `p_e = 2 * pi * (1 - pi)` where `pi` is the average positive
#' rate of the two raters. `AC1 = (p_a - p_e) / (1 - p_e)`. Averaged over
#' pairs. Pairwise reduction kept symmetric with PABAK.
#'
#' @param Y N x k integer matrix in \{0L, 1L\}.
#' @return Single numeric.
#' @keywords internal
#' @noRd
obs_mean_pairwise_ac1 <- function(Y) {
  .check_Y_multi(Y)
  k <- ncol(Y)
  marginals <- colMeans(Y)
  pairs <- utils::combn(k, 2)
  vals <- apply(pairs, 2, function(jj) {
    p_a <- mean(Y[, jj[1]] == Y[, jj[2]])
    pi_bar <- mean(marginals[jj])
    p_e <- 2 * pi_bar * (1 - pi_bar)
    (p_a - p_e) / (1 - p_e)
  })
  mean(vals)
}

#' Krippendorff's alpha (nominal, binary, fully crossed)
#'
#' Fully-crossed design (`m_i = k` for all i). Hayes-Krippendorff coincidence-
#' matrix reduction:
#' `o_{01} = 2 * sum_i r_i * (k - r_i) / (k - 1)`,
#' `n_0 = N*k - sum(r_i)`, `n_1 = sum(r_i)`, `n_total = N*k`,
#' `alpha = 1 - (n_total - 1) * o_{01} / (2 * n_0 * n_1)`.
#'
#' Returns `NA_real_` (with a `note` attribute) if the design is degenerate
#' (no positives or no negatives in `Y`).
#'
#' @param Y N x k integer matrix in \{0L, 1L\}.
#' @return Single numeric.
#' @keywords internal
#' @noRd
obs_krippendorff_alpha <- function(Y) {
  .check_Y_multi(Y)
  k <- ncol(Y)
  N <- nrow(Y)
  r <- rowSums(Y)
  n_total <- N * k
  n1 <- sum(r)
  n0 <- n_total - n1
  if (n0 == 0 || n1 == 0) {
    out <- NA_real_
    attr(out, "note") <- "Krippendorff alpha undefined: all-0 or all-1 ratings."
    return(out)
  }
  o_01 <- 2 * sum(r * (k - r)) / (k - 1)
  1 - (n_total - 1) * o_01 / (2 * n0 * n1)
}

#' Logit-mixed ICC(1,1) via lme4::glmer
#'
#' Random-intercept logistic GLMM:
#' `logit(P(Y_ij = 1)) = beta_0 + u_i`, `u_i ~ N(0, tau^2)`.
#' Standard logistic-latent convention (Nakagawa-Schielzeth 2013):
#' `ICC = tau^2 / (tau^2 + pi^2 / 3)`.
#'
#' If `lme4` is not installed, returns `NA_real_` with
#' `attr(., "note") = "lme4 not available"`. If `glmer` fails to converge,
#' returns `NA_real_` with a "glmer failed" note.
#'
#' @param Y N x k integer matrix in \{0L, 1L\}.
#' @return Single numeric or `NA_real_` with a `note` attribute.
#' @keywords internal
#' @noRd
obs_icc_glmer <- function(Y) {
  .check_Y_multi(Y)
  if (!requireNamespace("lme4", quietly = TRUE)) {
    out <- NA_real_
    attr(out, "note") <- "lme4 not available"
    return(out)
  }
  N <- nrow(Y)
  k <- ncol(Y)
  df <- data.frame(
    y = as.integer(as.vector(Y)),               # column-major: subjects rater1, then rater2, ...
    subject = factor(rep(seq_len(N), times = k))
  )
  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      lme4::glmer(y ~ 1 + (1 | subject), data = df, family = stats::binomial)
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    out <- NA_real_
    attr(out, "note") <- "glmer failed"
    return(out)
  }
  vc <- as.data.frame(lme4::VarCorr(fit))
  tau2 <- vc$vcov[vc$grp == "subject"]
  if (length(tau2) == 0L || !is.finite(tau2)) {
    out <- NA_real_
    attr(out, "note") <- "glmer returned no subject variance"
    return(out)
  }
  tau2 / (tau2 + pi^2 / 3)
}

#' Dispatch a single observed metric on an N x k matrix
#'
#' Convenience switch over the panel. `metric` is one of
#' `"pabak"`, `"fleiss_kappa"`, `"mean_ac1"`, `"krippendorff_a"`, `"icc"`.
#' At `k = 2`, `"pabak"` calls `obs_mean_pairwise_pabak()`, which reduces
#' correctly to the two-rater PABAK.
#'
#' @param metric Single string, one of the supported metric names.
#' @param Y N x k integer matrix in \{0L, 1L\}.
#' @return Single numeric.
#' @keywords internal
#' @noRd
compute_observed <- function(metric, Y) {
  .check_Y_multi(Y)
  metric <- match.arg(metric,
                      c("pabak", "fleiss_kappa", "mean_ac1",
                        "krippendorff_a", "icc"))
  switch(metric,
    pabak          = obs_mean_pairwise_pabak(Y),
    fleiss_kappa   = obs_fleiss_kappa(Y),
    mean_ac1       = obs_mean_pairwise_ac1(Y),
    krippendorff_a = obs_krippendorff_alpha(Y),
    icc            = obs_icc_glmer(Y)
  )
}
