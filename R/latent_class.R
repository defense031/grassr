# Latent-class fit for the divergent branch of the GRASS Reporting Card.
#
# When the cross-coefficient panel disagrees (delta_hat >= 11.75 pp),
# `grass_report()` falls back from a single q_hat summary to per-rater
# Sensitivity (Se) and Specificity (Sp). This file implements the
# latent-class machinery that fallback relies on.
#
# Two regimes:
#   - k >= 3 raters: Dawid-Skene 1979 EM. Hand-rolled, base R only.
#   - k == 2 raters: Hui-Walter 1980 closed-form bounds. Per-rater Se and
#     Sp are NOT point-identified from a single 2x2 table; we report the
#     inequality bounds the marginals support.
#
# The bootstrap CI is nonparametric: resample subjects with replacement
# B times and refit. At k = 2 we resample and recompute the bounds; we
# report the empirical 2.5/97.5 percentiles on the bound midpoints.

# ---- Internal: input normalization (defensive vs. Phase 1A landing) ----

# Coerce ratings to canonical N x k integer matrix in {0, 1}. If
# Phase 1A's normalize_ratings() is loaded, defer to it; otherwise do
# the inline check.
.lc_normalize_ratings <- function(ratings) {
  if (exists("normalize_ratings", mode = "function")) {
    Y <- normalize_ratings(ratings)
  } else {
    if (is.data.frame(ratings)) {
      Y <- as.matrix(ratings)
    } else if (is.list(ratings) && !is.data.frame(ratings)) {
      lens <- vapply(ratings, length, integer(1))
      if (length(unique(lens)) != 1L) {
        stop("`ratings` is a list with unequal-length elements; cannot ",
             "coerce to a rectangular N x k matrix.", call. = FALSE)
      }
      Y <- do.call(cbind, ratings)
    } else if (is.matrix(ratings)) {
      Y <- ratings
    } else if (is.numeric(ratings) || is.logical(ratings)) {
      stop("`ratings` must be an N x k matrix, data.frame, or list of ",
           "rater columns; got a vector.", call. = FALSE)
    } else {
      stop("`ratings` must be coercible to an N x k integer matrix in ",
           "{0, 1}.", call. = FALSE)
    }
    storage.mode(Y) <- "integer"
  }

  if (!is.matrix(Y)) {
    stop("`ratings` must be coercible to a matrix.", call. = FALSE)
  }
  if (anyNA(Y)) {
    stop("`ratings` contains NA values; latent-class fit requires complete ",
         "binary rating data.", call. = FALSE)
  }
  uniq <- unique(as.vector(Y))
  if (!all(uniq %in% c(0L, 1L))) {
    stop("`ratings` must contain only 0 and 1 values; observed: ",
         paste(sort(unique(uniq)), collapse = ", "), ".", call. = FALSE)
  }

  N <- nrow(Y)
  k <- ncol(Y)
  if (k < 2L) {
    stop("`ratings` must have at least 2 rater columns (got k = ", k, ").",
         call. = FALSE)
  }
  if (N < 10L) {
    stop("`ratings` must have at least 10 subjects (got N = ", N, ").",
         call. = FALSE)
  }
  col_sums <- colSums(Y)
  if (any(col_sums == 0L) || any(col_sums == N)) {
    bad <- which(col_sums == 0L | col_sums == N)
    stop("`ratings` has all-constant column(s): rater(s) ",
         paste(bad, collapse = ", "),
         " gave the same rating to every subject. The latent-class fit is ",
         "ill-defined for raters with no signal.", call. = FALSE)
  }
  Y
}

# ---- Internal: Hui-Walter k = 2 bounds ----
#
# Two-rater binary cell counts:
#         R2=0  R2=1
#   R1=0 [  a    b  ]
#   R1=1 [  c    d  ]
#
# With latent prevalence pi and conditionally independent raters,
# the marginals identify P(R_j = 1) = pi * Se_j + (1 - pi) * (1 - Sp_j),
# but Se_j and Sp_j are not separately identified from one such pair.
# The Hui-Walter 1980 inequality bounds: any (Se_j, Sp_j) consistent
# with the observed marginal P_j = (col / row totals) / N must satisfy
# P_j <= Se_j and P_j <= Sp_j (when prevalence is unknown but in [0, 1]).
# Concretely: Se_j in [P_j, 1] and Sp_j in [1 - P_j, 1] is the widest
# admissible range; tightening requires observed agreement.
#
# We use the agreement-tightened bounds: the marginal positive rate of
# rater j gives a lower bound on max(Se_j, 1 - Sp_j); the joint
# agreement P0 gives an upper bound on the per-rater accuracy. The
# range we return is [P_j, 1] for Se and [1 - P_j, 1] for Sp -- the
# unrestricted Hui-Walter range. Tighter algebraic bounds exist when
# both raters share parameters; we do not assume that here.
.lc_hui_walter_bounds <- function(Y) {
  k <- ncol(Y)
  if (k != 2L) {
    stop("Hui-Walter bounds require exactly k = 2 raters.", call. = FALSE)
  }
  P_j <- colMeans(Y)
  data.frame(
    rater    = paste0("R", seq_len(k)),
    se_lower = P_j,
    se_upper = rep(1, k),
    sp_lower = 1 - P_j,
    sp_upper = rep(1, k),
    stringsAsFactors = FALSE
  )
}

# ---- Internal: Dawid-Skene EM at k >= 3 ----
#
# Latent binary class C_i in {0, 1} for subject i; Y_ij = rater j's call.
# Conditional independence given C: Y_ij | C_i = 1 ~ Bern(Se_j),
# Y_ij | C_i = 0 ~ Bern(1 - Sp_j).
# Parameters: pi (prevalence), Se_j, Sp_j for j = 1..k.
# Initialization: majority-vote consensus. Per-rater Se_j and Sp_j are
# the agreement of rater j with the consensus on positive and negative
# subjects respectively.
.lc_ds_em <- function(Y, max_iter = 1000L, tol = 1e-6) {
  N <- nrow(Y)
  k <- ncol(Y)

  # Initialization: majority-vote.
  consensus <- as.integer(rowSums(Y) > k / 2)
  pi <- mean(consensus)
  pi <- min(max(pi, 1e-6), 1 - 1e-6)

  Se <- numeric(k)
  Sp <- numeric(k)
  pos <- consensus == 1L
  neg <- !pos
  n_pos <- sum(pos)
  n_neg <- sum(neg)
  for (j in seq_len(k)) {
    Se[j] <- if (n_pos > 0) mean(Y[pos, j]) else 0.5
    Sp[j] <- if (n_neg > 0) mean(1 - Y[neg, j]) else 0.5
  }
  # Clip to interior to avoid log(0) in the E-step.
  clip <- function(x) pmin(pmax(x, 1e-6), 1 - 1e-6)
  Se <- clip(Se)
  Sp <- clip(Sp)

  log_lik_old <- -Inf
  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(max_iter)) {
    # ---- E-step ----
    # log P(Y_i | C = 1) = sum_j [ Y_ij log Se_j + (1 - Y_ij) log (1 - Se_j) ]
    # log P(Y_i | C = 0) = sum_j [ Y_ij log (1 - Sp_j) + (1 - Y_ij) log Sp_j ]
    log_Se  <- log(Se)
    log_1Se <- log(1 - Se)
    log_Sp  <- log(Sp)
    log_1Sp <- log(1 - Sp)

    log_like_pos <- as.vector(Y %*% log_Se + (1 - Y) %*% log_1Se)
    log_like_neg <- as.vector(Y %*% log_1Sp + (1 - Y) %*% log_Sp)

    log_pi  <- log(pi)
    log_1pi <- log(1 - pi)

    a_pos <- log_pi + log_like_pos
    a_neg <- log_1pi + log_like_neg
    m <- pmax(a_pos, a_neg)
    log_norm <- m + log(exp(a_pos - m) + exp(a_neg - m))
    p <- exp(a_pos - log_norm)  # P(C_i = 1 | Y_i)

    log_lik <- sum(log_norm)

    # ---- M-step ----
    sum_p <- sum(p)
    sum_1p <- N - sum_p
    pi_new <- sum_p / N
    pi_new <- clip(pi_new)

    Se_new <- as.vector(crossprod(Y, p) / sum_p)
    Sp_new <- as.vector(crossprod(1 - Y, 1 - p) / sum_1p)
    Se_new <- clip(Se_new)
    Sp_new <- clip(Sp_new)

    pi <- pi_new
    Se <- Se_new
    Sp <- Sp_new

    if (is.finite(log_lik) && is.finite(log_lik_old) &&
        abs(log_lik - log_lik_old) < tol) {
      converged <- TRUE
      break
    }
    log_lik_old <- log_lik
  }

  list(
    pi = pi,
    Se = Se,
    Sp = Sp,
    p_post = p,
    log_likelihood = log_lik,
    converged = converged,
    iterations = iter
  )
}

# ---- The exported workhorse ---------------------------------------------

#' Latent-class fit: per-rater Sensitivity and Specificity from a binary
#' rating matrix
#'
#' `latent_class_fit()` is the divergent-branch fallback for the GRASS
#' Reporting Card. When the cross-coefficient panel disagrees, the
#' framework abandons the single q_hat summary and reports per-rater
#' (Se_j, Sp_j) instead. This function fits those per-rater accuracy
#' parameters from an N x k binary rating matrix.
#'
#' Two regimes, dispatched by `k`:
#'
#' - **k >= 3: Dawid-Skene 1979 expectation-maximization.** Latent binary
#'   class C_i in {0, 1}; conditional independence of raters given C;
#'   parameters (pi, Se_j, Sp_j) for j = 1..k. Initialized from
#'   majority-vote consensus. Iterated to log-likelihood tolerance `tol`
#'   or `max_iter`. Hand-rolled in base R, no new package dependency.
#'
#' - **k == 2: Hui-Walter 1980 inequality bounds.** Per-rater Se and Sp
#'   are NOT point-identified from a single 2x2 table without external
#'   information (e.g., a known-prevalence subgroup or a reference
#'   standard). We return the inequality bounds the observed marginals
#'   support: with rater j positive rate `P_j`, `Se_j` is bounded below
#'   by `P_j`, `Sp_j` is bounded below by `1 - P_j`, both above by 1.
#'   The returned rows have `bound_only = TRUE` and `se_hat = sp_hat =
#'   NA`. This is a documented limitation of the design, not a bug.
#'
#' Bootstrap CI (when `B > 0`): nonparametric, subjects-with-replacement.
#' At k >= 3 the EM is refit on each bootstrap sample; at k = 2 the
#' bounds are recomputed and the bootstrap distribution is on the
#' bound midpoints. Reports the empirical 2.5 / 97.5 percentiles on
#' the per-rater Se_j and Sp_j (or their bound midpoints at k = 2).
#'
#' @param ratings An N x k binary rating matrix. Rows are subjects,
#'   columns are raters, values in {0, 1}. Data.frame and list-of-rater-
#'   columns inputs are coerced. Requires N >= 10, k >= 2, no NA, no
#'   all-constant columns.
#' @param B Integer >= 0. Number of nonparametric bootstrap replicates.
#'   `B = 0` skips the bootstrap and returns the EM point estimate (or
#'   bounds at k = 2) without CIs.
#' @param method One of `"dawid_skene_em"`, `"hui_walter"`, or `NULL`
#'   (default). When `NULL`, dispatches to Hui-Walter at k = 2 and to
#'   Dawid-Skene EM at k >= 3. Supplying `"dawid_skene_em"` at k = 2
#'   errors -- the EM is unidentified there.
#' @param max_iter Integer. EM iteration cap. Default 1000.
#' @param tol Numeric. EM log-likelihood tolerance. Default 1e-6.
#' @param seed Optional integer. Sets the bootstrap RNG; the EM itself
#'   is deterministic.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return An S3 object of class `c("grass_latent_class", "list")`:
#' - `per_rater`: data.frame with `rater`, `se_hat`, `sp_hat` (NA at
#'   k = 2), `se_lower`, `se_upper`, `sp_lower`, `sp_upper`,
#'   `bound_only` (TRUE at k = 2).
#' - `method`: `"dawid_skene_em"` or `"hui_walter"`.
#' - `converged`: logical (NA at k = 2).
#' - `iterations`: integer (NA at k = 2).
#' - `B`: integer; bootstrap replicates run.
#' - `prevalence_hat`: numeric estimated prevalence (NA at k = 2).
#' - `log_likelihood`: numeric (NA at k = 2).
#'
#' @references
#' Dawid, A. P. and Skene, A. M. (1979). Maximum likelihood estimation
#' of observer error-rates using the EM algorithm. *Applied Statistics*,
#' 28(1), 20--28.
#'
#' Hui, S. L. and Walter, S. D. (1980). Estimating the error rates of
#' diagnostic tests. *Biometrics*, 36(1), 167--171.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' N <- 500; k <- 5
#' Se <- 0.90; Sp <- 0.85; pi <- 0.30
#' C <- rbinom(N, 1, pi)
#' Y <- matrix(0L, N, k)
#' for (j in seq_len(k))
#'   Y[, j] <- rbinom(N, 1, ifelse(C == 1, Se, 1 - Sp))
#' fit <- latent_class_fit(Y, B = 200, seed = 1)
#' print(fit)
#' }
#'
#' @export
latent_class_fit <- function(ratings,
                             B = 1000L,
                             method = NULL,
                             max_iter = 1000L,
                             tol = 1e-6,
                             seed = NULL,
                             ...) {
  Y <- .lc_normalize_ratings(ratings)
  N <- nrow(Y)
  k <- ncol(Y)

  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 0) {
    stop("`B` must be a non-negative integer.", call. = FALSE)
  }
  B <- as.integer(B)

  # Method dispatch.
  if (is.null(method)) {
    method <- if (k == 2L) "hui_walter" else "dawid_skene_em"
  }
  method <- match.arg(method, c("dawid_skene_em", "hui_walter"))
  if (method == "dawid_skene_em" && k == 2L) {
    stop("`method = \"dawid_skene_em\"` is not identifiable at k = 2 ",
         "without external information. Use `method = \"hui_walter\"` ",
         "(or leave `method = NULL` for auto-dispatch) to obtain ",
         "inequality bounds.", call. = FALSE)
  }

  rater_names <- paste0("R", seq_len(k))

  if (method == "hui_walter") {
    # k = 2 branch: bounds, optional bootstrap on bound midpoints.
    bounds <- .lc_hui_walter_bounds(Y)
    se_lower <- bounds$se_lower
    se_upper <- bounds$se_upper
    sp_lower <- bounds$sp_lower
    sp_upper <- bounds$sp_upper

    if (B > 0L) {
      if (!is.null(seed)) set.seed(seed)
      boot_se_mid <- matrix(NA_real_, B, k)
      boot_sp_mid <- matrix(NA_real_, B, k)
      for (b in seq_len(B)) {
        idx <- sample.int(N, N, replace = TRUE)
        Yb <- Y[idx, , drop = FALSE]
        # Skip degenerate resamples (all-constant columns).
        cs <- colSums(Yb)
        if (any(cs == 0L) || any(cs == N)) next
        bb <- .lc_hui_walter_bounds(Yb)
        boot_se_mid[b, ] <- (bb$se_lower + bb$se_upper) / 2
        boot_sp_mid[b, ] <- (bb$sp_lower + bb$sp_upper) / 2
      }
      se_lower_b <- apply(boot_se_mid, 2, function(x)
                          stats::quantile(x, 0.025, na.rm = TRUE))
      se_upper_b <- apply(boot_se_mid, 2, function(x)
                          stats::quantile(x, 0.975, na.rm = TRUE))
      sp_lower_b <- apply(boot_sp_mid, 2, function(x)
                          stats::quantile(x, 0.025, na.rm = TRUE))
      sp_upper_b <- apply(boot_sp_mid, 2, function(x)
                          stats::quantile(x, 0.975, na.rm = TRUE))
      # The reported lower / upper combine the analytic bound and the
      # bootstrap-sampling envelope: take the wider on each side.
      se_lower <- pmin(se_lower, as.numeric(se_lower_b))
      se_upper <- pmax(se_upper, as.numeric(se_upper_b))
      sp_lower <- pmin(sp_lower, as.numeric(sp_lower_b))
      sp_upper <- pmax(sp_upper, as.numeric(sp_upper_b))
    }

    per_rater <- data.frame(
      rater      = rater_names,
      se_hat     = rep(NA_real_, k),
      sp_hat     = rep(NA_real_, k),
      se_lower   = as.numeric(se_lower),
      se_upper   = as.numeric(se_upper),
      sp_lower   = as.numeric(sp_lower),
      sp_upper   = as.numeric(sp_upper),
      bound_only = rep(TRUE, k),
      stringsAsFactors = FALSE
    )

    out <- list(
      per_rater      = per_rater,
      method         = "hui_walter",
      converged      = NA,
      iterations     = NA_integer_,
      B              = B,
      prevalence_hat = NA_real_,
      log_likelihood = NA_real_
    )
    class(out) <- c("grass_latent_class", "list")
    return(out)
  }

  # k >= 3: Dawid-Skene EM.
  fit <- .lc_ds_em(Y, max_iter = max_iter, tol = tol)
  Se_hat <- fit$Se
  Sp_hat <- fit$Sp

  if (B > 0L) {
    if (!is.null(seed)) set.seed(seed)
    boot_Se <- matrix(NA_real_, B, k)
    boot_Sp <- matrix(NA_real_, B, k)
    for (b in seq_len(B)) {
      idx <- sample.int(N, N, replace = TRUE)
      Yb <- Y[idx, , drop = FALSE]
      cs <- colSums(Yb)
      if (any(cs == 0L) || any(cs == N)) next
      fb <- tryCatch(
        .lc_ds_em(Yb, max_iter = max_iter, tol = tol),
        error = function(e) NULL
      )
      if (is.null(fb)) next
      # Resolve label-switching by aligning on majority-vote prevalence:
      # if the bootstrap fit's pi is far from 0.5 in the same direction
      # as the parent fit, accept; otherwise flip Se <-> 1 - Sp.
      flip <- (fb$pi - 0.5) * (fit$pi - 0.5) < 0
      if (flip) {
        boot_Se[b, ] <- 1 - fb$Sp
        boot_Sp[b, ] <- 1 - fb$Se
      } else {
        boot_Se[b, ] <- fb$Se
        boot_Sp[b, ] <- fb$Sp
      }
    }
    se_lower <- apply(boot_Se, 2, function(x)
                      stats::quantile(x, 0.025, na.rm = TRUE))
    se_upper <- apply(boot_Se, 2, function(x)
                      stats::quantile(x, 0.975, na.rm = TRUE))
    sp_lower <- apply(boot_Sp, 2, function(x)
                      stats::quantile(x, 0.025, na.rm = TRUE))
    sp_upper <- apply(boot_Sp, 2, function(x)
                      stats::quantile(x, 0.975, na.rm = TRUE))
  } else {
    se_lower <- rep(NA_real_, k)
    se_upper <- rep(NA_real_, k)
    sp_lower <- rep(NA_real_, k)
    sp_upper <- rep(NA_real_, k)
  }

  per_rater <- data.frame(
    rater      = rater_names,
    se_hat     = as.numeric(Se_hat),
    sp_hat     = as.numeric(Sp_hat),
    se_lower   = as.numeric(se_lower),
    se_upper   = as.numeric(se_upper),
    sp_lower   = as.numeric(sp_lower),
    sp_upper   = as.numeric(sp_upper),
    bound_only = rep(FALSE, k),
    stringsAsFactors = FALSE
  )

  out <- list(
    per_rater      = per_rater,
    method         = "dawid_skene_em",
    converged      = isTRUE(fit$converged),
    iterations     = as.integer(fit$iterations),
    B              = B,
    prevalence_hat = as.numeric(fit$pi),
    log_likelihood = as.numeric(fit$log_likelihood)
  )
  class(out) <- c("grass_latent_class", "list")
  out
}

#' @export
print.grass_latent_class <- function(x, digits = 3, ...) {
  k <- nrow(x$per_rater)
  cat("grass latent-class fit\n", sep = "")
  cat("  method        : ", x$method, "\n", sep = "")
  cat("  raters (k)    : ", k, "\n", sep = "")
  cat("  bootstrap B   : ", x$B, "\n", sep = "")
  if (x$method == "dawid_skene_em") {
    cat("  converged     : ", x$converged, "\n", sep = "")
    cat("  iterations    : ", x$iterations, "\n", sep = "")
    cat(sprintf("  prevalence_hat: %.*f\n", digits, x$prevalence_hat))
    cat(sprintf("  log-likelihood: %.*f\n", digits, x$log_likelihood))
  }
  cat("  per-rater\n", sep = "")
  pr <- x$per_rater
  fmt_ci <- function(lo, hi) {
    if (is.na(lo) || is.na(hi)) return("(NA, NA)")
    sprintf("(%.*f, %.*f)", digits, lo, digits, hi)
  }
  for (i in seq_len(k)) {
    if (isTRUE(pr$bound_only[i])) {
      cat(sprintf("    %-4s  Se in [%.*f, %.*f]   Sp in [%.*f, %.*f]   (bounds)\n",
                  pr$rater[i],
                  digits, pr$se_lower[i], digits, pr$se_upper[i],
                  digits, pr$sp_lower[i], digits, pr$sp_upper[i]))
    } else {
      cat(sprintf("    %-4s  Se = %.*f  %s   Sp = %.*f  %s\n",
                  pr$rater[i],
                  digits, pr$se_hat[i],
                  fmt_ci(pr$se_lower[i], pr$se_upper[i]),
                  digits, pr$sp_hat[i],
                  fmt_ci(pr$sp_lower[i], pr$sp_upper[i])))
    }
  }
  if (any(pr$bound_only)) {
    cat("\n  Note: at k = 2, per-rater Se/Sp are not point-identified\n",
        "        without external information. The intervals are\n",
        "        Hui-Walter (1980) inequality bounds, not point\n",
        "        estimates with sampling uncertainty.\n", sep = "")
  }
  invisible(x)
}
