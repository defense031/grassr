# =====================================================================
# pairwise_agreement.R -- Pairwise reliability for divergent panels
# =====================================================================
#
# When grass_report() flags a panel as divergent, the framework's
# prescription (paper §3.3) is to abandon the panel-aggregate summary
# and instead report:
#
#   (1) a k x k pairwise PABAK matrix, with each entry placed on the
#       k = 2 reference surface at the pair's observed marginal; and
#   (2) per-rater (Se_tilde_j, Sp_tilde_j) against the panel-majority
#       of the OTHER k - 1 raters as a proxy reference.
#
# (1) preserves identifiability (PABAK_ij is well-defined regardless
# of which rater is "more accurate"); (2) uses the larger panel's
# information about each rater rather than discarding it to strictly
# pairwise comparisons. Together they answer the question the
# divergent flag actually raises: which raters are out of step with
# the panel and in which direction.

#' Pairwise reliability for a divergent panel
#'
#' Computes pairwise PABAK between every pair of raters in an
#' \eqn{N \times k} binary rating matrix and places each entry on the
#' \eqn{k = 2} reference surface at the pair's observed marginal. Also
#' returns per-rater pooled-reference sensitivity and specificity --- each
#' rater's call rate against the panel-majority of the other
#' \eqn{k - 1} raters --- which uses the larger panel's information
#' about rater behavior rather than discarding it to a strictly
#' pairwise comparison.
#'
#' This function is the recommended primary deliverable when
#' \code{grass_report()} flags the panel as \emph{divergent}; the
#' panel-aggregate coefficients no longer summarize the panel
#' adequately, but the pairwise matrix exposes the panel's structure
#' directly (uniform inconsistency, sub-group clustering, or single-
#' rater outliers).
#'
#' @param ratings An \eqn{N \times k} binary rating matrix
#'   (rows = subjects, columns = raters), or a data.frame whose
#'   columns are 0/1 / logical / two-level factor. Same input
#'   conventions as \code{\link{grass_report}}.
#' @param axis Character; \code{"inter"} (default) or \code{"intra"}.
#'   The intra-axis path treats the columns of \code{ratings} as
#'   per-rater viewing pairs; see paper §3.3 intra-rater section.
#'
#' @return An object of class \code{c("grass_pairwise", "list")} with
#'   fields:
#' \itemize{
#'   \item \code{pabak_matrix} --- \eqn{k \times k} symmetric numeric
#'     matrix; \eqn{[i, j]} is \eqn{\PABAK_{ij}}, the diagonal is 1.
#'   \item \code{percentile_matrix} --- \eqn{k \times k} symmetric
#'     numeric (0--100); \eqn{[i, j]} is the surface percentile of
#'     \eqn{\PABAK_{ij}} on the \eqn{k = 2} reference at the pair's
#'     observed marginal. Diagonal is \code{NA}.
#'   \item \code{marginal_matrix} --- \eqn{k \times k} symmetric
#'     numeric; \eqn{[i, j]} is \eqn{\hat\pi_{+, ij}}. Diagonal is
#'     \code{NA}.
#'   \item \code{band_matrix} --- \eqn{k \times k} character matrix;
#'     four-band label for each pairwise percentile. Diagonal is
#'     \code{NA}.
#'   \item \code{qualifier_matrix} --- \eqn{k \times k} character;
#'     decisive / moderate / weak per pairwise percentile.
#'   \item \code{pooled_per_rater} --- data frame with \eqn{k} rows,
#'     one per rater, columns \code{rater}, \code{se_tilde},
#'     \code{sp_tilde}, \code{n_pool_pos}, \code{n_pool_neg},
#'     \code{n_pool_excluded} (subjects with tied panel majority).
#'   \item \code{sample} --- list with \code{k}, \code{N},
#'     \code{pi_hat}, \code{tau2_hat}, \code{axis}.
#'   \item \code{notes} --- character vector with caveats, if any
#'     (e.g., undefined pooled-reference at \eqn{k = 2}).
#'   \item \code{call} --- the matched call.
#' }
#' @export
#' @examples
#' set.seed(6)
#' Se <- c(0.95, 0.75, 0.95, 0.75, 0.95)
#' Sp <- c(0.75, 0.95, 0.75, 0.95, 0.75)
#' truth <- rbinom(200, 1, 0.5)
#' Y <- sapply(seq_along(Se),
#'             function(j) ifelse(truth == 1,
#'                                rbinom(200, 1, Se[j]),
#'                                rbinom(200, 1, 1 - Sp[j])))
#' pw <- pairwise_agreement(Y)
#' pw
pairwise_agreement <- function(ratings, axis = c("inter", "intra")) {
  call <- match.call()
  axis <- match.arg(axis)
  Y <- normalize_ratings(ratings)
  k <- ncol(Y); N <- nrow(Y)
  if (k < 2L) {
    stop("`ratings` must have at least 2 rater columns (k >= 2).",
         call. = FALSE)
  }

  pi_hat   <- mean(Y)
  tau2_hat <- compute_tau2_hat(Y)

  notes <- character(0L)

  # ---- (1) Pairwise PABAK matrix + surface placement ------------------
  pabak_mat   <- matrix(NA_real_, nrow = k, ncol = k)
  pct_mat     <- matrix(NA_real_, nrow = k, ncol = k)
  marg_mat    <- matrix(NA_real_, nrow = k, ncol = k)
  band_mat    <- matrix(NA_character_, nrow = k, ncol = k)
  qual_mat    <- matrix(NA_character_, nrow = k, ncol = k)
  diag(pabak_mat) <- 1
  rater_names <- colnames(Y) %||% paste0("R", seq_len(k))
  dimnames(pabak_mat) <- dimnames(pct_mat) <- dimnames(marg_mat) <-
    dimnames(band_mat) <- dimnames(qual_mat) <-
      list(rater_names, rater_names)

  for (i in seq_len(k - 1L)) {
    for (j in seq.int(i + 1L, k)) {
      Yi <- Y[, i]; Yj <- Y[, j]
      P_o_ij     <- mean(Yi == Yj)
      pabak_ij   <- 2 * P_o_ij - 1
      pi_hat_ij  <- mean(c(Yi, Yj))

      pos <- tryCatch(
        position_on_surface(obs_value = pabak_ij,
                            pi_hat    = pi_hat_ij,
                            k         = 2L,
                            N         = N,
                            metric    = "pabak"),
        error = function(e) NULL)

      pabak_mat[i, j] <- pabak_mat[j, i] <- pabak_ij
      marg_mat[i, j]  <- marg_mat[j, i]  <- pi_hat_ij

      if (!is.null(pos)) {
        pct_pp <- as.numeric(pos$percentile) * 100
        pct_mat[i, j]  <- pct_mat[j, i]  <- pct_pp
        band_mat[i, j] <- band_mat[j, i] <-
          pos$modal_band_label %||% NA_character_
        qual_mat[i, j] <- qual_mat[j, i] <-
          pos$confidence %||% NA_character_
      }
    }
  }

  # ---- (2) Per-rater pooled-reference Se/Sp ---------------------------
  pooled <- data.frame(
    rater          = rater_names,
    se_tilde       = NA_real_,
    sp_tilde       = NA_real_,
    n_pool_pos     = NA_integer_,
    n_pool_neg     = NA_integer_,
    n_pool_excluded = NA_integer_,
    stringsAsFactors = FALSE
  )
  if (k >= 3L) {
    for (j in seq_len(k)) {
      others <- Y[, -j, drop = FALSE]
      others_mean <- rowMeans(others)
      pooled_call <- ifelse(others_mean > 0.5, 1L,
                     ifelse(others_mean < 0.5, 0L, NA_integer_))
      excluded <- is.na(pooled_call)
      n_excl <- sum(excluded)
      keep <- !excluded
      yj <- Y[keep, j]
      pc <- pooled_call[keep]
      n_pos <- sum(pc == 1L)
      n_neg <- sum(pc == 0L)
      pooled$se_tilde[j]       <- if (n_pos > 0) sum(yj == 1L & pc == 1L) / n_pos else NA_real_
      pooled$sp_tilde[j]       <- if (n_neg > 0) sum(yj == 0L & pc == 0L) / n_neg else NA_real_
      pooled$n_pool_pos[j]     <- n_pos
      pooled$n_pool_neg[j]     <- n_neg
      pooled$n_pool_excluded[j] <- n_excl
    }
  } else {
    notes <- c(notes,
      "Pooled-reference per-rater (Se_tilde, Sp_tilde) is undefined at k < 3 (no `rest of the panel` to pool against); the per-rater table is omitted.")
  }

  out <- list(
    pabak_matrix      = pabak_mat,
    percentile_matrix = pct_mat,
    marginal_matrix   = marg_mat,
    band_matrix       = band_mat,
    qualifier_matrix  = qual_mat,
    pooled_per_rater  = pooled,
    sample = list(
      k        = as.integer(k),
      N        = as.integer(N),
      pi_hat   = as.numeric(pi_hat),
      tau2_hat = as.numeric(tau2_hat),
      axis     = axis
    ),
    notes = notes,
    call  = call
  )
  class(out) <- c("grass_pairwise", "list")
  out
}

#' @export
print.grass_pairwise <- function(x, digits = 2, ...) {
  cat("GRASS Pairwise Reliability\n\n")
  cat(sprintf("  sample      = %d raters, N = %d, pi_hat = %.*f, tau2_hat = %.*f, axis = %s\n",
              x$sample$k, x$sample$N,
              digits, x$sample$pi_hat,
              max(digits + 1L, 3L), x$sample$tau2_hat,
              x$sample$axis))
  cat("\n  Pairwise PABAK (lower triangle = PABAK_ij; upper triangle = surface percentile):\n\n")
  k <- x$sample$k
  M <- matrix("", nrow = k, ncol = k)
  rn <- rownames(x$pabak_matrix)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i == j) {
        M[i, j] <- "    --"
      } else if (i > j) {
        M[i, j] <- formatC(x$pabak_matrix[i, j], format = "f",
                           digits = digits, width = 6)
      } else {
        pct <- x$percentile_matrix[i, j]
        M[i, j] <- if (is.finite(pct))
          sprintf("%4.0f%%", pct) else "    NA"
      }
    }
  }
  colnames(M) <- rn; rownames(M) <- rn
  print(noquote(M))

  if (k >= 3L) {
    cat("\n  Per-rater behavior against pooled panel-majority:\n\n")
    pp <- x$pooled_per_rater
    pp_show <- data.frame(
      rater    = pp$rater,
      Se_tilde = formatC(pp$se_tilde, format = "f", digits = digits),
      Sp_tilde = formatC(pp$sp_tilde, format = "f", digits = digits),
      n_pos    = pp$n_pool_pos,
      n_neg    = pp$n_pool_neg,
      n_excl   = pp$n_pool_excluded
    )
    print(pp_show, row.names = FALSE)
    cat("  (Se_tilde, Sp_tilde are calls vs panel-majority of OTHER raters,\n",
        "   not against external truth. n_excl: subjects with tied majority,\n",
        "   excluded from the per-rater pool.)\n", sep = "")
  }

  if (length(x$notes) > 0L) {
    cat("\n  Notes:\n")
    for (note in x$notes) cat("    -", note, "\n")
  }

  invisible(x)
}
