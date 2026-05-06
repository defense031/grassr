# build_icc_reference_curves.R -- Precompute ICC reference curves for
# grass::position_on_surface(metric = "icc").
#
# Unlike the four agreement-family metrics, whose closed-form E[metric]
# reduces to a function of (q, pi_+) under the symmetric DGP, ICC's closed
# form depends on sigma^2_subject, which is a function of the FULL F-shape
# (mu, tau2 for logit-normal; p_vec, w_vec for discrete-mixture), not just
# the marginal pi_+. The runtime path in position_on_surface() therefore
# cannot reconstruct E[ICC](q | F) from pi_hat alone; a bundled lookup is
# required if practitioners are to use the function without manually
# supplying surface_data$reference_curve.
#
# This script evaluates the ICC closed form at the 501-point q-grid for
# every F_key in paper2/code/06_grid.R's F_params_registry (48 logit-normal
# keys + 4 discrete-mixture keys = 52 total) and bundles the result into
# grass/R/sysdata.rda alongside the existing `reference_binary` and
# `empirical_q_hat_surface` objects.
#
# The closed form is k-independent in the large-N limit, so we do NOT key
# by k. At runtime, nearest-neighbor F_key lookup is by M1 (= E_F[p])
# matched against the user's observed pi_hat.
#
# Size: 52 F_keys x 501 q-points x 8 bytes ~ 210 KB uncompressed.
# Compressed under xz probably ~80 KB. Negligible sysdata delta.
#
# Run from PABAK_Investigation project root:
#   Rscript grass/data-raw/build_icc_reference_curves.R

suppressMessages({
  source("paper2/code/04_reference_closed_form.R")
  source("paper2/code/06_grid.R")
})

SYSDATA <- "grass/R/sysdata.rda"
Q_GRID  <- seq(0.5, 1.0, length.out = 501L)

# ---- Pull the F_params registry ------------------------------------------
f_keys <- names(F_params_registry)
stopifnot(length(f_keys) > 0)
cat(sprintf("F_params registry has %d keys.\n", length(f_keys)))

# ---- Evaluate E[ICC](q | F) at Q_GRID for every F_key --------------------
icc_matrix <- matrix(NA_real_,
                     nrow = length(f_keys), ncol = length(Q_GRID),
                     dimnames = list(F_key = f_keys, q = NULL))

F_family_by_key <- character(length(f_keys))
M1_by_key       <- numeric(length(f_keys))

cat("Evaluating cf_logit_mixed_icc_* at 501-point q-grid for every F_key ...\n")
for (i in seq_along(f_keys)) {
  key <- f_keys[i]
  params <- F_params_registry[[key]]
  if (startsWith(key, "LN_")) {
    F_family_by_key[i] <- "logit_normal"
    M1_by_key[i] <- e_p_logit_normal(params$mu, params$tau2)
    icc_matrix[i, ] <- cf_logit_mixed_icc_logit_normal(
      q = Q_GRID, mu = params$mu, tau2 = params$tau2
    )
  } else if (startsWith(key, "DM_")) {
    F_family_by_key[i] <- "discrete_mixture"
    M1_by_key[i] <- e_p_discrete_mixture(params$p, params$w)
    icc_matrix[i, ] <- cf_logit_mixed_icc_discrete_mixture(
      q = Q_GRID, p_vec = params$p, w_vec = params$w
    )
  } else {
    stop("Unrecognized F_key prefix: ", key)
  }
}

cat(sprintf("  done (%d F_keys).\n", length(f_keys)))
cat(sprintf("  ICC range across grid: [%.4f, %.4f]\n",
            min(icc_matrix, na.rm = TRUE),
            max(icc_matrix, na.rm = TRUE)))
cat(sprintf("  M1 range across F_keys: [%.4f, %.4f]\n",
            min(M1_by_key), max(M1_by_key)))

# ---- Index data.frame for nearest-neighbor lookup ------------------------
index <- data.frame(
  F_key    = f_keys,
  F_family = F_family_by_key,
  M1       = M1_by_key,
  stringsAsFactors = FALSE
)

# ---- Assemble the bundled object -----------------------------------------
icc_reference_curves <- list(
  curves   = icc_matrix,                 # [n_F_keys x 501]
  q_grid   = Q_GRID,
  index    = index,                      # columns: F_key, F_family, M1
  source   = "paper2/code/04_reference_closed_form.R :: cf_logit_mixed_icc_*",
  built_on = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

cat(sprintf("\nAssembled icc_reference_curves:\n"))
cat(sprintf("  curves   : %d F_keys x %d q-points\n",
            nrow(icc_matrix), ncol(icc_matrix)))
cat(sprintf("  families : logit_normal=%d, discrete_mixture=%d\n",
            sum(F_family_by_key == "logit_normal"),
            sum(F_family_by_key == "discrete_mixture")))

# ---- Merge into existing sysdata.rda (preserve other objects) ------------
e <- new.env(parent = emptyenv())
if (file.exists(SYSDATA)) {
  load(SYSDATA, envir = e)
  cat(sprintf("\nExisting sysdata objects: %s\n",
              paste(ls(e), collapse = ", ")))
}

assign("icc_reference_curves", icc_reference_curves, envir = e)

keep <- ls(e)
cat(sprintf("Saving objects to %s: %s\n",
            SYSDATA, paste(keep, collapse = ", ")))
save(list = keep, file = SYSDATA, envir = e, compress = "xz")

cat(sprintf("\nsysdata.rda size after save: %.1f KB\n",
            file.info(SYSDATA)$size / 1024))
cat("DONE.\n")
