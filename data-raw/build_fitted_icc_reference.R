# build_fitted_icc_reference.R -- Precompute GLMM-gap-corrected ICC reference
# curves for grass::position_on_surface(metric = "icc", ratings = ...).
#
# Background. The bundled `icc_reference_curves` (see build_icc_reference_curves.R)
# holds ORACLE E[ICC](q | F) -- the closed-form
# sigma^2_subject / (sigma^2_subject + pi^2/3) with sigma^2_subject known from
# F's parameters. Practitioners, however, compute obs_ICC via glmer -- a FITTED
# ICC that systematically inflates sigma^2_subject under the clustered-latent-
# class DGP because glmer absorbs rater classification noise into the variance
# components. This is the documented GLMM gap (framework_notes.md §0.4.iii +
# cont_sim probe #5: +0.013 at q=0.55 growing to +0.647 at q=0.99).
#
# Inverting a practitioner's obs_ICC against the ORACLE reference therefore
# clamps routinely. Inverting against the FITTED reference (= oracle + bias_emp
# correction from paper2's q_recovery_fitted_icc.rds sim) gives a recovery
# that matches the scale the practitioner actually computes.
#
# Construction:
#   For each (F_key, k, N) cell in the sim grid
#   (52 F_keys x {3,5,8,15} k x {50, 200} N = 416 cells):
#     bias_emp(q) is given empirically at 10 q values in [0.55, 0.99]
#     interpolate linearly across q to produce a 501-point bias-correction curve
#   fitted_ref(q | F_key, k, N) = oracle_ref(q | F_key) + bias_corr(q | F_key, k, N)
#
# Size: 52 F_keys x 4 k x 2 N x 501 q-grid = 208,416 numbers ~ 1.7 MB uncompressed,
# ~400 KB xz-compressed. Stored as 4D array [F_key, k, N, q] with dimnames.
#
# Run from PABAK_Investigation project root:
#   Rscript grass/data-raw/build_fitted_icc_reference.R

suppressMessages({
  source("paper2/code/04_reference_closed_form.R")
  source("paper2/code/06_grid.R")
})

FITTED_RDS <- "paper2/simulation_output/multirater_sim_v3/q_recovery_fitted_icc.rds"
SYSDATA    <- "grass/R/sysdata.rda"
Q_GRID     <- seq(0.5, 1.0, length.out = 501L)

# ---- Load fitted-ICC bias table ------------------------------------------
cat(sprintf("Loading fitted-ICC sim recovery table from %s ...\n", FITTED_RDS))
fit_tbl <- readRDS(FITTED_RDS)
cat(sprintf("  %d rows (%d F_keys, k in {%s}, N in {%s}, q_true in {%s})\n",
            nrow(fit_tbl),
            length(unique(fit_tbl$F_key)),
            paste(sort(unique(fit_tbl$k)), collapse = ","),
            paste(sort(unique(fit_tbl$N)), collapse = ","),
            paste(sort(unique(fit_tbl$q_true)), collapse = ",")))

F_KEYS <- sort(unique(fit_tbl$F_key))
K_VALUES <- sort(unique(fit_tbl$k))
N_VALUES <- sort(unique(fit_tbl$N))
Q_SIM    <- sort(unique(fit_tbl$q_true))

# ---- Load existing sysdata (for oracle reference curves) -----------------
e <- new.env(parent = emptyenv())
load(SYSDATA, envir = e)
stopifnot("icc_reference_curves" %in% ls(e))
oracle <- e$icc_reference_curves
stopifnot(all(F_KEYS %in% rownames(oracle$curves)))
cat(sprintf("Oracle reference loaded: %d F_keys x %d q-points.\n",
            nrow(oracle$curves), ncol(oracle$curves)))

# ---- Build 4D fitted reference array --------------------------------------
cat("Building fitted reference curves by cell ...\n")

fitted_arr <- array(
  NA_real_,
  dim = c(length(F_KEYS), length(K_VALUES), length(N_VALUES), length(Q_GRID)),
  dimnames = list(
    F_key = F_KEYS,
    k     = as.character(K_VALUES),
    N     = as.character(N_VALUES),
    q     = NULL
  )
)

# Index for (F_key, F_family, M1, mu, tau2) metadata per F_key
F_family_by_key <- character(length(F_KEYS))
M1_by_key       <- numeric(length(F_KEYS))
mu_by_key       <- numeric(length(F_KEYS))
tau2_by_key     <- numeric(length(F_KEYS))

n_cells <- length(F_KEYS) * length(K_VALUES) * length(N_VALUES)
cell_i  <- 0L
n_skipped_cells <- 0L

for (fi in seq_along(F_KEYS)) {
  key <- F_KEYS[fi]
  params <- F_params_registry[[key]]

  if (startsWith(key, "LN_")) {
    F_family_by_key[fi] <- "logit_normal"
    M1_by_key[fi]       <- e_p_logit_normal(params$mu, params$tau2)
    mu_by_key[fi]       <- params$mu
    tau2_by_key[fi]     <- params$tau2
  } else if (startsWith(key, "DM_")) {
    F_family_by_key[fi] <- "discrete_mixture"
    M1_by_key[fi]       <- e_p_discrete_mixture(params$p, params$w)
    mu_by_key[fi]       <- NA_real_
    tau2_by_key[fi]     <- NA_real_
  }

  oracle_row_idx <- which(rownames(oracle$curves) == key)
  oracle_curve   <- as.numeric(oracle$curves[oracle_row_idx, ])

  for (ki in seq_along(K_VALUES)) {
    for (ni in seq_along(N_VALUES)) {
      cell_i <- cell_i + 1L
      kv <- K_VALUES[ki]
      nv <- N_VALUES[ni]

      cell <- fit_tbl[fit_tbl$F_key == key &
                      fit_tbl$k == kv &
                      fit_tbl$N == nv, ]
      if (nrow(cell) == 0L) {
        n_skipped_cells <- n_skipped_cells + 1L
        next
      }

      # Sort by q_true, interpolate bias_emp across q_grid (linear, flat extrap)
      cell <- cell[order(cell$q_true), ]
      bias_finite <- is.finite(cell$bias_emp)
      if (sum(bias_finite) < 2L) {
        n_skipped_cells <- n_skipped_cells + 1L
        next
      }
      bias_curve <- stats::approx(
        x = cell$q_true[bias_finite],
        y = cell$bias_emp[bias_finite],
        xout = Q_GRID, rule = 2
      )$y

      # Fitted reference = oracle + bias (bias is positive; glmer inflates)
      fitted_arr[fi, ki, ni, ] <- oracle_curve + bias_curve
    }
  }
}

cat(sprintf("  done. %d of %d cells populated (%d skipped).\n",
            n_cells - n_skipped_cells, n_cells, n_skipped_cells))

# ---- Index frame for nearest-neighbor lookup -----------------------------
index <- data.frame(
  F_key    = F_KEYS,
  F_family = F_family_by_key,
  M1       = M1_by_key,
  mu       = mu_by_key,
  tau2     = tau2_by_key,
  stringsAsFactors = FALSE
)

fitted_icc_reference_curves <- list(
  curves   = fitted_arr,         # [F_key, k, N, q]
  q_grid   = Q_GRID,
  k_grid   = K_VALUES,
  N_grid   = N_VALUES,
  index    = index,              # F_key, F_family, M1, mu, tau2
  source   = paste0(
    "paper2/simulation_output/multirater_sim_v3/q_recovery_fitted_icc.rds + ",
    "paper2/code/04_reference_closed_form.R :: cf_logit_mixed_icc_* (oracle base)"
  ),
  construction = "fitted_ref(q | F_key, k, N) = oracle_ref(q | F_key) + bias_emp(q | F_key, k, N) interpolated linearly across 10 simulated q-values.",
  N_limit_note = paste0(
    "Fitted reference is available for N in {", paste(N_VALUES, collapse = ", "),
    "}. For N > ", max(N_VALUES),
    " use oracle reference (asymptotic bias is ~0)."
  ),
  built_on = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

cat(sprintf("\nAssembled fitted_icc_reference_curves:\n"))
cat(sprintf("  curves   : %d F_keys x %d k x %d N x %d q-points\n",
            dim(fitted_arr)[1], dim(fitted_arr)[2], dim(fitted_arr)[3],
            dim(fitted_arr)[4]))
cat(sprintf("  finite   : %d / %d cells populated\n",
            sum(!is.na(fitted_arr[, , , 251L])),
            prod(dim(fitted_arr)[1:3])))

# ---- Merge into sysdata.rda (preserve existing objects) ------------------
assign("fitted_icc_reference_curves", fitted_icc_reference_curves, envir = e)
keep <- ls(e)
cat(sprintf("Saving objects to %s: %s\n",
            SYSDATA, paste(keep, collapse = ", ")))
save(list = keep, file = SYSDATA, envir = e, compress = "xz")

cat(sprintf("\nsysdata.rda size after save: %.1f KB\n",
            file.info(SYSDATA)$size / 1024))
cat("DONE.\n")
