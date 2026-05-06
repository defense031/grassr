#!/usr/bin/env Rscript
# build_fitted_icc_reference_v0.4.R -- Bundle the unified fitted-ICC
# reference into grass::sysdata as `fitted_icc_reference_curves`.
#
# Background. v0.3 built the fitted reference as
#   fitted_ref(q | F, k, N) = oracle_ref(q | F) + bias_emp(q | F, k, N)
# from `q_recovery_fitted_icc.rds`. That coverage stopped at
# k in {3, 5, 8, 15} and N in {50, 200}; queries at small k=2 / k=25 or
# N in {30, 75, 100, 300, 500, 1000} routed silently to the oracle path
# and clamped at a 98.2% rate against common applied designs.
#
# v0.4 replaces the construction. The unified Julia sim
# (paper2/simulation_output/fitted_icc_unified.rds) directly simulates
# glmer-fitted ICC at every (F_key, k, N, q) cell on the extended grid:
#   F_keys : 52   (LN logit-normal grid + 4 DM presets)
#   k      : {3, 5, 8, 15, 25}
#   N      : {30, 50, 75, 100, 200, 300, 500, 1000}
#   q      : {0.35, 0.40, ..., 0.95, 0.99}   (14 points)
# Each cell averages 500-1000 fitted-ICC reps from MixedModels.jl GLMM
# fits. The unified `mean_icc` per cell IS the fitted reference value
# at that q -- no oracle base / bias_emp split.
#
# Construction:
#   For each (F_key, k, N) group, take the (q, mean_icc) pairs and
#   linearly interpolate mean_icc onto the 501-point q-grid
#   seq(0.5, 1.0, length.out = 501). cummax-enforce monotonicity along q
#   so the downstream inverter (invert_metric_to_q in
#   grass/R/position_on_surface.R) gets a clean monotone curve. Above
#   the highest sim q (0.99), extrapolate flat via approx(rule = 2).
#
# Size: 52 F_keys x 5 k x 8 N x 501 q = 1,042,080 numbers ~ 8.3 MB
# uncompressed, ~3-4 MB xz-compressed in sysdata.rda.
#
# Run from PABAK_Investigation project root:
#   Rscript grass/data-raw/build_fitted_icc_reference_v0.4.R

UNIFIED_RDS <- "paper2/simulation_output/fitted_icc_unified.rds"
SYSDATA     <- "grass/R/sysdata.rda"
Q_GRID      <- seq(0.5, 1.0, length.out = 501L)

# ---- Load unified sim ----------------------------------------------------
cat(sprintf("Loading unified fitted-ICC sim from %s ...\n", UNIFIED_RDS))
u <- readRDS(UNIFIED_RDS)
tbl <- u$index
cat(sprintf("  %s rows; F_keys=%d; k in {%s}; N in {%s}; q in {%s}\n",
            format(nrow(tbl), big.mark = ","),
            length(unique(tbl$F_key)),
            paste(sort(unique(tbl$k)), collapse = ","),
            paste(sort(unique(tbl$N)), collapse = ","),
            paste(sort(unique(tbl$q)), collapse = ",")))

F_KEYS   <- sort(unique(tbl$F_key))
K_VALUES <- as.integer(sort(unique(tbl$k)))
N_VALUES <- as.integer(sort(unique(tbl$N)))

# ---- Load existing sysdata (for F_key index metadata + preserve other slots) ----
e <- new.env(parent = emptyenv())
load(SYSDATA, envir = e)
stopifnot("fitted_icc_reference_curves" %in% ls(e))
old_index <- e$fitted_icc_reference_curves$index
stopifnot(nrow(old_index) == 52L,
          all(c("F_key", "F_family", "M1", "mu", "tau2") %in% colnames(old_index)))
stopifnot(setequal(old_index$F_key, F_KEYS))

# Order index rows to match F_KEYS so dimnames(curves)$F_key is consistent
index <- old_index[match(F_KEYS, old_index$F_key), , drop = FALSE]
rownames(index) <- NULL
cat(sprintf("Reused F_key index from existing sysdata (%d F_keys; %d LN, %d DM).\n",
            nrow(index), sum(index$F_family == "logit_normal"),
            sum(index$F_family == "discrete_mixture")))

# ---- Build 4D array via per-group monotone interpolation ----------------
cat("\nInterpolating mean_icc onto 501-point q-grid per (F_key, k, N) cell ...\n")

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

# Index the unified table for fast subset
tbl_split <- split(
  tbl,
  f = paste(tbl$F_key, tbl$k, tbl$N, sep = "|")
)

n_cells <- length(F_KEYS) * length(K_VALUES) * length(N_VALUES)
n_done    <- 0L
n_skipped <- 0L
n_nonmono <- 0L

for (fi in seq_along(F_KEYS)) {
  key <- F_KEYS[fi]
  for (ki in seq_along(K_VALUES)) {
    kv <- K_VALUES[ki]
    for (ni in seq_along(N_VALUES)) {
      nv <- N_VALUES[ni]
      cell_id <- paste(key, kv, nv, sep = "|")
      cell <- tbl_split[[cell_id]]
      if (is.null(cell) || nrow(cell) == 0L) {
        n_skipped <- n_skipped + 1L
        next
      }
      cell <- cell[order(cell$q), ]
      y <- cell$mean_icc
      x <- cell$q
      finite <- is.finite(y)
      if (sum(finite) < 2L) {
        n_skipped <- n_skipped + 1L
        next
      }
      x <- x[finite]; y <- y[finite]
      # Enforce monotone non-decreasing along q so the downstream inverter
      # in invert_metric_to_q gets a clean curve. Tiny epsilon ladder
      # breaks ties without producing a duplicate-x error in approx().
      if (any(diff(y) < 0)) {
        n_nonmono <- n_nonmono + 1L
        y <- cummax(y)
      }
      ladder <- seq(0, 1e-12, length.out = length(y))
      y <- y + ladder
      ref <- stats::approx(x = x, y = y, xout = Q_GRID, rule = 2,
                           ties = "ordered")$y
      fitted_arr[fi, ki, ni, ] <- ref
      n_done <- n_done + 1L
    }
  }
}

cat(sprintf("  %d / %d cells populated; %d skipped (no data); %d cells monotone-corrected.\n",
            n_done, n_cells, n_skipped, n_nonmono))

# ---- Sanity: §4.1 design point should not clamp -------------------------
cat("\n=== §4.1 sanity check ===\n")
target_F <- "LN_mu=-1.386_tau2=1.0000"
target_k <- 5L
target_N <- 1000L
target_obs <- 0.332
ref_curve <- fitted_arr[target_F, as.character(target_k), as.character(target_N), ]
rc_finite <- is.finite(ref_curve)
if (sum(rc_finite) >= 2L) {
  rng <- range(ref_curve[rc_finite])
  cat(sprintf("  F_key=%s, k=%d, N=%d\n", target_F, target_k, target_N))
  cat(sprintf("  ICC range on q-grid [0.5, 1.0]: [%.4f, %.4f]\n", rng[1], rng[2]))
  if (target_obs >= rng[1] && target_obs <= rng[2]) {
    qg <- Q_GRID[rc_finite]; rc <- ref_curve[rc_finite]
    ord <- order(rc)
    q_hat <- stats::approx(x = rc[ord], y = qg[ord],
                           xout = target_obs, rule = 2, ties = "ordered")$y
    cat(sprintf("  obs ICC = %.3f -> q_hat = %.3f (no clamp)\n",
                target_obs, q_hat))
  } else {
    cat(sprintf("  obs ICC = %.3f outside envelope [%.4f, %.4f] (CLAMPS)\n",
                target_obs, rng[1], rng[2]))
  }
} else {
  cat("  cell empty (sanity check failed)\n")
}

# ---- Assemble bundle -----------------------------------------------------
fitted_icc_reference_curves <- list(
  curves    = fitted_arr,
  q_grid    = Q_GRID,
  k_grid    = K_VALUES,
  N_grid    = N_VALUES,
  index     = index,
  source    = paste0(
    "paper2/simulation_output/fitted_icc_unified.rds (built ",
    u$metadata$built_at, "); ",
    "Julia ", u$metadata$julia_version, " + MixedModels.jl ",
    u$metadata$mixedmodels
  ),
  construction = paste0(
    "Per-cell glmer-fitted ICC mean from a unified Julia/MixedModels.jl sim ",
    "across 52 F_keys x 5 k {3,5,8,15,25} x 8 N {30,50,75,100,200,300,500,1000} ",
    "x 14 q [0.35, 0.99]. Per (F_key, k, N) the (q, mean_icc) points are ",
    "monotone-enforced (cummax) and linearly interpolated onto the 501-point ",
    "q-grid on [0.5, 1.0] with flat extrapolation above q=0.99 (approx rule=2)."
  ),
  N_limit_note = paste0(
    "Fitted reference is available for k in {", paste(K_VALUES, collapse = ", "),
    "} and N in {", paste(N_VALUES, collapse = ", "),
    "}. Out-of-range queries snap to nearest sim cell with a flag in `notes`."
  ),
  built_on = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

cat(sprintf("\nAssembled fitted_icc_reference_curves:\n"))
cat(sprintf("  curves   : %d F_keys x %d k x %d N x %d q-points\n",
            dim(fitted_arr)[1], dim(fitted_arr)[2], dim(fitted_arr)[3],
            dim(fitted_arr)[4]))
cat(sprintf("  finite   : %d / %d cells populated (slice at q-index 251)\n",
            sum(!is.na(fitted_arr[, , , 251L])),
            prod(dim(fitted_arr)[1:3])))
cat(sprintf("  size     : %.1f MB uncompressed in memory\n",
            as.numeric(object.size(fitted_icc_reference_curves)) / 1024^2))

# ---- Merge into sysdata.rda (preserve existing objects) -----------------
assign("fitted_icc_reference_curves", fitted_icc_reference_curves, envir = e)
keep <- ls(e)
cat(sprintf("\nSaving objects to %s: %s\n",
            SYSDATA, paste(keep, collapse = ", ")))
save(list = keep, file = SYSDATA, envir = e, compress = "xz")

cat(sprintf("\nsysdata.rda size after save: %.2f MB\n",
            file.info(SYSDATA)$size / 1024^2))
cat("DONE.\n")
