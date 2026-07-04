#!/usr/bin/env Rscript
# Tier 2 Arm B — correlated errors (TIER2_DESIGN.md).
# Exchangeable error correlation via Gaussian-threshold copula: shared
# item shock z_i; rater j errs when sqrt(rho) z_i + sqrt(1-rho) w_ij
# < qnorm(1 - q). rho = 0 recovers the calibration DGP.
source("grassr/simulation/v070_program/tier2/arm_common.R")
GRID <- expand.grid(
  rho  = c(0, 0.1, 0.25, 0.5),
  kN   = c("3_50", "5_200", "8_100", "10_500", "25_1000"),
  q    = c(0.75, 0.85, 0.92),
  prev = c(0.10, 0.30, 0.50),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
GRID$k <- as.integer(sub("_.*", "", GRID$kN))
GRID$N <- as.integer(sub(".*_", "", GRID$kN))
GRID$kN <- NULL
GRID$cell_id <- seq_len(nrow(GRID))
gen_panel <- function(row) {
  C <- rbinom(row$N, 1L, row$prev)
  z <- rnorm(row$N)
  thr <- qnorm(1 - row$q)   # error prob = 1 - q on both margins
  Y <- matrix(0L, row$N, row$k)
  for (j in seq_len(row$k)) {
    err <- (sqrt(row$rho) * z + sqrt(1 - row$rho) * rnorm(row$N)) < thr
    Y[, j] <- ifelse(err, 1L - C, C)
  }
  Y
}
run_arm("arm_b_correlated_errors", GRID, gen_panel)
