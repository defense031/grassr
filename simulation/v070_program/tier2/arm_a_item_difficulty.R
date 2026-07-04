#!/usr/bin/env Rscript
# Tier 2 Arm A — item difficulty (TIER2_DESIGN.md).
# d_i ~ N(0, sd_d^2); Se_ij = plogis(qlogis(q) - d_i), same for Sp.
# sd_d = 0 recovers the calibration DGP (null anchor).
source("grassr/simulation/v070_program/tier2/arm_common.R")
GRID <- expand.grid(
  sd_d = c(0, 0.5, 1.0, 1.5),
  kN   = c("3_50", "5_200", "8_100", "10_500", "25_1000"),
  q    = c(0.75, 0.85, 0.92),
  prev = c(0.10, 0.30, 0.50),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
GRID$k <- as.integer(sub("_.*", "", GRID$kN))
GRID$N <- as.integer(sub(".*_", "", GRID$kN))
GRID$kN <- NULL
GRID$cell_id <- seq_len(nrow(GRID))
gen_panel <- function(row) {
  d  <- rnorm(row$N, 0, row$sd_d)
  Se <- plogis(qlogis(row$q) - d)   # per-item accuracy, all raters
  C  <- rbinom(row$N, 1L, row$prev)
  p  <- ifelse(C == 1L, Se, 1 - Se)  # Sp_i = Se_i under symmetric raters
  Y  <- matrix(0L, row$N, row$k)
  for (j in seq_len(row$k)) Y[, j] <- rbinom(row$N, 1L, ifelse(C == 1L, Se, 1 - Se))
  Y
}
run_arm("arm_a_item_difficulty", GRID, gen_panel)
