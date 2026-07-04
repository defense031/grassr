#!/usr/bin/env Rscript
# Tier 2 Arm C — asymmetry patterns at fixed mean-norm A = 0.20.
# Patterns: half (current default: random +/- A/2 split), single (all
# asymmetry mass on rater 1, clipped; realized mean-norm recorded),
# ramp (linear gradient across raters, scaled to mean|offset| = A/2).
source("grassr/simulation/v070_program/tier2/arm_common.R")
A_FIX <- 0.20
GRID <- expand.grid(
  pattern = c("half", "single", "ramp"),
  k    = c(3L, 5L, 8L, 10L, 15L, 25L),
  N    = c(50L, 200L, 1000L),
  prev = c(0.20, 0.50),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
GRID$q <- 0.85
GRID$cell_id <- seq_len(nrow(GRID))
gen_panel <- function(row) {
  k <- row$k; q <- row$q
  off <- switch(row$pattern,
    half   = sample(c(-1, 1), k, replace = TRUE) * A_FIX / 2,
    single = c(k * A_FIX / 2, rep(0, k - 1L)),
    ramp   = { v <- seq(-1, 1, length.out = k); v * (A_FIX / 2) / mean(abs(v)) })
  Se <- pmin(pmax(q + off, 0.001), 0.999)
  Sp <- pmin(pmax(q - off, 0.001), 0.999)
  C  <- rbinom(row$N, 1L, row$prev)
  Y  <- matrix(0L, row$N, k)
  for (j in seq_len(k)) Y[, j] <- rbinom(row$N, 1L, ifelse(C == 1L, Se[j], 1 - Sp[j]))
  Y
}
run_arm("arm_c_asymmetry_patterns", GRID, gen_panel)
