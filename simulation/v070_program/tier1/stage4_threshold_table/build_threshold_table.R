#!/usr/bin/env Rscript
# Stage 4 — threshold table + q-flatness verdict + LOO snapping analysis.
#
# Consumes stage-3 per_cell draws (13,200 cells; A=0 rows are the null).
# Derivation matches build_smallN_threshold_table.R exactly: smallest t on
# a 0.25-pp grid with null exceedance <= alpha (0.05 caution / 0.01
# divergent), marginal over prevalence.
#
# Outputs (stage4 dir):
#   threshold_table_by_q.rds     t(k, N, q) — the bias check's raw answer
#   threshold_table_pooled.rds   t(k, N) pooling null draws across q
#   tpr_by_q.rds                 TPR at A = 0.20 by (k, N, q)
#   loo_analysis.rds             interpolation-vs-snap on the dense grid
#   analysis_summary.txt         human-readable verdicts
# DONE is written by this script only if the q-flatness verdict is clean
# (spread <= 1 pp everywhere); otherwise NEEDS_REVIEW — the t(k,N,q)
# lookup decision belongs to Austin.

S3 <- "grassr/simulation/v070_program/tier1/stage3_threshold_grid/per_cell"
S4 <- "grassr/simulation/v070_program/tier1/stage4_threshold_table"
dir.create(S4, showWarnings = FALSE, recursive = TRUE)

files <- list.files(S3, pattern = "^cell_[0-9]+\\.rds$", full.names = TRUE)
cat(sprintf("per_cell files: %d (expect 13200)\n", length(files)))
stopifnot(length(files) == 13200L)
dat <- do.call(rbind, lapply(files, readRDS))
dat <- dat[is.finite(dat$delta), ]
cat(sprintf("finite delta rows: %s\n", format(nrow(dat), big.mark = ",")))

np_threshold <- function(deltas, alpha) {
  thr_grid <- seq(0, 50, by = 0.25)
  fpr <- vapply(thr_grid, function(t) mean(deltas >= t), numeric(1))
  ok <- fpr <= alpha
  if (!any(ok)) return(NA_real_)
  min(thr_grid[ok])
}

null <- dat[dat$A == 0, ]
by_q <- do.call(rbind, lapply(split(null, list(null$k, null$N, null$q), drop = TRUE),
  function(g) data.frame(k = g$k[1], N = g$N[1], q = g$q[1],
                         t_caution = np_threshold(g$delta, 0.05),
                         t_divergent = np_threshold(g$delta, 0.01),
                         n_null = nrow(g))))
by_q <- by_q[order(by_q$k, by_q$N, by_q$q), ]
saveRDS(by_q, file.path(S4, "threshold_table_by_q.rds"))

pooled <- do.call(rbind, lapply(split(null, list(null$k, null$N), drop = TRUE),
  function(g) data.frame(k = g$k[1], N = g$N[1],
                         t_caution = np_threshold(g$delta, 0.05),
                         t_divergent = np_threshold(g$delta, 0.01),
                         n_null = nrow(g))))
pooled <- pooled[order(pooled$k, pooled$N), ]
saveRDS(pooled, file.path(S4, "threshold_table_pooled.rds"))

a20 <- dat[dat$A == 0.20, ]
tpr <- do.call(rbind, lapply(split(a20, list(a20$k, a20$N, a20$q), drop = TRUE),
  function(g) {
    row <- by_q[by_q$k == g$k[1] & by_q$N == g$N[1] & by_q$q == g$q[1], ]
    data.frame(k = g$k[1], N = g$N[1], q = g$q[1],
               tpr_caution = mean(g$delta >= row$t_caution),
               tpr_divergent = mean(g$delta >= row$t_divergent))
  }))
saveRDS(tpr, file.path(S4, "tpr_by_q.rds"))

# ---- q-flatness verdict -----------------------------------------------
flat <- do.call(rbind, lapply(split(by_q, list(by_q$k, by_q$N), drop = TRUE),
  function(g) data.frame(k = g$k[1], N = g$N[1],
    spread_caution = diff(range(g$t_caution, na.rm = TRUE)),
    spread_divergent = diff(range(g$t_divergent, na.rm = TRUE)),
    n_calibrated_q = sum(is.finite(g$t_divergent)))))
flat <- flat[order(-pmax(flat$spread_caution, flat$spread_divergent)), ]
max_spread <- max(pmax(flat$spread_caution, flat$spread_divergent), na.rm = TRUE)

# ---- LOO: interpolation vs snap on the dense pooled grid ---------------
cal <- pooled[is.finite(pooled$t_caution) & is.finite(pooled$t_divergent), ]
loo <- do.call(rbind, lapply(seq_len(nrow(cal)), function(i) {
  tgt <- cal[i, ]; rest <- cal[-i, ]
  # snap: nearest by normalized (k, log10 N)
  dk <- (rest$k - tgt$k) / diff(range(cal$k))
  dN <- (log10(rest$N) - log10(tgt$N)) / diff(range(log10(cal$N)))
  snap <- rest[which.min(sqrt(dk^2 + dN^2)), ]
  # interp: inverse-distance weight of <=4 nearest neighbours
  d <- sqrt(dk^2 + dN^2); w <- 1 / pmax(d, 1e-6)
  nb <- order(d)[seq_len(min(4L, nrow(rest)))]
  interp_c <- sum(rest$t_caution[nb] * w[nb]) / sum(w[nb])
  interp_d <- sum(rest$t_divergent[nb] * w[nb]) / sum(w[nb])
  data.frame(k = tgt$k, N = tgt$N,
             snap_err_c = abs(snap$t_caution - tgt$t_caution),
             snap_err_d = abs(snap$t_divergent - tgt$t_divergent),
             interp_err_c = abs(interp_c - tgt$t_caution),
             interp_err_d = abs(interp_d - tgt$t_divergent))
}))
saveRDS(loo, file.path(S4, "loo_analysis.rds"))

lines <- c(
  sprintf("stage4 summary (%s)", Sys.time()),
  sprintf("calibrated (k,N,q) cells: %d of %d (finite t_divergent)",
          sum(is.finite(by_q$t_divergent)), nrow(by_q)),
  sprintf("calibrated pooled (k,N) cells: %d of %d",
          sum(is.finite(pooled$t_divergent)), nrow(pooled)),
  "",
  sprintf("Q-FLATNESS: max threshold spread across q = %.2f pp", max_spread),
  "worst 8 cells:",
  capture.output(print(head(flat, 8), row.names = FALSE)),
  "",
  "LOO (pooled grid): mean abs err, snap vs interp:",
  sprintf("  caution:   snap %.2f  interp %.2f",
          mean(loo$snap_err_c), mean(loo$interp_err_c)),
  sprintf("  divergent: snap %.2f  interp %.2f",
          mean(loo$snap_err_d), mean(loo$interp_err_d)),
  sprintf("  worst-case: snap %.2f  interp %.2f",
          max(pmax(loo$snap_err_c, loo$snap_err_d)),
          max(pmax(loo$interp_err_c, loo$interp_err_d))),
  "",
  "shipped-vs-new modal cell (k=5, N=200, pooled):",
  capture.output(print(pooled[pooled$k == 5 & pooled$N == 200, ], row.names = FALSE))
)
writeLines(lines, file.path(S4, "analysis_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")

if (is.finite(max_spread) && max_spread <= 1.0) {
  writeLines(sprintf("stage4 %s: q-flat (max spread %.2f pp) — pooled table adopted",
                     Sys.time(), max_spread), file.path(S4, "DONE"))
  cat("DONE: thresholds flat in q; pooled t(k,N) table stands.\n")
} else {
  writeLines(sprintf("stage4 %s: thresholds vary with q (max spread %.2f pp) — t(k,N,q) lookup decision needed",
                     Sys.time(), max_spread), file.path(S4, "NEEDS_REVIEW"))
  cat("NEEDS_REVIEW: q-conditioning decision required.\n")
}
