#!/usr/bin/env Rscript
# TPR of the percentile-flag convention: position each A>0 delta-hat draw
# on the MATCHED null ECDF (stage-6 fine quantile grid), flag at >=95th
# (caution) / >=99th (divergent), report rates per (k, N, q, A) pooled
# over prevalence. Also realized false-positive rates from A=0 nulls are
# by construction ~5%/1% — the size guarantee is definitional here.
setwd("/Users/austinsemmel/Desktop/PABAK_Investigation")
S7 <- "grassr/simulation/v070_program/tier1/stage7_power"
nulls <- readRDS("grassr/simulation/v070_program/tier1/stage6_production_null/null_ecdf_cells.rds")
FINE <- sort(unique(c(seq(0.01, 0.99, by = 0.01), 0.995)))
nl <- new.env()
for (c in nulls) assign(paste(c$k, c$N, c$q), c$fine, nl)

files <- list.files(file.path(S7, "per_cell"), full.names = TRUE)
res <- vector("list", length(files))
for (i in seq_along(files)) {
  x <- readRDS(files[i]); cell <- x$cell
  fine <- get0(paste(cell$k, cell$N, cell$q), envir = nl)
  d <- x$delta[is.finite(x$delta)]
  if (is.null(fine) || length(d) < 100) next
  # percentile via ECDF interpolation on the fine grid
  pct <- approx(x = fine, y = FINE, xout = d, rule = 2, ties = "ordered")$y
  res[[i]] <- data.frame(k = cell$k, N = cell$N, q = cell$q, A = cell$A,
                         prev = cell$prev, n = length(d),
                         tpr95 = mean(pct >= 0.95), tpr99 = mean(pct >= 0.99))
}
res <- do.call(rbind, res)
agg <- aggregate(cbind(tpr95, tpr99) ~ k + N + q + A, data = res, FUN = mean)
saveRDS(agg, file.path(S7, "tpr_percentile_convention.rds"))
cat("== divergent-flag TPR (>=99th pct of matched null), A = 0.20 ==\n")
sub <- agg[agg$A == 0.20 & agg$q == 0.85 & agg$N %in% c(50, 200, 500, 1000), ]
print(sub[order(sub$k, sub$N), c("k","N","tpr95","tpr99")], row.names = FALSE, digits = 2)
cat("\n== TPR by q at (k=8, N=200), A = 0.20 ==\n")
print(agg[agg$k==8 & agg$N==200 & agg$A==0.20, c("q","tpr95","tpr99")], row.names = FALSE, digits = 2)
cat("\n== TPR by A at modal (k=5, N=200, q=0.85) ==\n")
print(agg[agg$k==5 & agg$N==200 & agg$q==0.85, c("A","tpr95","tpr99")], row.names = FALSE, digits = 2)
