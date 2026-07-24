#!/usr/bin/env Rscript
# Stage 7B analysis -- TPR of the percentile-flag convention under Option B.
# Position each A>0 delta_hat draw on the MATCHED Option-B null ridge
# (stage6_null_deltaB fine quantile grid), flag at >= 95th (caution) /
# >= 99th (divergent), report rates per (k, N, q, A) pooled over prevalence.
# Mirrors stage7_power/analyze_power.R; only the null source (Option-B) and
# the per-cell draw field ($draws$delta) differ.
#
# Requires extract_nulls_deltaB.R to have produced the Option-B ridge file.
# Run: Rscript analyze_power_deltaB.R

get_script_path <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[[1]])))
  of <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(of)) return(normalizePath(of))
  stop("cannot locate script; run via Rscript")
}
S7    <- dirname(get_script_path())
TIER1 <- dirname(S7)
NULLF <- file.path(TIER1, "stage6_null_deltaB", "null_ecdf_cells.rds")
if (!file.exists(NULLF))
  stop("Option-B null not found: ", NULLF, " -- run extract_nulls_deltaB.R first.")

nulls <- readRDS(NULLF)
# k = 2 is not calibrated (the two-coefficient family implies identical
# quality by construction; delta_hat reports not_applicable there) -- drop
# the degenerate all-zero ridges so no TPR row is computed at k = 2.
nulls <- Filter(function(c) c$k >= 3, nulls)
FINE  <- sort(unique(c(seq(0.01, 0.99, by = 0.01), 0.995)))
nl <- new.env()
for (c in nulls) assign(paste(c$k, c$N, c$q), c$fine, nl)

files <- list.files(file.path(S7, "per_cell"), pattern = "^cell_[0-9]+\\.rds$",
                    full.names = TRUE)
res <- vector("list", length(files))
for (i in seq_along(files)) {
  x <- tryCatch(readRDS(files[i]), error = function(e) NULL)
  if (is.null(x)) next
  cell <- x$cell
  fine <- get0(paste(cell$k, cell$N, cell$q), envir = nl)
  d <- x$draws$delta[is.finite(x$draws$delta)]
  if (is.null(fine) || length(d) < 100) next
  # percentile via ECDF interpolation on the matched null's fine grid
  pct <- approx(x = fine, y = FINE, xout = d, rule = 2, ties = "ordered")$y
  res[[i]] <- data.frame(k = cell$k, N = cell$N, q = cell$q, A = cell$A,
                         prev = cell$prev, n = length(d),
                         tpr95 = mean(pct >= 0.95), tpr99 = mean(pct >= 0.99))
}
res <- do.call(rbind, res)
agg <- aggregate(cbind(tpr95, tpr99) ~ k + N + q + A, data = res, FUN = mean)
saveRDS(agg, file.path(S7, "tpr_percentile_convention.rds"))

cat("== Option-B divergent-flag TPR (>=99th pct of matched null), A = 0.20 ==\n")
sub <- agg[agg$A == 0.20 & agg$q == 0.85 & agg$N %in% c(50, 200, 500, 1000), ]
print(sub[order(sub$k, sub$N), c("k", "N", "tpr95", "tpr99")], row.names = FALSE, digits = 2)
cat("\n== TPR by q at (k=8, N=200), A = 0.20 ==\n")
print(agg[agg$k == 8 & agg$N == 200 & agg$A == 0.20, c("q", "tpr95", "tpr99")],
      row.names = FALSE, digits = 2)
cat("\n== TPR by A at modal (k=5, N=200, q=0.85) ==\n")
print(agg[agg$k == 5 & agg$N == 200 & agg$q == 0.85, c("A", "tpr95", "tpr99")],
      row.names = FALSE, digits = 2)
