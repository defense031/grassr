#!/usr/bin/env Rscript
# Pool stage-6 production null draws over prevalence and emit the artifact
# data (nulls.json format: 56-bin histogram to 70pp, 50-point quantile
# grid, median/p95/p99 + bootstrap CI) plus the shipped-object precursor
# (fine tail grid) per (k, N, q).
S6 <- "grassr/simulation/v070_program/tier1/stage6_production_null"
files <- list.files(file.path(S6, "per_cell"), full.names = TRUE)
cells <- lapply(files, readRDS)
key <- sapply(cells, function(x) paste(x$cell$k, x$cell$N, x$cell$q))
BR <- seq(0, 70, length.out = 57)
PROBS <- seq(0.01, 0.99, by = 0.02)
FINE <- sort(unique(c(seq(0.01, 0.99, by = 0.01), 0.995)))  # shipped grid
out <- lapply(split(seq_along(cells), key), function(idx) {
  cell <- cells[[idx[1]]]$cell
  d <- unlist(lapply(idx, function(i) cells[[i]]$delta))
  d <- d[is.finite(d)]
  dd <- pmin(d, 69.99)
  bt <- replicate(200, quantile(sample(d, replace = TRUE), 0.99))
  list(k = cell$k, N = cell$N, q = cell$q, n = length(d),
       med = round(median(d), 2), p95 = round(quantile(d, .95), 2),
       p99 = round(quantile(d, .99), 2),
       p99lo = round(quantile(bt, .025), 2), p99hi = round(quantile(bt, .975), 2),
       hist = round(hist(dd, breaks = BR, plot = FALSE)$counts / length(dd), 5),
       qs = round(quantile(d, PROBS, names = FALSE), 2),
       fine = round(quantile(d, FINE, names = FALSE), 3))
})
json <- jsonlite::toJSON(unname(out), auto_unbox = TRUE)
writeLines(json, file.path(S6, "nulls_production.json"))
saveRDS(out, file.path(S6, "null_ecdf_cells.rds"))
cat(sprintf("extracted %d ridges; total finite draws %s\n",
    length(out), format(sum(sapply(out, `[[`, "n")), big.mark = ",")))
