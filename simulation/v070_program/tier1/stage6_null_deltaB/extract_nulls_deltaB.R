#!/usr/bin/env Rscript
# Stage 6B extract -- pool the Option-B null draws into the ridge ECDF object
# and JSON artifact. Mirrors stage6_production_null/extract_nulls_json.R and
# emits null_ecdf_cells.rds in the SAME per-ridge shape that
# grassr/data-raw/build_delta_null_ecdf.R consumes ($k,$N,$q,$n,$fine plus
# $p99lo/$p99hi for the unstable_tail flag; $med/$p95/$p99/$hist/$h2/$qs feed
# the JSON display artifact).
#
# Difference from 0.7.0 extract: prevalence is already pooled INSIDE each
# (k, N, q) cell (Option-B single-program regen), so each per_cell file is
# one ridge -- no 6b/6c top-up dirs to merge. Pooling by (k,N,q) key is kept
# for robustness (idempotent under any accidental duplicate cell files).
#
# The build script's input path is hardcoded to stage6_production_null (the
# 0.7.0 record, which we do NOT touch). To ship the Option-B null, repoint
# build_delta_null_ecdf.R at THIS file's output, or copy it across, as the
# operator's release step -- see README.md.
#
# Run: Rscript extract_nulls_deltaB.R   (WORKERS respected for the read/pool)

get_script_path <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[[1]])))
  of <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(of)) return(normalizePath(of))
  stop("cannot locate script; run via Rscript")
}
S6 <- dirname(get_script_path())
files <- list.files(file.path(S6, "per_cell"), pattern = "^cell_[0-9]+\\.rds$",
                    full.names = TRUE)
if (!length(files)) stop("no per_cell RDS found under ", file.path(S6, "per_cell"))

read_one <- function(f) tryCatch(readRDS(f), error = function(e) NULL)
cells <- Filter(Negate(is.null), lapply(files, read_one))
key   <- vapply(cells, function(x) paste(x$cell$k, x$cell$N, x$cell$q), character(1L))

BR    <- seq(0, 70, length.out = 57)                          # JSON histogram
PROBS <- seq(0.01, 0.99, by = 0.02)                           # JSON quantile grid
FINE  <- sort(unique(c(seq(0.01, 0.99, by = 0.01), 0.995)))   # shipped fine grid

out <- lapply(split(seq_along(cells), key), function(idx) {
  cell <- cells[[idx[1]]]$cell
  d <- unlist(lapply(idx, function(i) cells[[i]]$draws$delta))
  d <- d[is.finite(d)]
  dd <- pmin(d, 69.99)
  bt <- replicate(200, quantile(sample(d, replace = TRUE), 0.99))
  list(k = cell$k, N = cell$N, q = cell$q, n = length(d),
       med = round(median(d), 2), p95 = round(quantile(d, .95), 2),
       p99 = round(quantile(d, .99), 2),
       p99lo = round(quantile(bt, .025), 2), p99hi = round(quantile(bt, .975), 2),
       hist = round(hist(dd, breaks = BR, plot = FALSE)$counts / length(dd), 5),
       h2max = { hm <- min(70, max(10, quantile(d, .995) * 1.3)); round(hm, 1) },
       h2 = { hm <- min(70, max(10, quantile(d, .995) * 1.3))
              round(hist(pmin(d, hm - 1e-6), breaks = seq(0, hm, length.out = 81),
                         plot = FALSE)$counts / length(d), 5) },
       qs = round(quantile(d, PROBS, names = FALSE), 2),
       fine = round(quantile(d, FINE, names = FALSE), 3))
})

json <- jsonlite::toJSON(unname(out), auto_unbox = TRUE)
writeLines(json, file.path(S6, "nulls_production.json"))
saveRDS(out, file.path(S6, "null_ecdf_cells.rds"))
cat(sprintf("extracted %d ridges (Option B); total finite draws %s; min ridge n = %s\n",
            length(out), format(sum(sapply(out, `[[`, "n")), big.mark = ","),
            format(min(sapply(out, `[[`, "n")), big.mark = ",")))
