# Build the delta_null_ecdf sysdata object (v0.7.0) from the Stage-6
# production-faithful null program. Replaces delta_thresholds_lookup as
# the primary flag object: delta-hat is reported as its percentile on
# the matched (k, N, q) null ECDF; flags are conventions on that
# percentile (>= 95th caution, >= 99th divergent).
#
# Structure:
#   $probs   fine probability grid (1..99% by 1% + 99.5%)
#   $values  matrix [cells x probs] of delta-hat quantiles (pp)
#   $index   data.frame(k, N, q, n_draws, unstable_tail) — unstable_tail
#            marks ridges whose p99 bootstrap CI exceeded 1pp at the
#            100k-draw cap (quantile-INVERSION instability; the ECDF
#            percentile itself is MC-stable to <= 0.5pp everywhere)
#
# Run from repo root after the stage-6 extraction:
#   Rscript grassr/data-raw/build_delta_null_ecdf.R

cells <- readRDS("grassr/simulation/v070_program/tier1/stage6_production_null/null_ecdf_cells.rds")
FINE <- sort(unique(c(seq(0.01, 0.99, by = 0.01), 0.995)))
idx <- do.call(rbind, lapply(cells, function(c)
  data.frame(k = c$k, N = c$N, q = c$q, n_draws = c$n,
             unstable_tail = (c$p99hi - c$p99lo) > 1.0)))
vals <- do.call(rbind, lapply(cells, function(c) round(c$fine, 3)))
ord <- order(idx$k, idx$N, idx$q)
idx <- idx[ord, ]; rownames(idx) <- NULL
vals <- vals[ord, , drop = FALSE]
dimnames(vals) <- list(NULL, as.character(FINE))

delta_null_ecdf <- list(
  probs = FINE, values = vals, index = idx,
  flag_conventions = c(caution = 0.95, divergent = 0.99),
  source = "v0.7.0 stage-6 production-faithful null program (check_asymmetry pipeline)",
  built_on = as.character(Sys.Date()))

e <- new.env()
load("grassr/R/sysdata.rda", envir = e)  # candidate assembly happens at Stage 5;
assign("delta_null_ecdf", delta_null_ecdf, e)
save(list = ls(e), envir = e, file = "grassr/R/sysdata.rda", compress = "xz")
cat(sprintf("delta_null_ecdf: %d cells x %d probs; unstable_tail: %d; sysdata now %.2f MB\n",
    nrow(idx), length(FINE), sum(idx$unstable_tail),
    file.size("grassr/R/sysdata.rda") / 1e6))
