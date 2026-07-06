# Build the delta_null_ecdf sysdata object from the Stage-6 null program.
# Replaces delta_thresholds_lookup as the primary flag object: delta-hat
# is reported as its percentile on the matched (k, N, q) null ECDF; flags
# are conventions on that percentile (>= 95th caution, >= 99th divergent).
#
# v0.7.1 (Option B): delta-hat is the implied-quality spread in QUALITY
# pp; the source is the stage6_null_deltaB regeneration. Override with
# env NULL_CELLS to build from another extraction.
#
# Structure:
#   $probs   fine probability grid (1..99% by 1% + 99.5%)
#   $values  matrix [cells x probs] of delta-hat quantiles (quality pp)
#   $index   data.frame(k, N, q, n_draws, unstable_tail) — unstable_tail
#            marks ridges whose p99 bootstrap CI exceeded the tail-
#            stability tolerance (quantile-INVERSION instability; the
#            ECDF percentile itself is MC-stable)
#
# Run from repo root after the stage-6B extraction:
#   Rscript grassr/data-raw/build_delta_null_ecdf.R

src <- Sys.getenv("NULL_CELLS",
  "grassr/simulation/v070_program/tier1/stage6_null_deltaB/null_ecdf_cells.rds")
if (!file.exists(src))
  stop("null cells not found at ", src,
       " — run the stage-6B extraction first (or set NULL_CELLS).")
cells <- readRDS(src)
# k = 2 ridges are excluded from the shipped object: the two-coefficient
# family implies identical quality by construction, so the k = 2 null is
# a point mass at zero (verified on 2.75M Option-B draws) and delta_hat
# is reported not_applicable there (lookup_delta_null guards k < 3).
cells <- Filter(function(c) c$k >= 3, cells)
FINE <- sort(unique(c(seq(0.01, 0.99, by = 0.01), 0.995)))
# Tail-stability tolerance in the OBJECT'S OWN UNITS: quality-pp nulls are
# ~3 orders of magnitude smaller than the retired percentile-pp nulls, so
# the old absolute 1.0pp CI tolerance would never flag; scale it to the
# ridge's p99 magnitude (10% relative, floored at MC noise).
idx <- do.call(rbind, lapply(cells, function(c)
  data.frame(k = c$k, N = c$N, q = c$q, n_draws = c$n,
             unstable_tail = (c$p99hi - c$p99lo) >
               pmax(0.10 * c$p99, 1e-4))))
# Rounding must preserve relative precision of TINY quality-pp quantiles
# (tight ridges have medians ~0.001-0.01 pp): signif, not fixed decimals
# (the old round(., 3) was calibrated to percentile-pp magnitudes and
# would collapse tight ridges to a handful of distinct values).
vals <- do.call(rbind, lapply(cells, function(c) signif(c$fine, 5)))
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
