#!/usr/bin/env Rscript
# Tier 2 delta-B analysis -- reproduces the 0.7.0 arm deliverables in NEW units.
# =====================================================================
# Reads tier2_deltaB/per_cell (armA/armB/armC cell files from
# run_tier2_deltaB.R) and produces, restated for the Option-B pipeline
# (TIER2_DESIGN.md "Arms A/B/C under v0.7.1"):
#
#   (A/B) POOLED-PERCENTILE DRIFT vs the null-anchor level of each arm's
#         nuisance parameter (sd_d = 0 for A, rho = 0 for B): for each design
#         cell (k, N, q, prev) and coefficient, the shift in the median pooled
#         percentile relative to the anchor level; summarized as median and p95
#         |drift| by nuisance level. Plus DELTA-B FALSE-FLAG RATES by nuisance
#         level READ AGAINST THE NEW NULL -- the >=95th (caution) / >=99th
#         (divergent) rate of each rep's delta_hat positioned on the matched
#         Option-B null ridge (does item difficulty / correlated error fake
#         divergence? the anchor level is the ~5% / ~1% baseline).
#
#         HEADLINE null-cell selection replicates PRODUCTION semantics
#         (check_asymmetry -> lookup_delta_null): the practitioner never knows
#         true q, so per rep q_hat_panel = median of the finite implied q_hats
#         (q_pabak / q_fleiss_kappa / q_mean_ac1, stored by the driver), snapped
#         to the NEAREST calibrated q level (which.min over the null's global q
#         levels, same as lookup_delta_null). k and N match exactly -- every
#         tier2 (k, N) sits on the null grid, so production's normalized
#         (k, log10 N) snap is the identity here. A design-q-matched rate
#         (*_design columns) is kept as a comparison read: the gap between the
#         two isolates the q-hat-selection-noise component of the flag rate.
#
#   (C)   CAUTION / DIVERGENT TPR by asymmetry pattern (single / half / graded)
#         at A = 0.20 -- does concentrating the asymmetry help or hide it?
#
# The null ridge object is the Option-B extract (stage6_null_deltaB/
# null_ecdf_cells.rds). Point NULL_CELLS at it, or pass its path as arg 1.
# If it is absent (stage 6B not yet extracted), the drift and coefficient
# summaries still compute; the flag-rate deliverables degrade with a clear
# message.
#
# Run: Rscript analyze_tier2_deltaB.R  [null_ecdf_cells.rds]

get_script_path <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[[1]])))
  of <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(of)) return(normalizePath(of))
  stop("cannot locate script; run via Rscript")
}
SCRIPT_DIR <- dirname(get_script_path())              # .../tier2_deltaB
PER_CELL   <- file.path(SCRIPT_DIR, "per_cell")
RESULTS    <- file.path(SCRIPT_DIR, "results")
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
V070       <- dirname(dirname(SCRIPT_DIR))             # .../v070_program

# ---- locate the Option-B null ridge object (env NULL_CELLS or arg 1) ------
args    <- commandArgs(trailingOnly = TRUE)
NULLF   <- if (length(args) >= 1L && nzchar(args[1L])) args[1L] else
           Sys.getenv("NULL_CELLS",
             file.path(V070, "tier1", "stage6_null_deltaB", "null_ecdf_cells.rds"))
HAVE_NULL <- file.exists(NULLF)

FINE <- sort(unique(c(seq(0.01, 0.99, by = 0.01), 0.995)))   # shipped fine grid
nl   <- new.env()          # "k N q" -> fine quantile grid
NULL_QS <- numeric(0)      # global calibrated q levels (production snap domain)
if (HAVE_NULL) {
  nulls <- readRDS(NULLF)
  for (cl in nulls) assign(paste(cl$k, cl$N, cl$q), cl$fine, nl)
  NULL_QS <- sort(unique(vapply(nulls, `[[`, numeric(1L), "q")))
  cat(sprintf("null ridges loaded: %d (%s); calibrated q levels: %s\n",
              length(nulls), NULLF, paste(NULL_QS, collapse = ", ")))
} else {
  cat(sprintf("[warn] Option-B null not found at %s -- flag-rate deliverables ",
              NULLF),
      "will be skipped (run extract_nulls_deltaB.R first). Drift deliverables ",
      "still computed.\n", sep = "")
}

# ---- load one arm's per_cell files into a long per-rep data.frame ----------
load_arm <- function(prefix) {
  files <- list.files(PER_CELL, pattern = sprintf("^%s_cell_[0-9]+\\.rds$", prefix),
                      full.names = TRUE)
  if (!length(files)) return(NULL)
  parts <- lapply(files, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x) || is.null(x$draws)) return(NULL)
    cbind(x$cell[rep(1L, nrow(x$draws)), , drop = FALSE], x$draws, row.names = NULL)
  })
  do.call(rbind, Filter(Negate(is.null), parts))
}

# ---- position delta_hat on the matched null; add flags --------------------
# HEADLINE (flag95/flag99): PRODUCTION cell selection. Per rep,
# q_hat_panel = median of the finite implied q_hats, snapped to the nearest
# calibrated q (which.min over the null's global q levels, mirroring
# lookup_delta_null; ties resolve to the lower q, same as which.min). k and N
# match exactly (all tier2 (k, N) are on the null grid, so production's
# normalized (k, log10 N) snap is the identity).
# COMPARISON (flag95_design/flag99_design): oracle design-q matching (true q of
# the DGP). The headline-minus-design gap isolates q-hat-selection noise.
percentile_on <- function(fine, d) {
  ok  <- is.finite(d)
  pct <- rep(NA_real_, length(d))
  if (any(ok))
    pct[ok] <- approx(x = fine, y = FINE, xout = d[ok],
                      rule = 2, ties = "ordered")$y
  pct
}
add_flags <- function(df) {
  df$flag95 <- NA_real_; df$flag99 <- NA_real_               # headline (q-hat)
  df$flag95_design <- NA_real_; df$flag99_design <- NA_real_ # comparison (true q)
  if (!HAVE_NULL) return(df)

  # -- headline: per-rep q_hat_panel -> nearest calibrated q (production) ----
  qmat <- cbind(df$q_pabak, df$q_fleiss_kappa, df$q_mean_ac1)
  q_hat_panel <- apply(qmat, 1L, function(z) {
    z <- z[is.finite(z)]; if (length(z)) median(z) else NA_real_ })
  q_snap <- rep(NA_real_, nrow(df))
  okq <- is.finite(q_hat_panel)
  q_snap[okq] <- vapply(q_hat_panel[okq],
                        function(z) NULL_QS[which.min(abs(NULL_QS - z))],
                        numeric(1L))
  keys_hat <- paste(df$k, df$N, q_snap)
  for (ky in unique(keys_hat[okq])) {
    fine <- get0(ky, envir = nl, inherits = FALSE)
    if (is.null(fine)) next                     # ridge absent -> leave NA
    idx <- which(keys_hat == ky)
    pct <- percentile_on(fine, df$delta[idx])
    df$flag95[idx] <- as.numeric(pct >= 0.95)
    df$flag99[idx] <- as.numeric(pct >= 0.99)
  }

  # -- comparison: oracle true-q matching ------------------------------------
  keys_des <- paste(df$k, df$N, df$q)
  for (ky in unique(keys_des)) {
    fine <- get0(ky, envir = nl, inherits = FALSE)
    if (is.null(fine)) next
    idx <- which(keys_des == ky)
    pct <- percentile_on(fine, df$delta[idx])
    df$flag95_design[idx] <- as.numeric(pct >= 0.95)
    df$flag99_design[idx] <- as.numeric(pct >= 0.99)
  }
  df
}

PP_COLS <- c(pabak = "pp_pabak", fleiss_kappa = "pp_fleiss_kappa",
             mean_ac1 = "pp_mean_ac1")

# ---- (A/B) pooled-percentile drift vs the anchor (nuis == 0) --------------
# For each coefficient and design cell (k,N,q,prev): median pooled percentile
# at each nuisance level minus the same at the anchor level; then median and
# p95 of |drift| across design cells, by nuisance level.
drift_summary <- function(df, nuis) {
  design <- c("k", "N", "q", "prev")
  out <- list()
  for (co in names(PP_COLS)) {
    v <- df[[PP_COLS[[co]]]]
    grp <- df[c(design, nuis)]
    med <- aggregate(list(med_pp = v), by = grp,
                     FUN = function(z) median(z[is.finite(z)]))
    anch <- med[med[[nuis]] == 0, c(design, "med_pp")]
    names(anch)[names(anch) == "med_pp"] <- "anchor_pp"
    m <- merge(med, anch, by = design)
    m$abs_drift <- abs(m$med_pp - m$anchor_pp)
    s <- aggregate(list(abs_drift = m$abs_drift), by = m[nuis],
                   FUN = function(z) c(median = median(z, na.rm = TRUE),
                                       p95 = as.numeric(quantile(z, 0.95, na.rm = TRUE))))
    res <- data.frame(coefficient = co, level = s[[nuis]],
                      median_abs_drift = s$abs_drift[, "median"],
                      p95_abs_drift    = s$abs_drift[, "p95"],
                      stringsAsFactors = FALSE)
    names(res)[names(res) == "level"] <- nuis
    out[[co]] <- res
  }
  do.call(rbind, out)
}

# ---- (A/B) false-flag rate by nuisance level (against the new null) -------
# rate_caution / rate_divergent are the HEADLINE production-semantics rates
# (q-hat-matched null cell); *_design are the oracle true-q comparison.
flag_by_level <- function(df, nuis) {
  if (!HAVE_NULL) return(NULL)
  lv <- sort(unique(df[[nuis]]))
  res <- do.call(rbind, lapply(lv, function(L) {
    s <- df[df[[nuis]] == L, ]
    data.frame(level = L, n = sum(is.finite(s$flag95)),
               rate_caution          = mean(s$flag95, na.rm = TRUE),
               rate_divergent        = mean(s$flag99, na.rm = TRUE),
               rate_caution_design   = mean(s$flag95_design, na.rm = TRUE),
               rate_divergent_design = mean(s$flag99_design, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  names(res)[names(res) == "level"] <- nuis
  res
}

# ---- (C) caution / divergent TPR by asymmetry pattern ---------------------
# tpr_caution / tpr_divergent: headline production-semantics (q-hat-matched);
# *_design: oracle true-q comparison.
tpr_by_pattern <- function(df) {
  if (!HAVE_NULL) return(NULL)
  agg_fun <- function(z) mean(z, na.rm = TRUE)
  byp <- aggregate(cbind(tpr_caution = flag95, tpr_divergent = flag99,
                         tpr_caution_design = flag95_design,
                         tpr_divergent_design = flag99_design) ~ pattern,
                   data = df, FUN = agg_fun, na.action = na.pass)
  byp_kn <- aggregate(cbind(tpr_caution = flag95, tpr_divergent = flag99) ~
                        pattern + k + N, data = df,
                      FUN = agg_fun, na.action = na.pass)
  list(overall = byp, by_kN = byp_kn)
}

# =========================================================================
armA <- load_arm("armA"); armB <- load_arm("armB"); armC <- load_arm("armC")
results <- list()
lines   <- c("== Tier 2 delta-B analysis (Option-B pipeline) ==",
             sprintf("null ridges: %s", if (HAVE_NULL) NULLF else "ABSENT -- flag rates skipped"),
             "")

summarize_ab <- function(df, arm, nuis, results, lines) {
  if (is.null(df)) {
    lines <- c(lines, sprintf("Arm %s: no per_cell files found -- skipped.", arm), "")
    return(list(results = results, lines = lines))
  }
  df <- add_flags(df)
  dr <- drift_summary(df, nuis)
  fl <- flag_by_level(df, nuis)
  results[[sprintf("arm%s_drift", arm)]]      <- dr
  results[[sprintf("arm%s_flag_rate", arm)]]  <- fl
  lines <- c(lines,
    sprintf("-- Arm %s (%s): pooled-percentile drift vs %s=0 anchor --", arm, nuis, nuis),
    capture.output(print(dr, row.names = FALSE, digits = 3)), "")
  if (!is.null(fl))
    lines <- c(lines,
      sprintf("-- Arm %s: delta-B false-flag rate by %s (q-hat-matched null, production semantics; *_design = oracle true-q comparison) --", arm, nuis),
      capture.output(print(fl, row.names = FALSE, digits = 3)), "")
  list(results = results, lines = lines)
}

o <- summarize_ab(armA, "A", "sd_d", results, lines); results <- o$results; lines <- o$lines
o <- summarize_ab(armB, "B", "rho",  results, lines); results <- o$results; lines <- o$lines

if (is.null(armC)) {
  lines <- c(lines, "Arm C: no per_cell files found -- skipped.", "")
} else {
  armC <- add_flags(armC)
  tp <- tpr_by_pattern(armC)
  results[["armC_tpr"]] <- tp
  if (!is.null(tp)) {
    lines <- c(lines, "-- Arm C: caution/divergent TPR by asymmetry pattern (A=0.20; q-hat-matched null, production semantics; *_design = oracle true-q) --",
      capture.output(print(tp$overall, row.names = FALSE, digits = 3)), "",
      "-- Arm C: TPR by pattern x (k, N) --",
      capture.output(print(tp$by_kN[order(tp$by_kN$k, tp$by_kN$N), ],
                           row.names = FALSE, digits = 3)), "")
  } else {
    lines <- c(lines, "Arm C: null absent -- pattern TPR skipped.", "")
  }
}

saveRDS(results, file.path(RESULTS, "tier2_deltaB_results.rds"))
writeLines(lines, file.path(RESULTS, "tier2_deltaB_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
cat(sprintf("\nwrote %s and %s\n",
            file.path(RESULTS, "tier2_deltaB_results.rds"),
            file.path(RESULTS, "tier2_deltaB_summary.txt")))
