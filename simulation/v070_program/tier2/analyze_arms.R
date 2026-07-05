#!/usr/bin/env Rscript
# Tier 2 arms A-C analysis: percentile drift vs each arm's null anchor and
# delta-hat false-flag inflation at the stage-4 own-(k,N,q) divergent cut.
T2 <- "grassr/simulation/v070_program/tier2"
by_q <- readRDS("grassr/simulation/v070_program/tier1/stage4_threshold_table/threshold_table_by_q.rds")
cut_for <- function(k, N, q) {
  r <- by_q[by_q$k==k & by_q$N==N & by_q$q==q, "t_divergent"]
  if (length(r) && is.finite(r)) r else NA_real_
}
load_arm <- function(a) do.call(rbind, lapply(
  list.files(file.path(T2, a, "per_cell"), full.names = TRUE), readRDS))

summarize_null_arm <- function(dat, par) {
  agg <- do.call(rbind, lapply(split(dat, dat[, c(par,"k","N","q"), drop=FALSE], drop=TRUE),
    function(g) {
      cut <- cut_for(g$k[1], g$N[1], g$q[1])
      data.frame(v = g[[par]][1], k = g$k[1], N = g$N[1], q = g$q[1],
                 med_pct_pabak = median(g$pabak.1, na.rm=TRUE),
                 flag_rate = if (is.finite(cut)) mean(g$delta >= cut, na.rm=TRUE) else NA)
    }))
  # marginal over (k,N,q,prev): median |percentile - anchor median| and flag rate by parameter
  out <- do.call(rbind, lapply(split(agg, agg$v), function(g)
    data.frame(value = g$v[1],
               med_pabak_pct = round(median(g$med_pct_pabak, na.rm=TRUE),1),
               mean_false_flag = round(mean(g$flag_rate, na.rm=TRUE),3),
               max_false_flag  = round(max(g$flag_rate,  na.rm=TRUE),3))))
  out
}

for (spec in list(c("arm_a_item_difficulty","sd_d"), c("arm_b_correlated_errors","rho"))) {
  dat <- load_arm(spec[1])
  names(dat) <- sub("^pabak$","pabak.0",names(dat))  # obs cols vs pct cols disambiguation
  # run_arm cbinds t(obs) then t(pcts): columns pabak, fleiss_kappa, mean_ac1 (obs), then pabak.1, fleiss_kappa.1, mean_ac1.1 (pcts)
  cat(sprintf("\n===== %s (parameter: %s; symmetric raters => every flag is FALSE) =====\n", spec[1], spec[2]))
  print(summarize_null_arm(dat, spec[2]), row.names = FALSE)
}

dat <- load_arm("arm_c_asymmetry_patterns")
cat("\n===== arm_c: divergent TPR by asymmetry pattern (A=0.20, q=0.85 cut) =====\n")
agg <- do.call(rbind, lapply(split(dat, dat[, c("pattern","k","N")], drop=TRUE), function(g) {
  cut <- cut_for(g$k[1], g$N[1], 0.85)
  data.frame(pattern=g$pattern[1], k=g$k[1], N=g$N[1],
             tpr = if (is.finite(cut)) mean(g$delta >= cut, na.rm=TRUE) else NA)
}))
out <- do.call(rbind, lapply(split(agg, agg$pattern), function(g)
  data.frame(pattern=g$pattern[1], mean_tpr=round(mean(g$tpr,na.rm=TRUE),3),
             min=round(min(g$tpr,na.rm=TRUE),3), max=round(max(g$tpr,na.rm=TRUE),3))))
print(out, row.names = FALSE)
