# Tier 2 delta-B re-run (Option-B / v0.7.1 pipeline)
WHAT/WHY: re-runs Tier-2 Arms A/B/C through Option-B `check_asymmetry()` — the
0.7.0 outputs lack per-rep pi-hat, so recompute is impossible (TIER2_DESIGN.md,
"Arms A/B/C under v0.7.1"). DGPs copied VERBATIM from `tier2/arm_*.R`; only
pipeline + storage (per-rep implied q_hats, pooled percentiles) modernize.
GRIDS (exact 0.7.0 counts): A 180 + B 180 + C 108 = 468 cells x 1,000 reps.
Seed = 20260708 + arm_offset(0/1000/2000) + cell_id.

LAUNCH (per machine; chained behind stage-7B by `../chain_tier2_deltaB.sh`):
  SPLIT_MOD=2 SPLIT_REM={0|1} WORKERS=10 CANDIDATE_LIB=/path/to/optionB/lib \
    Rscript run_tier2_deltaB.R          # laptop REM=0, sebastian REM=1
ANALYZE (per_cell merged + extract_nulls_deltaB.R done):
  NULL_CELLS=../tier1/stage6_null_deltaB/null_ecdf_cells.rds Rscript analyze_tier2_deltaB.R
Headline flag/TPR rates use PRODUCTION null-cell selection (per-rep q_hat_panel =
median implied q_hat -> nearest calibrated q); *_design = oracle true-q comparison.
