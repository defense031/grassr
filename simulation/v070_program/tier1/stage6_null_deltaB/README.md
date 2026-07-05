# Stage 6B -- Option-B delta_hat null regeneration

**What.** Regenerates the delta_hat null through the v0.7.1 pipeline. delta_hat
is now the implied-quality spread (PABAK/mean-AC1/Fleiss, quality pp) from
`check_asymmetry()`. Each draw stores delta_hat + the three implied q_hats.

**Why.** Option B (design/v0.7.1_position_redesign.md, ratified 2026-07-05)
redefined delta_hat; the 0.7.0 null stored only old-units delta, so it cannot
be converted -- full re-sim. 440 (k,N,q) ridges x 50,000 draws (10k at each of
5 prevalences, swept inside the cell) = 22,000,000 draws. seed = 20260706 + cell_id.

**Launch** (per machine; build the Option-B candidate lib first, park lme4):
- laptop:    `SPLIT_MOD=2 SPLIT_REM=0 WORKERS=10 CANDIDATE_LIB=/path/lib ../../chain_deltaB.sh`
- Sebastian: `SPLIT_MOD=2 SPLIT_REM=1 WORKERS=8  CANDIDATE_LIB=/path/lib ../../chain_deltaB.sh`
- or run just this stage: `Rscript run_null_deltaB.R` with the same env vars.

**Done.** Each split writes `split_<rem>_of_<mod>.DONE` and a "split R/M COMPLETE"
line in PROGRESS.txt. After merging both machines' `per_cell/`, run
`Rscript extract_nulls_deltaB.R` -> `null_ecdf_cells.rds` (repoint
data-raw/build_delta_null_ecdf.R at it for the release build).
