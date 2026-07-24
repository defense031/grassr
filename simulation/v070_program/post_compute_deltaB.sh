#!/usr/bin/env bash
# post_compute_deltaB.sh -- operator's one-command finisher for the delta-B
# regen. Run on the LAPTOP after monitor_deltaB.sh reports BOTH machines
# COMPLETE (stage6B + stage7B splits DONE on laptop and sebastian).
#
# Safe to re-run at any point: partial rsync pulls resume; each compute step
# skips when its output exists and is newer than its inputs. Everything is
# tee'd to post_compute_deltaB.log.
#
# Steps:
#   1. Pull sebastian's stage6B + stage7B per_cell halves (rsync, NO --delete).
#   2. Verify 440 stage6 + 11,000 stage7 cell files (abort listing missing ids).
#   3. Extract the Option-B null (extract_nulls_deltaB.R -> null_ecdf_cells.rds).
#   4. Analyze power (analyze_power_deltaB.R -> tpr_percentile_convention.rds).
#   5. Pull sebastian's tier2_deltaB per_cell if present; verify 468.
#   6. Print a DONE summary + next manual steps.
#
# STEP-3 HANDOFF (verified by reading both scripts): extract_nulls_deltaB.R
# writes stage6_null_deltaB/null_ecdf_cells.rds, and analyze_power_deltaB.R
# reads THAT file directly -- so extract already produces the object this
# finisher's pipeline consumes; the data-raw builder is NOT run here.
# grassr/data-raw/build_delta_null_ecdf.R is the PACKAGE sysdata builder: it is
# hardcoded to the untouched 0.7.0 stage6_production_null and would need
# repointing to stage6_null_deltaB. That is a manual RELEASE step (rebuild
# sysdata), listed under "next manual steps" below -- deliberately not invoked
# here so this finisher never mutates grassr/R/sysdata.rda.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # v070_program root
LOG="$SCRIPT_DIR/post_compute_deltaB.log"
log() { echo "[$(date '+%F %T')] post_compute: $*" | tee -a "$LOG"; }

REMOTE_ROOT="sebastian:/Users/sebastian/grassr_calc/grassr/simulation/v070_program"

S6="$SCRIPT_DIR/tier1/stage6_null_deltaB"
S7="$SCRIPT_DIR/tier1/stage7_power_deltaB"
T2="$SCRIPT_DIR/tier2_deltaB"
S6_CELLS="$S6/per_cell"; S7_CELLS="$S7/per_cell"; T2_CELLS="$T2/per_cell"
NULL_OBJ="$S6/null_ecdf_cells.rds"
POWER_OBJ="$S7/tpr_percentile_convention.rds"

RSCRIPT="$(command -v Rscript)"
[ -n "$RSCRIPT" ] || { log "FATAL: Rscript not on PATH"; exit 127; }

# needs_run OUTPUT INPUT... : true (0) if OUTPUT missing or any INPUT newer.
needs_run() {
  local out="$1"; shift; local in
  [ -e "$out" ] || return 0
  for in in "$@"; do
    [ -e "$in" ] || continue
    [ -n "$(find "$in" -newer "$out" 2>/dev/null | head -1)" ] && return 0
  done
  return 1
}

# list missing sequential cell ids in DIR for printf-fmt FMT, ids 1..N
missing_ids() {
  local dir="$1" fmt="$2" n="$3" i f miss=""
  for ((i = 1; i <= n; i++)); do
    printf -v f "$fmt" "$i"
    [ -f "$dir/$f" ] || miss+="$i "
  done
  echo "$miss"
}

# ---- 1. pull sebastian's stage6B + stage7B halves ------------------------
log "step 1: rsync sebastian stage6B + stage7B per_cell (no --delete)"
mkdir -p "$S6_CELLS" "$S7_CELLS"
rsync -a --partial "$REMOTE_ROOT/tier1/stage6_null_deltaB/per_cell/" "$S6_CELLS/" \
  2>&1 | tee -a "$LOG" || { log "FATAL: rsync stage6B failed"; exit 1; }
rsync -a --partial "$REMOTE_ROOT/tier1/stage7_power_deltaB/per_cell/" "$S7_CELLS/" \
  2>&1 | tee -a "$LOG" || { log "FATAL: rsync stage7B failed"; exit 1; }

# ---- 2. verify counts ----------------------------------------------------
log "step 2: verify 440 stage6B + 11,000 stage7B cell files"
n6=$(find "$S6_CELLS" -name 'cell_[0-9]*.rds' 2>/dev/null | wc -l | tr -d ' ')
n7=$(find "$S7_CELLS" -name 'cell_[0-9]*.rds' 2>/dev/null | wc -l | tr -d ' ')
log "  stage6B: $n6/440   stage7B: $n7/11000"
abort=0
if [ "$n6" -ne 440 ]; then
  m6="$(missing_ids "$S6_CELLS" 'cell_%04d.rds' 440)"
  log "  ABORT: stage6B missing cell ids: $m6"; abort=1
fi
if [ "$n7" -ne 11000 ]; then
  m7="$(missing_ids "$S7_CELLS" 'cell_%05d.rds' 11000)"
  log "  ABORT: stage7B missing cell ids: $m7"; abort=1
fi
[ "$abort" -eq 0 ] || { log "FATAL: incomplete per_cell -- re-run after both machines COMPLETE."; exit 1; }

# ---- 3. extract the Option-B null ----------------------------------------
if needs_run "$NULL_OBJ" "$S6_CELLS"; then
  log "step 3: extract_nulls_deltaB.R (build null_ecdf_cells.rds)"
  "$RSCRIPT" "$S6/extract_nulls_deltaB.R" 2>&1 | tee -a "$LOG" \
    || { log "FATAL: extract_nulls_deltaB.R failed"; exit 1; }
else
  log "step 3: null_ecdf_cells.rds up to date -- skipping extract"
fi

# ---- 4. analyze power ----------------------------------------------------
if needs_run "$POWER_OBJ" "$NULL_OBJ" "$S7_CELLS"; then
  log "step 4: analyze_power_deltaB.R (TPR of the percentile convention)"
  "$RSCRIPT" "$S7/analyze_power_deltaB.R" 2>&1 | tee -a "$LOG" \
    || { log "FATAL: analyze_power_deltaB.R failed"; exit 1; }
else
  log "step 4: tpr_percentile_convention.rds up to date -- skipping analyze"
fi

# ---- 5. pull + verify tier2_deltaB (skip gracefully if not run there) -----
log "step 5: pull sebastian tier2_deltaB per_cell (if present)"
mkdir -p "$T2_CELLS"
if rsync -a --partial "$REMOTE_ROOT/tier2_deltaB/per_cell/" "$T2_CELLS/" 2>>"$LOG"; then
  :
else
  log "  tier2_deltaB not present on sebastian (or rsync skipped) -- continuing"
fi
nt=$(find "$T2_CELLS" -name 'arm[ABC]_cell_[0-9]*.rds' 2>/dev/null | wc -l | tr -d ' ')
if [ "$nt" -gt 0 ]; then
  na=$(find "$T2_CELLS" -name 'armA_cell_[0-9]*.rds' | wc -l | tr -d ' ')
  nb=$(find "$T2_CELLS" -name 'armB_cell_[0-9]*.rds' | wc -l | tr -d ' ')
  nc=$(find "$T2_CELLS" -name 'armC_cell_[0-9]*.rds' | wc -l | tr -d ' ')
  log "  tier2_deltaB: $nt/468 (armA $na/180, armB $nb/180, armC $nc/108)"
  [ "$nt" -ne 468 ] && log "  WARN: tier2_deltaB incomplete -- run analyze_tier2_deltaB.R only after 468."
else
  log "  tier2_deltaB per_cell empty -- tier2 has not run yet (fine; run it later)."
fi

# ---- 6. DONE summary -----------------------------------------------------
min_ndraws="n/a"
if [ -f "$NULL_OBJ" ]; then
  min_ndraws=$("$RSCRIPT" -e 'x<-readRDS(commandArgs(TRUE)[1]); cat(format(min(vapply(x,function(cl) cl$n, numeric(1))), big.mark=","))' "$NULL_OBJ" 2>/dev/null) || min_ndraws="n/a"
fi

echo    | tee -a "$LOG"
log "================= delta-B post-compute DONE ================="
log "null stage6B: $NULL_OBJ ; min n_draws across ridges = $min_ndraws"
log "power table:    $POWER_OBJ"
log "tier2_deltaB:   $nt/468 cell files present"
log "next manual steps:"
log "  1. rebuild package sysdata: repoint grassr/data-raw/build_delta_null_ecdf.R"
log "     from stage6_production_null to stage6_null_deltaB, then Rscript it"
log "     (rebuilds grassr/R/sysdata.rda delta_null_ecdf; this finisher does NOT)."
log "  2. if tier2_deltaB is 468/468: NULL_CELLS=$NULL_OBJ Rscript $T2/analyze_tier2_deltaB.R"
log "  3. un-skip the delta-hat tests; re-knit the vignettes on the new numbers."
log "  4. build the 0.7.1 tarball."
log "============================================================"
