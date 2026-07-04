#!/bin/bash
# Waterfall runner for the grassr v0.7.0 calibration program.
#
# Runs tier stages in order; a stage is skipped when its DONE marker
# exists and the waterfall halts at the first gated tier lacking an
# APPROVED file (Austin's explicit go). Designed to be safe to re-run:
# every stage script is internally resumable (per-scenario RDS).
#
# Usage:  ./runner.sh            # run/resume everything permitted
#         WORKERS=8 ./runner.sh  # cap workers
# Cross-machine splits are launched manually per stage (SPLIT_MOD/REM);
# the DONE marker for split stages is written by the merge step, not
# the runner.

set -euo pipefail
PROG="$(cd "$(dirname "$0")" && pwd)"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$PROG/runner.log"; }

run_stage() {  # $1 = stage dir, $2 = command
  local dir="$PROG/$1"
  if [ -f "$dir/DONE" ]; then log "skip $1 (DONE)"; return 0; fi
  log "start $1"
  ( cd "$dir" && eval "$2" ) 2>&1 | tee -a "$dir/stage.log"
  log "finished command for $1 (DONE marker is written by the stage's own verify step where applicable)"
}

# ---------------- Tier 1 ----------------
run_stage "tier1/stage1_surface_ndense"  "Rscript run_surface_ndense.R"
[ -f "$PROG/tier1/stage1_surface_ndense/DONE" ] || { log "halt: stage1 split incomplete (merge + DONE pending)"; exit 0; }
run_stage "tier1/stage2_rebuild_surface" "Rscript rebuild_surface.R"
[ -f "$PROG/tier1/stage2_rebuild_surface/DONE" ] || { log "halt: stage2 verification pending"; exit 0; }
run_stage "tier1/stage3_threshold_grid"  "Rscript run_threshold_grid.R"
[ -f "$PROG/tier1/stage3_threshold_grid/DONE" ] || { log "halt: stage3 split incomplete"; exit 0; }
run_stage "tier1/stage4_threshold_table" "Rscript build_threshold_table.R"
run_stage "tier1/stage5_assemble"        "Rscript assemble_and_verify.R"

# ---------------- Tier 2 (gated) ----------------
if [ ! -f "$PROG/tier2/APPROVED" ]; then
  log "halt: tier2 gated (create tier2/APPROVED after ratifying its scripts)"
  exit 0
fi
log "tier2 approved - add stage invocations here once scripts are ratified"

# ---------------- Tier 3 (gated) ----------------
if [ ! -f "$PROG/tier3/APPROVED" ]; then
  log "halt: tier3 gated"
  exit 0
fi
log "tier3 approved - add stage invocations here once scripts are ratified"
