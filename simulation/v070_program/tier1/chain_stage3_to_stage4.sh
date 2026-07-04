#!/bin/bash
# Waits for both stage-3 halves, merges Sebastian's per_cell, runs stage 4.
set -uo pipefail
cd /Users/austinsemmel/Desktop/PABAK_Investigation
S3="grassr/simulation/v070_program/tier1/stage3_threshold_grid"
log() { echo "[$(date '+%F %T')] chain34: $*"; }

log "waiting for both stage-3 halves"
while true; do
  l=$(grep -c COMPLETE "$S3/PROGRESS.txt" 2>/dev/null || true); [ -z "$l" ] && l=0
  s=$(ssh -o ConnectTimeout=20 sebastian "grep -c COMPLETE ~/grassr_calc/$S3/PROGRESS.txt 2>/dev/null" 2>/dev/null | tr -d ' '); [ -z "$s" ] && s=0
  [ "$l" -ge 1 ] && [ "$s" -ge 1 ] && break
  sleep 300
done
log "both halves COMPLETE — merging per_cell"
rsync -a "sebastian:~/grassr_calc/$S3/per_cell/" "$S3/per_cell/"
n=$(ls "$S3/per_cell" | wc -l | tr -d ' ')
log "merged: $n cells (expect 13200)"
[ "$n" -ne 13200 ] && { log "COUNT MISMATCH — halting"; exit 1; }
log "stage 4 build starting"
Rscript grassr/simulation/v070_program/tier1/stage4_threshold_table/build_threshold_table.R \
  || { log "STAGE 4 FAILED"; exit 1; }
log "stage 4 complete — see analysis_summary.txt + DONE/NEEDS_REVIEW"
