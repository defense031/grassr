#!/bin/bash
set -uo pipefail
cd /Users/austinsemmel/Desktop/PABAK_Investigation
S7="grassr/simulation/v070_program/tier1/stage7_power"
log() { echo "[$(date '+%F %T')] chain7: $*"; }
log "waiting for both power halves"
while true; do
  l=$(grep -c COMPLETE "$S7/PROGRESS.txt" 2>/dev/null || true); [ -z "$l" ] && l=0
  s=$(ssh -o ConnectTimeout=20 sebastian "grep -c COMPLETE ~/grassr_calc/$S7/PROGRESS.txt 2>/dev/null" 2>/dev/null | tr -d ' '); [ -z "$s" ] && s=0
  [ "$l" -ge 1 ] && [ "$s" -ge 1 ] && break
  sleep 600
done
log "both halves COMPLETE — merging"
rsync -a "sebastian:~/grassr_calc/$S7/per_cell/" "$S7/per_cell/"
n=$(ls "$S7/per_cell" | wc -l | tr -d ' ')
log "merged: $n cells (expect 11000)"
log "restoring Sebastian's lme4"
ssh sebastian 'LIB="/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"; [ -d "$LIB/lme4.stage6-parked" ] && mv "$LIB/lme4.stage6-parked" "$LIB/lme4" && echo restored'
log "STAGE 7 COMPLETE — TPR analysis under percentile convention next"
