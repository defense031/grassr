#!/bin/bash
set -uo pipefail
cd /Users/austinsemmel/Desktop/PABAK_Investigation
S6B="grassr/simulation/v070_program/tier1/stage6b_topup"
log() { echo "[$(date '+%F %T')] chain6b: $*"; }
log "waiting for both top-up halves"
while true; do
  l=$(grep -c COMPLETE "$S6B/PROGRESS.txt" 2>/dev/null || true); [ -z "$l" ] && l=0
  s=$(ssh -o ConnectTimeout=20 sebastian "grep -c COMPLETE ~/grassr_calc/$S6B/PROGRESS.txt 2>/dev/null" 2>/dev/null | tr -d ' '); [ -z "$s" ] && s=0
  [ "$l" -ge 1 ] && [ "$s" -ge 1 ] && break
  sleep 300
done
log "both halves COMPLETE — merging"
rsync -a "sebastian:~/grassr_calc/$S6B/per_cell/" "$S6B/per_cell/"
n=$(ls "$S6B/per_cell" | wc -l | tr -d ' ')
log "merged: $n cells (expect 285)"
log "restoring Sebastian's lme4"
ssh sebastian 'LIB="/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"; [ -d "$LIB/lme4.stage6-parked" ] && mv "$LIB/lme4.stage6-parked" "$LIB/lme4" && echo restored'
log "re-extracting pooled nulls"
Rscript grassr/simulation/v070_program/tier1/stage6_production_null/extract_nulls_json.R || { log "EXTRACT FAILED"; exit 1; }
log "STAGE 6B COMPLETE — pooled nulls ready for artifact regen + power sweep launch"
