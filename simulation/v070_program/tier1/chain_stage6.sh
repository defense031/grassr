#!/bin/bash
# Stage-6 completion chain: wait both halves -> merge -> extract null
# summaries JSON -> restore Sebastian's parked lme4.
set -uo pipefail
cd /Users/austinsemmel/Desktop/PABAK_Investigation
S6="grassr/simulation/v070_program/tier1/stage6_production_null"
log() { echo "[$(date '+%F %T')] chain6: $*"; }
log "waiting for both stage-6 halves"
while true; do
  l=$(grep -c COMPLETE "$S6/PROGRESS.txt" 2>/dev/null || true); [ -z "$l" ] && l=0
  s=$(ssh -o ConnectTimeout=20 sebastian "grep -c COMPLETE ~/grassr_calc/$S6/PROGRESS.txt 2>/dev/null" 2>/dev/null | tr -d ' '); [ -z "$s" ] && s=0
  [ "$l" -ge 1 ] && [ "$s" -ge 1 ] && break
  sleep 600
done
log "both halves COMPLETE — merging"
rsync -a "sebastian:~/grassr_calc/$S6/per_cell/" "$S6/per_cell/"
n=$(ls "$S6/per_cell" | wc -l | tr -d ' ')
log "merged: $n cells (expect 2200)"
[ "$n" -ne 2200 ] && { log "COUNT MISMATCH — halting (lme4 NOT restored)"; exit 1; }
log "restoring Sebastian's lme4"
ssh sebastian 'LIB="/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"; [ -d "$LIB/lme4.stage6-parked" ] && mv "$LIB/lme4.stage6-parked" "$LIB/lme4" && echo restored'
log "extracting production null summaries JSON"
Rscript "$S6/extract_nulls_json.R" || { log "EXTRACTION FAILED"; exit 1; }
log "STAGE 6 COMPLETE — production nulls.json ready; regenerate artifact + build ECDF sysdata object next"
