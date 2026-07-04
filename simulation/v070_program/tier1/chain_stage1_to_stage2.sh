#!/bin/bash
# Waits for Sebastian's stage-1 half, merges, runs stage 2 build + verify.
# Detached via nohup; progress in chain.log (gitignored).
set -uo pipefail
cd /Users/austinsemmel/Desktop/PABAK_Investigation
S1="grassr/simulation/v070_program/tier1/stage1_surface_ndense"
log() { echo "[$(date '+%F %T')] $*"; }

log "chain: waiting for Sebastian's stage-1 COMPLETE"
while ! ssh -o ConnectTimeout=20 sebastian \
    "grep -q COMPLETE ~/grassr_calc/$S1/PROGRESS.txt 2>/dev/null"; do
  sleep 120
done
log "chain: Sebastian COMPLETE — merging per_scenario"

rsync -a "sebastian:~/grassr_calc/$S1/per_scenario/" "$S1/per_scenario/"
n=$(ls "$S1/per_scenario" | wc -l | tr -d ' ')
log "chain: merged, $n scenario files (expect 34476)"
if [ "$n" -ne 34476 ]; then log "chain: COUNT MISMATCH — halting"; exit 1; fi

log "chain: stage 2 rebuild starting"
Rscript grassr/simulation/v070_program/tier1/stage2_rebuild_surface/rebuild_surface.R \
  || { log "chain: REBUILD FAILED"; exit 1; }
log "chain: stage 2 verify starting"
Rscript grassr/simulation/v070_program/tier1/stage2_rebuild_surface/verify_candidate.R \
  || { log "chain: VERIFY SCRIPT FAILED"; exit 1; }
log "chain: stage 2 complete — see DONE or NEEDS_REVIEW + size_report.txt"
