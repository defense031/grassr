#!/bin/bash
# Fires after stage 7 completes: run 6c high-q top-up on both machines,
# merge, re-extract pooled nulls (6 + 6b + 6c).
set -uo pipefail
cd /Users/austinsemmel/Desktop/PABAK_Investigation
S6C="grassr/simulation/v070_program/tier1/stage6c_highq"
log() { echo "[$(date '+%F %T')] chain6c: $*"; }
log "waiting for stage 7 to finish"
while ! grep -q "STAGE 7 COMPLETE" grassr/simulation/v070_program/tier1/chain7.log 2>/dev/null; do sleep 300; done
log "stage 7 done — launching 6c on both machines"
rsync -az "$S6C/" "sebastian:~/grassr_calc/$S6C/"
ssh sebastian 'export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
LIB="/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"
[ -d "$LIB/lme4" ] && mv "$LIB/lme4" "$LIB/lme4.stage6-parked"
tmux kill-session -t grasscalc 2>/dev/null || true
tmux new-session -d -s grasscalc
tmux send-keys -t grasscalc "cd ~/grassr_calc && SPLIT_MOD=2 SPLIT_REM=0 WORKERS=8 GRASS_SIM_ROOT=\$HOME/grassr_calc Rscript grassr/simulation/v070_program/tier1/stage6c_highq/run_topup.R 2>&1 | tee stage6c_run.log" Enter'
SPLIT_MOD=2 SPLIT_REM=1 WORKERS=10 Rscript "$S6C/run_topup.R" >> "$S6C/laptop_run.log" 2>&1
log "laptop half done — waiting for sebastian"
while ! ssh -o ConnectTimeout=20 sebastian "grep -q COMPLETE ~/grassr_calc/$S6C/PROGRESS.txt 2>/dev/null"; do sleep 120; done
rsync -a "sebastian:~/grassr_calc/$S6C/per_cell/" "$S6C/per_cell/"
ssh sebastian 'LIB="/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"; [ -d "$LIB/lme4.stage6-parked" ] && mv "$LIB/lme4.stage6-parked" "$LIB/lme4" && echo restored'
log "re-extracting pooled nulls (6 + 6b + 6c)"
Rscript grassr/simulation/v070_program/tier1/stage6_production_null/extract_nulls_json.R || { log "EXTRACT FAILED"; exit 1; }
log "STAGE 6C COMPLETE — high-q ridges at 50k; regenerate artifact"
