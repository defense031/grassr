#!/bin/bash
# Waits for stage-4 verdict (DONE or NEEDS_REVIEW — arms A-C don't depend
# on the threshold-lookup decision), then fans arms A-C across both
# machines. Arm D (coverage audit) stays gated on the stage-4 decision.
set -uo pipefail
cd /Users/austinsemmel/Desktop/PABAK_Investigation
S4="grassr/simulation/v070_program/tier1/stage4_threshold_table"
T2="grassr/simulation/v070_program/tier2"
log() { echo "[$(date '+%F %T')] chain-t2: $*"; }

log "waiting for stage-4 verdict"
while [ ! -f "$S4/DONE" ] && [ ! -f "$S4/NEEDS_REVIEW" ]; do sleep 120; done
log "stage-4 verdict present — deploying tier 2 arms"

rsync -az "$T2/" "sebastian:~/grassr_calc/$T2/"
ssh sebastian 'export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
tmux kill-session -t grasscalc 2>/dev/null || true
tmux new-session -d -s grasscalc
tmux send-keys -t grasscalc "cd ~/grassr_calc && for a in a_item_difficulty b_correlated_errors c_asymmetry_patterns; do SPLIT_MOD=2 SPLIT_REM=0 WORKERS=8 GRASS_SIM_ROOT=\$HOME/grassr_calc Rscript grassr/simulation/v070_program/tier2/arm_\$a.R; done 2>&1 | tee tier2_run.log" Enter'
log "sebastian tier-2 arms launched"
for a in a_item_difficulty b_correlated_errors c_asymmetry_patterns; do
  SPLIT_MOD=2 SPLIT_REM=1 WORKERS=10 Rscript "$T2/arm_${a}.R" >> "$T2/laptop_tier2.log" 2>&1
done
log "laptop tier-2 arms complete — waiting for sebastian, then merging"
for a in a_item_difficulty b_correlated_errors c_asymmetry_patterns; do
  while ! ssh -o ConnectTimeout=20 sebastian "grep -q 'arm_${a} split 0/2 COMPLETE' ~/grassr_calc/$T2/arm_${a}/PROGRESS.txt 2>/dev/null"; do sleep 120; done
  rsync -a "sebastian:~/grassr_calc/$T2/arm_${a}/per_cell/" "$T2/arm_${a}/per_cell/"
  n=$(ls "$T2/arm_${a}/per_cell" | wc -l | tr -d ' ')
  log "arm_${a} merged: $n cells"
done
log "TIER 2 ARMS A-C COMPLETE AND MERGED — ready for drift analysis + Arm D"
