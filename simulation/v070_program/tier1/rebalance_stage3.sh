#!/bin/bash
# v2: laptop joins Sebastian's half when its own split completes; after
# 15 min of laptop help, sync laptop's rem=0 cells to Sebastian and
# restart his run once (resume-skip re-lists done cells at startup).
# Deterministic per-cell seeds make racing/duplication bit-safe.
set -uo pipefail
cd /Users/austinsemmel/Desktop/PABAK_Investigation
S3="grassr/simulation/v070_program/tier1/stage3_threshold_grid"
log() { echo "[$(date '+%F %T')] rebal: $*"; }

while ! grep -q "split 1/2 COMPLETE" "$S3/PROGRESS.txt" 2>/dev/null; do sleep 60; done
log "laptop split done — syncing Sebastian's cells, joining his half"
rsync -a "sebastian:~/grassr_calc/$S3/per_cell/" "$S3/per_cell/"
nohup env SPLIT_MOD=2 SPLIT_REM=0 WORKERS=10 Rscript "$S3/run_threshold_grid.R" >> "$S3/laptop_run.log" 2>&1 &
HELPER=$!
log "laptop helper running (PID $HELPER); waiting 15 min before Sebastian restart"
sleep 900
log "syncing laptop cells to Sebastian and restarting his run (resume-safe)"
rsync -a "$S3/per_cell/" "sebastian:~/grassr_calc/$S3/per_cell/"
ssh sebastian 'export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
pkill -f run_threshold_grid || true
sleep 3
tmux send-keys -t grasscalc "cd ~/grassr_calc && SPLIT_MOD=2 SPLIT_REM=0 WORKERS=8 GRASS_SIM_ROOT=\$HOME/grassr_calc Rscript grassr/simulation/v070_program/tier1/stage3_threshold_grid/run_threshold_grid.R 2>&1 | tee -a stage3_run.log" Enter'
log "Sebastian restarted with shared progress; both machines racing remainder"
wait $HELPER || true
log "laptop helper finished"
