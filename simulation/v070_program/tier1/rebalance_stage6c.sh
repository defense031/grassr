#!/bin/bash
set -uo pipefail
cd /Users/austinsemmel/Desktop/PABAK_Investigation
S6C="grassr/simulation/v070_program/tier1/stage6c_highq"
while ! grep -q "split 1/2 COMPLETE" "$S6C/PROGRESS.txt" 2>/dev/null; do sleep 60; done
echo "[$(date '+%F %T')] rebal6c: laptop done — joining sebastian's half"
rsync -a "sebastian:~/grassr_calc/$S6C/per_cell/" "$S6C/per_cell/"
nohup env SPLIT_MOD=2 SPLIT_REM=0 WORKERS=10 GRASS_SIM_ROOT=/Users/austinsemmel/Desktop/PABAK_Investigation Rscript "$S6C/run_topup.R" >> "$S6C/laptop_run.log" 2>&1 &
sleep 600
rsync -a "$S6C/per_cell/" "sebastian:~/grassr_calc/$S6C/per_cell/"
ssh sebastian 'export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
pkill -f run_topup || true; sleep 3
tmux send-keys -t grasscalc "cd ~/grassr_calc && SPLIT_MOD=2 SPLIT_REM=0 WORKERS=8 GRASS_SIM_ROOT=\$HOME/grassr_calc Rscript grassr/simulation/v070_program/tier1/stage6c_highq/run_topup.R 2>&1 | tee -a stage6c_run.log" Enter'
echo "[$(date '+%F %T')] rebal6c: sebastian restarted with shared progress"
