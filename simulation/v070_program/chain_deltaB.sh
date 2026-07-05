#!/usr/bin/env bash
# chain_deltaB.sh -- per-machine driver for the Option-B delta_hat regen.
# Runs THIS machine's stage-6B null split to completion, then its stage-7B
# power split. Per-machine only: it does NOT rsync/merge across machines
# (the operator merges per_cell dirs and runs extract/analyze afterward).
# Idempotent: a completed stage is detected by its split DONE marker and
# skipped; an interrupted stage resumes (the R scripts skip readable cells).
#
# Launch (per machine), e.g. laptop:
#   SPLIT_MOD=2 SPLIT_REM=0 WORKERS=10 \
#   CANDIDATE_LIB=/path/to/optionB/lib ./chain_deltaB.sh
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # v070_program root
S6="$SCRIPT_DIR/tier1/stage6_null_deltaB"
S7="$SCRIPT_DIR/tier1/stage7_power_deltaB"

: "${CANDIDATE_LIB:?set CANDIDATE_LIB to the Option-B candidate grassr library}"
: "${SPLIT_MOD:=1}"
: "${SPLIT_REM:=0}"
: "${WORKERS:=8}"
export CANDIDATE_LIB SPLIT_MOD SPLIT_REM WORKERS

HOST="$(hostname -s)"
LOG="$SCRIPT_DIR/chain_deltaB_${HOST}.log"
log() { echo "[$(date '+%F %T')] chain_deltaB[$HOST $SPLIT_REM/$SPLIT_MOD]: $*" | tee -a "$LOG"; }

RSCRIPT="$(command -v Rscript)"
[ -n "$RSCRIPT" ] || { log "FATAL: Rscript not on PATH"; exit 127; }

S6_DONE="$S6/split_${SPLIT_REM}_of_${SPLIT_MOD}.DONE"
S7_DONE="$S7/split_${SPLIT_REM}_of_${SPLIT_MOD}.DONE"

log "start (CANDIDATE_LIB=$CANDIDATE_LIB, WORKERS=$WORKERS)"

# ---- stage 6B null ----
if [ -f "$S6_DONE" ]; then
  log "stage6B split already DONE -- skipping"
else
  log "stage6B: launching null split"
  if "$RSCRIPT" "$S6/run_null_deltaB.R" >>"$LOG" 2>&1; then
    log "stage6B split COMPLETE"
  else
    log "FATAL: stage6B split failed (see $LOG)"; exit 1
  fi
fi

# ---- stage 7B power ----
if [ -f "$S7_DONE" ]; then
  log "stage7B split already DONE -- skipping"
else
  log "stage7B: launching power split"
  if "$RSCRIPT" "$S7/run_power_deltaB.R" >>"$LOG" 2>&1; then
    log "stage7B split COMPLETE"
  else
    log "FATAL: stage7B split failed (see $LOG)"; exit 1
  fi
fi

log "chain_deltaB DONE for split ${SPLIT_REM}/${SPLIT_MOD} (both stages). "\
"Operator step next: merge per_cell across machines, then extract_nulls_deltaB.R "\
"+ analyze_power_deltaB.R."
