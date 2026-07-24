#!/usr/bin/env bash
# chain_tier2_deltaB.sh -- per-machine driver for the Tier-2 delta-B re-run.
# Runs THIS machine's tier2_deltaB split to completion, then writes a
# chain-level DONE marker. Per-machine only: it does NOT rsync/merge across
# machines (post_compute_deltaB.sh merges per_cell and runs the analysis).
# Idempotent: a completed split is detected by the R driver's split DONE marker
# and skipped; an interrupted split resumes (the R script skips readable cells).
#
# Launch it AFTER this machine's stage-7B split has finished. chain_deltaB.sh is
# already running, so the two are not chained in one command here; launch this
# standalone once chain_deltaB.sh exits, with the SAME split assignment, e.g.:
#   SPLIT_MOD=2 SPLIT_REM=0 WORKERS=10 \
#   CANDIDATE_LIB=/path/to/optionB/lib ./chain_tier2_deltaB.sh   # laptop
#   SPLIT_MOD=2 SPLIT_REM=1 WORKERS=10 \
#   CANDIDATE_LIB=/path/to/optionB/lib ./chain_tier2_deltaB.sh   # sebastian
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # v070_program root
T2="$SCRIPT_DIR/tier2_deltaB"

: "${CANDIDATE_LIB:?set CANDIDATE_LIB to the Option-B candidate grassr library}"
: "${SPLIT_MOD:=1}"
: "${SPLIT_REM:=0}"
: "${WORKERS:=8}"
export CANDIDATE_LIB SPLIT_MOD SPLIT_REM WORKERS

HOST="$(hostname -s)"
LOG="$SCRIPT_DIR/chain_tier2_deltaB_${HOST}.log"
log() { echo "[$(date '+%F %T')] chain_tier2_deltaB[$HOST $SPLIT_REM/$SPLIT_MOD]: $*" | tee -a "$LOG"; }

RSCRIPT="$(command -v Rscript)"
[ -n "$RSCRIPT" ] || { log "FATAL: Rscript not on PATH"; exit 127; }

SPLIT_DONE="$T2/split_${SPLIT_REM}_of_${SPLIT_MOD}.DONE"    # written by the R driver
CHAIN_DONE="$T2/chain_tier2_${HOST}_${SPLIT_REM}_of_${SPLIT_MOD}.DONE"

log "start (CANDIDATE_LIB=$CANDIDATE_LIB, WORKERS=$WORKERS)"

if [ -f "$SPLIT_DONE" ]; then
  log "tier2B split already DONE -- skipping compute"
else
  log "tier2B: launching split"
  if "$RSCRIPT" "$T2/run_tier2_deltaB.R" >>"$LOG" 2>&1; then
    log "tier2B split COMPLETE"
  else
    log "FATAL: tier2B split failed (see $LOG)"; exit 1
  fi
fi

# chain-level marker (this wrapper's own DONE, distinct from the R driver's)
writeln="tier2B chain done on $HOST split ${SPLIT_REM}/${SPLIT_MOD} $(date '+%F %T')"
echo "$writeln" > "$CHAIN_DONE"
log "chain_tier2_deltaB DONE for split ${SPLIT_REM}/${SPLIT_MOD}. "\
"Operator step next: on the laptop, ./post_compute_deltaB.sh (merges per_cell "\
"across machines, extracts the null, runs the power + tier2 analyses)."
