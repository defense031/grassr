#!/usr/bin/env bash
# monitor_deltaB.sh -- print ONE status line per stage and exit. No sleep
# loop; intended to be invoked periodically by an external scheduler.
# Reports, per stage/split: cells done / split-total, and a draws/sec estimate
# derived from the first vs last timestamped PROGRESS.txt cadence lines.
#
# Honors SPLIT_MOD / SPLIT_REM (default 1/0) to compute the split's cell total.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S6="$SCRIPT_DIR/tier1/stage6_null_deltaB"
S7="$SCRIPT_DIR/tier1/stage7_power_deltaB"
MOD="${SPLIT_MOD:-1}"
REM="${SPLIT_REM:-0}"

# cells assigned to this split out of a full grid of size $1
split_total() { awk -v n="$1" -v m="$MOD" -v r="$REM" \
  'BEGIN{c=0; for(i=1;i<=n;i++) if(i%m==r) c++; print c}'; }

# epoch (macOS `date -j`) from a "YYYY-MM-DD HH:MM:SS" string; empty on failure
epoch() { date -j -f "%Y-%m-%d %H:%M:%S" "$1" "+%s" 2>/dev/null || true; }

# draws/sec across cadence lines "[ts] i/N ..." in a PROGRESS file, given
# draws-per-cell. Prints "n/a" if fewer than two parsable cadence lines.
rate() {
  local prog="$1" per="$2"
  [ -f "$prog" ] || { echo "n/a"; return; }
  local lines first last t1 t2 i1 i2 e1 e2 dc dt
  lines="$(grep -E '^\[[0-9-]+ [0-9:]+\] [0-9]+/[0-9]+' "$prog" 2>/dev/null || true)"
  [ -z "$lines" ] && { echo "n/a"; return; }
  first="$(echo "$lines" | head -1)"; last="$(echo "$lines" | tail -1)"
  t1="$(echo "$first" | sed -E 's/^\[([0-9-]+ [0-9:]+)\].*/\1/')"
  t2="$(echo "$last"  | sed -E 's/^\[([0-9-]+ [0-9:]+)\].*/\1/')"
  i1="$(echo "$first" | sed -E 's/^\[[^]]*\] ([0-9]+)\/.*/\1/')"
  i2="$(echo "$last"  | sed -E 's/^\[[^]]*\] ([0-9]+)\/.*/\1/')"
  e1="$(epoch "$t1")"; e2="$(epoch "$t2")"
  [ -z "$e1" ] || [ -z "$e2" ] && { echo "n/a"; return; }
  dc=$(( (i2 - i1) * per )); dt=$(( e2 - e1 ))
  [ "$dt" -le 0 ] && { echo "n/a"; return; }
  awk -v dc="$dc" -v dt="$dt" 'BEGIN{ printf "%.0f", dc/dt }'
}

stage_line() {
  local tag="$1" dir="$2" grid="$3" per="$4"
  local total done complete r
  total="$(split_total "$grid")"
  done=0
  [ -d "$dir/per_cell" ] && done="$(find "$dir/per_cell" -name 'cell_[0-9]*.rds' 2>/dev/null | wc -l | tr -d ' ')"
  complete="no"; [ -f "$dir/split_${REM}_of_${MOD}.DONE" ] && complete="DONE"
  r="$(rate "$dir/PROGRESS.txt" "$per")"
  printf "%s split %s/%s: %s/%s cells  ~%s draws/sec  [%s]\n" \
    "$tag" "$REM" "$MOD" "$done" "$total" "$r" "$complete"
}

stage_line "stage6B" "$S6" 440   50000
stage_line "stage7B" "$S7" 11000 2000
