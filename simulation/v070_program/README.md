# grassr v0.7.0 calibration program

Design + ratified decisions: `grassr/design/v0.7_calibration_program.md`.
Waterfall: `./runner.sh` (resumable; halts at gated tiers and at
split-stage merge points).

- tier1/ — grid completeness (stages 1–5, ACTIVE)
- tier2/ — DGP realism (GATED: scripts written after Tier-1 review;
  requires APPROVED file)
- tier3/ — scope growth (GATED)

Cross-machine: Stage 1 and Stage 3 split by scenario parity —
laptop runs SPLIT_MOD=2 SPLIT_REM=1, Sebastian runs SPLIT_MOD=2
SPLIT_REM=0 in tmux session `grasscalc` (bundle at
/Users/sebastian/grassr_calc, mirroring repo-relative layout).
Merge = rsync Sebastian's per_scenario/ into the laptop's, then the
stage's verify step writes DONE.
