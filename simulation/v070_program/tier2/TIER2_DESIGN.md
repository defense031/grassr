# Tier 2 design — DGP realism arms (approved to execute 2026-07-04)

Approval: Austin, 2026-07-04 ("keep going on executing tier 2 — fan out
the jobs"). Discipline: every arm MEASURES misspecification cost against
the Tier-1 candidate surface; nothing recalibrates the estimand.
All arms run after Stage 3 (they share the candidate lib and the
machines). Split/resume mechanics identical to stage 3.

## Arm A — item difficulty
DGP extension: latent difficulty d_i ~ N(0, sd_d^2) per subject;
rater j's per-item accuracy becomes
  Se_ij = plogis(qlogis(Se_j) - d_i),  Sp_ij = plogis(qlogis(Sp_j) - d_i)
(hard items degrade everyone symmetrically; sd_d = 0 recovers the
calibration DGP exactly — the null anchor).
Grid: sd_d in {0, 0.5, 1.0, 1.5}; (k, N) in {(3,50), (5,200), (8,100),
(10,500), (25,1000)}; q in {0.75, 0.85, 0.92}; prev in {0.10, 0.30, 0.50};
1,000 reps. 4 x 5 x 3 x 3 = 180 cells.
Outputs per rep: 3 agreement coefficients + their candidate-surface
percentiles + delta-hat. Deliverables: (1) percentile drift vs the sd_d=0
anchor (median + tail), (2) delta-hat null-size inflation by sd_d —
i.e., does item difficulty fake divergence?

## Arm B — correlated errors
DGP extension: rater error indicators on item i are exchangeably
correlated via a Gaussian-threshold copula with correlation rho
(shared latent shock z_i; rater j errs when
Phi-mixed draw crosses its Se_j/Sp_j threshold). rho = 0 recovers the
calibration DGP.
Grid: rho in {0, 0.1, 0.25, 0.5}; same (k, N) x q x prev frame as Arm A.
180 cells x 1,000 reps.
Deliverables: percentile drift + delta-hat size inflation by rho, plus
the qualifier's behavior (does bootstrap dispersion absorb it?).

## Arm C — TPR breadth
Already half-covered: Stage 3's A-sweep runs at 5 q levels (the shipped
table had q = 0.85 only). Remaining: asymmetry PATTERNS at fixed
mean-norm A = 0.20 — (i) single bad rater (all asymmetry concentrated on
one j), (ii) half-panel split (current default), (iii) graded (linear
ramp across raters). k in {3, 5, 8, 10, 15, 25}, N in {50, 200, 1000},
prev {0.20, 0.50}, q 0.85, 1,000 reps. 3 x 6 x 3 x 2 = 108 cells.
Deliverable: divergent/caution TPR by pattern — does concentrating the
asymmetry help or hide it?

## Arm D — bootstrap-CI coverage + MC-SE audit
Runs LAST (needs Stage-4 thresholds for flag-flip accounting).
Sampled audit, not exhaustive: stratified sample of ~600 surface cells
(all k, all N bands, q in {0.75, 0.85, 0.95}, 6 F-keys spanning M1);
500 panels per cell, bootstrap_B = 200 as shipped. Deliverables:
qualifier calibration (decisive/moderate/weak vs realized percentile
dispersion), CI coverage vs nominal, MC-SE on every published threshold
(extends App D.7 to the full table).

### Arms A/B/C under v0.7.1 (2026-07-05, delta-B)
The stored per-rep files carry raw coefficients but not each rep's
realized pi-hat, which the AC1/Fleiss quality inversions need, so the
"arithmetic recompute" is replaced by a RE-RUN of the three arms
through the v0.7.1 pipeline (468 cells x 1,000 reps ~ 0.5M panels,
sub-hour) chained behind the stage 6/7 delta-B regeneration.
Deliverables restate in new units: pooled-percentile drift bounds and
delta-B flag rates. Per-draw q-hats stored this time.

### Arm D redefinition (2026-07-05 — thresholds table retired)
The Stage-4 threshold table no longer ships; the flag is delta-hat's
percentile on the matched (k, N, q-hat) null ECDF (ratified
2026-07-04). Arm D's deliverables restate as:
(i)   qualifier calibration — bootstrap p* bands vs realized surface-
      percentile dispersion (unchanged in concept);
(ii)  bootstrap CI coverage of the surface percentile vs nominal;
(iii) realized flag size — false caution/divergent rates of the
      >=95th / >=99th percentile convention on null panels drawn
      through the production pipeline at the sampled cells, vs the
      nominal 5% / 1%;
(iv)  MC-SE table for every published quantity: surface quantiles
      (from per-cell rep counts) and null-ECDF percentile positions
      (convergence program already bounds these at <= 0.22pp; Arm D
      folds them into one published table).
Frame unchanged: ~600 stratified cells, 500 panels per cell,
bootstrap_B = 200 as shipped, lme4 parked. Runs on the 0.7.0 package
(tarball build of 2026-07-05).

## Arm E — build pipeline consolidation (code, no sims)
One end-to-end `data-raw/build_sysdata.R` that regenerates every sysdata
object from per-scenario sources + rounding in one call; retires the
augment-script sprawl. Written during Arm A-D compute.
