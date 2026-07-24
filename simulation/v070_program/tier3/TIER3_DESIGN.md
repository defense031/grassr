# Tier 3 design — intra-rater program (DRAFT, awaiting Austin's ratification)

Status: DRAFT 2026-07-05. No APPROVED file until Austin ratifies. The
runner requires APPROVED before entering this tier (same gate as Tier 2).

## Decision 0 — scope of the writing gate (amends 2026-07-04 log entry)

The 2026-07-04 decision log recorded "intra cube and incomplete designs
join the pre-paper critical path." Verified against the manuscript
2026-07-05: the paper prints NO incomplete-design number and its
discussion limitation list does not name incomplete designs. Only the
intra axis is paper-coupled (limitation four: "intra-rater percentiles
should be treated as provisional until [the dedicated cube] ships").

PROPOSED: Tier 3 for paper-gating purposes = the intra program below.
Incomplete/unbalanced designs, grass_design(), gold-anchor API, and the
F-key compendium move to the post-paper 0.9/1.0 package roadmap. This
turns the remaining pre-writing compute from month-scale to day-scale.

## The reframe (why this is not a "cube" anymore)

The shipped calibration surface's DGP is the deterministic diagonal:
every rater has Se = Sp = q exactly, conditionally independent given
C_i (paper2/code/01_dgp.R). Under the idealized intra-rater model —
one rater, W occasions, no drift, no memory — the N x W rating matrix
has the SAME joint distribution as an inter-rater panel of k = W
raters on the diagonal: identical (Se, Sp) = (q, q) per column,
conditional independence given C_i. The two models coincide by
construction, not by approximation.

Consequence: the package's "intra percentile is approximate" disclaimer
is wrong in the flattering direction — under the reference model the
inter diagonal surface at k = W is EXACT for the intra axis. What
distinguishes intra from inter is two realism deviations the reference
model excludes:

1. **Occasion drift** — the rater's (Se, Sp) shifts between viewings
   (fatigue, learning, criterion shift). Already spec'd in
   paper2/code/31_intra_dgp.R (independent logit-scale Gaussian drift
   on Se and Sp per occasion; ratified spec (b), 2026-04-21).
2. **Within-rater memory** — the rater's calls on the SAME subject
   correlate across occasions (anchoring/recall), violating
   conditional independence. Not yet in any DGP script; mechanism
   proposed in Arm H below. Memory inflates observed intra agreement,
   so it biases the surface percentile UP — the direction that
   flatters the practitioner.

So Tier 3 ships **no new surface and no sysdata growth**. Like Tier 2,
it measures misspecification cost against the existing calibration and
publishes the bounds. delta-hat never fires on the intra axis (paper
S3: a W = 2 intra design "is already a pair of viewings, so there is
no coefficient summary to abandon"), so no intra null program exists.

## Arm F — equivalence verification (drift 0, memory 0)

Purpose: convert the equivalence argument from a derivation into a
verified claim; becomes a short proposition + simulation check in the
paper's intra appendix.

Grid: W in {2, 3, 5} (all exact surface k slices); N in {15, 50, 200,
1000}; q in {0.65, 0.85, 0.97}; F: 5 logit-normal keys, mu in
{-2.197, -0.847, 0, +0.847, +2.197} at tau2 = 1.0. 180 cells x 2,000
reps via simulate_dataset_intra(tau2_drift = 0) with the memory hook
off. Pass criterion: per-cell coefficient quantiles (PABAK, AC1;
Fleiss kappa at W >= 3) match the shipped surface quantiles within MC
error. All cells sit exactly on the calibrated grid — no snapping.

## Arm G — occasion-drift bounds

Grid: tau_drift (logit-scale SD) in {0.25, 0.5, 1.0}; W in {2, 3, 5};
N in {15, 50, 200, 1000}; q in {0.65, 0.75, 0.85, 0.92, 0.97}; F: 3
keys (mu in {-2.197, 0, +2.197}, tau2 = 1.0). 540 cells x 1,000 reps.
Per rep: agreement-family coefficients + their surface percentiles at
the matched (k = W, N, q-hat) cell. Deliverable: percentile
displacement vs the drift-0 anchor (median + p95 of |displacement|),
by (W, N, q, tau_drift) — the published drift bound.

## Arm H — within-subject dependence bounds (reframed 2026-07-05)

Austin's identifiability point (2026-07-05, adopted as the framing):
from a single rater's N x W matrix you cannot distinguish MEMORY (the
rater recalls their prior call) from IDIOSYNCRATIC RUBRIC (the rater
consistently applies their own private criterion to the same subject).
Both violate conditional independence given C_i the same way: they
raise the within-subject cross-occasion match probability above what
the reference model predicts. They are observationally equivalent, so
no mechanism-specific model is claimed and none is needed.

The knob: carryover parameterization. At occasion w > 1, with
probability rho the call on subject i repeats the occasion-(w-1) call;
with probability 1 - rho it is drawn fresh. rho = 0 recovers the
reference model exactly. This is NOT a psychological model — it is the
simplest parameterization of within-subject dependence. At W = 2 the
match probability becomes rho + (1 - rho) * m0 (m0 = reference match
probability), so for match-rate-driven coefficients ANY dependence
source maps to an effective rho: the bound is mechanism-agnostic by
construction, covering memory and rubric idiosyncrasy alike. (One
sentence of this argument goes in the paper's intra appendix.)

Grid: rho in {0.1, 0.25, 0.5}; same (W, N, q, F) frame as Arm G.
540 cells x 1,000 reps. Deliverable: percentile inflation vs the
rho = 0 anchor — the "your intra percentile may flatter you" bound,
signed (inflation expected upward). Why we keep the arm: dependence is
the deviation that inflates the percentile (drift deflates), it is the
classic reviewer objection to test-retest designs, and a published
bound is the only honest answer available given the identifiability
limit above.

### The auto-correlation check (Austin's suggestion -> Arm I diagnostic)

With a single rater the dependence is not identifiable — hence the
bound. But in a k-raters x W-occasions design (HRPU-style tensor),
it IS partially separable: between-rater agreement identifies the
C_i-driven component; a rater's within-self cross-occasion agreement
in EXCESS of the between-rater level estimates rater-specific
dependence (memory or rubric — still indistinguishable from each
other, but jointly distinguishable from subject signal). Arm I adds
this as a small package diagnostic on tensor input: report the
within-minus-between excess per rater, flag raters whose excess is
large, cite the Arm H bound for the implied percentile inflation.

## Arm I — package + paper integration (code/prose, no sims)

- Package: rewrite the intra-axis card note — the surface is exact
  under the reference model; drift/memory sensitivity context printed
  (bounds cited, not conditioned on: neither drift nor rho is
  identifiable from a single rater's W x N matrix, so no new lookup —
  consistent with the no-per-trio-lookup rule).
- Paper: intra appendix rebuilt around the equivalence proposition +
  Arm F verification + Arm G/H bound tables; S3 intra subsection
  sharpened; discussion limitation four rewritten from "provisional
  pending cube" to "exact reference, bounded deviations."

## Metric scope

Agreement family only (PABAK primary per the paper's intra path, AC1,
Fleiss kappa at W >= 3). No glmer/ICC in the arms: ICC(3,1) stays an
observed-value-only report on the intra card with no percentile claim,
matching existing scope discipline. lme4 parked on both machines for
the runs (restore after — see ops notes).

## Compute + mechanics

~1.44M small panels (N x W, W <= 5, no glmer): hours on the laptop,
comfortably under a day with Sebastian splitting. Same split/resume
mechanics as Tier 2 (SPLIT_MOD/SPLIT_REM, per-cell RDS, deterministic
seeds: master_seed 20260706 + cell_id). Tier 2 Arm D (bootstrap
coverage audit, redefined to the percentile-flag convention) can run
on the machines in parallel — it is the other pre-writing compute
leftover.

## Decisions

1. Decision 0: paper gate = intra program only; rest of Tier 3 moves
   post-paper (amends the 2026-07-04 log entry). RATIFIED 2026-07-05.
2. The reframe: no new intra surface; equivalence + bounds instead of
   a shipped cube. (The paper currently promises "a dedicated
   intra-rater calibration cube" — the writing pass replaces that
   promise with the stronger exactness claim.) RATIFIED 2026-07-05.
3. Arm H: reframed 2026-07-05 as mechanism-agnostic within-subject
   dependence (Austin's memory-vs-rubric identifiability point);
   carryover rho is the parameterization, not a psychological claim;
   within-vs-between diagnostic added to Arm I for tensor designs.
   AWAITING Austin's go on the reframed arm.
4. Tier 2 Arm D runs alongside. RATIFIED 2026-07-05 ("let's hit tier
   2 Arm D if not done yet"); redefinition under the percentile
   convention appended to TIER2_DESIGN.md.
