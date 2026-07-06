# grassr 0.7.1

The card-redesign release. The headline reading changes from a
nearest-cell surface percentile with a four-band adjective to a pooled
percentile plus a consistency band on panel quality, and `delta_hat`
becomes an implied-quality spread. Supersedes the 0.7.0 tarball.

## Pooled percentile + consistency band (sweep redesign)

`position_on_surface()` now evaluates the observed coefficient against
every calibrated quality level at the matched design -- a sweep -- and
reports three read-outs of one object:

* **Pooled percentile** -- the trapezoid-weighted average of
  `p(q) = P(coefficient <= obs | quality q, this design)` over the
  calibrated quality axis: the coefficient's position within the
  design's achievable agreement range, monotone in the coefficient by
  construction. This retires the nearest-q_hat-cell percentile, whose
  cohort selection by a statistic derived from the coefficient made it a
  non-monotone sawtooth (7-9 full 0-100 cycles across the PABAK range at
  k=8-10, N=200; five-agent panel review 2026-07-05).
* **Consistency band** -- the 95% test-inversion band on panel quality
  `q_hat`: the quality levels whose sampling distributions are
  consistent with the observed value at this design. `q_hat` is promoted
  to the card via this band (previously internal scaffolding, never
  printed).
* **`sweep`** -- the full `p(q)` profile (the sweep-ridgeline graphic).

The stipulated four-band partition (Poor / Moderate / Strong /
Excellent) and the modal-band bootstrap qualifier (decisive / moderate /
weak), with their 0.60 / 0.90 cuts, are retired. The `bands` and
`band_labels` arguments to `position_on_surface()` and `grass_report()`
are deprecated (non-default values warn and are ignored).

## delta-hat is now an implied-quality spread (Option B)

`delta_hat` becomes the max-min spread of the agreement family's
*implied panel qualities* (PABAK, mean AC1, Fleiss kappa), in quality
percentage points: each coefficient inverts to its own `q_hat` on the
shared (q, pi_+) reference, and their spread measures cross-coefficient
model discordance in interpretable units (Option B, ratified
2026-07-05). The previous definition -- spread of surface percentiles
across four coefficients including Krippendorff alpha -- ran through the
retired sawtooth machinery. Krippendorff alpha left the panel at 0.6.0
and does not re-enter; ICC is never in `delta_hat`. The flag convention
(percentile on the matched (k, N, q_hat) null, >= 95th caution,
>= 99th divergent) is unchanged in form.

* The delta null was REGENERATED under the Option-B implied-quality
  definition (2026-07-06): 440 (k, N, q) ridges, ~22M production-pipeline
  draws, at least 46,002 per ridge. The shipped `delta_null_ecdf`
  carries the 385 k >= 3 ridges on a fine 1% + 99.5% probability grid.
  At k = 2 the two-coefficient family implies identical quality by
  construction (the null is a point mass at zero on all draws), so
  `delta_hat` reports `not_applicable` there and asymmetry assessment
  routes to the pairwise/bounds path. Tie runs (point masses) use the
  mid-p convention. Divergent-flag power at the same asymmetry roughly
  doubles relative to the retired percentile-spread definition (e.g.
  TPR 0.82 at mean-norm A = 0.20 pooled over designs, vs 0.38 at the
  modal design previously); realized null flag rates are slightly
  conservative (caution 3.7% vs nominal 5%, divergent 0.6% vs 1%).

## Bug fixes

* **Zero-plateau percentile.** `delta_null_percentile()` now resolves a
  value at or below the null's lowest stored quantile to the TOP of the
  zero/floor plateau. At small N the null often has a point mass at
  `delta_hat = 0` reaching well past the 1st percentile; 0.7.0 reported a
  value at the plateau's lowest prob, so identical raters printed
  "1.0 percentile" on cells where `P(D <= 0)` was ~35%.
* **`grass_report()` lme4 NA-skip.** A missing lme4 (Suggests) or a
  failed glmer no longer hard-errors `grass_report()`: non-finite panel
  entries (e.g. ICC when lme4 is unavailable) are dropped before
  positioning and the coefficient degrades to "absent from the panel".
  0.7.0's NEWS claimed this fix, but the 0.7.0 code landed it in
  `check_asymmetry()` only; the `grass_report()` positioning loop still
  errored. It is genuinely fixed here.

## Documentation layer

* DESCRIPTION, the roxygen sources for `check_asymmetry()`,
  `grass_report()`, `position_on_surface()`, and `?grass_roadmap`, and
  the README were still describing the 0.6.x conventions (per-(k, N)
  size-alpha thresholds `c(9.25, 11.75)`, the four-band adjective, the
  four-coefficient percentile spread). All corrected to the 0.7.1 card
  convention. Suggests-dependent examples wrapped in `\donttest{}`.

## Hygiene

* Removed macOS duplicate `" 2"` file strays from the source tree.
* `inst/preview/` and dev-scratch example scripts excluded from the
  build via `.Rbuildignore`.
* `data-raw/round_sysdata_for_release.R` skips sysdata objects that are
  absent (the retired `delta_thresholds_lookup`) instead of erroring.

# grassr 0.7.0

The calibration release. Two structural changes, both from the v0.7.0
calibration program (July 2026).

## delta-hat flag redesigned: percentile on the matched null

The per-(k, N) threshold table is retired. `check_asymmetry()` and
`grass_report()` now report delta-hat as its **percentile on the null
distribution of delta-hat at the matched (k, N, q-hat) cell**, with
flags as conventions on that percentile (>= 95th caution, >= 99th
divergent). Motivation: recalibration across panel quality showed the
old single-quality (q = 0.85) thresholds were mis-sized off their
calibration point, and at small N the null's extreme quantiles are not
stably invertible; the ECDF-position reading is stable, holds its
false-positive rate at every quality by construction, and matches how
every other quantity on the card is read.

* New sysdata object `delta_null_ecdf`: 440 (k, N, q) null ECDFs
  generated through the production pipeline (24.85M panels, >= 50k
  draws per cell), with per-cell instability metadata.
* New card fields: `delta$delta_percentile`, `delta$matched_null`;
  `delta$thresholds` now carries the implied pp cuts (95th/99th of the
  matched null) as context.
* `delta_thresholds` argument deprecated (honored with legacy
  semantics plus a warning); `delta_thresholds_lookup` removed.
* Bug fix: a missing lme4 (Suggests) no longer hard-errors
  `check_asymmetry()`; non-finite panel entries are dropped before
  positioning.

## Reference surface densified

`empirical_q_hat_surface` grows from 10,140 to 44,616 cells:
N in {15, 20, 30, 50, 75, 100, 150, 200, 300, 500, 1000} (was
{50, 200, 1000}) and a unified 13-point q-grid at every k. Lookups
that previously clamped (e.g. N = 22 -> 50) now resolve to near
cells. sysdata is 6.91 MB (storage-precision reduction to follow
before any CRAN submission of this line).

# grassr 0.6.2

CRAN resubmission release. No changes to package functionality.

## CRAN check surface reduced

The 0.6.1 incoming pretest failed on CRAN's Windows R-devel machine
with silent process deaths in the test and vignette-rebuild steps (no
diagnostics in either log; the same tarball checks clean on the same
R revision -- 2026-07-03 r90206 ucrt -- on independent Windows
infrastructure, and on the Debian pretest). 0.6.2 minimizes what CRAN
machines execute, which is recommended practice regardless:

* Tests: on CRAN, a fast deterministic smoke subset now runs (~190
  assertions: metric arithmetic, input normalization, pairwise
  agreement, print methods, and the sysdata regression anchors). The
  full ~650-test suite runs on every push via the GitHub Actions
  matrix (windows/macos/ubuntu x devel/release).
* Vignettes: both vignettes are precomputed from `.Rmd.orig` sources
  (kept in the repository, regenerated with `knitr::knit()`); CRAN
  machines render static markdown only.
* `delta_thresholds_lookup` in sysdata is now a plain data.frame; it
  previously carried a stray `data.table` class from its build script
  although the package never depended on data.table.

# grassr 0.6.1

CRAN resubmission release. No functional changes.

## Package renamed: grass -> grassr

CRAN declined to release the archived name 'GRASS' (package names are
persistent, case-insensitively). The package is now 'grassr'. All
function names, arguments, and behavior are unchanged -- the API keeps
the `grass_*` prefix, matching the GRASS framework the package
implements.

## Bundled reference surfaces rounded to 5 decimals

The reference surfaces in `R/sysdata.rda` are now stored at 5-decimal
precision (previously full 64-bit doubles), shrinking the file from
9.55 MB to 2.83 MB and the tarball to under 5 MB. The discarded digits
are Monte Carlo noise: the surfaces' simulation error is on the order
of 1e-2, three orders of magnitude above the rounding. No lookup,
threshold, test, or reported number changes at reported precision. The
full-precision surfaces and the rounding script are archived in
`data-raw/`.

# grassr 0.6.0

First CRAN release. Two behavioral changes and a legacy-API removal,
plus CRAN compliance fixes.

## Behavioral change — Krippendorff's alpha removed from the Report Card panel and delta-hat

In the binary fully-crossed designs this package targets, Krippendorff's
alpha coincides with Fleiss' kappa asymptotically (median absolute
difference 0.00024 across the 10,140-cell calibration grid; the
asymptotic argument is in the accompanying paper's Appendix A.2). The
alpha row added no information to the Report Card and its
small-sample deviation from Fleiss' kappa could tip borderline
delta-hat readings across a calibrated threshold on noise alone. As of
0.6.0:

* `compute_panel()` no longer returns `krippendorff_a`; the Report Card
  panel is PABAK / AC1 / Cohen's kappa at k = 2 and PABAK / AC1 /
  Fleiss' kappa / ICC at k >= 3.
* `delta_hat` is the surface-percentile spread over PABAK, mean AC1,
  and Fleiss' kappa. The bundled per-(k, N) thresholds hold for the
  3-coefficient spread (sensitivity check: sub-pp shift, within MC-SE).
* `obs_krippendorff_alpha()` is newly exported for users who need the
  value for cross-study comparison, and
  `position_on_surface(metric = "krippendorff_a")` still positions it
  on its closed-form reference.

## Legacy API removed

The pre-0.2.0 `grass_result` object and its consumers were deprecated
at 0.2.0 and have been unreachable from the public API since then (no
exported constructor produced a `grass_result`). Removed:
`grass_format_report()`, `grass_methods()`, `grass_report_by()`,
`print.grass_result()`, `summary.grass_result()`,
`tidy.grass_result()`, `as.data.frame.grass_result()`. The modern equivalents are
`print(grass_report(Y))`, `summary()`, `as.data.frame()`, and the
Methods-paragraph templates in the reporting-card vignette.

## CRAN compliance

* All R sources ASCII-clean; Rd math fixed in `pairwise_agreement`;
  broken Rd cross-links resolved.
* `bootstrap_delta_B` argument documented.
* ggplot2 NSE column names declared via `globalVariables()`.
* DESCRIPTION rewritten (placeholder self-citation removed).

# grassr 0.5.2

Bug-fix and sysdata-correction release. The v0.5.0 / v0.5.1 bundled
`empirical_q_hat_surface` had a build-time mislabeling bug at k = 2
that paired empirical q-hat distributions with the wrong scenarios for
roughly half of all k = 2 sysdata cells; a corrected k = 2 sub-grid
plus a high-q augmentation are bundled here. Internal logic is
unchanged; the `0 FAIL` test posture is preserved and 12 new
regression-anchor tests pin known design points.

## Behavioral change (BUG FIX) — k = 2 surface percentile lookup corrected

- **Cause.** The script that originally built the k = 2 portion of
  `empirical_q_hat_surface` (`grass/data-raw/extend_empirical_bands_k2.R`,
  pre-2026-05-10) concatenated per-rep batches via
  `per_rep <- c(per_rep, batch)` over `list.files(...)` (alphabetical
  filename order), then assigned `names(per_rep) <- s$scenario_id`
  assuming the resulting list was in the same order as
  `summaries.rds`. The two source machines (laptop / sebastian)
  wrote scenarios in randomized scenario-id order within each batch,
  so the alphabetical concat misaligned 864 of 1728 k = 2 entries.
  The empirical q-hat distribution under `scenario_id = X` was, in
  half of all k = 2 cells, actually the q-hat distribution from a
  different scenario.

- **Effect.** `position_on_surface(metric, pi_hat, k = 2L, N, ...)`
  returned percentiles that were correct for ~50% of the (F-key,
  q_true, N) cells and silently wrong for the rest. The most visible
  manifestation was the published Soares (2021) HEART-score panel:
  the v0.5.0 surface placed it at the 96th percentile; the corrected
  surface places it at the ~9th percentile (a structural inversion,
  not a small drift).

- **Fix.** The replacement sim
  (`paper1_2_merged/output/multirater_sim_k2_extension/`) writes one
  RDS per scenario keyed by scenario_id in the filename, so the
  alphabetical-concat bug class is structurally impossible. The new
  `extend_empirical_bands_k2.R` reads per-scenario files directly.
  k = 2 q-true grid is densified from the previous 6-point set
  {0.60, 0.70, 0.80, 0.85, 0.90, 0.95} to the canonical 11-point set
  {0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97, 0.99}.
  F-key set densified from 16 to the canonical 52 (matches k >= 3).
  Total k = 2 scenario count: 1728 -> 1716 (slightly lower because the
  new sim does not include the 5 asymmetric profiles per cell that
  were inappropriately mixed into the symmetric-anchor surface in
  v0.5.0; those are preserved separately at
  `paper2/simulation_output/k2_asym_sim/` for non-surface uses such as
  the §4.2 op_strong worked example and pairwise-PABAK validation).

- **Migration impact.** Anyone who computed grass surface percentiles
  on a k = 2 panel with v0.5.0 or v0.5.1 should re-run under v0.5.2.
  Aggregate downstream: the merged GRASS paper §5 Soares anchor was
  written against the buggy v0.5.0 surface and reports 96th pct;
  corrected reading is ~9th pct, with the Klein 2018 panels stable
  (50th, 51st pct unchanged) and the Hinde 2025 k = 10 panels
  unaffected (canonical k >= 3 surface; built correctly in v0.5.0).

## Robustness checks — high-q surface coverage

Three targeted sub-sims added at v0.5.2 verify the surface lookup
behaves as designed at high-quality boundary cells where per-cohort
empirical distributions are tightest and grid coverage matters
most. These are robustness checks of the framework's reach into
the high-q corners of the design space, not bug fixes.

- **k = 2 q = 0.97 (`multirater_sim_k2_q097_aug/`, 156 scenarios).**
  At k = 2 with the canonical 10-point q-grid, panels with
  closed-form-inverted q_hat in (0.96, 0.98) snap to q_true_near = 0.99,
  where the empirical distribution is centered above q_hat. Adding
  q = 0.97 to the grid lets such panels snap to a cohort whose
  empirical distribution brackets q_hat, returning a centered
  percentile.

- **k = 25 q in {0.92, 0.97} (`multirater_sim_k25_q092_097_aug/`,
  312 scenarios).** At k = 25 with N = 1000, per-cohort q-hat sd is
  ~0.003 (panels are tightly informative), so grid spacing of 0.05
  leaves cells in the gap between adjacent cohorts. Adding q = 0.92
  and q = 0.97 closes the 0.90–0.95 and 0.95–0.99 gaps.

- **k = 2 + k = 25 q = 0.94 (`multirater_sim_q094_aug/`, 312
  scenarios).** A residual 0.92–0.95 midpoint gap surfaced when
  `tab:card_summary` row 6 (k = 25, PABAK = 0.78, q_hat = 0.9416)
  still snapped to q = 0.95 after the q = {0.92, 0.97} augmentation.
  Adding q = 0.94 to k = 2 and k = 25 lets such panels snap to
  q = 0.94 cohort and return a centered percentile.

**Net q-grids after all v0.5.2 augmentations:**
- **k = 2: 12 points** {0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
  0.90, 0.94, 0.95, 0.97, 0.99} — canonical 10-point grid plus
  0.94 and 0.97.
- **k = 25: 13 points** {0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
  0.90, 0.92, 0.94, 0.95, 0.97, 0.99} — canonical 10-point grid plus
  0.92, 0.94, and 0.97.
- **k in {3, 5, 8, 15}: 10 points** (canonical) — unchanged.

The augmentation sims use a per-scenario rds layout (kill-safe,
resume-safe) and skip glmer (the bundled surface uses only
`logit_mixed_icc_oracle`, not the glmer-fitted ICC; glmer was the
slow path that initially stalled the k = 25 run before being
disabled).

**Why three sub-sims rather than one denser grid up front?** Each
boundary gap surfaced through diagnostic work on a specific paper
anchor (post-hoc figure A at k = 2; tab:card_summary row 6 at k =
25). We added grid points where the diagnostics surfaced gaps,
keeping the augmentation tight to the actual coverage need rather
than running an order-of-magnitude denser global grid. The
structural fix — bilinear interpolation across two adjacent q_true
cohorts in `lookup_empirical_q_hat()` — is designed in
`grass/design/v0.5.3_lookup_interpolation.md` and queued for the
next minor release; it would resolve the boundary-snap pathology
*without* requiring ever-denser grids, but the augmentations bring
the v0.5.2 surface to a defensible coverage level for the paper's
specific anchor cells.

## Defensive guards in sysdata builders

- **Stronger preconditions.** `extend_empirical_bands_k2.R` and the new
  `augment_empirical_bands_k2_q097.R` both assert
  `setequal(names(per_rep), as.character(s$scenario_id))` and
  `!any(duplicated(names(per_rep)))` after labeling. These would
  have caught the original mislabeling at write time.
- **Recommended for v0.5.3.** Mirror the same preconditions in
  `build_empirical_bands.R` (canonical k >= 3 builder) — currently
  safe because `multirater_sim_v3_combined_full.rds` saves per_rep
  named by scenario_id, but the silent positional-fallback at L171-173
  could mask a future regression.

## New regression tests

- **`tests/testthat/test-sysdata-regression-anchors.R`** (12 tests, all
  PASS): pins surface-percentile lookups for the merged GRASS paper's
  §4 / §5 anchors. Soares now bracketed below 20 pct (was inflated
  to 96 in v0.5.0), Klein AE/PR ~50 pct, Hinde DSM-5-TR ~10 pct,
  Hinde ICD-11 PTSD ~95 pct, §4.1 PABAK 40-55, §4.2 PABAK 30-50 +
  AC1 50-65 (the divergent split), §4.3 unclamped 30-80 (the q = 0.97
  augmentation guard). Plus three structural invariants on the k = 2
  index. Regenerating sysdata in the future without producing these
  values will fail-fast.

## Sysdata size and shape

- `empirical_q_hat_surface`: 9,360 -> 10,140 scenarios after the
  fix and the three v0.5.2 augmentations (k = 2 went 1728 -> 1872;
  k = 25 went 1560 -> 2028; k in {3, 5, 8, 15} unchanged at 1560
  each).
- File size: 9.0 MB -> 9.1 MB on disk (xz-compressed).

## Report Card render fixes (`R/format.R`)

Two render-layer changes surfaced by the 2026-05-11 paper-vs-code
audit. Neither changes computation; both affect what
`print(grass_report(Y))` shows on screen.

- **Aligned and caution flags now render all panel coefficients,
  not just the primary.** Previously, under the aligned/caution
  branch, `format.grass_card()` emitted only the primary coefficient
  line; the other agreement-family coefficients (and ICC) were
  computed and stored in `card$panel` but suppressed in the print
  output. The new render lists every coefficient with its observed
  value and surface percentile; the primary carries the qualifier
  (`decisive` / `moderate` / `weak`) inline and a `<- primary`
  marker. This matches the §3.2 dictionary table in the paper
  (`coefficient` field "lists every coefficient with its own
  percentile") and the §3.1 Quick Start example output.

- **`[distribution-sensitive]` marker now rendered by
  `format.grass_card()` on the ICC line.** Previously the marker
  was emitted only inside `print.grass_asymmetry_panel`, so the
  Report Card a user saw via `print(grass_report(Y))` carried no
  marker on ICC. The marker now renders in both the aligned/caution
  branch and the divergent branch, alongside any
  `[oracle-fallback]` / `[oracle-icc]` reference-source markers.

- **`tau2_hat` removed from the inline sample header.** The sample
  line now reads `sample = k raters, N = ..., pi_hat = ...` only;
  the per-panel τ̂² lives in `card$sample$tau2_hat` and prints in
  the Notes section when a fitted-ICC F_key is picked via glmer.

## What did NOT change

- ICC reference surface (`fitted_icc_reference_curves`) is unchanged.
- δ̂ threshold lookup (`delta_thresholds_lookup`) is unchanged.
- Closed-form `reference_binary` lookup is unchanged.
- All exported function signatures and return objects are unchanged.
- Paper §4 worked examples (k = 5 aligned + divergent) reproduce
  exactly under the corrected sysdata; the fix touches k = 2 only.

# grassr 0.5.1

Cosmetic and documentation release. No internal-logic, sysdata, or
test changes; the `0 FAIL across 194 test files` posture from v0.5.0
is unchanged.

## Behavioral change

- **Print marker rename: `[F-shape sensitive]` → `[distribution-sensitive]`.**
  The earlier label collided with the F-statistic / ANOVA reading.
  The marker references the full subject-prevalence distribution F,
  not an F-statistic; the rename matches the merged GRASS paper's
  manuscript wording. Affects the `print` method on `grass_card` and
  `grass_asymmetry_panel` objects only. ICC panel rows now print as
  `no [distribution-sensitive]` in the `in_delta_hat` column.
  Downstream code that parses print output (rare) needs an update;
  return values, sysdata, and signatures are unchanged.

- **Print format simplified: separate `band = ...` line dropped.**
  The percentile and qualifier now appear inline as
  `Nth pct (decisive)` per the manifesto-locked Report Card format.
  The internal `card$coefficient$band` field is still computed and
  remains in the object slot for callers that need it; only the
  printed rendering changed. Divergent panels print
  `panel-agg. = suppressed (divergent)` instead of `band = suppressed`.

## Vignette

- **`methods-companion.Rmd` expanded** by ~120 lines absorbing
  Tier-1-cut content from the merged GRASS paper:
  - §3 *ICC and the subject-prevalence distribution* — full
    distribution-sensitivity discussion, bundled-grid coverage,
    Koo-Li variant decision tree.
  - §5 *δ̂ thresholds and operating characteristics* — Neyman-Pearson
    threshold derivation, snap-to-cell behavior, k = 2 collapse,
    shared-vs-split-bias detection comparison.
  - §5b *Pairwise PABAK and pooled-reference recovery* — structural
    why-PABAK-at-pair-level argument, identifiability conditional,
    54-cell recovery numbers (94% direction recovery / 90-97%
    rank correlation / 98% in-majority direction / Se-tilde mean
    bias 0.012).
  The vignette is the canonical landing for practitioners wanting
  depth beyond the published paper.

## Documentation

- README.md updated with v0.5.1 wording, refreshed example output,
  and a `ref = "v0.5.1"` install pin.
- `man/*.Rd` regenerated from roxygen blocks to reflect the marker
  text change throughout function-level documentation.

# grassr 0.5.0

Three substantive landings drive this release. (1) **Cross-coefficient
divergence `delta_hat` is now computed over the agreement family only**
(PABAK, mean AC1, Fleiss kappa, Krippendorff alpha). Each has a
closed-form reference depending on `(q, pi_+)` only, so the cross-family
spread is DGP-robust at the panel level. ICC remains on the Report
Card panel rows but does NOT enter `delta_hat`: ICC's reference
surface depends on the full F-shape variance structure rather than
`(q, pi_+)`, and a panel whose true F does not match the bundled
logit-normal reference can drift by ~20 pp at small designs. An
earlier delta_hat scope including ICC produced ~80% FPR at the
diagonal null under matched DGPs and was abandoned in favor of the
agreement-family-only construction adopted here (see § "Behavioral
change" below; paper §3 ICC scope and Appendix E for the full
quantitative case). (2) The bundled
`fitted_icc_reference_curves` sysdata object is rebuilt from a unified
Julia `MixedModels.jl` simulation extending coverage to
`k in {3, 5, 8, 15, 25} x N in {30, 50, 75, 100, 200, 300, 500, 1000}`
on a 14-point q-grid spanning `[0.35, 0.99]`. Across 1,000 random
common-applied-design panels the ICC clamp rate is 7.7% (residual
clamping concentrated at small N: 14.3% at N=30; or small k: 16.3%
at k=3). (3) The `delta_thresholds_lookup` table is freshly derived
under the agreement-family `delta_hat` definition (4-coefficient
max-min) for size-alpha calibration.

## Behavioral change

- **`delta_hat` over the agreement family only.** The panel-level
  `card$panel` data frame now carries an `in_delta_hat` logical
  column (TRUE for PABAK, AC1, Fleiss kappa, Krippendorff alpha;
  FALSE for ICC). The print method shows the column inline so users
  can see which coefficients enter the cross-family spread. The
  `notes` vector for any divergent / caution panel records the v0.5
  scope decision: ICC's surface depends on the full F-shape and is
  reported separately on the Report Card, but does not contribute
  to the asymmetry signal.

- **ICC marked `[F-shape sensitive]`.** The print method labels ICC's
  panel row with the `[F-shape sensitive]` marker, signaling to the
  reader that ICC's percentile reading carries a misspecification
  cost the agreement family does not. The cost is quantified in
  Appendix E of the merged GRASS paper.

## New

- **Per-(k, N) threshold lookup** -- `grass_report()` and
  `check_asymmetry()` now default `delta_thresholds = NULL`, which
  triggers a lookup against the `delta_thresholds_lookup` table
  (`k in {2, 3, 5, 6, 8, 10}, N in {15, 20, 30, 50, 75, 100, 200, 500}`).
  When the exact `(k, N)` has both thresholds calibrated, those are
  used directly. Otherwise the lookup searches all fully-calibrated
  cells and returns the nearest one in `(k, log10(N))` space (each
  axis range-normalized). Falling back to the modal default
  `c(9.25, 11.75)` is reserved for the impossible-with-bundled-sysdata
  case of an empty calibration table; in normal operation users always
  receive thresholds calibrated *somewhere* on the grid.

  The card field `card$delta$thresholds_source` records one of
  `"calibrated_at_k_N"`, `"snapped_to_nearest_calibrated"`,
  `"default_fallback"`, or `"user_supplied"`.

- **Print-method visibility** -- the Report Card now shows a
  `thresholds = (caution, divergent) [source]` line directly under the
  `delta = ... pp (flag)` line, so the user sees which pair was applied
  and why.

- **ICC reference-path transparency** -- `position_on_surface()` now
  records which path produced each coefficient's reference surface in
  a structured `reference_used` field (one of `"closed-form"`,
  `"fitted-icc"`, `"oracle-icc-fallback"`, `"oracle-icc-explicit"`,
  `"user-supplied"`), surfaced in `card$panel$reference_used`. When
  the fitted-ICC reference is unavailable at the user's `(k, N)` and
  the function falls back to the older oracle reference, the print
  output now marks the affected coefficient with `[oracle-fallback]`
  in addition to the existing notes-vector entry. Closes a quiet
  fallback that previously only showed up in the bottom-of-card
  notes block.

## Coverage extension

- **`fitted_icc_reference_curves` rebuilt for the v0.5 grid.** The
  bundled GLMM-gap-corrected ICC reference now covers
  `k in {3, 5, 8, 15, 25} x N in {30, 50, 75, 100, 200, 300, 500, 1000}`
  on a 14-point q-grid spanning `[0.35, 0.99]` -- 2,080 (F_key, k, N)
  lookups across 52 F_keys (logit-normal grid + 4 discrete-mixture
  presets). Built from a unified Julia + `MixedModels.jl` simulation
  averaging 500-1000 fitted-ICC reps per cell. Per-cell `mean_icc`
  values are monotone-enforced along q (cummax) and linearly
  interpolated onto the package's 501-point internal q-grid.
  Builder source: `data-raw/build_fitted_icc_reference_v0.4.R`.
  Sysdata size ~8.9 MB (xz-compressed).

- **Calibrated `delta_thresholds_lookup` for `delta_hat_4`.** The
  per-(k, N) threshold table is built against the agreement-family
  `delta_hat` definition over a 1,440-cell sim grid
  (`k in {2, 3, 5, 6, 8, 10} x N in {15, ..., 500} x A in {0, ..., 0.30}
  x mu in {0.05, 0.20, 0.50, 0.80, 0.95}`, 1,000 reps per cell). Cells
  where no threshold satisfies the size-alpha constraint at
  `t in [0, 50]` pp return `NA` and the Report Card defaults to
  caution-by-default per-rater routing per the paper's Appendix D
  specification.

- **§4.1 worked example clean under v0.5 sysdata.** The §4.1 design
  point `(k=5, N=1000, mu=-1.386, tau2=1.0)` resolves via the
  fitted-ICC path with no `[oracle-fallback]` marker; the previously
  load-bearing ICC clamp at this case is closed.


# grassr 0.1.2

Round 5 — spec-dispatch architecture + reviewer 5/5 items.

## Architecture

* **Spec-dispatch**. `grass_report()` now takes `spec = grass_spec_binary()` as its primary argument. Metric families (binary today; ordinal, multirater, continuous planned) are selected by the spec constructor rather than by family-specific function names. One package, one entry point, extensible to future families without API churn. See `?grass_roadmap`.
* **Stub constructors** `grass_spec_ordinal()`, `grass_spec_multirater()`, `grass_spec_continuous()` return valid placeholder specs so users can write release-ready code; passing one to `grass_report()` errors with a pointer to `?grass_roadmap`.
* **Result objects** now carry a `$spec` slot holding the spec used for the call. All downstream verbs (`print`, `plot`, `as.data.frame`, `tidy`, `grass_methods`, `grass_format_report`) dispatch on `result$spec$family`.

## New features

* **`reference_level` argument** on `grass_spec_binary()` accepts `0.70 / 0.80 / 0.85 / 0.90`. Each band is the analytical Youden-J-optimal expected metric value at Se = Sp = reference_level, in closed form under conditional independence of the two raters given the true class. Previously only `"high"` (= 0.85) and `"medium"` (= 0.70) were available.
* **`grass_methods(result)`** — GRRAS-compliant manuscript methods paragraph pre-filled with study numbers. Supports `format = c("markdown", "latex", "plain")`. Templated so that future metric families slot in without touching the dispatcher.
* **`.parallel = TRUE`** on `grass_report_by()` — backed by `future.apply::future_lapply()` with a `progressr::progressor()` handler. Respects the user's active `future::plan()`. Optional dependencies; clean install-command error if absent.
* **`?grass_roadmap`** help page — framework taxonomy, ecosystem architecture, what is implemented today, what is planned.
* **"The grass ecosystem"** vignette — one page on the one-package-many-families design.

## Breaking / estimand changes

* **Reference curve estimand changed.** Pre-0.1.2 shipped a simulation-derived ROC classification threshold (a classifier cutpoint separating rater scenarios of different ground-truth quality bands). 0.1.2 ships the expected metric value at the Se = Sp calibration point — a cleaner analytical estimand with a closed-form interpretation. Numeric values differ at `reference_level = 0.85` vs pre-0.1.2 `reference = "high"` for most prevalences. Full Monte-Carlo regeneration with empirical confidence bands is deferred to a future release.
* **Legacy `reference = "high" / "medium" / "none"` argument on `grass_report()`** is deprecated (soft — emits a one-time message). `"high"` maps to `reference_level = 0.85`, `"medium"` to `0.70`, `"none"` to `NULL`.
* **`grass_reference(quality = ...)`** is deprecated; use `grass_reference(reference_level = ...)`.
* **Internal `thresholds` sysdata renamed** to `reference_binary` (long-form: prevalence × reference_level × metric). `grass_reference_table()` now returns the long form by default and accepts an optional `reference_level` filter. Package-internal only; user impact limited to anyone previously calling `grass:::thresholds`.

---

# grassr 0.1.1

Reviewer-driven polish (round 4).

## New in 0.1.1 (round 4)

* `response =` argument — stats-modelling alias for `rater_cols =` (wide) and `rating =` (long). Conflicting canonical + alias errors with a "disagree" message.
* `grass_format_report(ci_width = TRUE)` — appends `"; kappa CI width = 0.26, wide"` with half-width cutpoints **calibrated for kappa**: `tight < 0.10`, `moderate < 0.20`, `wide >= 0.20`. Default `FALSE` (opt-in). Kappa Wilson-logit CIs run 2-2.5x wider than component Se/Sp CIs at the same N, so these cutpoints are higher than the Se/Sp conventions from diagnostic accuracy reporting.
* `broom::tidy.grass_result()` long-form — 9 rows (3 metrics × 3 quantities: `estimate`, `reference`, `distance`). Columns: `metric`, `quantity`, `value`, `conf.low`, `conf.high`, `n`, `prevalence`. `conf.low`/`conf.high` populated only on kappa estimate. `value` (not `estimate`) is the generic column so reference/distance rows do not read as point estimates without CIs. Plots via `ggplot(aes(metric, value, colour = quantity))` without a reshape.
* `grass_report_by(data, group, ...)` — splits `data` by `group` (bare symbol or string), runs `grass_report()` per subset, returns a tidy data.frame with `.cohort` column carrying the group value. Uses `as.data.frame(compact = TRUE)` internally.
* `as.data.frame(r, compact = TRUE)` argument (already in earlier 0.1.1 draft; surface again here).
* Wide-format error now adds a soft `id_col = ` hint when exactly one column is non-binary: *"Column 'xxx' looks like an identifier; you could try `id_col = 'xxx'`."* No hint when ambiguous (all columns binary-looking).
* `reset_grass_warnings()` — clears the once-per-session message cache so loops over many studies with consistently named columns re-emit the fallthrough warning.
* `?grass_report` now has a "Sample size caveats" section documenting the two N thresholds (`N < 10` warning, `N < 30` print note).

## Breaking

* Default `format` is now `"wide"` (was `"matrix"`) across `grass_report()`, `grass_compute()`, `grass_prevalence()`. Users passing a data.frame no longer need to name the format.
* `grass_format_report()` default is now `ascii = TRUE` (was `FALSE`). Safer for Slack, Markdown, and non-UTF-8 locales. Pass `ascii = FALSE` for a Unicode `kappa` in an RStudio console.
* `validate_prevalence()` now errors on exactly `0` or `1` (degenerate cases with no reference defined). Values in `(0, 0.01)` and `(0.99, 1)` still warn and clamp.

## New

* `id_col = NULL` argument on `format = "wide"` — drop a subject-identifier column before the two-rater detection. Example: `grass_report(df, id_col = "subject_id")`.
* `as.data.frame.grass_result()` now returns 19 columns including `kappa_ref`, `kappa_distance`, `PABAK_ref`, `PABAK_distance`, `AC1_ref`, `AC1_distance`, `reference_quality`, and `regime_note`. Pass `compact = TRUE` to drop `regime_note` when binding many rows.
* `print.grass_result()` prints a "sensitive to sampling noise" caveat in the reference-comparison section when `N < 30`. (This is distinct from the `N < 10` warning on `grass_compute()`, which flags the metrics themselves as unreliable at very small N.)

## Fixes

* When no keyword rule matches, `pick_positive()` now escalates from a `message()` to a one-time `warning()` with explicit override guidance. This prevents silent prevalence inversion when rater levels are clinical terms (e.g., `"abnormal"`, `"normal"`) that fall outside the `yes` / `1` / `true` / `positive` keyword list. Flagged by Dr. M. Chen.
* Three-column wide-format error message now names both fix options (`rater_cols = c(...)` and `id_col = "..."`) and lists the actual column names present, without guessing which are raters. Flagged by J. Okafor.

## Still v0.1.0 entry-points

`grass_compute()`, `grass_report()`, `grass_reference()`, `grass_reference_table()`, `grass_prevalence()`, `grass_format_report()`, `grass_plot()`.


# grassr 0.1.0

Initial local draft.

## Public API

* `grass_compute()` — raw metric panel from 2x2 or rater data
* `grass_report()` — contextual profile: metrics + PI + BI + regime + reference
* `grass_reference()` — look up reference curves at a given prevalence
* `grass_reference_table()` — full internal reference-curve table (23 prevalence points, High and Medium calibration levels)
* `grass_prevalence()` — estimate prevalence from rater marginals
* `grass_format_report()` — one-line paper-ready summary of a result
* `grass_plot()` — plot a grass object

## Design

* `grass_report()` returns a contextual profile (metrics + PI + BI + regime + signed distance from the reference curve). Regime is one of `balanced`, `prevalence-dominated`, `bias-dominated`, `mixed`; each carries a short structural-implication note describing what the algebra of that regime forces on the metrics.
* Plotting: `plot.grass_result()` landing plot and regime scatter, `plot.grass_reference()` curves, `theme_grass()`. Inline or legend curve labels via `labels = c("auto", "inline", "legend")`.
* Vignettes: `grass-intro`, `skewed-examples`.

## Notes for early readers

* In the paper, the reference values are called *prevalence-conditioned thresholds* (Youden-J-optimal metric values from the simulation). The package uses the name `reference` throughout its public API to signal that the values serve as comparison points reported alongside PI and BI, rather than as cutoffs for a binary verdict.
