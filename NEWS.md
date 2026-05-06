# grass 0.5.0

Three substantive landings drive this release. (1) **Cross-coefficient
divergence `delta_hat` is now computed over the agreement family only**
(PABAK, mean AC1, Fleiss kappa, Krippendorff alpha). Each has a
closed-form reference depending on `(q, pi_+)` only, so the cross-family
spread is DGP-robust at the panel level. ICC remains on the Report
Card panel rows but does NOT enter `delta_hat`: ICC's reference
surface depends on the full F-shape variance structure rather than
`(q, pi_+)`, and a panel whose true F does not match the bundled
logit-normal reference can drift by ~20 pp at small designs. The
v0.4.0-pre tag (which shipped a 5-coefficient `delta_hat`) was
withdrawn after a sanity probe showed FPR ~80% at the diagonal null
under matched DGPs (see § "Behavioral change" below; paper §3 ICC
scope and Appendix E for the full quantitative case). (2) The bundled
`fitted_icc_reference_curves` sysdata object is rebuilt from a unified
Julia `MixedModels.jl` simulation extending coverage to
`k in {3, 5, 8, 15, 25} x N in {30, 50, 75, 100, 200, 300, 500, 1000}`
on a 14-point q-grid spanning `[0.35, 0.99]`. Across 1,000 random
common-applied-design panels the ICC clamp rate drops from a 98.2%
v0.3 baseline to 7.7%; residual clamping is concentrated at small N
(14.3% at N=30) or small k (16.3% at k=3). (3) The
`delta_thresholds_lookup` table is re-derived under the new
`delta_hat` definition (4-coefficient max-min over the agreement
family) for size-alpha calibration; the new table replaces the v0.3
calibration which assumed an effectively-4-coefficient `delta_hat`
(via v0.2.2's silent ICC-clamp exclusion).

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
  averaging 500-1000 fitted-ICC reps per cell. Replaces the v0.3
  bundle's narrower coverage (4 k x 2 N). Per-cell `mean_icc` values
  are monotone-enforced along q (cummax) and linearly interpolated
  onto the package's 501-point internal q-grid. Builder source:
  `data-raw/build_fitted_icc_reference_v0.4.R`. Sysdata size grows
  from ~4.5 MB to ~8.9 MB (xz-compressed).

- **Re-calibrated `delta_thresholds_lookup` for `delta_hat_4`.** The
  per-(k, N) threshold table is freshly built against the v0.5
  4-coefficient `delta_hat` definition over a 1,440-cell sim grid
  (`k in {2, 3, 5, 6, 8, 10} x N in {15, ..., 500} x A in {0, ..., 0.30}
  x mu in {0.05, 0.20, 0.50, 0.80, 0.95}`, 1,000 reps per cell). The
  new lookup replaces the v0.3 table without API change; cells where
  no threshold satisfies the size-alpha constraint at `t in [0, 50]`
  pp return `NA` and the Report Card defaults to caution-by-default
  per-rater routing per the paper's Appendix D specification.

- **§4.1 worked example clean under v0.5 sysdata.** The §4.1 design
  point `(k=5, N=1000, mu=-1.386, tau2=1.0)` resolves via the
  fitted-ICC path with no `[oracle-fallback]` marker; the previously
  load-bearing ICC clamp at this case is closed.

# grass 0.4.0-pre (WITHDRAWN 2026-05-05)

The 0.4.0 tag was withdrawn after sanity probing revealed the
5-coefficient `delta_hat` had ~80% FPR at the diagonal null even under
matched logit-normal DGPs. Root cause: ICC's reference surface depends
on the full F-shape variance structure, while the agreement family's
references depend only on `(q, pi_+)`; under any practical
finite-sample condition glmer's tau-hat for individual panels drifts
enough that ICC's surface percentile lands ~20 pp away from the
agreement family even when the population truth matches the reference.
v0.3 hid this via 98% ICC clamping (which silently excluded ICC from
`delta_hat` per the v0.2.2 clamp policy); v0.4 closed the clamp gap
but the threshold table was still calibrated to v0.3's
effectively-4-coefficient behavior.

The fix in v0.5.0 is structural: `delta_hat` is now defined over the
agreement family by construction, and ICC stays on the Report Card as
a separate F-shape-sensitive reading. The two infrastructure landings
the v0.4.0 tag was supposed to ship (per-(k, N) threshold table +
extended fitted-ICC reference) are preserved and ship in v0.5.0; only
the `delta_hat` definition and the threshold-calibration were
corrected. See the v0.5.0 entry above for the full bullet list.

# grass 0.2.3

New exported `pairwise_agreement()` for the divergent-flag recovery
path, auto-triggered from `grass_report()`; bundled
`empirical_q_hat_surface` extended with `k = 2`; per-subject
between-cluster variance `tau2_hat` surfaced on the Report Card; and
a paper- and package-wide rename of the prior "Neyman-Pearson
optimum" thresholds to "Neyman-Pearson-motivated size-alpha
thresholds" with precision about the test class.

## New

- **`pairwise_agreement(ratings, axis = c("inter", "intra"))`** —
  exported. When the panel is flagged `divergent`, the panel-aggregate
  coefficients no longer summarize the panel adequately and the
  framework's prescription (paper §3.3) is to report (a) a `k x k`
  pairwise PABAK matrix with each entry placed on the `k = 2`
  reference surface at the pair's observed marginal, and (b) per-rater
  pooled-reference sensitivity and specificity `(Se_tilde,
  Sp_tilde)` against the panel-majority of the *other* `k - 1`
  raters. The pooled-reference recovery is identifiability-clean (no
  label-switching, since the "reference" is the observable panel
  majority, not a latent class). Returns a `grass_pairwise` S3 object
  with `pabak_matrix`, `percentile_matrix`, `marginal_matrix`,
  `band_matrix`, `qualifier_matrix`, `pooled_per_rater`, `sample`,
  `notes`, and `call`. Print method shows the pairwise PABAK matrix
  (lower triangle) + surface percentiles (upper triangle) + the
  pooled-reference per-rater table.

- **`grass_report()` auto-trigger.** When the asymmetry flag is
  `divergent`, `pairwise_agreement(ratings = Y, axis = axis)` is now
  computed and stored at `grass_card$pairwise`; the divergent-branch
  block of the printed Report Card reports the pairwise matrix and
  per-rater pooled-reference table alongside the (still-reported)
  `latent_class_fit()`. Under `caution`, the Report Card emits a
  one-line note suggesting the user run `pairwise_agreement()`
  manually. Under `aligned`, no pairwise computation is triggered.

- **`plot(card, type = "pairwise")`** — new view on
  [plot.grass_card()]. `k x k` tile heatmap of the pairwise surface
  percentile matrix with the percentile (pp) printed in each cell and
  `PABAK_ij` in parentheses below; viridis-coded fill on `[0, 100]`.
  Errors when `card$pairwise` is `NULL` (i.e., when the panel is not
  divergent). The `type = "diagnostic"` composite now substitutes the
  pairwise view for the `per_rater` view when the card is divergent;
  this surfaces panel structure (clustering / outlier raters) directly
  in the at-a-glance composite.

- **Sample-line `tau2_hat`.** A method-of-moments estimator for the
  between-subject variance of positive rates (`compute_tau2_hat()`,
  internal) is now surfaced as `grass_card$sample$tau2_hat` and
  printed in the Report Card sample line. Previously this slot was
  always `NA_real_`; the four-field summary `(pi_hat, tau2_hat, k, N)`
  the paper documents is now backed by an actual estimator.

## Coverage extension

- **`empirical_q_hat_surface` extended to `k = 2`.** Built via
  `data-raw/extend_empirical_bands_k2.R`, which inverts each per-rep
  observed metric value through the closed-form `k = 2` reference
  curves and summarizes the per-cell `q_hat` distribution at 13
  quantile probabilities. Source: `paper2/simulation_output/
  k2_asym_sim/per_rep` (36 batches, 1,728 scenarios x 2,000 reps x 6
  metric columns). The bundled lookup now covers
  `k in {2, 3, 5, 8, 15, 25}`, **9,528 scenarios** total (was 7,800
  at v0.2.2). The v0.2.0 `"k=2 clamped to nearest sim-grid k=3"`
  approximation note is retired for paper-2 `k = 2` entries; pairwise
  PABAK percentiles in `pairwise_agreement()` now read directly from
  the `k = 2` surface. `sysdata.rda` grows from 4.0 MB to 4.6 MB.

## Renamed (vocabulary; behavior unchanged)

- **"Neyman-Pearson optimum" -> "Neyman-Pearson-motivated size-alpha
  thresholds"** across the package (DESCRIPTION, README, `?check_asymmetry`,
  `?grass_report`, `?grass_roadmap`, `vignette("reporting-card")`).
  The construction *is* NP-motivated — size-controlled threshold
  selection, with the smallest threshold `t` such that
  `P(delta_hat > t | symmetric panel) <= alpha` chosen by calibration on
  the simulation grid and the TPR read at the asymmetric alternative.
  The new wording is precise: NP-optimal *within the test class of
  threshold rules on `delta_hat`*, not derived from a likelihood-ratio
  statistic on the joint distribution of the panel coefficients
  (`delta_hat` is a paired-margin difference). The numerical
  thresholds (9.25, 11.75 pp) and the `aligned` / `caution` /
  `divergent` tier vocabulary are unchanged.

## Documentation

- `?pairwise_agreement` and a new vignette section walk the divergent
  Case-2 worked example end to end (paper §4.2 op_strong DGP, seed=17).
- `?grass_report` `@return` documents `surface_percentile in [0, 100]`
  in the Report Card; underlying `position_on_surface()` returns
  `percentile in [0, 1]` and is multiplied by 100 by the headline.
- `vignette("reporting-card")`: Section 3 walks the §4.2 worked
  example (op_strong, pi=0.50) and, via the new auto-trigger, walks
  `pairwise_agreement()` output directly. The `latent_class_fit()`
  section is reframed as the alongside alternative ("primary value
  when external orientation is available").

## Tests

- `tests/testthat/test-pairwise_agreement.R` — new test file with
  seven `test_that` blocks covering input validation (`k < 2`),
  output structure (S3 class, documented fields, sample list keys),
  matrix shape and symmetry (`k x k`, diagonals, percentiles in
  `[0, 100]`), the pooled-reference per-rater data frame at `k >= 3`,
  the `k = 2` omit-pool branch with the explanatory note, and the
  print method headers.
- `tests/testthat/test-grass_report-card.R` — added two auto-trigger
  tests: the divergent panel populates `card$pairwise` with a
  `grass_pairwise` object; non-divergent panels leave it `NULL`.
- `tests/testthat/test-plot-card.R` — added two tests for the new
  `type = "pairwise"` view (errors under non-divergent; returns a
  ggplot under divergent) and extended the `.mock_card()` fixture to
  include a `pairwise` slot on divergent flag.
- Final tally: **624 PASS / 0 FAIL / 28 SKIP** (553 -> 624 across
  v0.2.2 -> v0.2.3, +71 expectations from `tau2_hat` surfacing,
  pairwise coverage, and plot-view coverage).

# grass 0.2.2

Behavioral fix to `delta_hat` calculation in `check_asymmetry()` and
`grass_report()` plus two additive items: a standalone
`plot_surface()` for prospective design and an intra-rater axis
approximation note. **The behavioral fix changes the divergent-flag
behavior on panels where a coefficient's observed value falls outside
its reference-surface envelope** (most commonly ICC at `N > 200`,
where the bundled fitted-ICC reference is unavailable and the oracle
reference over-predicts).

## New

- **`plot_surface(metric, pi_hat = NULL, observed = NULL, k = NULL,
  N = NULL, ...)`**: standalone heatmap of the closed-form expectation
  `E[metric](M_1, q)` over the reference surface, intended for
  prospective study design. Uses the same closed-form machinery as
  `plot.grass_card(type = "surface")`, but driven by scalar inputs
  rather than a fitted card. When `pi_hat` is supplied, draws a
  dashed vertical reference at the design's marginal; when `observed`
  is supplied with `pi_hat`, pins a marker at `(pi_hat, hat q)`
  with `hat q` recovered by closed-form inversion. ICC surfaces are
  not supported via this function (they require fitted reference
  curves; use `plot.grass_card(type = "surface")` for ICC).

## Behavior change

- `delta_hat` is now computed only over coefficients whose surface
  percentile is *not* clamped to the envelope boundary (`0` or `100`).
  If at least two unclamped coefficients remain, they alone determine
  the spread; if fewer than two remain, the function falls back to the
  raw spread including clamped values and surfaces a note. Previously,
  ICC clamping at boundary 0 or 100 could inflate `delta_hat` by 30+ pp
  and fire a false `divergent` flag on panels where the four agreement
  coefficients (PABAK, AC1, Fleiss kappa, Krippendorff alpha) actually
  agree on percentile.
- The `panel` data.frame on the return object now includes a
  `clamped` logical column flagging which coefficients were excluded
  from `delta_hat`.
- A `clamp_note` is added to the `notes` vector whenever clamping
  affects `delta_hat`, identifying the affected coefficients.
- The default `print.grass_card()` output for `divergent` panels now
  marks clamped coefficient lines with `[clamped — excluded from delta]`
  for transparency at output time.

## Why it matters

The previous behavior fired `divergent` flag artifacts on routine
panels where ICC reference coverage was missing. Several test fixtures
that asserted `flag == "divergent"` for the §4 paper case were relying
on this artifact (the four agreement coefficients were aligned within
4 pp; ICC's clamp to percentile 0 inflated `delta_hat` to ~32 pp).
Those fixtures have been updated to use a genuinely divergent DGP
(`op_strong` profile from `paper2/code/21_asym_grid.R` — alternating
Se-favoring and Sp-favoring raters, where AC1's surface position
genuinely separates from the kappa-family).

## Notes

- **Intra-rater axis approximation surfaced.** When
  `grass_report(axis = "intra")` is called, the Report Card now
  carries a note that the intra-axis surface percentile uses the
  bundled inter-rater diagonal calibration as an approximation; the
  dedicated intra-axis calibration cube (the 1,500-rep extension
  from `paper2/code/31_intra_dgp.R`) is queued for a follow-on data
  release. Behavior is otherwise unchanged from v0.2.1; the note
  surfaces a previously silent approximation.

## Documentation

- `?check_asymmetry` has a new "Surface-envelope clamp" section
  explaining when and why clamping fires and how `delta_hat` excludes
  clamped coefficients.
- `?plot_surface` documents the new standalone surface-plot entry
  point; its `\examples` block shows the bare-surface, design-context,
  and pinned-observation use cases.
- The `vignette("reporting-card")` walks an ICC-clamp example so
  practitioners can recognize the surfaced note (Section 7), and a
  separate section demonstrates `plot_surface()` for prospective
  study design.

# grass 0.2.1

Two patches flagged at the v0.2.0 release; both are non-breaking and
unblock the merged-paper line edits.

## Fixes

- `plot(card, type = "panel")` no longer errors with "Discrete values
  supplied to continuous scale" under ggplot2 (>= 3.5). The forest of
  per-coefficient surface percentiles now renders cleanly. The fix
  switches the y-axis from a factor mapping to integer positions with
  `scale_y_continuous(breaks, labels)` so the band-name annotations
  above the top row no longer collide with a discrete scale.
- ICC at `k > 15` (outside the `fitted_icc_reference_curves` bundle)
  now falls back to the oracle ICC reference and surfaces a prominent
  user-visible note that names the gap, e.g. `"Fitted-ICC reference
  unavailable at k=25 (fitted-ICC k_grid maxes at 15); using oracle
  ICC reference (GLMM-gap not corrected). Treat the surface position
  as an approximation."` Previously the function silently clamped to
  the fitted reference at k = 15, which over-corrects the GLMM gap at
  large k. The N > 200 fallback is unchanged.

## Tests

- New regression test in `test-position_on_surface.R` asserts the
  k = 25 fallback produces the explicit-gap note and does *not* touch
  the fitted-ICC reference at k = 15.

# grass 0.2.0

Headline API rewrite for paper alignment. The Report Card is now a
single-call workflow on the rating matrix.

## New

- New headline `grass_report(ratings = Y)` — single-call workflow that
  returns a four-field Report Card (sample summary, primary coefficient
  with surface percentile, cross-coefficient asymmetry diagnostic, and
  divergent-branch per-rater latent-class fit).
- New S3 methods on the `grass_card` class: `print`, `summary`,
  `as.data.frame`, `plot` (with six view types: `surface`, `panel`,
  `thermometer`, `intervals`, `per_rater`, `diagnostic`).
- New `latent_class_fit()` — Dawid-Skene EM at k >= 3 plus Hui-Walter
  bounds at k = 2, with nonparametric bootstrap CIs on per-rater Se / Sp.
- New `check_asymmetry(ratings, ...)` — cross-coefficient surface-percentile
  spread `delta_hat` with NP-motivated size-alpha thresholds (9.25, 11.75) pp and
  tiers `aligned` / `caution` / `divergent`.
- `position_on_surface()` now accepts `ratings = Y` directly; auto-derives
  `obs_value`, `pi_hat`, `k`, `N` from the rating matrix.

## Renamed

- `check_asymmetry(se, sp, ...)` -> `check_rater_asymmetry(se, sp, ...)`.
  The old name now serves the new ratings-input flow. Old calls
  soft-deprecate to the renamed function for one cycle.

## Retired

- `classify()`, `emr_panel()`, `grass_use_case_ladder()` — Column B and
  the use-case-ladder machinery are removed from the headline workflow
  (2026-04-22 framework decision: reliability is not quality control).
  Functions still callable for one cycle with `.Deprecated()` warnings;
  will be removed entirely in v0.3.0.

## Deprecated

- `grass_compute()`, `grass_reference()`, `grass_format_report()`,
  `grass_methods()`, `grass_report_by()` — superseded by
  `grass_report(ratings = Y)`. Soft deprecation; old user code continues
  to work for one cycle.

## Documentation

- New vignette: `reporting-card` (replaces the old `reporting-card.Rmd`
  content).
- Other vignettes (`grass-intro`, `grass-ecosystem`, `skewed-examples`)
  carry the old framework and will be retired in a future release. Each
  carries a banner pointing readers to `vignette("reporting-card")`.

## Coverage notes

Sysdata audit (2026-05-02) against design doc §6.4 expectations:

- `empirical_q_hat_surface` covers k in {3, 5, 8, 15, 25}, N in
  {50, 200, 1000}, all five metrics (`pabak`, `mean_ac1`, `fleiss_kappa`,
  `krippendorff_a`, `icc`), 520 scenarios per (k, N) cell.
  **Gap: k = 2 not present.** `grass_report(ratings = Y)` at k = 2 emits
  a `"k=2 clamped to nearest sim-grid k=3"` note; surface percentiles at
  k = 2 are interpolated from the k = 3 surface and should be read as
  approximate. A k = 2 simulation slice is queued for a follow-on data
  release.
- `fitted_icc_reference_curves` covers k in {3, 5, 8, 15}, N in {50, 200},
  52 F-keys per cell. **Gap: k = 25 not present** (design doc §6.4 listed
  k in {3, 5, 8, 15, 25}). At k = 25 the package falls back to the
  oracle ICC reference curve; the GLMM-gap correction for fitted ICC is
  unavailable. The fallback is silent in the printed Report Card; users
  who need fitted-ICC interpretation at k = 25 should consult
  `position_on_surface()$notes`.
- `icc_reference_curves` (oracle) covers all 52 F-keys at the closed-form
  resolution (501-point q grid). No gaps.

---

# grass 0.1.2

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

# grass 0.1.1

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


# grass 0.1.0

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
