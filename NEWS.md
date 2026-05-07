# grass 0.5.1

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

# grass 0.5.0

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
