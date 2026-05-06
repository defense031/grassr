# grass

**G**uide for **R**ater **A**greement under **S**tructural **S**kew.

Rater-reliability coefficients for binary outcomes are routinely
interpreted using fixed labels from published scales, but those labels
drift with prevalence, rater count, sample size, and latent rater
asymmetry. The same panel can land in different conventional categories
just because prevalence changes. `grass` replaces fixed cutoffs with
context-conditioned **surface-position reporting**: the rating matrix
goes in, and a four-field **Report Card** comes out — sample summary,
primary coefficient with calibrated percentile, cross-coefficient
asymmetry diagnostic, and (when the panel disagrees with itself) a
per-rater latent-class fit.

## A 30-second example

```r
library(grass)
set.seed(29)

# Five raters, 200 subjects, Se = Sp = 0.85, prevalence 0.30.
truth <- rbinom(200, 1, 0.30)
Y <- sapply(1:5, function(j) {
  ifelse(truth == 1, rbinom(200, 1, 0.85), rbinom(200, 1, 0.15))
})

grass_report(ratings = Y)
```

```
GRASS Report Card

  sample      = 5 raters, N = 200, pi_hat = 0.35, tau2_hat = 0.083
  PABAK  = 0.41  ->  95th percentile
  band        = Strong (decisive)
  delta       = 1.2 pp (aligned)

  Notes:
    - Fitted-ICC F_key picked via glmer: mu_hat=-0.877, tau2_hat=2.801
        -> F_key tau2=4.0000, mu=-0.847.
    - Fitted-ICC reference (GLMM-gap corrected) at F_key=LN_mu=-0.847_tau2=4.0000,
        k=5, N=200 (family=logit_normal, M1=0.374).
    - obs_value 0.4642 above achievable maximum (0.3381); q_hat clamped.
    - Delta-method SE undefined: dE/dq near zero at q_hat.

  See `summary(...)` for full panel and CI details.
  See `plot(...)` for a surface-position visualization.
```

The percentile is computed against a reference surface calibrated for
the study's `(k, N, pi_hat, tau2_hat)` — not against a fixed 1977
cutoff table. The `aligned` flag (`delta = 1.2 pp`) means the panel
of unclamped coefficients agrees about where the study lands, so a
single number is safe to cite. The clamp note in this example is for
ICC, whose observed value (`0.4642`) sits above the achievable maximum
of its reference surface at this `(k, N)`; ICC is therefore excluded
from `delta_hat` per the v0.2.2 clamp policy (see vignette section 7).
At the *divergent* tier the headline is replaced with the per-rater
sensitivity / specificity table from a latent-class fit (point
estimates at `k >= 3`, bounded estimates at `k = 2`).

The full walkthrough — the easy case, the killer case, layered access
via `summary()` / `as.data.frame()` / `plot()`, granular building blocks
(`position_on_surface()`, `check_asymmetry()`, `latent_class_fit()`),
and the two-rater branch — is in
[`vignette("reporting-card")`](vignettes/reporting-card.Rmd).

## Install (local dev)

```r
remotes::install_github("defense031/grass")

# from the package directory:
devtools::load_all()

# or install the built tarball:
# R CMD INSTALL grass_0.2.2.tar.gz
```

## What's new in 0.2.2

- **`delta_hat` clamp-aware.** Coefficients whose observed value
  falls outside the achievable range of their reference surface
  (most often ICC at `N > 200`, where the bundled fitted-ICC
  reference is unavailable) are now excluded from `delta_hat`
  rather than inflating it to 0 or 100 percentile points. The
  exclusion is surfaced in the printed Report Card (clamped rows
  marked `[clamped — excluded from delta]` under `divergent`),
  in the `panel` data.frame's new `clamped` column, and via a
  `clamp_note` in `$notes`. See `vignette("reporting-card")`,
  Section 7.
- **`plot_surface()` for prospective design.** New top-level
  `plot_surface(metric, pi_hat = NULL, observed = NULL, ...)`
  renders the closed-form `E[metric](M_1, q)` reference surface
  as a heatmap; supports `pi_hat`-only (vertical reference line)
  and `pi_hat + observed` (pin marker via closed-form q-hat
  inversion). Useful for "what does my reference surface look
  like at this design, before I collect data?" See
  `vignette("reporting-card")`, Section 8.
- **Intra-axis approximation note.** When `grass_report(axis =
  "intra")` is called, the Report Card surfaces a note that the
  intra-axis surface percentile uses the bundled inter-rater
  diagonal calibration as an approximation; the dedicated
  intra-axis calibration cube is queued for a follow-on data
  release.

## What's new in 0.2.1

- **`plot(card, type = "panel")` bug fix.** No longer errors with
  "Discrete values supplied to continuous scale" under ggplot2
  >= 3.5; the per-coefficient surface-percentile forest renders
  cleanly.
- **ICC fallback note at `k > 15`.** When the fitted-ICC reference
  is unavailable beyond the bundled grid maximum,
  `position_on_surface()` now surfaces an explicit user-visible
  note rather than silently clamping to a smaller-`k` reference.

## What's new in 0.2.0

- **Headline API.** New `grass_report(ratings = Y)` — single-call
  workflow returning a four-field Report Card object (`grass_card`)
  with `print` / `summary` / `as.data.frame` / `plot` methods.
- **Layered access.** Six `plot()` view types: `surface` (default),
  `panel`, `thermometer`, `intervals`, `per_rater`, `diagnostic`.
- **Latent-class fit.** New `latent_class_fit()` — point estimates
  at k >= 3, bounded estimates at k = 2, with bootstrap CIs.
- **Granular building blocks.** New `check_asymmetry(ratings = Y)`
  (cross-coefficient percentile spread, NP-motivated size-alpha thresholds 9.25 and
  11.75 pp). `position_on_surface()` now accepts `ratings = Y` directly.
- **Retired Column B.** The use-case-ladder machinery
  (`classify()`, `emr_panel()`, `grass_use_case_ladder()`) is removed
  from the headline workflow; reliability is not quality control.
- **Soft deprecations.** `grass_compute()`, `grass_reference()`,
  `grass_format_report()`, `grass_methods()`, and `grass_report_by()`
  are superseded by `grass_report(ratings = Y)`. Old calls work for
  one cycle.

See `NEWS.md` for the full changelog.

## Status

v0.2.2 (development). Binary inter-rater family fully implemented; the
intra-rater axis is supported via inter-rater approximation pending the
v0.3 calibration-cube release. See `?grass_roadmap` for planned
families (ordinal, multi-rater nominal, continuous).
