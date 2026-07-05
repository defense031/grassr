# grassr

**G**uide for **R**ater **A**greement under **S**tructural **S**kew.

Rater-reliability coefficients for binary outcomes are routinely
interpreted using fixed labels from published scales, but those labels
drift with prevalence, rater count, sample size, and latent rater
asymmetry. The same panel can land in different conventional categories
just because prevalence changes. `grassr` replaces fixed cutoffs with
context-conditioned **surface-position reporting**: the rating matrix
goes in, and a four-field **Report Card** comes out — sample summary,
primary coefficient with calibrated percentile, cross-coefficient
asymmetry diagnostic, and (when the panel disagrees with itself) a
per-rater latent-class fit.

## A 30-second example

```r
library(grassr)
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
  PABAK  = 0.41  ->  95th pct; consistent with panel quality [0.82, 0.90]
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
from `delta_hat` per the clamping policy (see vignette section 7).
ICC also carries a `[distribution-sensitive]` marker because its
reference surface depends on the full subject-prevalence distribution
F rather than on `(q, pi_+)` alone, and is excluded from `delta_hat`
by construction regardless of clamp status. At the *divergent* tier the headline is
replaced with the per-rater sensitivity / specificity table from a
latent-class fit (point estimates at `k >= 3`, bounded estimates at
`k = 2`).

The full walkthrough — the easy case, the killer case, layered access
via `summary()` / `as.data.frame()` / `plot()`, granular building blocks
(`position_on_surface()`, `check_asymmetry()`, `latent_class_fit()`),
and the two-rater branch — is in
[`vignette("reporting-card")`](vignettes/reporting-card.Rmd).

## Install (local dev)

```r
remotes::install_github("defense031/grassr")

# from the package directory:
devtools::load_all()

# the command above installs the repository default branch (development
# head); the CRAN submission line is the 0.6.x series.
```

## What you get in v0.5.x

- **Headline API.** `grass_report(ratings = Y)` — single-call
  workflow returning a four-field Report Card object (`grass_card`)
  with `print` / `summary` / `as.data.frame` / `plot` methods. Six
  `plot()` view types: `surface` (default), `panel`, `thermometer`,
  `intervals`, `per_rater`, `diagnostic`.
- **Surface-position primitive.** `position_on_surface(ratings = Y,
  metric = ...)` places one observed coefficient on its
  context-conditioned reference surface and returns the pooled
  percentile (the coefficient's position within the design's achievable
  agreement range) together with a 95% test-inversion consistency band
  on panel quality -- the quality levels whose sampling distributions
  are consistent with the observed value at this design. The pooled
  percentile is the categorical score; no four-way labeled band is
  interposed between the percentile and the reader. Closed-form
  references for PABAK, AC1, Fleiss kappa; bundled fitted-ICC reference
  for ICC across `k in {3, 5, 8, 15, 25} x N in {30, 50, 75, 100, 200,
  300, 500, 1000}` on a 14-point q-grid.
- **Cross-coefficient stability signal.** `check_asymmetry(ratings = Y)`
  returns the cross-coefficient implied-quality spread `delta_hat` (in
  pp of quality) computed over the three-coefficient agreement family
  (PABAK, AC1, Fleiss kappa) plus an aligned / caution / divergent flag
  set by `delta_hat`'s percentile on the matched (k, N, q_hat) null
  distribution (>= 95th caution, >= 99th divergent). ICC is reported on
  the panel with a `[distribution-sensitive]` marker but does not
  enter `delta_hat` by construction, because its reference surface
  depends on the full subject-prevalence distribution F rather than
  on `(q, pi_+)` alone.
- **Divergent-flag recovery path.** `pairwise_agreement()` returns a
  pairwise PABAK matrix on the `k = 2` reference surface plus
  per-rater pooled-reference `(Se_tilde, Sp_tilde)` against the
  panel-majority of the other `k - 1` raters. Auto-triggered from
  `grass_report()` under the divergent flag.
- **Latent-class fit.** `latent_class_fit()` returns per-rater
  Dawid-Skene `(Se_hat, Sp_hat)` point estimates at `k >= 3` and
  Hui-Walter bounds at `k = 2`, with bootstrap CIs.
- **Prospective design.** `plot_surface(metric, pi_hat, observed)`
  renders the closed-form reference surface as a heatmap with an
  optional pin marker at the observed value; useful for
  "what does my reference surface look like at this design, before
  I collect data?"

See `NEWS.md` for the release history and `?grass_roadmap`
for planned future families (ordinal, multi-rater nominal, continuous).

## Status

v0.7.1 (in development). The pooled-percentile + consistency-band card
redesign and the implied-quality `delta_hat` (see NEWS.md); the CRAN
submission line is the 0.6.2 series. Binary inter-rater and intra-rater
families fully implemented. For the intra-rater axis, the inter-rater
diagonal calibration is exact for the idealized intra model (an
equivalence proposition); drift and within-subject dependence bounds are
published rather than conditioned on. See `?grass_roadmap` for planned
families.
