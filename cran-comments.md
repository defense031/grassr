# cran-comments.md — grassr 0.7.2

## Submission notes

This is the first update to grassr (0.6.2 was accepted as a new
submission; 0.7.x was not submitted before this tarball). It revises
the package's two headline statistics, ships the regenerated
calibration data behind them, and adds a contributor pipeline for
extending the calibration:

* The reported percentile is now the pooled percentile (position of the
  observed coefficient within the design's full achievable range,
  monotone in the observed value), replacing the matched-cohort rank.
  A 95% consistency band on panel quality is reported alongside it.
* The cross-coefficient diagnostic delta_hat is now the implied-quality
  spread of the agreement family, flagged by its percentile on a
  bundled null distribution at the matched design (new sysdata object;
  385 cells, quantiles stored at 5 significant digits).
* The two previous vignettes are consolidated into a single package
  vignette (vignette("grassr")) that adds two synthetic worked panels.
  It is precomputed; CRAN machines render static markdown with no
  computation.
* At two raters the diagnostic reports not-applicable rather than a
  flag. Never-released pre-0.7.1 API (the four-band label arguments,
  the legacy threshold override, and an unwired spec-constructor
  family) is removed outright; the changes are documented in NEWS.md.
* New in 0.7.2: three exported functions (grass_calibration_manifest,
  grass_contribute, grass_verify_contribution) let users run open,
  seeded calibration blocks locally and build a submission bundle for
  the project repository. No network access; writes only to a
  user-supplied directory; examples are \donttest and run in seconds.

The CRAN test posture is unchanged from 0.6.2: a fast deterministic
smoke subset runs on CRAN (~2 s), and the full suite (703 tests, 0
failures) runs on every push on a five-platform CI matrix
(https://github.com/defense031/grassr/actions).

### "Possibly misspelled words in DESCRIPTION"

Dawid, Skene, Fleiss, and Gwet are surnames (Dawid-Skene latent-class
model; Fleiss' kappa; Gwet's AC1). PABAK is the standard acronym for
prevalence-adjusted bias-adjusted kappa (Byrt, Bishop, and Carlin,
1993). "Intraclass" and "intra" are standard reliability terminology.

## Test environments

* local: macOS (Darwin 24.6), R 4.3.1
* GitHub Actions: windows-latest R-devel and R-release, macos-latest
  R-release, ubuntu-latest R-devel and R-release — R CMD check
  --as-cran, Status: OK on all five
  (https://github.com/defense031/grassr/actions).

## R CMD check results

0 ERRORs, 0 WARNINGs. One NOTE:

* installed size — the bundled Monte Carlo calibration reference
  surfaces and the delta_hat null lookup in R/sysdata.rda (the
  package's core functionality) plus the precomputed vignette figures.

## Downstream dependencies

None known (first update; no packages import grassr).
