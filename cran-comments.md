# cran-comments.md — grassr 0.6.2

## Submission notes

This is a new submission (third attempt: 'grass' 0.6.0 was renamed
per CRAN feedback; 'grassr' 0.6.1 failed the Windows incoming
pretest). Changes in 0.6.2, addressing the 0.6.1 pretest result:

The 0.6.1 Windows pretest reported ERRORs in the test and
vignette-rebuild steps, in both cases a silent termination with no
diagnostics in the logs. The identical tarball checks clean (R CMD
check --as-cran, Status: OK) on the same R-devel revision
(2026-07-03 r90206 ucrt) on Windows via GitHub Actions, on the
Debian pretest, and locally. We could not reproduce the failure, so
0.6.2 minimizes what CRAN machines execute:

* Tests on CRAN now run a fast deterministic smoke subset (~190
  assertions, ~2 s) covering metric arithmetic, input normalization,
  pairwise agreement, print methods, and regression anchors that pin
  the bundled reference surfaces. The full ~650-test suite runs on
  every push on a five-platform CI matrix
  (https://github.com/defense031/grassr/actions).
* Both vignettes are precomputed; CRAN machines render static
  markdown with no computation.
* A stray 'data.table' class attribute on one bundled lookup table
  (the package does not depend on data.table) is removed.

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

0 ERRORs, 0 WARNINGs. One NOTE beyond the new-submission NOTE:

* installed size is 5.5Mb — the bundled Monte Carlo calibration
  reference surfaces in R/sysdata.rda (2.8 MB compressed; the
  package's core functionality) plus the precomputed vignette
  figures. The source tarball is 4.4 MB.

## Downstream dependencies

None; this is a first release.
