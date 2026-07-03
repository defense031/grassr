# cran-comments.md — grassr 0.6.1

## Submission notes

This is a new submission. It is a resubmission of the package
previously submitted as 'grass' 0.6.0, revised per CRAN review
feedback (Uwe Ligges, 2026-07-03):

* Package renamed 'grass' -> 'grassr' (the archived name 'GRASS'
  is persistent and was not released for reuse).
* Tarball size reduced from 10.7 MB to under 5 MB: the bundled
  calibration reference surfaces in R/sysdata.rda are now stored
  at 5-decimal precision (their Monte Carlo simulation error is
  on the order of 1e-2, so the discarded digits carried no
  information). No functional changes.

### "Possibly misspelled words in DESCRIPTION"

Dawid, Skene, Fleiss, and Gwet are surnames (Dawid-Skene latent-class
model; Fleiss' kappa; Gwet's AC1). PABAK is the standard acronym for
prevalence-adjusted bias-adjusted kappa (Byrt, Bishop, and Carlin,
1993). "Intraclass" and "intra" are standard reliability terminology.

## Test environments

* local: macOS (Darwin 24.6), R 4.3.1
* win-builder: R-devel and R-release

## R CMD check results

0 ERRORs, 0 WARNINGs. Remaining local NOTEs:

* "unable to verify current time" — local environment artifact.
* HTML manual tidy warnings (<table> lacks "summary" attribute,
  <script> type attributes) — produced by the R 4.3 Rd-to-HTML
  toolchain on every page, not by package markup.

## Downstream dependencies

None; this is a first release.
