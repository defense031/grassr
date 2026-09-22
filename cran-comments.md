# cran-comments.md — grassr 0.8.0

## Resubmission

The 2026-09-21 submission of this version was returned by the incoming
pre-test with one NOTE on r-devel-linux-x86_64-debian-gcc: the
`grass_power()` example took 5.1 s elapsed. Windows was OK. The example
ran three solves; the two secondary ones (the prevalence-range solve
and the `target =` mode) are removed from the example and pointed to
the vignette, leaving one solve and its plot, which runs in under 1 s
locally and about a third of the previous total. Nothing else changed.

## Update

This is an update to grassr 0.7.4, currently on CRAN. It is the first
update since that release.

The package's calibration data is regenerated at much finer resolution,
and the lookups that read it now interpolate rather than snapping to the
nearest calibrated cell. The user-facing functions and their arguments
are unchanged apart from one new argument (`pi_hat` on an internal
lookup's exported wrapper); no function is removed, renamed, or
deprecated in this release. Full detail is in NEWS.md:

* The bundled null distribution for the `delta_hat` diagnostic moves
  from 385 prevalence-pooled cells to 11,616 fully resolved cells
  (rater count, sample size, panel quality, prevalence), each at 50,000
  Monte Carlo draws. The pooled version's realized flag rate ran to
  roughly twice its nominal level at extreme prevalence; the resolved
  version removes that by construction.
* The null and reference-surface lookups interpolate between calibrated
  nodes. A pre-registered hold-out program (280 off-grid designs, 25,000
  draws each) scores the interpolated reading within 0.9 percentile
  points of each design's own simulated truth on average.
* Reported percentiles at off-grid designs therefore change relative to
  0.7.4, as intended. On-node readings are unchanged.

## Installed size is smaller than in 0.7.4

The size NOTE on the previous release prompted a storage-design pass.
The two large bundled arrays hold values that are exact multiples of
1e-5, so storing them as doubles spent 64 bits per value on at most 17
bits of content. They are now stored as integer-delta byte planes:

* R/sysdata.rda: 8.6 MB -> 4.9 MB
* tarball: 9.1 MB -> 5.4 MB
* installed size: about 6 Mb (was ~9 Mb)

Reconstruction is exact rather than approximate — the decoded arrays are
`identical()` to those in 0.7.4 — so no reported value changes as a
result of the storage change. The build step is in data-raw/ and the
measured design record is in design/v0.8.0_surface_encoding.md.

The remaining installed size is the Monte Carlo calibration reference
surfaces and the null lookup in R/sysdata.rda. These are the package's
core functionality: without them it cannot position a coefficient
against its design, which is what the package is for.

## Test posture

Unchanged from 0.7.4. A fast deterministic smoke subset runs on CRAN
(~2 s); the full suite (874 assertions, 0 failures) runs on every push
on a five-platform CI matrix
(https://github.com/defense031/grassr/actions). The vignette is
precomputed, so CRAN machines render static markdown with no
computation.

### "Possibly misspelled words in DESCRIPTION"

Dawid, Skene, Fleiss, Gwet, Byrt, Carlin, Hui, and Walter are surnames
(cited in the Description references). PABAK is the standard acronym
for the prevalence-adjusted bias-adjusted kappa, expanded in the text.
"Intraclass" and "intra" are standard reliability terminology.

## Test environments

* local: macOS (Darwin 24.6), R 4.3.1 — R CMD check --as-cran
  --run-donttest
* GitHub Actions: windows-latest R-devel and R-release, macos-latest
  R-release, ubuntu-latest R-devel and R-release — R CMD check
  --as-cran, Status: OK on all five. The windows R-devel job had run
  without lme4 (Suggests) since 2026-07-11 because of a confirmed
  upstream lme4 bug (https://github.com/lme4/lme4/issues/990). lme4
  2.0-6 (CRAN, 2026-07-16) fixes both defects; we verified the fix
  with 13 consecutive clean R-devel Windows check runs against lme4
  2.0.6 + Rcpp 1.1.2 and restored lme4 on that job, so all five
  platforms now check with the full Suggests set.

## R CMD check results

0 ERRORs, 0 WARNINGs. One NOTE:

* installed size — the bundled calibration reference surfaces and the
  delta_hat null lookup in R/sysdata.rda, reduced by 43% in this
  release as described above.

## Downstream dependencies

None known; no packages on CRAN import grassr.
