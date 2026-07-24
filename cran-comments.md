# cran-comments.md — grassr 0.7.4

## Resubmission

This resubmits grassr after the CRAN review of the 0.6.2 tarball
(Konstanze Lauseker, 2026-07-11) and the incoming pretest of 0.7.3
(2026-07-12). All requests are addressed:

* **checkRd "Lost braces" fixed.** The 0.7.3 pretest flagged
  grass_roadmap.Rd line 46 (a stray LaTeX-style thousands separator,
  "1{,}000"); now plain "1,000". This was the only new issue in that
  pretest; the remaining NOTE was the new-submission/surname check
  answered below.

* **References in the Description field.** The methods' sources are
  cited in the requested format: Byrt, Bishop, and Carlin (1993)
  <doi:10.1016/0895-4356(93)90018-V> for PABAK; Gwet (2008)
  <doi:10.1348/000711006X126600> for AC1; Fleiss (1971)
  <doi:10.1037/h0031619>; Dawid and Skene (1979) <doi:10.2307/2346806>
  and Hui and Walter (1980) <doi:10.2307/2530508> for the latent-class
  fit.
* **Acronyms explained.** PABAK is expanded in the Description text
  (prevalence-adjusted bias-adjusted kappa), as is AC1 (first-order
  agreement coefficient).
* **\dontrun removed and examples unwrapped.** The latent_class_fit()
  example runs in well under 5 seconds and is now unwrapped. The four
  \donttest wrappers elsewhere in the package dated from when those
  code paths were slower; each example now completes in under a
  second, so all are unwrapped as well. No \dontrun and no \donttest
  remain anywhere in the package; every example executes
  unconditionally during checks.

## Version change since the reviewed tarball

The reviewed tarball was 0.6.2; this submission is 0.7.4. While 0.6.2
waited in the review queue we completed a planned revision of the
package's two headline statistics and the calibration data behind
them, so the resubmission carries the current package rather than a
superseded one. The changes are documented in NEWS.md:

* The reported percentile is now the pooled percentile (position of the
  observed coefficient within the design's full achievable range,
  monotone in the observed value), replacing the matched-cohort rank.
  A 95% consistency band on panel quality is reported alongside it.
* The cross-coefficient diagnostic delta_hat is now the implied-quality
  spread of the agreement family, flagged by its percentile on a
  bundled null distribution at the matched design (new sysdata object;
  385 cells, quantiles stored at 5 significant digits).
* The two previous vignettes are consolidated into a single package
  vignette (vignette("grassr")) with two synthetic worked panels. It is
  precomputed; CRAN machines render static markdown with no
  computation.
* At two raters the diagnostic reports not-applicable rather than a
  flag. Never-released pre-0.7.1 API (the four-band label arguments,
  the legacy threshold override, and an unwired spec-constructor
  family) is removed outright; the changes are documented in NEWS.md.
* Three exported functions (grass_calibration_manifest,
  grass_contribute, grass_verify_contribution) let users run open,
  seeded calibration blocks locally and build a submission bundle for
  the project repository. No network access; writes only to a
  user-supplied directory; examples run in under a second.

The CRAN test posture is unchanged from 0.6.2: a fast deterministic
smoke subset runs on CRAN (~2 s), and the full suite runs on every
push on a five-platform CI matrix
(https://github.com/defense031/grassr/actions).

### "Possibly misspelled words in DESCRIPTION"

Dawid, Skene, Fleiss, Gwet, Byrt, Carlin, Hui, and Walter are surnames
(cited in the Description references). PABAK is the standard acronym
for the prevalence-adjusted bias-adjusted kappa, expanded in the text.
"Intraclass" and "intra" are standard reliability terminology.

## Test environments

* local: macOS (Darwin 24.6), R 4.3.1
* GitHub Actions: windows-latest R-devel and R-release, macos-latest
  R-release, ubuntu-latest R-devel and R-release — R CMD check
  --as-cran, Status: OK on all five
  (https://github.com/defense031/grassr/actions). The windows R-devel
  job currently runs without lme4 (Suggests) because of a confirmed
  upstream lme4 bug (https://github.com/lme4/lme4/issues/990: two
  heap-use-after-free defects exposed by Rcpp 1.1.2, causing
  GC-timing-dependent intermittent glmer failures; fixed in the lme4
  development version and expected on CRAN shortly). The package
  degrades gracefully without lme4 by design, and lme4-dependent
  paths are exercised on the other four platforms. We restore lme4 on
  that job as soon as the fixed lme4 reaches CRAN.

## R CMD check results

0 ERRORs, 0 WARNINGs. One NOTE:

* installed size — the bundled Monte Carlo calibration reference
  surfaces and the delta_hat null lookup in R/sysdata.rda (the
  package's core functionality) plus the precomputed vignette figures.

## Downstream dependencies

None known (grassr is not yet on CRAN; no packages import it).
