# cran-comments.md — grass 0.6.0

## Submission notes

This is a new submission.

### Package name: request to reuse the archived name 'grass'

The incoming check reports a case-insensitive conflict with the
archived package 'GRASS' (last released 2005-2006, versions
0.1-5 through 0.2-15). That package was an interface to the
GRASS GIS 5 system, was removed from CRAN many years ago, and
its role has long been superseded by the actively maintained
'rgrass' (previously 'rgrass7'/'spgrass6') family. The present
package is unrelated to GIS: 'grass' here is the acronym of the
statistical framework it implements (Guide for Rater Agreement
under Structural Skew), described in a methods manuscript by the
package authors. We respectfully request that the long-archived
name be released for reuse. If the CRAN team prefers to keep the
name reserved, we will resubmit under an alternative name.

### Installed size NOTE (~12 MB installed, ~10 MB tarball)

The package bundles pre-computed calibration reference surfaces
in R/sysdata.rda (9.5 MB, already xz-compressed). These surfaces
are the package's core functionality: every reported percentile
is a lookup against a Monte Carlo reference distribution built
from roughly 24 million simulated rating panels across 18,780
design scenarios. Recomputing them at load time or on demand is
not feasible (multi-day compute), and thinning the grid would
degrade the calibration the package exists to provide.

## Test environments

* local: macOS (Darwin 24.6), R 4.3.1
* R CMD check --as-cran: 0 ERRORs (other than the name-conflict
  feasibility item above), 0 WARNINGs
* win-builder / R-devel: to be run before submission

## R CMD check results

Remaining NOTEs:

* installed package size — justified above.
* "unable to verify current time" — local environment artifact.
* HTML manual tidy warnings (<table> lacks "summary" attribute,
  <script> type attributes) — produced by the R 4.3 Rd-to-HTML
  toolchain on every page, not by package markup.

## Downstream dependencies

None; this is a first release.
