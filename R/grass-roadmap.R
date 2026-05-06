#' The grass ecosystem roadmap
#'
#' `grass` is the binary categorical-agreement submodule of the MEADOW
#' framework. Its single headline entry point — [grass_report()] — takes
#' an `N x k` binary rating matrix and returns a `grass_card`: a four-field
#' Report Card carrying the sample summary, primary coefficient with its
#' surface-position percentile, the four-band qualitative label, and the
#' cross-coefficient asymmetry diagnostic `delta_hat` with its three-tier
#' stability flag. The full panel of coefficients, percentiles, bootstrap
#' CIs, and reference-surface artifacts ride along on the same object for
#' `summary()`, `as.data.frame()`, and `plot()` access.
#'
#' @section Framework foundation:
#'
#' A core idea animates the framework: **fixed interpretation bands
#' (Landis-Koch 1977 and its descendants) are mathematically invalid
#' across prevalences, sample sizes, and rater counts.** The Target-2
#' principle of MEADOW is *context-conditioned reporting*: the percentile
#' a coefficient lands at and the qualitative band that label maps to are
#' both conditioned on the actual study `(k, N, pi_hat)`, computed by
#' inverting a calibrated reference surface rather than read off a fixed
#' table. `grass` ships the binary submodule; FIELD is planned to ship
#' alongside Paper 3 with the same Target-2 contract for variance-component
#' reliability on continuous outcomes.
#'
#' @section MEADOW submodules:
#'
#' MEADOW is the umbrella; each submodule covers one scale type. The
#' user-facing API of each submodule is the same — a single rating
#' matrix in, a Report Card out — so user code that works on a binary
#' panel today will work on a continuous panel once FIELD ships.
#'
#' \tabular{lll}{
#'   **Submodule** \tab **Scope** \tab **Status** \cr
#'   GRASS \tab Binary categorical agreement (Cohen's kappa, PABAK, AC1,
#'     Fleiss kappa, Krippendorff alpha, observed ICC for binary). Surface
#'     positioning calibrated over `k in {2,3,5,8,15,25}` and
#'     `N in {50,200,1000}`. \tab **implemented** in grass v0.2.0 — see
#'     [grass_report()] \cr
#'   FIELD \tab Continuous variance-component reliability (Shrout-Fleiss
#'     ICC family, Lin's CCC, Bland-Altman bounds, generalisability-theory
#'     variance components). Same Target-2 surface-positioning contract as
#'     GRASS. \tab planned for v1.0.0 alongside Paper 3 — see
#'     [grass_spec_continuous()] placeholder \cr
#' }
#'
#' Earlier drafts of the roadmap referenced TURF and a separate MEADOW
#' submodule for nominal multi-rater agreement. Those are retired: GRASS
#' now covers the binary multi-rater case (Fleiss kappa, Krippendorff
#' alpha, observed ICC) directly, and the framework taxonomy collapses
#' to MEADOW = GRASS + FIELD.
#'
#' @section Stable API:
#'
#' The Target-2 contract is the same across submodules:
#'
#' * [grass_report()] — primary entry point; rating matrix in,
#'   `grass_card` out. The `grass_card` carries the four-field summary
#'   (sample, primary coefficient with surface percentile, four-band
#'   label, `delta_hat` with stability flag), the full panel of
#'   coefficients with their bootstrap CIs, and the reference-surface
#'   artifacts.
#' * [position_on_surface()] — granular access to a single coefficient's
#'   surface percentile, four-band label, and qualifier.
#' * [check_asymmetry()] — granular access to the cross-coefficient
#'   `delta_hat` and three-tier stability flag.
#' * [latent_class_fit()] — per-rater Se/Sp via Dawid-Skene EM (k >= 3)
#'   or Hui-Walter bounds (k = 2), populated automatically in the
#'   divergent branch of [grass_report()].
#' * `summary()`, `as.data.frame()`, `plot()` — layered access to the
#'   full underlying panel and the surface-position visualization.
#'
#' @section Target-2 vocabulary:
#'
#' These terms, used throughout the package documentation and printed
#' Report Card, come from the merged GRASS binary-rater-reliability
#' paper (§§3-4):
#'
#' * **context-conditioned reporting convention** — the principle that
#'   the percentile and band reported for a coefficient must condition on
#'   `(k, N, pi_hat)`, not on a fixed table.
#' * **surface-position percentile** — the empirical percentile of the
#'   observed coefficient against the calibrated reference surface at
#'   the study's `(k, N, pi_hat)`.
#' * **four-band label** — the qualitative tier (`Poor` / `Moderate` /
#'   `Strong` / `Excellent`) mapped from `q_hat` (the operating-quality
#'   projection onto the `Se = Sp` diagonal) via the partition
#'   `c(0.5, 0.625, 0.75, 0.875, 1.0)`.
#' * **delta-hat (`delta_hat`) stability flag** — the cross-coefficient
#'   percentile spread (in pp) with three tiers `aligned` / `caution` /
#'   `divergent` at NP-motivated size-alpha thresholds `c(9.25, 11.75)` (paper §3.2,
#'   App G operating characteristics). When `divergent`, the band is
#'   suppressed and per-rater Se/Sp from a latent-class fit are
#'   reported instead.
#'
#' @section What is not yet implemented:
#'
#' Constructing one of the placeholder specs ([grass_spec_continuous()],
#' [grass_spec_multirater()], [grass_spec_ordinal()]) is legal so users
#' can write code ready for FIELD or for future GRASS extensions.
#' Passing a placeholder spec to a dispatching function errors with a
#' pointer back here. The placeholders are kept so that a paper or a
#' user-side script that already names a future submodule by spec
#' continues to parse.
#'
#' @section Paper:
#'
#' The foundational paper for MEADOW and the GRASS submodule is in
#' review (Semmel 202X, *Context-Conditioned Reporting for Binary Rater
#' Reliability*). The FIELD paper will follow, citing the GRASS paper as
#' the methodological precedent. Each paper accompanies a minor or major
#' release of this package rather than a new package.
#'
#' @name grass_roadmap
#' @aliases grass-roadmap
NULL
