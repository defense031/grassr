# Round the bundled reference surfaces to 5 decimals for CRAN release.
#
# The full-precision surfaces store ~15 significant digits, but their
# Monte Carlo error is on the order of 1e-2; digits beyond the 5th are
# simulation noise that xz cannot compress. Rounding to 1e-5 (three
# orders of magnitude below the MC noise floor) shrinks R/sysdata.rda
# from 9.55 MB to 2.83 MB with no detectable effect on any lookup,
# test, or reported number (v0.6.1, 2026-07-03; full suite 647/0
# against the rounded surfaces).
#
# Input:  data-raw/sysdata_fullprecision_v0.6.0.rda  (archived original)
# Output: R/sysdata.rda                              (shipped)
#
# Rebuilding the full-precision file itself is done by the sibling
# build_* scripts in this directory from the simulation outputs.

DIGITS <- 5

e <- new.env()
load("data-raw/sysdata_fullprecision_v0.6.0.rda", envir = e)

# Round only the objects actually present. sysdata contents drift across
# releases (e.g. `delta_thresholds_lookup` was retired at 0.7.x, and
# `delta_null_ecdf` was added at 0.7.0), so guard every object rather than
# assuming a fixed set. `exists_in()` returns TRUE only for a non-NULL bind.
exists_in <- function(name) !is.null(e[[name]])

if (exists_in("empirical_q_hat_surface")) {
  s <- e$empirical_q_hat_surface
  s$quantiles <- round(s$quantiles, DIGITS)
  s$mean      <- round(s$mean,      DIGITS)
  s$sd        <- round(s$sd,        DIGITS)
  e$empirical_q_hat_surface <- s
}

if (exists_in("fitted_icc_reference_curves")) {
  f <- e$fitted_icc_reference_curves
  f$curves <- round(f$curves, DIGITS)
  e$fitted_icc_reference_curves <- f
}

# The thresholds table was built with data.table but the package does not
# depend on it; strip the class so S3 dispatch can never route into
# data.table methods on machines where that namespace happens to be loaded
# (v0.6.2). Retired at 0.7.x, so skip it when absent.
if (exists_in("delta_thresholds_lookup")) {
  d <- e$delta_thresholds_lookup
  attr(d, ".internal.selfref") <- NULL
  class(d) <- "data.frame"
  e$delta_thresholds_lookup <- d
}

save(list = ls(e), envir = e, file = "R/sysdata.rda", compress = "xz")
cat(sprintf("R/sysdata.rda written: %.2f MB\n", file.size("R/sysdata.rda") / 1e6))
