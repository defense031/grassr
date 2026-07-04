# =====================================================================
#  grassr playground
#
#  A hands-on tour of the package. Edit, rerun, and break things.
#  Every section is independent — you can skip around.
# =====================================================================

# ---------------------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------------------
# Pick ONE of the two blocks below. Re-run the block every time the
# package source has changed — otherwise you'll be using a stale version.
#
#  (a) Dev mode, no install. Fastest iteration.
#      Run this from the package root (the directory that holds DESCRIPTION):
#
#          devtools::load_all()
#
#      load_all() picks up every edit to R/ and tests/; no reinstall needed.
#
#  (b) Installed-package mode. Needed if you're running from RStudio's
#      console with the package already attached via library():
#
#          detach("package:grassr", unload = TRUE)           # if already loaded
#          devtools::install(quick = TRUE, upgrade = "never")
#          library(grassr)
#
# If you get "could not find function grass_report" after a rename, you're
# running stale code — redo one of the two blocks above.

library(grassr)
library(ggplot2)


# ---------------------------------------------------------------------
# 1. A 2x2 count matrix
# ---------------------------------------------------------------------
# Cell layout:   R1 rows, R2 columns. matrix(c(a, c, b, d), nrow = 2).
tab <- matrix(c(88, 10, 14, 88), nrow = 2,
              dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
tab

# Raw metric panel.
grass_compute(tab, format = "matrix")

# Full contextual report: metrics + prevalence + bias + regime + reference.
grass_report(tab, format = "matrix")


# ---------------------------------------------------------------------
# 2. Same data, four input shapes — identical results
# ---------------------------------------------------------------------
r1 <- c(rep(0, 88 + 14), rep(1, 10 + 88))
r2 <- c(rep(0, 88), rep(1, 14), rep(0, 10), rep(1, 88))

grass_compute(list(r1, r2), format = "paired")

grass_compute(data.frame(rater1 = r1, rater2 = r2), format = "wide")

long_df <- data.frame(
  subject = rep(seq_along(r1), 2),
  rater   = rep(c("rater1", "rater2"), each = length(r1)),
  rating  = c(r1, r2)
)
grass_compute(long_df, format = "long")


# ---------------------------------------------------------------------
# 3. Low-prevalence example (Feinstein & Cicchetti 1990)
# ---------------------------------------------------------------------
# Observed 80% agreement, negative kappa.
# PI is high, BI is zero — prevalence-dominated regime.
tab_skewed <- matrix(c(40, 5, 5, 0), nrow = 2,
                      dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
r_skewed <- grass_report(tab_skewed, format = "matrix")
r_skewed


# ---------------------------------------------------------------------
# 4. Extreme prevalence (Byrt-Bishop-Carlin regime)
# ---------------------------------------------------------------------
tab_rare <- matrix(c(118, 2, 5, 0), nrow = 2,
                   dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
r_rare <- grass_report(tab_rare, format = "matrix")
r_rare


# ---------------------------------------------------------------------
# 5. Quality of the reference curve
# ---------------------------------------------------------------------
# The reference curves come in two calibration levels: high and medium.
# Pick whichever matches the rater quality you are comparing against.
grass_report(tab, format = "matrix", reference = "high")
grass_report(tab, format = "matrix", reference = "medium")

# Or skip the reference entirely — metrics in context, no curves attached.
grass_report(tab, format = "matrix", reference = "none")


# ---------------------------------------------------------------------
# 6. Text-level inputs (factor / character)
# ---------------------------------------------------------------------
yes_no_df <- data.frame(
  rater1 = sample(c("yes", "no"), 60, replace = TRUE, prob = c(0.3, 0.7)),
  rater2 = sample(c("yes", "no"), 60, replace = TRUE, prob = c(0.3, 0.7))
)
grass_compute(yes_no_df, format = "wide")   # one-time message: positive = "yes"
grass_compute(yes_no_df, format = "wide", positive = "no")


# ---------------------------------------------------------------------
# 7. Prevalence override
# ---------------------------------------------------------------------
# Default: estimate prevalence from rater marginals.
# Override with a known population rate:
grass_report(tab_skewed, format = "matrix", prevalence = 0.05)
grass_report(tab_skewed, format = "matrix", prevalence = 0.50)


# ---------------------------------------------------------------------
# 8. The full reference table
# ---------------------------------------------------------------------
head(grass_reference_table(), 5)
# 23 prevalence points; grassr uses linear interpolation between them.


# ---------------------------------------------------------------------
# 9. Plotting
# ---------------------------------------------------------------------
# Landing plot: metrics on reference curves. Inline curve labels by
# default; pass labels = "legend" to use a top legend instead.
plot(r_skewed)
plot(r_skewed, labels = "legend")

# Turn on the medium-quality curves as a faint dotted band:
plot(r_skewed, show_medium = TRUE)

# Regime plot: PI^2 vs BI^2.
plot(r_skewed, type = "regime")

# Just the reference curves at a given prevalence.
plot(r_skewed$reference)

# Override title / subtitle for papers where the caption does the talking.
plot(r_skewed, title = NULL, subtitle = NULL)


# ---------------------------------------------------------------------
# 9b. Paper-ready one-line report
# ---------------------------------------------------------------------
# A single character string suitable for a manuscript sentence or a
# results-table cell. Unicode kappa by default; ascii = TRUE falls back
# to "kappa" for contexts that can't render it.
grass_format_report(r_skewed)
grass_format_report(r_skewed, ascii = TRUE)


# ---------------------------------------------------------------------
# 10. Downstream analysis — tidy output
# ---------------------------------------------------------------------
# Every grassr object supports as.data.frame().
as.data.frame(r_skewed)
as.data.frame(r_skewed$metrics)
as.data.frame(r_skewed$reference)

if (requireNamespace("broom", quietly = TRUE)) {
  broom::tidy(r_skewed)
}


# ---------------------------------------------------------------------
# 11. Synthetic data: drive prevalence and rater quality
# ---------------------------------------------------------------------
simulate_raters <- function(n, prevalence, Se = 0.9, Sp = 0.9) {
  truth <- rbinom(n, 1, prevalence)
  r1 <- ifelse(truth == 1, rbinom(n, 1, Se), rbinom(n, 1, 1 - Sp))
  r2 <- ifelse(truth == 1, rbinom(n, 1, Se), rbinom(n, 1, 1 - Sp))
  data.frame(rater1 = r1, rater2 = r2)
}

set.seed(42)

# Balanced prevalence, competent raters.
df_balanced <- simulate_raters(n = 200, prevalence = 0.50, Se = 0.90, Sp = 0.90)
plot(grass_report(df_balanced, format = "wide"))

# Rare positive class, same raters.
df_skewed <- simulate_raters(n = 200, prevalence = 0.05, Se = 0.90, Sp = 0.90)
plot(grass_report(df_skewed, format = "wide"))

# Introduce a bias mismatch — rater 2 is stricter.
df_biased <- simulate_raters(n = 200, prevalence = 0.30, Se = 0.90, Sp = 0.90)
mask <- df_biased$rater2 == 1 & rbinom(nrow(df_biased), 1, 0.5) == 1
df_biased$rater2[mask] <- 0
plot(grass_report(df_biased, format = "wide"), type = "regime")


# ---------------------------------------------------------------------
# 12. Things worth trying
# ---------------------------------------------------------------------
# * Vary Se/Sp in simulate_raters. Watch kappa drift while PABAK stays flat.
# * Loop prevalence from 0.01 to 0.99 and see when the three metrics diverge.
# * Pass a 2-level factor c("present", "absent"); confirm "present" is mapped to 1.
# * Introduce NAs and watch the pairwise-drop warning fire.
# * Try to break the package: 3-level factor, 3x3 matrix, long with 3 raters.

# End of playground.
