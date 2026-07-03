# =====================================================================
#  grassr — first tour
#
#  A short walkthrough of what the grassr package does. Runs top-to-bottom
#  in about 30 seconds and hits the headline features that support Paper 1.
# =====================================================================

# ---- 0. Install & load --------------------------------------------------
# If you don't have remotes yet:  install.packages("remotes")
# remotes::install_github("defense031/grassr")

library(grassr)
library(ggplot2)

# Optional companion vignettes — browse after install:
#   vignette("grass-intro", package = "grassr")
#   vignette("skewed-examples", package = "grassr")
#   vignette("grass-ecosystem", package = "grassr")


# ---- 1. A real example dataset ------------------------------------------
# 200 chest X-rays independently read by two radiologists as abnormal/normal.
# Intentional prevalence skew and rater bias so the framework earns its keep.
cxr <- read.csv(system.file("examples", "cxr_review.csv", package = "grassr"),
                stringsAsFactors = FALSE)
head(cxr)
table(cxr$rater_A, cxr$rater_B, dnn = c("rater_A", "rater_B"))


# ---- 2. The headline call -----------------------------------------------
# Pass the data frame, tell it which column is the subject ID, tell it
# which text value is the clinical positive. Everything else is figured out.
result <- grass_report(cxr,
                       id_col   = "xray_id",
                       positive = "abnormal")
result


# ---- 3. The landing plot ------------------------------------------------
# Observed kappa / PABAK / AC1 overlaid on their prevalence-conditioned
# reference curves at the study's observed prevalence. This is the figure
# we want in the paper.
plot(result)

# Skew regime diagnostic: PI^2 vs BI^2 tells you whether the metric gap
# is driven by prevalence, bias, or both.
plot(result, type = "regime")


# ---- 4. Manuscript methods paragraph (THE headline feature) -------------
# Drop-in for your methods section. Pre-filled with the study's numbers,
# GRRAS-compliant, available in three formats.
cat(grass_methods(result, format = "markdown"))
# cat(grass_methods(result, format = "latex"))
# cat(grass_methods(result, format = "plain"))


# ---- 5. One-line paper-ready summary ------------------------------------
# For results tables or Slack threads.
grass_format_report(result)
grass_format_report(result, ci_width = TRUE)


# ---- 6. Downstream analysis --------------------------------------------
# Every grassr result coerces cleanly to a one-row data.frame for dashboards
# or results-table pipelines.
as.data.frame(result)

# Long form for ggplot / gt:
broom::tidy(result)


# ---- 7. Cohort-split reporting -----------------------------------------
# One row per cohort, tidy, ready for bind_rows or write.csv.
set.seed(42)
cxr$site <- sample(c("site1", "site2", "site3"), nrow(cxr), replace = TRUE)
grass_report_by(cxr, site, id_col = "xray_id", positive = "abnormal")


# ---- 8. Where to go next ----------------------------------------------
# See the ROADMAP in the package for what's coming:
#   - TURF family (ordinal, weighted kappa)
#   - MEADOW family (more than two raters)
#   - FIELD family (continuous / ICC)
# All dispatch through the same grass_report() call via different spec
# constructors. The binary family above is implemented today.
# ?grass_roadmap
