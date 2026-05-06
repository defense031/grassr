## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4.5, dpi = 110)
library(grass)

## -----------------------------------------------------------------------------
se <- c(R1 = 0.88, R2 = 0.86, R3 = 0.89, R4 = 0.87, R5 = 0.85)
sp <- c(R1 = 0.86, R2 = 0.87, R3 = 0.88, R4 = 0.85, R5 = 0.86)

## -----------------------------------------------------------------------------
safety <- check_asymmetry(se = se, sp = sp, rater = names(se))
safety

## -----------------------------------------------------------------------------
q_hat <- 0.87
emr_point <- emr_panel(q = q_hat, k = 5)
# Rough 95% CI on EMR via propagating q_hat ± 1.96 * 0.02 through emr_panel()
emr_lo <- emr_panel(q = q_hat + 1.96 * 0.02, k = 5)  # lower EMR = upper q
emr_hi <- emr_panel(q = q_hat - 1.96 * 0.02, k = 5)  # upper EMR = lower q
c(point = emr_point, lower = emr_lo, upper = emr_hi)

