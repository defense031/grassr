## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 110
)
library(grass)

## -----------------------------------------------------------------------------
tab <- matrix(c(88, 10, 14, 88), nrow = 2,
              dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
tab

result <- grass_report(tab, format = "matrix")
result

## -----------------------------------------------------------------------------
tab_rare <- matrix(c(40, 5, 5, 0), nrow = 2,
                   dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
result_rare <- grass_report(tab_rare, format = "matrix")
result_rare

## ----eval = requireNamespace("ggplot2", quietly = TRUE)-----------------------
plot(result_rare)

## ----eval = requireNamespace("ggplot2", quietly = TRUE)-----------------------
plot(result_rare, type = "regime")

## -----------------------------------------------------------------------------
r1 <- c(rep(0, 40), rep(0, 5), rep(1, 5))
r2 <- c(rep(0, 40), rep(1, 5), rep(0, 5))

grass_compute(list(r1, r2), format = "paired")$values["kappa"]
grass_compute(data.frame(r1 = r1, r2 = r2), format = "wide")$values["kappa"]

long_df <- data.frame(
  subject = rep(seq_along(r1), 2),
  rater   = rep(c("r1", "r2"), each = length(r1)),
  rating  = c(r1, r2)
)
grass_compute(long_df, format = "long")$values["kappa"]

