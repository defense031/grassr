## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4.5, dpi = 110)
library(grass)

## -----------------------------------------------------------------------------
tab <- matrix(c(118, 2, 5, 0), nrow = 2,
              dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
result <- grass_report(tab, format = "matrix")
result

## -----------------------------------------------------------------------------
tab <- matrix(c(40, 5, 5, 0), nrow = 2,
              dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
result <- grass_report(tab, format = "matrix")
result

## -----------------------------------------------------------------------------
tab <- matrix(c(60, 10, 25, 5), nrow = 2,
              dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
result <- grass_report(tab, format = "matrix")
result

## ----eval = requireNamespace("ggplot2", quietly = TRUE)-----------------------
plot(result, type = "regime")

## -----------------------------------------------------------------------------
# Helper: simulate two raters at Se = Sp = q against a latent truth.
draw <- function(n, p, q, seed = 1) {
  set.seed(seed)
  truth <- rbinom(n, 1, p)
  r1 <- ifelse(truth == 1, rbinom(n, 1, q), rbinom(n, 1, 1 - q))
  r2 <- ifelse(truth == 1, rbinom(n, 1, q), rbinom(n, 1, 1 - q))
  data.frame(r1 = r1, r2 = r2)
}

## ----fig.show = "hold", out.width = "50%", fig.height = 4.2, eval = requireNamespace("ggplot2", quietly = TRUE)----
r_small <- grass_report(draw(n =  50, p = 0.10, q = 0.85), reference = "high")
r_large <- grass_report(draw(n = 500, p = 0.10, q = 0.85), reference = "high")
plot(r_small, labels = "legend") + ggplot2::labs(title = "N = 50")
plot(r_large, labels = "legend") + ggplot2::labs(title = "N = 500")

## ----fig.show = "hold", out.width = "50%", fig.height = 4.2, eval = requireNamespace("ggplot2", quietly = TRUE)----
r_skew <- grass_report(draw(n = 50, p = 0.10, q = 0.85), reference = "high")
r_bal  <- grass_report(draw(n = 50, p = 0.50, q = 0.85), reference = "high")
plot(r_skew, labels = "legend") + ggplot2::labs(title = "Prevalence = 0.10 (skewed)")
plot(r_bal,  labels = "legend") + ggplot2::labs(title = "Prevalence = 0.50 (balanced)")

