## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----eval=FALSE---------------------------------------------------------------
#  grass_report(data, spec = grass_spec_binary(), ...)

## -----------------------------------------------------------------------------
library(grass)
tab <- matrix(c(88, 10, 14, 88), nrow = 2,
              dimnames = list(R1 = c("0","1"), R2 = c("0","1")))
r <- grass_report(tab, format = "matrix",
                  spec = grass_spec_binary(reference_level = 0.85))
r

## -----------------------------------------------------------------------------
cat(grass_methods(r))

