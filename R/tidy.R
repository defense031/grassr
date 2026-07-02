# broom::tidy methods -- registered conditionally in zzz.R so `broom` is not
# a hard dependency.

tidy.grass_metrics <- function(x, ...) {
  v <- x$values
  data.frame(
    term     = names(v),
    estimate = unname(as.numeric(v)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

tidy.grass_reference <- function(x, ...) {
  x$reference
}
