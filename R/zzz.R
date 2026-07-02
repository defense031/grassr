# Package-level environment for one-time message tracking.
.grass_env <- new.env(parent = emptyenv())
.grass_env$msg_seen <- list()

.onLoad <- function(libname, pkgname) {
  # Reset one-time message state on load.
  .grass_env$msg_seen <- list()

  # Register broom::tidy methods only if broom is available. Avoids a hard
  # Imports on broom.
  if (requireNamespace("broom", quietly = TRUE)) {
    registerS3method("tidy", "grass_metrics", tidy.grass_metrics,
                     envir = asNamespace("broom"))
    registerS3method("tidy", "grass_reference", tidy.grass_reference,
                     envir = asNamespace("broom"))
  }
  invisible()
}

.onUnload <- function(libpath) {
  .grass_env$msg_seen <- list()
}
