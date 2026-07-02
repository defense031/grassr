# S3 constructors for grass objects.

new_grass_metrics <- function(values, n, table, positive_level,
                              n_dropped = 0L, call = NULL) {
  structure(
    list(values = values,
         n = n,
         table = table,
         positive_level = positive_level,
         n_dropped = n_dropped,
         call = call),
    class = "grass_metrics"
  )
}

new_grass_reference <- function(prevalence, quality, reference,
                                call = NULL) {
  structure(
    list(prevalence = prevalence,
         quality = quality,
         reference = reference,
         call = call),
    class = "grass_reference"
  )
}
