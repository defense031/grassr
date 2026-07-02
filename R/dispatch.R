# Spec-dispatch for the metric families.
#
# grass_report() drives every metric family through the same three
# operations:
#
#   compute_agreement(data, spec, ...)  -> grass_metrics
#   reference_for(spec, context)        -> grass_reference or NULL
#   classify_regime(metrics, spec)      -> list(regime, note)
#
# These are package-internal. Rather than use S3 method names (which
# roxygen flags as needing @export), each operation dispatches manually
# on spec$family via switch(). The intent is the same -- extensible per
# family -- without exporting an internal generic.

# ---- compute_agreement -------------------------------------------------

compute_agreement <- function(data, spec, ...) {
  switch(spec$family,
    binary     = grass_compute(data, ...),
    ordinal    = stop_family_unimplemented("ordinal"),
    multirater = stop_family_unimplemented("multirater"),
    continuous = stop_family_unimplemented("continuous"),
    stop("Unknown spec family: ", spec$family, call. = FALSE))
}

# ---- reference_for -----------------------------------------------------

reference_for <- function(spec, context, ...) {
  switch(spec$family,
    binary     = {
      level <- spec$reference_level
      if (is.null(level)) NULL
      else reference_for_binary(context$prevalence, level)
    },
    ordinal    = stop_family_unimplemented("ordinal"),
    multirater = stop_family_unimplemented("multirater"),
    continuous = stop_family_unimplemented("continuous"),
    stop("Unknown spec family: ", spec$family, call. = FALSE))
}

# ---- classify_regime ---------------------------------------------------
#
# The binary regime is driven by PI and BI; the free function
# classify_regime_binary() holds the existing logic. Future families
# whose structural regimes are not PI/BI-indexed (e.g., a continuous
# family by between/within variance ratio) slot in as new switch arms.

classify_regime <- function(metrics, spec, ...) {
  switch(spec$family,
    binary     = {
      v <- metrics$values
      classify_regime_binary(v["prevalence_index"], v["bias_index"])
    },
    ordinal    = stop_family_unimplemented("ordinal"),
    multirater = stop_family_unimplemented("multirater"),
    continuous = stop_family_unimplemented("continuous"),
    stop("Unknown spec family: ", spec$family, call. = FALSE))
}

# ---- Helpers -----------------------------------------------------------

stop_family_unimplemented <- function(family) {
  stop("The '", family, "' family is declared in the grass roadmap but not ",
       "yet implemented. See `?grass_roadmap` for the planned scope.",
       call. = FALSE)
}
