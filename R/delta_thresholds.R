# delta_thresholds.R
# Per-(k, N) threshold lookup for the cross-coefficient asymmetry flag.
#
# Background. The original (9.25, 11.75) pair was calibrated at the modal
# inter-rater design (k = 5, N = 1000). When applied at smaller N or k the
# size-alpha guarantee fails (FPR climbs above the nominal alpha = 0.05).
# `paper1_2_merged/scripts/run_delta_roc_smallN.R` re-derives the (caution,
# divergent) cuts cell-by-cell on a (k in {2,3,5,8,15,25}) x (N in
# {15,20,30,50,75,100,200,500}) grid. The resulting `delta_thresholds_lookup`
# table is bundled in sysdata.rda; this file gives the access path.
#
# Behavior policy (transparency over silence):
#   - exact (k, N) match with BOTH thresholds calibrated -> use directly
#   - exact match has NA, or (k, N) is off-grid -> search the table for
#     the nearest cell where both thresholds are calibrated (Euclidean
#     distance in (k, log10(N)) space, each axis range-normalized so
#     k and N contribute comparably). Surface the snap.
#   - table missing entirely -> hard fallback to modal default
#     (9.25, 11.75); should not happen with sysdata bundled correctly.
#   - user-supplied delta_thresholds -> always honored, no lookup
#
# Why search the calibrated set rather than fall back to a default at
# the first NA: the modal default (9.25, 11.75) was calibrated at
# (k = 5, N = 1000); at small (k, N) it produces FPR an order of
# magnitude above the nominal alpha. A pair calibrated somewhere on
# the small-design regime, even if not at the user's exact (k, N), is
# closer to size-controlling than the modal default.

#' Look up calibrated cross-coefficient asymmetry thresholds at (k, N)
#'
#' Returns the per-(k, N) calibrated `(caution, divergent)` threshold pair
#' from the bundled `delta_thresholds_lookup` table, or the default
#' modal-design pair `c(9.25, 11.75)` when no calibrated entry is available.
#' Always returns the source so callers can surface the choice in
#' `grass_report()` output.
#'
#' @param k Integer panel size (number of raters).
#' @param N Integer panel sample size.
#'
#' @return A list with
#' \describe{
#'   \item{thresholds}{Length-2 numeric in pp, named c("caution", "divergent").}
#'   \item{source}{One of `"calibrated_at_k_N"`,
#'     `"snapped_to_nearest_calibrated"`, or `"default_fallback"`. The
#'     last is reserved for the (impossible-with-bundled-sysdata) case
#'     of an empty table.}
#'   \item{note}{Human-readable explanation suitable for the card's notes.
#'     Empty string when source = `"calibrated_at_k_N"`.}
#'   \item{snapped_k}{The k value actually used (NA when no snap was needed).}
#'   \item{snapped_N}{The N value actually used (NA when no snap was needed).}
#' }
#'
#' @keywords internal
lookup_delta_thresholds <- function(k, N) {
  default <- c(caution = 9.25, divergent = 11.75)

  tbl <- tryCatch(
    get0("delta_thresholds_lookup", envir = asNamespace("grass"),
         inherits = FALSE),
    error = function(e) NULL
  )
  if (is.null(tbl) || !nrow(tbl)) {
    return(list(
      thresholds = default,
      source     = "default_fallback",
      note       = paste0("delta_thresholds_lookup table not bundled in ",
                          "sysdata.rda; using modal-design default (9.25, ",
                          "11.75)."),
      snapped_k  = NA_integer_,
      snapped_N  = NA_integer_
    ))
  }

  k <- as.integer(k)
  N <- as.integer(N)

  # 1. Exact (k, N) match with both thresholds calibrated -> use directly.
  exact_row <- tbl[tbl$k == k & tbl$N == N, , drop = FALSE]
  if (nrow(exact_row) > 0L &&
      is.finite(exact_row$t_caution[1L]) &&
      is.finite(exact_row$t_divergent[1L])) {
    return(list(
      thresholds = c(caution   = unname(exact_row$t_caution[1L]),
                     divergent = unname(exact_row$t_divergent[1L])),
      source     = "calibrated_at_k_N",
      note       = "",
      snapped_k  = NA_integer_,
      snapped_N  = NA_integer_
    ))
  }

  # 2. Restrict to fully-calibrated cells. Among those, find the cell
  # nearest in (k, log10(N)) space with each axis range-normalized.
  cal <- tbl[is.finite(tbl$t_caution) & is.finite(tbl$t_divergent), ,
             drop = FALSE]
  if (nrow(cal) == 0L) {
    return(list(
      thresholds = default,
      source     = "default_fallback",
      note       = paste0("no fully-calibrated cell in delta_thresholds_lookup; ",
                          "using modal-design default (9.25, 11.75)."),
      snapped_k  = NA_integer_,
      snapped_N  = NA_integer_
    ))
  }

  k_range <- max(diff(range(cal$k)), 1L)
  N_range <- max(diff(range(log10(cal$N))), 1e-6)
  k_norm  <- (cal$k - k) / k_range
  N_norm  <- (log10(cal$N) - log10(N)) / N_range
  d <- sqrt(k_norm^2 + N_norm^2)
  i_nearest <- which.min(d)
  row <- cal[i_nearest, , drop = FALSE]

  list(
    thresholds = c(caution   = unname(row$t_caution[1L]),
                   divergent = unname(row$t_divergent[1L])),
    source     = "snapped_to_nearest_calibrated",
    note       = sprintf(
      "(k = %d, N = %d) has no fully-calibrated threshold pair on the grid (or is off-grid); using nearest fully-calibrated cell (k = %d, N = %d).",
      k, N, as.integer(row$k[1L]), as.integer(row$N[1L])),
    snapped_k  = as.integer(row$k[1L]),
    snapped_N  = as.integer(row$N[1L])
  )
}
