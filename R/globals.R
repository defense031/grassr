# ggplot2 non-standard-evaluation column names used in plot-card.R /
# plot.R aesthetics. Declared so R CMD check does not flag them as
# undefined globals.
utils::globalVariables(c(
  "coefficient", "col_rater", "estimate", "fill", "M1",
  "observed_value", "pabak_label", "percentile", "pos", "rater",
  "row_rater", "surface_percentile", "xmax", "xmin", "y_pos"
))
