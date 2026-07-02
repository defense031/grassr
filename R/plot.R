# Plotting. All entry points require ggplot2; checked at call time so ggplot2
# stays in Suggests.

# Paper palette.
grass_colors <- c(kappa = "#E41A1C", PABAK = "#377EB8", AC1 = "#4DAF4A")

# Diverging gradient from the paper's accuracy heatmap.
grass_gradient_colors <- c(low = "#C1272D", mid = "#FFF3B0", high = "#006837")

check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the 'ggplot2' package. ",
         "Install it with install.packages(\"ggplot2\").",
         call. = FALSE)
  }
}

#' ggplot2 theme for grass plots
#'
#' Minimal publication-oriented theme. No chrome beyond what the data
#' require.
#'
#' @param base_size Base font size passed to `ggplot2::theme_minimal()`.
#' @return A ggplot2 theme object.
#' @export
theme_grass <- function(base_size = 12) {
  check_ggplot2()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      plot.title.position = "plot",
      plot.title    = ggplot2::element_text(face = "bold",
                                            size = ggplot2::rel(1.05),
                                            margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(color = "grey25",
                                            size = ggplot2::rel(0.90),
                                            margin = ggplot2::margin(b = 8)),
      plot.margin   = ggplot2::margin(10, 14, 8, 10),
      axis.title    = ggplot2::element_text(face = "plain"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(0, 0, 2, 0),
      legend.key.width = grid::unit(1.4, "lines")
    )
}

#' Color scale for grass metrics
#'
#' Maps `kappa`, `PABAK`, and `AC1` to the paper palette.
#'
#' @param ... Passed through to `ggplot2::scale_color_manual` /
#'   `ggplot2::scale_fill_manual`.
#' @return A ggplot2 scale.
#' @name scale_grass_metric
#' @export
scale_color_grass_metric <- function(...) {
  check_ggplot2()
  ggplot2::scale_color_manual(
    values = grass_colors,
    breaks = c("kappa", "PABAK", "AC1"),
    labels = c("kappa", "PABAK", "AC1"),
    guide  = ggplot2::guide_legend(override.aes = list(
      linewidth = 1.4, shape = 16, size = 3, label = NA, fill = NA
    )),
    ...
  )
}

#' @rdname scale_grass_metric
#' @export
scale_fill_grass_metric <- function(...) {
  check_ggplot2()
  ggplot2::scale_fill_manual(values = grass_colors, ...)
}

# Short interpretive caption for plot subtitles. The `regime` slot on a
# grass_result carries a bare label ("prevalence-dominated" etc.) that on
# its own reads as jargon. We append a one-clause interpretation so a
# reader who has never opened the package help can tell -- at a glance --
# why the three metrics disagree the way they do.
regime_caption <- function(regime) {
  if (is.null(regime) || is.na(regime)) return("regime = unknown")
  switch(regime,
    "balanced" =
      "balanced regime (kappa, PABAK, AC1 algebraically close)",
    "prevalence-dominated" =
      "prevalence-dominated (kappa suppressed; PABAK/AC1 insulated)",
    "bias-dominated" =
      "bias-dominated (AC1 trails kappa/PABAK)",
    "mixed" =
      "mixed regime (report kappa, PABAK, AC1 alongside PI and BI)",
    paste0("regime = ", regime))
}

# ---- Shared helpers ----------------------------------------------------

# Map legacy "high"/"medium" quality labels to their numeric reference_level.
# Pre-0.1.2 plot code requested curves by string; keep the alias working.
quality_to_level <- function(quality) {
  switch(quality,
         high   = 0.85,
         medium = 0.70,
         stop("Unknown quality label: ", quality, call. = FALSE))
}

# Curves for one level. Accepts either a numeric band (0.70/0.80/0.85/
# 0.90) OR a legacy string ("high" = 0.85, "medium" = 0.70, "both" =
# union of high + medium). The plot API calls this with either form
# depending on whether the result was produced by a spec with a numeric
# reference_level or by the legacy quality-string path.
threshold_curves_df <- function(ref = "high") {
  rb <- reference_binary
  pick <- function(lvl, label) {
    sub <- rb[abs(rb$reference_level - lvl) < 1e-8, , drop = FALSE]
    data.frame(
      prevalence = sub$prevalence,
      metric     = sub$metric,
      threshold  = sub$reference,
      quality    = label,
      stringsAsFactors = FALSE
    )
  }
  if (is.numeric(ref)) {
    lvl <- ref
    return(pick(lvl, format(lvl, nsmall = 2)))
  }
  if (ref == "both") {
    return(rbind(pick(0.85, "high"), pick(0.70, "medium")))
  }
  lvl <- quality_to_level(ref)
  pick(lvl, ref)
}

# Right-edge x-coordinate for inline curve labels. 0.99 is the densest grid
# point; we label just past it.
inline_x <- 1.005

inline_curve_labels <- function(ref = "high") {
  lvl <- if (is.numeric(ref)) ref else quality_to_level(ref)
  rb <- reference_binary
  sub <- rb[abs(rb$reference_level - lvl) < 1e-8, , drop = FALSE]
  metrics <- c("kappa", "PABAK", "AC1")
  tips <- vapply(metrics, function(m) {
    rows <- sub[sub$metric == m, , drop = FALSE]
    rows$reference[which.max(rows$prevalence)]
  }, numeric(1))
  data.frame(
    metric     = metrics,
    prevalence = inline_x,
    threshold  = tips,
    stringsAsFactors = FALSE
  )
}

# ---- Landing plot -----------------------------------------------------

#' Plot a grass report
#'
#' The default `"landing"` plot overlays observed kappa, PABAK, and AC1 on
#' the GRASS prevalence-conditioned reference curves. The `"regime"` plot
#' places the study on the PI^2 vs BI^2 plane.
#'
#' Reference curves are the Youden-J-optimal metric values from the GRASS
#' simulation across prevalence, for competent raters (Se, Sp at or above
#' 0.85 for `"high"` quality; 0.70 for `"medium"`).
#'
#' @param x A `grass_result` object.
#' @param type Either `"landing"` (default) or `"regime"`.
#' @param labels Where to identify the three curves: `"auto"` (default --
#'   inline at the right edge unless curves are too close together),
#'   `"inline"`, or `"legend"`.
#' @param show_medium If `TRUE`, also plot the medium-quality reference
#'   curves as a faint dotted band. Default `FALSE`.
#' @param title,subtitle Optional strings to override the default title and
#'   subtitle. Pass `NULL` to omit either.
#' @param ... Unused.
#'
#' @return A ggplot object.
#' @export
plot.grass_result <- function(x, type = c("landing", "regime"),
                              labels = c("auto", "inline", "legend"),
                              show_medium = FALSE,
                              title = NULL, subtitle = NULL, ...) {
  check_ggplot2()
  type   <- match.arg(type)
  labels <- match.arg(labels)

  if (type == "landing") {
    plot_landing(x, labels = labels, show_medium = show_medium,
                 title = title, subtitle = subtitle)
  } else {
    plot_regime(x, title = title, subtitle = subtitle)
  }
}

# Auto-decide inline vs legend. Inline works when curve endpoints at x=0.99
# are reasonably separated (>= 0.06 between nearest pair). Otherwise legend.
choose_label_placement <- function(ref = "high") {
  endpoints <- inline_curve_labels(ref)$threshold
  diffs <- diff(sort(endpoints))
  if (any(diffs < 0.06)) "legend" else "inline"
}

# Pull the active reference level from a grass_result. Spec-driven
# results store it on the spec; legacy results with a quality string
# fall back to that. NULL means "no reference attached".
active_reference_level <- function(x) {
  if (!is.null(x$spec) && !is.null(x$spec$reference_level)) {
    return(x$spec$reference_level)
  }
  if (!is.null(x$reference) && !is.null(x$reference$reference_level) &&
      !is.na(x$reference$reference_level)) {
    return(x$reference$reference_level)
  }
  if (!is.null(x$reference) && !is.null(x$reference$quality) &&
      !is.na(x$reference$quality)) {
    return(quality_to_level(x$reference$quality))
  }
  0.85
}

plot_landing <- function(x, labels = "auto", show_medium = FALSE,
                         title = NULL, subtitle = NULL) {
  ref_level <- active_reference_level(x)
  if (labels == "auto") labels <- choose_label_placement(ref_level)

  curves_primary <- threshold_curves_df(ref_level)
  obs <- data.frame(
    metric = c("kappa", "PABAK", "AC1"),
    value  = c(unname(x$metrics$values["kappa"]),
               unname(x$metrics$values["PABAK"]),
               unname(x$metrics$values["AC1"])),
    prevalence = x$prevalence,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot()

  # Vertical anchor at the observed prevalence.
  p <- p +
    ggplot2::geom_vline(xintercept = x$prevalence,
                        color = "grey80", linewidth = 0.4)

  # Optional medium-band backdrop.
  if (isTRUE(show_medium) && abs(ref_level - 0.85) < 1e-8) {
    curves_med <- threshold_curves_df("medium")
    p <- p +
      ggplot2::geom_line(
        data = curves_med,
        ggplot2::aes(x = prevalence, y = threshold, color = metric),
        linewidth = 0.4, linetype = "dotted", alpha = 0.6, show.legend = FALSE
      )
  }

  # Primary reference curves.
  p <- p +
    ggplot2::geom_line(
      data = curves_primary,
      ggplot2::aes(x = prevalence, y = threshold, color = metric),
      linewidth = 1, alpha = 0.9
    )

  # Observed-metric points with white ring for contrast.
  p <- p +
    ggplot2::geom_point(
      data = obs,
      ggplot2::aes(x = prevalence, y = value, color = metric),
      size = 4.5, shape = 21, fill = "white", stroke = 1.2, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = obs,
      ggplot2::aes(x = prevalence, y = value, color = metric),
      size = 2, show.legend = FALSE
    )

  # Inline curve labels at the right edge when requested.
  if (labels == "inline") {
    tip <- inline_curve_labels(ref_level)
    p <- p +
      ggplot2::geom_text(
        data = tip,
        ggplot2::aes(x = prevalence, y = threshold, color = metric, label = metric),
        hjust = 0, nudge_x = 0.005, size = 3.3, fontface = "bold",
        show.legend = FALSE
      )
  }

  # Scales, labels, theme.
  x_limits <- if (labels == "inline") c(-0.02, 1.12) else c(-0.02, 1.02)

  default_title <- "Agreement metrics on GRASS reference curves"
  default_subtitle <- sprintf("Prevalence = %.2f   |   N = %d   |   %s",
                              x$prevalence, x$metrics$n,
                              regime_caption(x$regime))
  use_title    <- if (is.null(title))    default_title    else title
  use_subtitle <- if (is.null(subtitle)) default_subtitle else subtitle

  p <- p +
    scale_color_grass_metric() +
    ggplot2::scale_x_continuous(limits = x_limits, breaks = seq(0, 1, 0.1),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-1.05, 1.05), breaks = seq(-1, 1, 0.25)) +
    ggplot2::labs(
      title = use_title,
      subtitle = use_subtitle,
      x = "Prevalence of positive class",
      y = "Metric value"
    ) +
    theme_grass()

  if (labels == "inline") {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p
}

plot_regime <- function(x, title = NULL, subtitle = NULL) {
  pi  <- unname(x$metrics$values["prevalence_index"])
  bi  <- unname(x$metrics$values["bias_index"])
  pi2 <- pi^2
  bi2 <- bi^2
  # Nudge the point off the axes when PI^2 or BI^2 is (near) zero so the
  # marker is not bisected by the plot edge. Cosmetic only -- the printed
  # PI / BI in the subtitle still report the true values.
  nudge <- 0.015
  df_point <- data.frame(PI2 = pmax(pi2, nudge),
                         BI2 = pmax(bi2, nudge))

  default_title    <- "Skew regime"
  default_subtitle <- sprintf("PI = %.3f   |   BI = %.3f   |   %s",
                              pi, bi, regime_caption(x$regime))
  use_title    <- if (is.null(title))    default_title    else title
  use_subtitle <- if (is.null(subtitle)) default_subtitle else subtitle

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = data.frame(x = c(0, 1, 1), y = c(0, 0, 1)),
      ggplot2::aes(x = x, y = y),
      fill = "#E41A1C", alpha = 0.07
    ) +
    ggplot2::geom_polygon(
      data = data.frame(x = c(0, 1, 0), y = c(0, 1, 1)),
      ggplot2::aes(x = x, y = y),
      fill = "#377EB8", alpha = 0.07
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey50",
                         linewidth = 0.4) +
    ggplot2::annotate("text", x = 0.78, y = 0.10,
                      label = "prevalence-dominated", color = "#C1272D",
                      size = 3.4, fontface = "plain") +
    ggplot2::annotate("text", x = 0.22, y = 0.90,
                      label = "bias-dominated", color = "#1B6DB8",
                      size = 3.4, fontface = "plain") +
    ggplot2::geom_point(data = df_point,
                        ggplot2::aes(x = PI2, y = BI2),
                        size = 5, color = "black", shape = 21,
                        stroke = 1.3, fill = "white") +
    ggplot2::geom_point(data = df_point,
                        ggplot2::aes(x = PI2, y = BI2),
                        size = 2, color = "black") +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(
      title = use_title,
      subtitle = use_subtitle,
      x = expression(PI^2),
      y = expression(BI^2)
    ) +
    theme_grass()
}

# ---- Standalone reference curves --------------------------------------

#' @export
plot.grass_reference <- function(x, labels = c("auto", "inline", "legend"),
                                 show_medium = FALSE, ...) {
  check_ggplot2()
  labels <- match.arg(labels)
  lvl <- x$reference$reference_level[1]
  if (is.null(lvl) || is.na(lvl)) lvl <- 0.85
  if (labels == "auto") labels <- choose_label_placement(lvl)

  show_secondary <- isTRUE(show_medium) && abs(lvl - 0.85) < 1e-8
  curves <- threshold_curves_df(if (show_secondary) "both" else lvl)
  current <- data.frame(
    metric     = x$reference$metric,
    threshold  = x$reference$reference,
    prevalence = x$prevalence,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = x$prevalence,
                        color = "grey80", linewidth = 0.4)

  curves_primary   <- if (show_secondary) curves[curves$quality == "high", ] else curves
  curves_secondary <- if (show_secondary) curves[curves$quality == "medium", ] else curves[integer(0), ]
  if (nrow(curves_secondary) > 0) {
    p <- p +
      ggplot2::geom_line(data = curves_secondary,
                         ggplot2::aes(x = prevalence, y = threshold, color = metric),
                         linewidth = 0.45, linetype = "dotted", alpha = 0.7,
                         show.legend = FALSE)
  }
  p <- p +
    ggplot2::geom_line(data = curves_primary,
                       ggplot2::aes(x = prevalence, y = threshold, color = metric),
                       linewidth = 1) +
    ggplot2::geom_point(data = current,
                        ggplot2::aes(x = prevalence, y = threshold, color = metric),
                        size = 4.5, shape = 21, fill = "white", stroke = 1.2,
                        show.legend = FALSE) +
    ggplot2::geom_point(data = current,
                        ggplot2::aes(x = prevalence, y = threshold, color = metric),
                        size = 2, show.legend = FALSE)

  if (labels == "inline") {
    tip <- inline_curve_labels(lvl)
    p <- p +
      ggplot2::geom_text(data = tip,
                         ggplot2::aes(x = prevalence, y = threshold, color = metric,
                                      label = metric),
                         hjust = 0, nudge_x = 0.005, size = 3.3,
                         fontface = "bold", show.legend = FALSE)
  }

  x_limits <- if (labels == "inline") c(-0.02, 1.12) else c(-0.02, 1.02)

  p <- p +
    scale_color_grass_metric() +
    ggplot2::scale_x_continuous(limits = x_limits, breaks = seq(0, 1, 0.1),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(
      title = sprintf("GRASS reference curves (Se = Sp = %.2f)", lvl),
      subtitle = sprintf("Evaluation prevalence = %.2f", x$prevalence),
      x = "Prevalence of positive class",
      y = "Reference value"
    ) +
    theme_grass()

  if (labels == "inline") {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}

# ---- Bare metrics plot ------------------------------------------------

# Dropped -- bare metrics have no prevalence or bias context. Guide the user
# toward the report workflow.
#' @export
plot.grass_metrics <- function(x, ...) {
  stop("Plotting raw metrics without prevalence context is not supported. ",
       "Use grass_report() and plot() the resulting grass_result.",
       call. = FALSE)
}

# ---- Generic entry point ----------------------------------------------

#' Plot a grass object
#'
#' Dispatches to `plot.grass_result` or `plot.grass_reference`.
#'
#' @param x A grass object.
#' @param ... Passed through to the class-specific plot method.
#' @return A ggplot object.
#' @export
grass_plot <- function(x, ...) {
  check_ggplot2()
  UseMethod("grass_plot")
}

#' @export
grass_plot.default <- function(x, ...) {
  plot(x, ...)
}
