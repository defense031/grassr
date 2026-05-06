# plot.grass_card — multi-view S3 method for the v0.2.0 Report Card.
#
# Six views, all returning ggplot objects (or a patchwork composite for
# `type = "diagnostic"`). The flagship "surface" view is a 2D heatmap of
# E[primary metric] over (M_1, q) at the study's (k, N), with the observation
# pinned and dotted band boundaries at q in {0.5, 0.625, 0.75, 0.875, 1.0}.
#
# Closed-form expectations for non-ICC metrics follow paper2 §6.1
# (paper2/code/04_reference_closed_form.R is the reference implementation;
# we reimplement the large-N closed forms inline here so the package has no
# project-path runtime dependency). For ICC we prefer
# `card$surface$reference_curves$icc` if available; otherwise we fall back to
# a non-ICC primary with a note.

# ---- package availability check ----------------------------------------

check_patchwork <- function() {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("The 'diagnostic' view requires the 'patchwork' package. ",
         "Install it with install.packages(\"patchwork\").",
         call. = FALSE)
  }
}

# ---- closed-form surface helpers ---------------------------------------
#
# Large-N closed forms for E[metric] under the clustered latent-class DGP
# with diagonal Se = Sp = q and overall positive marginal pi_+. We
# parameterize the surface by (M_1, q) where M_1 = E[p] (the mean prevalence
# under F). pi_+ = (1 - q) + (2 q - 1) * M_1.
#
# These reduce to F-only-through-pi_+ for the four agreement-family metrics
# (PABAK, Fleiss kappa, mean-pairwise AC1, Krippendorff alpha) — see paper2
# §6.1, App A.2-A.3.

.cf_pi_plus <- function(q, M1) (1 - q) + (2 * q - 1) * M1

.cf_pabak  <- function(q, M1) (2 * q - 1)^2

.cf_mean_pabak <- function(q, M1) (2 * q - 1)^2  # same form

.cf_fleiss_kappa <- function(q, M1) {
  pi_p <- .cf_pi_plus(q, M1)
  P_e  <- pi_p^2 + (1 - pi_p)^2
  P_a  <- 1 - 2 * q * (1 - q)
  ifelse(abs(1 - P_e) < 1e-12, NA_real_, (P_a - P_e) / (1 - P_e))
}

.cf_mean_ac1 <- function(q, M1) {
  pi_p <- .cf_pi_plus(q, M1)
  p_a  <- 1 - 2 * q * (1 - q)
  p_e  <- 2 * pi_p * (1 - pi_p)
  ifelse(abs(1 - p_e) < 1e-12, NA_real_, (p_a - p_e) / (1 - p_e))
}

.cf_krippendorff_alpha <- function(q, M1) {
  pi_p <- .cf_pi_plus(q, M1)
  denom <- pi_p * (1 - pi_p)
  ifelse(denom < 1e-12, NA_real_, 1 - q * (1 - q) / denom)
}

# Two-rater Cohen's kappa under symmetric Se = Sp = q, identical marginals:
# the per-pair p_a / p_e structure is the same as Fleiss in this special case.
# We treat "kappa" for plotting purposes as Fleiss' large-N form.
.cf_kappa <- .cf_fleiss_kappa

# Dispatcher for which closed form to use.
.cf_metric_dispatch <- function(metric) {
  switch(metric,
    "pabak"          = .cf_pabak,
    "mean_pabak"     = .cf_mean_pabak,
    "ac1"            = .cf_mean_ac1,
    "mean_ac1"       = .cf_mean_ac1,
    "fleiss_kappa"   = .cf_fleiss_kappa,
    "kappa"          = .cf_kappa,
    "krippendorff_a" = .cf_krippendorff_alpha,
    NULL)
}

# 2D grid of E[metric] over (M_1, q). Returns a long-format data frame for
# ggplot's geom_raster / geom_contour.
.surface_grid <- function(metric, n_grid_M1 = 60, n_grid_q = 60,
                          M1_range = c(0.02, 0.98), q_range = c(0.5, 1.0)) {
  fn <- .cf_metric_dispatch(metric)
  if (is.null(fn)) {
    return(NULL)  # caller decides how to handle (e.g. ICC fallback)
  }
  M1_seq <- seq(M1_range[1], M1_range[2], length.out = n_grid_M1)
  q_seq  <- seq(q_range[1],  q_range[2],  length.out = n_grid_q)
  grid <- expand.grid(M1 = M1_seq, q = q_seq, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  grid$value <- fn(grid$q, grid$M1)
  grid
}

# ---- standalone surface plot (prospective design) ----------------------

#' Standalone reference-surface plot for prospective study design
#'
#' Renders the same context-conditioned reference surface that
#' `plot.grass_card(type = "surface")` shows, but driven by scalar inputs
#' rather than a fitted [grass_report()] result. The intended use case is
#' prospective study design: a practitioner planning a multi-rater study at
#' design `(pi_hat, k, N)` can ask "what does my reference surface look
#' like for this metric, before I collect data?" — and `plot_surface()`
#' answers without requiring a rating matrix.
#'
#' The plot is a heatmap of the closed-form expectation
#' `E[metric](M_1, q)` over the surface, where `M_1` is the marginal
#' positive rate and `q in [0.5, 1]` is the diagonal rater operating
#' quality. Dotted contours overlay the four-band partition
#' (`bands` argument). When `pi_hat` is supplied alone, a dashed vertical
#' line marks the design's marginal. When both `pi_hat` and `observed`
#' are supplied, a filled marker pins the implied `(M_1, hat q)` point on
#' the surface; `hat q` is recovered by closed-form inversion of `observed`
#' against the reference curve at `pi_hat`.
#'
#' This function uses the closed-form non-ICC surfaces only (PABAK,
#' Fleiss kappa, mean-pairwise AC1, Krippendorff alpha). ICC surfaces
#' depend on the full F-shape and the GLMM-gap-corrected fitted reference;
#' those are accessed via [position_on_surface()] with `metric = "icc"`
#' and (if available) a fitted card via `plot.grass_card(type = "surface")`.
#'
#' @param metric Character scalar. One of `"pabak"`, `"mean_pabak"`,
#'   `"fleiss_kappa"`, `"kappa"`, `"mean_ac1"`, `"ac1"`,
#'   `"krippendorff_a"`. ICC is not supported (see Details).
#' @param pi_hat Optional numeric scalar in `(0, 1)`. The design's
#'   marginal positive rate. If supplied alone, draws a vertical
#'   reference line on the surface; if supplied together with
#'   `observed`, anchors the pin point.
#' @param observed Optional numeric scalar. An observed (or hypothetical)
#'   coefficient value. When supplied with `pi_hat`, the function
#'   inverts to `hat q` and pins a marker at `(pi_hat, hat q)`.
#' @param k Optional integer. Rater count. Used for the subtitle only;
#'   non-ICC closed-form surfaces are k-invariant.
#' @param N Optional integer. Sample size. Used for the subtitle only;
#'   non-ICC closed-form surfaces are N-invariant in the large-N limit
#'   the surface represents.
#' @param axis One of `"inter"` (default) or `"intra"`. Used for the
#'   subtitle only.
#' @param bands Numeric length-5 vector giving band boundaries on `q`.
#'   Default `c(0.5, 0.625, 0.75, 0.875, 1.0)`.
#' @param ... Reserved for future extension.
#'
#' @return A `ggplot` object.
#' @seealso [plot.grass_card()] for the card-driven surface view that
#'   adds card-specific annotations (band, qualifier, primary-coefficient
#'   pin); [position_on_surface()] for the underlying inversion machinery.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # Bare surface — what does PABAK's reference look like?
#'   plot_surface("pabak")
#'
#'   # Mark a design context (no observed value yet — pre-data).
#'   plot_surface("fleiss_kappa", pi_hat = 0.30, k = 5, N = 200)
#'
#'   # Pin a hypothetical observation on the surface.
#'   plot_surface("mean_ac1", pi_hat = 0.30, observed = 0.62,
#'                k = 5, N = 200)
#' }
#' }
plot_surface <- function(metric,
                         pi_hat = NULL,
                         observed = NULL,
                         k = NULL,
                         N = NULL,
                         axis = c("inter", "intra"),
                         bands = c(0.5, 0.625, 0.75, 0.875, 1.0),
                         ...) {
  check_ggplot2()
  axis <- match.arg(axis)

  allowed <- c("pabak", "mean_pabak", "fleiss_kappa", "kappa",
               "mean_ac1", "ac1", "krippendorff_a")
  if (!is.character(metric) || length(metric) != 1L || !metric %in% allowed) {
    stop("`metric` must be one of: ",
         paste(shQuote(allowed), collapse = ", "),
         ". ICC surfaces require fitted reference curves and are not ",
         "supported via plot_surface(); see ?position_on_surface and ",
         "?plot.grass_card.", call. = FALSE)
  }
  if (!is.numeric(bands) || length(bands) != 5L ||
      any(!is.finite(bands)) || any(diff(bands) <= 0)) {
    stop("`bands` must be a length-5 strictly increasing numeric vector.",
         call. = FALSE)
  }

  grid <- .surface_grid(metric)
  if (is.null(grid)) {
    stop("No closed-form surface available for metric = '", metric, "'.",
         call. = FALSE)
  }

  pin <- NULL
  if (!is.null(pi_hat)) {
    if (!is.numeric(pi_hat) || length(pi_hat) != 1L ||
        !is.finite(pi_hat) || pi_hat <= 0 || pi_hat >= 1) {
      stop("`pi_hat` must be a numeric scalar in (0, 1).", call. = FALSE)
    }
  }
  if (!is.null(observed)) {
    if (!is.numeric(observed) || length(observed) != 1L ||
        !is.finite(observed)) {
      stop("`observed` must be a finite numeric scalar.", call. = FALSE)
    }
    if (is.null(pi_hat)) {
      stop("`observed` requires `pi_hat` to invert to q_hat for the pin.",
           call. = FALSE)
    }
    q_grid_pin <- seq(0.5, 1.0, length.out = 501L)
    ref_curve  <- closed_form_reference_curve(
      metric  = if (metric == "kappa") "fleiss_kappa" else metric,
      pi_plus = pi_hat, q_grid = q_grid_pin
    )
    inv <- invert_metric_to_q(observed, ref_curve, q_grid_pin)
    if (is.finite(inv$q_hat)) {
      pin <- data.frame(M1 = pi_hat, q = inv$q_hat)
    }
  }

  title <- sprintf("E[%s] reference surface", .pretty_metric_name(metric))
  sub_parts <- character(0L)
  if (!is.null(pi_hat)) sub_parts <- c(sub_parts,
                                       sprintf("pi_hat = %.2f", pi_hat))
  if (!is.null(k)) sub_parts <- c(sub_parts,
                                  sprintf("k = %d", as.integer(k)))
  if (!is.null(N)) sub_parts <- c(sub_parts,
                                  sprintf("N = %d", as.integer(N)))
  sub_parts <- c(sub_parts, sprintf("axis = %s", axis))
  if (!is.null(observed)) sub_parts <- c(sub_parts,
                                         sprintf("observed = %.3f", observed))
  subtitle <- paste(sub_parts, collapse = ",  ")

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = M1, y = q, fill = value)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::geom_contour(ggplot2::aes(z = q),
                          breaks = bands,
                          color = "white", linetype = "dotted",
                          linewidth = 0.5)

  if (!is.null(pin) && is.finite(pin$M1) && is.finite(pin$q)) {
    p <- p + ggplot2::geom_point(data = pin, ggplot2::aes(x = M1, y = q),
                                  inherit.aes = FALSE,
                                  size = 4.5, shape = 21, fill = "white",
                                  color = "black", stroke = 1.4)
  } else if (!is.null(pi_hat)) {
    p <- p + ggplot2::geom_vline(xintercept = pi_hat, linetype = "dashed",
                                  color = "white", alpha = 0.7,
                                  linewidth = 0.6)
  }

  if (requireNamespace("viridisLite", quietly = TRUE)) {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = viridisLite::viridis(64),
      name = paste0("E[", .pretty_metric_name(metric), "]"))
  } else {
    p <- p + ggplot2::scale_fill_gradient(
      low = "#440154", high = "#FDE725",
      name = paste0("E[", .pretty_metric_name(metric), "]"))
  }

  p +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                                breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0),
                                breaks = seq(0.5, 1.0, 0.1)) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x = expression(M[1]~"(mean prevalence)"),
      y = expression(q~"(rater operating quality)")
    ) +
    theme_grass()
}

# ---- view 1: surface ---------------------------------------------------

.plot_card_surface <- function(x, ...) {
  check_ggplot2()
  primary  <- x$coefficient$primary %||% "pabak"
  obs_val  <- x$coefficient$observed_value
  pct      <- x$coefficient$surface_percentile
  band     <- x$coefficient$band
  qual     <- x$coefficient$qualifier
  k        <- x$sample$k
  N        <- x$sample$N
  pi_hat   <- x$sample$pi_hat
  axis     <- x$sample$axis %||% "inter"

  primary_for_grid <- primary
  fallback_note <- NULL
  if (identical(primary, "icc")) {
    # ICC's surface needs a fitted reference curve (the GLMM gap). Use the
    # bundled curve if the card carries it; otherwise fall back to a non-ICC
    # primary and note the substitution.
    has_icc_curve <- !is.null(x$surface$reference_curves) &&
                     !is.null(x$surface$reference_curves$icc)
    if (!has_icc_curve) {
      fallback_note <- paste0(
        "ICC surface requires a fitted reference curve; ",
        "displaying PABAK surface as a fallback.")
      primary_for_grid <- "pabak"
    }
  }

  grid <- .surface_grid(primary_for_grid)
  if (is.null(grid)) {
    # No closed form available for this metric — fall back to PABAK.
    fallback_note <- sprintf(
      "No closed-form surface for '%s'; displaying PABAK surface.", primary)
    primary_for_grid <- "pabak"
    grid <- .surface_grid(primary_for_grid)
  }

  # The pinned observation: the observed (pi_hat, q_hat) point.
  q_hat <- NA_real_
  if (!is.null(x$panel) && nrow(x$panel) > 0L) {
    row <- x$panel[x$panel$coefficient == primary, , drop = FALSE]
    if (nrow(row) >= 1L) q_hat <- row$q_hat[1]
  }
  pin <- data.frame(M1 = pi_hat, q = q_hat)

  # Top-line annotation in the four-field card style.
  pct_int <- if (is.numeric(pct) && is.finite(pct)) round(pct) else NA_integer_
  pct_lbl <- if (!is.na(pct_int)) {
    paste0(pct_int, .ord_suffix(pct_int))
  } else "--"
  band_qual <- if (!is.null(band) && !is.na(band) &&
                    !is.null(qual) && !is.na(qual))
    sprintf("(%s, %s)", band, qual) else ""
  card_summary <- sprintf("%s = %.2f -> %s percentile %s",
                          .pretty_metric_name(primary),
                          obs_val %||% NA_real_, pct_lbl, band_qual)
  if (!is.null(fallback_note)) {
    card_summary <- paste0(card_summary, "  |  ", fallback_note)
  }

  subtitle <- sprintf("k = %d, N = %d, pi_hat = %.2f, axis = %s",
                      as.integer(k), as.integer(N),
                      as.numeric(pi_hat), as.character(axis))

  bands <- x$inputs$bands %||% c(0.5, 0.625, 0.75, 0.875, 1.0)

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = M1, y = q, fill = value)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::geom_contour(ggplot2::aes(z = q),
                          breaks = bands,
                          color = "white", linetype = "dotted",
                          linewidth = 0.5)

  # The pinned observation: solid white-filled black-bordered point.
  if (is.finite(pin$M1) && is.finite(pin$q)) {
    p <- p +
      ggplot2::geom_point(data = pin, ggplot2::aes(x = M1, y = q),
                          inherit.aes = FALSE,
                          size = 4.5, shape = 21, fill = "white",
                          color = "black", stroke = 1.4)
  }

  # Color scale: viridis if available; otherwise a manual blue-yellow-red.
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = viridisLite::viridis(64),
      name = paste0("E[", .pretty_metric_name(primary_for_grid), "]"))
  } else {
    p <- p + ggplot2::scale_fill_gradient(
      low = "#440154", high = "#FDE725",
      name = paste0("E[", .pretty_metric_name(primary_for_grid), "]"))
  }

  p <- p +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                                breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0),
                                breaks = seq(0.5, 1.0, 0.1)) +
    ggplot2::labs(
      title    = card_summary,
      subtitle = subtitle,
      x = expression(M[1]~"(mean prevalence)"),
      y = expression(q~"(rater operating quality)")
    ) +
    theme_grass()

  p
}

# Ordinal suffix ("st", "nd", "rd", "th") for an integer percentile.
.ord_suffix <- function(n) {
  if (!is.finite(n)) return("")
  n <- as.integer(n)
  last_two <- n %% 100
  last_one <- n %% 10
  if (last_two %in% 11:13) return("th")
  switch(as.character(last_one),
         "1" = "st", "2" = "nd", "3" = "rd",
         "th")
}

# Pretty display names for the metric labels in plot titles / axis labels.
.pretty_metric_name <- function(metric) {
  switch(metric,
    "pabak"          = "PABAK",
    "mean_pabak"     = "mean-pairwise PABAK",
    "ac1"            = "AC1",
    "mean_ac1"       = "mean-pairwise AC1",
    "fleiss_kappa"   = "Fleiss kappa",
    "kappa"          = "Cohen's kappa",
    "krippendorff_a" = "Krippendorff alpha",
    "icc"            = "ICC",
    metric)
}

# ---- view 2: panel (forest on percentile axis) -------------------------

.plot_card_panel <- function(x, ...) {
  check_ggplot2()
  panel <- x$panel
  if (is.null(panel) || nrow(panel) == 0L) {
    stop("The card has no `panel` data frame to plot.", call. = FALSE)
  }

  # Coefficient names get factor levels for stable ordering, then we map to
  # integer y positions. The integer-y path lets us mix the discrete coefficient
  # labels with numeric annotations (band-name text above the top row) on the
  # same axis without the discrete/continuous scale conflict.
  coef_levels <- rev(panel$coefficient)
  n_rows <- nrow(panel)
  df <- data.frame(
    coefficient = factor(panel$coefficient, levels = coef_levels),
    surface_percentile = panel$surface_percentile,
    band = panel$band,
    qualifier = panel$qualifier,
    stringsAsFactors = FALSE
  )
  df$y_pos <- as.integer(df$coefficient)

  band_breaks <- data.frame(pos = c(25, 50, 75))

  ttl <- sprintf("Panel of coefficients: surface percentile (delta = %.1f pp, %s)",
                 x$delta$delta_hat %||% NA_real_,
                 x$delta$flag %||% "?")
  sub <- sprintf("k = %d, N = %d, pi_hat = %.2f",
                 as.integer(x$sample$k), as.integer(x$sample$N),
                 as.numeric(x$sample$pi_hat))

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = surface_percentile, y = y_pos)) +
    ggplot2::geom_vline(data = band_breaks,
                        ggplot2::aes(xintercept = pos),
                        linetype = "dotted", color = "grey50") +
    ggplot2::annotate("text", x = c(12.5, 37.5, 62.5, 87.5),
                      y = n_rows + 0.4,
                      label = c("Poor", "Moderate", "Strong", "Excellent"),
                      size = 3.0, color = "grey25", fontface = "italic") +
    ggplot2::geom_point(size = 3.5, color = "black") +
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                                expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
    ggplot2::scale_y_continuous(breaks = seq_len(n_rows),
                                labels = levels(df$coefficient),
                                limits = c(0.5, n_rows + 0.7),
                                expand = c(0, 0)) +
    ggplot2::labs(title = ttl, subtitle = sub,
                  x = "surface percentile (pp)",
                  y = NULL) +
    theme_grass()
  p
}

# ---- view 3: thermometer ----------------------------------------------

.plot_card_thermometer <- function(x, ...) {
  check_ggplot2()
  thr <- x$delta$thresholds %||% c(caution = 9.25, divergent = 11.75)
  caut <- unname(thr[1])
  div  <- unname(thr[2])
  delta_hat <- x$delta$delta_hat %||% NA_real_
  flag <- x$delta$flag %||% "?"

  # Three-band rectangle stack; pointer arrow at delta_hat.
  rng <- data.frame(
    xmin = c(0, caut, div),
    xmax = c(caut, div, 100),
    fill = c("aligned", "caution", "divergent")
  )
  fill_pal <- c(aligned = "#A6D96A", caution = "#FEE08B", divergent = "#F46D43")

  ttl <- expression(hat(delta) ~ "asymmetry gauge")
  sub <- sprintf("delta_hat = %.2f pp (%s); thresholds %.2f / %.2f",
                 delta_hat, flag, caut, div)

  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_rect(
    data = rng,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = fill),
    color = NA, alpha = 0.85
  )
  p <- p + ggplot2::scale_fill_manual(values = fill_pal,
                                       breaks = c("aligned", "caution",
                                                  "divergent"),
                                       name = NULL)
  # Pointer.
  if (is.finite(delta_hat)) {
    p <- p + ggplot2::geom_segment(
      data = data.frame(x = delta_hat),
      ggplot2::aes(x = x, xend = x, y = 1.05, yend = 1.5),
      arrow = grid::arrow(angle = 25, length = grid::unit(0.18, "inches"),
                          ends = "first", type = "closed"),
      linewidth = 1.2, color = "black"
    )
    p <- p + ggplot2::annotate("text", x = delta_hat, y = 1.65,
                                label = sprintf("%.1f pp", delta_hat),
                                fontface = "bold", size = 4)
  }
  p <- p +
    ggplot2::scale_x_continuous(limits = c(0, max(50, ceiling(delta_hat) + 5,
                                                   div + 10, na.rm = TRUE)),
                                expand = c(0, 0),
                                breaks = c(0, caut, div,
                                           seq(0, 100, 10))) +
    ggplot2::scale_y_continuous(limits = c(-0.1, 2.1), expand = c(0, 0),
                                breaks = NULL) +
    ggplot2::labs(title = ttl, subtitle = sub,
                  x = expression(hat(delta) ~ "(percentile-spread, pp)"),
                  y = NULL) +
    theme_grass() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  p
}

# ---- view 4: intervals (forest of observed coefficients with CIs) -------

.plot_card_intervals <- function(x, ...) {
  check_ggplot2()
  panel <- x$panel
  if (is.null(panel) || nrow(panel) == 0L) {
    stop("The card has no `panel` data frame to plot.", call. = FALSE)
  }

  # The panel may carry ci_lower / ci_upper from compute_panel; if not we
  # widen the point by se_q_hat as a fallback (better than nothing for tests).
  has_ci <- all(c("ci_lower", "ci_upper") %in% names(panel))
  if (has_ci) {
    df <- data.frame(
      coefficient = factor(panel$coefficient, levels = rev(panel$coefficient)),
      observed_value = panel$observed_value,
      lower = panel$ci_lower,
      upper = panel$ci_upper,
      stringsAsFactors = FALSE
    )
  } else {
    se <- panel$se_q_hat %||% rep(NA_real_, nrow(panel))
    df <- data.frame(
      coefficient = factor(panel$coefficient, levels = rev(panel$coefficient)),
      observed_value = panel$observed_value,
      lower = panel$observed_value - 1.96 * se,
      upper = panel$observed_value + 1.96 * se,
      stringsAsFactors = FALSE
    )
  }

  ttl <- "Observed coefficients with 95% CIs"
  sub <- sprintf("k = %d, N = %d, pi_hat = %.2f",
                 as.integer(x$sample$k), as.integer(x$sample$N),
                 as.numeric(x$sample$pi_hat))

  ymin <- if (any(is.finite(df$lower))) min(df$lower, na.rm = TRUE) else min(df$observed_value, na.rm = TRUE)
  ymax <- if (any(is.finite(df$upper))) max(df$upper, na.rm = TRUE) else max(df$observed_value, na.rm = TRUE)
  pad  <- 0.05
  xlim <- c(min(0, ymin - pad), max(1, ymax + pad))

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = observed_value, y = coefficient)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted",
                        color = "grey50") +
    ggplot2::geom_vline(xintercept = 1, linetype = "dotted",
                        color = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper,
                                        yend = coefficient),
                          linewidth = 0.8, color = "grey25",
                          na.rm = TRUE) +
    ggplot2::geom_point(size = 3.2, color = "black") +
    ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::labs(title = ttl, subtitle = sub,
                  x = "coefficient value", y = NULL) +
    theme_grass()
  p
}

# ---- view 5: per_rater forest ------------------------------------------

.plot_card_per_rater <- function(x, ...) {
  check_ggplot2()
  if (is.null(x$per_rater)) {
    stop("per_rater view is only available when card$delta$flag == 'divergent'.",
         call. = FALSE)
  }
  pr <- x$per_rater
  if (!is.data.frame(pr) || nrow(pr) == 0L) {
    stop("per_rater view requires a non-empty data frame in card$per_rater.",
         call. = FALSE)
  }

  bound_only <- isTRUE(any(pr$bound_only))

  # Long-format: one row per (rater, parameter in {Se, Sp}).
  long_se <- data.frame(
    rater = pr$rater,
    parameter = "Se",
    estimate = pr$se_hat,
    lower = pr$se_lower,
    upper = pr$se_upper,
    stringsAsFactors = FALSE
  )
  long_sp <- data.frame(
    rater = pr$rater,
    parameter = "Sp",
    estimate = pr$sp_hat,
    lower = pr$sp_lower,
    upper = pr$sp_upper,
    stringsAsFactors = FALSE
  )
  df <- rbind(long_se, long_sp)
  df$rater <- factor(df$rater, levels = rev(unique(pr$rater)))

  ttl <- "Per-rater Se/Sp (latent-class fit)"
  sub <- if (bound_only)
    sprintf("k = 2 -> Hui-Walter bounds; CIs shown without point estimates")
    else sprintf("Bootstrap CIs; k = %d, N = %d",
                  as.integer(x$sample$k), as.integer(x$sample$N))

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = estimate, y = rater)) +
    ggplot2::facet_wrap(~ parameter, ncol = 2) +
    ggplot2::geom_vline(xintercept = c(0, 1), linetype = "dotted",
                        color = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper,
                                        yend = rater),
                          linewidth = 0.8, color = "grey25",
                          na.rm = TRUE)

  if (!bound_only) {
    p <- p + ggplot2::geom_point(size = 3.2, color = "black", na.rm = TRUE)
  }

  p <- p +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(title = ttl, subtitle = sub,
                  x = "estimate", y = NULL) +
    theme_grass()
  p
}

# ---- view 6: pairwise (k x k surface-percentile heatmap) ---------------

.plot_card_pairwise <- function(x, ...) {
  check_ggplot2()
  if (is.null(x$pairwise)) {
    stop("pairwise view is only available when card$delta$flag == 'divergent' (card$pairwise is populated automatically by grass_report() under divergent).",
         call. = FALSE)
  }
  pw <- x$pairwise
  pct <- pw$percentile_matrix
  pab <- pw$pabak_matrix
  k   <- nrow(pct)
  rn  <- rownames(pct) %||% paste0("R", seq_len(k))

  # Long-form for ggplot tile.
  rows <- rep(rn, times = k)
  cols <- rep(rn, each  = k)
  df <- data.frame(
    row_rater = factor(rows, levels = rev(rn)),
    col_rater = factor(cols, levels = rn),
    percentile = as.numeric(pct),
    pabak      = as.numeric(pab),
    stringsAsFactors = FALSE
  )
  # Cell label: surface percentile (off-diagonal); "—" on the diagonal.
  df$label <- ifelse(rows == cols, "—",
                     sprintf("%.0f%%", df$percentile))
  # Secondary label: PABAK_ij in parentheses (off-diagonal).
  df$pabak_label <- ifelse(rows == cols, "",
                           sprintf("(%.2f)", df$pabak))

  ttl <- "Pairwise reliability (divergent panel)"
  sub <- sprintf("k = %d, N = %d, pi_hat = %.2f; cell = surface percentile, parentheses = PABAK_ij",
                 as.integer(x$sample$k), as.integer(x$sample$N),
                 as.numeric(x$sample$pi_hat))

  has_viridis <- requireNamespace("viridisLite", quietly = TRUE)
  fill_scale <- if (has_viridis) {
    ggplot2::scale_fill_gradientn(
      colours = viridisLite::viridis(7),
      limits = c(0, 100), na.value = "white",
      name = "percentile (pp)")
  } else {
    ggplot2::scale_fill_gradient(
      low = "#fde0dd", high = "#7a0177",
      limits = c(0, 100), na.value = "white",
      name = "percentile (pp)")
  }

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = col_rater, y = row_rater,
                                    fill = percentile)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       size = 3.4, color = "black",
                       fontface = "bold") +
    ggplot2::geom_text(ggplot2::aes(label = pabak_label),
                       size = 2.8, color = "black",
                       nudge_y = -0.22) +
    fill_scale +
    ggplot2::coord_equal() +
    ggplot2::labs(title = ttl, subtitle = sub,
                  x = NULL, y = NULL) +
    theme_grass() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  p
}

# ---- view 7: diagnostic composite --------------------------------------

.plot_card_diagnostic <- function(x, ...) {
  check_ggplot2()
  check_patchwork()
  p_panel <- .plot_card_panel(x)
  p_therm <- .plot_card_thermometer(x)
  # Under divergent: pairwise heatmap is the most informative bottom view
  # because it directly exposes panel structure (clustering / outliers).
  # Fall back to per_rater (latent-class fit), then intervals.
  if (!is.null(x$pairwise)) {
    p_lower <- .plot_card_pairwise(x)
  } else if (!is.null(x$per_rater)) {
    p_lower <- .plot_card_per_rater(x)
  } else {
    p_lower <- .plot_card_intervals(x)
  }
  composite <- patchwork::wrap_plots(p_panel, p_therm, p_lower, ncol = 1,
                                     heights = c(2, 1, 2))
  composite
}

# ---- dispatch entry point ----------------------------------------------

#' Plot a `grass_card` Report Card object
#'
#' Six views of a Report Card returned by [grass_report()]. The default
#' `"surface"` view places the study on its DGP-calibrated reference surface;
#' the other views are for cross-coefficient comparison, the asymmetry
#' diagnostic, observed-coefficient confidence intervals, per-rater
#' Se/Sp from the latent-class fit (only valid when `card$delta$flag ==
#' "divergent"`), and a `patchwork` composite of the most useful at-a-glance
#' panels.
#'
#' @param x A `grass_card` object.
#' @param type One of:
#'   - `"surface"` (default) — 2D heatmap of E[primary metric] over (M_1, q)
#'     at the study's (k, N), with the observation pinned and dotted band
#'     contours at q in {0.5, 0.625, 0.75, 0.875, 1.0}.
#'   - `"panel"` — forest plot of all panel coefficients on the percentile
#'     axis (0-100 pp), so the cross-coefficient spread is visible.
#'   - `"thermometer"` — colored gauge for `delta_hat` with the
#'     aligned/caution/divergent thresholds shown.
#'   - `"intervals"` — forest plot of observed coefficients with 95%
#'     Wilson-logit CIs.
#'   - `"per_rater"` — forest plot of per-rater Ŝe and Ŝp (latent-class
#'     bootstrap CIs). Errors when `card$per_rater` is `NULL`.
#'   - `"pairwise"` — k x k tile heatmap of pairwise surface percentiles
#'     (cell label = percentile in pp; parenthetical label = PABAK_ij). Only
#'     available when `card$delta$flag == "divergent"` (auto-populated by
#'     [grass_report()] via [pairwise_agreement()]).
#'   - `"diagnostic"` — `patchwork` composite of "panel" + "thermometer" +
#'     "pairwise" (under divergent) / "per_rater" / "intervals".
#' @param ... Currently unused.
#'
#' @return A `ggplot` object (or a `patchwork` composite for
#'   `type = "diagnostic"`).
#' @export
plot.grass_card <- function(x,
                            type = c("surface", "panel", "thermometer",
                                     "intervals", "per_rater", "pairwise",
                                     "diagnostic"),
                            ...) {
  type <- match.arg(type)
  switch(type,
    surface     = .plot_card_surface(x, ...),
    panel       = .plot_card_panel(x, ...),
    thermometer = .plot_card_thermometer(x, ...),
    intervals   = .plot_card_intervals(x, ...),
    per_rater   = .plot_card_per_rater(x, ...),
    pairwise    = .plot_card_pairwise(x, ...),
    diagnostic  = .plot_card_diagnostic(x, ...)
  )
}
