# Tests for plot.grass_card. Phase 4A's grass_card constructor is being
# implemented in parallel; we build a minimal viable mock here so the tests
# can run without depending on the parallel work.

.mock_card <- function(flag = "aligned", primary = "pabak",
                       k = 5, N = 200) {
  card <- list(
    sample = list(k = k, N = N, pi_hat = 0.30, tau2_hat = NA,
                  axis = "inter"),
    coefficient = list(
      primary = primary, observed_value = 0.55,
      surface_percentile = 65, band = "Strong",
      qualifier = "moderate"
    ),
    delta = list(
      delta_hat = if (flag == "divergent") 25 else 5,
      flag = flag,
      thresholds = c(caution = 9.25, divergent = 11.75)
    ),
    panel = data.frame(
      coefficient = c("pabak", "mean_ac1", "fleiss_kappa", "krippendorff_a"),
      observed_value = c(0.55, 0.59, 0.43, 0.42),
      surface_percentile = if (flag == "divergent")
        c(35, 70, 50, 50) else c(65, 68, 64, 63),
      band = "Strong",
      qualifier = "moderate",
      band_probability_modal = 0.7,
      q_hat = 0.83,
      se_q_hat = 0.02,
      stringsAsFactors = FALSE
    ),
    per_rater = if (flag == "divergent")
      data.frame(rater = paste0("R", 1:5),
                 se_hat = rep(0.93, 5), sp_hat = rep(0.67, 5),
                 se_lower = rep(0.88, 5), se_upper = rep(0.97, 5),
                 sp_lower = rep(0.60, 5), sp_upper = rep(0.74, 5),
                 bound_only = FALSE,
                 stringsAsFactors = FALSE)
      else NULL,
    pairwise = if (flag == "divergent") {
      rn <- paste0("R", 1:5); kk <- 5L
      pab <- matrix(0.30, kk, kk); diag(pab) <- 1
      pct <- matrix(45,   kk, kk); diag(pct) <- NA_real_
      mar <- matrix(0.30, kk, kk); diag(mar) <- NA_real_
      bnd <- matrix("Moderate", kk, kk); diag(bnd) <- NA_character_
      qul <- matrix("moderate", kk, kk); diag(qul) <- NA_character_
      dimnames(pab) <- dimnames(pct) <- dimnames(mar) <-
        dimnames(bnd) <- dimnames(qul) <- list(rn, rn)
      pw <- list(
        pabak_matrix      = pab,
        percentile_matrix = pct,
        marginal_matrix   = mar,
        band_matrix       = bnd,
        qualifier_matrix  = qul,
        pooled_per_rater  = data.frame(
          rater = rn,
          se_tilde = rep(0.85, kk), sp_tilde = rep(0.75, kk),
          n_pool_pos = rep(50L, kk), n_pool_neg = rep(50L, kk),
          n_pool_excluded = rep(0L, kk),
          stringsAsFactors = FALSE),
        sample = list(k = kk, N = N, pi_hat = 0.30,
                      tau2_hat = NA_real_, axis = "inter"),
        notes = character(0L),
        call  = NULL)
      class(pw) <- c("grass_pairwise", "list")
      pw
    } else NULL,
    surface = list(q_grid = seq(0.5, 1, length.out = 501),
                   reference_curves = list(),
                   reference_type = "closed-form"),
    inputs = list(ratings_dim = c(N = N, k = k), axis = "inter",
                  delta_thresholds = c(9.25, 11.75),
                  bands = c(0.5, 0.625, 0.75, 0.875, 1.0),
                  band_labels = c("Poor", "Moderate", "Strong", "Excellent")),
    notes = character(0)
  )
  class(card) <- c("grass_card", "list")
  card
}

# ---- dispatch ---------------------------------------------------------

test_that("default plot.grass_card returns a ggplot via the surface view", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card()
  p <- plot(card)
  expect_s3_class(p, "ggplot")
})

test_that("invalid `type` errors via match.arg", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card()
  expect_error(plot(card, type = "nonsense"))
})

# ---- view: surface ----------------------------------------------------

test_that("type = 'surface' renders for primary = 'pabak'", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card(primary = "pabak")
  p <- plot(card, type = "surface")
  expect_s3_class(p, "ggplot")
})

test_that("type = 'surface' renders for primary = 'fleiss_kappa'", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card(primary = "fleiss_kappa")
  p <- plot(card, type = "surface")
  expect_s3_class(p, "ggplot")
})

test_that("type = 'surface' falls back when primary = 'icc' lacks a curve", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card(primary = "icc")
  # ICC + no fitted curve -> falls back to PABAK; should still return ggplot.
  expect_s3_class(plot(card, type = "surface"), "ggplot")
})

# ---- view: panel ------------------------------------------------------

test_that("type = 'panel' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card()
  expect_s3_class(plot(card, type = "panel"), "ggplot")
})

# ---- view: thermometer ------------------------------------------------

test_that("type = 'thermometer' returns a ggplot for both flag states", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(.mock_card(flag = "aligned"), type = "thermometer"),
                  "ggplot")
  expect_s3_class(plot(.mock_card(flag = "divergent"), type = "thermometer"),
                  "ggplot")
})

# ---- view: intervals --------------------------------------------------

test_that("type = 'intervals' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card()
  expect_s3_class(plot(card, type = "intervals"), "ggplot")
})

# ---- view: per_rater --------------------------------------------------

test_that("type = 'per_rater' errors when card$per_rater is NULL", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card(flag = "aligned")
  expect_error(plot(card, type = "per_rater"),
               "per_rater view is only available")
})

test_that("type = 'per_rater' returns a ggplot when per_rater is a data.frame", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card(flag = "divergent")
  expect_s3_class(plot(card, type = "per_rater"), "ggplot")
})

test_that("type = 'per_rater' handles bound_only = TRUE (k = 2)", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card(flag = "divergent", k = 2)
  card$per_rater$bound_only <- TRUE
  expect_s3_class(plot(card, type = "per_rater"), "ggplot")
})

# ---- view: pairwise ---------------------------------------------------

test_that("type = 'pairwise' errors when card$pairwise is NULL", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card(flag = "aligned")
  expect_error(plot(card, type = "pairwise"),
               "pairwise view is only available")
})

test_that("type = 'pairwise' returns a ggplot when card$pairwise is set", {
  skip_if_not_installed("ggplot2")
  card <- .mock_card(flag = "divergent")
  p <- plot(card, type = "pairwise")
  expect_s3_class(p, "ggplot")
})

# ---- view: diagnostic -------------------------------------------------

test_that("type = 'diagnostic' returns a patchwork object (or ggplot)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  card <- .mock_card(flag = "divergent")
  p <- plot(card, type = "diagnostic")
  # patchwork::wrap_plots returns an object inheriting from "patchwork"
  # and "gg".
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("type = 'diagnostic' substitutes intervals when per_rater is NULL", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  card <- .mock_card(flag = "aligned")
  p <- plot(card, type = "diagnostic")
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

# ---- standalone plot_surface() (prospective design) ----------------------

test_that("plot_surface returns a ggplot for each non-ICC metric", {
  skip_if_not_installed("ggplot2")
  for (m in c("pabak", "fleiss_kappa", "mean_ac1", "krippendorff_a")) {
    p <- plot_surface(m)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot_surface rejects ICC and unknown metrics with informative errors", {
  expect_error(plot_surface("icc"), "ICC surfaces require fitted")
  expect_error(plot_surface("nonsense"), "must be one of")
})

test_that("plot_surface with pi_hat alone draws a vertical reference line", {
  skip_if_not_installed("ggplot2")
  p <- plot_surface("pabak", pi_hat = 0.30, k = 5, N = 200)
  built <- ggplot2::ggplot_build(p)
  # Layers: raster, contour, vline -> at least 3
  expect_gte(length(p$layers), 3L)
  # Subtitle should mention pi_hat = 0.30
  expect_match(p$labels$subtitle, "pi_hat = 0\\.30")
  expect_match(p$labels$subtitle, "k = 5")
  expect_match(p$labels$subtitle, "N = 200")
})

test_that("plot_surface with observed pins a marker via closed-form inversion", {
  skip_if_not_installed("ggplot2")
  # PABAK = 0.62 inverts algebraically to q = (1 + sqrt(0.62))/2 ~ 0.894
  p <- plot_surface("pabak", pi_hat = 0.30, observed = 0.62, k = 5, N = 200)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "observed = 0\\.620")
  # The pin layer (geom_point) carries the data.frame with M1 / q.
  point_layer_idx <- vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"),
                            logical(1L))
  expect_true(any(point_layer_idx))
  pin_data <- p$layers[[which(point_layer_idx)[1]]]$data
  expect_equal(pin_data$M1, 0.30)
  expect_true(abs(pin_data$q - 0.5 * (1 + sqrt(0.62))) < 1e-3)
})

test_that("plot_surface errors when observed supplied without pi_hat", {
  expect_error(plot_surface("pabak", observed = 0.62),
               "requires `pi_hat`")
})

test_that("plot_surface validates pi_hat and bands", {
  expect_error(plot_surface("pabak", pi_hat = 1.5), "in \\(0, 1\\)")
  expect_error(plot_surface("pabak", pi_hat = -0.1), "in \\(0, 1\\)")
  expect_error(plot_surface("pabak", bands = c(0.5, 0.7, 0.6, 0.8, 1.0)),
               "strictly increasing")
})
