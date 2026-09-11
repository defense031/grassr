#' Power analysis for a planned rating study
#'
#' `grass_power()` reads the calibrated reference surface forward. Given a
#' coefficient and a target value, it returns the probability that a panel
#' of quality `q` rating `N` subjects with `k` raters at observed positive
#' rate `pi_hat` produces a coefficient at or above `target`, or solves for
#' whichever one of `q`, `pi_hat`, `k`, `N`, `power` is left `NULL`. The
#' convention follows [stats::power.t.test()]: exactly one argument is
#' `NULL` and is solved for.
#'
#' The probability is a direct read of the quality sweep that
#' [position_on_surface()] returns, `P(coefficient >= target | q, design)
#' = 1 - p(q)`, with `p(q)` interpolated between the calibrated quality
#' levels. Nothing is simulated at call time.
#'
#' @section Solving for `N` or `k`:
#' The smallest value on the calibrated surface at which `power` is
#' reached. When the expected coefficient at the given quality and
#' prevalence sits below `target`, no sample size or rater count reaches
#' the requested power: larger designs concentrate the sampling
#' distribution around that expected value. The result then carries
#' `NA`, `feasible = FALSE`, and the reason.
#'
#' @section Solving for `pi_hat`:
#' The solution is the range of observed positive rates over which
#' `power` is reached, because feasibility in prevalence is an interval,
#' not a point. `solution` holds its two endpoints.
#'
#' @param metric One of `"pabak"`, `"fleiss_kappa"`, `"mean_ac1"`, `"icc"`.
#' @param target Coefficient value the study is planned to reach.
#' @param q Panel quality, the probability of a correct call on the
#'   `Se = Sp` diagonal, in `[0.55, 0.99]`.
#' @param pi_hat Observed positive rate in `[0.05, 0.95]`.
#' @param k Number of raters. Snaps to the nearest calibrated rater count
#'   (2, 3, 5, 8, 15, 25), as the surfaces do everywhere.
#' @param N Number of subjects in `[15, 1000]`.
#' @param power Probability of reaching `target`, in `(0, 1)`.
#'
#' @return An object of class `grass_power`: the six quantities with the
#'   solved one filled in, `solved` naming it, `feasible`, `reason` (when
#'   not feasible), `expected` (the median coefficient at the fixed
#'   design, when `q` is fixed), `curve` (power across the solved
#'   variable's range, the data `plot()` draws), and `notes` from the
#'   surface lookup.
#'
#' @examples
#' # What is the chance that five raters of quality 0.90 reach a "substantial"
#' # Fleiss' kappa on 200 subjects at 50% prevalence?
#' grass_power("fleiss_kappa", target = 0.61, q = 0.90, pi_hat = 0.50,
#'             k = 5, N = 200)
#'
#' # How many subjects for an 80% chance?
#' pw <- grass_power("fleiss_kappa", target = 0.61, q = 0.90, pi_hat = 0.50,
#'                   k = 5, power = 0.80)
#' pw
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(pw)
#'
#' # Over what prevalence range does the same study keep that chance?
#' grass_power("fleiss_kappa", target = 0.61, q = 0.90, k = 5, N = 200,
#'             power = 0.80)
#' @export
grass_power <- function(metric, target, q = NULL, pi_hat = NULL, k = NULL,
                        N = NULL, power = NULL) {
  allowed <- c("pabak", "fleiss_kappa", "mean_ac1", "icc")
  if (!is.character(metric) || length(metric) != 1L || !metric %in% allowed) {
    stop("`metric` must be one of: ", paste(shQuote(allowed), collapse = ", "),
         ".", call. = FALSE)
  }
  if (!is.numeric(target) || length(target) != 1L || !is.finite(target)) {
    stop("`target` must be a finite numeric scalar.", call. = FALSE)
  }
  args <- list(q = q, pi_hat = pi_hat, k = k, N = N, power = power)
  nulls <- names(args)[vapply(args, is.null, logical(1))]
  if (length(nulls) != 1L) {
    stop("Exactly one of `q`, `pi_hat`, `k`, `N`, `power` must be NULL ",
         "(the one to solve for); got ", length(nulls), ".", call. = FALSE)
  }
  solved <- nulls
  .chk <- function(x, nm, lo, hi) {
    if (is.null(x)) return(invisible())
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lo || x > hi)
      stop(sprintf("`%s` must be a single number in [%s, %s].", nm, lo, hi),
           call. = FALSE)
  }
  .chk(q, "q", .pw_q_range[1], .pw_q_range[2])
  .chk(pi_hat, "pi_hat", .pw_pi_range[1], .pw_pi_range[2])
  .chk(N, "N", .pw_n_range[1], .pw_n_range[2])
  if (!is.null(k) && (!is.numeric(k) || length(k) != 1L || k < 2 ||
                      k != as.integer(k)))
    stop("`k` must be an integer >= 2.", call. = FALSE)
  if (!is.null(power) && (!is.numeric(power) || length(power) != 1L ||
                          power <= 0 || power >= 1))
    stop("`power` must be a single number in (0, 1).", call. = FALSE)

  notes <- character()
  res <- list(metric = metric, target = target, q = q, pi_hat = pi_hat,
              k = k, N = N, power = power, solved = solved,
              solution = NA_real_, feasible = TRUE, reason = NULL,
              expected = NA_real_, curve = NULL, notes = notes)

  if (solved == "power") {
    ev <- .pw_eval(metric, target, q, pi_hat, k, N)
    res$power <- ev$power; res$solution <- ev$power; res$notes <- ev$notes
    res$expected <- .pw_expected(metric, q, pi_hat, k, N)
    res$curve <- .pw_curve(metric, target, q, pi_hat, k, N, over = "N")
    res$curve_var <- "N"
  } else if (solved == "N") {
    res$expected <- .pw_expected(metric, q, pi_hat, k, .pw_n_range[2])
    cv <- .pw_curve(metric, target, q, pi_hat, k, NULL, over = "N")
    res$curve <- cv; res$curve_var <- "N"
    ok <- which(cv$power >= power)
    if (length(ok)) {
      i <- ok[1L]
      n_hat <- if (i == 1L) cv$x[1L] else
        .pw_refine(function(n) .pw_eval(metric, target, q, pi_hat, k, n)$power - power,
                   cv$x[i - 1L], cv$x[i])
      res$N <- res$solution <- ceiling(n_hat)
    } else {
      res$N <- NA_real_; res$feasible <- FALSE
      res$reason <- .pw_reason(metric, target, q, pi_hat, res$expected, power, "N", cv)
    }
  } else if (solved == "k") {
    res$expected <- .pw_expected(metric, q, pi_hat, .pw_k_grid[length(.pw_k_grid)], N)
    cv <- .pw_curve(metric, target, q, pi_hat, NULL, N, over = "k")
    res$curve <- cv; res$curve_var <- "k"
    ok <- which(cv$power >= power)
    if (length(ok)) {
      res$k <- res$solution <- cv$x[ok[1L]]
    } else {
      res$k <- NA_real_; res$feasible <- FALSE
      res$reason <- .pw_reason(metric, target, q, pi_hat, res$expected, power, "k", cv)
    }
  } else if (solved == "q") {
    cv <- .pw_curve(metric, target, NULL, pi_hat, k, N, over = "q")
    res$curve <- cv; res$curve_var <- "q"
    f <- function(qq) .pw_eval(metric, target, qq, pi_hat, k, N)$power - power
    if (f(.pw_q_range[2]) < 0) {
      res$q <- NA_real_; res$feasible <- FALSE
      res$reason <- sprintf(
        "No calibrated panel quality (up to %.2f) reaches power %.2f for %s >= %.2f at pi_hat = %.2f, k = %d, N = %d.",
        .pw_q_range[2], power, .coef_label(metric), target, pi_hat, k, N)
    } else if (f(.pw_q_range[1]) >= 0) {
      res$q <- res$solution <- .pw_q_range[1]
      res$notes <- c(res$notes, sprintf(
        "Power %.2f is reached at the lowest calibrated quality %.2f; the solution is a floor.",
        power, .pw_q_range[1]))
    } else {
      res$q <- res$solution <- .pw_refine(f, .pw_q_range[1], .pw_q_range[2])
    }
  } else if (solved == "pi_hat") {
    cv <- .pw_curve(metric, target, q, NULL, k, N, over = "pi_hat")
    res$curve <- cv; res$curve_var <- "pi_hat"
    ok <- cv$x[cv$power >= power]
    if (length(ok)) {
      res$pi_hat <- res$solution <- range(ok)
      runs <- rle(cv$power >= power)
      if (sum(runs$values) > 1L)
        res$notes <- c(res$notes,
          "The feasible prevalence set is not one contiguous interval; `curve` holds the full profile.")
    } else {
      res$pi_hat <- NA_real_; res$feasible <- FALSE
      res$reason <- sprintf(
        "No observed positive rate on the calibrated surface (%.2f to %.2f) reaches power %.2f for %s >= %.2f at q = %.2f, k = %d, N = %d.",
        .pw_pi_range[1], .pw_pi_range[2], power, .coef_label(metric), target, q, k, N)
    }
  }
  if (is.null(res$curve_var)) res$curve_var <- solved
  # Lookup notes from one evaluation at the resolved design (k snap, N or
  # prevalence clamp, F-shape preset). Band notes describe the consistency
  # band on an observed value and do not apply to a power reading.
  if (res$feasible) {
    pi_eval <- if (length(res$pi_hat) == 2L) mean(res$pi_hat) else res$pi_hat
    ev <- .pw_eval(metric, target, res$q, pi_eval, res$k, res$N)
    keep <- ev$notes[!grepl("Consistency band|band", ev$notes)]
    res$notes <- unique(c(res$notes, keep))
  }
  class(res) <- "grass_power"
  res
}

# ---- internals -------------------------------------------------------------

.pw_q_range  <- c(0.55, 0.99)
.pw_pi_range <- c(0.05, 0.95)
.pw_n_range  <- c(15, 1000)
.pw_k_grid   <- c(2L, 3L, 5L, 8L, 15L, 25L)

# One evaluation: P(coefficient >= target | q, pi_hat, k, N) = 1 - p(q).
.pw_eval <- function(metric, target, q, pi_hat, k, N) {
  s <- suppressMessages(suppressWarnings(
    position_on_surface(obs_value = target, metric = metric,
                        pi_hat = pi_hat, k = k, N = N)))
  sw <- s$sweep
  if (is.null(sw) || !nrow(sw)) return(list(power = NA_real_, notes = s$notes))
  p <- stats::approx(sw$q, sw$p, xout = q, rule = 2)$y
  list(power = 1 - p, notes = s$notes)
}

# Median coefficient a quality-q panel produces at the design: the value c
# with P(coefficient <= c | q) = 0.5, found by inverting the sweep.
.pw_expected <- function(metric, q, pi_hat, k, N) {
  f <- function(cc) {
    sw <- suppressMessages(suppressWarnings(
      position_on_surface(obs_value = cc, metric = metric,
                          pi_hat = pi_hat, k = k, N = N)))$sweep
    if (is.null(sw)) return(NA_real_)
    stats::approx(sw$q, sw$p, xout = q, rule = 2)$y - 0.5
  }
  lo <- -0.99; hi <- 0.999
  flo <- f(lo); fhi <- f(hi)
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) return(NA_real_)
  tryCatch(stats::uniroot(f, c(lo, hi), tol = 1e-4)$root,
           error = function(e) NA_real_)
}

.pw_refine <- function(f, lo, hi) {
  tryCatch(stats::uniroot(f, c(lo, hi), tol = 1e-4)$root,
           error = function(e) hi)
}

.pw_reason <- function(metric, target, q, pi_hat, expected, power, var, cv) {
  what <- if (var == "N") "sample size (15 to 1,000)" else "rater count (2 to 25)"
  best <- cv[which.max(cv$power), ]
  reach <- sprintf("The largest power on the calibrated surface is %.2f, at %s = %s.",
                   best$power, var, format(best$x, big.mark = ","))
  if (is.finite(expected) && expected < target) {
    sprintf(paste0(
      "Expected %s at q = %.2f and pi_hat = %.2f is %.2f, below the target %.2f; ",
      "no %s reaches power %.2f. Larger designs concentrate the sampling ",
      "distribution around %.2f. %s"),
      .coef_label(metric), q, pi_hat, expected, target, what, power, expected, reach)
  } else {
    sprintf(paste0(
      "Expected %s at q = %.2f and pi_hat = %.2f is %.2f, near the target %.2f; ",
      "no %s reaches power %.2f. %s"),
      .coef_label(metric), q, pi_hat, expected, target, what, power, reach)
  }
}

# Power across one variable's range; the other four are fixed.
.pw_curve <- function(metric, target, q, pi_hat, k, N, over) {
  xs <- switch(over,
    N      = sort(unique(c(15L, 20L, 30L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L,
                          as.integer(round(exp(seq(log(15), log(1000), length.out = 40))))))),
    k      = .pw_k_grid,
    q      = seq(.pw_q_range[1], .pw_q_range[2], by = 0.01),
    pi_hat = seq(.pw_pi_range[1], .pw_pi_range[2], by = 0.01))
  pw <- vapply(xs, function(x) {
    .pw_eval(metric, target,
             q      = if (over == "q") x else q,
             pi_hat = if (over == "pi_hat") x else pi_hat,
             k      = if (over == "k") x else k,
             N      = if (over == "N") x else N)$power
  }, numeric(1))
  data.frame(x = xs, power = pw)
}

.pw_var_label <- function(v) {
  switch(v, N = "Number of subjects (N)", k = "Number of raters (k)",
         q = "Panel quality (q)", pi_hat = "Observed positive rate (pi_hat)",
         power = "Power")
}

#' @export
print.grass_power <- function(x, digits = 2, ...) {
  lab <- .coef_label(x$metric)
  cat(sprintf("\n     GRASS power analysis: %s >= %s\n\n", lab,
              formatC(x$target, digits = digits, format = "f")))
  fmt <- function(v, d = digits) {
    if (is.null(v) || all(is.na(v))) return("NA")
    if (length(v) == 2L) return(sprintf("%s to %s",
                                       formatC(v[1], digits = d, format = "f"),
                                       formatC(v[2], digits = d, format = "f")))
    if (v == round(v)) format(v, big.mark = ",") else
      formatC(v, digits = d, format = "f")
  }
  rows <- c(q = fmt(x$q), pi_hat = fmt(x$pi_hat), k = fmt(x$k, 0),
            N = fmt(x$N, 0), power = fmt(x$power))
  for (nm in names(rows)) {
    mark <- if (nm == x$solved) "  <- solved" else ""
    cat(sprintf("  %8s = %s%s\n", nm, rows[[nm]], mark))
  }
  if (is.finite(x$expected) && !is.null(x$q)) {
    cat(sprintf("\n  expected %s at this quality and prevalence: %s\n", lab,
                formatC(x$expected, digits = digits, format = "f")))
  }
  if (!x$feasible) {
    cat("\n"); cat(.wrap_note_lines(x$reason), sep = "\n")
  }
  if (length(x$notes)) {
    cat("\n  notes:\n"); for (n in x$notes) cat(.wrap_note_lines(n), sep = "\n")
  }
  cat("\n  Power is P(", lab, " >= target) on the calibrated reference surface;\n",
      "  see `plot()` for the curve over ", .pw_var_label(x$curve_var), ".\n",
      sep = "")
  invisible(x)
}

#' Plot a power curve from `grass_power()`
#'
#' Draws power across the solved variable's range, with the requested
#' power as a reference line and the solution marked when one exists.
#' When `power` itself was solved, the curve runs over `N`.
#'
#' @param x A `grass_power` object.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.grass_power <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot().", call. = FALSE)
  cv <- x$curve
  lab <- .coef_label(x$metric)
  p <- ggplot2::ggplot(cv, ggplot2::aes(x = x, y = power)) +
    ggplot2::geom_line(linewidth = 1, colour = "#1a1a1a") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(
      x = .pw_var_label(x$curve_var),
      y = sprintf("P(%s >= %.2f)", lab, x$target),
      title = sprintf("Power to reach %s >= %.2f", lab, x$target),
      subtitle = .pw_fixed_label(x)) +
    theme_grass()
  if (x$curve_var == "N") p <- p + ggplot2::scale_x_log10()
  if (x$solved != "power") {
    p <- p + ggplot2::geom_hline(yintercept = x$power, linetype = "dashed",
                                 colour = "#6b6b6b")
    if (x$feasible) {
      sol <- x$solution
      if (length(sol) == 2L) {
        p <- p + ggplot2::annotate("rect", xmin = sol[1], xmax = sol[2],
                                   ymin = 0, ymax = 1, alpha = 0.08,
                                   fill = "#377EB8")
      } else {
        xs <- if (x$curve_var == "N") ceiling(sol) else sol
        p <- p + ggplot2::annotate("point", x = xs, y = x$power, size = 3,
                                   colour = "#377EB8") +
          ggplot2::geom_vline(xintercept = xs, linetype = "dotted",
                              colour = "#377EB8")
      }
    }
  } else {
    p <- p + ggplot2::annotate("point", x = x$N, y = x$power, size = 3,
                               colour = "#377EB8")
  }
  p
}

.pw_fixed_label <- function(x) {
  parts <- character()
  if (x$curve_var != "q"      && !is.null(x$q))      parts <- c(parts, sprintf("q = %.2f", x$q))
  if (x$curve_var != "pi_hat" && !is.null(x$pi_hat) && length(x$pi_hat) == 1L)
    parts <- c(parts, sprintf("pi_hat = %.2f", x$pi_hat))
  if (x$curve_var != "k"      && !is.null(x$k))      parts <- c(parts, sprintf("k = %d", as.integer(x$k)))
  if (x$curve_var != "N"      && !is.null(x$N) && !is.na(x$N))
    parts <- c(parts, sprintf("N = %d", as.integer(x$N)))
  paste(parts, collapse = ", ")
}
