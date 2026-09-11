#' Plan a rating study: what a design can show
#'
#' `grass_power()` plans a study's size. Assume a prevalence and a rough
#' rater quality, and it says what a design of a given size can show
#' about the panel, and how many more subjects or raters it would take to
#' show more. The result describes the design, how much a study of that
#' size can learn about its raters.
#'
#' A study sized to tell a panel of quality `q` from one of quality `q0`
#' returns a 95% consistency band on quality narrow enough to separate
#' the two. `q0` sets the resolution of the plan. A panel weaker than
#' assumed returns a band that sits lower and is about as wide, so the
#' study reports the quality it finds at the precision it was planned
#' for. The convention follows [stats::power.t.test()]. Fix four of `q`,
#' `pi_hat`, `k`, `N`, `power`, leave one `NULL`, and the function solves
#' for it.
#'
#' Whether more subjects or more raters raises power depends on
#' prevalence. At balanced prevalence a few more raters do the work of
#' subjects. At a rare or very common finding they do not, because a
#' small sample holds only a handful of the minority class, and the
#' answer is more subjects.
#'
#' The function also accepts `target`, a fixed coefficient value, for a
#' threshold imposed from outside (a journal's or regulator's band). A
#' fixed coefficient value means something different at each design.
#' When the value a panel of quality `q` produces at the design sits
#' below `target`, power falls with `N`. Give `q0` or `target`, not both.
#'
#' @section How it is computed:
#' Every quantity is a direct read of the quality sweep that
#' [position_on_surface()] returns, `p(q) = P(coefficient <= c | q, design)`.
#' Nothing is simulated at call time. For `q0`, by test inversion the lower
#' end of the band is above `q0` exactly when the observed coefficient exceeds
#' the 97.5th percentile of the `q0` distribution, `c0`, so
#' `power = 1 - p_q(c0)`. For `target`, `power = 1 - p_q(target)`.
#'
#' @section Solving for `N` or `k`:
#' The smallest value on the calibrated surface at which `power` is
#' reached. When none reaches it the result carries `NA`,
#' `feasible = FALSE`, and the reason, including the best power any
#' design on the surface reaches.
#'
#' @section Solving for `pi_hat`:
#' The range of observed positive rates over which `power` is reached;
#' `solution` holds its two endpoints.
#'
#' @param metric One of `"pabak"`, `"fleiss_kappa"`, `"mean_ac1"`, `"icc"`.
#' @param q Panel quality, the probability of a correct call on the
#'   `Se = Sp` diagonal, in `[0.55, 0.99]`.
#' @param q0 The lower edge of the quality resolution the study is planned
#'   for, in `[0.55, 0.99]`: the study is sized to tell a panel of quality
#'   `q` from one of quality `q0`. Give `q0` or `target`, not both.
#' @param target A fixed coefficient value to reach. Give `q0` or `target`,
#'   not both.
#' @param pi_hat Observed positive rate in `[0.05, 0.95]`.
#' @param k Number of raters. Snaps to the nearest calibrated rater count
#'   (2, 3, 5, 8, 15, 25), as the surfaces do everywhere.
#' @param N Number of subjects in `[15, 1000]`.
#' @param power Probability, in `(0, 1)`.
#'
#' @return An object of class `grass_power`: the design quantities with the
#'   solved one filled in, `solved` naming it, `mode` (`"quality"` or
#'   `"value"`), `feasible`, `reason` (when not feasible), `expected` (the
#'   median coefficient a panel of quality `q` produces at the design, in
#'   `"value"` mode), `curve` (power across the solved variable's range,
#'   the data `plot()` draws), and `notes` from the surface lookup.
#'
#' @examples
#' # Raters assumed near quality 0.90, a 10% positive rate, three raters:
#' # how many subjects to tell a 0.90 panel from a 0.80 one, 80% power?
#' pw <- grass_power("fleiss_kappa", q = 0.90, q0 = 0.80, pi_hat = 0.10,
#'                   k = 3, power = 0.80)
#' pw
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(pw)
#'
#' # Over what positive rates does a fixed design keep that power?
#' grass_power("fleiss_kappa", q = 0.90, q0 = 0.80, k = 5, N = 100,
#'             power = 0.80)
#'
#' # Against a fixed coefficient value imposed from outside:
#' grass_power("fleiss_kappa", target = 0.61, q = 0.90, pi_hat = 0.50,
#'             k = 5, N = 200)
#' @export
grass_power <- function(metric, q = NULL, q0 = NULL, target = NULL,
                        pi_hat = NULL, k = NULL, N = NULL, power = NULL) {
  allowed <- c("pabak", "fleiss_kappa", "mean_ac1", "icc")
  if (!is.character(metric) || length(metric) != 1L || !metric %in% allowed) {
    stop("`metric` must be one of: ", paste(shQuote(allowed), collapse = ", "),
         ".", call. = FALSE)
  }
  if (is.null(q0) == is.null(target)) {
    stop("Give exactly one of `q0` (the lower edge of the quality resolution) or ",
         "`target` (a fixed coefficient value).", call. = FALSE)
  }
  mode <- if (is.null(target)) "quality" else "value"
  .chk <- function(x, nm, lo, hi) {
    if (is.null(x)) return(invisible())
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lo || x > hi)
      stop(sprintf("`%s` must be a single number in [%s, %s].", nm, lo, hi),
           call. = FALSE)
  }
  .chk(q0, "q0", .pw_q_range[1], .pw_q_range[2])
  if (!is.null(target) && (!is.numeric(target) || length(target) != 1L ||
                           !is.finite(target)))
    stop("`target` must be a finite numeric scalar.", call. = FALSE)
  args <- list(q = q, pi_hat = pi_hat, k = k, N = N, power = power)
  nulls <- names(args)[vapply(args, is.null, logical(1))]
  if (length(nulls) != 1L) {
    stop("Exactly one of `q`, `pi_hat`, `k`, `N`, `power` must be NULL ",
         "(the one to solve for); got ", length(nulls), ".", call. = FALSE)
  }
  solved <- nulls
  .chk(q, "q", .pw_q_range[1], .pw_q_range[2])
  .chk(pi_hat, "pi_hat", .pw_pi_range[1], .pw_pi_range[2])
  .chk(N, "N", .pw_n_range[1], .pw_n_range[2])
  if (!is.null(k) && (!is.numeric(k) || length(k) != 1L || k < 2 ||
                      k != as.integer(k)))
    stop("`k` must be an integer >= 2.", call. = FALSE)
  if (!is.null(power) && (!is.numeric(power) || length(power) != 1L ||
                          power <= 0 || power >= 1))
    stop("`power` must be a single number in (0, 1).", call. = FALSE)
  if (mode == "quality" && !is.null(q) && q <= q0)
    stop("`q` (the assumed panel quality) must exceed `q0` (the lower edge of the resolution).",
         call. = FALSE)

  tg <- list(mode = mode, q0 = q0, target = target)
  res <- list(metric = metric, mode = mode, q0 = q0, target = target,
              q = q, pi_hat = pi_hat, k = k, N = N, power = power,
              solved = solved, solution = NA_real_, feasible = TRUE,
              reason = NULL, expected = NA_real_, curve = NULL,
              curve_var = NULL, notes = character())

  if (solved == "power") {
    res$power <- res$solution <- .pw_eval(metric, tg, q, pi_hat, k, N)$power
    res$curve <- .pw_curve(metric, tg, q, pi_hat, k, N, over = "N")
    res$curve_var <- "N"
  } else if (solved %in% c("N", "k")) {
    cv <- .pw_curve(metric, tg, q, pi_hat,
                    if (solved == "k") NULL else k,
                    if (solved == "N") NULL else N, over = solved)
    res$curve <- cv; res$curve_var <- solved
    ok <- which(cv$power >= power)
    if (length(ok)) {
      i <- ok[1L]
      if (solved == "N") {
        n_hat <- if (i == 1L) cv$x[1L] else
          .pw_refine(function(n) .pw_eval(metric, tg, q, pi_hat, k, n)$power - power,
                     cv$x[i - 1L], cv$x[i])
        res$N <- res$solution <- ceiling(n_hat)
      } else {
        res$k <- res$solution <- cv$x[i]
      }
    } else {
      res[[solved]] <- NA_real_; res$feasible <- FALSE
      res$reason <- .pw_reason(metric, tg, q, pi_hat, k, N, power, solved, cv)
    }
  } else if (solved == "q") {
    lo <- if (mode == "quality") max(.pw_q_range[1], q0 + 0.005) else .pw_q_range[1]
    cv <- .pw_curve(metric, tg, NULL, pi_hat, k, N, over = "q", q_lo = lo)
    res$curve <- cv; res$curve_var <- "q"
    f <- function(qq) .pw_eval(metric, tg, qq, pi_hat, k, N)$power - power
    if (f(.pw_q_range[2]) < 0) {
      res$q <- NA_real_; res$feasible <- FALSE
      res$reason <- sprintf(
        "No calibrated panel quality (up to %.2f) reaches power %.2f %s at pi_hat = %.2f, k = %d, N = %d.",
        .pw_q_range[2], power, .pw_goal(metric, tg), pi_hat, k, N)
    } else if (f(lo) >= 0) {
      res$q <- res$solution <- lo
      res$notes <- c(res$notes, sprintf(
        "Power %.2f is reached at the lowest admissible quality %.3f; the solution is a floor.",
        power, lo))
    } else {
      res$q <- res$solution <- .pw_refine(f, lo, .pw_q_range[2])
    }
  } else if (solved == "pi_hat") {
    cv <- .pw_curve(metric, tg, q, NULL, k, N, over = "pi_hat")
    res$curve <- cv; res$curve_var <- "pi_hat"
    hit <- !is.na(cv$power) & cv$power >= power
    ok <- cv$x[hit]
    if (length(ok)) {
      res$pi_hat <- res$solution <- range(ok)
      if (sum(rle(hit)$values) > 1L)
        res$notes <- c(res$notes,
          "The feasible prevalence set is not one contiguous interval; `curve` holds the full profile.")
    } else {
      res$pi_hat <- NA_real_; res$feasible <- FALSE
      res$reason <- sprintf(
        "No observed positive rate on the calibrated surface (%.2f to %.2f) reaches power %.2f %s at q = %.2f, k = %d, N = %d.",
        .pw_pi_range[1], .pw_pi_range[2], power, .pw_goal(metric, tg), q, k, N)
    }
  }
  # Value mode: the median coefficient a quality-q panel produces at the
  # design (at the largest calibrated N or k when that axis is unsolved).
  if (mode == "value" && !is.null(res$q) && !is.na(res$q) &&
      length(res$pi_hat) == 1L && !is.na(res$pi_hat)) {
    N_e <- if (is.null(res$N) || is.na(res$N)) .pw_n_range[2] else res$N
    k_e <- if (is.null(res$k) || is.na(res$k)) .pw_k_grid[length(.pw_k_grid)] else res$k
    res$expected <- .pw_expected(metric, res$q, res$pi_hat, k_e, N_e)
  }
  # Lookup notes from one evaluation at the resolved design (k snap, N or
  # prevalence clamp, F-shape preset). Band notes describe the consistency
  # band on an observed value and do not apply to a power reading.
  if (res$feasible) {
    pi_eval <- if (length(res$pi_hat) == 2L) mean(res$pi_hat) else res$pi_hat
    ev <- .pw_eval(metric, tg, res$q, pi_eval, res$k, res$N)
    keep <- ev$notes[!grepl("band", ev$notes)]
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

.pw_sweep <- function(metric, cc, pi_hat, k, N) {
  s <- suppressMessages(suppressWarnings(
    position_on_surface(obs_value = cc, metric = metric,
                        pi_hat = pi_hat, k = k, N = N)))
  list(sweep = s$sweep, notes = s$notes)
}

# p_q(c): P(coefficient <= c | quality q, design), q interpolated between
# the calibrated levels.
.pw_p <- function(sw, qq) stats::approx(sw$q, sw$p, xout = qq, rule = 2)$y

# The coefficient value c with p_{qq}(c) = prob at the design.
.pw_quantile <- function(metric, qq, prob, pi_hat, k, N) {
  f <- function(cc) {
    sw <- .pw_sweep(metric, cc, pi_hat, k, N)$sweep
    if (is.null(sw)) return(NA_real_)
    .pw_p(sw, qq) - prob
  }
  flo <- f(-0.999); fhi <- f(0.9999)
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) return(NA_real_)
  tryCatch(stats::uniroot(f, c(-0.999, 0.9999), tol = 1e-4)$root,
           error = function(e) NA_real_)
}

# One evaluation of power at a design.
.pw_eval <- function(metric, tg, q, pi_hat, k, N) {
  cc <- if (tg$mode == "value") tg$target else
    .pw_quantile(metric, tg$q0, 0.975, pi_hat, k, N)
  if (is.na(cc)) return(list(power = NA_real_, notes = character()))
  s <- .pw_sweep(metric, cc, pi_hat, k, N)
  if (is.null(s$sweep) || !nrow(s$sweep)) return(list(power = NA_real_, notes = s$notes))
  list(power = 1 - .pw_p(s$sweep, q), notes = s$notes)
}

.pw_expected <- function(metric, q, pi_hat, k, N) .pw_quantile(metric, q, 0.5, pi_hat, k, N)

.pw_refine <- function(f, lo, hi) {
  tryCatch(stats::uniroot(f, c(lo, hi), tol = 1e-4)$root,
           error = function(e) hi)
}

.pw_goal <- function(metric, tg) {
  if (tg$mode == "quality") sprintf("to resolve panel quality from %.2f", tg$q0)
  else sprintf("for %s >= %.2f", .coef_label(metric), tg$target)
}

.pw_reason <- function(metric, tg, q, pi_hat, k, N, power, var, cv) {
  what <- if (var == "N") "sample size (15 to 1,000)" else "rater count (2 to 25)"
  best <- cv[which.max(cv$power), ]
  reach <- sprintf("The largest power on the calibrated surface is %.2f, at %s = %s.",
                   best$power, var, format(best$x, big.mark = ","))
  if (tg$mode == "quality") {
    return(sprintf("No %s reaches power %.2f to separate panel quality %.2f from %.2f at pi_hat = %.2f. %s",
                   what, power, q, tg$q0, pi_hat, reach))
  }
  expected <- .pw_expected(metric, q, pi_hat,
                           if (var == "k") .pw_k_grid[length(.pw_k_grid)] else k,
                           if (var == "N") .pw_n_range[2] else N)
  if (is.finite(expected) && expected < tg$target) {
    sprintf(paste0(
      "Expected %s at q = %.2f and pi_hat = %.2f is %.2f, below the target %.2f; ",
      "no %s reaches power %.2f. Larger designs concentrate the sampling ",
      "distribution around %.2f. %s"),
      .coef_label(metric), q, pi_hat, expected, tg$target, what, power, expected, reach)
  } else {
    sprintf(paste0(
      "Expected %s at q = %.2f and pi_hat = %.2f is %.2f, near the target %.2f; ",
      "no %s reaches power %.2f. %s"),
      .coef_label(metric), q, pi_hat, expected, tg$target, what, power, reach)
  }
}

# Power across one variable's range; the other four are fixed.
.pw_curve <- function(metric, tg, q, pi_hat, k, N, over, q_lo = .pw_q_range[1]) {
  xs <- switch(over,
    N      = sort(unique(c(15L, 20L, 30L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L,
                          as.integer(round(exp(seq(log(15), log(1000), length.out = 40))))))),
    k      = .pw_k_grid,
    q      = seq(q_lo, .pw_q_range[2], length.out = 45),
    pi_hat = seq(.pw_pi_range[1], .pw_pi_range[2], by = 0.01))
  pw <- vapply(xs, function(x) {
    .pw_eval(metric, tg,
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

.pw_ylab <- function(x) {
  if (x$mode == "quality") sprintf("P(band separates quality from %.2f)", x$q0)
  else sprintf("P(%s >= %.2f)", .coef_label(x$metric), x$target)
}

#' @export
print.grass_power <- function(x, digits = 2, ...) {
  lab <- .coef_label(x$metric)
  hdr <- if (x$mode == "quality")
    sprintf("resolve panel quality %s from %s", if (is.null(x$q) || is.na(x$q)) "(solved)" else formatC(x$q, digits = digits, format = "f"), formatC(x$q0, digits = digits, format = "f"))
  else
    sprintf("reach %s >= %s (fixed value)", lab, formatC(x$target, digits = digits, format = "f"))
  cat(sprintf("\n     GRASS power analysis: %s\n", hdr))
  cat(sprintf("     coefficient: %s\n\n", lab))
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
  if (x$mode == "value" && is.finite(x$expected)) {
    cat(sprintf("\n  expected %s at this quality and prevalence: %s\n", lab,
                formatC(x$expected, digits = digits, format = "f")))
  }
  if (!x$feasible) {
    cat("\n"); cat(.wrap_note_lines(x$reason), sep = "\n")
  }
  if (length(x$notes)) {
    cat("\n  notes:\n"); for (n in x$notes) cat(.wrap_note_lines(n), sep = "\n")
  }
  cat("\n")
  if (x$mode == "quality") {
    cat(.wrap_note_lines(sprintf(
      "Power is the probability that a study of this size tells panel quality %s from %.2f.",
      if (is.null(x$q) || is.na(x$q)) "the solved quality" else formatC(x$q, digits = digits, format = "f"), x$q0),
      indent = "  "), sep = "\n")
  } else {
    cat(.wrap_note_lines(sprintf(
      "Power is P(%s >= %.2f) at this design. A fixed coefficient value means something different at each design; `q0 =` sizes the study on panel quality instead.",
      lab, x$target), indent = "  "), sep = "\n")
  }
  cat("  See `plot()` for the curve over ", .pw_var_label(x$curve_var), ".\n", sep = "")
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
  ttl <- if (x$mode == "quality")
    sprintf("Power to resolve panel quality %s from %.2f", if (is.null(x$q) || is.na(x$q)) "(solved)" else sprintf("%.2f", x$q), x$q0)
  else sprintf("Power to reach %s >= %.2f", .coef_label(x$metric), x$target)
  p <- ggplot2::ggplot(cv, ggplot2::aes(x = x, y = power)) +
    ggplot2::geom_line(linewidth = 1, colour = "#1a1a1a") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(x = .pw_var_label(x$curve_var), y = .pw_ylab(x),
                  title = ttl, subtitle = .pw_fixed_label(x)) +
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
        p <- p + ggplot2::annotate("point", x = sol, y = x$power, size = 3,
                                   colour = "#377EB8") +
          ggplot2::geom_vline(xintercept = sol, linetype = "dotted",
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
  parts <- sprintf("%s", .coef_label(x$metric))
  if (x$curve_var != "q"      && !is.null(x$q) && !is.na(x$q))
    parts <- c(parts, sprintf("q = %.2f", x$q))
  if (x$curve_var != "pi_hat" && !is.null(x$pi_hat) && length(x$pi_hat) == 1L && !is.na(x$pi_hat))
    parts <- c(parts, sprintf("pi_hat = %.2f", x$pi_hat))
  if (x$curve_var != "k"      && !is.null(x$k) && !is.na(x$k))
    parts <- c(parts, sprintf("k = %d", as.integer(x$k)))
  if (x$curve_var != "N"      && !is.null(x$N) && !is.na(x$N))
    parts <- c(parts, sprintf("N = %d", as.integer(x$N)))
  paste(parts, collapse = ", ")
}
