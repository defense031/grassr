# Mechanical fixes from the 2026-09-13 five-agent first-use trial.

test_that("input error names the offending column and caps the value list", {
  d <- data.frame(id = sprintf("P%03d", 1:30), r1 = rbinom(30, 1, .5), r2 = rbinom(30, 1, .5))
  expect_error(grass_report(d), "Column `id` has 30 distinct values")
  expect_error(grass_report(d), "\\.\\.\\.")
})

test_that("two-rater card exports NA delta_hat, not zero", {
  set.seed(3); Y <- cbind(rbinom(80, 1, .3), rbinom(80, 1, .3))
  d <- as.data.frame(grass_report(Y))
  expect_true(all(is.na(d$delta_hat)))
  expect_equal(unique(d$delta_flag), "not_applicable")
})

test_that("surface refuses coefficients that are not on the two-rater card", {
  expect_error(position_on_surface(0.3, "fleiss_kappa", pi_hat = .3, k = 2, N = 100),
               "not on the two-rater card")
  expect_error(position_on_surface(0.3, "icc", pi_hat = .3, k = 2, N = 100),
               "not on the two-rater card")
  expect_s3_class(position_on_surface(0.3, "pabak", pi_hat = .3, k = 2, N = 100),
                  "grass_surface_position")
})

test_that("agreement-family lookups never carry a discrete-mixture preset note", {
  for (m in c("pabak", "mean_ac1", "fleiss_kappa")) {
    r <- position_on_surface(0.4, m, pi_hat = 0.50, k = 5, N = 1000)
    expect_false(any(grepl("preset", r$notes)), info = m)
  }
})

test_that("grass_power prints the calibrated rater count when k snaps", {
  pw <- grass_power("pabak", q = .88, q0 = .80, prevalence = .2, k = 4, N = 60)
  expect_equal(pw$k_ref, 3L)
  expect_output(print(pw), "k = 4 \\(calibrated 3\\)")
  pw3 <- grass_power("pabak", q = .88, q0 = .80, prevalence = .2, k = 3, N = 60)
  expect_output(print(pw3), "k = 3\n")
})

test_that("power's k solve tolerates coefficients absent at k = 2", {
  r <- grass_power("fleiss_kappa", q = .90, q0 = .80, prevalence = .3, N = 200, power = .8)
  expect_true(is.na(r$curve$power[r$curve$x == 2]))
  expect_true(is.finite(r$k) && r$k >= 3)
})

test_that("summary withholds q_hat on a divergent card", {
  gen <- function(seed, k, N, Se, Sp, pi) {
    set.seed(seed); C <- rbinom(N, 1L, pi)
    sapply(seq_len(k), function(j) rbinom(N, 1L, ifelse(C == 1L, Se[j], 1 - Sp[j])))
  }
  Y <- gen(51L, 5L, 1000L, Se = rep(0.92, 5), Sp = c(0.92, 0.92, 0.92, 0.51, 0.51), pi = 0.30)
  card <- suppressWarnings(grass_report(Y, bootstrap_B = 0, verbose = FALSE))
  skip_if_not(identical(card$delta$flag, "divergent"))
  expect_output(print(summary(card)), "q_hat        : withheld")
})

test_that("grass_report and check_asymmetry no longer take `occasion`", {
  expect_false("occasion" %in% names(formals(grass_report)))
  expect_false("occasion" %in% names(formals(check_asymmetry)))
})
