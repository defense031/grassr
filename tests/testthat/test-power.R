# grass_power(): the reference surface read forward (2026-09-11).
# Anchors are the k = 5, N = 200, q = 0.90 values recorded in
# paper1_2_merged/decisions.md (2026-09-10 entry).

test_that("exactly one of q/pi_hat/k/N/power must be NULL", {
  expect_error(grass_power("pabak", 0.61, q = 0.9, pi_hat = 0.5, k = 5, N = 200, power = 0.8),
               "Exactly one")
  expect_error(grass_power("pabak", 0.61, q = 0.9, pi_hat = 0.5),
               "Exactly one")
  expect_error(grass_power("nope", 0.61, q = 0.9, pi_hat = 0.5, k = 5, N = 200),
               "`metric`")
})

test_that("solving for power reproduces the recorded anchors", {
  a <- grass_power("fleiss_kappa", 0.61, q = 0.90, pi_hat = 0.50, k = 5, N = 200)
  expect_s3_class(a, "grass_power")
  expect_equal(a$solved, "power")
  expect_equal(a$power, 0.80, tolerance = 0.03)
  b <- grass_power("fleiss_kappa", 0.61, q = 0.90, pi_hat = 0.85, k = 5, N = 200)
  expect_lt(b$power, 0.01)
  p1 <- grass_power("pabak", 0.61, q = 0.90, pi_hat = 0.50, k = 5, N = 200)$power
  p2 <- grass_power("pabak", 0.61, q = 0.90, pi_hat = 0.85, k = 5, N = 200)$power
  expect_lt(abs(p1 - p2), 0.05)   # PABAK flat across prevalence
  expect_true(is.finite(a$expected) && a$expected > 0.61)
})

test_that("solving for N returns the smallest feasible N, or NA with a reason", {
  s <- grass_power("fleiss_kappa", 0.61, q = 0.90, pi_hat = 0.50, k = 5, power = 0.80)
  expect_equal(s$solved, "N")
  expect_true(s$feasible)
  expect_true(is.finite(s$N) && s$N >= 15 && s$N <= 1000)
  at <- grass_power("fleiss_kappa", 0.61, q = 0.90, pi_hat = 0.50, k = 5, N = s$N)$power
  expect_gte(at, 0.80 - 0.02)
  expect_true(is.data.frame(s$curve) && all(c("x", "power") %in% names(s$curve)))

  f <- grass_power("fleiss_kappa", 0.61, q = 0.90, pi_hat = 0.70, k = 5, power = 0.80)
  expect_false(f$feasible)
  expect_true(is.na(f$N))
  expect_match(f$reason, "below the target")
  expect_true(is.finite(f$expected) && f$expected < 0.61)
})

test_that("solving for k, q, and pi_hat behave", {
  kk <- grass_power("fleiss_kappa", 0.61, q = 0.90, pi_hat = 0.50, N = 200, power = 0.80)
  expect_true(kk$feasible); expect_true(kk$k %in% c(2, 3, 5, 8, 15, 25))

  qq <- grass_power("fleiss_kappa", 0.61, pi_hat = 0.85, k = 5, N = 200, power = 0.50)
  expect_true(qq$feasible)
  expect_true(qq$q > 0.90 && qq$q <= 0.99)
  expect_equal(grass_power("fleiss_kappa", 0.61, q = qq$q, pi_hat = 0.85, k = 5, N = 200)$power,
               0.50, tolerance = 0.02)

  pp <- grass_power("fleiss_kappa", 0.61, q = 0.90, k = 5, N = 200, power = 0.50)
  expect_true(pp$feasible)
  expect_length(pp$solution, 2L)
  expect_true(pp$solution[1] <= 0.50 && pp$solution[2] >= 0.50)
  expect_lt(pp$solution[2], 0.85)
})

test_that("print and plot methods run", {
  pw <- grass_power("pabak", 0.61, q = 0.90, pi_hat = 0.50, k = 5, power = 0.80)
  expect_output(print(pw), "solved")
  expect_output(print(pw), "power")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(pw), "ggplot")
  expect_s3_class(plot(grass_power("pabak", 0.61, q = 0.90, k = 5, N = 200, power = 0.5)), "ggplot")
  expect_s3_class(plot(grass_power("pabak", 0.61, q = 0.90, pi_hat = 0.5, k = 5, N = 200)), "ggplot")
})
