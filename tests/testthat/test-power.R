# grass_power(): the reference surface read forward (2026-09-11).
# Primary target is q0, the lower edge of the quality resolution (ratified 2026-09-11:
# a fixed coefficient value is the external-threshold mode only).
# Value-mode anchors are the k = 5, N = 200, q = 0.90 numbers in
# paper1_2_merged/decisions.md (2026-09-10 entry).

test_that("argument contract: one target, one unknown", {
  expect_error(grass_power("pabak", q = 0.9, pi_hat = 0.5, k = 5, N = 200), "exactly one of `q0`")
  expect_error(grass_power("pabak", q0 = 0.8, target = 0.6, q = 0.9, pi_hat = 0.5, k = 5, N = 200),
               "exactly one of `q0`")
  expect_error(grass_power("pabak", q0 = 0.8, q = 0.9, pi_hat = 0.5, k = 5, N = 200, power = 0.8),
               "Exactly one")
  expect_error(grass_power("pabak", q0 = 0.8, q = 0.9, pi_hat = 0.5), "Exactly one")
  expect_error(grass_power("pabak", q0 = 0.9, q = 0.85, pi_hat = 0.5, k = 5, N = 200), "must exceed")
  expect_error(grass_power("nope", q0 = 0.8, q = 0.9, pi_hat = 0.5, k = 5, N = 200), "`metric`")
})

test_that("quality mode: power rises with N and reproduces the recorded anchors", {
  # screening study k=3, pi=0.10: N for 0.80 power is ~120 (record 2026-09-11)
  s <- grass_power("fleiss_kappa", q = 0.90, q0 = 0.80, pi_hat = 0.10, k = 3, power = 0.80)
  expect_equal(s$mode, "quality"); expect_equal(s$solved, "N")
  expect_true(s$feasible)
  expect_true(s$N >= 90 && s$N <= 150)
  expect_true(all(diff(s$curve$power) >= -0.02))   # monotone in N (MC slack)
  at <- grass_power("fleiss_kappa", q = 0.90, q0 = 0.80, pi_hat = 0.10, k = 3, N = s$N)$power
  expect_gte(at, 0.78)
  # modal design k=5, pi=0.50: ~25 subjects
  m <- grass_power("fleiss_kappa", q = 0.90, q0 = 0.80, pi_hat = 0.50, k = 5, power = 0.80)
  expect_true(m$feasible && m$N <= 40)
  # a hopeless gap: q barely above q0 at a small design -> infeasible with a reason
  h <- grass_power("fleiss_kappa", q = 0.81, q0 = 0.80, pi_hat = 0.50, k = 3, power = 0.90)
  expect_false(h$feasible); expect_true(is.na(h$N)); expect_match(h$reason, "largest power")
})

test_that("quality mode: k, q, and pi_hat solve", {
  kk <- grass_power("fleiss_kappa", q = 0.90, q0 = 0.80, pi_hat = 0.10, N = 80, power = 0.80)
  expect_true(kk$feasible); expect_true(kk$k %in% c(2, 3, 5, 8, 15, 25))
  qq <- grass_power("fleiss_kappa", q0 = 0.80, pi_hat = 0.10, k = 3, N = 80, power = 0.80)
  expect_true(qq$feasible); expect_true(qq$q > 0.80 && qq$q <= 0.99)
  expect_equal(grass_power("fleiss_kappa", q = qq$q, q0 = 0.80, pi_hat = 0.10, k = 3, N = 80)$power,
               0.80, tolerance = 0.03)
  pp <- grass_power("fleiss_kappa", q = 0.90, q0 = 0.80, k = 5, N = 100, power = 0.80)
  expect_true(pp$feasible); expect_length(pp$solution, 2L)
  expect_true(pp$solution[1] <= 0.50 && pp$solution[2] >= 0.50)
})

test_that("value mode reproduces the 2026-09-10 anchors, including the falling curve", {
  a <- grass_power("fleiss_kappa", target = 0.61, q = 0.90, pi_hat = 0.50, k = 5, N = 200)
  expect_equal(a$mode, "value"); expect_equal(a$power, 0.80, tolerance = 0.03)
  expect_true(is.finite(a$expected) && a$expected > 0.61)
  b <- grass_power("fleiss_kappa", target = 0.61, q = 0.90, pi_hat = 0.85, k = 5, N = 200)
  expect_lt(b$power, 0.01)
  p1 <- grass_power("pabak", target = 0.61, q = 0.90, pi_hat = 0.50, k = 5, N = 200)$power
  p2 <- grass_power("pabak", target = 0.61, q = 0.90, pi_hat = 0.85, k = 5, N = 200)$power
  expect_lt(abs(p1 - p2), 0.05)
  f <- grass_power("fleiss_kappa", target = 0.61, q = 0.90, pi_hat = 0.70, k = 5, power = 0.80)
  expect_false(f$feasible); expect_true(is.na(f$N)); expect_match(f$reason, "below the target")
  expect_lt(f$curve$power[nrow(f$curve)], f$curve$power[1])   # the fixed-value pathology
})

test_that("print and plot methods run in both modes", {
  pw <- grass_power("pabak", q = 0.90, q0 = 0.80, pi_hat = 0.50, k = 5, power = 0.80)
  expect_output(print(pw), "resolve panel quality")
  expect_output(print(pw), "solved")
  vv <- grass_power("pabak", target = 0.61, q = 0.90, pi_hat = 0.50, k = 5, N = 200)
  expect_output(print(vv), "fixed value")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(pw), "ggplot")
  expect_s3_class(plot(vv), "ggplot")
  expect_s3_class(plot(grass_power("pabak", q = 0.9, q0 = 0.8, k = 5, N = 100, power = 0.8)), "ggplot")
})
