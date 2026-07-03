# Regression tests pinning surface-percentile lookups for paper §4 / §5
# anchors against the bundled sysdata. These would have caught the
# 2026-05-10 k=2 per_rep ordering bug (which had inflated the Soares §5
# percentile from the correct ~9 to 96 by mis-pairing per-rep matrices
# with scenarios at sysdata build time).
#
# The values are tolerance-bracketed (not exact) because the underlying
# empirical reference is a sampling distribution; small (<2 pp) drift
# is acceptable across sysdata rebuilds. Hard mismatches indicate
# either a structural bug (like the per_rep ordering bug) or an
# intentional surface change that should be reflected in updated test
# expectations.

test_that("Soares 2021 §5 anchor: PABAK lookup stays in lower tail (~9 pct)", {
  # Soares ED HEART concordance, k=2, N=336, pi_hat=0.71, observed
  # PABAK = 0.5654762. Paper §5 was originally written against a
  # corrupted k=2 surface that gave 96th pct; the corrected reading
  # is ~9th pct.
  res <- position_on_surface(
    obs_value = 0.5654762, metric = "pabak",
    pi_hat = 0.71, k = 2L, N = 336L,
    method = "empirical"
  )
  expect_lt(res$percentile * 100, 20)   # NOT in upper tail
  expect_gt(res$percentile * 100, 1)    # NOT clamped
})

test_that("Klein 2018 AE §5 anchor: PABAK lookup near 50th pct", {
  # Klein adverse-event panel, k=2, N=50, pi_hat=0.53, PABAK=0.64.
  res <- position_on_surface(
    obs_value = 0.64, metric = "pabak",
    pi_hat = 0.53, k = 2L, N = 50L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 40)
  expect_lt(res$percentile * 100, 60)
})

test_that("Klein 2018 PR §5 anchor: PABAK lookup near 50th pct", {
  # Klein preventability panel, k=2, N=22, pi_hat=0.57, PABAK=0.3636.
  res <- position_on_surface(
    obs_value = 0.3636364, metric = "pabak",
    pi_hat = 0.57, k = 2L, N = 22L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 40)
  expect_lt(res$percentile * 100, 60)
})

test_that("Hinde DSM-5-TR §5 anchor: PABAK lookup near 10th pct", {
  # Hinde DSM-5-TR PTSD panel, k=10, N=200, pi_hat=0.663, PABAK=0.4556.
  res <- position_on_surface(
    obs_value = 0.4555556, metric = "pabak",
    pi_hat = 0.663, k = 10L, N = 200L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 5)
  expect_lt(res$percentile * 100, 15)
})

test_that("Hinde ICD-11 PTSD §5 anchor: PABAK lookup near 95th pct", {
  # Hinde ICD-11 PTSD panel, k=10, N=200, pi_hat=0.634, PABAK=0.5324.
  res <- position_on_surface(
    obs_value = 0.5324444, metric = "pabak",
    pi_hat = 0.634, k = 10L, N = 200L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 90)
  expect_lt(res$percentile * 100, 99)
})

test_that("§4.1 aligned-case PABAK: pct near 46 under matched DGP", {
  # §4.1 aligned: k=5, Se=Sp=0.85, pi=0.08, N=1000, seed=5.
  # Realized: pi_hat=0.212, PABAK=0.4884.
  res <- position_on_surface(
    obs_value = 0.4884, metric = "pabak",
    pi_hat = 0.212, k = 5L, N = 1000L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 40)
  expect_lt(res$percentile * 100, 55)
})

test_that("§4.2 divergent-case PABAK: pct near 40 with op_strong DGP", {
  # §4.2 divergent: k=5, op_strong (alternating Se/Sp asymmetry),
  # pi=0.50, N=1000, seed=6. Realized: pi_hat=0.499, PABAK=0.4868.
  res <- position_on_surface(
    obs_value = 0.4868, metric = "pabak",
    pi_hat = 0.499, k = 5L, N = 1000L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 30)
  expect_lt(res$percentile * 100, 50)
})

test_that("§4.2 divergent-case AC1: pct near 58 with op_strong DGP", {
  # AC1 separates from kappa-family in op_strong; this is the
  # divergent signal that drives the §4.2 narrative.
  res <- position_on_surface(
    obs_value = 0.4937208, metric = "mean_ac1",
    pi_hat = 0.499, k = 5L, N = 1000L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 50)
  expect_lt(res$percentile * 100, 65)
})

test_that("§4.3 k=2 easy case: q=0.97 augmentation prevents clamp-to-0", {
  # §4.3 high-q k=2 case: pi=0.077, PABAK=0.89 inverts to q_hat~0.972.
  # Without the q=0.97 grid point this snaps to q_true=0.99 (tight
  # distribution centered above q_hat) and clamps to pct=0. With
  # the augmentation it snaps to q_true=0.97 and reads near median.
  res <- position_on_surface(
    obs_value = 0.89, metric = "pabak",
    pi_hat = 0.077, k = 2L, N = 200L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 30)
  expect_lt(res$percentile * 100, 80)
})

test_that("k=2 surface index has the 12-point q-grid (canonical + 0.94, 0.97)", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  k2_qs <- sort(unique(surf$index$q_true[surf$index$k == 2L]))
  expect_setequal(
    k2_qs,
    c(0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
      0.94, 0.95, 0.97, 0.99)
  )
})

test_that("k=25 surface index has the 13-point q-grid (canonical + 0.92, 0.94, 0.97)", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  k25_qs <- sort(unique(surf$index$q_true[surf$index$k == 25L]))
  expect_setequal(
    k25_qs,
    c(0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
      0.92, 0.94, 0.95, 0.97, 0.99)
  )
})

test_that("Row 6 (k=25 PABAK=0.78 pi=0.30 N=1000) unclamps post-q094 aug", {
  # tab:card_summary row 6 panel. q_hat = 0.9416 falls in the
  # 0.92-0.95 midpoint gap before q=0.94 augmentation; clamps to 0.
  # After q=0.94 aug it snaps to q=0.94 cohort and reads centered
  # in the upper portion of that cohort's distribution (since
  # q_hat is slightly above 0.94).
  res <- position_on_surface(
    obs_value = 0.78, metric = "pabak",
    pi_hat = 0.30, k = 25L, N = 1000L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 30)   # NOT clamped at 0
  expect_lt(res$percentile * 100, 99)   # NOT clamped at 100
})

test_that("k=2 surface F-key set matches canonical k>=3 set (52 keys)", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  k2_F <- sort(unique(surf$index$F_key[surf$index$k == 2L]))
  k3_F <- sort(unique(surf$index$F_key[surf$index$k == 3L]))
  expect_equal(length(k2_F), 52L)
  expect_setequal(k2_F, k3_F)
})

test_that("k=2 surface scenario count = 1872 (52 F x 12 q x 3 N)", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  expect_equal(sum(surf$index$k == 2L), 1872L)
})

test_that("k=25 surface scenario count = 2028 (52 F x 13 q x 3 N)", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  expect_equal(sum(surf$index$k == 25L), 2028L)
})

test_that("total surface = 10,140 scenarios across 6 k values", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  expect_equal(nrow(surf$index), 10140L)
  expect_setequal(unique(surf$index$k), c(2L, 3L, 5L, 8L, 15L, 25L))
})
