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
#
# v0.7.1 sweep redesign (2026-07-05): position_on_surface() now reports
# the POOLED percentile -- the observed value's position within the
# design's full achievable agreement range (trapezoid-weighted mixture
# over every calibrated quality level), monotone in obs_value -- instead
# of the retired nearest-q_hat-cell percentile (the sawtooth). The
# brackets below were re-anchored to the pooled convention; they still
# pin the deterministic sysdata lookup and would catch a per_rep
# ordering bug (which shifts the pooled percentile just as it shifted
# the old nearest-cell one). The narrative percentiles in the original
# §4/§5 comments referred to the retired convention.

test_that("Soares 2021 §5 anchor: PABAK pooled lookup is deterministic", {
  # Soares ED HEART concordance, k=2, N=336, pi_hat=0.71, observed
  # PABAK = 0.5654762. Pooled percentile ~73 under the v0.7.1 convention.
  res <- position_on_surface(
    obs_value = 0.5654762, metric = "pabak",
    pi_hat = 0.71, k = 2L, N = 336L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 60)
  expect_lt(res$percentile * 100, 85)
})

test_that("Klein 2018 AE §5 anchor: PABAK pooled lookup", {
  # Klein adverse-event panel, k=2, N=50, pi_hat=0.53, PABAK=0.64.
  # Pooled percentile ~78 under the v0.7.1 convention.
  res <- position_on_surface(
    obs_value = 0.64, metric = "pabak",
    pi_hat = 0.53, k = 2L, N = 50L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 65)
  expect_lt(res$percentile * 100, 90)
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

test_that("Hinde DSM-5-TR §5 anchor: PABAK pooled lookup", {
  # Hinde DSM-5-TR PTSD panel, k=10, N=200, pi_hat=0.663, PABAK=0.4556.
  # Pooled percentile ~64 under the v0.7.1 convention.
  res <- position_on_surface(
    obs_value = 0.4555556, metric = "pabak",
    pi_hat = 0.663, k = 10L, N = 200L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 50)
  expect_lt(res$percentile * 100, 75)
})

test_that("Hinde ICD-11 PTSD §5 anchor: PABAK pooled lookup", {
  # Hinde ICD-11 PTSD panel, k=10, N=200, pi_hat=0.634, PABAK=0.5324.
  # Pooled percentile ~73 under the v0.7.1 convention.
  res <- position_on_surface(
    obs_value = 0.5324444, metric = "pabak",
    pi_hat = 0.634, k = 10L, N = 200L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 60)
  expect_lt(res$percentile * 100, 85)
})

test_that("§4.1 aligned-case PABAK: pooled pct under matched DGP", {
  # §4.1 aligned: k=5, Se=Sp=0.85, pi=0.08, N=1000, seed=5.
  # Realized: pi_hat=0.212, PABAK=0.4884. Pooled percentile ~68.
  res <- position_on_surface(
    obs_value = 0.4884, metric = "pabak",
    pi_hat = 0.212, k = 5L, N = 1000L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 55)
  expect_lt(res$percentile * 100, 80)
})

test_that("§4.2 divergent-case PABAK: pooled pct with op_strong DGP", {
  # §4.2 divergent: k=5, op_strong (alternating Se/Sp asymmetry),
  # pi=0.50, N=1000, seed=6. Realized: pi_hat=0.499, PABAK=0.4868.
  # Pooled percentile ~67.
  res <- position_on_surface(
    obs_value = 0.4868, metric = "pabak",
    pi_hat = 0.499, k = 5L, N = 1000L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 55)
  expect_lt(res$percentile * 100, 80)
})

test_that("§4.2 divergent-case AC1: pooled pct with op_strong DGP", {
  # AC1 separates from kappa-family in op_strong; pooled percentile ~69.
  res <- position_on_surface(
    obs_value = 0.4937208, metric = "mean_ac1",
    pi_hat = 0.499, k = 5L, N = 1000L,
    method = "empirical"
  )
  expect_gt(res$percentile * 100, 58)
  expect_lt(res$percentile * 100, 80)
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
  # High-q observation: pooled percentile ~95 under the v0.7.1 convention.
  expect_gt(res$percentile * 100, 30)   # NOT clamped at 0
  expect_lt(res$percentile * 100, 99)   # NOT clamped at 100
})

test_that("k=2 surface index has the unified 13-point q-grid (v0.7.0)", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  k2_qs <- sort(unique(surf$index$q_true[surf$index$k == 2L]))
  expect_setequal(
    k2_qs,
    c(0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.92, 0.94, 0.95, 0.97, 0.99)
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

test_that("k=2 surface scenario count = 7436 (52 F x 13 q x 11 N)", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  expect_equal(sum(surf$index$k == 2L), 7436L)
})

test_that("k=25 surface scenario count = 7436 (52 F x 13 q x 11 N)", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  expect_equal(sum(surf$index$k == 25L), 7436L)
})

test_that("total surface = 44,616 scenarios across 6 k values (v0.7.0 densified)", {
  surf <- get("empirical_q_hat_surface",
              envir = asNamespace("grassr"), inherits = FALSE)
  expect_equal(nrow(surf$index), 44616L)
  expect_setequal(unique(surf$index$k), c(2L, 3L, 5L, 8L, 15L, 25L))
})
