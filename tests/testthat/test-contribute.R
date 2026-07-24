# Contributor pipeline: manifest integrity, seeded reproducibility, and
# bundle verification. Real blocks run 50,000 draws; tests override to a
# handful so the whole file stays fast.

test_that("bundled manifest is intact and seeds are disjoint from shipped programs", {
  m <- grass_calibration_manifest()
  expect_equal(nrow(m), 2166L)
  expect_setequal(unique(m$program), c("tail_topup", "prev_strata", "lattice_k"))
  expect_equal(sum(m$program == "tail_topup"), 21L)
  expect_equal(sum(m$program == "prev_strata"), 1925L)
  expect_equal(sum(m$program == "lattice_k"), 220L)
  # explicit per-block seeds, one per block, all above the shipped
  # master-seed range (20260706 + cell_id and 20260707 offsets)
  expect_equal(m$seed, 20300000L + m$block_id)
  expect_false(any(duplicated(m$seed)))
  expect_true(all(m$seed > 20270000L))
  # every block is a k >= 3 cell: delta_hat is not_applicable at k = 2
  expect_true(all(m$k >= 3L))
  # prev_strata blocks pin one calibration prevalence; pooled blocks are NA
  expect_true(all(is.na(m$prev[m$program != "prev_strata"])))
  expect_true(all(m$prev[m$program == "prev_strata"] %in%
                    c(0.05, 0.20, 0.50, 0.80, 0.95)))
})

test_that("manifest program filter works", {
  m <- grass_calibration_manifest("tail_topup")
  expect_equal(unique(m$program), "tail_topup")
  expect_equal(nrow(m), 21L)
})

test_that("dry_run returns a plan inside the budget without writing anything", {
  d <- file.path(tempdir(), "contrib_dryrun_test")
  plan <- suppressMessages(
    grass_contribute(dir = d, hours = 2, dry_run = TRUE))
  expect_s3_class(plan, "data.frame")
  expect_true(nrow(plan) >= 1L)
  expect_true(all(c("block_id", "seed", "est_secs") %in% names(plan)))
  expect_false(dir.exists(d))
})

test_that("a block run is bit-reproducible from its manifest seed", {
  d1 <- file.path(tempdir(), "contrib_repro_1")
  d2 <- file.path(tempdir(), "contrib_repro_2")
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  suppressMessages(grass_contribute(dir = d1, blocks = 30L, draws = 25L))
  suppressMessages(grass_contribute(dir = d2, blocks = 30L, draws = 25L))
  b1 <- readRDS(file.path(d1, "block_0030.rds"))
  b2 <- readRDS(file.path(d2, "block_0030.rds"))
  expect_identical(b1$draws, b2$draws)
  expect_true(b1$demo)
  expect_equal(b1$target_draws, 25L)
})

test_that("prev_strata blocks run every draw at the pinned prevalence", {
  m <- grass_calibration_manifest("prev_strata")
  blk <- m$block_id[m$k == 5L & m$N == 30L & m$q == 0.85 & m$prev == 0.20][1L]
  d <- file.path(tempdir(), "contrib_strata_test")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  suppressMessages(grass_contribute(dir = d, blocks = blk, draws = 10L))
  x <- readRDS(file.path(d, sprintf("block_%04d.rds", blk)))
  expect_true(all(x$draws$prev == 0.20))
})

test_that("verification passes a clean bundle and catches tampering", {
  d <- file.path(tempdir(), "contrib_verify_test")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  suppressMessages(grass_contribute(dir = d, blocks = c(30L, 2100L),
                                    draws = 25L))
  v <- suppressMessages(grass_verify_contribution(d, draws = 10L))
  expect_true(all(v$checksum_ok))
  expect_true(all(v$complete))
  expect_true(all(v$replay_ok))
  # tamper with one block file: checksum must fail
  f <- file.path(d, "block_0030.rds")
  x <- readRDS(f); x$draws$delta[1L] <- x$draws$delta[1L] + 1; saveRDS(x, f)
  expect_warning(v2 <- grass_verify_contribution(d, draws = 10L),
                 "did not fully verify")
  expect_false(v2$checksum_ok[v2$block_id == 30L])
})

test_that("replay tolerates last-bit BLAS drift but not fabricated draws", {
  a <- matrix(c(0.0037574410853693507, 0.97486828240252910,
                0.97483070799167537, 0.97486828240252910), 1L, 4L)
  expect_true(replay_matches(a, a, 1e-8))

  # Rachel's block 10, PR #4: one ULP in the AC1 mean, scaled by 100 into
  # `delta` because delta is carried in quality percentage points.
  drift <- a
  drift[1L, 1L] <- 0.0037574410853551399
  drift[1L, 3L] <- 0.97483070799167548
  expect_true(max(abs(drift - a)) < 1e-13)
  expect_false(identical(drift, a))
  expect_true(replay_matches(drift, a, 1e-8))

  # A block that did not come from the pipeline cannot land within tol.
  fake <- a
  fake[1L, 1L] <- a[1L, 1L] + 1e-6
  expect_false(replay_matches(fake, a, 1e-8))

  # Non-finite draws must agree position for position, not just in count.
  na1 <- matrix(c(NA_real_, 1, 2, 3), 1L, 4L)
  na2 <- matrix(c(1, NA_real_, 2, 3), 1L, 4L)
  expect_false(replay_matches(na1, na2, 1e-8))
  expect_true(replay_matches(na1, na1, 1e-8))

  all_na <- matrix(NA_real_, 1L, 4L)
  expect_true(replay_matches(all_na, all_na, 1e-8))
  expect_false(replay_matches(a, a[, 1:3, drop = FALSE], 1e-8))
})

test_that("grass_verify_contribution accepts drift and rejects fabrication", {
  d <- file.path(tempdir(), "contrib_tol_test")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  suppressMessages(grass_contribute(dir = d, blocks = 30L, draws = 25L))
  f <- file.path(d, "block_0030.rds")
  mf <- file.path(d, "bundle_manifest.csv")

  perturb <- function(by) {
    x <- readRDS(f)
    x$draws$delta[1L] <- x$draws$delta[1L] + by
    saveRDS(x, f)
    b <- utils::read.csv(mf, stringsAsFactors = FALSE)
    b$md5 <- unname(tools::md5sum(f))
    utils::write.csv(b, mf, row.names = FALSE)
  }

  perturb(1e-14)
  v <- suppressMessages(grass_verify_contribution(d, draws = 10L))
  expect_true(v$checksum_ok)
  expect_true(v$replay_ok)

  perturb(1e-3)
  expect_warning(v2 <- grass_verify_contribution(d, draws = 10L),
                 "did not fully verify")
  expect_true(v2$checksum_ok)
  expect_false(v2$replay_ok)
})

test_that("skipping the discarded ICC fit leaves the draw stream untouched", {
  skip_if_not_installed("lme4")

  # one_null_draw() reports delta_hat and the three implied qualities, and
  # never ICC. The glmer fit behind ICC therefore only costs time -- provided
  # it draws no random numbers. If it ever did, a contributor without lme4
  # would produce a different stream than a maintainer with it, and every
  # seed-replay across the open calibration program would break. Pin it.
  draw <- function(N, k, prev, q, fit) {
    res <- suppressWarnings(suppressMessages(
      check_asymmetry(gen_null_panel(N, k, prev, q), fit_icc = fit)))
    p <- res$panel
    q_of <- function(nm) {
      i <- which(p$coefficient == nm)
      if (length(i) == 1L) p$implied_q[i] else NA_real_
    }
    c(as.numeric(res$delta_hat), q_of("pabak"), q_of("mean_ac1"),
      q_of("fleiss_kappa"))
  }

  for (cell in list(c(40, 5, 0.20, 0.85), c(30, 3, 0.50, 0.75))) {
    N <- cell[1]; k <- cell[2]; prev <- cell[3]; q <- cell[4]
    set.seed(4242)
    with_icc <- t(vapply(1:5, function(i) draw(N, k, prev, q, TRUE), numeric(4)))
    set.seed(4242)
    no_icc <- t(vapply(1:5, function(i) draw(N, k, prev, q, FALSE), numeric(4)))
    expect_identical(no_icc, with_icc)
  }

  # The panel still carries ICC when asked, and drops it when not.
  set.seed(7)
  Y <- gen_null_panel(60L, 4L, 0.4, 0.9)
  full <- suppressWarnings(suppressMessages(check_asymmetry(Y)))
  lean <- suppressWarnings(suppressMessages(check_asymmetry(Y, fit_icc = FALSE)))
  expect_true("icc" %in% full$panel$coefficient)
  expect_false("icc" %in% lean$panel$coefficient)
  expect_identical(as.numeric(lean$delta_hat), as.numeric(full$delta_hat))
})
