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
