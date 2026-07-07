# Contributor pipeline: run open calibration blocks from the installed
# package and produce a submission bundle. The draw machinery mirrors the
# stage-6B null program (simulation/v070_program/tier1/stage6_null_deltaB/):
# one set.seed() per block before any RNG, prevalence-major sweep, and a
# draw order independent of batching, so a block is bit-reproducible from
# its manifest seed at a pinned package version.

CALIB_PREVS <- c(0.05, 0.20, 0.50, 0.80, 0.95)

#' Open calibration blocks
#'
#' The bundled reference surfaces and null distributions are finite, and the
#' bounds they carry are compute bounds rather than method bounds. This
#' manifest lists the open calibration blocks: seeded, bite-sized simulation
#' cells that any machine running this package can compute and submit. See
#' the CONTRIBUTING file in the source repository
#' (\url{https://github.com/defense031/grassr}) for the submission protocol,
#' and \code{\link{grass_contribute}} to run blocks.
#'
#' A block is one cell of a null-calibration program. Blocks with
#' \code{prev = NA} sweep the five calibration prevalences inside the cell
#' (one fifth of the draws each), exactly like the shipped null. Blocks with
#' a \code{prev} value run every draw at that single prevalence.
#'
#' The bundled manifest is a snapshot from this package version. The live
#' manifest, including which blocks are already claimed, is the copy on the
#' repository's \code{calibration-contrib} branch.
#'
#' @param program Optional program name to filter to. One of
#'   \code{"tail_topup"} (deeper draws for shipped null cells whose extreme
#'   tails carry a transparency flag), \code{"prev_strata"} (a
#'   prevalence-stratified null), or \code{"lattice_k"} (intermediate rater
#'   counts between the shipped grid points).
#' @return A data frame with one row per block: \code{block_id},
#'   \code{program}, the cell's \code{k}, \code{N}, \code{q}, \code{prev},
#'   the required \code{draws}, the block's \code{seed}, and \code{status}.
#' @examples
#' m <- grass_calibration_manifest()
#' table(m$program)
#' @export
grass_calibration_manifest <- function(program = NULL) {
  path <- system.file("extdata", "calibration_manifest.csv",
                      package = "grassr", mustWork = TRUE)
  m <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!is.null(program)) {
    program <- match.arg(program, unique(m$program))
    m <- m[m$program == program, , drop = FALSE]
    rownames(m) <- NULL
  }
  m
}

# Null-program panel generator. Identical draw order to the stage-6B script:
# subject classes first, then one rbinom per rater column.
gen_null_panel <- function(N, k, prev, q) {
  C <- stats::rbinom(N, 1L, prev)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) Y[, j] <- stats::rbinom(N, 1L, ifelse(C == 1L, q, 1 - q))
  Y
}

# One null draw -> c(delta_hat, implied q for PABAK / mean AC1 / Fleiss).
one_null_draw <- function(N, k, prev, q) {
  res <- tryCatch(suppressWarnings(suppressMessages(
    check_asymmetry(gen_null_panel(N, k, prev, q)))),
    error = function(e) NULL)
  if (is.null(res)) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  p <- res$panel
  q_of <- function(nm) {
    i <- which(p$coefficient == nm)
    if (length(i) == 1L) p$implied_q[i] else NA_real_
  }
  c(as.numeric(res$delta_hat), q_of("pabak"), q_of("mean_ac1"),
    q_of("fleiss_kappa"))
}

# Run one manifest block. Seeded once, prevalence-major, sequential; the
# draw order never depends on how the loop is chunked.
run_calibration_block <- function(block, draws = NULL) {
  target <- if (is.null(draws)) block$draws else as.integer(draws)
  prevs <- if (is.na(block$prev)) CALIB_PREVS else block$prev
  per_prev <- target %/% length(prevs)
  target <- per_prev * length(prevs)
  set.seed(block$seed)
  delta <- numeric(target); qp <- numeric(target); qa <- numeric(target)
  qf <- numeric(target); pv <- numeric(target)
  pos <- 0L
  t0 <- proc.time()[["elapsed"]]
  for (prev in prevs) {
    for (r in seq_len(per_prev)) {
      d <- one_null_draw(block$N, block$k, prev, block$q)
      pos <- pos + 1L
      delta[pos] <- d[1L]; qp[pos] <- d[2L]; qa[pos] <- d[3L]; qf[pos] <- d[4L]
      pv[pos] <- prev
    }
  }
  list(block = block[, c("block_id", "program", "k", "N", "q", "prev",
                         "draws", "seed")],
       draws = data.frame(delta = delta, q_pabak = qp, q_mean_ac1 = qa,
                          q_fleiss_kappa = qf, prev = pv),
       target_draws = target,
       n_finite = sum(is.finite(delta)),
       demo = !is.null(draws) && !identical(target, as.integer(block$draws)),
       package_version = as.character(utils::packageVersion("grassr")),
       r_version = R.version.string,
       wall_secs = proc.time()[["elapsed"]] - t0)
}

# Per-draw cost model, fit on two anchor cells: cost = a + b * (N * k).
benchmark_draw_cost <- function(bench_draws = 8L) {
  anchors <- list(c(k = 3L, N = 30L), c(k = 8L, N = 300L))
  secs <- vapply(anchors, function(a) {
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(bench_draws)) one_null_draw(a[["N"]], a[["k"]], 0.5, 0.85)
    (proc.time()[["elapsed"]] - t0) / bench_draws
  }, numeric(1L))
  nk <- vapply(anchors, function(a) a[["N"]] * a[["k"]], numeric(1L))
  b <- (secs[2L] - secs[1L]) / (nk[2L] - nk[1L])
  a <- secs[1L] - b * nk[1L]
  function(k, N) pmax(a + b * (N * k), min(secs) / 4)
}

#' Run open calibration blocks and build a submission bundle
#'
#' Runs one or more open calibration blocks from
#' \code{\link{grass_calibration_manifest}} on this machine and writes a
#' submission bundle to \code{dir}. Each block is seeded from the manifest,
#' so the run is bit-reproducible at this package version and maintainers
#' verify a submission by re-executing its seeds. Everything runs locally;
#' the function makes no network connections.
#'
#' With \code{dry_run = TRUE} the function benchmarks this machine on two
#' small anchor cells, estimates the wall time of every candidate block, and
#' returns the plan that fits the time budget without running it. This is
#' the answer to "I have five hours to give: what would that buy?".
#'
#' Blocks are drawn at random from the open candidates so that uncoordinated
#' contributors rarely collide. Duplicate runs of a block are not wasted;
#' they cross-verify. For a large budget, claim a specific block range
#' through a repository issue and pass it as \code{blocks}.
#'
#' @param dir Directory to write the bundle to. Created if missing. Must be
#'   supplied; nothing is written anywhere else.
#' @param hours Time budget in hours. Blocks are selected so the estimated
#'   total stays inside it. Ignored when \code{blocks} is supplied.
#' @param program Optional program filter, as in
#'   \code{\link{grass_calibration_manifest}}.
#' @param blocks Optional integer vector of \code{block_id}s to run,
#'   typically a range claimed through a repository issue.
#' @param dry_run If \code{TRUE}, benchmark and return the plan only.
#' @param draws Override the per-block draw count, for demonstrations and
#'   tests. A bundle built with overridden draws is marked \code{demo} and
#'   is not submittable.
#' @return Invisibly, a data frame with one row per selected block and its
#'   estimated (and, after a real run, actual) wall time.
#' @examples
#' \donttest{
#' # What would five hours on this machine buy?
#' plan <- grass_contribute(dir = tempdir(), hours = 5, dry_run = TRUE)
#' head(plan)
#' }
#' @export
grass_contribute <- function(dir, hours = 1, program = NULL, blocks = NULL,
                             dry_run = FALSE, draws = NULL) {
  if (missing(dir) || !is.character(dir) || length(dir) != 1L)
    stop("`dir` must be a single path; grass_contribute() writes only there.")
  m <- grass_calibration_manifest(program)
  m <- m[m$status == "open", , drop = FALSE]
  if (!is.null(blocks)) {
    m <- m[m$block_id %in% blocks, , drop = FALSE]
    if (nrow(m) == 0L) stop("none of the requested block_ids are open blocks")
  }
  if (nrow(m) == 0L) stop("no open blocks match the request")

  message("Benchmarking this machine on two small anchor cells...")
  cost <- benchmark_draw_cost()
  m$est_secs <- cost(m$k, m$N) * (if (is.null(draws)) m$draws else draws)

  if (is.null(blocks)) {
    m <- m[sample.int(nrow(m)), , drop = FALSE]     # collision avoidance
    fit <- cumsum(m$est_secs) <= hours * 3600
    if (!any(fit)) fit[1L] <- TRUE                  # always offer one block
    m <- m[fit, , drop = FALSE]
  }
  m <- m[order(m$block_id), , drop = FALSE]
  rownames(m) <- NULL

  message(sprintf(
    "%d block(s), %s draws, estimated %.1f hours on this machine.",
    nrow(m), format(sum(if (is.null(draws)) m$draws else
      rep(draws, nrow(m))), big.mark = ","), sum(m$est_secs) / 3600))
  if (dry_run) return(invisible(m))

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  m$wall_secs <- NA_real_
  for (i in seq_len(nrow(m))) {
    res <- run_calibration_block(m[i, ], draws = draws)
    f <- file.path(dir, sprintf("block_%04d.rds", m$block_id[i]))
    saveRDS(res, f)
    m$wall_secs[i] <- res$wall_secs
    message(sprintf("block %d/%d done (block_id %d, %.1f min): %s",
                    i, nrow(m), m$block_id[i], res$wall_secs / 60, f))
  }

  files <- file.path(dir, sprintf("block_%04d.rds", m$block_id))
  bundle <- data.frame(m[, c("block_id", "program", "k", "N", "q", "prev",
                             "draws", "seed", "wall_secs")],
                       file = basename(files),
                       md5 = unname(tools::md5sum(files)),
                       demo = !is.null(draws),
                       package_version =
                         as.character(utils::packageVersion("grassr")),
                       r_version = R.version.string)
  utils::write.csv(bundle, file.path(dir, "bundle_manifest.csv"),
                   row.names = FALSE)
  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(dir, "SESSIONINFO.txt"))
  message("Bundle written. To submit: fork https://github.com/defense031/grassr, ",
          "add this directory under contrib/contributions/<your-name>-<date>/ ",
          "on the `calibration-contrib` branch, and open a pull request. ",
          "See CONTRIBUTING.md in the repository.")
  if (!is.null(draws))
    message("NOTE: draw count was overridden; this bundle is a demo and ",
            "is not submittable.")
  invisible(m)
}

#' Verify a contribution bundle
#'
#' Checks a bundle written by \code{\link{grass_contribute}}: file checksums
#' against the bundle manifest, package-version match, completeness of each
#' block, and bit-reproducibility of the first draws of every block from its
#' manifest seed. Contributors can run it before submitting; maintainers run
#' it (and a full re-execution of sampled blocks) at intake.
#'
#' @param dir The bundle directory.
#' @param draws Number of leading draws per block to replay from the seed.
#'   The draw order is sequential from one seeding, so a leading replay is
#'   exact regardless of how the original run was chunked.
#' @return Invisibly, a data frame with one row per block and logical
#'   columns \code{checksum_ok}, \code{complete}, \code{replay_ok}.
#' @export
grass_verify_contribution <- function(dir, draws = 200L) {
  mf <- file.path(dir, "bundle_manifest.csv")
  if (!file.exists(mf)) stop("no bundle_manifest.csv in `dir`")
  bundle <- utils::read.csv(mf, stringsAsFactors = FALSE)
  pkg_now <- as.character(utils::packageVersion("grassr"))
  out <- bundle[, c("block_id", "program", "k", "N", "q", "prev", "draws",
                    "seed")]
  out$checksum_ok <- out$complete <- out$replay_ok <- NA
  for (i in seq_len(nrow(bundle))) {
    f <- file.path(dir, bundle$file[i])
    out$checksum_ok[i] <- file.exists(f) &&
      identical(unname(tools::md5sum(f)), bundle$md5[i])
    if (!out$checksum_ok[i]) next
    x <- readRDS(f)
    out$complete[i] <- isTRUE(nrow(x$draws) == x$target_draws) &&
      (isTRUE(bundle$demo[i]) || isTRUE(x$target_draws == bundle$draws[i]))
    if (!identical(x$package_version, pkg_now)) {
      warning(sprintf(
        "block %d was run under grassr %s; this is %s. Replay needs the same version.",
        bundle$block_id[i], x$package_version, pkg_now))
      next
    }
    n_replay <- min(as.integer(draws), x$target_draws)
    prevs <- if (is.na(bundle$prev[i])) CALIB_PREVS else bundle$prev[i]
    per_prev <- x$target_draws %/% length(prevs)
    set.seed(bundle$seed[i])
    replay <- matrix(NA_real_, n_replay, 4L)
    pos <- 0L
    for (prev in prevs) {
      for (r in seq_len(per_prev)) {
        pos <- pos + 1L
        if (pos > n_replay) break
        replay[pos, ] <- one_null_draw(bundle$N[i], bundle$k[i], prev,
                                       bundle$q[i])
      }
      if (pos >= n_replay) break
    }
    got <- as.matrix(x$draws[seq_len(n_replay),
                             c("delta", "q_pabak", "q_mean_ac1",
                               "q_fleiss_kappa")])
    dimnames(got) <- NULL
    out$replay_ok[i] <- identical(got, replay)
  }
  ok <- all(out$checksum_ok %in% TRUE) && all(out$complete %in% TRUE) &&
    all(out$replay_ok[!is.na(out$replay_ok)] %in% TRUE)
  message(sprintf("%d block(s): checksums %s, completeness %s, replay %s.",
                  nrow(out),
                  if (all(out$checksum_ok %in% TRUE)) "ok" else "FAILED",
                  if (all(out$complete %in% TRUE)) "ok" else "FAILED",
                  if (all(out$replay_ok[!is.na(out$replay_ok)] %in% TRUE))
                    "ok" else "FAILED"))
  if (!ok) warning("bundle did not fully verify; see the returned table")
  invisible(out)
}
