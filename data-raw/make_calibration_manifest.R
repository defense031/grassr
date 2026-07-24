# Build inst/extdata/calibration_manifest.csv -- the open calibration blocks
# that grass_contribute() can run.
#
# A block is one null-program cell: a (k, N, q) ridge, optionally pinned to a
# single calibration prevalence. Blocks with prev = NA sweep all five
# calibration prevalences inside the cell (10,000 draws each), exactly like
# the shipped stage-6B null. Blocks with a prev value run all draws at that
# one prevalence (the prevalence-stratified program).
#
# SEEDS. Every block carries an explicit seed, 20300000 + block_id. The
# shipped programs used master seeds 20260706 (null) and 20260707 (power)
# offset by cell_id, so contributed draws can never replay shipped draws.
# Block ids are assigned in the fixed program order below. DO NOT reorder
# programs or grids: reordering reassigns seeds and orphans claimed blocks.
#
# Run from the package root: Rscript data-raw/make_calibration_manifest.R

K_GRID <- c(3L, 5L, 6L, 8L, 10L, 15L, 25L)             # shipped null, k >= 3
N_GRID <- c(15L, 20L, 30L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L)
Q_GRID <- c(0.65, 0.75, 0.85, 0.92, 0.97)
PREVS  <- c(0.05, 0.20, 0.50, 0.80, 0.95)
DRAWS  <- 50000L
SEED_BASE <- 20300000L

# Program 1: tail_topup -- deeper draws for the 21 shipped cells whose
# extreme tails carry the unstable-tail transparency flag. Read the cell
# list from the shipped null itself so the manifest cannot drift from it.
e <- new.env(); load("R/sysdata.rda", envir = e)
idx <- e$delta_null_ecdf$index
tail_cells <- idx[idx$unstable_tail, c("k", "N", "q")]
tail_cells <- tail_cells[order(tail_cells$k, tail_cells$N, tail_cells$q), ]
p1 <- data.frame(program = "tail_topup", k = tail_cells$k, N = tail_cells$N,
                 q = tail_cells$q, prev = NA_real_)

# Program 2: prev_strata -- the shipped null pools five prevalences per
# ridge; the realized flag size runs to roughly twice nominal at extreme
# prevalence. One block per (ridge, prevalence) rebuilds the null with a
# full 50,000 draws in every stratum.
g2 <- expand.grid(q = Q_GRID, N = N_GRID, k = K_GRID, prev = PREVS,
                  KEEP.OUT.ATTRS = FALSE)
g2 <- g2[order(g2$k, g2$N, g2$q, g2$prev), ]
p2 <- data.frame(program = "prev_strata", k = g2$k, N = g2$N, q = g2$q,
                 prev = g2$prev)

# Program 3: lattice_k -- intermediate rater counts. Off-grid designs snap
# to the nearest calibrated cell; these ridges cut the snap distance.
g3 <- expand.grid(q = Q_GRID, N = N_GRID, k = c(4L, 7L, 12L, 20L),
                  KEEP.OUT.ATTRS = FALSE)
g3 <- g3[order(g3$k, g3$N, g3$q), ]
p3 <- data.frame(program = "lattice_k", k = g3$k, N = g3$N, q = g3$q,
                 prev = NA_real_)

manifest <- rbind(p1, p2, p3)
manifest$draws <- DRAWS
manifest$block_id <- seq_len(nrow(manifest))
manifest$seed <- SEED_BASE + manifest$block_id
manifest$status <- "open"
manifest <- manifest[, c("block_id", "program", "k", "N", "q", "prev",
                         "draws", "seed", "status")]

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
write.csv(manifest, "inst/extdata/calibration_manifest.csv", row.names = FALSE)
cat(sprintf("wrote %d blocks (%s)\n", nrow(manifest),
            paste(sprintf("%s: %d", names(table(manifest$program)),
                          table(manifest$program)), collapse = ", ")))
