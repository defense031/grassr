#!/usr/bin/env Rscript
# Build the grassr candidate (0.6.2 code + stage-2 candidate sysdata) into
# stage3_threshold_grid/lib/ as an isolated library, so stage-3 percentile
# inversions run against the densified surface without touching the
# machine's installed grassr. Run from repo/bundle root on each machine.

ROOT <- Sys.getenv("GRASS_SIM_ROOT",
                   "/Users/austinsemmel/Desktop/PABAK_Investigation")
setwd(ROOT)
STAGE2 <- "grassr/simulation/v070_program/tier1/stage2_rebuild_surface"
STAGE3 <- "grassr/simulation/v070_program/tier1/stage3_threshold_grid"
CANDIDATE <- file.path(STAGE2, "candidate_sysdata.rda")
stopifnot(file.exists(CANDIDATE))

src <- file.path(tempdir(), "grassr")
unlink(src, recursive = TRUE)
dir.create(src)
# copy package sources only (skip simulation/, docs/, .git)
for (d in c("R", "man", "inst", "tests", "vignettes", "DESCRIPTION",
            "NAMESPACE", "LICENSE", "NEWS.md", ".Rbuildignore")) {
  p <- file.path("grassr", d)
  if (file.exists(p)) file.copy(p, src, recursive = TRUE)
}
file.copy(CANDIDATE, file.path(src, "R", "sysdata.rda"), overwrite = TRUE)

lib <- file.path(ROOT, STAGE3, "lib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
r <- system2(file.path(R.home("bin"), "R"),
             c("CMD", "INSTALL", "--no-docs", "--no-multiarch",
               paste0("--library=", lib), src))
stopifnot(r == 0L)
.libPaths(c(lib, .libPaths()))
library(grassr, lib.loc = lib)
n <- nrow(get("empirical_q_hat_surface", envir = asNamespace("grassr"))$index)
cat(sprintf("candidate lib ready: grassr %s, surface cells %d\n",
            as.character(packageVersion("grassr")), n))
stopifnot(n > 40000L)   # densified surface really is in
