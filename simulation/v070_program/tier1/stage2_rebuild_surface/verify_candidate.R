#!/usr/bin/env Rscript
# Stage 2 verification — run the package test suite against the candidate
# sysdata in a scratch copy. Zero failures -> DONE. Any failures ->
# NEEDS_REVIEW with the list (expected class: tests pinning old-grid snap
# behavior, since denser N/q grids change which cell a lookup resolves to).

STAGE2 <- "grassr/simulation/v070_program/tier1/stage2_rebuild_surface"
CANDIDATE <- file.path(STAGE2, "candidate_sysdata.rda")
stopifnot(file.exists(CANDIDATE))

scratch <- file.path(tempdir(), "grassr_stage2_verify")
unlink(scratch, recursive = TRUE)
dir.create(scratch, recursive = TRUE)
file.copy("grassr", scratch, recursive = TRUE)
pkg <- file.path(scratch, "grassr")
unlink(file.path(pkg, "simulation"), recursive = TRUE)  # keep the copy light
file.copy(CANDIDATE, file.path(pkg, "R", "sysdata.rda"), overwrite = TRUE)

res <- as.data.frame(testthat::test_local(pkg, reporter = "silent"))
fails <- res[res$failed > 0 | res$error, c("file", "test")]
cat(sprintf("suite vs candidate: %d passed, %d failed/errored\n",
            sum(res$passed), nrow(fails)))

if (nrow(fails) == 0L) {
  writeLines(sprintf("stage2 verified %s: 0 failures", Sys.time()),
             file.path(STAGE2, "DONE"))
  cat("DONE written.\n")
} else {
  writeLines(c(sprintf("stage2 %s: %d failures need adjudication", Sys.time(), nrow(fails)),
               apply(fails, 1, paste, collapse = " :: ")),
             file.path(STAGE2, "NEEDS_REVIEW"))
  cat("NEEDS_REVIEW written:\n")
  print(fails, row.names = FALSE)
}
