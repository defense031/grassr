library(testthat)
library(grassr)

# On CRAN, run a fast deterministic smoke subset: core metric arithmetic,
# input normalization, pairwise agreement, print methods, and the sysdata
# regression anchors that pin the bundled reference surfaces at known
# design points. The full suite (~650 tests, incl. bootstrap, plotting,
# and lme4-dependent paths) runs on every push via the package's
# GitHub Actions matrix (windows/macos/ubuntu x devel/release):
# https://github.com/defense031/grassr/actions
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  test_check("grassr")
} else {
  message("grassr: running CRAN smoke subset (full suite runs on CI)")
  test_check(
    "grassr",
    filter = "^(metrics-core|metrics-multi-rater|normalize|pairwise_agreement|sysdata-regression-anchors|print)$"
  )
  message("grassr: CRAN smoke subset complete")
}
