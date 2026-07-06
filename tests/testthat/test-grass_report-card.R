# Tests for v0.2.0 Target-2 grass_report() / grass_card class.
#
# Phase 4A scope: the new headline `grass_report(ratings = Y)` flow plus the
# non-plot S3 methods (print, summary, as.data.frame, format). Plotting is
# tested separately by Phase 4B.

# ---- Shared fixtures ------------------------------------------------------
# `.simulate_op_strong_panel()` builds the §4 divergent worked-example case
# from the merged paper: k=5 raters with alternating bias direction
# (R1, R3, R5 favor sensitivity at Se=0.95/Sp=0.75; R2, R4 favor specificity
# at Se=0.75/Sp=0.95), F = logit-normal(pi_target=0.50, tau2=0.25), N=1000.
# This DGP is the `op_strong` profile from paper2/code/21_asym_grid.R at
# q = 0.85. seed = 6 produces the divergent panel cited in §4
# (delta_hat ~ 18 pp; AC1 separates from kappa-family on the surface).
.simulate_op_strong_panel <- function(seed = 6L) {
  set.seed(seed)
  N <- 1000L; k <- 5L
  Se <- c(0.95, 0.75, 0.95, 0.75, 0.95)
  Sp <- c(0.75, 0.95, 0.75, 0.95, 0.75)
  mu_logit <- qlogis(0.50)
  p_i <- plogis(rnorm(N, mu_logit, sqrt(0.25)))
  C <- rbinom(N, 1, p_i)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) {
    Y[, j] <- rbinom(N, 1, ifelse(C == 1L, Se[j], 1 - Sp[j]))
  }
  Y
}

# ---- Test 1 ---------------------------------------------------------------
# Basic call returns a grass_card with the documented top-level structure.

test_that("grass_report returns a grass_card with all top-level fields", {
  set.seed(1)
  Y <- matrix(rbinom(500, 1, 0.3), nrow = 100, ncol = 5)
  card <- grass_report(Y, bootstrap_B = 50)

  expect_s3_class(card, "grass_card")
  expect_type(card, "list")
  expected_fields <- c("sample", "coefficient", "delta", "panel", "per_rater",
                       "surface", "call", "grass_version", "timestamp",
                       "inputs", "notes")
  for (f in expected_fields) {
    expect_true(f %in% names(card),
                info = paste("missing field:", f))
  }

  # sample fields
  expect_true(all(c("k", "N", "pi_hat", "tau2_hat", "axis") %in%
                  names(card$sample)))
  expect_equal(card$sample$k, 5L)
  expect_equal(card$sample$N, 100L)
  expect_true(card$sample$pi_hat > 0 && card$sample$pi_hat < 1)
  expect_equal(card$sample$axis, "inter")

  # coefficient fields (v0.7.1: `qualifier` retired; `band` is a rendered
  # consistency-band string, plus percentile_basis / consistency_band / sweep).
  expect_true(all(c("primary", "observed_value", "surface_percentile",
                    "band", "q_hat", "consistency_band") %in%
                  names(card$coefficient)))
  expect_false("qualifier" %in% names(card$coefficient))

  # delta fields
  expect_true(all(c("delta_hat", "flag", "thresholds") %in%
                  names(card$delta)))
  expect_true(card$delta$flag %in% c("aligned", "caution", "divergent"))
  expect_named(card$delta$thresholds, c("caution", "divergent"))

  # panel data.frame structure (v0.7.1 columns: consistency-band endpoints
  # instead of a modal band string / qualifier).
  expect_s3_class(card$panel, "data.frame")
  expect_true(nrow(card$panel) >= 2)
  expect_true(all(c("coefficient", "observed_value", "surface_percentile",
                    "band_lo", "band_hi", "band_open_low", "band_open_high",
                    "q_hat", "se_q_hat", "clamped", "reference_used",
                    "in_delta_hat") %in% names(card$panel)))
  expect_false("qualifier" %in% names(card$panel))
})

# ---- Test 2 ---------------------------------------------------------------
# Paper §4 divergent worked-example reproduction (load-bearing integration
# test). Uses the op_strong asym-grid profile at q=0.85, pi_target=0.50,
# N=1000, seed=6 — see .simulate_op_strong_panel() at top of file.
# Heterogeneous per-rater bias produces genuine cross-coefficient
# discordance; the divergent flag is real (not an ICC-clamp artifact).

test_that("grass_report reproduces the paper §4 divergent worked example", {
  skip_if_not_installed("lme4")  # the fixture pins the 4-row panel (ICC row)
  Y <- .simulate_op_strong_panel(seed = 6L)
  card <- grass_report(Y, bootstrap_B = 200)

  # Re-pinned 2026-07-06 against the stage-6B delta-B null: divergent at
  # the >= 99th percentile of the matched null (delta ~0.25 quality pp).
  expect_equal(card$delta$flag, "divergent")
  expect_gte(card$delta$delta_percentile, 99)

  # The divergent flag must come from real cross-coefficient disagreement,
  # not an ICC-clamp artifact: the non-ICC pooled percentiles also spread
  # (observed ~2 pp at re-pin; aligned panels sit well under 0.5 pp).
  nonicc <- card$panel$surface_percentile[card$panel$coefficient != "icc"]
  expect_true(diff(range(nonicc, na.rm = TRUE)) >= 1,
              info = sprintf("non-ICC pooled-pct range = %.2f pp; expected >= 1",
                             diff(range(nonicc, na.rm = TRUE))))

  # Four-coefficient panel (PABAK, AC1, Fleiss kappa, ICC; alpha removed
  # from the panel at v0.6.0).
  expect_equal(nrow(card$panel), 4L)
  expect_true("clamped" %in% names(card$panel))

  # Per-rater data.frame with five rows of finite Se/Sp + CIs.
  expect_s3_class(card$per_rater, "data.frame")
  expect_equal(nrow(card$per_rater), 5L)
  expect_true(all(is.finite(card$per_rater$se_hat)))
  expect_true(all(is.finite(card$per_rater$sp_hat)))
  expect_true(all(is.finite(card$per_rater$se_lower)))
  expect_true(all(is.finite(card$per_rater$se_upper)))
  expect_true(all(is.finite(card$per_rater$sp_lower)))
  expect_true(all(is.finite(card$per_rater$sp_upper)))

  # Per-rater point estimates separate Se-favoring vs Sp-favoring raters
  # (R1, R3, R5 should have higher Se than Sp; R2, R4 the reverse).
  pr <- card$per_rater
  for (j in c(1L, 3L, 5L)) {
    expect_gt(pr$se_hat[j], pr$sp_hat[j],
              label = sprintf("R%d should be Se-favoring", j))
  }
  for (j in c(2L, 4L)) {
    expect_gt(pr$sp_hat[j], pr$se_hat[j],
              label = sprintf("R%d should be Sp-favoring", j))
  }

  # Headline coefficient band suppressed under the divergent flag.
  expect_equal(card$coefficient$band, "suppressed")
})

# ---- Test 3 ---------------------------------------------------------------
# Paper §4 k=2 symmetric counter-example: 200 subjects, prevalence 0.05,
# Se = Sp = 0.94. Should resolve as `aligned`, primary band Strong /
# Excellent, no per-rater table.

test_that("grass_report handles the §4 k=2 symmetric counter-example", {
  set.seed(1)
  N <- 200; k <- 2
  pi <- 0.05
  Se <- 0.94; Sp <- 0.94
  C <- rbinom(N, 1, pi)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) {
    Y[, j] <- rbinom(N, 1, ifelse(C == 1, Se, 1 - Sp))
  }

  card <- grass_report(Y, bootstrap_B = 50)

  # v0.7.1 delta-B: at k = 2 the two-coefficient family implies identical
  # quality by construction, so delta_hat is structurally uninformative
  # and the flag reports not_applicable (never aligned/caution/divergent).
  expect_equal(card$delta$flag, "not_applicable")
  expect_equal(card$delta$thresholds_source, "not_applicable_k2")
  # v0.7.1: the primary `band` is a rendered consistency-band string on
  # quality (not the retired Strong/Excellent adjective), and is only
  # "suppressed" under the divergent flag.
  expect_type(card$coefficient$band, "character")
  expect_false(identical(card$coefficient$band, "suppressed"))
  expect_match(card$coefficient$band, "quality")
  expect_true(is.list(card$coefficient$consistency_band))
  expect_null(card$per_rater)
})

# ---- Test 4 ---------------------------------------------------------------
# print() emits "GRASS Report Card" header and shape varies with flag.

.make_align_card_k2 <- function() {
  set.seed(1)
  N <- 200; k <- 2
  pi <- 0.05
  Se <- 0.94; Sp <- 0.94
  C <- rbinom(N, 1, pi)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) {
    Y[, j] <- rbinom(N, 1, ifelse(C == 1, Se, 1 - Sp))
  }
  grass_report(Y, bootstrap_B = 0)
}

test_that("print.grass_card renders the aligned card header, gloss and delta", {
  card_align <- .make_align_card_k2()

  expect_output(print(card_align), "GRASS Report Card")

  out_align <- capture.output(print(card_align))
  txt_align <- paste(out_align, collapse = "\n")
  # k = 2 card: delta is not_applicable (v0.7.1) with the n/a matched-null
  # line routing to the pairwise/bounds path.
  expect_match(txt_align, "not_applicable", fixed = TRUE)
  expect_match(txt_align, "n/a at k = 2", fixed = TRUE)
  # v0.7.1: the aligned card carries a plain-language "read:" gloss line.
  expect_match(txt_align, "read:", fixed = TRUE)
  # The primary coefficient (PABAK) and its consistency band on quality show.
  expect_match(txt_align, "PABAK", fixed = TRUE)
  expect_match(txt_align, "percentile", fixed = TRUE)
  # Debug-grade notes (glmer / F_key provenance) stay off the headline card.
  expect_false(grepl("F_key picked via glmer", txt_align, fixed = TRUE))
})

test_that("print.grass_card renders the divergent card (all coefficients + suppressed)", {
  # Re-pinned 2026-07-06: op_strong fires divergent on the delta-B null.
  Y_div <- .simulate_op_strong_panel(seed = 6L)
  card_div <- grass_report(Y_div, bootstrap_B = 50)

  expect_output(print(card_div), "GRASS Report Card")
  out_div <- capture.output(print(card_div))
  txt_div <- paste(out_div, collapse = "\n")
  expect_match(txt_div, "PABAK", fixed = TRUE)
  expect_match(txt_div, "AC1", fixed = TRUE)
  expect_match(txt_div, "Fleiss kappa", fixed = TRUE)
  expect_match(txt_div, "suppressed", fixed = TRUE)
  expect_match(txt_div, "divergent", fixed = TRUE)
})

# ---- Test 5 ---------------------------------------------------------------
# summary() returns a summary.grass_card object, prints multi-section block.

test_that("summary.grass_card returns the documented structure", {
  Y <- .simulate_op_strong_panel(seed = 6L)
  card <- grass_report(Y, bootstrap_B = 50)

  s <- summary(card)
  expect_s3_class(s, "summary.grass_card")
  expect_true("panel" %in% names(s))
  expect_s3_class(s$panel, "data.frame")
  expect_output(print(s), "panel")
})

# ---- Test 6 ---------------------------------------------------------------
# as.data.frame() returns a data.frame for non-divergent cards, list for
# divergent. Either way, the panel rows have an `is_primary` logical column.

test_that("as.data.frame.grass_card has an is_primary column", {
  set.seed(1)
  N <- 200; k <- 2
  pi <- 0.05
  Se <- 0.94; Sp <- 0.94
  C <- rbinom(N, 1, pi)
  Y <- matrix(0L, N, k)
  for (j in seq_len(k)) {
    Y[, j] <- rbinom(N, 1, ifelse(C == 1, Se, 1 - Sp))
  }
  card_align <- grass_report(Y, bootstrap_B = 50)

  df_align <- as.data.frame(card_align)
  expect_s3_class(df_align, "data.frame")
  expect_true("is_primary" %in% names(df_align))
  expect_true(any(df_align$is_primary))
})

test_that("as.data.frame.grass_card returns a list(panel, per_rater) on divergent cards", {
  # Re-pinned 2026-07-06: op_strong fires divergent on the delta-B null.
  Y <- .simulate_op_strong_panel(seed = 6L)
  card_div <- grass_report(Y, bootstrap_B = 50)
  df_div <- as.data.frame(card_div)
  expect_type(df_div, "list")
  expect_true("panel" %in% names(df_div))
  expect_s3_class(df_div$panel, "data.frame")
  expect_true("is_primary" %in% names(df_div$panel))
})

# ---- Test 7 ---------------------------------------------------------------
# verbose = TRUE emits at least one progress message.

test_that("grass_report with verbose = TRUE emits a progress message", {
  set.seed(1)
  Y <- matrix(rbinom(500, 1, 0.3), nrow = 100, ncol = 5)
  expect_message(grass_report(Y, bootstrap_B = 25, verbose = TRUE),
                 regexp = NULL)
})

# ---- v0.2.3 auto-trigger of pairwise_agreement under divergent flag -------
# Per paper §3.3, divergent panels populate `card$pairwise` with the
# pairwise PABAK matrix + per-rater pooled-reference table; aligned panels
# do not. This is the load-bearing wiring for the divergent recovery path.

test_that("grass_report auto-populates card$pairwise on divergent panels", {
  # Re-pinned 2026-07-06: op_strong fires divergent on the delta-B null.
  Y <- .simulate_op_strong_panel(seed = 6L)
  card <- grass_report(Y, bootstrap_B = 50)

  expect_equal(card$delta$flag, "divergent")
  expect_s3_class(card$pairwise, "grass_pairwise")
  expect_equal(dim(card$pairwise$pabak_matrix), c(5L, 5L))
  expect_equal(nrow(card$pairwise$pooled_per_rater), 5L)
})

test_that("grass_report leaves card$pairwise NULL when flag is not divergent", {
  set.seed(2)
  Y <- matrix(rbinom(500, 1, 0.3), nrow = 100, ncol = 5)
  card <- grass_report(Y, bootstrap_B = 25)

  expect_false(card$delta$flag == "divergent")
  expect_null(card$pairwise)
})
