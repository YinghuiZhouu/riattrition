source(testthat::test_path("../../R/0_function_sharp_null.R"))
source(testthat::test_path("../../R/0_function_sharp_null_twostep.R"))
source(testthat::test_path("../../R/1_result_output.R"))
source(testthat::test_path("../../R/2_visualization.R"))

test_that("ri_sharp_attrition() returns a printable applied-user result", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  set.seed(123)
  out <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    c = 0,
    missing = "general",
    class = "RS",
    method.list = list(name = "Wilcoxon"),
    nperm = 200,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  expect_s3_class(out, "riattrition_result")
  expect_equal(out$assumption, "general")
  expect_true("results" %in% names(out))
  expect_equal(
    colnames(out$results),
    c("Missingness assumption", "Method", "Test stat.", "Test stat. value", "P-value", "CI lower", "CI upper")
  )
  expect_equal(out$results$`Missingness assumption`[[1]], "general")
  expect_equal(out$results$Method[[1]], "Sharp")
  expect_equal(out$results$`Test stat.`[[1]], "Wilcoxon rank-sum")
  expect_true(is.numeric(out$results$`Test stat. value`[[1]]))
  expect_true(is.numeric(out$results$`P-value`[[1]]))
  expect_invisible(print(out))

  printed <- paste(capture.output(print(out)), collapse = "\n")
  expect_false(grepl("CI for constant effect", printed, fixed = TRUE))
})

test_that("ri_sharp_attrition() supports blocked randomization summaries", {
  Z <- c(1, 0, 1, 0, 1, 0)
  Y <- c(2, NA, 3, 1, NA, 0)
  block <- c("A", "A", "B", "B", "C", "C")

  out <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    c = 0,
    missing = "general",
    class = "RS",
    method.list = list(name = "Wilcoxon"),
    block = block,
    nperm = 100,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  expect_s3_class(out, "riattrition_result")
  expect_equal(out$design, "blocked randomization")
  expect_true(is.numeric(out$counts$attrition_rate_total))
  expect_true(is.numeric(out$counts$attrition_rate_treat))
  expect_true(is.numeric(out$counts$attrition_rate_control))
  expect_true(is.character(out$assumption_hint))
  expect_equal(out$block_diagnostics$n_blocks, 3)
  expect_true(is.numeric(out$block_diagnostics$weighted_mean_within_block_attrition_diff))
  expect_true(is.numeric(out$block_diagnostics$median_within_block_attrition_diff))
  expect_true(is.numeric(out$block_diagnostics$min_within_block_attrition_diff))
  expect_true(is.numeric(out$block_diagnostics$max_within_block_attrition_diff))
  expect_equal(nrow(out$block_diagnostics$block_table), 3)
})

test_that("within-block attrition diff is averaged with block-size weights", {
  # Block A (4 units): treated attrition 1/2, control 0/2   -> diff  0.50
  # Block B (8 units): treated attrition 0/4, control 1/4   -> diff -0.25
  # Block C (2 units): treated only, so no within-block diff -> dropped
  # Weighted: (4 * 0.50 + 8 * -0.25) / 12 = 0; unweighted would be 0.125.
  Z <- c(1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
  Y <- c(NA, 1, 2, 3, 4, 5, 6, 7, NA, 8, 9, 10, 11, NA)
  block <- c(rep("A", 4), rep("B", 8), rep("C", 2))

  diag <- .ri_block_diagnostics(Z = Z, Y = Y, block = block)

  expect_equal(diag$n_blocks, 3)
  expect_equal(diag$n_blocks_used, 2)
  expect_equal(diag$block_table$attrition_diff, c(0.5, -0.25, NA))
  expect_equal(diag$weighted_mean_within_block_attrition_diff, 0)
  expect_equal(diag$median_within_block_attrition_diff, 0.125)

  out <- ri_sharp_attrition(
    Z = Z, Y = Y, block = block, nperm = 100,
    include_ci = FALSE, include_twostep = FALSE
  )
  printed <- paste(capture.output(print(out)), collapse = "\n")
  expect_true(grepl("Blocks: 3 total; 2 contain both treatment arms.", printed, fixed = TRUE))
})

test_that("missingness is printed as counts and percentages", {
  # Same blocks as above: 3 of 14 missing overall, 2 of 8 treated, 1 of 6 control;
  # within-block diffs are +50% (A) and -25% (B), weighted mean 0%.
  Z <- c(1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
  Y <- c(NA, 1, 2, 3, 4, 5, 6, 7, NA, 8, 9, 10, 11, NA)
  block <- c(rep("A", 4), rep("B", 8), rep("C", 2))

  out <- ri_sharp_attrition(
    Z = Z, Y = Y, block = block, nperm = 100,
    include_ci = FALSE, include_twostep = FALSE
  )
  printed <- paste(capture.output(print(out)), collapse = "\n")

  expect_true(grepl("Overall +3 / 14 +21\\.4%", printed))
  expect_true(grepl("Treatment +2 / 8 +25\\.0%", printed))
  expect_true(grepl("Control +1 / 6 +16\\.7%", printed))
  expect_true(grepl("Unadjusted: +\\+8\\.3%", printed))
  expect_true(grepl("Block-weighted: +0\\.0%", printed))
  expect_true(grepl("Median: +12\\.5%", printed))
  expect_true(grepl("Range: +[[]\u221225\\.0%, 50\\.0%[]]", printed))

  # Stored values stay proportions; only the display changes.
  expect_equal(out$counts$attrition_rate_treat, 0.25)
  expect_equal(out$block_diagnostics$max_within_block_attrition_diff, 0.5)

  # Cluster designs: clusters observe 1/2, 2/2, 2/2, 0/2 -> response rate 62.5%.
  Zc <- c(1, 1, 1, 1, 0, 0, 0, 0)
  Yc <- c(10, NA, 8, 6, 2, 4, NA, NA)
  cluster <- c("A", "A", "B", "B", "C", "C", "D", "D")
  out_cluster <- ri_sharp_attrition(
    Z = Zc, Y = Yc, cluster = cluster, nperm = 100,
    include_ci = FALSE, include_twostep = FALSE
  )
  printed_cluster <- paste(capture.output(print(out_cluster)), collapse = "\n")
  expect_true(grepl("Mean cluster response rate +62\\.5%", printed_cluster))
  expect_true(grepl("Overall +1 / 4 +25\\.0%", printed_cluster))
  expect_true(grepl("Control +1 / 2 +50\\.0%", printed_cluster))
})

test_that("assumption hint compares arms within blocks in blocked designs", {
  # Block A (high attrition): 8 treated (3 missing), 2 control (1 missing)
  # Block B (low attrition):  2 treated (0 missing), 8 control (1 missing)
  # Within each block treated attrition is lower (diff -0.125 in both), but
  # pooled treated attrition is higher (0.30 vs 0.20): Simpson's paradox.
  Z <- c(rep(1, 8), rep(0, 2), rep(1, 2), rep(0, 8))
  M <- c(0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1)
  Y <- ifelse(M == 1, seq_along(M), NA)
  block <- rep(c("A", "B"), each = 10)

  out <- ri_sharp_attrition(
    Z = Z, Y = Y, block = block, nperm = 100,
    include_ci = FALSE, include_twostep = FALSE
  )

  expect_equal(out$counts$attrition_rate_treat - out$counts$attrition_rate_control, 0.1)
  expect_equal(out$block_diagnostics$weighted_mean_within_block_attrition_diff, -0.125)
  expect_true(startsWith(out$assumption_hint, "Within blocks (block-size weighted), treated attrition is lower"))
  expect_true(grepl("more consistent with mp than mn", out$assumption_hint, fixed = TRUE))
  expect_true(grepl("pooled attrition rates point the other way", out$assumption_hint, fixed = TRUE))

  printed <- paste(capture.output(print(out)), collapse = "\n")
  expect_true(grepl("Unadjusted: +\\+10\\.0%", printed))
  expect_true(grepl("Block-weighted: +\u221212\\.5%", printed))

  # The mn note uses the block-weighted difference: treated units drop out less
  # within blocks, so the rates point against mn and agree with mp.
  out_mn <- ri_sharp_attrition(Z = Z, Y = Y, block = block, missing = "mn", nperm = 100, include_ci = FALSE)
  out_mp <- ri_sharp_attrition(Z = Z, Y = Y, block = block, missing = "mp", nperm = 100, include_ci = FALSE)
  expect_true(grepl("point the other way", .ri_notes_text(out_mn)[2], fixed = TRUE))
  expect_true(grepl("are compatible with this assumption", .ri_notes_text(out_mp)[2], fixed = TRUE))
})

test_that("assumption hint adds no reversal note when pooled and within-block agree", {
  # Same treated share in both blocks, so pooled and within-block diffs match:
  # block A diff = 0 - 1/2 = -0.5, block B diff = 0, weighted = pooled = -0.25.
  Z <- rep(c(1, 1, 0, 0), 2)
  Y <- c(1, 2, NA, 3, 4, 5, 6, 7)
  block <- rep(c("A", "B"), each = 4)

  counts <- .ri_missing_counts(Z, Y)
  diag <- .ri_block_diagnostics(Z = Z, Y = Y, block = block)
  hint <- .ri_assumption_hint(counts, diag)

  expect_equal(counts$attrition_rate_treat - counts$attrition_rate_control, -0.25)
  expect_equal(diag$weighted_mean_within_block_attrition_diff, -0.25)
  expect_equal(
    hint,
    "Within blocks (block-size weighted), treated attrition is lower than control attrition; this is more consistent with mp than mn, although rates alone do not identify the missingness mechanism."
  )
  expect_false(grepl("point the other way", hint, fixed = TRUE))
})

test_that("assumption hint without blocks still uses pooled attrition rates", {
  Z <- c(1, 1, 0, 0)
  Y <- c(NA, 1, 2, 3)

  expect_equal(
    .ri_assumption_hint(.ri_missing_counts(Z, Y)),
    "Treated attrition is higher than control attrition; this is more consistent with mn than mp, although rates alone do not identify the missingness mechanism."
  )
})

test_that("general is always reported first when user requests another assumption", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  out <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    c = 0,
    missing = "mp",
    class = "RS",
    method.list = list(name = "Wilcoxon"),
    nperm = 200,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  expect_equal(out$results$`Missingness assumption`, c("general", "mp"))
})

test_that("block_missing_sum controls expanded block printing", {
  Z <- c(1, 0, 1, 0, 1, 0)
  Y <- c(2, NA, 3, 1, NA, 0)
  block <- c("A", "A", "B", "B", "C", "C")

  out_hidden <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    block = block,
    block_missing_sum = FALSE,
    nperm = 100,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  out_expanded <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    block = block,
    block_missing_sum = TRUE,
    nperm = 100,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  printed_hidden <- paste(capture.output(print(out_hidden)), collapse = "\n")
  printed_expanded <- paste(capture.output(print(out_expanded)), collapse = "\n")

  expect_true(grepl("Blocks: 3 total", printed_hidden, fixed = TRUE))
  expect_false(grepl("Within-block treatment", printed_hidden, fixed = TRUE))
  expect_true(grepl("Within-block treatment \u2212 control differences:", printed_expanded, fixed = TRUE))
  expect_true(grepl("IQR:", printed_expanded, fixed = TRUE))
})

test_that("CI column is shown when include_ci = TRUE", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  out <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    nperm = 200,
    include_ci = TRUE,
    include_twostep = FALSE
  )

  printed <- paste(capture.output(print(out)), collapse = "\n")
  expect_true(grepl("95% CI for constant effect", printed, fixed = TRUE))
})

test_that("ri_sharp_attrition() supports cluster-level analysis", {
  Z <- c(1, 1, 1, 1, 0, 0, 0, 0)
  Y <- c(10, NA, 8, 6, 2, 4, NA, NA)
  cluster <- c("A", "A", "B", "B", "C", "C", "D", "D")

  out <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    cluster = cluster,
    nperm = 100,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  expect_equal(out$inference_unit, "cluster")
  expect_equal(out$design, "cluster randomization")
  expect_equal(out$counts$n_total, 4)
  expect_equal(out$counts$n_miss, 1)
  expect_equal(out$cluster_diagnostics$n_individuals, 8)
  expect_equal(out$cluster_diagnostics$n_clusters, 4)
  expect_equal(out$cluster_diagnostics$clusters_with_no_observed_outcome, 1)
  expect_equal(out$cluster_diagnostics$clusters_with_any_missing, 2)

  printed <- paste(capture.output(print(out)), collapse = "\n")
  expect_true(grepl("Number of clusters", printed, fixed = TRUE))
  expect_true(grepl("Cluster summary", printed, fixed = TRUE))
  expect_true(grepl("Design: cluster randomization", printed, fixed = TRUE))
})

test_that("ri_sharp_attrition() supports blocked cluster randomization", {
  Z <- c(1, 1, 0, 0, 1, 1, 0, 0)
  Y <- c(10, 9, 3, 4, 8, NA, 2, 1)
  cluster <- c("A", "A", "B", "B", "C", "C", "D", "D")
  block <- c("X", "X", "X", "X", "Y", "Y", "Y", "Y")

  out <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    cluster = cluster,
    block = block,
    nperm = 100,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  expect_equal(out$design, "blocked cluster randomization")
  expect_equal(out$block_diagnostics$n_blocks, 2)
  expect_equal(out$counts$n_treat, 2)
  expect_equal(out$counts$n_control, 2)
})

test_that("cluster_missing = any_missing is available as a sensitivity option", {
  Z <- c(1, 1, 1, 1, 0, 0, 0, 0)
  Y <- c(10, NA, 8, 6, 2, 4, NA, NA)
  cluster <- c("A", "A", "B", "B", "C", "C", "D", "D")

  out_all <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    cluster = cluster,
    cluster_missing = "all_missing",
    nperm = 100,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  out_any <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    cluster = cluster,
    cluster_missing = "any_missing",
    nperm = 100,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  expect_equal(out_all$counts$n_miss, 1)
  expect_equal(out_any$counts$n_miss, 2)
})

test_that("cluster-level analysis requires constant treatment within cluster", {
  Z <- c(1, 0, 1, 1)
  Y <- c(10, 9, 8, 6)
  cluster <- c("A", "A", "B", "B")

  expect_error(
    ri_sharp_attrition(
      Z = Z,
      Y = Y,
      cluster = cluster,
      nperm = 100,
      include_ci = FALSE,
      include_twostep = FALSE
    ),
    "Treatment assignment must be constant within each cluster"
  )
})

test_that("plot_ri_test_stat() returns visualization data aligned with ri_sharp_attrition()", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  set.seed(123)
  out <- ri_sharp_attrition(
    Z = Z,
    Y = Y,
    c = 0,
    missing = "general",
    class = "RS",
    method.list = list(name = "Wilcoxon"),
    nperm = 200,
    include_ci = FALSE,
    include_twostep = FALSE
  )

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  set.seed(123)
  viz <- plot_ri_test_stat(
    Z = Z,
    Y = Y,
    c = 0,
    missing = "general",
    class = "RS",
    method.list = list(name = "Wilcoxon"),
    nperm = 200
  )

  expect_s3_class(viz, "riattrition_plot")
  expect_equal(viz$assumption, "general")
  expect_equal(viz$test_name, "Wilcoxon rank-sum")
  expect_equal(viz$stat_obs, out$results$`Test stat. value`[[1]])
  expect_equal(viz$p_value, out$results$`P-value`[[1]])
  expect_length(viz$stat_null, 200)
  expect_equal(viz$off_support_mode, "inside")
  expect_equal(viz$obs_position, "inside")
  expect_equal(viz$mc_replications, 200)
  expect_true(grepl("MC replications = 200", viz$subtitle[[2]], fixed = TRUE))
})

test_that("the header states the outcome, design, test and sharp null", {
  my_data <- data.frame(Z = c(1, 1, 1, 0, 0, 0), y = c(NA, 8, 6, 3, 4, NA))

  out <- ri_sharp_attrition(Z = my_data$Z, Y = my_data$y, nperm = 100, include_ci = FALSE)
  lines <- capture.output(print(out))
  expect_equal(lines[1], "Randomization inference with attrition")
  expect_equal(lines[2], strrep("=", nchar(lines[1])))
  for (title in c("Outcome missingness", "Inference", "Notes")) {
    at <- which(lines == title)
    expect_length(at, 1)
    expect_equal(lines[at + 1], strrep("-", nchar(title)))
  }
  expect_true("Outcome: y" %in% lines)
  expect_true("Design: individual-level complete randomization" %in% lines)
  expect_true("Test: Wilcoxon rank-sum with 100 Monte Carlo permutations" %in% lines)
  expect_true("Sharp null: Y_i(1) \u2212 Y_i(0) = 0 for every unit" %in% lines)

  out_named <- ri_sharp_attrition(Z = my_data$Z, Y = my_data$y, nperm = 100, include_ci = FALSE, outcome_name = "Media")
  expect_true("Outcome: Media" %in% capture.output(print(out_named)))

  out_shift <- ri_sharp_attrition(Z = my_data$Z, Y = my_data$y, c = -2, nperm = 100, include_ci = FALSE)
  lines_shift <- capture.output(print(out_shift))
  expect_true("Sharp null: Y_i(1) \u2212 Y_i(0) = \u22122 for every unit" %in% lines_shift)
  expect_true(grepl("treatment effect of exactly \u22122.", paste(lines_shift, collapse = " "), fixed = TRUE))
})

test_that("p-values and confidence intervals are formatted as in the meeting notes", {
  expect_equal(.ri_pvalue(0), "<0.001")
  expect_equal(.ri_pvalue(0.0004), "<0.001")
  expect_equal(.ri_pvalue(1), "1.000")
  expect_equal(.ri_ci_text(-Inf, Inf), "(\u2212\u221e, \u221e)")
  expect_equal(.ri_ci_text(0.011, 0.066), "[0.011, 0.066]")
  expect_equal(.ri_ci_text(-0.5, Inf), "[\u22120.500, \u221e)")
  expect_equal(.ri_pct(-0.127), "\u221212.7%")
  expect_equal(.ri_pct(0.011, plus = TRUE), "+1.1%")
  expect_equal(.ri_pct(-0.0004, plus = TRUE), "0.0%")
})

test_that("by default the table has only the unrestricted row and its notes", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  out <- ri_sharp_attrition(Z = Z, Y = Y, nperm = 100, include_ci = FALSE)
  notes <- .ri_notes_text(out)

  expect_equal(out$results$Method, "Sharp")
  expect_length(notes, 2)
  expect_true(startsWith(notes[1], "Unrestricted missingness uses least-favorable imputation"))
  expect_equal(notes[2], "P-values test the sharp null that every unit has zero treatment effect.")

  lines <- capture.output(print(out))
  expect_true(any(grepl("^Missingness assumption +MC p-value$", lines)))
  expect_true(any(grepl("^Unrestricted \\(worst case\\) +[0-9.<]+$", lines)))
})

test_that("two-step rows are added only when asked, with their own note", {
  Z <- c(rep(1, 8), rep(0, 2), rep(1, 2), rep(0, 8))
  M <- c(0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1)
  Y <- ifelse(M == 1, seq_along(M), NA)
  block <- rep(c("A", "B"), each = 10)

  out_default <- ri_sharp_attrition(Z = Z, Y = Y, block = block, missing = "mn", nperm = 100, include_ci = FALSE)
  expect_equal(out_default$results$Method, c("Sharp", "Sharp"))

  out <- ri_sharp_attrition(
    Z = Z, Y = Y, block = block, missing = "mn", nperm = 100,
    include_ci = TRUE, include_twostep = TRUE
  )
  lines <- capture.output(print(out))
  notes <- .ri_notes_text(out)

  expect_true(any(startsWith(lines, "Unrestricted (worst case), two-step")))
  expect_true(any(startsWith(lines, "Monotone negative, M(1) \u2264 M(0), two-step")))
  expect_true(any(grepl("\u03b2 = 0.005", notes, fixed = TRUE)))
  expect_true(grepl("The 95% confidence intervals are obtained by inverting tests of constant additive effects, Y_i(1) \u2212 Y_i(0) = \u03c4.", notes[length(notes)], fixed = TRUE))
  expect_true(all(nchar(lines) <= 90))
})

test_that("notes explain sharp and random only when they are in the table", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  out <- ri_sharp_attrition(Z = Z, Y = Y, missing = "random", nperm = 100, include_ci = FALSE)
  notes <- .ri_notes_text(out)

  expect_length(notes, 3)
  expect_true(startsWith(notes[2], "Missing at random means"))
  expect_false(any(grepl("M(1) \u2264 M(0)", notes, fixed = TRUE)))
})

test_that("MC replications are stored and printed with a uniform label", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  out <- ri_sharp_attrition(
    Z = Z, Y = Y, missing = "mp", nperm = 1500,
    include_ci = FALSE, include_twostep = FALSE
  )
  printed_lines <- capture.output(print(out))
  printed <- paste(printed_lines, collapse = "\n")
  expect_equal(out$mc_replications, 1500)
  expect_true(grepl("Test: Wilcoxon rank-sum with 1,500 Monte Carlo permutations", printed, fixed = TRUE))
  expect_false(grepl("nperm", printed, fixed = TRUE))

  # A user-supplied assignment matrix determines the count actually used.
  set.seed(1)
  Z.perm <- assign_CRE(n = length(Z), m = sum(Z), nperm = 50)
  out_custom <- ri_sharp_attrition(
    Z = Z, Y = Y, Z.perm = Z.perm, nperm = 10^4,
    include_ci = FALSE, include_twostep = FALSE
  )
  expect_equal(out_custom$mc_replications, 50)

  # Large counts are written out in full rather than as 1e+05.
  expect_equal(.ri_format_count(10^5), "100,000")
})

test_that("plot_ri_test_stat() supports density style", {
  Z <- c(1, 0, 1, 0, 1, 0)
  Y <- c(2, NA, 3, 1, NA, 0)
  block <- c("A", "A", "B", "B", "C", "C")

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  expect_no_error(
    plot_ri_test_stat(
      Z = Z,
      Y = Y,
      block = block,
      nperm = 100,
      style = "density"
    )
  )
})

test_that("plot_ri_test_stat() supports complete-case visualization", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  expect_no_error(
    viz <- plot_ri_test_stat(
      Z = Z,
      Y = Y,
      analysis = "complete-case",
      nperm = 100
    )
  )

  expect_equal(viz$analysis, "complete-case")
  expect_equal(viz$assumption, "complete-case")
  expect_equal(viz$mode_label, "complete-case analysis")
  expect_length(viz$stat_null, 100)
})

find_off_support_example <- function(target_mode = c("extend", "split"),
                                     nperm = 50,
                                     seed = 1,
                                     max_iter = 1000) {
  target_mode <- match.arg(target_mode)
  set.seed(seed)

  for (iter in seq_len(max_iter)) {
    Z <- c(rep(1, 6), rep(0, 6))
    Y <- sample(c(NA, -2:8), 12, replace = TRUE)

    tf <- tempfile(fileext = ".pdf")
    grDevices::pdf(tf)
    viz <- try(
      plot_ri_test_stat(
        Z = Z,
        Y = Y,
        nperm = nperm,
        off_support = "auto"
      ),
      silent = TRUE
    )
    grDevices::dev.off()
    unlink(tf)

    if (!inherits(viz, "try-error") && identical(viz$off_support_mode, target_mode)) {
      return(list(Z = Z, Y = Y, viz = viz))
    }
  }

  stop("Could not find a deterministic off-support example for testing.")
}

test_that("plot_ri_test_stat() can extend the x-axis for modest off-support gaps", {
  example <- find_off_support_example(target_mode = "extend")
  Z <- example$Z
  Y <- example$Y

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  expect_no_error(
    viz <- plot_ri_test_stat(
      Z = Z,
      Y = Y,
      nperm = 100,
      off_support = "auto"
    )
  )

  expect_equal(viz$obs_position, "left")
  expect_equal(viz$off_support_mode, "extend")
})

test_that("plot_ri_test_stat() supports split-axis displays for off-support statistics", {
  example <- find_off_support_example(target_mode = "extend")
  Z <- example$Z
  Y <- example$Y

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  expect_no_error(
    viz <- plot_ri_test_stat(
      Z = Z,
      Y = Y,
      nperm = 100,
      off_support = "split"
    )
  )

  expect_equal(viz$obs_position, "left")
  expect_equal(viz$off_support_mode, "split")
})

test_that("plot_ri_test_stat() supports blocked complete-case visualization", {
  Z <- c(1, 0, 1, 0, 1, 0)
  Y <- c(2, NA, 3, 1, NA, 0)
  block <- c("A", "A", "B", "B", "C", "C")

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  expect_no_error(
    plot_ri_test_stat(
      Z = Z,
      Y = Y,
      block = block,
      analysis = "complete-case",
      nperm = 100
    )
  )
})

test_that("plot_ri_test_stat() restricts non-general assumptions in general mode", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  expect_error(
    plot_ri_test_stat(
      Z = Z,
      Y = Y,
      missing = "mn",
      nperm = 100
    ),
    "`missing` should match the requested analysis"
  )
})
