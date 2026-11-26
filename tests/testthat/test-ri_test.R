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

test_that("missingness is printed as counts and percentages", {
  # 3 of 14 missing overall, 2 of 8 treated, 1 of 6 control.
  Z <- c(1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
  Y <- c(NA, 1, 2, 3, 4, 5, 6, 7, NA, 8, 9, 10, 11, NA)

  out <- ri_sharp_attrition(
    Z = Z, Y = Y, nperm = 100,
    include_ci = FALSE, include_twostep = FALSE
  )
  printed <- paste(capture.output(print(out)), collapse = "\n")

  expect_true(grepl("Overall +3 / 14 +21\\.4%", printed))
  expect_true(grepl("Treatment +2 / 8 +25\\.0%", printed))
  expect_true(grepl("Control +1 / 6 +16\\.7%", printed))
  expect_true(grepl("Unadjusted: +\\+8\\.3%", printed))
  expect_true(grepl("Design: individual-level complete randomization", printed, fixed = TRUE))

  # Stored values stay proportions; only the display changes.
  expect_equal(out$counts$attrition_rate_treat, 0.25)

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

test_that("assumption hint compares treated and control attrition rates", {
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
  # The figure names the observed statistic and its value in full.
  expect_true(startsWith(viz$subtitle[[1]], paste0("Observed statistic = ", .ri_format_num(viz$stat_obs))))
  expect_true(grepl("p-value", viz$subtitle[[1]], fixed = TRUE))
  expect_true(startsWith(viz$subtitle[[2]], "Complete randomization"))
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

  out_default <- ri_sharp_attrition(Z = Z, Y = Y, missing = "mn", nperm = 100, include_ci = FALSE)
  expect_equal(out_default$results$Method, c("Sharp", "Sharp"))

  out <- ri_sharp_attrition(
    Z = Z, Y = Y, missing = "mn", nperm = 100,
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
      nperm = 100,
      style = "density"
    )
  )
})

test_that("plot_ri_test_stat() supports sharp missingness", {
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
      analysis = "sharp",
      nperm = 100
    )
  )

  expect_equal(viz$analysis, "sharp")
  expect_equal(viz$assumption, "sharp")
  expect_equal(viz$mode_label, "sharp missingness")
  expect_length(viz$stat_null, 100)
  expect_equal(viz$counts$n_obs, 4)   # only the observed outcomes are used
})

test_that("the figure and the table give the same p-value with the same seed", {
  set.seed(3)
  Z <- sample(rep(c(1, 0), 30))
  Y <- rnorm(60) + 0.4 * Z
  Y[sample(60, 12)] <- NA

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  set.seed(10)
  table_rows <- ri_sharp_attrition(Z = Z, Y = Y, missing = "mp", nperm = 500,
                                   include_ci = FALSE)$results
  for (assumption in c("general", "mp")) {
    set.seed(10)
    viz <- plot_ri_test_stat(Z = Z, Y = Y, analysis = assumption,
                             missing = assumption, nperm = 500)
    row <- table_rows$`Missingness assumption` == assumption
    expect_equal(viz$p_value, table_rows$`P-value`[row])
    expect_equal(viz$stat_obs, table_rows$`Test stat. value`[row])
  }
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

test_that("units are put in a random order before ranking", {
  # .ri_random_order() returns a permutation and leaves the random-number state unchanged.
  set.seed(3)
  new_order <- .ri_random_order(10)
  after <- runif(1)
  set.seed(3)
  expect_equal(sort(new_order), 1:10)
  expect_equal(after, runif(1))

  # Without ties between the two arms, the order of the rows does not matter.
  set.seed(4)
  Z <- rep(c(1, 0), each = 30)
  Y <- rnorm(60) + 0.5 * Z
  Y[c(1:5, 31:40)] <- NA
  p_for_order <- function(o) {
    set.seed(9)
    ri_sharp_attrition(Z = Z[o], Y = Y[o], missing = "general", nperm = 500,
                       include_ci = FALSE)$results$`P-value`
  }
  expect_equal(p_for_order(1:60), p_for_order(60:1))
})

