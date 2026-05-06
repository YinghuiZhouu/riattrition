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
  expect_false(grepl("\\[95% C\\.I\\.\\]", printed))
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
  expect_true(is.numeric(out$block_diagnostics$mean_within_block_attrition_diff))
  expect_true(is.numeric(out$block_diagnostics$median_within_block_attrition_diff))
  expect_true(is.numeric(out$block_diagnostics$min_within_block_attrition_diff))
  expect_true(is.numeric(out$block_diagnostics$max_within_block_attrition_diff))
  expect_equal(nrow(out$block_diagnostics$block_table), 3)
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

  expect_false(grepl("Block missingness", printed_hidden, fixed = TRUE))
  expect_true(grepl("Within-block attrition diff (T - C)", printed_expanded, fixed = TRUE))
  expect_true(grepl("Mean       Min        Q1         Median     Q3         Max", printed_expanded, fixed = TRUE))
  expect_false(grepl("Median within-block attrition diff", printed_expanded, fixed = TRUE))
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
  expect_true(grepl("\\[95% C\\.I\\.\\]", printed))
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
  expect_true(grepl("Inference unit", printed, fixed = TRUE))
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
    "When analysis = \"general\""
  )
})
