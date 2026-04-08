test_that("ri_test() returns a printable applied-user result", {
  Z <- c(1, 1, 1, 0, 0, 0)
  Y <- c(NA, 8, 6, 3, 4, NA)

  out <- ri_test(
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
  expect_equal(colnames(out$results), c("Method", "Test stat.", "P-value", "CI lower", "CI upper"))
  expect_equal(out$results$Method[[1]], "Sharp-null")
  expect_true(is.numeric(out$results$`Test stat.`[[1]]))
  expect_true(is.numeric(out$results$`P-value`[[1]]))
  expect_invisible(print(out))
})

test_that("ri_test() supports blocked randomization summaries", {
  Z <- c(1, 0, 1, 0, 1, 0)
  Y <- c(2, NA, 3, 1, NA, 0)
  block <- c("A", "A", "B", "B", "C", "C")

  out <- ri_test(
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
  expect_equal(nrow(out$block_diagnostics$block_table), 3)
})
