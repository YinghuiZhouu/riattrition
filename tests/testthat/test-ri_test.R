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
