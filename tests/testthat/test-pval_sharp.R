test_that("pval_sharp() and pval_sharp_twostep() works", {
  Z <- chen_and_yang_2019$treatment_vpnnl
  Y <- chen_and_yang_2019$az_overall_a_w3
  result <- application_wrapper(Z = Z, Y = Y, seed = 1)
  expect_equal(result$pval.g, 0.47185)
  expect_equal(result$pval.mx, 0)
  expect_equal(result$pval.s, 0)
  expect_equal(result$pval.r, 0)

  expect_equal(result$pval.g.twostep, 0.31598)
  expect_equal(result$pval.mx.twostep, 0.01)
})
