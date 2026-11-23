simulate_missing <- function(rate_treat, rate_control, n = 200) {
  set.seed(3)
  Z <- rep(c(1, 0), each = n / 2)
  Y <- rnorm(n)
  observed <- ifelse(Z == 1, runif(n) < rate_treat, runif(n) < rate_control)
  Y[!observed] <- NA
  list(Z = Z, Y = Y)
}

test_that("under complete randomization the p-value equals the one-sided Fisher exact test", {
  d <- simulate_missing(0.7, 0.85)
  tab <- table(factor(d$Z, levels = c(1, 0)), factor(!is.na(d$Y), levels = c(TRUE, FALSE)))

  mp <- test_monotone_missingness(d$Z, d$Y, assumption = "mp")
  mn <- test_monotone_missingness(d$Z, d$Y, assumption = "mn")
  expect_equal(mp$p_value, fisher.test(tab, alternative = "less")$p.value)
  expect_equal(mn$p_value, fisher.test(tab, alternative = "greater")$p.value)
  expect_equal(mp$method, "exact")

  # Fewer treated units observed: evidence against mp, not against mn.
  expect_lt(mp$p_value, 0.05)
  expect_gt(mn$p_value, 0.5)
  expect_equal(mp$statistic, mean(!is.na(d$Y[d$Z == 1])) - mean(!is.na(d$Y[d$Z == 0])))
})

test_that("the blocked Monte Carlo p-value matches full enumeration in a small example", {
  # Two blocks (4 and 6 units), 2 and 3 treated; enumerate all assignments.
  block <- c(rep("a", 4), rep("b", 6))
  Z <- c(1, 1, 0, 0, 1, 1, 1, 0, 0, 0)
  Y <- c(NA, 1, 1, 1, NA, NA, 1, 1, 1, NA)
  M <- as.numeric(!is.na(Y))

  weighted_diff <- function(z) {
    diff_a <- mean(M[1:4][z[1:4] == 1]) - mean(M[1:4][z[1:4] == 0])
    diff_b <- mean(M[5:10][z[5:10] == 1]) - mean(M[5:10][z[5:10] == 0])
    (4 * diff_a + 6 * diff_b) / 10
  }
  assign_a <- utils::combn(4, 2)
  assign_b <- utils::combn(6, 3)
  stats_all <- c()
  for (i in seq_len(ncol(assign_a))) {
    for (j in seq_len(ncol(assign_b))) {
      z <- numeric(10)
      z[assign_a[, i]] <- 1
      z[4 + assign_b[, j]] <- 1
      stats_all <- c(stats_all, weighted_diff(z))
    }
  }
  exact_mp <- mean(stats_all <= weighted_diff(Z) + 1e-12)

  set.seed(10)
  mc <- test_monotone_missingness(Z, Y, assumption = "mp", block = block, nperm = 20000)
  expect_equal(mc$statistic, weighted_diff(Z))
  expect_equal(mc$p_value, exact_mp, tolerance = 0.02)
  expect_equal(mc$method, "Monte Carlo")
})

test_that("ri_sharp_attrition() reports the test for mp and mn only", {
  d <- simulate_missing(0.9, 0.7)

  out_mp <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "mp", nperm = 200, include_ci = FALSE)
  expect_equal(out_mp$monotone_test$p_value,
               test_monotone_missingness(d$Z, d$Y, assumption = "mp")$p_value)
  printed <- paste(capture.output(print(out_mp)), collapse = "\n")
  expect_true(grepl("Test of M(1) \u2265 M(0) for every unit: p = ", printed, fixed = TRUE))

  out_mn <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "mn", nperm = 200, include_ci = FALSE)
  expect_equal(out_mn$monotone_test$assumption, "mn")
  expect_true(out_mn$monotone_test$p_value < 0.05)

  out_general <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "general", nperm = 200, include_ci = FALSE)
  expect_null(out_general$monotone_test)
  expect_false(grepl("Test of M(1)", paste(capture.output(print(out_general)), collapse = "\n"), fixed = TRUE))
})

test_that("the test does not change the random-number state", {
  set.seed(1)
  invisible(.ri_keep_seed(runif(5)))
  after <- runif(1)
  set.seed(1)
  expect_equal(after, runif(1))
})

test_that("inputs are checked", {
  expect_error(test_monotone_missingness(c(1, 2, 0), c(1, NA, 2)), "0/1")
  expect_error(test_monotone_missingness(c(1, 0), c(1, NA, 2)), "same length")
})
