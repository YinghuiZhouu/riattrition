test_that("pval_sharp() works", {
  set.seed(1)
  N <- 500
  Y0 <- rnorm(N, 0, 1)
  Y1 <- Y0

  # (1) Threshold Missingness
  set.seed(1)
  M0 <- as.numeric(Y0 <= qnorm(0.95))
  M1 <- as.numeric(Y1 >= qnorm(0.05))
  Z <- sample(c(rep(1, N/2), rep(0, N/2)))
  M <- Z * M1 + (1 - Z) * M0
  Y <- Z * Y1 + (1 - Z) * Y0 # No missing value
  Y_observed <- ifelse(M == 0, NA, Y)
  p_value_g <- pval_sharp(Z, Y_observed, c = 0, missing = "general", class = "RS",
                          method.list = list(name = "Wilcoxon"),
                          stat.null = NULL, Z.perm = NULL, nperm = 10^4)
  p_value_g_twostep <- pval_sharp_twostep(Z, Y_observed, c = 0, missing = "general",
                                          method.list = list(name = "Wilcoxon"),
                                          stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta = 0.1 * 0.1)
  expect_equal(p_value_g, 0.7241)
  expect_equal(p_value_g_twostep, 0.7383)

  # (2) Monotone Positive Missingness
  set.seed(1)
  M0 <- as.numeric(Y0 <= qnorm(0.92))
  M1 <- as.numeric(Y0 <= qnorm(0.98))
  Z <- sample(c(rep(1, N/2), rep(0, N/2)))
  M <- Z * M1 + (1 - Z) * M0
  Y <- Z * Y1 + (1 - Z) * Y0 # No missing value
  Y_observed <- ifelse(M == 0, NA, Y)
  p_value_mp <- pval_sharp(Z, Y_observed, c = 0, missing = "mp", class = "RS",
                           method.list = list(name = "Wilcoxon"),
                           stat.null = NULL, Z.perm = NULL, nperm = 10^4)
  p_value_mp_twostep <- pval_sharp_twostep(Z, Y_observed, c = 0, missing = "mp",
                                          method.list = list(name = "Wilcoxon"),
                                          stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta = 0.1 * 0.1)
  expect_equal(p_value_mp, 0.7294)
  expect_equal(p_value_mp_twostep, 0.7436)

  # (3) Monotone Negative Missingness
  set.seed(1)
  M1 <- as.numeric(Y0 >= qnorm(0.08))
  M0 <- as.numeric(Y0 >= qnorm(0.02))
  Z <- sample(c(rep(1, N/2), rep(0, N/2)))
  M <- Z * M1 + (1 - Z) * M0
  Y <- Z * Y1 + (1 - Z) * Y0 # No missing value
  Y_observed <- ifelse(M == 0, NA, Y)
  p_value_mn <- pval_sharp(Z, Y_observed, c = 0, missing = "mn", class = "RS",
                           method.list = list(name = "Wilcoxon"),
                           stat.null = NULL, Z.perm = NULL, nperm = 10^4)
  p_value_mn_twostep <- pval_sharp_twostep(Z, Y_observed, c = 0, missing = "mn",
                                           method.list = list(name = "Wilcoxon"),
                                           stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta = 0.1 * 0.1)
  expect_equal(p_value_mn, 0.7216)
  expect_equal(p_value_mn_twostep, 0.736)

  # (4) Sharp Missingness
  set.seed(1)
  M0 <- as.numeric(Y0 <= qnorm(0.95))
  M1 <- M0
  Z <- sample(c(rep(1, N/2), rep(0, N/2)))
  M <- Z * M1 + (1 - Z) * M0
  Y <- Z * Y1 + (1 - Z) * Y0 # No missing value
  Y_observed <- ifelse(M == 0, NA, Y)
  p_value_s <- pval_sharp(Z, Y_observed, c = 0, missing = "sharp", class = "RS",
                          method.list = list(name = "Wilcoxon"),
                          stat.null = NULL, Z.perm = NULL, nperm = 10^4)
  expect_equal(p_value_s, 0.755)

  # (5) Missing at random
  set.seed(1)
  p <- 0.95
  M0 <- rbinom(N, 1, p)
  M1 <- rbinom(N, 1, p)
  Z <- sample(c(rep(1, N/2), rep(0, N/2)))
  M <- Z * M1 + (1 - Z) * M0
  Y <- Z * Y1 + (1 - Z) * Y0 # No missing value
  Y_observed <- ifelse(M == 0, NA, Y)
  p_value_r <- pval_sharp(Z, Y_observed, c = 0, missing = "random", class = "RS",
                          method.list = list(name = "Wilcoxon"),
                          stat.null = NULL, Z.perm = NULL, nperm = 10^4)
  expect_equal(p_value_r, 0.1044)

})
