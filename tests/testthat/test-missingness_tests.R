simulate_missing <- function(rate_treat, rate_control, n = 200) {
  set.seed(3)
  Z <- rep(c(1, 0), each = n / 2)
  Y <- rnorm(n)
  observed <- ifelse(Z == 1, runif(n) < rate_treat, runif(n) < rate_control)
  Y[!observed] <- NA
  list(Z = Z, Y = Y)
}

test_that("with all possible assignments, the p-value is exact", {
  # 10 units, 4 treated: all choose(10, 4) = 210 assignments as the columns of
  # Z.perm. The p-value is then exact and equals the one-sided Fisher exact test.
  Z <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  Y <- c(NA, NA, 1, 1, 1, 1, 1, NA, 1, 1)
  M <- as.numeric(!is.na(Y))
  Z.all <- apply(utils::combn(10, 4), 2, function(treated) as.numeric(seq_len(10) %in% treated))

  tab <- table(factor(Z, levels = c(1, 0)), factor(M, levels = c(1, 0)))
  expect_equal(test_monotone_missingness(Z, Y, "mp", Z.perm = Z.all),
               fisher.test(tab, alternative = "less")$p.value)
  expect_equal(test_monotone_missingness(Z, Y, "mn", Z.perm = Z.all),
               fisher.test(tab, alternative = "greater")$p.value)
})

test_that("for sharp, with all possible assignments, the p-value is the exact two-sided one", {
  # Same example. Under sharp every unit keeps its M, so the number of treated
  # units with an observed outcome, k, is hypergeometric. The p-value adds the
  # probabilities of all k whose difference in response rates is at least as
  # far from 0 as the observed one.
  Z <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  Y <- c(NA, NA, 1, 1, 1, 1, 1, NA, 1, 1)
  M <- as.numeric(!is.na(Y))
  Z.all <- apply(utils::combn(10, 4), 2, function(treated) as.numeric(seq_len(10) %in% treated))

  n_observed <- sum(M)
  k <- 0:4
  stat_k <- k / 4 - (n_observed - k) / 6
  stat_obs <- mean(M[Z == 1]) - mean(M[Z == 0])
  exact <- sum(dhyper(k, n_observed, 10 - n_observed, 4)[abs(stat_k) >= abs(stat_obs) - 1e-9])

  expect_equal(test_monotone_missingness(Z, Y, "sharp", Z.perm = Z.all), exact)
})

test_that("for sharp, a gap in response rates in either direction is evidence", {
  # Lower response rate among treated units.
  d <- simulate_missing(0.7, 0.85)
  set.seed(1)
  expect_lt(test_monotone_missingness(d$Z, d$Y, "sharp", nperm = 10^4), 0.05)

  # Higher response rate among treated units.
  d <- simulate_missing(0.9, 0.7)
  set.seed(1)
  expect_lt(test_monotone_missingness(d$Z, d$Y, "sharp", nperm = 10^4), 0.05)

  # Same response rate in both groups.
  d <- simulate_missing(0.8, 0.8)
  set.seed(1)
  expect_gt(test_monotone_missingness(d$Z, d$Y, "sharp", nperm = 10^4), 0.1)
})

test_that("with random assignments, the p-value is close to the exact one", {
  d <- simulate_missing(0.7, 0.85)
  tab <- table(factor(d$Z, levels = c(1, 0)), factor(!is.na(d$Y), levels = c(TRUE, FALSE)))

  set.seed(1)
  p_mp <- test_monotone_missingness(d$Z, d$Y, assumption = "mp", nperm = 10^4)
  set.seed(1)
  p_mn <- test_monotone_missingness(d$Z, d$Y, assumption = "mn", nperm = 10^4)
  # Monte Carlo standard error is at most 0.005 with 10,000 assignments.
  expect_lt(abs(p_mp - fisher.test(tab, alternative = "less")$p.value), 0.02)
  expect_lt(abs(p_mn - fisher.test(tab, alternative = "greater")$p.value), 0.02)

  # Fewer treated units observed: evidence against mp, not against mn.
  expect_lt(p_mp, 0.05)
  expect_gt(p_mn, 0.5)
})

test_that("a supplied Z.perm is used as given and draws no random numbers", {
  d <- simulate_missing(0.7, 0.85)
  set.seed(2)
  Z.perm <- assign_CRE(length(d$Z), sum(d$Z), 500)

  before <- .Random.seed
  p1 <- test_monotone_missingness(d$Z, d$Y, "mp", Z.perm = Z.perm)
  expect_identical(.Random.seed, before)
  expect_identical(test_monotone_missingness(d$Z, d$Y, "mp", Z.perm = Z.perm), p1)
})

test_that("ri_sharp_attrition() reports the test in the notes for mp, mn and sharp only", {
  d <- simulate_missing(0.9, 0.7)

  set.seed(1)
  out_mp <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "mp", nperm = 2000, include_ci = FALSE)
  expect_true(out_mp$monotone_pval > 0.5)
  printed <- paste(capture.output(print(out_mp)), collapse = " ")
  expect_true(grepl("Testing this assumption gives p", printed, fixed = TRUE))

  set.seed(1)
  out_mn <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "mn", nperm = 2000, include_ci = FALSE)
  expect_true(out_mn$monotone_pval < 0.05)

  set.seed(1)
  out_sharp <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "sharp", nperm = 2000, include_ci = FALSE)
  expect_true(out_sharp$monotone_pval < 0.05)
  printed <- paste(capture.output(print(out_sharp)), collapse = " ")
  expect_true(grepl("Testing this assumption gives p", printed, fixed = TRUE))

  for (missing in c("general", "random")) {
    set.seed(1)
    out_other <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = missing, nperm = 2000, include_ci = FALSE)
    expect_null(out_other$monotone_pval)
    expect_false(grepl("Testing this assumption",
                       paste(capture.output(print(out_other)), collapse = " "), fixed = TRUE))
  }
})

test_that("the test does not change the other results of ri_sharp_attrition()", {
  # The test reuses the re-assignments of the general and mp rows, so the
  # general row is the same with missing = "mp" as with missing = "general",
  # and so is the random-number state after the call.
  d <- simulate_missing(0.9, 0.7)

  set.seed(4)
  out_mp <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "mp", nperm = 500, include_ci = FALSE)
  seed_after_mp <- .Random.seed
  set.seed(4)
  out_general <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "general", nperm = 500, include_ci = FALSE)
  seed_after_general <- .Random.seed

  expect_equal(out_mp$results[1, "P-value"], out_general$results[1, "P-value"])
  expect_identical(seed_after_mp, seed_after_general)
})

test_that("an unknown assumption is an error", {
  expect_error(test_monotone_missingness(c(1, 0), c(1, NA), "general", nperm = 10),
               "Unrecognized assumption")
  # random is tested by test_missing_at_random(), which needs covariates.
  expect_error(test_monotone_missingness(c(1, 0), c(1, NA), "random", nperm = 10),
               "Unrecognized assumption")
})


#----------- test_missing_at_random() -----------------------------------------#

# 1,000 units, half treated, two covariates. `age_effect` says how missingness
# depends on age: "none", response rates of 60% among treated and 85% among
# control units, unrelated to age; "same", older units less often observed,
# with the same probability at each age in both groups; or "opposite", rates
# of 60% and 85% at average age, with older units less often observed among
# treated units and more often among control units.
simulate_covariates <- function(age_effect, n = 1000) {
  set.seed(5)
  Z <- rep(c(1, 0), each = n / 2)
  X <- data.frame(age = rnorm(n), female = rbinom(n, 1, 0.5))
  Y <- rnorm(n)
  logit <- qlogis(ifelse(Z == 1, 0.6, 0.85))
  if (age_effect == "same") {
    logit <- qlogis(0.75) - X$age
  }
  if (age_effect == "opposite") {
    logit <- logit + ifelse(Z == 1, -1, 1) * X$age
  }
  Y[runif(n) > plogis(logit)] <- NA
  list(Z = Z, Y = Y, X = X)
}

test_that("the statistic is the F statistic of the regression of Z on X among observed units", {
  d <- simulate_covariates("opposite", n = 300)
  observed <- !is.na(d$Y)
  fit <- lm(Z ~ age + female, data = data.frame(Z = d$Z, d$X)[observed, ])
  expect_equal(.balance_F_stat(model.matrix(fit), d$Z[observed]),
               unname(summary(fit)$fstatistic["value"]))
})

test_that("the p-value is the share of re-randomizations with at least the observed statistic", {
  d <- simulate_covariates("none", n = 200)
  observed <- !is.na(d$Y)
  dat <- data.frame(Z = d$Z, d$X)[observed, ]

  # The F statistic computed with lm().
  F_by_lm <- function(z) {
    dat_z <- dat
    dat_z$Z <- z
    unname(summary(lm(Z ~ age + female, data = dat_z))$fstatistic["value"])
  }

  # 200 re-randomizations of Z among the observed units.
  set.seed(9)
  F_drawn <- rep(NA, 200)
  for (iter in 1:200) {
    F_drawn[iter] <- F_by_lm(sample(dat$Z))
  }
  p_by_hand <- mean(F_drawn >= F_by_lm(dat$Z) - 1e-9)

  set.seed(9)
  expect_equal(test_missing_at_random(d$Z, d$Y, d$X, nperm = 200), p_by_hand)
})

test_that("only imbalance among the observed units is evidence", {
  # Different response rates (60% vs 85%) unrelated to the covariates.
  d <- simulate_covariates("none")
  set.seed(1)
  expect_gt(test_missing_at_random(d$Z, d$Y, d$X, nperm = 2000), 0.1)

  # Missingness related to age but not to the group: the observed units are
  # still balanced, so the test does not detect it.
  d <- simulate_covariates("same")
  set.seed(1)
  expect_gt(test_missing_at_random(d$Z, d$Y, d$X, nperm = 2000), 0.05)

  # Missingness related to age in opposite directions: the observed treated
  # units are younger than the observed control units.
  d <- simulate_covariates("opposite")
  set.seed(1)
  expect_lt(test_missing_at_random(d$Z, d$Y, d$X, nperm = 2000), 0.001)
})

test_that("with a rare 0/1 covariate in a small sample, the test keeps its level", {
  # 500 experiments with 100 units in which missing at random holds (80%
  # observed in both groups); one covariate equals 1 for about 10% of units.
  # With 200 re-randomizations the test rejects at level 0.05 with probability
  # 11/201 = 0.055; 0.09 leaves room for simulation error (standard error 0.01).
  set.seed(6)
  p <- rep(NA, 500)
  for (r in 1:500) {
    Z <- sample(rep(c(1, 0), each = 50))
    X <- data.frame(age = rnorm(100), rare = rbinom(100, 1, 0.1), female = rbinom(100, 1, 0.5))
    Y <- ifelse(runif(100) < 0.8, rnorm(100), NA)
    p[r] <- test_missing_at_random(Z, Y, X, nperm = 200)
  }
  expect_lt(mean(p <= 0.05), 0.09)
})

test_that("units with a missing covariate are left out", {
  d <- simulate_covariates("opposite", n = 300)
  d$X$age[1:10] <- NA
  set.seed(2)
  p_with_na <- test_missing_at_random(d$Z, d$Y, d$X, nperm = 500)
  set.seed(2)
  p_removed <- test_missing_at_random(d$Z[-(1:10)], d$Y[-(1:10)], d$X[-(1:10), ], nperm = 500)
  expect_equal(p_with_na, p_removed)
})

test_that("unusable inputs are errors", {
  d <- simulate_covariates("none", n = 100)
  expect_error(test_missing_at_random(d$Z, d$Y, d$X[1:50, ]), "one row per unit")

  X_text <- d$X
  X_text$region <- rep(c("north", "south"), 50)
  expect_error(test_missing_at_random(d$Z, d$Y, X_text), "numeric covariates only")

  Y_no_treated <- d$Y
  Y_no_treated[d$Z == 1] <- NA
  expect_error(test_missing_at_random(d$Z, Y_no_treated, d$X), "in both groups")

  expect_error(test_missing_at_random(d$Z, d$Y, data.frame(constant = rep(1, 100))),
               "No covariate varies")
})

test_that("ri_sharp_attrition(missing = \"random\", X = X) reports the test in the notes", {
  d <- simulate_covariates("opposite", n = 300)

  set.seed(1)
  out <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "random", X = d$X,
                            nperm = 500, include_ci = FALSE)
  # The observed units are unbalanced on age here, so the p-value is small.
  expect_lt(out$random_pval, 0.01)
  notes <- .ri_notes_text(out)
  expect_true(any(grepl("A check of this assumption", notes, fixed = TRUE)))
})

test_that("X does not change the table, and is ignored for the other assumptions and with clusters", {
  d <- simulate_covariates("opposite", n = 300)

  # The table is the same with or without X: the re-randomizations of the test
  # are drawn after the table is computed.
  set.seed(1)
  out_x <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "random", X = d$X,
                              nperm = 500, include_ci = FALSE)
  set.seed(1)
  out_no_x <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "random",
                                 nperm = 500, include_ci = FALSE)
  expect_equal(out_x$results, out_no_x$results)
  expect_null(out_no_x$random_pval)
  expect_false(any(grepl("A check of this assumption", .ri_notes_text(out_no_x), fixed = TRUE)))

  # Other assumptions: no test of missing at random.
  for (missing in c("general", "mp", "mn", "sharp")) {
    set.seed(1)
    out_other <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = missing, X = d$X,
                                    nperm = 200, include_ci = FALSE)
    expect_null(out_other$random_pval)
  }

  # Clusters of 10 units (each cluster is in one group): no test either.
  set.seed(1)
  out_cluster <- ri_sharp_attrition(Z = d$Z, Y = d$Y, missing = "random", X = d$X,
                                    cluster = rep(1:30, each = 10),
                                    nperm = 200, include_ci = FALSE)
  expect_null(out_cluster$random_pval)
})
