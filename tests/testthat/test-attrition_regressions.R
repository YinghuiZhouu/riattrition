helper_path <- testthat::test_path("../../analysis/attrition_regressions.R")
skip_if_not(file.exists(helper_path), "analysis/attrition_regressions.R is not available")
source(helper_path)

simulate_attrition <- function(n = 400) {
  set.seed(42)
  X <- data.frame(female = rbinom(n, 1, 0.5), age = rnorm(n))
  Z <- rbinom(n, 1, 0.6)
  M <- rbinom(n, 1, plogis(0.5 + 0.4 * X$female - 0.3 * Z * X$age))
  list(M = M, Z = Z, X = X)
}

test_that("group differences equal treatment minus control coefficients", {
  d <- simulate_attrition()
  ar <- attrition_regressions(d$M, d$Z, d$X, labels = c(female = "Female", age = "Age"))
  co <- attrition_group_coefs(ar)

  expect_equal(co$label, c("Female", "Age", "Constant"))
  expect_equal(co$diff_est, co$treatment_est - co$control_est)
})

test_that("HC1 standard errors and joint tests match fixest", {
  skip_if_not_installed("fixest")
  d <- simulate_attrition()
  ar <- attrition_regressions(d$M, d$Z, d$X)
  dat <- data.frame(M = d$M, Z = d$Z, d$X)
  ref <- fixest::feols(M ~ Z * (female + age), data = dat, vcov = "hetero")

  terms <- names(ar$fits$reg_b$coef)
  expect_equal(unname(ar$fits$reg_b$coef), unname(stats::coef(ref)[terms]), tolerance = 1e-8)
  expect_equal(unname(ar$fits$reg_b$se), unname(fixest::se(ref)[terms]), tolerance = 1e-8)
  expect_equal(ar$tests$slopes$p, fixest::wald(ref, keep = "^Z:", print = FALSE)$p, tolerance = 1e-8)
  expect_equal(ar$tests$chow$p, fixest::wald(ref, keep = "^Z", print = FALSE)$p, tolerance = 1e-8)
  expect_equal(ar$tests$chow$df1, 3)
})

test_that("tables report adjusted R-squared, the joint tests and optional stars", {
  d <- simulate_attrition()
  ar <- attrition_regressions(d$M, d$Z, d$X)

  group <- attrition_table(ar, type = "group")
  expect_equal(names(group), c("Covariate", "Control group", "Treatment group", "Difference (T - C)", "p-value"))
  expect_true(any(grepl("^Adjusted R", group$Covariate)))
  expect_true(any(startsWith(group$Covariate, "Chow test")))
  expect_false(any(grepl("\\*", group[["Control group"]])))

  starred <- attrition_table(ar, type = "pooled", stars = TRUE)
  expect_true(any(grepl("\\*", starred[[2]])))
  expect_true(grepl("p < 0.10", attrition_table_note(ar, "pooled", stars = TRUE), fixed = TRUE))

  expect_true(startsWith(attrition_table(ar, "group", format = "markdown"), "| Covariate |"))
  expect_true(grepl("\\begin{tabular}", attrition_table(ar, "pooled", format = "latex"), fixed = TRUE))
  expect_error(attrition_regressions(d$M + 1, d$Z, d$X), "0/1")
})

test_that("prognostic_check() fits Y on X among observed controls and compares predicted Y*(0)", {
  d <- simulate_attrition()
  set.seed(7)
  Y <- 1 + 2 * d$X$female + d$X$age + 0.5 * d$Z + rnorm(length(d$Z))
  Y[d$M == 0] <- NA

  pc <- prognostic_check(Y = Y, Z = d$Z, X = d$X)

  ctrl_obs <- d$Z == 0 & !is.na(Y)
  ref <- stats::lm(Y ~ female + age, data = data.frame(Y = Y, d$X)[ctrl_obs, ])
  expect_equal(unname(pc$prognostic$coef), unname(stats::coef(ref)))

  score <- as.numeric(cbind(1, as.matrix(d$X)) %*% stats::coef(ref))
  expect_equal(pc$score$all$est, mean(score[d$Z == 1]) - mean(score[d$Z == 0]))
  obs <- !is.na(Y)
  expect_equal(pc$score$respondents$est, mean(score[obs & d$Z == 1]) - mean(score[obs & d$Z == 0]))
  expect_equal(pc$score$observed_y$est, mean(Y[obs & d$Z == 1]) - mean(Y[obs & d$Z == 0]))

  side <- prognostic_table(pc, type = "side_by_side")
  expect_equal(side$Covariate, c("female", "age", "N", "Adjusted R²", "Joint test"))
  expect_equal(ncol(side), 5)
  score_tab <- prognostic_table(pc, type = "score")
  expect_equal(nrow(score_tab), 3)
  expect_true(grepl("treat b as fixed", prognostic_table_note(pc, "score"), fixed = TRUE))
})

test_that("sparse 0/1 covariates trigger a warning and a caution in the notes", {
  d <- simulate_attrition()
  rare <- integer(length(d$Z))
  rare[which(d$Z == 0)[1:2]] <- 1L          # equals 1 for only 2 control units
  rare[which(d$Z == 1)[1:30]] <- 1L
  X <- cbind(d$X, rare = rare)

  expect_warning(
    ar <- attrition_regressions(d$M, d$Z, X, labels = c(rare = "Rare group")),
    "Rare group equals 1 for only 2 units in the control group"
  )
  expect_equal(ar$sparse$covariate, "rare")
  expect_true(grepl("Caution: Rare group equals 1 for only 2 units", attrition_table_note(ar, "group"), fixed = TRUE))
  expect_true(grepl("Caution:", attrition_table_note(ar, "pooled"), fixed = TRUE))

  # Without sparse covariates there is no warning and no caution.
  expect_silent(ar_ok <- attrition_regressions(d$M, d$Z, d$X))
  expect_equal(nrow(ar_ok$sparse), 0)
  expect_false(grepl("Caution:", attrition_table_note(ar_ok, "group"), fixed = TRUE))

  # The prognostic sample (observed control units) is checked too.
  Y <- rnorm(length(d$Z))
  Y[d$M == 0] <- NA
  expect_warning(pc <- prognostic_check(Y, d$Z, X), "prognostic sample")
  expect_true(grepl("Caution:", prognostic_table_note(pc, "score"), fixed = TRUE))
})
