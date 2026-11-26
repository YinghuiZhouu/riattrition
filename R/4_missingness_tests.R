#----------- Test of monotone missingness ------------------------------------#
#
# M_i(z) = 1 if unit i's outcome is observed under treatment z, and 0 if it is
# missing (Li, Sheng & Yu 2025). Their monotone missingness assumptions are
#   * mp: M_i(1) >= M_i(0) for every unit;
#   * mn: M_i(1) <= M_i(0) for every unit.
# M is observed for every unit, so each assumption is a bounded null hypothesis
# about the effect of Z on M, tau_i = M_i(1) - M_i(0):
#   * mp: tau_i >= 0 for every unit;
#   * mn: tau_i <= 0 for every unit.
#
# Caughey, Dafoe, Li & Miratrix (2023, "Randomization inference beyond the
# sharp null", arXiv:2101.09195) show that such a bounded null can be tested
# with the randomization test of the sharp null tau_i = 0 for every unit, as long
# as the test statistic increases with the individual effects, as the
# difference in means does. The function below follows Xinran Li's RIQITE code
# (pval_imp_control() in R/RI_bound_20220919.R). With M as the outcome, it gives
# the same p-value as
#   pval_bound(Z, M, c = 0, method.list = list(name = "DIM"), alternative = "less")
# for mp ("greater" for mn) when both use the same random assignments Z.perm.
#
# The function also tests sharp missingness, M_i(1) = M_i(0) for every unit
# (Li, Sheng & Yu 2025). This is the sharp null tau_i = 0 for every unit
# itself, and it holds exactly when mp and mn both hold. Evidence against it
# can go either way, so its p-value is two-sided: the share of random
# assignments whose statistic is at least as far from 0 as the observed one.
# RIQITE's pval_bound() has no two-sided version.
#
# A small p-value is evidence against the assumption. A large p-value does not
# show that the assumption holds: the test only sees the net difference in
# response rates, so units that violate the assumption can be offset by units
# with an effect in the other direction.

#' Test of monotone missingness
#'
#' Randomization p-value for H0: M_i(1) >= M_i(0) for every unit (`"mp"`),
#' H0: M_i(1) <= M_i(0) for every unit (`"mn"`) or H0: M_i(1) = M_i(0) for
#' every unit (`"sharp"`), where M_i = 1 if unit i's outcome is observed, under
#' complete randomization.
#'
#' @param Z Treatment assignment (0/1 vector).
#' @param Y Outcome vector; `NA` marks a missing outcome.
#' @param assumption `"mp"`, `"mn"` or `"sharp"`.
#' @param Z.perm Optional matrix of random assignments with one row per unit
#'   and one column per assignment. If `NULL`, `nperm` complete
#'   randomizations are drawn.
#' @param nperm Number of random assignments drawn when `Z.perm` is `NULL`.
#'
#' @return The p-value.
#' @export
test_monotone_missingness <- function(Z, Y, assumption = "mp", Z.perm = NULL, nperm = 10^4){
  M = as.numeric(!is.na(Y))   # outcome of this test: 1 = observed, 0 = missing
  m = sum(Z)                  # number of treated units
  n = length(Z)

  # Random assignments, one per column: each treats m of the n units at random
  # (complete randomization). assign_CRE() is in 0_function_sharp_null.R.
  if(is.null(Z.perm)){
    Z.perm = assign_CRE(n, m, nperm)
  }

  # mp and mn are bounded nulls, but we compute the p-value as if tau_i = 0
  # for every unit, which is valid for both (see the top of this file). For
  # sharp, tau_i = 0 for every unit is the hypothesis itself. In all cases
  # every unit keeps its M whatever its treatment, so the loop uses the same M.
  nperm = ncol(Z.perm)
  stat.perm = rep(NA, nperm)
  for(iter in 1:nperm){
    Z.iter = Z.perm[, iter]
    stat.perm[iter] = mean(M[Z.iter == 1]) - mean(M[Z.iter == 0])
  }
  stat.obs = mean(M[Z == 1]) - mean(M[Z == 0])

  # p-value: the share of random assignments whose statistic is at least as
  # extreme as the observed one. Adding 1e-9 makes ties count despite rounding.
  if(assumption == "mp"){
    # Evidence against mp: a low response rate among treated units.
    pval = mean(stat.perm <= stat.obs + 1e-9)
  } else if(assumption == "mn"){
    # Evidence against mn: a high response rate among treated units.
    pval = mean(stat.perm >= stat.obs - 1e-9)
  } else if(assumption == "sharp"){
    # Evidence against sharp: a gap in response rates in either direction.
    pval = mean(abs(stat.perm) >= abs(stat.obs) - 1e-9)
  } else {
    stop(paste0("Unrecognized assumption '", assumption, "'"))
  }
  return(pval)
}


#----------- Test of missing at random ---------------------------------------#
#
# Missing at random (missing = "random" in ri_sharp_attrition()) is Assumption 4
# of Li, Sheng & Yu (2025). Under complete randomization it implies that the
# units with an observed outcome (M = 1) form a completely randomized
# experiment among themselves (their equation (19) and Theorem 5), so their
# treated and control groups should be balanced on pretreatment covariates X.
# The function keeps these units and tests this balance: the statistic is the
# F statistic of the regression of Z on X, and the p-value comes from
# re-randomizing Z among these units. A small p-value is evidence against
# missing at random; a large one means only that we fail to reject it. The test
# cannot detect missingness that depends on the outcome, on characteristics
# that are not in X, or on X but not on the group.

# F statistic of the OLS regression of z on the columns of D (an intercept and
# the covariates), for H0: all coefficients on the covariates are 0. It
# compares the residual sum of squares with and without the covariates.
.balance_F_stat <- function(D, z){
  n = nrow(D)
  k = ncol(D)
  b = solve(crossprod(D), crossprod(D, z))   # OLS coefficients
  rss = sum((z - D %*% b)^2)                  # residual sum of squares
  tss = sum((z - mean(z))^2)                  # the same with the intercept only
  stat = ((tss - rss) / (k - 1)) / (rss / (n - k))
  return(stat)
}

#' Test of missing at random
#'
#' Randomization p-value of a balance test that checks missing at random
#' (Assumption 4 of Li, Sheng and Yu 2025): among the units with an observed
#' outcome, are the covariates `X` balanced between the treatment and the
#' control group?
#'
#' @param Z Treatment assignment (0/1 vector).
#' @param Y Outcome vector; `NA` marks a missing outcome.
#' @param X Data frame (or matrix) of numeric pretreatment covariates, one row
#'   per unit, in the same order as `Z` and `Y`. Units with a missing covariate
#'   are left out.
#' @param nperm Number of re-randomizations of Z among the units with an
#'   observed outcome.
#'
#' @return The p-value.
#' @export
test_missing_at_random <- function(Z, Y, X, nperm = 10^4){
  M = as.numeric(!is.na(Y))   # 1 = observed, 0 = missing
  X = as.data.frame(X)
  if(nrow(X) != length(Z)){
    stop("X must have one row per unit, in the same order as Z and Y.")
  }
  if(!all(sapply(X, is.numeric))){
    stop("X must contain numeric covariates only (code a categorical covariate as 0/1 columns).")
  }

  # Keep the units with an observed outcome (M = 1). Units with a missing
  # covariate are also left out.
  dat = stats::na.omit(data.frame(Z = Z, X)[M == 1, , drop = FALSE])
  if(length(unique(dat$Z)) < 2){
    stop("There must be units with an observed outcome in both groups.")
  }

  # Balance test. Columns of the regression of Z on the covariates: an
  # intercept and X ("." in the formula stands for all the covariates). As in
  # .ar_fit(), columns that lm() cannot estimate (NA coefficient, e.g. a
  # covariate that does not vary) are left out.
  fit = stats::lm(Z ~ ., data = dat)
  D = stats::model.matrix(fit)[, !is.na(stats::coef(fit)), drop = FALSE]
  if(ncol(D) < 2){
    stop("No covariate varies among the units with an observed outcome.")
  }
  stat.obs = .balance_F_stat(D, dat$Z)

  # Re-randomize Z among these units: sample() puts their treatment labels in
  # a random order, so each draw treats the same number of them, chosen at
  # random (complete randomization). The columns D stay the same.
  stat.perm = rep(NA, nperm)
  for(iter in 1:nperm){
    Z.iter = sample(dat$Z)
    stat.perm[iter] = .balance_F_stat(D, Z.iter)
  }

  # p-value: the share of re-randomizations whose statistic is at least the
  # observed one. Adding 1e-9 makes ties count despite rounding.
  pval = mean(stat.perm >= stat.obs - 1e-9)
  return(pval)
}
