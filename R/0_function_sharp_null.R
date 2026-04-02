#----------- Functions for Testing Fisher Sharp Null Hypothesis ---------------#

#----------------- Helper Functions for pval_sharp() --------------------------#

#' Calculate the rank score
#'
#' `rank_score()` calculates the scores of `n` units (ordered from low to high)
#' given a specified test statistic, not allowing for ties.
#' @param n A positive integer that specifies the number of units.
#' @param method.list A list that specifies the choice of test statistic:
#' * `list(name = "Wilcoxon")` (the default): the Wilcoxon rank-sum statistic.
#' * `list(name = "Stephenson", s = 10)`: the Stephenson rank-sum statistic with
#' parameter `s = 10`, where `s` is a positive integer.
#'
#' @return An \eqn{n \times 1} vector.
#' @examples
#' rank_score(n = 5, method.list = list(name = "Wilcoxon"))
#' # c(0.2, 0.4, 0.6, 0.8, 1.0)
#' rank_score(n = 5, method.list = list(name = "Stephenson", s = 2))
#' # c(0.00, 0.25, 0.50, 0.75, 1.00)
#' @export
rank_score <- function(n, method.list = list(name = "Wilcoxon")) {
  if (method.list$name == "Wilcoxon") {
    score <- c(1:n)
    score <- score / max(score)
    return(score)
  }

  if (method.list$name == "Stephenson") {
    score <- choose(c(1:n) - 1, method.list$s - 1)
    score <- score / max(score)
    return(score)
  }
}

#' Calculate the rank of treated units relative to control units
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Realized outcome (\eqn{n \times 1} vector, no `NA`s).
#'
#' @return An \eqn{n_1 \times 1} vector.
#' @examples
#' Z <- c(1, 0, 0, 1)
#' Y <- c(1:4)
#' rank_treat_relative_to_control(Z, Y)
#' # c(0, 2)
#' @export
rank_treat_relative_to_control <- function(Z, Y) {
  r <- rank(Y, ties.method = "first")
  r1 <- rank(r[Z == 1])
  r10 <- r[Z == 1] - r1 # n1 x 1 vector
  return(r10)
}

#' Calculate the rank of control units relative to treated units
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Realized outcome (\eqn{n \times 1} vector, no `NA`s).
#'
#' @return A \eqn{n_0 \times 1} vector
#' @examples
#' Z <- c(1, 0, 0, 1)
#' Y <- c(1:4)
#' rank_control_relative_to_treat(Z, Y)
#' # c(1, 1)
#' @export
rank_control_relative_to_treat <- function(Z, Y) {
  r <- rank(Y, ties.method = "first")
  r0 <- rank(r[Z == 0])
  r01 <- r[Z == 0] - r0 # n0 x 1 vector
  return(r01)
}

#' Calculate the rank-based test statistic
#'
#' `test_stat()` calculates two classes of rank-based test statistics.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Realized outcome (\eqn{n \times 1} vector, no `NA`s).
#' @param class A string that specifies the class of test statistic:
#' * `class = "RS"` (first class, default): rank-sum statistic.
#' * `class = "MWU+"` (second class): generalized Mann-Whitney U statistic,
#' defined as the sum of relative ranks for treated units.
#' * `class = "MWU-"` (second class): generalized Mann-Whitney U statistic,
#' defined as the negative sum of relative ranks for control units.
#' @param method.list A list that specifies the choice of test statistic:
#' * if `class = "RS"`, `method.list` can be `list(name = "Wilcoxon")` or
#' `list(name = "Stephenson", s = 10)`.
#' * if `class = "MWU+"` or `class = "MWU-"`, `method.list` can be
#' `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#'
#' @return A scalar.
#' @examples
#' Z <- c(1, 1, 1, 0, 0)
#' Y <- c(2, 4, 3, 1, 0)
#' test_stat(Z, Y, class = "RS", list(name = "Wilcoxon")) # 2.4
#' test_stat(Z, Y, class = "RS", list(name = "Stephenson", s = 2)) # 2.25
#' test_stat(Z, Y, class = "MWU+", list(name = "Wilcoxon")) # 6
#' test_stat(Z, Y, class = "MWU+", list(name = "Polynomial", s = 2)) # 6
#' test_stat(Z, Y, class = "MWU-", list(name = "Wilcoxon")) # 0
#' test_stat(Z, Y, class = "MWU-", list(name = "Polynomial", s = 2)) # 0
#' @export
test_stat <- function(Z, Y, class = "RS", method.list = list(name = "Wilcoxon")) {
  # Class 1: rank-sum statistics
  if (class == "RS") {
    if (method.list$name %in% c("Wilcoxon", "Stephenson")) {
      n <- length(Y)
      score <- rank_score(n, method.list = method.list)
      stat <- sum(score[rank(Y, ties.method = "first")[Z == 1]])
      return(stat)
    }
  }

  # Class 2: generalized Mann-Whitney U statistics
  if (class %in% c("MWU+", "MWU-")) {
    # rank transformation
    if (method.list$name == "Wilcoxon") {
      phi <- function(x) {
        x
      }
    }
    if (method.list$name == "Polynomial") {
      phi <- function(x) {
        x^(method.list$s - 1)
      }
    }

    # sum of relative ranks for treated units
    if (class == "MWU+") {
      r10 <- rank_treat_relative_to_control(Z, Y)
      return(sum(phi(r10)))
    }

    # negative sum of relative ranks for control units
    if (class == "MWU-") {
      r01 <- rank_control_relative_to_treat(Z, Y)
      return(-sum(phi(r01)))
    }
  }
}

#' Generate treatment assignments from a completely randomized experiment
#'
#' `assign_CRE()` draws `nperm` independent assignments, each assigning
#' exactly `m` of `n` units to treatment.
#'
#' @param n A positive integer that specifies the number of units.
#' @param m A positive integer that specifies the number of treated units.
#' @param nperm A positive integer that specifies the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return An \eqn{n \times nperm} matrix.
#' @examples
#' assign_CRE(n = 5, m = 3, nperm = 10^4)
#' @export
assign_CRE <- function(n, m, nperm) {
  Z.perm <- matrix(0, nrow = n, ncol = nperm)
  for (iter in 1:nperm) {
    Z.perm[sample(c(1:n), m, replace = FALSE), iter] <- 1
  }
  return(Z.perm)
}

#' Generate randomization distribution of the rank-based test statistic
#'
#' `null_dist()` generates the null distribution of the given rank-based test
#' statistic for an experiment with `m` treated out of `n` units.
#'
#' @param n A positive integer that specifies the number of units.
#' @param m A positive integer that specifies the number of treated units.
#' @param class A string that specifies the class of test statistic:
#' * `class = "RS"` (first class, default): rank-sum statistic.
#' * `class = "MWU+"` (second class): generalized Mann-Whitney U statistic,
#' defined as the sum of relative ranks for treated units.
#' * `class = "MWU-"` (second class): generalized Mann-Whitney U statistic,
#' defined as the negative sum of relative ranks for control units.
#' @param method.list A list that specifies the choice of test statistic:
#' * if `class = "RS"`, `method.list` can be `list(name = "Wilcoxon")` or
#' `list(name = "Stephenson", s = 10)`.
#' * if `class = "MWU+"` or `class = "MWU-"`, `method.list` can be
#' `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permuted assignments
#' for approximating the null distribution of the test statistic
#' (\eqn{n} choose \eqn{m}).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer that specifies the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return An \eqn{nperm \times 1} vector.
#' @examples
#' null_dist(
#'   n = 5, m = 3, class = "RS", method.list = list(name = "Wilcoxon"),
#'   Z.perm = NULL, nperm = 10^4
#' )
#' null_dist(
#'   n = 5, m = 3, class = "MWU+", method.list = list(name = "Wilcoxon"),
#'   Z.perm = NULL, nperm = 10^4
#' )
#' null_dist(
#'   n = 5, m = 3, class = "MWU-", method.list = list(name = "Wilcoxon"),
#'   Z.perm = NULL, nperm = 10^4
#' )
#' @export
null_dist <- function(n, m, class = "RS", method.list = list(name = "Wilcoxon"),
                      Z.perm = NULL, nperm = 10^4) {
  # Generate the Z.perm matrix
  if (is.null(Z.perm)) {
    Z.perm <- assign_CRE(n, m, nperm)
  }
  nperm <- ncol(Z.perm)
  stat.null <- rep(NA, nperm)

  # Class 1: rank-sum statistics
  if (class == "RS") {
    # Calculate the score
    score <- rank_score(n, method.list)
    for (iter in 1:nperm) {
      stat.null[iter] <- sum(score[Z.perm[, iter] == 1])
    }

    return(stat.null)
  }

  # Class 2: generalized Mann-Whitney U statistics
  if (class %in% c("MWU+", "MWU-")) {
    if (method.list$name == "Wilcoxon") {
      phi <- function(x) {
        x
      }
    }
    if (method.list$name == "Polynomial") {
      phi <- function(x) {
        x^(method.list$s - 1)
      }
    }

    r <- c(1:n)

    # sum of relative ranks for treated units
    if (class == "MWU+") {
      for (iter in 1:nperm) {
        Z <- Z.perm[, iter]
        r1 <- rank(r[Z == 1])
        # rank of treated units relative to control units
        r10 <- r[Z == 1] - r1
        stat.null[iter] <- sum(phi(r10))
      }
      return(stat.null)
    }

    # negative sum of relative ranks for control units
    if (class == "MWU-") {
      for (iter in 1:nperm) {
        Z <- Z.perm[, iter]
        r0 <- rank(r[Z == 0])
        # rank of control units relative to treated units
        r01 <- r[Z == 0] - r0
        stat.null[iter] <- -sum(phi(r01))
      }
      return(stat.null)
    }
  }
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' `pval_sharp_control()` obtains the p-value for testing the sharp null hypothesis
#' \eqn{H_0: \tau = c}, using the worst-case imputed control potential outcomes \eqn{Y(0)}.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param class A string that specifies the class of test statistic:
#' * `class = "RS"` (first class, default): rank-sum statistic.
#' * `class = "MWU+"` (second class): generalized Mann-Whitney U statistic,
#' defined as the sum of relative ranks for treated units.
#' * `class = "MWU-"` (second class): generalized Mann-Whitney U statistic,
#' defined as the negative sum of relative ranks for control units.
#' @param method.list A list that specifies the choice of test statistic:
#' * if `class = "RS"`, `method.list` can be `list(name = "Wilcoxon")` or
#' `list(name = "Stephenson", s = 10)`.
#' * if `class = "MWU+"` or `class = "MWU-"`, `method.list` can be
#' `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-based statistic
#' (\eqn{n} choose \eqn{n_1}).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permuted assignments
#' for approximating the null distribution of the test statistic
#' (\eqn{n} choose \eqn{n_1}).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer that specifies the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(NA, 8, 6, 3, 4, NA)
#' pval_sharp_control(Z, Y,
#'   c = 0, class = "RS", method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4
#' )
#' @export
pval_sharp_control <- function(Z, Y, c = 0, class = "RS", method.list = list(name = "Wilcoxon"),
                               stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  # Impute the potential outcome if null hypothesis H0: tau = c is true
  Y1.imp <- Y + (1 - Z) * c
  Y0.imp <- Y - Z * c

  # Impute the worst-case control potential outcome Y0
  Y0.imp[Z == 1 & M == 0] <- -Inf
  Y0.imp[Z == 0 & M == 0] <- Inf

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()

  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, class = class, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z, Y = Y0.imp, class = class, method.list = method.list)

  pval <- mean(stat.null >= stat.obs)
  return(pval)
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' `pval_sharp_composite()` obtains the p-value for testing the sharp null hypothesis
#' \eqn{H_0: \tau = c}, using the worst-case imputed composite control potential outcome \eqn{Y_b(0)}.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param b00,b01,b10 Three scalars that specifies the composite control potential outcome.
#' @param missing A string that specifies the missing mechanism:
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 \geq M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 \leq M_0}).
#' * `missing = "sharp"`: sharp missing mechanism (\eqn{M_1 = M_0}).
#' @param class A string that specifies the class of test statistic:
#' * `class = "RS"` (first class, default): rank-sum statistic.
#' * `class = "MWU+"` (second class): generalized Mann-Whitney U statistic,
#' defined as the sum of relative ranks for treated units.
#' * `class = "MWU-"` (second class): generalized Mann-Whitney U statistic,
#' defined as the negative sum of relative ranks for control units.
#' @param method.list A list that specifies the choice of test statistic:
#' * if `class = "RS"`, `method.list` can be `list(name = "Wilcoxon")` or
#' `list(name = "Stephenson", s = 10)`.
#' * if `class = "MWU+"` or `class = "MWU-"`, `method.list` can be
#' `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-based statistic
#' (\eqn{n} choose \eqn{n_1}).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permuted assignments
#' for approximating the null distribution of the test statistic
#' (\eqn{n} choose \eqn{n_1}).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer that specifies the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(NA, 8, 6, 3, 4, NA)
#' pval_sharp_composite(Z, Y,
#'   c = 0, b00 = 0, b01 = Inf, b10 = -Inf, missing = "general",
#'   class = "RS", method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4
#' )
#' @export
pval_sharp_composite <- function(Z, Y, c = 0, b00, b01, b10, missing = "general",
                                 class = "RS", method.list = list(name = "Wilcoxon"),
                                 stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  # Impute the potential outcome if null hypothesis H0: tau = c is true
  Y1.imp <- Y + (1 - Z) * c
  Y0.imp <- Y - Z * c

  # Impute the worst-case composite control potential outcome
  # Y0 * M0 * M1 + b00 * (1 - M0) * (1 - M1) + b01 * (1 - M0) * M1 + b10 * M0 * (1 - M1)
  Y0.com.imp <- rep(NA, n)
  # General missing mechanism (including random, threshold missingness pattern)
  if (missing == "general") {
    Y0.com.imp[Z == 1 & M == 1] <- pmin(Y0.imp[Z == 1 & M == 1], b01) # min{Yi - c, b01}
    Y0.com.imp[Z == 1 & M == 0] <- min(b00, b10)
    Y0.com.imp[Z == 0 & M == 1] <- pmax(Y0.imp[Z == 0 & M == 1], b10) # max{Yi, b10}
    Y0.com.imp[Z == 0 & M == 0] <- max(b00, b01)
  }
  # Monotone missing mechanism (M1 >= M0)
  else if (missing == "mp") {
    Y0.com.imp[Z == 1 & M == 1] <- pmin(Y0.imp[Z == 1 & M == 1], b01) # min{Yi - c, b01}
    Y0.com.imp[Z == 1 & M == 0] <- b00
    Y0.com.imp[Z == 0 & M == 1] <- Y0.imp[Z == 0 & M == 1]
    Y0.com.imp[Z == 0 & M == 0] <- max(b00, b01)
  }
  # Monotone missing mechanism (M1 <= M0)
  else if (missing == "mn") {
    Y0.com.imp[Z == 1 & M == 1] <- Y0.imp[Z == 1 & M == 1]
    Y0.com.imp[Z == 1 & M == 0] <- min(b00, b10)
    Y0.com.imp[Z == 0 & M == 1] <- pmax(Y0.imp[Z == 0 & M == 1], b10) # max{Yi, b10}
    Y0.com.imp[Z == 0 & M == 0] <- b00
  }
  # Sharp missing mechanism (M1 == M0)
  else if (missing == "sharp") {
    Y0.com.imp[Z == 1 & M == 1] <- Y0.imp[Z == 1 & M == 1]
    Y0.com.imp[Z == 1 & M == 0] <- b00
    Y0.com.imp[Z == 0 & M == 1] <- Y0.imp[Z == 0 & M == 1]
    Y0.com.imp[Z == 0 & M == 0] <- b00
  }

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, class = class, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z, Y = Y0.com.imp, class = class, method.list = method.list)

  pval <- mean(stat.null >= stat.obs)
  return(pval)
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' `pval_sharp_naive()` obtains the p-value for testing the sharp null hypothesis
#' \eqn{H_0: \tau = c} using naive approach. Implement CRE for only units without
#' missing outcomes, and perform usual randomization test with \eqn{n_{11}} treated
#' and \eqn{n_{01}} control.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param class A string that specifies the class of test statistic:
#' * `class = "RS"` (first class, default): rank-sum statistic.
#' * `class = "MWU+"` (second class): generalized Mann-Whitney U statistic,
#' defined as the sum of relative ranks for treated units.
#' * `class = "MWU-"` (second class): generalized Mann-Whitney U statistic,
#' defined as the negative sum of relative ranks for control units.
#' @param method.list A list that specifies the choice of test statistic:
#' * if `class = "RS"`, `method.list` can be `list(name = "Wilcoxon")` or
#' `list(name = "Stephenson", s = 10)`.
#' * if `class = "MWU+"` or `class = "MWU-"`, `method.list` can be
#' `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-based statistic
#' (\eqn{n_{11} + n_{01}} choose \eqn{n_{11}}).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{(n_{11} + n_{01}) \times nperm} matrix that specifies the permuted assignments
#' for approximating the null distribution of the test statistic
#' (\eqn{n_{11} + n_{01}} choose \eqn{n_{11}}).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer that specifies the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(NA, 8, 6, 3, 4, NA)
#' pval_sharp_naive(Z, Y,
#'   c = 0, class = "RS", method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4
#' )
#' @export
pval_sharp_naive <- function(Z, Y, c = 0, class = "RS", method.list = list(name = "Wilcoxon"),
                             stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  # Impute the potential outcome if null hypothesis H0: tau = c is true
  Y1.imp <- Y + (1 - Z) * c
  Y0.imp <- Y - Z * c

  # Drop missing data
  Y1.imp.obs <- Y1.imp[M == 1]
  Y0.imp.obs <- Y0.imp[M == 1]
  Z.obs <- Z[M == 1]

  n.obs <- length(Z.obs) # n11 + n01
  n1.obs <- sum(Z.obs) # n11

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n.obs, n1.obs, class = class, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z.obs, Y = Y0.imp.obs, class = class, method.list = method.list)

  pval <- mean(stat.null >= stat.obs)
  return(pval)
}

#--------------------- Helper Functions for ci_sharp() ------------------------#

#' Confidence interval assuming constant treatment effect
#'
#' `ci_sharp_greater()` obtains one-sided confidence interval for the maximum individual effect.
#'
#' @inheritParams pval_sharp
#'
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A numerical object that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return The one-sided confidence interval for the maximum individual effect.
#' @export
ci_sharp_greater <- function(Z, Y, missing = "general", class = "RS", method.list = list(name = "Wilcoxon"),
                             stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  Z.obs <- Z[M == 1]
  n.obs <- length(Z.obs) # n11 + n01
  n1.obs <- sum(Z.obs) # n11

  # Generate Z.perm as matrix of possible treatment assignments with nperm combinations
  if (is.null(Z.perm)) {
    if (missing %in% c("general", "mp", "mn")) {
      Z.perm <- assign_CRE(n, n1, nperm)
    }
    if (missing %in% c("sharp", "random")) {
      Z.perm <- assign_CRE(n.obs, n1.obs, nperm)
    }
  }
  nperm <- ncol(Z.perm)
  # Z.perm: n x nperm or n.obs x nperm matrix in which each column is one permutation of treatment assignment.

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    if (missing %in% c("general", "mp", "mn")) {
      stat.null <- null_dist(n, n1, class = class, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
    }
    if (missing %in% c("sharp", "random")) {
      stat.null <- null_dist(n.obs, n1.obs, class = class, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
    }
  }

  f <- function(c) {
    pval <- pval_sharp(
      Z = Z, Y = Y, c = c, missing = missing, class = class,
      method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm
    )
    return(pval - alpha)
  }

  c.max <- max(Y, na.rm = TRUE) - min(Y, na.rm = TRUE) + tol
  c.min <- -max(Y, na.rm = TRUE) + min(Y, na.rm = TRUE) - tol

  if (f(c.min) >= 0) {
    return(-Inf)
  }

  if (f(c.max) < 0) {
    return(Inf)
  }

  c_sol <- stats::uniroot(f, interval = c(c.min, c.max), extendInt = "upX", tol = tol)$root
  c_sol <- round(c_sol, digits = -log10(tol))

  if (f(c_sol) > 0) {
    while (f(c_sol) > 0) {
      c_sol <- c_sol - tol
    }
    c_sol <- c_sol + tol
  } else {
    while (f(c_sol) <= 0) {
      c_sol <- c_sol + tol
    }
  }
  return(c_sol)
}

#------------------------------ Main Functions --------------------------------#

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' `pval_sharp()` obtains the p-value for testing the sharp null hypothesis \eqn{H_0: \tau = c}.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param missing A string that specifies the missing mechanism:
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 \geq M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 \leq M_0}).
#' * `missing = "sharp"`: sharp missing mechanism (\eqn{M_1 = M_0}).
#' * `missing = "random"`: missing at random mechanism.
#' @param class A string that specifies the class of test statistic:
#' * `class = "RS"` (first class, default): rank-sum statistic.
#' * `class = "MWU+"` (second class): generalized Mann-Whitney U statistic,
#' defined as the sum of relative ranks for treated units.
#' * `class = "MWU-"` (second class): generalized Mann-Whitney U statistic,
#' defined as the negative sum of relative ranks for control units.
#' @param method.list A list that specifies the choice of test statistic:
#' * if `class = "RS"`, `method.list` can be `list(name = "Wilcoxon")` or
#' `list(name = "Stephenson", s = 10)`.
#' * if `class = "MWU+"` or `class = "MWU-"`, `method.list` can be
#' `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-based statistic.
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' * if `missing = "general"`, `mp`, `mn`, then \eqn{n} choose \eqn{n_1}.
#' * if `missing = "sharp"`, `random`, then \eqn{n_{11} + n_{01}} choose \eqn{n_{11}}.
#' @param Z.perm A matrix that specifies the permuted assignments
#' for approximating the null distribution of the test statistic.
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' * if `missing = "general"`, `mp`, `mn`, then `Z.perm` is an \eqn{n \times nperm} matrix
#' (\eqn{n} choose \eqn{n_1}).
#' * if `missing = "sharp"`, `random`, then `Z.perm` is an \eqn{(n_{11} + n_{01}) \times nperm} matrix
#' (\eqn{n_{11} + n_{01}} choose \eqn{n_{11}}).
#' @param nperm A positive integer that specifies the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
pval_sharp <- function(Z, Y, c = 0, missing = "general", class = "RS", method.list = list(name = "Wilcoxon"),
                       stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  if (missing == "general") {
    pval <- pval_sharp_composite(
      Z = Z, Y = Y, c = c, b00 = 0, b01 = Inf, b10 = -Inf,
      missing = "general", class = class, method.list = method.list,
      stat.null = stat.null, Z.perm = Z.perm, nperm = nperm
    )
  }
  if (missing == "mp") {
    pval <- pval_sharp_composite(
      Z = Z, Y = Y, c = c, b00 = Inf, b01 = Inf, b10 = 0,
      missing = "mp", class = class, method.list = method.list,
      stat.null = stat.null, Z.perm = Z.perm, nperm = nperm
    )
  }
  if (missing == "mn") {
    pval <- pval_sharp_composite(
      Z = Z, Y = Y, c = c, b00 = -Inf, b01 = 0, b10 = -Inf,
      missing = "mn", class = class, method.list = method.list,
      stat.null = stat.null, Z.perm = Z.perm, nperm = nperm
    )
  }
  if (missing %in% c("sharp", "random")) {
    pval <- pval_sharp_naive(
      Z = Z, Y = Y, c = c, class = class, method.list = method.list,
      stat.null = stat.null, Z.perm = Z.perm, nperm = nperm
    )
  }
  return(pval)
}

#' Confidence interval assuming constant treatment effect
#'
#' `ci_sharp()` obtains one-sided or two-sided confidence interval for the maximum individual effect.
#'
#' @inheritParams pval_sharp
#'
#' @param alternative A string that specifies the direction of the alternative hypothesis:
#' * `alternative = "greater"`: lower bound.
#' * `alternative = "less"`: upper bound.
#' * `alternative = "two.sided"`: (lower bound, upper bound).
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A numerical object that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return The one-sided or two-sided confidence interval for the maximum individual effect.
#' @export
ci_sharp <- function(Z, Y, alternative, missing = "general", class = "RS", method.list = list(name = "Wilcoxon"),
                     stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  if (alternative == "greater") {
    ci.lower <- ci_sharp_greater(
      Z = Z, Y = Y, missing = missing, class = class,
      method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha, tol = tol
    )
    return(ci.lower)
  }
  if (alternative == "less") {
    ci.upper <- -1 * ci_sharp_greater(
      Z = Z, Y = -Y, missing = missing, class = class,
      method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha, tol = tol
    )
    return(ci.upper)
  }
  if (alternative == "two.sided") {
    ci.lower <- ci_sharp_greater(
      Z = Z, Y = Y, missing = missing, class = class,
      method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha / 2, tol = tol
    )
    ci.upper <- -1 * ci_sharp_greater(
      Z = Z, Y = -Y, missing = missing, class = class,
      method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha / 2, tol = tol
    )
    return(c(ci.lower, ci.upper))
  }
}
