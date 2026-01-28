#----------- Functions for Testing Fisher Sharp Null Hypothesis ---------------#

#----------------- Helper Functions for pval_sharp() --------------------------#

#' Calculate the rank score
#'
#' `rank_score()` calculates the scores of n units (ordered from low to high score)
#' given a specified test statistic, not allowing for ties.
#' @param n A positive integer that specifies the number of units.
#' @param method.list A list that specifies the choice of the rank-sum test statistic.
#' * `list(name = "Wilcoxon")` (the default): the Wilcoxon rank-sum statistic.
#' * `list(name = "Stephenson", s = 10)`: the Stephenson rank-sum statistic with parameter s = 10.
#'
#' @return An \eqn{n \times 1} vector.
#' @export
#'
#' @examples
#' rank_score(n = 5, method.list = list(name = "Wilcoxon"))
#' # c(0.2, 0.4, 0.6, 0.8, 1.0)
#' rank_score(n = 5, method.list = list(name = "Stephenson", s = 2))
#' # c(0.00, 0.25, 0.50, 0.75, 1.00)
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
#' @return A scalar.
#' @export
#'
#' @examples
#' Z <- c(1, 0, 0, 1)
#' Y <- c(1:4)
#' rank_treat_relative_to_control(Z, Y)
#' # c(0, 2)
rank_treat_relative_to_control <- function(Z, Y) {
  r <- rank(Y, ties.method = "first")
  r1 <- rank(r[Z == 1])
  r10 <- r[Z == 1] - r1 # n1 x 1 vector
  return(r10)
}

#' Calculate the rank-based test statistic
#'
#' `test_stat()` calculates two classes of rank-based test statistics.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Realized outcome (\eqn{n \times 1} vector, no `NA`s).
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @return A scalar.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0)
#' Y <- c(2, 4, 3, 1, 0)
#' test_stat(Z, Y, MWU = FALSE, list(name = "Wilcoxon")) # 2.4
#' test_stat(Z, Y, MWU = FALSE, list(name = "Stephenson", s = 2)) # 2.25
#' test_stat(Z, Y, MWU = TRUE, list(name = "Wilcoxon")) # 6
#' test_stat(Z, Y, MWU = TRUE, list(name = "Polynomial", s = 2)) # 6
test_stat <- function(Z, Y, MWU = FALSE, method.list = list(name = "Wilcoxon")) {
  # Class 1: rank-sum statistics
  if (!MWU) {
    if (method.list$name %in% c("Wilcoxon", "Stephenson")) {
      n <- length(Y)
      score <- rank_score(n, method.list = method.list)
      stat <- sum(score[rank(Y, ties.method = "first")[Z == 1]])
      return(stat)
    }
  }

  # Class 2: generalized Mann-Whitney U statistics
  if (MWU) {
    r10 <- rank_treat_relative_to_control(Z, Y)

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

    return(sum(phi(r10)))
  }
}

#' Generate matrix of complete randomization assignments
#'
#' Draw multiple assignments from a completely randomized experiment
#'
#' @param n A positive integer that specifies the number of units.
#' @param m A positive integer that specifies the number of treated units.
#' @param nperm A positive integer representing the number of permutations.
#'
#' @return An \eqn{n \times nperm} matrix.
#' @export
#'
#' @examples
#' assign_CRE(n = 5, m = 3, nperm = 10^4)
assign_CRE <- function(n, m, nperm) {
  Z.perm <- matrix(0, nrow = n, ncol = nperm)
  for (iter in 1:nperm) {
    Z.perm[sample(c(1:n), m, replace = FALSE), iter] <- 1
  }
  return(Z.perm)
}

#' Generate randomization distribution of the rank-based test statistic
#'
#' Generate the null distribution of the given rank-based test statistic
#' for an experiment with m treated out of n units.
#' @param n A positive integer that specifies the number of units.
#' @param m A positive integer that specifies the number of treated units.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose m).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return An \eqn{nperm \times 1} vector.
#' @export
#'
#' @examples
#' null_dist(
#'   n = 5, m = 3, MWU = FALSE, method.list = list(name = "Wilcoxon"),
#'   Z.perm = NULL, nperm = 10^4
#' )
#' null_dist(
#'   n = 5, m = 3, MWU = TRUE, method.list = list(name = "Wilcoxon"),
#'   Z.perm = NULL, nperm = 10^4
#' )
null_dist <- function(n, m, MWU = FALSE, method.list = list(name = "Wilcoxon"),
                      Z.perm = NULL, nperm = 10^4) {
  # Class 1: rank-sum statistics
  if (!MWU) {
    # Calculate the score
    score <- rank_score(n, method.list)

    # Generate the Z.perm matrix
    if (is.null(Z.perm)) {
      Z.perm <- assign_CRE(n, m, nperm)
    }
    nperm <- ncol(Z.perm)

    stat.null <- rep(NA, nperm)
    for (iter in 1:nperm) {
      stat.null[iter] <- sum(score[Z.perm[, iter] == 1])
    }

    return(stat.null)
  }

  # Class 2: generalized Mann-Whitney U statistics
  if (MWU) {
    r <- c(1:n)

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

    # Generate the Z.perm matrix
    if (is.null(Z.perm)) {
      Z.perm <- assign_CRE(n, m, nperm)
    }
    nperm <- ncol(Z.perm)

    stat.null <- rep(NA, nperm)
    for (iter in 1:nperm) {
      Z <- Z.perm[, iter]
      r1 <- rank(r[Z == 1])
      # rank of treated units relative to control units
      r10 <- r[Z == 1] - r1
      stat.null[iter] <- sum(phi(r10))
    }

    return(stat.null)
  }
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' @description `pval_sharp_control()` obtains the p-value for testing the sharp null hypothesis \eqn{H_0: \tau = c}.
#'
#' Approach 1: use the worst-case imputed control potential outcomes \eqn{Y(0)} to calculate the test statistic.
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-sum statistic (n choose n1).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose n1).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(NA, 8, 6, 3, 4, NA)
#' pval_sharp_control(
#'   Z, Y, c = 0, MWU = FALSE,
#'   method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4
#' )
pval_sharp_control <- function(Z, Y, c = 0, MWU = FALSE,
                               method.list = list(name = "Wilcoxon"),
                               stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  # Impute the potential outcome if null hypothesis H0: tau = c is true
  Y1.imp <- Y + (1 - Z) * c
  Y0.imp <- Y - Z * c

  # Impute the worst-case control potential outcome Y0
  Y0.imp[Z == 1 & is.na(Y)] <- -Inf
  Y0.imp[Z == 0 & is.na(Y)] <- Inf

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()

  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, MWU = MWU, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z, Y = Y0.imp, MWU = MWU, method.list = method.list)

  pval <- mean(stat.null >= stat.obs)
  return(pval)
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' @description `pval_sharp_composite()` obtains the p-value for testing the sharp null hypothesis \eqn{H_0: \tau = c}.
#'
#' Approach 2: use the worst-case imputed composite control potential outcome
#' \eqn{Y_0M_0 + b(1-M_0)} to calculate the test statistic.
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param b A scalar specifying the composite control potential outcome.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' * `missing = "sharp"`: sharp missing mechanism (\eqn{M_1 = M_0}).
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-sum statistic (n choose n1).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose n1).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(NA, 8, 6, 3, 4, NA)
#' pval_sharp_composite(
#'   Z, Y, c = 0, b = Inf, missing = "general", MWU = FALSE,
#'   method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4
#' )
pval_sharp_composite <- function(Z, Y, c = 0, b, missing = "general", MWU = FALSE,
                                 method.list = list(name = "Wilcoxon"),
                                 stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  # Impute the potential outcome if null hypothesis H0: tau = c is true
  Y1.imp <- Y + (1 - Z) * c
  Y0.imp <- Y - Z * c

  # Impute the worst-case composite control potential outcome Y0 * M0 + b * (1 - M0)
  Y0.com.imp <- rep(NA, n)
  # General missing mechanism (including random, threshold missingness pattern)
  if (missing == "general") {
    Y0.com.imp[Z == 1 & !is.na(Y)] <- pmin(Y0.imp[Z == 1 & !is.na(Y)], b) # min{Yi - c, b}
    Y0.com.imp[Z == 1 & is.na(Y)] <- -Inf
    Y0.com.imp[Z == 0 & !is.na(Y)] <- Y0.imp[Z == 0 & !is.na(Y)]
    Y0.com.imp[Z == 0 & is.na(Y)] <- b
  }
  # Monotone missing mechanism (M1 >= M0)
  else if (missing == "mp") {
    Y0.com.imp[Z == 1 & !is.na(Y)] <- pmin(Y0.imp[Z == 1 & !is.na(Y)], b) # min{Yi - c, b}
    Y0.com.imp[Z == 1 & is.na(Y)] <- b
    Y0.com.imp[Z == 0 & !is.na(Y)] <- Y0.imp[Z == 0 & !is.na(Y)]
    Y0.com.imp[Z == 0 & is.na(Y)] <- b
  }
  # Monotone missing mechanism (M1 <= M0)
  else if (missing == "mn") {
    Y0.com.imp[Z == 1 & !is.na(Y)] <- Y0.imp[Z == 1 & !is.na(Y)]
    Y0.com.imp[Z == 1 & is.na(Y)] <- -Inf
    Y0.com.imp[Z == 0 & !is.na(Y)] <- Y0.imp[Z == 0 & !is.na(Y)]
    Y0.com.imp[Z == 0 & is.na(Y)] <- b
  }
  # Sharp missing mechanism (M1 == M0)
  else if (missing == "sharp") {
    Y0.com.imp[Z == 1 & !is.na(Y)] <- Y0.imp[Z == 1 & !is.na(Y)]
    Y0.com.imp[Z == 1 & is.na(Y)] <- b
    Y0.com.imp[Z == 0 & !is.na(Y)] <- Y0.imp[Z == 0 & !is.na(Y)]
    Y0.com.imp[Z == 0 & is.na(Y)] <- b
  }

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, MWU = MWU, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z, Y = Y0.com.imp, MWU = MWU, method.list = method.list)

  pval <- mean(stat.null >= stat.obs)
  return(pval)
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' @description `pval_sharp_naive()` obtains the p-value for testing the sharp null hypothesis \eqn{H_0: \tau = c}.
#'
#' Naive approach: drop missing data and set \eqn{n_{11}=\sum_{i=1}^n Z_i*M_i} treated
#' and \eqn{n_{01}=\sum_{i=1}^n (1-Z_i)*M_i} control.
#' Implement CRE for only units without missing outcomes,
#' and perform usual randomization test with \eqn{n_{11}} treated and \eqn{n_{01}} control.
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-sum statistic (n.obs choose n1.obs).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n.obs \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n.obs choose n1.obs).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(NA, 8, 6, 3, 4, NA)
#' pval_sharp_naive(
#'   Z, Y, c = 0, MWU = FALSE,
#'   method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4
#' )
pval_sharp_naive <- function(Z, Y, c = 0, MWU = FALSE,
                             method.list = list(name = "Wilcoxon"),
                             stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  # Impute the potential outcome if null hypothesis H0: tau = c is true
  Y1.imp <- Y + (1 - Z) * c
  Y0.imp <- Y - Z * c

  # Drop missing data
  Y1.imp.obs <- Y1.imp[!is.na(Y)]
  Y0.imp.obs <- Y0.imp[!is.na(Y)]
  Z.obs <- Z[!is.na(Y)]

  n.obs <- length(Z.obs) # n11 + n01
  n1.obs <- sum(Z.obs) # n11

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n.obs, n1.obs, MWU = MWU, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z.obs, Y = Y0.imp.obs, MWU = MWU, method.list = method.list)

  pval <- mean(stat.null >= stat.obs)
  return(pval)
}

#------------- Helper Functions for pval_sharp_twostep() ----------------------#

#' Step 1 of two-step method: construct confidence intervals of hypergeometric parameter
#'
#' @description Generate the lower limit of \eqn{1 - \alpha} CP interval.
#' Consider hypergeometric with \eqn{N} units, \eqn{M} of them are good (\eqn{N - M} are bad),
#' select \eqn{n} from \eqn{N} units, there are \eqn{x} good units in the sample.
#' P(X = x) = dhyper(x, m, N - m, n) for max(0, n + M - N) <= x <= min(n, M)
#' Here, \eqn{n} and \eqn{N} are known, and we want to estimate \eqn{M} based on \eqn{x}.
#' Source: "Exact Optimal Confidence Intervals for Hypergeometric Parameters"
#' https://doi.org/10.1080/01621459.2014.966191
#' @param N Number of total units.
#' @param x Number of observed good units.
#' @param n Number of draws.
#' @param alpha Confidence parameter.
#'
#' @return The lower limit of \eqn{1 - \alpha} CP interval.
#' @export
#'
#' @examples
#' lci(N = 100, x = 50, n = 50, alpha = 0.1)
lci <- function(N, x, n, alpha) {
  kk <- 1:length(x)
  for (i in kk) {
    if (x[i] < 0.5) {
      kk[i] <- 0
    } else {
      aa <- 0:N
      bb <- aa + 1
      bb[2:(N + 1)] <- stats::phyper(x[i] - 1, aa[2:(N + 1)] - 1, N - aa[2:(N + 1)] + 1, n)
      dd <- cbind(aa, bb)
      dd <- dd[which(dd[, 2] >= 1 - alpha), ]
      if (length(dd) == 2) {
        kk[i] <- dd[1]
      } else {
        kk[i] <- max(dd[, 1])
      }
    }
  }
  return(kk)
}

#' Step 1 of two-step method: construct confidence intervals of hypergeometric parameter
#'
#' @description Generate the upper limit of \eqn{1 - \alpha} CP interval.
#' Consider hypergeometric with \eqn{N} units, \eqn{M} of them are good (\eqn{N - M} are bad),
#' select \eqn{n} from \eqn{N} units, there are \eqn{x} good units in the sample.
#' P(X = x) = dhyper(x, m, N - m, n) for max(0, n + M - N) <= x <= min(n, M)
#' Here, \eqn{n} and \eqn{N} are known, and we want to estimate \eqn{M} based on \eqn{x}.
#' Source: "Exact Optimal Confidence Intervals for Hypergeometric Parameters"
#' https://doi.org/10.1080/01621459.2014.966191
#' @param N Number of total units.
#' @param x Number of observed good units.
#' @param n Number of draws.
#' @param alpha Confidence parameter.
#'
#' @return The upper limit of \eqn{1 - \alpha} CP interval.
#' @export
#'
#' @examples
#' uci(N = 100, x = 50, n = 50, alpha = 0.1)
uci <- function(N, x, n, alpha) {
  return(N - lci(N, n - x, n, alpha))
}

#' Indicator function for pairwise comparison
#' @param i Index for x.
#' @param j Index for y.
#' @param x A scalar.
#' @param y A scalar.
#'
#' @return The indicator function for pairwise comparison.
#' @export
psi <- function(i, j, x, y) {
  return(as.integer(x > y) + as.integer(x == y) * as.integer(i >= j))
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' @description `pval_sharp_g_twostep()` obtains the p-value for testing the sharp null hypothesis \eqn{H_0: \tau = c}.
#'
#' Approach 3: use the two-step procedure and the worst-case imputed composite control potential outcome
#' to calculate the test statistic (b = Inf).
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-sum statistic (n choose n1).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose n1).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param beta A real number that belongs to \eqn{[0, \alpha]}.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
pval_sharp_g_twostep <- function(Z, Y, c = 0, method.list = list(name = "Wilcoxon"),
                                 stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta) {
  n <- length(Z) # number of observations
  n.obs <- length(Z[!is.na(Y)]) # number of observed units (n11 + n01)
  n1 <- sum(Z) # number of treated units
  n11 <- sum(Z[!is.na(Y)]) # number of observed treated units
  n01 <- n.obs - n11 # number of observed control units

  # Step 1: construct confidence intervals of hypergeometric parameter
  Mhat <- uci(N = n, x = n01, n = n - n1, alpha = beta)
  mbar <- max(Mhat - n01, 0)

  # Step 2
  M <- as.numeric(!is.na(Y))
  Y.comp <- ifelse(M == 1, Y, Inf) # Y_j * M_j + Inf * (1 - M_j)
  ind.treat.obs <- which(Z == 1 & M == 1)
  ind.treat.miss <- which(Z == 1 & M == 0)

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

  A <- rep(NA, n)
  B <- rep(NA, n)
  for (i in ind.treat.obs) {
    A_i <- 0
    B_i <- 0
    for (j in 1:n) {
      if (Z[j] == 0) {
        A_i <- A_i + psi(i, j, Inf, Y.comp[j])
        B_i <- B_i + psi(i, j, Y[i] - c, Y.comp[j])
      }
    }
    A[i] <- A_i
    B[i] <- B_i
  }
  D <- phi(A) - phi(B) # n x 1 vector (where n11 elements are non NAs)

  A <- rep(NA, n)
  C <- rep(NA, n)
  for (i in ind.treat.miss) {
    A_i <- 0
    C_i <- 0
    for (j in 1:n) {
      if (Z[j] == 0) {
        A_i <- A_i + psi(i, j, Inf, Y.comp[j])
        C_i <- C_i + psi(i, j, -Inf, Y.comp[j])
      }
    }
    A[i] <- A_i
    C[i] <- C_i
  }
  E <- phi(A) - phi(C) # n x 1 vector (where n10 elements are non NAs)
  F <- ifelse(is.na(D), E, D) # n x 1 vector (where n1 elements are non NAs)

  r <- rank(F, ties.method = "first") # n x 1 vector
  ind.sort <- sort.int(r, index.return = TRUE)$ix # n x 1 vector
  ind.sort.treat <- ind.sort[Z[ind.sort] == 1] # n1 x 1 vector

  # Impute the worst-case composite control potential outcome
  Y0.com.imp <- rep(NA, n)
  Y0.com.imp[intersect(which((Z == 1) & (M == 1)), ind.sort.treat[seq_len(n1 - mbar)])] <- Inf
  Y0.com.imp[intersect(which((Z == 1) & (M == 1)), ind.sort.treat[(n1 - mbar + 1):n1])] <- Y[intersect(which((Z == 1) & (M == 1)), ind.sort.treat[(n1 - mbar + 1):n1])] - c
  Y0.com.imp[intersect(which((Z == 1) & (M == 0)), ind.sort.treat[seq_len(n1 - mbar)])] <- Inf
  Y0.com.imp[intersect(which((Z == 1) & (M == 0)), ind.sort.treat[(n1 - mbar + 1):n1])] <- -Inf
  Y0.com.imp[Z == 0 & !is.na(Y)] <- Y[Z == 0 & !is.na(Y)]
  Y0.com.imp[Z == 0 & is.na(Y)] <- Inf

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, MWU = TRUE, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z, Y = Y0.com.imp, MWU = TRUE, method.list = method.list)

  pval <- mean(stat.null >= stat.obs) + beta
  return(pval)
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' @description `pval_sharp_mp_twostep()` obtains the p-value for testing the sharp null hypothesis \eqn{H_0: \tau = c}.
#'
#' Approach 3: use the two-step procedure and the worst-case imputed composite control potential outcome
#' to calculate the test statistic (b = Inf).
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-sum statistic (n choose n1).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose n1).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param beta A real number that belongs to \eqn{[0, \alpha]}.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
pval_sharp_mp_twostep <- function(Z, Y, c = 0, method.list = list(name = "Wilcoxon"),
                                  stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta) {
  n <- length(Z) # number of observations
  n.obs <- length(Z[!is.na(Y)]) # number of observed units (n11 + n01)
  n1 <- sum(Z) # number of treated units
  n11 <- sum(Z[!is.na(Y)]) # number of observed treated units
  n01 <- n.obs - n11 # number of observed control units

  # Step 1: construct confidence intervals of hypergeometric parameter
  Mhat <- uci(N = n, x = n01, n = n - n1, alpha = beta)
  mbar <- max(n11 + n01 - Mhat, 0)

  # Step 2
  M <- as.numeric(!is.na(Y))
  Y.comp <- ifelse(M == 1, Y, Inf) # Y_j * M_j + Inf * (1 - M_j)
  ind.treat.obs <- which(Z == 1 & M == 1)

  A <- rep(NA, n)
  B <- rep(NA, n)

  for (i in ind.treat.obs) {
    A_i <- 0
    B_i <- 0
    for (j in 1:n) {
      if (Z[j] == 0) {
        A_i <- A_i + psi(i, j, Inf, Y.comp[j])
        B_i <- B_i + psi(i, j, Y[i] - c, Y.comp[j])
      }
    }
    A[i] <- A_i
    B[i] <- B_i
  }

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

  D <- phi(A) - phi(B) # n x 1 vector (where n11 elements are non NAs)
  r <- rank(D, ties.method = "first") # n x 1 vector
  ind.sort <- sort.int(r, index.return = TRUE)$ix # n x 1 vector
  ind.sort.treat.obs <- ind.sort[Z[ind.sort] == 1 & !is.na(Y[ind.sort])] # n11 x 1 vector

  # Impute the worst-case composite control potential outcome
  Y0.com.imp <- rep(NA, n)
  Y0.com.imp[ind.sort.treat.obs[seq_len(mbar)]] <- Inf
  Y0.com.imp[ind.sort.treat.obs[(mbar + 1):n11]] <- Y[ind.sort.treat.obs[(mbar + 1):n11]] - c
  Y0.com.imp[Z == 1 & is.na(Y)] <- Inf
  Y0.com.imp[Z == 0 & !is.na(Y)] <- Y[Z == 0 & !is.na(Y)]
  Y0.com.imp[Z == 0 & is.na(Y)] <- Inf

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, MWU = TRUE, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z, Y = Y0.com.imp, MWU = TRUE, method.list = method.list)

  pval <- mean(stat.null >= stat.obs) + beta
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
ci_sharp_greater <- function(Z, Y, missing = "general", MWU = FALSE,
                             method.list = list(name = "Wilcoxon"),
                             stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  Z.obs <- Z[!is.na(Y)]
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
      stat.null <- null_dist(n, n1, MWU = MWU,
                             method.list = method.list, Z.perm = Z.perm, nperm = nperm)
    }
    if (missing %in% c("sharp", "random")) {
      stat.null <- null_dist(n.obs, n1.obs, MWU = MWU,
                             method.list = method.list, Z.perm = Z.perm, nperm = nperm)
    }
  }

  f <- function(c) {
    pval <- pval_sharp(Z = Z, Y = Y, c = c, missing = missing, MWU = MWU,
                       method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm)
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
#' @description `pval_sharp()` obtains the p-value for testing the sharp null hypothesis \eqn{H_0: \tau = c}.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' * `missing = "sharp"`: sharp missing mechanism (\eqn{M_1 = M_0}).
#' * `missing = "random"`: missing at random mechanism.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-sum statistic.
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' * if `missing = "general"`, `mp`, `mn`, then n choose n1.
#' * if `missing = "sharp"`, `random`, then n.obs choose n1.obs.
#' @param Z.perm A matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic.
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' * if `missing = "general"`, `mp`, `mn`, then `Z.perm` is an \eqn{n \times nperm} matrix (n choose n1).
#' * if `missing = "sharp"`, `random`, then `Z.perm` is an \eqn{n.obs \times nperm} matrix (n.obs choose n1.obs).
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
pval_sharp <- function(Z, Y, c = 0, missing = "general", MWU = FALSE,
                       method.list = list(name = "Wilcoxon"),
                       stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  if (missing == "general") {
    pval <- pval_sharp_composite(Z = Z, Y = Y, c = c, b = Inf, missing = "general", MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm)
  }
  if (missing == "mp") {
    pval <- pval_sharp_composite(Z = Z, Y = Y, c = c, b = Inf, missing = "mp", MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm)
  }
  if (missing == "mn") {
    pval <- pval_sharp_composite(Z = Z, Y = Y, c = c, b = -Inf, missing = "mn", MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm)
  }
  if (missing %in% c("sharp", "random")) {
    pval <- pval_sharp_naive(Z = Z, Y = Y, c = c, MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm)
  }
  return(pval)
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' @description `pval_sharp_twostep()` obtains the p-value for testing the sharp null hypothesis \eqn{H_0: \tau = c}.
#'
#' Approach 3: use the two-step procedure and the worst-case imputed composite control potential outcome
#' to calculate the test statistic (b = Inf).
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank-sum statistic (n choose n1).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose n1).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param beta A real number that belongs to \eqn{[0, \alpha]}.
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
pval_sharp_twostep <- function(Z, Y, c = 0, missing,
                               method.list = list(name = "Wilcoxon"),
                               stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta) {
  if (missing == "general") {
    return(pval_sharp_g_twostep(Z, Y, c, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, beta = beta))
  }
  if (missing == "mp") {
    return(pval_sharp_mp_twostep(Z, Y, c, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, beta = beta))
  }
  if (missing == "mn") {
    return(pval_sharp_mp_twostep(1 - Z, -Y, c, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, beta = beta))
  }
}

#' Confidence interval assuming constant treatment effect
#'
#' `ci_sharp()` obtains one-sided or two-sided confidence interval for the maximum individual effect.
#'
#' @inheritParams pval_sharp
#'
#' @param alternative A character that specifies the direction of the alternative hypothesis.
#' * `alternative = "greater"`: lower bound
#' * `alternative = "less"`: upper bound
#' * `alternative = "two.sided"`: (lower bound, upper bound)
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A numerical object that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return The one-sided or two-sided confidence interval for the maximum individual effect.
#' @export
ci_sharp <- function(Z, Y, alternative, missing = "general", MWU = FALSE,
                     method.list = list(name = "Wilcoxon"),
                     stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  if (alternative == "greater") {
    ci.lower <- ci_sharp_greater(Z = Z, Y = Y, missing = missing, MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha, tol = tol)
    return(ci.lower)
  }
  if (alternative == "less") {
    ci.upper <- -1 * ci_sharp_greater(Z = Z, Y = -Y, missing = missing, MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha, tol = tol)
    return(ci.upper)
  }
  if (alternative == "two.sided") {
    ci.lower <- ci_sharp_greater(Z = Z, Y = Y, missing = missing, MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha / 2, tol = tol)
    ci.upper <- -1 * ci_sharp_greater(Z = Z, Y = -Y, missing = missing, MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha / 2, tol = tol)
    return(c(ci.lower, ci.upper))
  }
}
