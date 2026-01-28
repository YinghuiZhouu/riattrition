#----- Functions for Testing Quantiles of Individual Treatment Effects --------#

#------------------------------------ Basics ----------------------------------#

#' Calculate minimum test statistic value in \eqn{H_{k,c}} with missing outcomes for general/mp/mn missing
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k A positive integer that specifies which quantile of individual effect is of interest (\eqn{1 \leq k \leq n_1}).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#'
#' @return A scalar.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(NA, 8, 6, 3, 4, NA)
#' min_stat_g_mp_mn(
#'   Z, Y, k = 2, c = 0, missing = "general",
#'   MWU = FALSE, method.list = list(name = "Wilcoxon")
#' )
min_stat_g_mp_mn <- function(Z, Y, k, c, missing = "general", MWU = FALSE, method.list = list(name = "Wilcoxon")) {
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units
  n11 <- sum(Z[!is.na(Y)]) # number of observed treated units

  # Sort the observed treated units using the "first" tie method
  r <- rank(Y, ties.method = "first", na.last = FALSE) # n x 1 vector
  ind.sort <- sort.int(r, index.return = TRUE)$ix # n x 1 vector
  ind.sort.treat.obs <- ind.sort[Z[ind.sort] == 1 & !is.na(Y[ind.sort])] # n11 x 1 vector

  # Generate xi vector
  xi <- rep(c, n)
  if (k < n1) {
    xi[ind.sort.treat.obs[(n11 + 1 - min(n11, n1 - k)):n11]] <- max(Y[Z == 1 & !is.na(Y)]) - min(Y[Z == 0 & !is.na(Y)]) + 1
  }

  # Impute the worst-case composite control potential outcome Y_{t,k,c}(0)
  Y0.com.imp <- rep(NA, n)

  # General missing mechanism (including threshold missing mechanism)
  if (missing == "general") {
    b <- Inf
    Y0.com.imp[Z == 1 & !is.na(Y)] <- pmin(Y[Z == 1 & !is.na(Y)] - xi[Z == 1 & !is.na(Y)], b) # min{Y_i - xi_i, b}
    Y0.com.imp[Z == 1 & is.na(Y)] <- -Inf
    Y0.com.imp[Z == 0 & !is.na(Y)] <- Y[Z == 0 & !is.na(Y)]
    Y0.com.imp[Z == 0 & is.na(Y)] <- b
  }
  # Monotone missing mechanism (M1 >= M0)
  else if (missing == "mp") {
    b <- Inf
    Y0.com.imp[Z == 1 & !is.na(Y)] <- pmin(Y[Z == 1 & !is.na(Y)] - xi[Z == 1 & !is.na(Y)], b) # min{Y_i - xi_i, b}
    Y0.com.imp[Z == 1 & is.na(Y)] <- b
    Y0.com.imp[Z == 0 & !is.na(Y)] <- Y[Z == 0 & !is.na(Y)]
    Y0.com.imp[Z == 0 & is.na(Y)] <- b
  }
  # Monotone missing mechanism (M1 <= M0)
  else if (missing == "mn") {
    b <- -Inf
    Y0.com.imp[Z == 1 & !is.na(Y)] <- Y[Z == 1 & !is.na(Y)] - xi[Z == 1 & !is.na(Y)] # Y_i - xi_i
    Y0.com.imp[Z == 1 & is.na(Y)] <- -Inf
    Y0.com.imp[Z == 0 & !is.na(Y)] <- Y[Z == 0 & !is.na(Y)]
    Y0.com.imp[Z == 0 & is.na(Y)] <- b
  }

  stat.min <- test_stat(Z = Z, Y = Y0.com.imp, MWU = MWU, method.list = method.list)
  return(stat.min)
}

#' Calculate minimum test statistic value in \eqn{H_{k,c}} with no missing outcomes
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, no `NA`s).
#' @param k A positive integer that specifies which quantile of individual effect is of interest (\eqn{1 \leq k \leq n_1}).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#'
#' @return A scalar.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(7, 8, 6, 3, 4, 1)
#' min_stat(Z, Y, k = 2, c = 0, MWU = FALSE, method.list = list(name = "Wilcoxon"))
min_stat <- function(Z, Y, k, c, MWU = FALSE, method.list = list(name = "Wilcoxon")) {
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  # Sort the treated units using the "first" tie method
  r <- rank(Y, ties.method = "first") # n x 1 vector
  ind.sort <- sort.int(r, index.return = TRUE)$ix # n x 1 vector
  ind.sort.treat <- ind.sort[Z[ind.sort] == 1] # n1 x 1 vector

  # Generate xi vector
  xi <- rep(c, n)
  if (k < n1) {
    xi[ind.sort.treat[(k + 1):n1]] <- Inf
  }

  stat.min <- test_stat(Z = Z, Y = Y - Z * xi, MWU = MWU, method.list = method.list)
  return(stat.min)
}

#------------------ p-value for testing tau_{(k)} <= c ------------------------#

#' Randomization test for quantiles of individual treatment effects with missing outcomes for general/mp/mn missing
#'
#' @description `pval_quantile_g_mp_mn()` obtains the p-value for testing the null hypothesis
#' \eqn{H_{k,c}: \tau_{(k)} \leq c}, where \eqn{\tau_{(k)}} denotes individual treatment effect at rank \eqn{k}.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k A positive integer that specifies which quantile of individual effect is of interest (\eqn{1 \leq k \leq n_1}).
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n choose n1).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose n1).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified null hypothesis of interest (a scalar).
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(8, 9, NA, 3, 4, NA)
#' pval <- pval_quantile_g_mp_mn(
#'   Z, Y, k = 2, c = 0, missing = "general",
#'   MWU = FALSE, method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4
#' )
pval_quantile_g_mp_mn <- function(Z, Y, k, c, missing = "general", MWU = FALSE,
                                  method.list = list(name = "Wilcoxon"), stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units
  n11 <- sum(Z[!is.na(Y)]) # number of observed treated units

  # Generate null distribution
  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, MWU = MWU, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Generate minimum test statistic value under H_{k,c}
  stat.min <- min_stat_g_mp_mn(Z, Y, k, c, missing = missing, MWU = MWU, method.list = method.list)

  # Calculate p-value
  pval <- mean(stat.null >= stat.min)

  return(pval)
}

#' Randomization test for quantiles of individual treatment effects with missing outcomes for sharp/random missing
#'
#' @description `pval_quantile_s_r()` obtains the p-value for testing the null hypothesis
#' \eqn{H_{k,c}: \tau_{(k)} \leq c}, where \eqn{\tau_{(k)}} denotes individual treatment effect at rank \eqn{k}.
#' Naive approach: drop missing data and set \eqn{n_{11}=\sum_{i=1}^n Z_i*M_i} treated
#' and \eqn{n_{01}=\sum_{i=1}^n (1-Z_i)*M_i} control.
#' Implement CRE for only units without missing outcomes,
#' and perform usual randomization test with \eqn{n_{11}} treated and \eqn{n_{01}} control.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k A positive integer that specifies which quantile of individual effect is of interest (\eqn{1 \leq k \leq n_{11}})
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n.obs choose n1.obs).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n.obs \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n.obs choose n1.obs).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#'
#' @return The p-value for testing the specified null hypothesis of interest (a scalar).
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(8, 9, NA, 3, 4, NA)
#' pval <- pval_quantile_s_r(
#'   Z, Y, k = 2, c = 0, MWU = FALSE, method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4
#' )
pval_quantile_s_r <- function(Z, Y, k, c, MWU = FALSE, method.list = list(name = "Wilcoxon"),
                              stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  # Drop missing data
  Y.obs <- Y[!is.na(Y)]
  Z.obs <- Z[!is.na(Y)]

  n.obs <- length(Z.obs) # n_s = n11 + n01
  n1.obs <- sum(Z.obs) # n11

  # Generate null distribution
  if (is.null(stat.null)) {
    stat.null <- null_dist(n.obs, n1.obs, MWU = MWU, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Generate minimum test statistic value under H_{k,c}
  stat.min <- min_stat(Z = Z.obs, Y = Y.obs, k, c, MWU = MWU, method.list = method.list)

  # p-value
  pval <- mean(stat.null >= stat.min)

  return(pval)
}

#' Randomization test for quantiles of individual treatment effects with missing outcomes
#'
#' @description `pval_quantile()` obtains the p-value for testing the null hypothesis
#' \eqn{H_{k,c}: \tau_{(k)} \leq c}, where \eqn{\tau_{(k)}} denotes individual treatment effect at rank \eqn{k}.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k A positive integer that specifies which quantile of individual effect is of interest
#' * if `missing = "general"` or `missing = "mp"` or `missing = "mn"`, \eqn{1 \leq k \leq n_1}.
#' * if `missing = "sharp"` or `missing = "random"`, \eqn{1 \leq k \leq n_{11}}.
#' @param c A scalar or vector specifying the sharp null hypothesis.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' * `missing = "sharp"`: sharp missing mechanism (\eqn{M_1 = M_0}).
#' * `missing = "random"`: random missing mechanism.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistics.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic.
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
#' @return The p-value for testing the specified null hypothesis of interest (a scalar).
#' @export
pval_quantile <- function(Z, Y, k, c, missing = "general", MWU = FALSE,
                          method.list = list(name = "Wilcoxon"), stat.null = NULL, Z.perm = NULL, nperm = 10^4) {
  if (missing %in% c("general", "mp", "mn")) {
    pval <- pval_quantile_g_mp_mn(Z, Y, k, c, missing = missing, MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm)
  }

  if (missing %in% c("sharp", "random")) {
    pval <- pval_quantile_s_r(Z, Y, k, c, MWU = MWU, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm)
  }

  return(pval)
}

#--------------------------- Confidence Interval ------------------------------#

#--- general/mp/mn missing ---#

#' Calculate one-sided confidence interval (c, Inf) for all quantiles among treated units with missing outcomes for general/mp/mn missing
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k.vec A vector that specifies the quantiles of individual effects of interest.
#' If `k.vec = NULL`, then we consider all quantiles of individual effects \eqn{1 \leq k \leq n_1}.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n choose n1).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose n1).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A scalar that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return An \eqn{L \times 1} vector corresponding to the lower confidence limits of the given \eqn{\tau_{(k)}}.
#' * if `k.vec = NULL`, \eqn{L = n_{1}}.
#' * if `k.vec` is not `NULL`, \eqn{L} is the length of `k.vec`.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(8, 9, NA, 3, 4, NA)
#' c.limit <- conf_quant_treated_g_mp_mn(Z, Y,
#'   k.vec = NULL, missing = "general",
#'   MWU = FALSE, method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)
#' )
conf_quant_treated_g_mp_mn <- function(Z, Y, k.vec = NULL, missing = "general", MWU = FALSE, method.list = list(name = "Wilcoxon"),
                                       stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units
  n11 <- sum(Z[!is.na(Y)]) # number of observed treated units

  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, MWU = MWU, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Find threshold such that p-value <= alpha <=> test statistic > threshold
  thres <- sort(stat.null, decreasing = TRUE)[floor(nperm * alpha) + 1]

  # Calculate range of c
  Y1.max <- max(Y[Z == 1], na.rm = T)
  Y1.min <- min(Y[Z == 1], na.rm = T)
  Y0.max <- max(Y[Z == 0], na.rm = T)
  Y0.min <- min(Y[Z == 0], na.rm = T)
  c.max <- Y1.max - Y0.min + tol
  c.min <- Y1.min - Y0.max - tol

  if (is.null(k.vec)) {
    # Initialize conf interval for all quantiles tau_{(k)}, 1 <= k <= n1
    c.limit <- rep(NA, n1)

    # When k <= n1 - n11, conf interval for tau_{(k)} is often (-Inf, Inf)
    for (k in n1:(n1 - n11)) {
      if (k < n1) {
        c.max <- c.limit[k + 1]
      }
      # define the target fun
      # f > 0 <=> p-value <= alpha
      # f decreases in c, p value increases in c
      f <- function(c) {
        stat.min <- min_stat_g_mp_mn(Z, Y, k, c, missing = missing, MWU = MWU, method.list = method.list)
        return(stat.min - thres)
      }

      if (f(c.min) <= 0) {
        c.sol <- -Inf
      }
      if (f(c.max) > 0) {
        c.sol <- c.max
      }
      if (f(c.min) > 0 & f(c.max) <= 0) {
        c.sol <- stats::uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        # find the min c s.t. p-value > alpha <=> f <= 0
        c.sol <- round(c.sol, digits = -log10(tol))
        if (f(c.sol) <= 0) {
          while (f(c.sol) <= 0) {
            c.sol <- c.sol - tol
          }
          c.sol <- c.sol + tol
        } else {
          while (f(c.sol) > 0) {
            c.sol <- c.sol + tol
          }
        }
      }
      c.limit[k] <- c.sol
    }

    if (n1 - n11 > 1) {
      c.limit[1:(n1 - n11 - 1)] <- c.limit[n1 - n11]
    }

    c.limit[c.limit > (Y1.max - Y0.min) + tol / 2] <- Inf
  } else {
    k.vec.sort <- sort(k.vec, decreasing = FALSE)
    j.max <- length(k.vec.sort)
    j.min <- max(sum(k.vec <= (n1 - n11)), 1)

    # conf interval for tau_{(k)} with k in k.vec
    c.limit <- rep(NA, j.max)

    # when k <= j.min, conf interval for tau_{(k)} is often (-Inf, Inf)
    for (j in j.max:j.min) {
      k <- k.vec.sort[j]
      if (j < j.max) {
        c.max <- c.limit[j + 1]
      }
      # define the target fun
      # f > 0 <=> p-value <= alpha
      # f decreases in c, p value increases in c
      f <- function(c) {
        stat.min <- min_stat_g_mp_mn(Z, Y, k, c, missing = missing, MWU = MWU, method.list = method.list)
        return(stat.min - thres)
      }
      # check whether f(-Inf) = f(c.min) > 0
      if (f(c.min) <= 0) {
        c.sol <- -Inf
      }
      if (f(c.max) > 0) {
        c.sol <- c.max
      }
      if (f(c.min) > 0 & f(c.max) <= 0) {
        c.sol <- stats::uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        # find the min c st p-value > alpha <=> f <= 0
        c.sol <- round(c.sol, digits = -log10(tol))
        if (f(c.sol) <= 0) {
          while (f(c.sol) <= 0) {
            c.sol <- c.sol - tol
          }
          c.sol <- c.sol + tol
        } else {
          while (f(c.sol) > 0) {
            c.sol <- c.sol + tol
          }
        }
      }
      c.limit[j] <- c.sol
    }

    if (j.min > 1) {
      c.limit[1:(j.min - 1)] <- c.limit[j.min]
    }

    c.limit[c.limit > (Y1.max - Y0.min) + tol / 2] <- Inf
  }

  return(c.limit)
}

#' Randomization inference for quantiles of individual treatment effects under general/mp/mn missing mechanism
#'
#' @description `ci_quantile_g_mp_mn()` obtains one-sided confidence intervals
#' for all quantiles of individual treatment effects.
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k.vec A vector that specifies the quantiles of individual effects of interest.
#' If `k.vec = NULL`, then we consider all quantiles of individual effects \eqn{1 \leq k \leq n}.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list.treated A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param method.list.control A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null.treated An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n choose n1).
#' * if `stat.null.treated = NULL` (default), the function will calculate it internally.
#' @param Z.perm.treated An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose n1).
#' * if `Z.perm.treated = NULL` (default), the function will calculate it internally.
#' @param stat.null.control An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n choose n0).
#' * if `stat.null.control = NULL` (default), the function will calculate it internally.
#' @param Z.perm.control An \eqn{n \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n choose n0).
#' * if `Z.perm.control = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A scalar that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return A dataframe with \eqn{L} rows and 3 columns (\eqn{k, lower, upper}).
#' Each row corresponds to the lower and upper confidence limits of the given \eqn{\tau_{(k)}}.
#' * if `k.vec = NULL`, \eqn{L = n}.
#' * if `k.vec` is not `NULL`, \eqn{L} is the length of `k.vec`.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(8, 9, NA, 3, 4, NA)
#' conf.int <- ci_quantile_g_mp_mn(
#'   Z, Y, k.vec = NULL, missing = "general", MWU = FALSE,
#'   method.list.treated = list(name = "Wilcoxon"),
#'   method.list.control = list(name = "Wilcoxon"),
#'   stat.null.treated = NULL, Z.perm.treated = NULL,
#'   stat.null.control = NULL, Z.perm.control = NULL,
#'   nperm = 10^4, alpha = 0.05, tol = 10^(-3)
#' )
ci_quantile_g_mp_mn <- function(Z, Y, k.vec = NULL, missing = "general", MWU = FALSE,
                                method.list.treated = list(name = "Wilcoxon"), method.list.control = list(name = "Wilcoxon"),
                                stat.null.treated = NULL, Z.perm.treated = NULL, stat.null.control = NULL, Z.perm.control = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units
  n0 <- n - n1 # number of control units

  # confidence interval for treated units
  conf.int.treated <- data.frame(k = c(1:n1), lower = rep(NA, n1), upper = rep(NA, n1))
  missing.treated <- missing
  c.lower.treated <- conf_quant_treated_g_mp_mn(
    Z = Z, Y = Y, k.vec = NULL, missing = missing.treated, MWU = MWU, method.list = method.list.treated,
    stat.null = stat.null.treated, Z.perm = Z.perm.treated, nperm = nperm, alpha = alpha / 2, tol = tol
  )
  conf.int.treated$upper <- Inf
  conf.int.treated$lower <- c.lower.treated

  # confidence interval for control units
  conf.int.control <- data.frame(k = c(1:n0), lower = rep(NA, n0), upper = rep(NA, n0))

  if (missing == "general") {
    missing.control <- "general"
  } else if (missing == "mp") {
    missing.control <- "mn"
  } else if (missing == "mn") {
    missing.control <- "mp"
  }

  c.lower.control <- conf_quant_treated_g_mp_mn(
    Z = 1 - Z, Y = -Y, k.vec = NULL, missing = missing.control, MWU = MWU, method.list = method.list.control,
    stat.null = stat.null.control, Z.perm = Z.perm.control, nperm = nperm, alpha = alpha / 2, tol = tol
  )

  conf.int.control$upper <- Inf
  conf.int.control$lower <- c.lower.control

  # Combine prediction intervals for treated and control units
  conf.int <- rbind(conf.int.treated, conf.int.control)
  conf.int$k <- NULL
  conf.int <- conf.int[order(conf.int$lower), ]
  conf.int$k <- c(1:n)

  if (!is.null(k.vec)) {
    conf.int <- conf.int[k.vec, ]
  }

  return(conf.int)
}

#--- sharp missing ---#

#' Calculate one-sided confidence interval (c, Inf) for all quantiles among observed treated units with missing outcomes for sharp/random missing
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k.vec A vector that specifies the quantiles of individual effects of interest.
#' If `k.vec = NULL`, then we consider all quantiles of individual effects \eqn{1 \leq k \leq n_{11}}.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistics.
#' @param method.list A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n.obs choose n1.obs).
#' * if `stat.null = NULL` (default), the function will calculate it internally.
#' @param Z.perm An \eqn{n.obs \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n.obs choose n1.obs).
#' * if `Z.perm = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A scalar that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return An \eqn{L \times 1} vector corresponding to the lower confidence limits of the given \eqn{\tau_{(k)}}.
#' * if `k.vec = NULL`, \eqn{L = n_{11}}.
#' * if `k.vec` is not `NULL`, \eqn{L} is the length of `k.vec`.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(8, 9, NA, 3, 4, NA)
#' c.limit <- conf_quant_treated_obs_s_r(
#'   Z, Y, k.vec = NULL,
#'   MWU = FALSE, method.list = list(name = "Wilcoxon"),
#'   stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)
#' )
conf_quant_treated_obs_s_r <- function(Z, Y, k.vec = NULL, MWU = FALSE, method.list = list(name = "Wilcoxon"),
                                       stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  # Drop missing data
  Y.obs <- Y[!is.na(Y)]
  Z.obs <- Z[!is.na(Y)]

  # Approach 1:
  # c.limit <- conf_quant_treated_g_mp_mn(Z = Z.obs, Y = Y.obs, k.vec = k.vec, missing = "general", MWU = MWU, method.list = method.list,
  #                                       stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha, tol = tol)
  # return (c.limit)

  # Approach 2:
  n.obs <- length(Z.obs) # n11 + n01
  n1.obs <- sum(Z.obs) # n11
  if (is.null(stat.null)) {
    stat.null <- null_dist(n.obs, n1.obs, MWU = MWU, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Find threshold such that p-value <= alpha <=> test statistic > threshold
  thres <- sort(stat.null, decreasing = TRUE)[floor(nperm * alpha) + 1]

  # Calculate range of c
  Y1.max <- max(Y.obs[Z.obs == 1])
  Y1.min <- min(Y.obs[Z.obs == 1])
  Y0.max <- max(Y.obs[Z.obs == 0])
  Y0.min <- min(Y.obs[Z.obs == 0])
  c.max <- Y1.max - Y0.min + tol
  c.min <- Y1.min - Y0.max - tol

  if (is.null(k.vec)) {
    # Initialize conf interval for all quantiles tau_{(k)}, 1<=k<=n1.obs
    c.limit <- rep(NA, n1.obs)

    for (k in n1.obs:1) {
      if (k < n1.obs) {
        c.max <- c.limit[k + 1]
      }
      # define the target fun
      # f > 0 <=> p-value <= alpha
      # f decreases in c, p value increases in c
      f <- function(c) {
        stat.min <- min_stat(Z.obs, Y.obs, k, c, MWU = MWU, method.list = method.list)
        return(stat.min - thres)
      }

      if (f(c.min) <= 0) {
        c.sol <- -Inf
      }
      if (f(c.max) > 0) {
        c.sol <- c.max
      }
      if (f(c.min) > 0 & f(c.max) <= 0) {
        c.sol <- stats::uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        # find the min c s.t. p-value > alpha <=> f <= 0
        c.sol <- round(c.sol, digits = -log10(tol))
        if (f(c.sol) <= 0) {
          while (f(c.sol) <= 0) {
            c.sol <- c.sol - tol
          }
          c.sol <- c.sol + tol
        } else {
          while (f(c.sol) > 0) {
            c.sol <- c.sol + tol
          }
        }
      }
      c.limit[k] <- c.sol
    }

    c.limit[c.limit > (Y1.max - Y0.min) + tol / 2] <- Inf
  } else {
    k.vec.sort <- sort(k.vec, decreasing = FALSE)
    j.max <- length(k.vec.sort)
    j.min <- 1

    # conf interval for tau_{(k)} with k in k.vec
    c.limit <- rep(NA, j.max)

    for (j in j.max:j.min) {
      k <- k.vec.sort[j]
      if (j < j.max) {
        c.max <- c.limit[j + 1]
      }
      # define the target fun
      # f > 0 <=> p-value <= alpha
      # f decreases in c, p value increases in c
      f <- function(c) {
        stat.min <- min_stat(Z.obs, Y.obs, k, c, MWU = MWU, method.list = method.list)
        return(stat.min - thres)
      }
      # check whether f(-Inf) = f(c.min) > 0
      if (f(c.min) <= 0) {
        c.sol <- -Inf
      }
      if (f(c.max) > 0) {
        c.sol <- c.max
      }
      if (f(c.min) > 0 & f(c.max) <= 0) {
        c.sol <- stats::uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        # find the min c st p-value > alpha <=> f <= 0
        c.sol <- round(c.sol, digits = -log10(tol))
        if (f(c.sol) <= 0) {
          while (f(c.sol) <= 0) {
            c.sol <- c.sol - tol
          }
          c.sol <- c.sol + tol
        } else {
          while (f(c.sol) > 0) {
            c.sol <- c.sol + tol
          }
        }
      }
      c.limit[j] <- c.sol
    }

    c.limit[c.limit > (Y1.max - Y0.min) + tol / 2] <- Inf
  }

  return(c.limit)
}

#' Calculate one-sided confidence interval (c, Inf) for all quantiles among observed units with missing outcomes for sharp/random missing
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k.vec A vector that specifies the quantiles of individual effects of interest.
#' If `k.vec = NULL`, then we consider all quantiles of individual effects \eqn{1 \leq k \leq n_{s}}.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistics.
#' @param method.list.treated A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param method.list.control A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null.treated An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n.obs choose n1.obs).
#' * if `stat.null.treated = NULL` (default), the function will calculate it internally.
#' @param Z.perm.treated An \eqn{n.obs \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n.obs choose n1.obs).
#' * if `Z.perm.treated = NULL` (default), the function will calculate it internally.
#' @param stat.null.control An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n.obs choose n0.obs).
#' * if `stat.null.control = NULL` (default), the function will calculate it internally.
#' @param Z.perm.control An \eqn{n.obs \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n.obs choose n0.obs).
#' * if `Z.perm.control = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A scalar that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return An \eqn{L \times 1} vector corresponding to the lower confidence limits of the given \eqn{\tau_{(k)}}.
#' * if `k.vec = NULL`, \eqn{L = n_{s}}.
#' * if `k.vec` is not `NULL`, \eqn{L} is the length of `k.vec`.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(8, 9, NA, 3, 4, NA)
#' c.limit <- conf_quant_obs_s_r(
#'   Z, Y, k.vec = NULL, MWU = FALSE,
#'   method.list.treated = list(name = "Wilcoxon"),
#'   method.list.control = list(name = "Wilcoxon"),
#'   stat.null.treated = NULL, Z.perm.treated = NULL,
#'   stat.null.control = NULL, Z.perm.control = NULL,
#'   nperm = 10^4, alpha = 0.05, tol = 10^(-3)
#' )
conf_quant_obs_s_r <- function(Z, Y, k.vec = NULL, MWU = FALSE, method.list.treated = list(name = "Wilcoxon"), method.list.control = list(name = "Wilcoxon"),
                               stat.null.treated = NULL, Z.perm.treated = NULL, stat.null.control = NULL, Z.perm.control = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  n_s <- length(Z[!is.na(Y)]) # n_s = n11 + n01
  n11 <- sum(Z[!is.na(Y)]) # n11
  n01 <- n_s - n11 # n01

  # confidence interval for observed treated units
  conf.int.treated.obs <- data.frame(k = c(1:n11), lower = rep(NA, n11), upper = rep(NA, n11))
  c.lower.treated.obs <- conf_quant_treated_obs_s_r(
    Z = Z, Y = Y, k.vec = NULL, MWU = MWU, method.list = method.list.treated,
    stat.null = stat.null.treated, Z.perm = Z.perm.treated, nperm = nperm, alpha = alpha / 2, tol = tol
  )
  conf.int.treated.obs$upper <- Inf
  conf.int.treated.obs$lower <- c.lower.treated.obs

  # confidence interval for observed control units
  conf.int.control.obs <- data.frame(k = c(1:n01), lower = rep(NA, n01), upper = rep(NA, n01))
  c.lower.control.obs <- conf_quant_treated_obs_s_r(
    Z = 1 - Z, Y = -Y, k.vec = NULL, MWU = MWU, method.list = method.list.control,
    stat.null = stat.null.control, Z.perm = Z.perm.control, nperm = nperm, alpha = alpha / 2, tol = tol
  )
  conf.int.control.obs$upper <- Inf
  conf.int.control.obs$lower <- c.lower.control.obs

  # Combine prediction intervals for observed treated and control units
  conf.int.obs <- rbind(conf.int.treated.obs, conf.int.control.obs)
  conf.int.obs$k <- NULL
  conf.int.obs <- conf.int.obs[order(conf.int.obs$lower), ]
  conf.int.obs$k <- c(1:n_s)

  if (!is.null(k.vec)) {
    conf.int.obs <- conf.int.obs[k.vec, ]
  }

  return(conf.int.obs)
}

#' Randomization inference for quantiles of individual treatment effects under sharp missing mechanism
#'
#' @description `ci_quantile_s()` obtains one-sided confidence intervals
#' for all quantiles of individual treatment effects.
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k.vec A vector that specifies the quantiles of individual effects of interest.
#' If `k.vec = NULL`, then we consider all quantiles of individual effects \eqn{1 \leq k \leq n}.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistics.
#' @param method.list.treated A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param method.list.control A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null.treated An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n.obs choose n1.obs).
#' * if `stat.null.treated = NULL` (default), the function will calculate it internally.
#' @param Z.perm.treated An \eqn{n.obs \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n.obs choose n1.obs).
#' * if `Z.perm.treated = NULL` (default), the function will calculate it internally.
#' @param stat.null.control An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n.obs choose n0.obs).
#' * if `stat.null.control = NULL` (default), the function will calculate it internally.
#' @param Z.perm.control An \eqn{n.obs \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n.obs choose n0.obs).
#' * if `Z.perm.control = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A scalar that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return A dataframe with \eqn{L} rows and 3 columns (\eqn{k, lower, upper}).
#' Each row corresponds to the lower and upper confidence limits of the given \eqn{\tau_{(k)}}.
#' * if `k.vec = NULL`, \eqn{L = n}.
#' * if `k.vec` is not `NULL`, \eqn{L} is the length of `k.vec`.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(8, 9, NA, 3, 4, NA)
#' conf.int <- ci_quantile_s(
#'   Z, Y, k.vec = NULL, MWU = FALSE,
#'   method.list.treated = list(name = "Wilcoxon"),
#'   method.list.control = list(name = "Wilcoxon"),
#'   stat.null.treated = NULL, Z.perm.treated = NULL,
#'   stat.null.control = NULL, Z.perm.control = NULL,
#'   nperm = 10^4, alpha = 0.05, tol = 10^(-3)
#' )
ci_quantile_s <- function(Z, Y, k.vec = NULL, MWU = FALSE,
                          method.list.treated = list(name = "Wilcoxon"), method.list.control = list(name = "Wilcoxon"),
                          stat.null.treated = NULL, Z.perm.treated = NULL, stat.null.control = NULL, Z.perm.control = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  conf.int.obs <- conf_quant_obs_s_r(Z, Y,
    k.vec = NULL, MWU = MWU, method.list.treated = method.list.treated, method.list.control = method.list.control,
    stat.null.treated = stat.null.treated, Z.perm.treated = Z.perm.treated, stat.null.control = stat.null.control, Z.perm.control = Z.perm.control, nperm = nperm, alpha = alpha, tol = tol
  )

  n <- length(Z)
  n_s <- length(Z[!is.na(Y)]) # n_s = n11 + n01
  n11 <- sum(Z[!is.na(Y)]) # n11
  n01 <- n_s - n11 # n01

  # confidence interval for all units
  conf.int <- data.frame(k = c(1:n), lower = rep(NA, n), upper = rep(NA, n))
  conf.int$upper <- Inf
  conf.int$lower <- -Inf
  conf.int$lower[(n - n_s + 1):n] <- conf.int.obs$lower

  if (!is.null(k.vec)) {
    conf.int <- conf.int[k.vec, ]
  }

  return(conf.int)
}

#--- missing at random ---#

correct_hypergeom <- function(N = 100, n = 50, K.vec = c(60, 80), kk.vec = c(30, 40), ndraw = 10^4, hg.draw = NULL) {
  if (length(K.vec) > 1) {
    if (is.null(hg.draw)) {
      hg.draw <- extraDistr::rmvhyper(ndraw, diff(c(0, K.vec, N)), n)
      hg.draw <- hg.draw[, -1, drop = FALSE]
    }
    H.sum <- apply(hg.draw[, ncol(hg.draw):1, drop = FALSE], 1, cumsum)
    H.sum <- matrix(H.sum, nrow = ncol(hg.draw))
    H.sum <- H.sum[nrow(H.sum):1, , drop = FALSE]
    prob <- mean(apply(H.sum - (n - kk.vec) > 0, 2, max))
    return(prob)
  }

  if (length(K.vec) == 1) {
    prob <- 1 - stats::phyper(n - kk.vec, N - K.vec, K.vec, n)
    return(prob)
  }
}

#' Calculate the correction term Delta using the proposed choice of kk.vec (k_prime)
threshold_correct_hypergeom <- function(N = 100, n = 50, K.vec = c(60, 80), alpha = 0.05, ndraw = 10^4, hg.draw = NULL, tol = 10^(-2)) {
  if (length(K.vec) > 1) {
    if (is.null(hg.draw)) {
      hg.draw <- extraDistr::rmvhyper(ndraw, diff(c(0, K.vec, N)), n)
      hg.draw <- hg.draw[, -1, drop = FALSE]
    }

    kappa <- NULL

    kappa.lower <- 1 / length(K.vec)
    kappa.upper <- 1
    prob.lower <- correct_hypergeom(N = N, n = n, K.vec = K.vec, kk.vec = n - stats::qhyper(1 - kappa.lower * alpha, N - K.vec, K.vec, n), hg.draw = hg.draw)
    prob.upper <- correct_hypergeom(N = N, n = n, K.vec = K.vec, kk.vec = n - stats::qhyper(1 - kappa.upper * alpha, N - K.vec, K.vec, n), hg.draw = hg.draw)

    if (prob.upper <= alpha) {
      kappa <- kappa.upper
    }

    if (prob.lower > alpha) {
      kappa <- kappa.lower
    }

    while (is.null(kappa)) {
      kappa.mid <- (kappa.upper + kappa.lower) / 2
      prob.mid <- correct_hypergeom(N = N, n = n, K.vec = K.vec, kk.vec = n - stats::qhyper(1 - kappa.mid * alpha, N - K.vec, K.vec, n), hg.draw = hg.draw)

      if (prob.mid <= alpha) {
        kappa.lower <- kappa.mid
      } else {
        kappa.upper <- kappa.mid
      }

      if (kappa.upper - kappa.lower <= tol) {
        kappa <- kappa.lower
      }
    }

    kk.vec <- n - stats::qhyper(1 - kappa * alpha, N - K.vec, K.vec, n)
    prob <- correct_hypergeom(N = N, n = n, K.vec = K.vec, kk.vec = kk.vec, hg.draw = hg.draw)
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }

  if (length(K.vec) == 1) {
    kappa <- 1
    temp <- stats::qhyper(1 - alpha, N - K.vec, K.vec, n)
    kk.vec <- n - stats::qhyper(1 - alpha, N - K.vec, K.vec, n)
    prob <- correct_hypergeom(N = N, n = n, K.vec = K.vec, kk.vec = kk.vec)
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }
}

#' Randomization inference for quantiles of individual treatment effects under random missing mechanism
#'
#' @description `ci_quantile_r()` obtains one-sided confidence intervals
#' for all quantiles of individual treatment effects.
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k.vec A vector that specifies the quantiles of individual effects of interest.
#' If `k.vec = NULL`, then we consider all quantiles of individual effects \eqn{1 \leq k \leq n}.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistic.
#' @param method.list.treated A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param method.list.control A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null.treated An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n.obs choose n1.obs).
#' * if `stat.null.treated = NULL` (default), the function will calculate it internally.
#' @param Z.perm.treated An \eqn{n.obs \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n.obs choose n1.obs).
#' * if `Z.perm.treated = NULL` (default), the function will calculate it internally.
#' @param stat.null.control An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n.obs choose n0.obs).
#' * if `stat.null.control = NULL` (default), the function will calculate it internally.
#' @param Z.perm.control An \eqn{n.obs \times nperm} matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic (n.obs choose n0.obs).
#' * if `Z.perm.control = NULL` (default), the function will calculate it internally.
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A scalar that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return A dataframe with \eqn{L} rows and 3 columns (\eqn{k, lower, upper}).
#' Each row corresponds to the lower and upper confidence limits of the given \eqn{\tau_{(k)}}.
#' * if `k.vec = NULL`, \eqn{L = n}.
#' * if `k.vec` is not `NULL`, \eqn{L} is the length of `k.vec`.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(8, 9, NA, 3, 4, NA)
#' conf.int <- ci_quantile_r(
#'   Z, Y, k.vec = NULL, MWU = FALSE,
#'   method.list.treated = list(name = "Wilcoxon"),
#'   method.list.control = list(name = "Wilcoxon"),
#'   stat.null.treated = NULL, Z.perm.treated = NULL,
#'   stat.null.control = NULL, Z.perm.control = NULL,
#'   nperm = 10^4, alpha = 0.05, tol = 10^(-3)
#' )
ci_quantile_r <- function(Z, Y, k.vec = NULL, MWU = FALSE,
                          method.list.treated = list(name = "Wilcoxon"), method.list.control = list(name = "Wilcoxon"),
                          stat.null.treated = NULL, Z.perm.treated = NULL, stat.null.control = NULL, Z.perm.control = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  n <- length(Z)
  n_s <- length(Z[!is.na(Y)]) # n_s = n11 + n01

  # confidence interval for all units
  conf.int <- data.frame(k = c(1:n), lower = rep(NA, n), upper = rep(NA, n))
  conf.int$upper <- Inf
  conf.int$lower <- -Inf

  conf.int.obs <- conf_quant_obs_s_r(Z, Y,
    k.vec = NULL, MWU = MWU, method.list.treated = method.list.treated, method.list.control = method.list.control,
    stat.null.treated = stat.null.treated, Z.perm.treated = Z.perm.treated, stat.null.control = stat.null.control, Z.perm.control = Z.perm.control, nperm = nperm, alpha = alpha / 2, tol = tol
  )

  for (k in 1:n) {
    step1 <- threshold_correct_hypergeom(N = n, n = n_s, K.vec = c(k), alpha = alpha / 2, ndraw = 10^4, hg.draw = NULL, tol = 10^(-2))

    # conf.int.obs <- conf_quant_obs_s_r(Z, Y, k.vec = NULL, MWU = MWU, method.list.treated = method.list.treated, method.list.control = method.list.control,
    #                                    stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha - step1$prob, tol = tol)
    if (step1$kk.vec == 0) {
      conf.int$lower[conf.int$k == k] <- -Inf
      conf.int$upper[conf.int$k == k] <- Inf
    } else {
      conf.int$lower[conf.int$k == k] <- conf.int.obs$lower[conf.int.obs$k == step1$kk.vec]
      conf.int$upper[conf.int$k == k] <- conf.int.obs$upper[conf.int.obs$k == step1$kk.vec]
    }
  }

  if (!is.null(k.vec)) {
    conf.int <- conf.int[k.vec, ]
  }

  return(conf.int)
}

#' Randomization inference for quantiles of individual treatment effects with missing outcomes
#'
#' @description `ci_quantile()` obtains one-sided confidence intervals
#' for all quantiles of individual treatment effects.
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param k.vec A vector that specifies the quantiles of individual effects of interest.
#' If `k.vec = NULL`, then we consider all quantiles of individual effects \eqn{1 \leq k \leq n}.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' * `missing = "sharp"`: sharp missing mechanism (\eqn{M_1 = M_0}).
#' * `missing = "random"`: random missing mechanism.
#' @param MWU A logical value that specifies the class of test statistic.
#' * `MWU = FALSE` (first class statistic, default): rank-sum statistic.
#' * `MWU = TRUE` (second class statistic): generalized Mann-Whitney U statistics.
#' @param method.list.treated A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.treated` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param method.list.control A list that specifies the choice of the rank-based test statistic.
#' * if `MWU = FALSE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Stephenson", s = 10)`.
#' * if `MWU = TRUE`, `method.list.control` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
#' @param stat.null.treated An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic.
#' * if `stat.null.treated = NULL` (default), the function will calculate it internally.
#' * if `missing = "general"`, `mp`, `mn`, then n choose n1.
#' * if `missing = "sharp"`, `random`, then n.obs choose n1.obs.
#' @param Z.perm.treated A matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic.
#' * if `Z.perm.treated = NULL` (default), the function will calculate it internally.
#' * if `missing = "general"`, `mp`, `mn`, then `Z.perm.treated` is an \eqn{n \times nperm} matrix (n choose n1).
#' * if `missing = "sharp"`, `random`, then `Z.perm.treated` is an \eqn{n.obs \times nperm} matrix (n.obs choose n1.obs).
#' @param stat.null.control An \eqn{nperm \times 1} vector whose empirical distribution
#' approximates the randomization distribution of the rank sum statistic (n.obs choose n0.obs).
#' * if `stat.null.control = NULL` (default), the function will calculate it internally.
#' * if `missing = "general"`, `mp`, `mn`, then n choose n0.
#' * if `missing = "sharp"`, `random`, then n.obs choose n0.obs.
#' @param Z.perm.control A matrix that specifies the permutated assignments
#' for approximating the null distribution of the test statistic.
#' * if `Z.perm.control = NULL` (default), the function will calculate it internally.
#' * if `missing = "general"`, `mp`, `mn`, then `Z.perm.control` is an \eqn{n \times nperm} matrix (n choose n0).
#' * if `missing = "sharp"`, `random`, then `Z.perm.control` is an \eqn{n.obs \times nperm} matrix (n.obs choose n0.obs).
#' @param nperm A positive integer representing the number of permutations to
#' approximate the randomization distribution of the test statistic.
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A scalar that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return A dataframe with \eqn{L} rows and 3 columns (\eqn{k, lower, upper}).
#' Each row corresponds to the lower and upper confidence limits of the given \eqn{\tau_{(k)}}.
#' * if `k.vec = NULL`, \eqn{L = n}.
#' * if `k.vec` is not `NULL`, \eqn{L} is the length of `k.vec`.
#' @export
#'
#' @examples
#' Z <- c(1, 1, 1, 0, 0, 0)
#' Y <- c(8, 9, NA, 3, 4, NA)
#' conf.int <- ci_quantile(
#'   Z, Y, k.vec = NULL, missing = "general", MWU = FALSE,
#'   method.list.treated = list(name = "Wilcoxon"),
#'   method.list.control = list(name = "Wilcoxon"),
#'   stat.null.treated = NULL, Z.perm.treated = NULL,
#'   stat.null.control = NULL, Z.perm.control = NULL,
#'   nperm = 10^4, alpha = 0.05, tol = 10^(-3)
#' )
ci_quantile <- function(Z, Y, k.vec = NULL, missing = "general", MWU = FALSE,
                        method.list.treated = list(name = "Wilcoxon"), method.list.control = list(name = "Wilcoxon"),
                        stat.null.treated = NULL, Z.perm.treated = NULL, stat.null.control = NULL, Z.perm.control = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3)) {
  if (missing %in% c("general", "mp", "mn")) {
    conf.int <- ci_quantile_g_mp_mn(Z, Y,
      k.vec = k.vec, missing = missing, MWU = MWU,
      method.list.treated = method.list.treated, method.list.control = method.list.control,
      stat.null.treated = stat.null.treated, Z.perm.treated = Z.perm.treated, stat.null.control = stat.null.control, Z.perm.control = Z.perm.control, nperm = nperm, alpha = alpha, tol = tol
    )
  }

  if (missing == "sharp") {
    conf.int <- ci_quantile_s(Z, Y,
      k.vec = k.vec, MWU = MWU,
      method.list.treated = method.list.treated, method.list.control = method.list.control,
      stat.null.treated = stat.null.treated, Z.perm.treated = Z.perm.treated, stat.null.control = stat.null.control, Z.perm.control = Z.perm.control, nperm = nperm, alpha = alpha, tol = tol
    )
  }

  if (missing == "random") {
    conf.int <- ci_quantile_r(Z, Y,
      k.vec = k.vec, MWU = MWU,
      method.list.treated = method.list.treated, method.list.control = method.list.control,
      stat.null.treated = stat.null.treated, Z.perm.treated = Z.perm.treated, stat.null.control = stat.null.control, Z.perm.control = Z.perm.control, nperm = nperm, alpha = alpha, tol = tol
    )
  }

  return(conf.int)
}
