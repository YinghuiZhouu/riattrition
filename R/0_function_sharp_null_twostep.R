#--- Functions for Testing Fisher Sharp Null Hypothesis with Two-Step Method --#

#------------- Helper Functions for pval_sharp_twostep() ----------------------#

#' Construct confidence intervals of hypergeometric parameter
#'
#' `lci()` generates the lower limit of \eqn{1 - \alpha} CP interval.
#' Consider hypergeometric with \eqn{N} units, \eqn{M} of them are good (\eqn{N - M} are bad),
#' select \eqn{n} from \eqn{N} units, there are \eqn{x} good units in the sample.
#' P(X = x) = dhyper(x, m, N - m, n) for max(0, n + M - N) <= x <= min(n, M)
#' Here, \eqn{n} and \eqn{N} are known, and we want to estimate \eqn{M} based on \eqn{x}.
#' Source: "Exact Optimal Confidence Intervals for Hypergeometric Parameters"
#' https://doi.org/10.1080/01621459.2014.966191
#'
#' @param N Number of total units.
#' @param x Number of observed good units.
#' @param n Number of draws.
#' @param alpha Confidence parameter.
#'
#' @return The lower limit of \eqn{1 - \alpha} CP interval.
#' @examples
#' lci(N = 100, x = 50, n = 50, alpha = 0.1)
#' @export
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

#' Construct confidence intervals of hypergeometric parameter
#'
#' `uci()` generates the upper limit of \eqn{1 - \alpha} CP interval.
#' Consider hypergeometric with \eqn{N} units, \eqn{M} of them are good (\eqn{N - M} are bad),
#' select \eqn{n} from \eqn{N} units, there are \eqn{x} good units in the sample.
#' P(X = x) = dhyper(x, m, N - m, n) for max(0, n + M - N) <= x <= min(n, M)
#' Here, \eqn{n} and \eqn{N} are known, and we want to estimate \eqn{M} based on \eqn{x}.
#' Source: "Exact Optimal Confidence Intervals for Hypergeometric Parameters"
#' https://doi.org/10.1080/01621459.2014.966191
#'
#' @param N Number of total units.
#' @param x Number of observed good units.
#' @param n Number of draws.
#' @param alpha Confidence parameter.
#'
#' @return The upper limit of \eqn{1 - \alpha} CP interval.
#' @examples
#' uci(N = 100, x = 50, n = 50, alpha = 0.1)
#' @export
uci <- function(N, x, n, alpha) {
  return(N - lci(N, n - x, n, alpha))
}

#' Construct confidence intervals of hypergeometric parameter
#'
#' `exactci()` generates the two-sided of \eqn{1 - \alpha} CP interval.
#' Consider hypergeometric with \eqn{N} units, \eqn{M} of them are good (\eqn{N - M} are bad),
#' select \eqn{n} from \eqn{N} units, there are \eqn{x} good units in the sample.
#' P(X = x) = dhyper(x, m, N - m, n) for max(0, n + M - N) <= x <= min(n, M)
#' Here, \eqn{n} and \eqn{N} are known, and we want to estimate \eqn{M} based on \eqn{x}.
#' Source: "Exact Optimal Confidence Intervals for Hypergeometric Parameters"
#' https://doi.org/10.1080/01621459.2014.966191
#'
#' @param N Number of total units.
#' @param x Number of observed good units.
#' @param n Number of draws.
#' @param alpha Confidence parameter.
#'
#' @return A list of lower and upper limit of \eqn{1 - \alpha} CP interval.
#' @examples
#' exactci(N = 100, x = 50, n = 50, alpha = 0.1)
#' @export
exactci <- function(N, x, n, alpha) {
  # indicator fucntion of interval [a, b]
  ind <- function(x, a, b) {
    (x >= a) * (x <= b)
  }

  xx <- 0:n
  lcin1 <- lci(N, xx, n, alpha / 2)
  ucin1 <- uci(N, xx, n, alpha / 2) # two-sided 1-alpha interval
  lcin2 <- lci(N, xx, n, alpha)
  ucin2 <- uci(N, xx, n, alpha) # two-sided 1-2alpha interval

  lciw <- lcin1 #  lcin1<=lciw<=lcin2
  uciw <- ucin1 #  ucin2<=uciw<=ucin1

  # for an even number n
  xvalue <- n / 2 + 1 # start from the center
  aa <- lciw[xvalue]:floor(N / 2)

  ii <- 1
  while (ii < length(aa) + 0.5) {
    lciw[xvalue] <- aa[ii]
    uciw[xvalue] <- N - aa[ii]

    # the coverage probability function for the combined interval.
    cpci <- function(M) {
      kk <- 1:length(M)
      for (i in kk) {
        xx <- 0:n
        indp <- xx
        uu <- 0
        while (uu < n + 0.5) {
          indp[uu + 1] <- ind(M[i], lciw[uu + 1], uciw[uu + 1]) * stats::dhyper(uu, M[i], N - M[i], n)
          uu <- uu + 1
        }
        kk[i] <- sum(indp)
      }
      return(kk)
    }
    M <- 0:N
    bb <- min(cpci(M))
    if (bb >= 1 - alpha) {
      ii1 <- ii
      ii <- ii + 1
    } else {
      ii <- length(aa) + 1
    }
  }
  lciw[xvalue] <- aa[ii1]
  uciw[xvalue] <- N - lciw[xvalue]


  xvalue <- n / 2 # xvalue >=1
  while (xvalue > 0.5) {
    al <- lcin2[xvalue] - lciw[xvalue] + 1
    au <- uciw[xvalue] - ucin2[xvalue] + 1

    if (al * au > 1) {
      ff <- array(, dim = c(al * au, 4))
      for (i in 1:al) {
        ff[((i - 1) * au + 1):(i * au), 1] <- lciw[xvalue] + i - 1
        ff[((i - 1) * au + 1):(i * au), 2] <- (ucin2[xvalue]):(uciw[xvalue])
        ff[((i - 1) * au + 1):(i * au), 3] <- ff[((i - 1) * au + 1):(i * au), 2] - ff[((i - 1) * au + 1):(i * au), 1]
      }

      for (ii in 1:dim(ff)[1]) {
        lciw[xvalue] <- ff[ii, 1]
        uciw[xvalue] <- ff[ii, 2]
        lciw[n + 2 - xvalue] <- N - uciw[xvalue]
        uciw[n + 2 - xvalue] <- N - lciw[xvalue]

        # the coverage probability function for the combined interval.
        cpci <- function(M) {
          kk <- 1:length(M)
          for (i in kk) {
            xx <- 0:n
            indp <- xx
            uu <- 0
            while (uu < n + 0.5) {
              indp[uu + 1] <- ind(M[i], lciw[uu + 1], uciw[uu + 1]) * stats::dhyper(uu, M[i], N - M[i], n)
              uu <- uu + 1
            }
            kk[i] <- sum(indp)
          }
          return(kk)
        }
        M <- 0:N
        ff[ii, 4] <- min(cpci(M))
      }

      ff <- ff[which(ff[, 4] >= 1 - alpha), ]
      if (length(ff) > 4) {
        ff <- ff[order(ff[, 3]), ]
        lciw[xvalue] <- ff[1, 1]
        uciw[xvalue] <- ff[1, 2]
      } else {
        lciw[xvalue] <- ff[1]
        uciw[xvalue] <- ff[2]
      }
      lciw[n + 2 - xvalue] <- N - uciw[xvalue]
      uciw[n + 2 - xvalue] <- N - lciw[xvalue]
    }
    xvalue <- xvalue - 1
  }
  cbind(xx, lciw, uciw, lcin1, ucin1) # the improved and original intervals
  lower <- lciw[xx == x]
  upper <- uciw[xx == x]
  return(list(lower = lower, upper = upper))
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
#' `pval_sharp_g_twostep()` obtains the p-value for testing the sharp null hypothesis
#' \eqn{H_0: \tau = c}, using the two-step procedure under general missing mechanism.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param method.list A list that specifies the choice of test statistic:
#' `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
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
#' @param beta A real number that belongs to \eqn{[0, \alpha]}.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
pval_sharp_g_twostep <- function(Z, Y, c = 0, method.list = list(name = "Wilcoxon"),
                                 stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z) # number of observations
  n.obs <- length(Z[M == 1]) # number of observed units (n11 + n01)
  n1 <- sum(Z) # number of treated units
  n0 <- n - n1 # number of control units
  n11 <- sum(Z[M == 1]) # number of observed treated units
  n10 <- sum(Z[M == 0]) # number of missing treated units
  n01 <- n.obs - n11 # number of observed control units
  n00 <- n0 - n01 # number of missing control units

  # Step 1: construct confidence intervals of hypergeometric parameter
  beta1 <- beta / 2
  beta2 <- beta / 2

  if (beta > 0) {
    M_mat <- ExactCIone::WhyperCI_M(x = n01, n = n0, N = n, conf.level = 1 - beta1)
    M1hat <- M_mat$CI[1, 2]
    M2hat <- M_mat$CI[1, 3]
    # Approach 2 (run slowly)
    # Mhat_list <- exactci(N = n, x = n01, n = n0, alpha = beta1)
    # M1hat <- Mhat_list$lower
    # M2hat <- Mhat_list$upper

    # Approach 3
    # M1hat <- lci(N = n, x = n01, n = n0, alpha = beta1 / 2)
    # M2hat <- uci(N = n, x = n01, n = n0, alpha = beta1 / 2)
  }
  if (beta == 0) {
    M1hat <- n01
    M2hat <- n1 + n01
  }
  mlbar <- M1hat - n01
  mubar <- M2hat - n01

  # Calculate dubar by looping over all candidate m11 values in [0, n.obs]
  m11_grid <- 0:n.obs
  q_hg_vec <- stats::qhyper(p = 1 - beta2, m = m11_grid, n = n - m11_grid, k = n0)
  dubar <- max(n / (n1 * n0) * q_hg_vec - m11_grid / n1)

  # Step 2
  ind.treat <- which(Z == 1) # n1 x 1 vector
  ind.control <- which(Z == 0) # n0 x 1 vector
  ind.control.obs <- which(Z == 0 & M == 1) # n01 x 1 vector
  ind.control.miss <- which(Z == 0 & M == 0) # n00 x 1 vector
  ind.treat.obs <- which(Z == 1 & M == 1) # n11 x 1 vector
  ind.treat.miss <- which(Z == 1 & M == 0) # n10 x 1 vector

  r <- rank(Y, ties.method = "first") # n x 1 vector
  ind.sort <- sort.int(r, index.return = TRUE)$ix # n x 1 vector
  ind.sort.control.obs <- ind.sort[Z[ind.sort] == 0 & M[ind.sort] == 1] # n01 x 1 vector

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

  # Calculate A_i for treated units
  A <- rep(NA, n) # n x 1 vector (where n1 elements are non NAs)
  A[ind.treat] <- n01 + findInterval(ind.treat, ind.control.miss)

  # Precompute B_J for missing treated units (n10 x 1 vector)
  # sum_{j:Z_j=0,M_j=1} 1(i>j)
  B_J_precompute_10 <- findInterval(ind.treat.miss, ind.control.obs)

  # Precompute B_J for observed treated units (n11 x (n01 + 1) matrix)
  # sum_{j in F_J} psi(i,j,Y_i-c,Y_j) + (n01 - J)
  if (n11 > 0 & n01 > 0) {
    # psi_mat stores all pairwise psi-values between treated observed units and control observed units
    psi_mat <- matrix(0, nrow = n11, ncol = n01) # n11 x n01 matrix

    for (a in 1:n11) {
      i <- ind.treat.obs[a]
      for (b in 1:n01) {
        j <- ind.sort.control.obs[b]
        psi_mat[a, b] <- psi(i, j, Y[i] - c, Y[j])
      }
    }

    # B_J_precompute_11[a, J+1] = sum of psi_mat[a, last J columns]
    B_J_precompute_11 <- matrix(0, nrow = n11, ncol = n01 + 1)
    psi_rev <- psi_mat[, n01:1, drop = FALSE]
    B_J_precompute_11[, 2:(n01 + 1)] <- t(apply(psi_rev, 1, cumsum))
  } else {
    B_J_precompute_11 <- matrix(0, nrow = n11, ncol = n01 + 1)
  }

  Klbar <- max(0, mlbar - n10)
  Kubar <- min(n11, mubar)

  TK_vec <- rep(NA, Kubar - Klbar + 1)
  for (K in Klbar:Kubar) {
    J <- min(floor(n0 * (dubar + K / n1)), n01)
    L <- min(mubar - K, n10)

    # Calculate B_Ji for treated units
    B_J <- rep(NA, n) # n x 1 vector (where n1 elements are non NAs)
    if (n11 > 0) {
      B_J[ind.treat.obs] <- B_J_precompute_11[, J + 1] + (n01 - J)
    }
    if (n10 > 0) {
      B_J[ind.treat.miss] <- pmax(0, B_J_precompute_10 - J)
    }
    C_J <- phi(A) - phi(B_J) # n x 1 vector (where n1 elements are non NAs)
    r_C <- rank(C_J, ties.method = "first") # n x 1 vector
    ind.sort.C <- sort.int(r_C, index.return = TRUE)$ix # n x 1 vector
    ind.sort.C.treat.obs <- ind.sort.C[Z[ind.sort.C] == 1 & M[ind.sort.C] == 1] # n11 x 1 vector
    ind.sort.C.treat.miss <- ind.sort.C[Z[ind.sort.C] == 1 & M[ind.sort.C] == 0] # n10 x 1 vector
    TK_vec[K - Klbar + 1] <- sum(phi(B_J)[Z == 1]) + sum(C_J[ind.sort.C.treat.obs[seq_len(n11 - K)]]) + sum(C_J[ind.sort.C.treat.miss[seq_len(n10 - L)]])
  }
  TK_min <- min(TK_vec)

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, class = "MWU+", method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- TK_min

  pval <- mean(stat.null >= stat.obs) + beta
  return(min(pval, 1))
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' `pval_sharp_mp_twostep()` obtains the p-value for testing the sharp null hypothesis
#' \eqn{H_0: \tau = c}, using the two-step procedure under monotone positive missing mechanism.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param method.list A list that specifies the choice of test statistic:
#' `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
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
#' @param beta A real number that belongs to \eqn{[0, \alpha]}.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
pval_sharp_mp_twostep <- function(Z, Y, c = 0, method.list = list(name = "Wilcoxon"),
                                  stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z) # number of observations
  n.obs <- length(Z[M == 1]) # number of observed units (n11 + n01)
  n1 <- sum(Z) # number of treated units
  n0 <- n - n1 # number of control units
  n11 <- sum(Z[M == 1]) # number of observed treated units
  n01 <- n.obs - n11 # number of observed control units

  # Step 1: construct confidence intervals of hypergeometric parameter
  Mhat <- uci(N = n, x = n01, n = n0, alpha = beta)
  mlbar <- max(n11 + n01 - Mhat, 0)

  # Step 2
  Y.comp <- ifelse(M == 1, Y, Inf) # Y_j * M_j + Inf * (1 - M_j)
  ind.treat.obs <- which(Z == 1 & M == 1)
  ind.control <- which(Z == 0)

  A <- rep(NA, n)
  B <- rep(NA, n)
  for (i in ind.treat.obs) {
    A[i] <- sum(psi(i = i, j = ind.control, x = Inf, y = Y.comp[ind.control]))
    B[i] <- sum(psi(i = i, j = ind.control, x = Y[i] - c, y = Y.comp[ind.control]))
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

  C <- phi(A) - phi(B) # n x 1 vector (where n11 elements are non NAs)
  r <- rank(C, ties.method = "first") # n x 1 vector
  ind.sort <- sort.int(r, index.return = TRUE)$ix # n x 1 vector
  ind.sort.treat.obs <- ind.sort[Z[ind.sort] == 1 & M[ind.sort] == 1] # n11 x 1 vector

  # Impute the worst-case composite control potential outcome
  Y0.com.imp <- rep(NA, n)
  Y0.com.imp[ind.sort.treat.obs[seq_len(mlbar)]] <- Inf
  Y0.com.imp[ind.sort.treat.obs[(mlbar + 1):n11]] <- Y[ind.sort.treat.obs[(mlbar + 1):n11]] - c
  Y0.com.imp[Z == 1 & M == 0] <- Inf
  Y0.com.imp[Z == 0 & M == 1] <- Y[Z == 0 & M == 1]
  Y0.com.imp[Z == 0 & M == 0] <- Inf

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, class = "MWU+", method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z, Y = Y0.com.imp, class = "MWU+", method.list = method.list)

  pval <- mean(stat.null >= stat.obs) + beta
  return(min(pval, 1))
}

#' Randomization test for sharp null hypotheses with missing outcomes
#'
#' `pval_sharp_mn_twostep()` obtains the p-value for testing the sharp null hypothesis
#' \eqn{H_0: \tau = c}, using the two-step procedure under monotone negative missing mechanism.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param method.list A list that specifies the choice of test statistic:
#' `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
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
#' @param beta A real number that belongs to \eqn{[0, \alpha]}.
#'
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
pval_sharp_mn_twostep <- function(Z, Y, c = 0, method.list = list(name = "Wilcoxon"),
                                  stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z) # number of observations
  n.obs <- length(Z[M == 1]) # number of observed units (n11 + n01)
  n1 <- sum(Z) # number of treated units
  n11 <- sum(Z[M == 1]) # number of observed treated units
  n01 <- n.obs - n11 # number of observed control units

  # Step 1: construct confidence intervals of hypergeometric parameter
  Mhat <- uci(N = n, x = n11, n = n1, alpha = beta)
  mlbar <- max(n11 + n01 - Mhat, 0)

  # Step 2
  Y.comp <- ifelse(M == 1, Y - c, -Inf) # (Y_j-c) * M_j + (-Inf) * (1 - M_j)
  ind.control.obs <- which(Z == 0 & M == 1)
  ind.treat <- which(Z == 1)

  A <- rep(NA, n)
  B <- rep(NA, n)
  for (i in ind.control.obs) {
    A[i] <- sum(psi(i = i, j = ind.treat, x = Y[i], y = Y.comp[ind.treat]))
    B[i] <- sum(psi(i = i, j = ind.treat, x = -Inf, y = Y.comp[ind.treat]))
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

  C <- phi(A) - phi(B) # n x 1 vector (where n01 elements are non NAs)
  r <- rank(C, ties.method = "first") # n x 1 vector
  ind.sort <- sort.int(r, index.return = TRUE)$ix # n x 1 vector
  ind.sort.control.obs <- ind.sort[Z[ind.sort] == 0 & M[ind.sort] == 1] # n01 x 1 vector

  # Impute the worst-case composite control potential outcome
  Y0.com.imp <- rep(NA, n)
  Y0.com.imp[Z == 1 & M == 1] <- Y[Z == 1 & M == 1] - c
  Y0.com.imp[Z == 1 & M == 0] <- -Inf
  Y0.com.imp[ind.sort.control.obs[(mlbar + 1):n01]] <- Y[ind.sort.control.obs[(mlbar + 1):n01]]
  Y0.com.imp[ind.sort.control.obs[seq_len(mlbar)]] <- -Inf
  Y0.com.imp[Z == 0 & M == 0] <- -Inf

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, class = "MWU-", method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  # Calculate the observed test statistic
  stat.obs <- test_stat(Z = Z, Y = Y0.com.imp, class = "MWU-", method.list = method.list)

  pval <- mean(stat.null >= stat.obs) + beta
  return(min(pval, 1))
}

#----------------- Helper Functions for ci_sharp_twostep() --------------------#

#' Confidence interval assuming constant treatment effect
#'
#' `ci_sharp_greater_twostep()` obtains one-sided confidence interval for the maximum individual effect.
#'
#' @inheritParams pval_sharp_twostep
#'
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A numerical object that specifies the precision of the obtained confidence intervals.
#' For example, if `tol = 10^(-3)`, then the confidence limits are precise up to 3 digits.
#'
#' @return The one-sided confidence interval for the maximum individual effect.
#' @export
ci_sharp_greater_twostep <- function(Z, Y, missing = "general", method.list = list(name = "Wilcoxon"),
                                     stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, beta = 0.1 * 0.05, tol = 10^(-3)) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z) # number of observations
  n1 <- sum(Z) # number of treated units

  if (missing %in% c("general", "mp")) {
    class <- "MWU+"
  }
  if (missing == "mn") {
    class <- "MWU-"
  }

  # Generate Z.perm as matrix of possible treatment assignments with nperm combinations
  if (is.null(Z.perm)) {
    Z.perm <- assign_CRE(n, n1, nperm)
  }
  nperm <- ncol(Z.perm)
  # Z.perm: n x nperm matrix in which each column is one permutation of treatment assignment.

  # Calculate the test statistic for each treatment assignment permutation
  # As nperm -> infinity, stat.null can approximate to the null distribution G0()
  if (is.null(stat.null)) {
    stat.null <- null_dist(n, n1, class = class, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  f <- function(c) {
    pval <- pval_sharp_twostep(Z = Z, Y = Y, c = c, missing = missing, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, beta = beta)
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
#' `pval_sharp_twostep()` obtains the p-value for testing the sharp null hypothesis
#' \eqn{H_0: \tau = c}, using the two-step procedure.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param missing A string specifying the missing mechanism.
#' * `missing = "general"`: general missing mechanism.
#' * `missing = "mp"`: monotone positive missing mechanism (\eqn{M_1 >= M_0}).
#' * `missing = "mn"`: monotone negative missing mechanism (\eqn{M_1 <= M_0}).
#' @param method.list A list that specifies the choice of test statistic:
#' `method.list` can be `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
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
#' @param beta A real number that belongs to \eqn{[0, \alpha]}.
#' @return The p-value for testing the specified sharp null hypothesis of interest.
#' @export
pval_sharp_twostep <- function(Z, Y, c = 0, missing, method.list = list(name = "Wilcoxon"),
                               stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta) {
  if (missing == "general") {
    return(pval_sharp_g_twostep(Z, Y, c, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, beta = beta))
  }
  if (missing == "mp") {
    return(pval_sharp_mp_twostep(Z, Y, c, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, beta = beta))
  }
  if (missing == "mn") {
    return(pval_sharp_mn_twostep(Z, Y, c, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, beta = beta))
    # return(pval_sharp_mp_twostep(1 - Z, -Y, c, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, beta = beta))
  }
}

#' Confidence interval assuming constant treatment effect
#'
#' `ci_sharp_twostep()` obtains one-sided or two-sided confidence interval for the maximum individual effect.
#'
#' @inheritParams pval_sharp_twostep
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
ci_sharp_twostep <- function(Z, Y, alternative, missing = "general", method.list = list(name = "Wilcoxon"),
                             stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, beta = 0.1 * 0.05, tol = 10^(-3)) {
  if (alternative == "greater") {
    ci.lower <- ci_sharp_greater_twostep(Z = Z, Y = Y, missing = missing, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha, beta = beta, tol = tol)
    return(ci.lower)
  }
  if (alternative == "less") {
    ci.upper <- -1 * ci_sharp_greater_twostep(Z = Z, Y = -Y, missing = missing, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha, beta = beta, tol = tol)
    return(ci.upper)
  }
  if (alternative == "two.sided") {
    ci.lower <- ci_sharp_greater_twostep(Z = Z, Y = Y, missing = missing, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha / 2, beta = beta / 2, tol = tol)
    ci.upper <- -1 * ci_sharp_greater_twostep(Z = Z, Y = -Y, missing = missing, method.list = method.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm, alpha = alpha / 2, beta = beta / 2, tol = tol)
    return(c(ci.lower, ci.upper))
  }
}
