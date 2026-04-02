
<!-- README.md is generated from README.Rmd. Please edit that file -->

# riattrition

<!-- badges: start -->

<!-- badges: end -->

## Overview

The `riattrition` package provides a set of functions for conducting
randomization inference in randomized experiments with sample attrition.
It constructs valid p-values and confidence intervals for testing sharp
null hypotheses when outcomes are missing for some units. The package
accommodates a range of missing mechanisms, including general, monotone
positive, monotone negative, sharp, and missing-at-random mechanisms.

## Installation

`riattrition` can be installed from its GitHub repository via

``` r
devtools::install_github("peizansheng/riattrition")
```

## Main Functions

This package provides the following functions to conduct randomization
inference with sample attrition.

`pval_sharp()` obtains the p-value for testing the sharp null hypothesis
$H_0: \tau = c$.

``` r
pval_sharp(Z, Y, c = 0, missing = "general", class = "RS", method.list = list(name = "Wilcoxon"),
           stat.null = NULL, Z.perm = NULL, nperm = 10^4)
ci_sharp(Z, Y, alternative, missing = "general", class = "RS", method.list = list(name = "Wilcoxon"),
         stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3))
```

where:

- `Z`: Treatment assignment ($n \times 1$ vector).
- `Y`: Observed outcome ($n \times 1$ vector, including `NA`s).
- `c`: A scalar that specifies the sharp null hypothesis.
- `missing`: A string that specifies the missing mechanism:
  - `missing = "general"`: general missing mechanism.
  - `missing = "mp"`: monotone positive missing mechanism
    ($M_1 \geq M_0$).
  - `missing = "mn"`: monotone negative missing mechanism
    ($M_1 \leq M_0$).
  - `missing = "sharp"`: sharp missing mechanism ($M_1 = M_0$).
  - `missing = "random"`: missing at random mechanism.
- `class`: A string that specifies the class of test statistic:
  - `class = "RS"` (first class, default): rank-sum statistic.
  - `class = "MWU+"` (second class): generalized Mann-Whitney U
    statistic, defined as the sum of relative ranks for treated units.
  - `class = "MWU-"` (second class): generalized Mann-Whitney U
    statistic, defined as the negative sum of relative ranks for control
    units.
- `method.list`: A list that specifies the choice of test statistic:
  - if `class = "RS"`, `method.list` can be `list(name = "Wilcoxon")` or
    `list(name = "Stephenson", s = 10)`.
  - if `class = "MWU+"` or `class = "MWU-"`, `method.list` can be
    `list(name = "Wilcoxon")` or `list(name = "Polynomial", s = 10)`.
- `stat.null`: An $nperm \times 1$ vector whose empirical distribution
  approximates the randomization distribution of the rank-based
  statistic.
  - if `stat.null = NULL` (default), the function will calculate it
    internally.
  - if `missing = "general"`, `mp`, `mn`, then $n$ choose $n_1$.
  - if `missing = "sharp"`, `random`, then $n_{11} + n_{01}$ choose
    $n_{11}$.
- Z.perm: A matrix that specifies the permuted assignments for
  approximating the null distribution of the test statistic.
  - if `Z.perm = NULL` (default), the function will calculate it
    internally.
  - if `missing = "general"`, `mp`, `mn`, then `Z.perm` is an
    $n \times nperm$ matrix ($n$ choose $n_1$).
  - if `missing = "sharp"`, `random`, then `Z.perm` is an
    $(n_{11} + n_{01}) \times nperm$ matrix ($n_{11} + n_{01}$ choose
    $n_{11}$).
- `nperm`: A positive integer that specifies the number of permutations
  to approximate the randomization distribution of the test statistic.

`ci_sharp()` obtains one-sided or two-sided confidence interval assuming
constant treatment effect.

``` r
ci_sharp(Z, Y, alternative, missing = "general", class = "RS", method.list = list(name = "Wilcoxon"),
         stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, tol = 10^(-3))
```

`pval_sharp_twostep()` obtains the p-value for testing the sharp null
hypothesis $H_0: \tau = c$ using the two-step procedure.

``` r
pval_sharp_twostep(Z, Y, c = 0, missing, method.list = list(name = "Wilcoxon"),
                   stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta)
```

`ci_sharp_twostep()` obtains one-sided or two-sided confidence interval
assuming constant treatment effect using the two-step procedure.

``` r
ci_sharp_twostep(Z, Y, alternative, missing = "general", method.list = list(name = "Wilcoxon"),
                 stat.null = NULL, Z.perm = NULL, nperm = 10^4, alpha = 0.05, beta = 0.1 * 0.05, tol = 10^(-3))
```

## Usage

We will use a simulated data to illustrate the usage of the main
functions.

``` r
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
p_value_g <- pval_sharp(Z, Y_observed, c = 0, missing = "general", MWU = FALSE,
                        method.list = list(name = "Wilcoxon"),
                        stat.null = NULL, Z.perm = NULL, nperm = 10^4)
p_value_g_twostep <- pval_sharp_twostep(Z, Y_observed, c = 0, missing = "general",
                                        method.list = list(name = "Wilcoxon"),
                                        stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta = 0.1 * 0.1)

# (2) Monotone Positive Missingness
set.seed(1)
M0 <- as.numeric(Y0 <= qnorm(0.92))
M1 <- as.numeric(Y0 <= qnorm(0.98))
Z <- sample(c(rep(1, N/2), rep(0, N/2)))
M <- Z * M1 + (1 - Z) * M0
Y <- Z * Y1 + (1 - Z) * Y0 # No missing value
Y_observed <- ifelse(M == 0, NA, Y)
p_value_mp <- pval_sharp(Z, Y_observed, c = 0, missing = "mp", MWU = FALSE,
                         method.list = list(name = "Wilcoxon"),
                         stat.null = NULL, Z.perm = NULL, nperm = 10^4)
p_value_mp_twostep <- pval_sharp_twostep(Z, Y_observed, c = 0, missing = "mp",
                                        method.list = list(name = "Wilcoxon"),
                                        stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta = 0.1 * 0.1)

# (3) Monotone Negative Missingness
set.seed(1)
M1 <- as.numeric(Y0 >= qnorm(0.08))
M0 <- as.numeric(Y0 >= qnorm(0.02))
Z <- sample(c(rep(1, N/2), rep(0, N/2)))
M <- Z * M1 + (1 - Z) * M0
Y <- Z * Y1 + (1 - Z) * Y0 # No missing value
Y_observed <- ifelse(M == 0, NA, Y)
p_value_mn <- pval_sharp(Z, Y_observed, c = 0, missing = "mn", MWU = FALSE,
                         method.list = list(name = "Wilcoxon"),
                         stat.null = NULL, Z.perm = NULL, nperm = 10^4)
p_value_mn_twostep <- pval_sharp_twostep(Z, Y_observed, c = 0, missing = "mn",
                                         method.list = list(name = "Wilcoxon"),
                                         stat.null = NULL, Z.perm = NULL, nperm = 10^4, beta = 0.1 * 0.1)
                                         
# (4) Sharp Missingness
set.seed(1)
M0 <- as.numeric(Y0 <= qnorm(0.95))
M1 <- M0
Z <- sample(c(rep(1, N/2), rep(0, N/2)))
M <- Z * M1 + (1 - Z) * M0
Y <- Z * Y1 + (1 - Z) * Y0 # No missing value
Y_observed <- ifelse(M == 0, NA, Y)
p_value_s <- pval_sharp(Z, Y_observed, c = 0, missing = "sharp", MWU = FALSE,
                        method.list = list(name = "Wilcoxon"),
                        stat.null = NULL, Z.perm = NULL, nperm = 10^4)

# (5) Missing at random
set.seed(1)
p <- 0.95
M0 <- rbinom(N, 1, p)
M1 <- rbinom(N, 1, p)
Z <- sample(c(rep(1, N/2), rep(0, N/2)))
M <- Z * M1 + (1 - Z) * M0
Y <- Z * Y1 + (1 - Z) * Y0 # No missing value
Y_observed <- ifelse(M == 0, NA, Y)
p_value_r <- pval_sharp(Z, Y_observed, c = 0, missing = "random", MWU = FALSE,
                        method.list = list(name = "Wilcoxon"),
                        stat.null = NULL, Z.perm = NULL, nperm = 10^4)
```
