% Generated manually for the test of missing at random
\name{test_missing_at_random}
\alias{test_missing_at_random}
\title{Test of missing at random}
\usage{
test_missing_at_random(Z, Y, X, nperm = 10^4)
}
\arguments{
\item{Z}{Treatment assignment (0/1 vector).}

\item{Y}{Outcome vector; \code{NA} marks a missing outcome.}

\item{X}{Data frame (or matrix) of numeric pretreatment covariates, one row
per unit, in the same order as \code{Z} and \code{Y}. Units with a missing
covariate are left out.}

\item{nperm}{Number of re-randomizations of Z among the units with an observed
outcome (default 10^4).}
}
\value{
The p-value.
}
\description{
Randomization p-value of a balance test that checks missing at random: among
the units with an observed outcome, are the covariates \code{X} balanced
between the treatment and the control group?
}
\details{
Missing at random (\code{missing = "random"} in
\code{\link{ri_sharp_attrition}()}) is Assumption 4 of Li, Sheng and Yu (2025):
the pairs (M_i(1), M_i(0)), where M_i(z) = 1 if unit i's outcome is observed
under treatment z, are independent across units and all have the same
distribution, which does not depend on the potential outcomes. Treatment may
still change how many outcomes are missing, so the response rates of the two
groups may differ.

Under complete randomization, the assumption implies that the units with an
observed outcome form a completely randomized experiment among themselves
(their equation (19) and Theorem 5); \code{\link{ri_sharp_attrition}()} relies
on this when it drops the missing outcomes. The missing outcomes are not seen,
but the covariates are, so the function checks this implication with a balance
test:
\enumerate{
\item keep the units with an observed outcome (M = 1);
\item test whether the covariates are balanced between their treatment and
control groups;
\item if the balance test rejects, this is evidence against missing at random;
if it does not reject, we fail to reject missing at random.
}

\strong{Test statistic.} The F statistic of the OLS regression of Z on X among
the units with an observed outcome, for H0: all coefficients on X are zero.

\strong{p-value.} Under the assumption, the treatment assignment of the units
with an observed outcome is a complete randomization. The function therefore
re-randomizes Z among them \code{nperm} times, each time treating the same
number of them, and returns the share of re-randomizations whose statistic is
at least the observed one. Like \code{\link{test_monotone_missingness}()}, this
is a randomization test and needs no large-sample approximation. The p-value
changes a little with the random seed and with \code{nperm}.

\code{ri_sharp_attrition(Z, Y, missing = "random", X = X)} reports this
p-value in its note on the assumption.

A small p-value is evidence against the assumption. A large p-value does not
show that the assumption holds. The test does not detect missingness that
depends on the outcome itself or on characteristics that are not in \code{X},
nor missingness that depends on \code{X} but not on the group (the same
probability of being observed at each value of \code{X} in both groups): the
units with an observed outcome are then still balanced between the groups.
}
\examples{
set.seed(1)
n <- 400
Z <- rep(c(1, 0), each = n / 2)
X <- data.frame(age = rnorm(n), female = rbinom(n, 1, 0.5))
Y <- rnorm(n)

# Missing at random, with different response rates in the two groups.
Y_random <- Y
Y_random[runif(n) < ifelse(Z == 1, 0.1, 0.3)] <- NA
test_missing_at_random(Z, Y_random, X, nperm = 2000)

# Older units are less often observed among treated units, and more often
# among control units: the observed units are no longer balanced on age.
Y_age <- Y
Y_age[runif(n) > plogis(qlogis(0.8) + ifelse(Z == 1, -1, 1) * X$age)] <- NA
test_missing_at_random(Z, Y_age, X, nperm = 2000)
}
\seealso{
\code{\link{ri_sharp_attrition}()}, \code{\link{test_monotone_missingness}()},
\code{\link{attrition_regressions}()}
}
