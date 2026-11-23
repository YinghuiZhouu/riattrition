% Generated manually for the test of monotone missingness
\name{test_monotone_missingness}
\alias{test_monotone_missingness}
\alias{print.monotone_missingness_test}
\title{Test of monotone missingness}
\usage{
test_monotone_missingness(
  Z,
  Y,
  assumption = c("mp", "mn"),
  block = NULL,
  nperm = 10^4
)

\method{print}{monotone_missingness_test}(x, ...)
}
\arguments{
\item{Z}{Treatment assignment (0/1 vector).}

\item{Y}{Outcome vector; \code{NA} marks a missing outcome. Only whether each
outcome is missing is used.}

\item{assumption}{\code{"mp"} tests H0: M_i(1) >= M_i(0) for every unit;
\code{"mn"} tests H0: M_i(1) <= M_i(0) for every unit. M_i = 1 if unit i's
outcome is observed (Li, Sheng and Yu 2025).}

\item{block}{Optional block identifier vector. If supplied, treatment is
reshuffled within blocks.}

\item{nperm}{Number of Monte Carlo replications. Used only when \code{block} is
supplied; under complete randomization the p-value is exact.}

\item{x}{An object returned by \code{test_monotone_missingness()}.}

\item{...}{Not used.}
}
\value{
An object of class \code{monotone_missingness_test}, a list with
\code{assumption}, \code{hypothesis}, \code{design}, \code{method}
(\code{"exact"} or \code{"Monte Carlo"}), \code{mc_replications},
\code{response_rate_treat}, \code{response_rate_control}, \code{statistic}
(the treatment-minus-control difference in response rates; block-size
weighted under blocked randomization) and \code{p_value}.
}
\description{
Tests the monotone missingness assumption used by
\code{ri_sharp_attrition(missing = "mp")} or \code{missing = "mn"}.
\code{ri_sharp_attrition()} runs this test automatically for these two
assumptions and prints its p-value.
}
\details{
The missingness indicator M is observed for every unit, so the effect of the
treatment on M can be tested like the effect on any other outcome. The
assumption \code{"mp"} states that every unit-level effect M_i(1) - M_i(0) is
non-negative, which is a bounded null hypothesis. A randomization test of no
effect on M with an effect-increasing statistic, here the difference in
response rates, is valid for this bounded null (Caughey, Dafoe, Li and
Miratrix 2023). The test is one-sided: for \code{"mp"}, a lower response rate
in the treatment group is evidence against the assumption.

Under complete randomization, the p-value is exact (a one-sided Fisher exact
test). Under blocked randomization, the statistic is the block-size weighted
difference in response rates over the blocks that contain both arms, and the
p-value is computed by Monte Carlo, reshuffling treatment within blocks.

A small p-value is evidence against the assumption. A large p-value does not
show that the assumption holds: the test only sees the net difference in
response rates, so units that violate the assumption can be offset by units
with an effect in the other direction.
}
\references{
Caughey, D., Dafoe, A., Li, X. and Miratrix, L. (2023). Randomization
inference beyond the sharp null: bounded null hypotheses and quantiles of
individual treatment effects. Journal of the Royal Statistical Society Series
B, 85(5), 1471-1491.

Li, X., Sheng, P. and Yu, Z. (2025). Randomization inference with sample
attrition. arXiv:2507.00795.
}
\examples{
Z <- c(rep(1, 50), rep(0, 50))
Y <- rnorm(100)
Y[c(1:5, 51:60)] <- NA   # 10\% missing under treatment, 20\% under control
test_monotone_missingness(Z, Y, assumption = "mp")
test_monotone_missingness(Z, Y, assumption = "mn")
}
