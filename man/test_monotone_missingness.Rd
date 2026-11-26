% Generated manually for the test of monotone missingness
\name{test_monotone_missingness}
\alias{test_monotone_missingness}
\title{Test of monotone missingness}
\usage{
test_monotone_missingness(Z, Y, assumption = "mp", Z.perm = NULL, nperm = 10^4)
}
\arguments{
\item{Z}{Treatment assignment (0/1 vector).}

\item{Y}{Outcome vector; \code{NA} marks a missing outcome.}

\item{assumption}{\code{"mp"}, \code{"mn"} or \code{"sharp"}.}

\item{Z.perm}{Optional matrix of random assignments with one row per unit and
one column per assignment. If \code{NULL}, \code{nperm} complete
randomizations are drawn.}

\item{nperm}{Number of random assignments drawn when \code{Z.perm} is
\code{NULL}.}
}
\value{
The p-value.
}
\description{
Randomization p-value for H0: M_i(1) >= M_i(0) for every unit (\code{"mp"}),
H0: M_i(1) <= M_i(0) for every unit (\code{"mn"}) or H0: M_i(1) = M_i(0) for
every unit (\code{"sharp"}), where M_i = 1 if unit i's outcome is observed,
under complete randomization.
}
\details{
mp and mn are bounded null hypotheses about the effect of treatment on M.
Following Caughey, Dafoe, Li and Miratrix (2023), each is tested with the
one-sided randomization test of "no effect on M", using the difference in
means of M (treated minus control) as the test statistic. The p-value is the
share of the random assignments in \code{Z.perm} whose statistic is at least as
extreme as the observed one. The code follows \code{pval_imp_control()} in
Xinran Li's RIQITE package: with the same \code{Z.perm}, the p-value equals
that of RIQITE's
\code{pval_bound(Z, M, c = 0, method.list = list(name = "DIM"))} with
\code{alternative = "less"} for mp and \code{"greater"} for mn.

sharp is the sharp null of no effect on M itself, and it holds exactly when mp
and mn both hold. Evidence against it can go either way, so its p-value is
two-sided: the share of the random assignments whose statistic is at least as
far from 0 as the observed one. RIQITE's \code{pval_bound()} has no two-sided
version.

A small p-value is evidence against the assumption. A large p-value does not
show that the assumption holds: units that violate it can be offset by units
with an effect in the other direction.

\code{\link{ri_sharp_attrition}()} reports this p-value in its notes when
\code{missing = "mp"}, \code{"mn"} or \code{"sharp"} and the design is complete randomization
of individuals. It passes the random assignments used for its own tests, so
the test draws no new random numbers.
}
\references{
Caughey, D., Dafoe, A., Li, X. and Miratrix, L. (2023). Randomization inference
beyond the sharp null: bounded null hypotheses and quantiles of individual
treatment effects. arXiv:2101.09195.
}
\examples{
set.seed(1)
Z <- rep(c(1, 0), each = 50)
Y <- rnorm(100)
Y[runif(100) < ifelse(Z == 1, 0.1, 0.3)] <- NA
test_monotone_missingness(Z, Y, assumption = "mp")
test_monotone_missingness(Z, Y, assumption = "mn")
test_monotone_missingness(Z, Y, assumption = "sharp")
}
