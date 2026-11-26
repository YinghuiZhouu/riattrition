% Generated manually for applied-user output wrapper
\name{ri_sharp_attrition}
\alias{ri_sharp_attrition}
\alias{ri_test}
\title{Applied-user wrapper for attrition inference output}
\usage{
ri_sharp_attrition(
  Z,
  Y,
  c = 0,
  missing = "general",
  X = NULL,
  class = "RS",
  method.list = list(name = "Wilcoxon"),
  cluster = NULL,
  cluster_outcome = "mean_observed",
  cluster_missing = "all_missing",
  alternative = "two.sided",
  stat.null = NULL,
  Z.perm = NULL,
  nperm = 10^4,
  alpha = 0.05,
  tol = 10^(-3),
  include_ci = TRUE,
  include_twostep = FALSE,
  beta = 0.1 * alpha,
  outcome_name = NULL
)
}
\arguments{
\item{Z}{Treatment assignment (an \eqn{n \times 1} vector).}

\item{Y}{Observed outcome (an \eqn{n \times 1} vector, including \code{NA}s).}

\item{c}{A scalar that specifies the sharp null hypothesis.}

\item{missing}{A string that specifies the missing mechanism.}

\item{X}{Optional data frame (or matrix) of numeric pretreatment covariates,
one row per unit, in the same order as \code{Z} and \code{Y}. It is used only
to check missing at random (see Details); it is not used in the main
analysis.}

\item{class}{A string that specifies the class of test statistic for the sharp-null procedure.}

\item{method.list}{A list that specifies the choice of test statistic.}

\item{cluster}{An optional cluster identifier vector. If supplied, the analysis
treats clusters as the unit of analysis.}

\item{cluster_outcome}{A string that specifies how individual outcomes are
aggregated to the cluster level. Currently only \code{"mean_observed"} is
supported.}

\item{cluster_missing}{A string that specifies when the cluster-level outcome is
missing. \code{"all_missing"} treats a cluster as missing only if all
individual outcomes are missing. \code{"any_missing"} treats a cluster as
missing if any individual outcome is missing.}

\item{alternative}{A string that specifies the direction of the confidence interval.}

\item{stat.null}{An optional null distribution for the sharp-null procedure.}

\item{Z.perm}{An optional assignment matrix for the sharp-null procedure.}

\item{nperm}{A positive integer that specifies the number of Monte Carlo (MC) replications, i.e. random re-assignments drawn to approximate the randomization distribution (default 10^4).}

\item{alpha}{A scalar, where \eqn{1-\alpha} indicates the confidence level.}

\item{tol}{A numerical object that specifies the precision of confidence intervals.}

\item{include_ci}{A logical value indicating whether confidence intervals should be computed.}

\item{include_twostep}{A logical value indicating whether two-step rows should be added for \code{general}, \code{mp}, and \code{mn} (default \code{FALSE}).}

\item{beta}{A real number in \eqn{[0, \alpha]} for the two-step procedure.}

\item{outcome_name}{Name of the outcome shown in the printed output. By default it is taken from the \code{Y} argument, e.g. \code{data$y} gives \code{y}.}
}
\value{
An object of class \code{riattrition_result}. When \code{missing} is
\code{"mp"}, \code{"mn"} or \code{"sharp"} and the design is complete randomization of
individuals, its component \code{monotone_pval} holds the p-value of the test
of that assumption from \code{\link{test_monotone_missingness}()}, and the
note on the assumption reports it. When \code{missing = "random"} and
\code{X} is given, its component \code{random_pval} holds the p-value of
\code{\link{test_missing_at_random}()} in the same way. Otherwise these
components are \code{NULL}.
}
\description{
\code{ri_sharp_attrition()} wraps the existing sharp-null attrition inference
functions and returns a richer result object with a default printed summary.
The printed output is a compact plain-text report: a header (outcome, design,
test, sharp null), outcome missingness overall and by arm, one row of Monte
Carlo p-value and confidence interval per missingness assumption, and notes
explaining the assumptions shown. Individual-level data are analyzed as a
completely randomized experiment.
If \code{cluster} is supplied, individual-level data are first aggregated to
cluster-level outcomes, and inference treats clusters as the unit of analysis.
}
\details{
The covariates \code{X} are used only to check the assumption of missing at
random. When \code{missing = "random"} and \code{X} is given,
\code{ri_sharp_attrition()} calls \code{test_missing_at_random(Z, Y, X)},
which tests whether, among the units with an observed outcome, the covariates
are balanced between the treatment and the control group (a consequence of
missing at random that can be checked), and adds its p-value to the note on the
assumption. Its p-value comes from \code{nperm} re-randomizations of the
treatment among those units, drawn after the table is computed. A small p-value
is evidence against missing at random; a large one means that we fail to reject
it, not that it holds.

\code{X} is not used in the main analysis: it does not adjust for the
covariates, and the test statistics, p-values and confidence intervals in the
table are the same with or without \code{X}. \code{X} is ignored when
\code{missing} is \code{"general"}, \code{"mp"}, \code{"mn"} or
\code{"sharp"}, whose checks need no covariates, and when \code{cluster} is
supplied, because \code{test_missing_at_random()} treats individuals as the
units.
}
\examples{
Z <- c(1, 1, 1, 0, 0, 0)
Y <- c(NA, 8, 6, 3, 4, NA)
out <- ri_sharp_attrition(
  Z = Z, Y = Y, c = 0, missing = "general",
  class = "RS", method.list = list(name = "Wilcoxon"),
  nperm = 200, include_ci = FALSE, include_twostep = FALSE
)
print(out)

# Missing at random, with two covariates to check this assumption
set.seed(1)
n <- 200
Z <- rep(c(1, 0), each = n / 2)
X <- data.frame(age = rnorm(n), female = rbinom(n, 1, 0.5))
Y <- rnorm(n) + Z
Y[runif(n) < 0.2] <- NA
out_random <- ri_sharp_attrition(
  Z = Z, Y = Y, missing = "random", X = X,
  nperm = 200, include_ci = FALSE
)
print(out_random)
out_random$random_pval
}
