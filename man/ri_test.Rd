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
  class = "RS",
  method.list = list(name = "Wilcoxon"),
  block = NULL,
  block_missing_sum = FALSE,
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
  include_twostep = NULL,
  beta = 0.1 * alpha
)
}
\arguments{
\item{Z}{Treatment assignment (an \eqn{n \times 1} vector).}

\item{Y}{Observed outcome (an \eqn{n \times 1} vector, including \code{NA}s).}

\item{c}{A scalar that specifies the sharp null hypothesis.}

\item{missing}{A string that specifies the missing mechanism.}

\item{class}{A string that specifies the class of test statistic for the sharp-null procedure.}

\item{method.list}{A list that specifies the choice of test statistic.}

\item{block}{An optional block identifier vector. If supplied, permutations are
generated within blocks while preserving the treated count in each block.}

\item{block_missing_sum}{A logical value indicating whether the printed block
missingness section should include additional summaries beyond the default
two-line display.}

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

\item{nperm}{A positive integer that specifies the number of permutations.}

\item{alpha}{A scalar, where \eqn{1-\alpha} indicates the confidence level.}

\item{tol}{A numerical object that specifies the precision of confidence intervals.}

\item{include_ci}{A logical value indicating whether confidence intervals should be computed.}

\item{include_twostep}{A logical value indicating whether the two-step procedure should be included.}

\item{beta}{A real number in \eqn{[0, \alpha]} for the two-step procedure.}
}
\value{
An object of class \code{riattrition_result}.
}
\description{
\code{ri_sharp_attrition()} wraps the existing sharp-null attrition inference
functions and returns a richer result object with a default printed summary.
The printed output emphasizes missingness counts and attrition rates, reports
the observed test statistic, and omits regression-style columns such as
standard errors and z-statistics. If \code{block} is supplied, the wrapper
uses block-respecting permutation draws for the randomization distribution and
reports simple block-level attrition diagnostics in the printed output.
If \code{cluster} is supplied, individual-level data are first aggregated to
cluster-level outcomes, and inference treats clusters as the unit of analysis.
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
}
