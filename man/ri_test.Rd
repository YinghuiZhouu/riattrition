% Generated manually for applied-user output wrapper
\name{ri_test}
\alias{ri_test}
\title{Applied-user wrapper for attrition inference output}
\usage{
ri_test(
  Z,
  Y,
  c = 0,
  missing = "general",
  class = "RS",
  method.list = list(name = "Wilcoxon"),
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
\code{ri_test()} wraps the existing sharp-null attrition inference functions and
returns a richer result object with a default printed summary. The printed
output emphasizes missingness counts, reports the observed test statistic, and
omits regression-style columns such as standard errors and z-statistics.
}
\examples{
Z <- c(1, 1, 1, 0, 0, 0)
Y <- c(NA, 8, 6, 3, 4, NA)
out <- ri_test(
  Z = Z, Y = Y, c = 0, missing = "general",
  class = "RS", method.list = list(name = "Wilcoxon"),
  nperm = 200, include_ci = FALSE, include_twostep = FALSE
)
print(out)
}
