% Generated manually for null-distribution visualization
\name{plot_ri_test_stat}
\alias{plot_ri_test_stat}
\title{Visualize the null distribution under general missingness}
\usage{
plot_ri_test_stat(
  Z,
  Y,
  c = 0,
  missing = "general",
  class = "RS",
  analysis = c("general", "complete-case"),
  method.list = list(name = "Wilcoxon"),
  block = NULL,
  cluster = NULL,
  cluster_outcome = "mean_observed",
  cluster_missing = "all_missing",
  nperm = 10^4,
  style = c("histogram", "density"),
  bins = "FD",
  show_density = TRUE,
  off_support = c("auto", "arrow", "extend", "split"),
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  null_col = "#D7DDE3",
  tail_col = "#C9B8B4",
  density_col = "#3E5F7E",
  stat_col = "#8F3B3B",
  border = "white",
  lwd = 2,
  ...
)
}
\arguments{
\item{Z}{Treatment assignment (an \eqn{n \times 1} vector).}

\item{Y}{Observed outcome (an \eqn{n \times 1} vector, including \code{NA}s).}

\item{c}{A scalar that specifies the sharp null hypothesis.}

\item{analysis}{A string that specifies which inferential path to visualize:
\code{"general"} uses the package's general-missingness procedure, while
\code{"complete-case"} restricts the analysis to observed outcomes only.}

\item{missing}{A string that specifies the missing mechanism. The current
visualization method uses \code{"general"} when
\code{analysis = "general"}. It is ignored when
\code{analysis = "complete-case"}.}

\item{class}{A string that specifies the class of test statistic for the
sharp-null procedure.}

\item{method.list}{A list that specifies the choice of test statistic.}

\item{block}{An optional block identifier vector. If supplied, permutations are
generated within blocks while preserving the treated count in each block.}

\item{cluster}{An optional cluster identifier vector. If supplied, the analysis
treats clusters as the unit of analysis.}

\item{cluster_outcome}{A string that specifies how individual outcomes are
aggregated to the cluster level. Currently only \code{"mean_observed"} is
supported.}

\item{cluster_missing}{A string that specifies when the cluster-level outcome is
missing.}

\item{nperm}{A positive integer that specifies the number of permutations.}

\item{style}{Plot style. \code{"histogram"} is the default because it directly
represents the empirical permutation distribution. \code{"density"} is
available as a smoother alternative.}

\item{bins}{Histogram break rule passed to \code{graphics::hist()} when
\code{style = "histogram"}.}

\item{show_density}{A logical value indicating whether a density curve should
be overlaid when \code{style = "histogram"}.}

\item{off_support}{A string that specifies how to display the observed test
statistic when it falls outside the support of the permutation null
distribution. \code{"auto"} extends the axis for modest gaps and uses a split
axis for large gaps. \code{"arrow"} keeps the current support and points
toward the off-support value. \code{"extend"} enlarges the x-axis to include
the observed statistic. \code{"split"} uses a split-axis display so that the
observed statistic remains visible without compressing the null distribution.}

\item{main}{Optional main title. If \code{NULL}, the function supplies a
default.}

\item{xlab}{Optional x-axis label.}

\item{ylab}{Optional y-axis label.}

\item{null_col}{Fill color for the main null distribution.}

\item{tail_col}{Fill color for the right-tail region used in the p-value.}

\item{density_col}{Line color for the optional density overlay.}

\item{stat_col}{Color of the observed-statistic reference line.}

\item{border}{Border color for histogram bars.}

\item{lwd}{Width of the observed-statistic reference line.}

\item{...}{Additional graphical arguments passed to the underlying plotting
function.}
}
\value{
Invisibly returns an object of class \code{riattrition_plot} containing the
null distribution, observed statistic, p-value, and plotting metadata.
}
\description{
\code{plot_ri_test_stat()} draws the randomization null distribution, marks the
observed test statistic, and highlights the right tail that determines the
sharp-null p-value. It supports the package's general-missingness procedure
and an observed-only complete-case visualization mode.
}
\examples{
Z <- c(1, 1, 1, 0, 0, 0)
Y <- c(NA, 8, 6, 3, 4, NA)
plot_ri_test_stat(
  Z = Z,
  Y = Y,
  nperm = 200
)
}
