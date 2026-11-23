% Generated manually for the missingness regression diagnostics
\name{attrition_regressions}
\alias{attrition_regressions}
\alias{print.attrition_regressions}
\title{Missingness regression diagnostics}
\usage{
attrition_regressions(M, Z, X, labels = NULL)

\method{print}{attrition_regressions}(x, ...)
}
\arguments{
\item{M}{A 0/1 vector; \code{M = 1} if the unit's outcome is observed and
\code{M = 0} if it is missing (Li, Sheng and Yu 2025).}

\item{Z}{A 0/1 treatment vector (\code{Z = 1} for treatment).}

\item{X}{A data frame of numeric pretreatment covariates.}

\item{labels}{Optional named character vector of display labels for the
covariates, e.g. \code{c(female = "Female")}. Covariates without a label are
shown under their own name.}

\item{x}{An object returned by \code{attrition_regressions()}.}

\item{...}{Not used.}
}
\value{
An object of class \code{attrition_regressions}: a list whose main
components are four data frames of numbers.
\describe{
\item{\code{coefficients}}{One row per model and coefficient: \code{model},
\code{term}, \code{label}, \code{estimate}, \code{std_error} (HC1) and
\code{p_value}. The models are \code{control} and \code{treatment} (M on X
within each group), \code{reg_a} (M on X, all units) and \code{reg_b} (M on X,
Z and Z x X, all units).}
\item{\code{group_effects}}{One row per covariate (and the constant) and
\code{effect}: the control-group coefficient, the treatment-group coefficient,
and their \code{difference} (treatment minus control). The difference equals the
Z x X coefficient of \code{reg_b} (Z for the constant) and its standard error
and p-value come from \code{reg_b}.}
\item{\code{model_statistics}}{Sample size and adjusted R-squared of each
model.}
\item{\code{joint_tests}}{HC1-robust Wald tests in F form: \code{all_x} (all X
coefficients = 0 in \code{reg_a}), \code{slopes} (all Z x X = 0 in
\code{reg_b}) and \code{chow} (Z = 0 and all Z x X = 0 in \code{reg_b}).}
}
It also contains \code{covariates}, \code{labels}, \code{sparse} (0/1
covariates with fewer than 10 ones or zeros in a group) and \code{fits} (the
fitted models, including their HC1 covariance matrices).
}
\description{
Fits the two attrition diagnostic tables. Table 1 regresses the missingness
indicator \code{M} on the covariates \code{X} separately in the control group
and in the treatment group, and reports the treatment-minus-control difference
of each coefficient. Table 2 reports reg. (a), \code{M} on \code{X}, and reg.
(b), \code{M} on \code{X}, \code{Z} and \code{Z x X}.
}
\details{
All regressions are linear probability models estimated by OLS with HC1
(heteroskedasticity-robust) standard errors. Units with a missing value of
\code{M}, \code{Z} or any covariate are dropped. If a 0/1 covariate equals 1
(or 0) for fewer than 10 units in the control or treatment group, a warning is
given, because its coefficients and joint tests that include it are unreliable.

Use \code{\link{attrition_table}()} to build the display tables and
\code{\link{export_attrition_tables}()} to write them to a Markdown or LaTeX
file.
}
\examples{
set.seed(1)
n <- 400
X <- data.frame(female = rbinom(n, 1, 0.5), age = rnorm(n))
Z <- rbinom(n, 1, 0.5)
M <- rbinom(n, 1, 0.8)

ar <- attrition_regressions(M, Z, X, labels = c(female = "Female", age = "Age"))
ar
ar$group_effects
}
