% Generated manually for the missingness regression diagnostics
\name{attrition_table}
\alias{attrition_table}
\alias{attrition_table_note}
\title{Display tables for the missingness regressions}
\usage{
attrition_table(
  x,
  type = c("group", "pooled"),
  format = c("data.frame", "markdown", "latex"),
  stars = FALSE,
  digits = 3,
  caption = NULL
)

attrition_table_note(x, type = c("group", "pooled"), stars = FALSE)
}
\arguments{
\item{x}{An object returned by \code{\link{attrition_regressions}()}.}

\item{type}{\code{"group"} for Table 1 (regressions by experimental group,
with the treatment-minus-control difference, the slopes test and the Chow
test) or \code{"pooled"} for Table 2 (reg. (a) and reg. (b)).}

\item{format}{\code{"data.frame"} (for example for \code{knitr::kable()}),
\code{"markdown"} or \code{"latex"}. The last two return a character string
that includes the table note.}

\item{stars}{Add significance stars to the estimates (* p < 0.10,
** p < 0.05, *** p < 0.01).}

\item{digits}{Number of decimals.}

\item{caption}{Optional note for the Markdown or LaTeX output. By default the
note from \code{attrition_table_note()}.}
}
\value{
\code{attrition_table()} returns a data frame of formatted cells or a character
string. \code{attrition_table_note()} returns the table note as a character
string.
}
\description{
\code{attrition_table()} builds one of the two diagnostic tables from the
numbers stored by \code{\link{attrition_regressions}()}. Each cell shows an
estimate with its HC1 standard error in parentheses.
\code{attrition_table_note()} returns the matching table note, which defines
\code{M}, explains how to read the coefficients and reports the joint tests.
}
\examples{
set.seed(1)
n <- 400
X <- data.frame(female = rbinom(n, 1, 0.5), age = rnorm(n))
Z <- rbinom(n, 1, 0.5)
M <- rbinom(n, 1, 0.8)
ar <- attrition_regressions(M, Z, X, labels = c(female = "Female", age = "Age"))

attrition_table(ar, type = "group")
attrition_table(ar, type = "pooled")
cat(attrition_table(ar, type = "pooled", format = "latex"), "\n")
attrition_table_note(ar, type = "group")
}
