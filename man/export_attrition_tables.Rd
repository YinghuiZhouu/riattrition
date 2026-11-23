% Generated manually for the missingness regression diagnostics
\name{export_attrition_tables}
\alias{export_attrition_tables}
\title{Write the missingness diagnostics as a Markdown or LaTeX file}
\usage{
export_attrition_tables(
  x,
  format = c("markdown", "latex"),
  file = NULL,
  stars = FALSE,
  digits = 3
)
}
\arguments{
\item{x}{An object returned by \code{\link{attrition_regressions}()}.}

\item{format}{\code{"markdown"} or \code{"latex"}.}

\item{file}{Path of the file to write. If \code{NULL}, nothing is written and
the text is returned.}

\item{stars}{Add significance stars to the estimates in the display tables.}

\item{digits}{Number of decimals in the display tables.}
}
\value{
The text of the file as a character string (invisibly when \code{file} is
given).
}
\description{
Writes both diagnostic tables (Table 1, by experimental group, and Table 2,
pooled) with their notes.

The LaTeX output contains two \code{table} environments that can be pasted
into a paper.

The Markdown output is meant to be read by people and by AI tools. It contains
definitions of all quantities, the two display tables with their notes, and the
underlying numbers in long format: coefficients, group effects, model
statistics and joint tests, with estimates, standard errors and p-values in
separate columns.
}
\examples{
set.seed(1)
n <- 400
X <- data.frame(female = rbinom(n, 1, 0.5), age = rnorm(n))
Z <- rbinom(n, 1, 0.5)
M <- rbinom(n, 1, 0.8)
ar <- attrition_regressions(M, Z, X, labels = c(female = "Female", age = "Age"))

md_file <- tempfile(fileext = ".md")
export_attrition_tables(ar, format = "markdown", file = md_file)

tex_file <- tempfile(fileext = ".tex")
export_attrition_tables(ar, format = "latex", file = tex_file)
}
