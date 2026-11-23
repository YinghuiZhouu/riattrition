#----------- Missingness regression diagnostics -------------------------------#
#
# Two diagnostic tables for attrition. Notation follows Li, Sheng & Yu (2025):
#   * M_i = 1 if unit i's outcome is observed, M_i = 0 if it is missing;
#   * Z_i = 1 if unit i is assigned to treatment, Z_i = 0 for control;
#   * X_i are pretreatment covariates.
#
# Table 1 (by experimental group): regress M on X separately in the control
#   group and in the treatment group, and report the treatment-minus-control
#   difference of each coefficient.
# Table 2 (pooled): reg. (a) M on X, and reg. (b) M on X, Z and Z x X.
#
# All regressions are linear probability models fitted with lm(), with HC1
# (heteroskedasticity-robust) standard errors computed below. Only base R is
# used.
#
# The results are stored as numbers in four data frames (coefficients,
# group_effects, model_statistics, joint_tests). The display tables, the LaTeX
# code and the Markdown file are all built from these data frames.


#----------- Step 1: OLS with HC1 standard errors -----------------------------#

# Fits `formula` by OLS and returns coefficients, HC1 standard errors,
# two-sided p-values, the HC1 covariance matrix, n and the adjusted R-squared.
# HC1 = sandwich covariance multiplied by n / (n - k).
.ar_fit <- function(formula, data) {
  fit <- stats::lm(formula, data = data)

  # Drop coefficients that lm() could not estimate (NA), e.g. a constant column.
  b <- stats::coef(fit)
  keep <- !is.na(b)
  b <- b[keep]
  X <- stats::model.matrix(fit)[, keep, drop = FALSE]
  e <- stats::residuals(fit)

  n <- nrow(X)
  k <- ncol(X)
  bread <- solve(crossprod(X))
  meat <- crossprod(X * e)
  V <- n / (n - k) * bread %*% meat %*% bread
  se <- sqrt(diag(V))

  list(
    coef = b,
    se = se,
    p = 2 * stats::pt(-abs(b / se), df = n - k),
    vcov = V,
    n = n,
    df_resid = n - k,
    adj_r2 = summary(fit)$adj.r.squared
  )
}


#----------- Step 2: joint Wald test ------------------------------------------#

# HC1-robust Wald test of H0: all coefficients in `terms` are zero.
# Reported in F form: F = W / q with q and n - k degrees of freedom.
.ar_wald <- function(fit, terms) {
  terms <- intersect(terms, names(fit$coef))
  b <- fit$coef[terms]
  V <- fit$vcov[terms, terms, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(V, b))
  q <- length(terms)
  list(
    F = W / q,
    df1 = q,
    df2 = fit$df_resid,
    p = stats::pf(W / q, q, fit$df_resid, lower.tail = FALSE)
  )
}


#----------- Step 3: sparse 0/1 covariates ------------------------------------#

# A 0/1 covariate that equals 1 (or 0) for only a few units in a regression
# sample is fit almost exactly for those units, so its HC1 standard error and
# any joint test that includes it are unreliable. Such covariates are listed in
# the table notes and reported as an R warning.
.ar_min_cell <- 10

.ar_sparse <- function(X_sample, labels, sample_name) {
  out <- data.frame(sample = character(0), covariate = character(0), label = character(0),
                    n_ones = integer(0), n_zeros = integer(0), stringsAsFactors = FALSE)
  for (v in names(X_sample)) {
    x <- X_sample[[v]]
    if (!all(x %in% c(0, 1))) {
      next
    }
    n_ones <- sum(x == 1)
    n_zeros <- sum(x == 0)
    if (min(n_ones, n_zeros) < .ar_min_cell) {
      out <- rbind(out, data.frame(sample = sample_name, covariate = v, label = unname(labels[v]),
                                   n_ones = n_ones, n_zeros = n_zeros, stringsAsFactors = FALSE))
    }
  }
  if (nrow(out) > 0) {
    warning(.ar_sparse_text(out), call. = FALSE)
  }
  out
}

.ar_sparse_text <- function(sparse) {
  if (is.null(sparse) || nrow(sparse) == 0) {
    return("")
  }
  cells <- ifelse(
    sparse$n_ones < sparse$n_zeros,
    paste0(sparse$label, " equals 1 for only ", sparse$n_ones, " units in the ", sparse$sample),
    paste0(sparse$label, " equals 0 for only ", sparse$n_zeros, " units in the ", sparse$sample)
  )
  paste0(
    " Caution: ", paste(cells, collapse = "; "), " (fewer than ", .ar_min_cell,
    "), so the corresponding coefficients, their HC1 standard errors, and joint tests that include them are unreliable."
  )
}


#----------- Step 4: the main function ----------------------------------------#

# Display labels for covariates; names without a label keep their own name.
.ar_labels <- function(covars, labels = NULL) {
  out <- stats::setNames(covars, covars)
  if (!is.null(labels)) {
    hit <- intersect(covars, names(labels))
    out[hit] <- labels[hit]
  }
  out
}

# One row per coefficient of one fitted model, as numbers.
.ar_coef_rows <- function(fit, model, terms, labels) {
  data.frame(
    model = model,
    term = terms,
    label = labels,
    estimate = unname(fit$coef[terms]),
    std_error = unname(fit$se[terms]),
    p_value = unname(fit$p[terms]),
    stringsAsFactors = FALSE
  )
}

#' Missingness regression diagnostics
#'
#' Fits the two attrition diagnostic tables: regressions of the missingness
#' indicator `M` on covariates `X` separately for the control and treatment
#' groups, and the pooled regressions of `M` on `X` (reg. (a)) and of `M` on
#' `X`, `Z` and `Z x X` (reg. (b)).
#'
#' @param M A 0/1 vector; `M = 1` if the unit's outcome is observed.
#' @param Z A 0/1 treatment vector.
#' @param X A data frame of numeric pretreatment covariates.
#' @param labels Optional named character vector of display labels for the
#'   covariates, e.g. `c(female = "Female")`.
#'
#' @return An object of class `attrition_regressions`.
#' @export
attrition_regressions <- function(M, Z, X, labels = NULL) {
  # 1. Check the inputs.
  X <- as.data.frame(X)
  covars <- names(X)
  if (!all(M %in% c(0, 1, NA)) || !all(Z %in% c(0, 1, NA))) {
    stop("`M` and `Z` must be 0/1 vectors (M = 1 means the outcome is observed).")
  }
  labels <- .ar_labels(covars, labels)

  # 2. Keep units with M, Z and every covariate recorded.
  dat <- stats::na.omit(data.frame(M = M, Z = Z, X, check.names = FALSE))
  control <- dat[dat$Z == 0, , drop = FALSE]
  treated <- dat[dat$Z == 1, , drop = FALSE]

  # 3. Fit the four regressions.
  rhs <- paste(sprintf("`%s`", covars), collapse = " + ")
  m_on_x <- stats::as.formula(paste("M ~", rhs))
  m_on_xz <- stats::as.formula(paste0("M ~ Z * (", rhs, ")"))
  fits <- list(
    control = .ar_fit(m_on_x, control),
    treatment = .ar_fit(m_on_x, treated),
    reg_a = .ar_fit(m_on_x, dat),
    reg_b = .ar_fit(m_on_xz, dat)
  )

  # 4. Joint tests.
  interactions <- paste0("Z:", covars)
  all_x <- .ar_wald(fits$reg_a, covars)
  slopes <- .ar_wald(fits$reg_b, interactions)
  chow <- .ar_wald(fits$reg_b, c("Z", interactions))

  # 5. Collect the results as numbers.
  x_terms <- c(covars, "(Intercept)")
  x_labels <- c(unname(labels), "Constant")
  b_terms <- c(covars, "Z", interactions, "(Intercept)")
  b_labels <- c(unname(labels), "Treatment", paste("Treatment \u00d7", unname(labels)), "Constant")

  coefficients <- rbind(
    .ar_coef_rows(fits$control, "control", x_terms, x_labels),
    .ar_coef_rows(fits$treatment, "treatment", x_terms, x_labels),
    .ar_coef_rows(fits$reg_a, "reg_a", x_terms, x_labels),
    .ar_coef_rows(fits$reg_b, "reg_b", b_terms, b_labels)
  )

  # Group effects: the control and treatment coefficients come from the two
  # separate regressions; the difference (treatment - control) is the Z x X
  # coefficient of reg. (b) (Z for the constant), with its HC1 standard error.
  difference <- .ar_coef_rows(fits$reg_b, "reg_b", c(interactions, "Z"), x_labels)
  difference$term <- x_terms
  group_effects <- rbind(
    cbind(effect = "control", .ar_coef_rows(fits$control, "control", x_terms, x_labels)),
    cbind(effect = "treatment", .ar_coef_rows(fits$treatment, "treatment", x_terms, x_labels)),
    cbind(effect = "difference", difference)
  )
  group_effects$model <- NULL

  model_statistics <- data.frame(
    model = c("control", "treatment", "reg_a", "reg_b"),
    description = c(
      "M on X, control group (Z = 0)",
      "M on X, treatment group (Z = 1)",
      "reg. (a): M on X, all units",
      "reg. (b): M on X, Z and Z x X, all units"
    ),
    n = c(fits$control$n, fits$treatment$n, fits$reg_a$n, fits$reg_b$n),
    adj_r_squared = c(fits$control$adj_r2, fits$treatment$adj_r2, fits$reg_a$adj_r2, fits$reg_b$adj_r2),
    stringsAsFactors = FALSE
  )

  joint_tests <- data.frame(
    test = c("all_x", "slopes", "chow"),
    model = c("reg_a", "reg_b", "reg_b"),
    hypothesis = c(
      "all X coefficients = 0",
      "all Z x X coefficients = 0",
      "Z = 0 and all Z x X coefficients = 0"
    ),
    F = c(all_x$F, slopes$F, chow$F),
    df1 = c(all_x$df1, slopes$df1, chow$df1),
    df2 = c(all_x$df2, slopes$df2, chow$df2),
    p_value = c(all_x$p, slopes$p, chow$p),
    stringsAsFactors = FALSE
  )

  sparse <- rbind(
    .ar_sparse(control[covars], labels, "control group"),
    .ar_sparse(treated[covars], labels, "treatment group")
  )

  structure(
    list(
      covariates = covars,
      labels = labels,
      coefficients = coefficients,
      group_effects = group_effects,
      model_statistics = model_statistics,
      joint_tests = joint_tests,
      sparse = sparse,
      fits = fits
    ),
    class = "attrition_regressions"
  )
}


#----------- Step 5: formatting helpers ---------------------------------------#

.ar_adj_r2_label <- "Adjusted R\u00b2"

.ar_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

.ar_count <- function(n) {
  format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
}

.ar_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", ""))))
}

# "0.052 (0.044)", with optional significance stars after the estimate.
.ar_est_se <- function(est, se, p, stars = FALSE, digits = 3) {
  star_text <- if (stars) .ar_stars(p) else ""
  out <- paste0(.ar_num(est, digits), star_text, " (", .ar_num(se, digits), ")")
  out[is.na(est)] <- ""
  out
}

# "F(6, 771) = 0.41" with or without ", p = 0.875".
.ar_f_text <- function(test, with_p = FALSE) {
  out <- sprintf("F(%d, %s) = %.2f", as.integer(test$df1), .ar_count(test$df2), test$F)
  if (with_p) {
    out <- paste0(out, sprintf(", p = %.3f", test$p_value))
  }
  out
}

# Rows of a coefficient data frame for one model, in the order of `terms`.
.ar_pick <- function(coefficients, model, terms) {
  rows <- coefficients[coefficients$model == model, , drop = FALSE]
  rows[match(terms, rows$term), , drop = FALSE]
}


#----------- Step 6: display tables -------------------------------------------#

#' Display tables for the missingness regressions
#'
#' `attrition_table()` builds Table 1 (`type = "group"`: regressions by
#' experimental group) or Table 2 (`type = "pooled"`: reg. (a) and reg. (b)).
#' `attrition_table_note()` returns the matching table note.
#'
#' @param x An object returned by [attrition_regressions()].
#' @param type `"group"` or `"pooled"`.
#' @param format `"data.frame"` (e.g. for `knitr::kable()`), `"markdown"` or
#'   `"latex"`. The last two return a character string that includes the note.
#' @param stars Add significance stars to the estimates.
#' @param digits Number of decimals.
#' @param caption Optional note for the Markdown or LaTeX output; by default
#'   `attrition_table_note()`.
#'
#' @return A data frame or a character string.
#' @export
attrition_table <- function(x, type = c("group", "pooled"),
                            format = c("data.frame", "markdown", "latex"),
                            stars = FALSE, digits = 3, caption = NULL) {
  type <- match.arg(type)
  format <- match.arg(format)
  stats_row <- function(model) x$model_statistics[x$model_statistics$model == model, ]
  test_row <- function(test) x$joint_tests[x$joint_tests$test == test, ]

  if (type == "group") {
    ge <- x$group_effects
    control <- ge[ge$effect == "control", ]
    treatment <- ge[ge$effect == "treatment", ]
    difference <- ge[ge$effect == "difference", ]

    tab <- data.frame(
      Covariate = control$label,
      `Control group` = .ar_est_se(control$estimate, control$std_error, control$p_value, stars, digits),
      `Treatment group` = .ar_est_se(treatment$estimate, treatment$std_error, treatment$p_value, stars, digits),
      `Difference (T - C)` = .ar_est_se(difference$estimate, difference$std_error, difference$p_value, stars, digits),
      `p-value` = .ar_num(difference$p_value, digits),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    footer <- data.frame(
      Covariate = c("N", .ar_adj_r2_label, "Slopes test: all Z \u00d7 X = 0", "Chow test: Z = 0 and all Z \u00d7 X = 0"),
      `Control group` = c(.ar_count(stats_row("control")$n), .ar_num(stats_row("control")$adj_r_squared, digits), "", ""),
      `Treatment group` = c(.ar_count(stats_row("treatment")$n), .ar_num(stats_row("treatment")$adj_r_squared, digits), "", ""),
      `Difference (T - C)` = c("", "", .ar_f_text(test_row("slopes")), .ar_f_text(test_row("chow"))),
      `p-value` = c("", "", .ar_num(test_row("slopes")$p_value, digits), .ar_num(test_row("chow")$p_value, digits)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    tab <- rbind(tab, footer)
  } else {
    terms <- c(x$covariates, "Z", paste0("Z:", x$covariates), "(Intercept)")
    reg_a <- .ar_pick(x$coefficients, "reg_a", terms)
    reg_b <- .ar_pick(x$coefficients, "reg_b", terms)

    tab <- data.frame(
      Covariate = c(reg_b$label, "N", .ar_adj_r2_label,
                    "p-value: all Z \u00d7 X = 0", "p-value: Z and all Z \u00d7 X = 0"),
      reg_a = c(.ar_est_se(reg_a$estimate, reg_a$std_error, reg_a$p_value, stars, digits),
                .ar_count(stats_row("reg_a")$n), .ar_num(stats_row("reg_a")$adj_r_squared, digits),
                "", ""),
      reg_b = c(.ar_est_se(reg_b$estimate, reg_b$std_error, reg_b$p_value, stars, digits),
                .ar_count(stats_row("reg_b")$n), .ar_num(stats_row("reg_b")$adj_r_squared, digits),
                .ar_num(test_row("slopes")$p_value, digits), .ar_num(test_row("chow")$p_value, digits)),
      stringsAsFactors = FALSE
    )
    names(tab) <- c("Covariate", "Reg. (a): M on X", "Reg. (b): M on X, Z, and Z \u00d7 X")
  }
  rownames(tab) <- NULL

  if (format == "data.frame") {
    return(tab)
  }
  if (is.null(caption)) {
    caption <- attrition_table_note(x, type = type, stars = stars)
  }
  if (format == "markdown") {
    return(.ar_markdown(tab, caption))
  }
  .ar_latex(tab, caption)
}

#' @rdname attrition_table
#' @export
attrition_table_note <- function(x, type = c("group", "pooled"), stars = FALSE) {
  type <- match.arg(type)
  test_row <- function(test) x$joint_tests[x$joint_tests$test == test, ]
  dv <- paste(
    "Dependent variable: M_i = 1 if unit i's outcome is observed and 0 if it is missing (Li, Sheng & Yu 2025).",
    "Entries are OLS (linear probability) coefficients with HC1 standard errors in parentheses;",
    "a positive coefficient means a higher probability of being observed, i.e. less attrition."
  )
  star_note <- if (stars) " * p < 0.10, ** p < 0.05, *** p < 0.01." else ""

  if (type == "group") {
    return(paste0(
      dv,
      " Columns 1-2 estimate M on X separately within the control (Z = 0) and treatment (Z = 1) groups.",
      " Difference (T - C) is the treatment-group coefficient minus the control-group coefficient;",
      " it equals the Z \u00d7 X coefficient (Z for the constant) in the pooled regression of M on X, Z, and Z \u00d7 X,",
      " and its p-value is from an HC1-robust t test.",
      " The slopes test asks whether the covariates predict M in the same way in both groups (H0: all Z \u00d7 X = 0);",
      " the Chow test asks whether the two groups share the whole regression, constant included (H0: Z = 0 and all Z \u00d7 X = 0).",
      " Both are HC1-robust Wald tests reported in F form.",
      .ar_sparse_text(x$sparse),
      star_note
    ))
  }
  paste0(
    dv,
    " In reg. (b), the Treatment coefficient is the treatment-minus-control difference in Pr(M = 1) for units in the",
    " covariates' reference categories (all X = 0), and each Treatment \u00d7 X coefficient is how much the association",
    " between that covariate and M differs between the treatment and control groups.",
    " Joint HC1-robust Wald tests (F form): all X coefficients in reg. (a) = 0, ", .ar_f_text(test_row("all_x"), with_p = TRUE),
    "; all Z \u00d7 X coefficients in reg. (b) = 0 (slopes test), ", .ar_f_text(test_row("slopes"), with_p = TRUE),
    "; Z = 0 and all Z \u00d7 X = 0 in reg. (b) (Chow test), ", .ar_f_text(test_row("chow"), with_p = TRUE), ".",
    .ar_sparse_text(x$sparse),
    star_note
  )
}


#----------- Step 7: Markdown and LaTeX text ----------------------------------#

# Pipe table. By default the first column is left-aligned and the others are
# right-aligned; `left` gives the columns to left-align instead.
.ar_pipe_table <- function(tab, left = 1) {
  align <- rep("--:", ncol(tab))
  align[left] <- ":--"
  header <- paste0("| ", paste(names(tab), collapse = " | "), " |")
  rule <- paste0("|", paste(align, collapse = "|"), "|")
  rows <- apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(header, rule, rows)
}

.ar_markdown <- function(tab, caption) {
  paste(c(.ar_pipe_table(tab), "", paste0("Note: ", caption)), collapse = "\n")
}

.ar_latex_escape <- function(s) {
  s <- gsub("([&%$#_{}])", "\\\\\\1", s)
  s <- gsub("\u00d7", "$\\\\times$", s)
  s <- gsub("\u00b2", "$^2$", s)
  s
}

.ar_latex <- function(tab, caption) {
  cells <- as.data.frame(lapply(tab, .ar_latex_escape), check.names = FALSE, stringsAsFactors = FALSE)
  header <- paste(.ar_latex_escape(names(tab)), collapse = " & ")
  rows <- apply(cells, 1, function(r) paste(r, collapse = " & "))
  paste(c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\begin{tabular}{l", strrep("c", ncol(tab) - 1), "}"),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline",
    paste0(rows, " \\\\"),
    "\\hline",
    "\\end{tabular}",
    paste0("\\par\\vspace{4pt}\\parbox{\\linewidth}{\\footnotesize Note: ", .ar_latex_escape(caption), "}"),
    "\\end{table}"
  ), collapse = "\n")
}

# Data sections of the Markdown file: numbers with 6 significant digits, NA
# written as "NA"; text columns left-aligned, number columns right-aligned.
.ar_data_table <- function(df) {
  is_text <- !sapply(df, is.numeric)
  for (v in names(df)[!is_text]) {
    df[[v]] <- ifelse(is.na(df[[v]]), "NA", trimws(formatC(df[[v]], digits = 6, format = "g")))
  }
  .ar_pipe_table(df, left = which(is_text))
}

#' Write the missingness diagnostics as a Markdown or LaTeX file
#'
#' Writes both diagnostic tables (by experimental group and pooled) with their
#' notes. The Markdown file also contains the definitions of all quantities and
#' the underlying numbers (coefficients, group effects, model statistics and
#' joint tests) in long format, so that it can be read by people and by AI
#' tools.
#'
#' @param x An object returned by [attrition_regressions()].
#' @param format `"markdown"` or `"latex"`.
#' @param file Path of the file to write. If `NULL`, nothing is written.
#' @param stars Add significance stars to the estimates in the display tables.
#' @param digits Number of decimals in the display tables.
#'
#' @return The text of the file (invisibly when `file` is given).
#' @export
export_attrition_tables <- function(x, format = c("markdown", "latex"), file = NULL,
                                    stars = FALSE, digits = 3) {
  format <- match.arg(format)

  if (format == "latex") {
    text <- paste(
      "% Missingness regression diagnostics (riattrition)",
      "% Table 1: regressions by experimental group",
      attrition_table(x, type = "group", format = "latex", stars = stars, digits = digits),
      "",
      "% Table 2: missingness regressions with and without interactions",
      attrition_table(x, type = "pooled", format = "latex", stars = stars, digits = digits),
      sep = "\n"
    )
  } else {
    text <- paste(c(
      "# Missingness regression diagnostics",
      "",
      "## Definitions",
      "",
      "- `M` = 1 if a unit's outcome is observed and 0 if it is missing (Li, Sheng & Yu 2025).",
      "- `Z` = 1 for the treatment group and 0 for the control group. `X` are pretreatment covariates.",
      "- Models: `control` = M on X in the control group; `treatment` = M on X in the treatment group;",
      "  `reg_a` = M on X, all units; `reg_b` = M on X, Z and Z x X, all units.",
      "- All models are OLS linear probability models. `std_error` is the HC1 (heteroskedasticity-robust)",
      "  standard error and `p_value` is from a two-sided t test.",
      "- In `group_effects`, `difference` = treatment coefficient - control coefficient. It equals the",
      "  Z x X coefficient of `reg_b` (Z for the constant), and its standard error comes from `reg_b`.",
      "- Joint tests are HC1-robust Wald tests in F form: `slopes` tests all Z x X = 0 and `chow` tests",
      "  Z = 0 and all Z x X = 0, both in `reg_b`; `all_x` tests all X = 0 in `reg_a`.",
      "",
      "## Table 1. Regressions by experimental group",
      "",
      attrition_table(x, type = "group", format = "markdown", stars = stars, digits = digits),
      "",
      "## Table 2. Missingness regressions with and without interactions",
      "",
      attrition_table(x, type = "pooled", format = "markdown", stars = stars, digits = digits),
      "",
      "## Data: coefficients",
      "",
      .ar_data_table(x$coefficients),
      "",
      "## Data: group effects",
      "",
      .ar_data_table(x$group_effects),
      "",
      "## Data: model statistics",
      "",
      .ar_data_table(x$model_statistics),
      "",
      "## Data: joint tests",
      "",
      .ar_data_table(x$joint_tests)
    ), collapse = "\n")
  }

  if (is.null(file)) {
    return(text)
  }
  writeLines(enc2utf8(text), file, useBytes = TRUE)
  invisible(text)
}


#----------- Step 8: print method ---------------------------------------------#

#' @rdname attrition_regressions
#' @param ... Not used.
#' @export
print.attrition_regressions <- function(x, ...) {
  n <- x$model_statistics$n
  names(n) <- x$model_statistics$model
  cat("Missingness regression diagnostics (M = 1 means the outcome is observed)\n")
  cat("Units used: ", .ar_count(n[["reg_b"]]), " (control ", .ar_count(n[["control"]]),
      ", treatment ", .ar_count(n[["treatment"]]), "); covariates: ", length(x$covariates), "\n\n", sep = "")

  cat("Joint tests (HC1-robust Wald tests, F form):\n")
  for (i in seq_len(nrow(x$joint_tests))) {
    test <- x$joint_tests[i, ]
    cat("  ", test$model, ": ", test$hypothesis, "  ", .ar_f_text(test, with_p = TRUE), "\n", sep = "")
  }
  if (nrow(x$sparse) > 0) {
    cat("\n", trimws(.ar_sparse_text(x$sparse)), "\n", sep = "")
  }

  cat("\nTables: attrition_table(x, type = \"group\") or attrition_table(x, type = \"pooled\").\n")
  cat("Files:  export_attrition_tables(x, format = \"markdown\" or \"latex\", file = ...).\n")
  invisible(x)
}
