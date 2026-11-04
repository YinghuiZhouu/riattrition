#----------- Missingness regression diagnostics for the application notes -----#
#
# Shared by the application notes in analysis/ (Bowles et al. 2025, Chen and
# Yang 2019). Notation follows Li, Sheng & Yu (2025, Section 2.1):
#   * M_i = 1 if unit i's outcome is observed, M_i = 0 if it is missing;
#   * Z_i = 1 if unit i is assigned to treatment, Z_i = 0 for control;
#   * X_i are pretreatment covariates.
#
# The diagnostics are linear probability models of M on X:
#   * by experimental group: M ~ X within control and within treatment, plus the
#     treatment-minus-control difference in each coefficient;
#   * pooled: reg. (a) M ~ X and reg. (b) M ~ X + Z + Z x X.
# The group differences equal the Z and Z x X coefficients of reg. (b), so they
# are taken from reg. (b) together with their HC1 standard errors.
#
# Only base R is used (lm() plus HC1 standard errors computed here), so the file
# can later move into R/ unchanged if the team decides to add it to the package.

# OLS fit with HC1 (heteroskedasticity-robust, n / (n - k)) standard errors.
.ar_fit <- function(formula, data) {
  fit <- stats::lm(formula, data = data)
  b <- stats::coef(fit)
  keep <- !is.na(b)
  X <- stats::model.matrix(fit)[, keep, drop = FALSE]
  e <- stats::residuals(fit)
  n <- nrow(X)
  k <- ncol(X)
  bread <- solve(crossprod(X))
  V <- n / (n - k) * bread %*% crossprod(X * e) %*% bread
  b <- b[keep]
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

# HC1-robust Wald test of H0: all listed coefficients = 0, in chi-square and F
# form (F = W / q with q and n - k degrees of freedom).
.ar_wald <- function(fit, terms) {
  terms <- intersect(terms, names(fit$coef))
  b <- fit$coef[terms]
  V <- fit$vcov[terms, terms, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(V, b))
  q <- length(terms)
  list(
    terms = terms,
    chi2 = W,
    F = W / q,
    df1 = q,
    df2 = fit$df_resid,
    p_chi2 = stats::pchisq(W, df = q, lower.tail = FALSE),
    p = stats::pf(W / q, q, fit$df_resid, lower.tail = FALSE)
  )
}

# Display labels for covariates; names without a label keep their own name.
.ar_labels <- function(covars, labels = NULL) {
  out <- stats::setNames(covars, covars)
  if (!is.null(labels)) {
    hit <- intersect(covars, names(labels))
    out[hit] <- labels[hit]
  }
  out
}

# Main entry point.
#   M      : 0/1 vector, 1 = outcome observed (Li, Sheng & Yu 2025).
#   Z      : 0/1 treatment vector.
#   X      : data frame of numeric pretreatment covariates.
#   labels : optional named character vector, covariate name -> display label.
attrition_regressions <- function(M, Z, X, labels = NULL) {
  X <- as.data.frame(X)
  covars <- names(X)
  if (!all(M %in% c(0, 1, NA)) || !all(Z %in% c(0, 1, NA))) {
    stop("`M` and `Z` must be 0/1 vectors (M = 1 means the outcome is observed).")
  }
  labels <- .ar_labels(covars, labels)

  dat <- stats::na.omit(data.frame(M = M, Z = Z, X, check.names = FALSE))
  rhs <- paste(sprintf("`%s`", covars), collapse = " + ")
  f_x <- stats::as.formula(paste("M ~", rhs))
  f_b <- stats::as.formula(paste0("M ~ Z * (", rhs, ")"))

  fits <- list(
    control = .ar_fit(f_x, dat[dat$Z == 0, , drop = FALSE]),
    treatment = .ar_fit(f_x, dat[dat$Z == 1, , drop = FALSE]),
    reg_a = .ar_fit(f_x, dat),
    reg_b = .ar_fit(f_b, dat)
  )

  slopes <- paste0("Z:", covars)
  structure(
    list(
      covariates = covars,
      labels = labels,
      fits = fits,
      tests = list(
        all_x = .ar_wald(fits$reg_a, covars),
        slopes = .ar_wald(fits$reg_b, slopes),
        chow = .ar_wald(fits$reg_b, c("Z", slopes))
      ),
      sparse = rbind(
        .ar_sparse(dat[dat$Z == 0, covars, drop = FALSE], labels, "control group"),
        .ar_sparse(dat[dat$Z == 1, covars, drop = FALSE], labels, "treatment group")
      )
    ),
    class = "attrition_regressions"
  )
}

#----------- Sparse 0/1 covariates --------------------------------------------#
#
# A 0/1 covariate that equals 1 (or 0) for only a handful of units in a
# regression sample is fit almost exactly for those units: their residuals are
# near zero, so the HC1 standard error of that coefficient, and any joint test
# that includes it, is unreliable. Each regression sample is checked and any
# such covariate is reported in the table note (and as an R warning).
.ar_min_cell <- 10

.ar_sparse <- function(X_sample, labels, sample_name, min_n = .ar_min_cell) {
  rows <- lapply(names(X_sample), function(v) {
    x <- X_sample[[v]]
    if (!all(x %in% c(0, 1))) {
      return(NULL)
    }
    n1 <- sum(x == 1)
    n0 <- sum(x == 0)
    if (min(n1, n0) >= min_n) {
      return(NULL)
    }
    data.frame(sample = sample_name, covariate = v, label = unname(labels[v]),
               n_ones = n1, n_zeros = n0, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(sample = character(0), covariate = character(0), label = character(0),
                      n_ones = integer(0), n_zeros = integer(0), stringsAsFactors = FALSE)
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

# Structured coefficient tables (numbers, not strings).
attrition_group_coefs <- function(x) {
  terms <- c(x$covariates, "(Intercept)")
  diff_terms <- c(paste0("Z:", x$covariates), "Z")
  get <- function(fit, what, tt) unname(fit[[what]][tt])
  data.frame(
    term = terms,
    label = c(unname(x$labels), "Constant"),
    control_est = get(x$fits$control, "coef", terms),
    control_se = get(x$fits$control, "se", terms),
    treatment_est = get(x$fits$treatment, "coef", terms),
    treatment_se = get(x$fits$treatment, "se", terms),
    diff_est = get(x$fits$reg_b, "coef", diff_terms),
    diff_se = get(x$fits$reg_b, "se", diff_terms),
    diff_p = get(x$fits$reg_b, "p", diff_terms),
    stringsAsFactors = FALSE
  )
}

attrition_pooled_coefs <- function(x) {
  terms <- c(x$covariates, "Z", paste0("Z:", x$covariates), "(Intercept)")
  labels <- c(
    unname(x$labels), "Treatment",
    paste("Treatment \u00d7", unname(x$labels)), "Constant"
  )
  get <- function(fit, what) unname(fit[[what]][terms])
  data.frame(
    term = terms,
    label = labels,
    a_est = get(x$fits$reg_a, "coef"),
    a_se = get(x$fits$reg_a, "se"),
    a_p = get(x$fits$reg_a, "p"),
    b_est = get(x$fits$reg_b, "coef"),
    b_se = get(x$fits$reg_b, "se"),
    b_p = get(x$fits$reg_b, "p"),
    stringsAsFactors = FALSE
  )
}

# Formatting helpers. Non-ASCII symbols are written as \u escapes so the file
# sources the same way in any locale.
.ar_adj_r2_label <- "Adjusted R\u00b2"

.ar_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

.ar_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", ""))))
}

.ar_est_se <- function(est, se, p, stars = FALSE, digits = 3) {
  out <- paste0(.ar_num(est, digits), if (stars) .ar_stars(p) else "", " (", .ar_num(se, digits), ")")
  out[is.na(est)] <- ""
  out
}

.ar_test_text <- function(test) {
  sprintf("F(%d, %s) = %.2f, p = %.3f", test$df1, format(test$df2, big.mark = ","), test$F, test$p)
}

.ar_count <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

# Display tables. type = "group" (by experimental group, with the T - C
# difference column and the Chow test) or "pooled" (reg. (a) and reg. (b)).
# format = "data.frame" (for knitr::kable), "markdown" or "latex".
attrition_table <- function(x, type = c("group", "pooled"),
                            format = c("data.frame", "markdown", "latex"),
                            stars = FALSE, digits = 3, caption = NULL) {
  type <- match.arg(type)
  format <- match.arg(format)

  if (type == "group") {
    co <- attrition_group_coefs(x)
    tab <- data.frame(
      Covariate = co$label,
      `Control group` = .ar_est_se(co$control_est, co$control_se, x$fits$control$p[co$term], stars, digits),
      `Treatment group` = .ar_est_se(co$treatment_est, co$treatment_se, x$fits$treatment$p[co$term], stars, digits),
      `Difference (T - C)` = .ar_est_se(co$diff_est, co$diff_se, co$diff_p, stars, digits),
      `p-value` = .ar_num(co$diff_p, digits),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    tab <- rbind(
      tab,
      data.frame(
        Covariate = c("N", .ar_adj_r2_label, "Slopes test: all Z \u00d7 X = 0", "Chow test: Z = 0 and all Z \u00d7 X = 0"),
        `Control group` = c(.ar_count(x$fits$control$n), .ar_num(x$fits$control$adj_r2, digits), "", ""),
        `Treatment group` = c(.ar_count(x$fits$treatment$n), .ar_num(x$fits$treatment$adj_r2, digits), "", ""),
        `Difference (T - C)` = c("", "", sprintf("F(%d, %s) = %.2f", x$tests$slopes$df1, .ar_count(x$tests$slopes$df2), x$tests$slopes$F),
                                 sprintf("F(%d, %s) = %.2f", x$tests$chow$df1, .ar_count(x$tests$chow$df2), x$tests$chow$F)),
        `p-value` = c("", "", .ar_num(x$tests$slopes$p, digits), .ar_num(x$tests$chow$p, digits)),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    )
  } else {
    co <- attrition_pooled_coefs(x)
    # The last two rows give the p-values of the joint tests on reg. (b):
    # the slopes test (all Z x X = 0) and the Chow test (Z and all Z x X = 0).
    tab <- data.frame(
      Covariate = c(
        co$label, "N", .ar_adj_r2_label,
        "p-value: all Z \u00d7 X = 0", "p-value: Z and all Z \u00d7 X = 0"
      ),
      reg_a = c(
        .ar_est_se(co$a_est, co$a_se, co$a_p, stars, digits),
        .ar_count(x$fits$reg_a$n), .ar_num(x$fits$reg_a$adj_r2, digits),
        "", ""
      ),
      reg_b = c(
        .ar_est_se(co$b_est, co$b_se, co$b_p, stars, digits),
        .ar_count(x$fits$reg_b$n), .ar_num(x$fits$reg_b$adj_r2, digits),
        .ar_num(x$tests$slopes$p, digits), .ar_num(x$tests$chow$p, digits)
      ),
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

# Table note: what M is, how to read the coefficients, and the joint tests.
attrition_table_note <- function(x, type = c("group", "pooled"), stars = FALSE) {
  type <- match.arg(type)
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
    " Joint HC1-robust Wald tests (F form): all X coefficients in reg. (a) = 0, ", .ar_test_text(x$tests$all_x),
    "; all Z \u00d7 X coefficients in reg. (b) = 0 (slopes test), ", .ar_test_text(x$tests$slopes),
    "; Z = 0 and all Z \u00d7 X = 0 in reg. (b) (Chow test), ", .ar_test_text(x$tests$chow), ".",
    .ar_sparse_text(x$sparse),
    star_note
  )
}

#----------- Are covariates linked to differential attrition also prognostic? -#
#
# Differential attrition on a covariate can bias the observed-outcome
# comparison only if that covariate is also prognostic of the outcome. Three
# pieces:
#   1. prognostic regression: Y on X among control-group units with observed
#      outcomes, so it describes how X predicts Y*_i(0);
#   2. side by side, per covariate: the Z x X coefficient of the missingness
#      regression (differential attrition) and the prognostic coefficient;
#   3. prognostic score: predicted Y*_i(0) = X_i' beta-hat for every unit, compared
#      between treated and control among all units (chance imbalance) and among
#      units with observed outcomes (composition after attrition).
#
#   Y : outcome, NA when missing.  M : 1 = observed (default: !is.na(Y)).
prognostic_check <- function(Y, Z, X, M = as.integer(!is.na(Y)), labels = NULL) {
  X <- as.data.frame(X)
  covars <- names(X)
  labels <- .ar_labels(covars, labels)
  ar <- attrition_regressions(M, Z, X, labels = labels)

  dat <- data.frame(Y = Y, M = M, Z = Z, X, check.names = FALSE)
  dat <- dat[stats::complete.cases(dat[, c("M", "Z", covars)]), , drop = FALSE]
  rhs <- paste(sprintf("`%s`", covars), collapse = " + ")

  prognostic_sample <- dat[dat$Z == 0 & dat$M == 1 & !is.na(dat$Y), , drop = FALSE]
  prognostic <- .ar_fit(stats::as.formula(paste("Y ~", rhs)), prognostic_sample)
  b <- prognostic$coef
  X_all <- stats::model.matrix(stats::as.formula(paste("~", rhs)), dat)
  dat$score <- as.numeric(X_all[, names(b), drop = FALSE] %*% b)

  compare <- function(formula, d) {
    fit <- .ar_fit(formula, d)
    list(est = unname(fit$coef["Z"]), se = unname(fit$se["Z"]), p = unname(fit$p["Z"]), n = fit$n)
  }
  respondents <- dat[dat$M == 1 & !is.na(dat$Y), , drop = FALSE]

  structure(
    list(
      covariates = covars,
      labels = labels,
      attrition = ar,
      prognostic = prognostic,
      prognostic_test = .ar_wald(prognostic, covars),
      sparse_prognostic = .ar_sparse(
        prognostic_sample[, covars, drop = FALSE], labels,
        "prognostic sample (control units with observed outcomes)"
      ),
      score = list(
        all = compare(score ~ Z, dat),
        respondents = compare(score ~ Z, respondents),
        observed_y = compare(Y ~ Z, respondents)
      )
    ),
    class = "prognostic_check"
  )
}

# type = "side_by_side" (steps 1-2) or "score" (step 3).
prognostic_table <- function(x, type = c("side_by_side", "score"),
                             format = c("data.frame", "markdown", "latex"),
                             stars = FALSE, digits = 3, caption = NULL) {
  type <- match.arg(type)
  format <- match.arg(format)

  if (type == "side_by_side") {
    rb <- x$attrition$fits$reg_b
    pr <- x$prognostic
    inter <- paste0("Z:", x$covariates)
    tab <- data.frame(
      Covariate = c(unname(x$labels), "N", .ar_adj_r2_label, "Joint test"),
      attr = c(
        .ar_est_se(rb$coef[inter], rb$se[inter], rb$p[inter], stars, digits),
        .ar_count(rb$n), .ar_num(rb$adj_r2, digits),
        sprintf("F(%d, %s) = %.2f", x$attrition$tests$slopes$df1, .ar_count(x$attrition$tests$slopes$df2), x$attrition$tests$slopes$F)
      ),
      attr_p = c(.ar_num(rb$p[inter], digits), "", "", .ar_num(x$attrition$tests$slopes$p, digits)),
      prog = c(
        .ar_est_se(pr$coef[x$covariates], pr$se[x$covariates], pr$p[x$covariates], stars, digits),
        .ar_count(pr$n), .ar_num(pr$adj_r2, digits),
        sprintf("F(%d, %s) = %.2f", x$prognostic_test$df1, .ar_count(x$prognostic_test$df2), x$prognostic_test$F)
      ),
      prog_p = c(.ar_num(pr$p[x$covariates], digits), "", "", .ar_num(x$prognostic_test$p, digits)),
      stringsAsFactors = FALSE
    )
    names(tab) <- c(
      "Covariate", "Differential attrition: Z \u00d7 X", "p-value",
      "Prognostic for Y: coefficient", "p-value"
    )
  } else {
    s <- x$score
    row <- function(r) c(.ar_est_se(r$est, r$se, r$p, stars, digits), .ar_num(r$p, digits), .ar_count(r$n))
    tab <- as.data.frame(rbind(
      c("Predicted Y*(0), all randomized units", row(s$all)),
      c("Predicted Y*(0), units with observed Y", row(s$respondents)),
      c("Observed Y, units with observed Y", row(s$observed_y))
    ), stringsAsFactors = FALSE)
    names(tab) <- c("Comparison", "Treated - control (SE)", "p-value", "N")
  }
  rownames(tab) <- NULL

  if (format == "data.frame") {
    return(tab)
  }
  if (is.null(caption)) {
    caption <- prognostic_table_note(x, type = type, stars = stars)
  }
  if (format == "markdown") {
    return(.ar_markdown(tab, caption))
  }
  .ar_latex(tab, caption)
}

prognostic_table_note <- function(x, type = c("side_by_side", "score"), stars = FALSE) {
  type <- match.arg(type)
  star_note <- if (stars) " * p < 0.10, ** p < 0.05, *** p < 0.01." else ""
  if (type == "side_by_side") {
    return(paste0(
      "Differential attrition: the Z \u00d7 X coefficients of the missingness regression of M on X, Z, and Z \u00d7 X",
      " (M_i = 1 if unit i's outcome is observed), i.e. how differently each covariate predicts being observed",
      " in the treatment and control groups; the joint test is the slopes test (all Z \u00d7 X = 0).",
      " Prognostic for Y: OLS regression of Y on X among control-group units with observed outcomes, so it describes",
      " how the covariates predict Y*_i(0); the joint test is all X = 0.",
      " HC1 standard errors in parentheses; joint tests are HC1-robust Wald tests in F form.",
      " A covariate can bias the comparison of observed outcomes only if it is both associated with differential",
      " attrition and prognostic of Y.",
      .ar_sparse_text(rbind(x$attrition$sparse, x$sparse_prognostic)),
      star_note
    ))
  }
  paste0(
    "Predicted Y*_i(0) = X_i' b, with b from the prognostic regression (Y on X among control-group units with",
    " observed outcomes), computed for every unit, including units whose outcome is missing.",
    " Among all randomized units, a treated-minus-control difference in predicted Y*(0) reflects chance imbalance only.",
    " Among units with observed outcomes, it shows how differently the covariate mix of respondents is related to Y",
    " in the two groups: it is the part of the observed-outcome difference (last row) that differential attrition on",
    " these covariates would produce even with no treatment effect.",
    " Differences are OLS coefficients on Z with HC1 standard errors; they treat b as fixed.",
    " Unobserved factors can still drive differential attrition.",
    .ar_sparse_text(x$sparse_prognostic),
    star_note
  )
}

.ar_markdown <- function(tab, caption) {
  header <- paste0("| ", paste(names(tab), collapse = " | "), " |")
  rule <- paste0("|", paste(c(":--", rep("--:", ncol(tab) - 1)), collapse = "|"), "|")
  rows <- apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, rule, rows, "", paste0("Note: ", caption)), collapse = "\n")
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
