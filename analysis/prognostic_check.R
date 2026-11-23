#----------- Are covariates linked to differential attrition also prognostic? -#
#
# Used by the application notes in analysis/. It is not part of the package
# yet and relies on the helpers in R/3_attrition_diagnostics.R, so source that
# file first.
#
# Differential attrition on a covariate can bias the comparison of observed
# outcomes only if that covariate also predicts the outcome. Three pieces:
#   1. prognostic regression: Y on X among control-group units with observed
#      outcomes, so it describes how X predicts Y*_i(0);
#   2. side by side, per covariate: the Z x X coefficient of the missingness
#      regression (differential attrition) and the prognostic coefficient;
#   3. prognostic score: predicted Y*_i(0) = X_i' b for every unit, compared
#      between treated and control units among all units (chance imbalance)
#      and among units with observed outcomes (composition after attrition).

# Y : outcome, NA when missing.  M : 1 = observed (default: !is.na(Y)).
prognostic_check <- function(Y, Z, X, M = as.integer(!is.na(Y)), labels = NULL) {
  X <- as.data.frame(X)
  covars <- names(X)
  labels <- .ar_labels(covars, labels)

  # Missingness regressions (for the Z x X coefficients and the slopes test).
  attrition <- attrition_regressions(M, Z, X, labels = labels)

  # Units with M, Z and every covariate recorded.
  dat <- data.frame(Y = Y, M = M, Z = Z, X, check.names = FALSE)
  dat <- dat[stats::complete.cases(dat[, c("M", "Z", covars)]), , drop = FALSE]
  rhs <- paste(sprintf("`%s`", covars), collapse = " + ")

  # 1. Prognostic regression among control units with observed outcomes.
  prognostic_sample <- dat[dat$Z == 0 & dat$M == 1 & !is.na(dat$Y), , drop = FALSE]
  prognostic <- .ar_fit(stats::as.formula(paste("Y ~", rhs)), prognostic_sample)

  # 3. Predicted Y*(0) for every unit.
  b <- prognostic$coef
  X_all <- stats::model.matrix(stats::as.formula(paste("~", rhs)), dat)
  dat$score <- as.numeric(X_all[, names(b), drop = FALSE] %*% b)

  # Treated-minus-control difference (OLS coefficient on Z, HC1 SE).
  compare <- function(formula, d) {
    fit <- .ar_fit(formula, d)
    list(est = unname(fit$coef["Z"]), se = unname(fit$se["Z"]), p = unname(fit$p["Z"]), n = fit$n)
  }
  respondents <- dat[dat$M == 1 & !is.na(dat$Y), , drop = FALSE]

  structure(
    list(
      covariates = covars,
      labels = labels,
      attrition = attrition,
      prognostic = prognostic,
      prognostic_test = .ar_wald(prognostic, covars),
      sparse_prognostic = .ar_sparse(
        prognostic_sample[covars], labels,
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
    inter <- .ar_pick(x$attrition$coefficients, "reg_b", paste0("Z:", x$covariates))
    reg_b <- x$attrition$model_statistics[x$attrition$model_statistics$model == "reg_b", ]
    slopes <- x$attrition$joint_tests[x$attrition$joint_tests$test == "slopes", ]
    pr <- x$prognostic
    v <- x$covariates

    tab <- data.frame(
      Covariate = c(unname(x$labels), "N", .ar_adj_r2_label, "Joint test"),
      attr = c(.ar_est_se(inter$estimate, inter$std_error, inter$p_value, stars, digits),
               .ar_count(reg_b$n), .ar_num(reg_b$adj_r_squared, digits), .ar_f_text(slopes)),
      attr_p = c(.ar_num(inter$p_value, digits), "", "", .ar_num(slopes$p_value, digits)),
      prog = c(.ar_est_se(pr$coef[v], pr$se[v], pr$p[v], stars, digits),
               .ar_count(pr$n), .ar_num(pr$adj_r2, digits), .ar_f_text(x$prognostic_test)),
      prog_p = c(.ar_num(pr$p[v], digits), "", "", .ar_num(x$prognostic_test$p, digits)),
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
