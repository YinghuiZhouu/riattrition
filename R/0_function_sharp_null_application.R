#------------------ Functions for Empirical Applications ----------------------#

#' Wrapper function for empirical application
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param seed Random seed.
#'
#' @return A list of randomization inference results.
#' @export
application_wrapper <- function(Z, Y, seed = 1) {
  set.seed(seed)
  # reshuffle the data
  idx <- sample(seq_along(Z))
  Z <- Z[idx]
  Y <- Y[idx]

  # Set up parameters
  alpha <- 0.1
  nperm <- 10^5
  tol <- 10^(-3)

  # One-step: precompute
  Z.perm.g.mp.mn <- assign_CRE(n = length(Z), m = sum(Z), nperm = nperm)
  Z.perm.s.r <- assign_CRE(n = length(Z[!is.na(Y)]), m = sum(Z[!is.na(Y)]), nperm = nperm)
  stat.null.g.mp.mn <- null_dist(n = length(Z), m = sum(Z), class = "RS", method.list = list(name = "Wilcoxon"), Z.perm = Z.perm.g.mp.mn, nperm = nperm)
  stat.null.s.r <- null_dist(n = length(Z[!is.na(Y)]), m = sum(Z[!is.na(Y)]), class = "RS", method.list = list(name = "Wilcoxon"), Z.perm = Z.perm.s.r, nperm = nperm)

  # Two-step: precompute
  stat.null.g.mp.twostep <- null_dist(n = length(Z), m = sum(Z), class = "MWU+", method.list = list(name = "Wilcoxon"), Z.perm = Z.perm.g.mp.mn, nperm = nperm)
  stat.null.mn.twostep <- null_dist(n = length(Z), m = sum(Z), class = "MWU-", method.list = list(name = "Wilcoxon"), Z.perm = Z.perm.g.mp.mn, nperm = nperm)

  # switch treatment label and outcome sign: prepare for sensitivity analysis
  switched <- mean(is.na(Y[Z == 1])) > mean(is.na(Y[Z == 0]))
  if (switched) {
    mx <- "mn"
    Z.switch <- 1 - Z
    Y.switch <- - Y
    Z.perm.switch <- assign_CRE(n = length(Z.switch), m = sum(Z.switch), nperm = nperm)
    stat.null.switch <- null_dist(n = length(Z.switch), m = sum(Z.switch), class = "MWU+", method.list = list(name = "Wilcoxon"), Z.perm = Z.perm.switch, nperm = nperm)
    # Z.switch, Y.switch, Z.perm.switch, stat.null.switch are used only for sensitivity analysis
  } else {
    mx <- "mp"
  }

  # Sharp null: p-value (one-step)
  pval.g <- pval_sharp(Z = Z, Y = Y, c = 0, missing = "general", class = "RS", method.list = list(name = "Wilcoxon"), stat.null = stat.null.g.mp.mn, Z.perm = Z.perm.g.mp.mn, nperm = nperm)
  pval.mx <- pval_sharp(Z = Z, Y = Y, c = 0, missing = mx, class = "RS", method.list = list(name = "Wilcoxon"), stat.null = stat.null.g.mp.mn, Z.perm = Z.perm.g.mp.mn, nperm = nperm)
  pval.s  <- pval_sharp(Z = Z, Y = Y, c = 0, missing = "sharp", class = "RS", method.list = list(name = "Wilcoxon"), stat.null = stat.null.s.r, Z.perm = Z.perm.s.r, nperm = nperm)
  pval.r  <- pval_sharp(Z = Z, Y = Y, c = 0, missing = "random", class = "RS", method.list = list(name = "Wilcoxon"), stat.null = stat.null.s.r, Z.perm = Z.perm.s.r, nperm = nperm)

  # Sharp null: CI (one-step)
  ci.sharp.g  <- ci_sharp(Z = Z, Y = Y, alternative = "two.sided", missing = "general", class = "RS", method.list = list(name = "Wilcoxon"), stat.null = stat.null.g.mp.mn, Z.perm = Z.perm.g.mp.mn, nperm = nperm, alpha = alpha, tol = tol)
  ci.sharp.mx <- ci_sharp(Z = Z, Y = Y, alternative = "two.sided", missing = mx, class = "RS", method.list = list(name = "Wilcoxon"), stat.null = stat.null.g.mp.mn, Z.perm = Z.perm.g.mp.mn, nperm = nperm, alpha = alpha, tol = tol)
  ci.sharp.s  <- ci_sharp(Z = Z, Y = Y, alternative = "two.sided", missing = "sharp", class = "RS", method.list = list(name = "Wilcoxon"), stat.null = stat.null.s.r, Z.perm = Z.perm.s.r, nperm = nperm, alpha = alpha, tol = tol)
  ci.sharp.r  <- ci_sharp(Z = Z, Y = Y, alternative = "two.sided", missing = "random", class = "RS", method.list = list(name = "Wilcoxon"), stat.null = stat.null.s.r, Z.perm = Z.perm.s.r, nperm = nperm, alpha = alpha, tol = tol)

  # Sharp null: p-value (two-step)
  pval.g.twostep <- pval_sharp_twostep(Z = Z, Y = Y, c = 0, missing = "general", method.list = list(name = "Wilcoxon"), stat.null = stat.null.g.mp.twostep, Z.perm = Z.perm.g.mp.mn, nperm = nperm, beta = 0.1 * alpha)
  if (mx == "mp") {
    pval.mx.twostep <- pval_sharp_twostep(Z = Z, Y = Y, c = 0, missing = mx, method.list = list(name = "Wilcoxon"), stat.null = stat.null.g.mp.twostep, Z.perm = Z.perm.g.mp.mn, nperm = nperm, beta = 0.1 * alpha)
  } else {
    pval.mx.twostep <- pval_sharp_twostep(Z = Z, Y = Y, c = 0, missing = mx, method.list = list(name = "Wilcoxon"), stat.null = stat.null.mn.twostep, Z.perm = Z.perm.g.mp.mn, nperm = nperm, beta = 0.1 * alpha)
  }

  # Sharp null: CI (two-step)
  ci.sharp.g.twostep  <- ci_sharp_twostep(Z = Z, Y = Y, alternative = "two.sided", missing = "general", method.list = list(name = "Wilcoxon"), stat.null = stat.null.g.mp.twostep, Z.perm = Z.perm.g.mp.mn, nperm = nperm, alpha = alpha, beta = 0.1 * alpha, tol = tol)
  if (mx == "mp") {
    ci.sharp.mx.twostep  <- ci_sharp_twostep(Z = Z, Y = Y, alternative = "two.sided", missing = mx, method.list = list(name = "Wilcoxon"), stat.null = stat.null.g.mp.twostep, Z.perm = Z.perm.g.mp.mn, nperm = nperm, alpha = alpha, beta = 0.1 * alpha, tol = tol)
  } else {
    ci.sharp.mx.twostep  <- ci_sharp_twostep(Z = Z, Y = Y, alternative = "two.sided", missing = mx, method.list = list(name = "Wilcoxon"), stat.null = stat.null.mn.twostep, Z.perm = Z.perm.g.mp.mn, nperm = nperm, alpha = alpha, beta = 0.1 * alpha, tol = tol)
  }

  # Sensitivity analysis
  if (mx == "mp") {
    sens_pval <- sens_pval_sharp_mp_twostep(
      Z = Z, Y = Y, c = 0, method.list = list(name = "Wilcoxon"),
      stat.null = stat.null.g.mp.twostep, Z.perm = Z.perm.g.mp.mn, nperm = nperm,
      beta = 0.1 * alpha
    )
  } else {
    sens_pval <- sens_pval_sharp_mp_twostep(
      Z = Z.switch, Y = Y.switch, c = 0, method.list = list(name = "Wilcoxon"),
      stat.null = stat.null.switch, Z.perm = Z.perm.switch, nperm = nperm,
      beta = 0.1 * alpha
    )
  }

  return(list(
    Z = Z,
    Y = Y,
    mx = mx,
    pval.g = pval.g,
    pval.mx = pval.mx,
    pval.s = pval.s,
    pval.r = pval.r,
    ci.sharp.g = ci.sharp.g,
    ci.sharp.mx = ci.sharp.mx,
    ci.sharp.s = ci.sharp.s,
    ci.sharp.r = ci.sharp.r,
    pval.g.twostep = pval.g.twostep,
    pval.mx.twostep = pval.mx.twostep,
    ci.sharp.g.twostep = ci.sharp.g.twostep,
    ci.sharp.mx.twostep = ci.sharp.mx.twostep,
    sens_pval = sens_pval
  ))
}

#' Generate LaTex table
#'
#' @param result.list A list of outputs from `application_wrapper()`, one per table row.
#'
#' @return A list of two LaTex tables.
#' @export
make_latex_table <- function(result.list) {
  mx <- result.list[[1]]$mx
  if (any(vapply(result.list, function(result) result$mx, character(1)) != mx)) {
    stop("All rows in one table must share the same `mx`.")
  }
  if (mx == "mp") {
    coltitle <- textab::TexRow(c("Missing Proportion Total", "Missing Proportion Treated", "Missing Proportion Control", "General", "MP", "Sharp", "Random", "Two-Step General", "Two-Step MP"))
  } else {
    coltitle <- textab::TexRow(c("Missing Proportion Total", "Missing Proportion Treated", "Missing Proportion Control", "General", "MN", "Sharp", "Random", "Two-Step General", "Two-Step MN"))
  }
  table_pval <- coltitle
  table_ci <- coltitle
  for (result in result.list) {
    row_pval <- textab::TexRow(c(mean(is.na(result$Y)), mean(is.na(result$Y)[result$Z == 1]), mean(is.na(result$Y)[result$Z == 0]), result$pval.g, result$pval.mx, result$pval.s, result$pval.r, result$pval.g.twostep, result$pval.mx.twostep))
    row_ci <- textab::TexRow(c(mean(is.na(result$Y)), mean(is.na(result$Y)[result$Z == 1]), mean(is.na(result$Y)[result$Z == 0])), dec = 3) / textab::TexRow(c(paste0("[", paste(sprintf("%.3f", result$ci.sharp.g), collapse = ", "), "]"), paste0("[", paste(sprintf("%.3f", result$ci.sharp.mx), collapse = ", "), "]"), paste0("[", paste(sprintf("%.3f", result$ci.sharp.s), collapse = ", "), "]"), paste0("[", paste(sprintf("%.3f", result$ci.sharp.r), collapse = ", "), "]"), paste0("[", paste(sprintf("%.3f", result$ci.sharp.g.twostep), collapse = ", "), "]"), paste0("[", paste(sprintf("%.3f", result$ci.sharp.mx.twostep), collapse = ", "), "]")))
    table_pval <- table_pval + row_pval
    table_ci <- table_ci + row_ci
  }
  return(list(table_pval = table_pval, table_ci = table_ci))
}

#' Draw sensitivity analysis plot
#'
#' @param sens_pval A dataframe of sensitivity analysis.
#' @param output_file A string of output file path for the plot.
#'
#' @return A ggplot object.
#' @export
plot_sens <- function(sens_pval, output_file) {
  max_pval <- max(sens_pval$pval, na.rm = TRUE)
  x_position <- 0.8 * max(sens_pval$m10)
  label_size <- 10
  alpha <- 0.1 # nominal level

  p <- ggplot(sens_pval, aes(x = .data$m10, y = .data$pval)) +
    geom_line(linewidth = 0.8) +
    # nominal level
    geom_hline(yintercept = alpha, linetype = "dashed") +
    annotate("text", x = x_position, y = alpha, vjust = 1.5,
             label = paste("nominal level =", alpha),
             size = label_size, family = "CMU Serif") +
    # maximum p-value over m10
    geom_hline(yintercept = max_pval, linetype = "dashed") +
    annotate("text", x = x_position, y = max_pval, vjust = 1.5,
             label = sprintf("maximum p-value = %.3f", max_pval),
             size = label_size, family = "CMU Serif") +
    labs(x = expression(m[10]), y = "p-value") +
    theme_minimal() +
    theme(
      text = element_text(family = "CMU Serif"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey", fill = NA, linewidth = 1)
    )

  ggsave(filename = output_file, plot = p, width = 30, height = 20, units = "cm", device = CairoPDF)
  return(p)
}
