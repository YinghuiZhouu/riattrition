#----------- Visualization helper for null distributions ---------------------#
#
# This file adds a lightweight plotting layer for applied users.  It does not
# change any inferential logic.  Instead, it reuses the same ingredients as the
# printed wrapper:
#   * .ri_prepare_analysis_data() from R/1_result_output.R
#   * .ri_default_assignments()   from R/1_result_output.R
#   * .ri_sharp_details()         from R/1_result_output.R
#   * null_dist() / test_stat()   from R/0_function_sharp_null.R
#
# The function supports two visualization modes:
#   * analysis = "general": use the package's general-missingness procedure
#   * analysis = "complete-case": restrict to observed outcomes only
#
# This keeps the plotting layer close to applied practice while still reusing
# the underlying randomization-inference engine.

#' Visualize the null distribution under general missingness
#'
#' `plot_ri_test_stat()` draws the randomization null distribution, marks the
#' observed test statistic, and highlights the right tail that determines the
#' sharp-null p-value. It supports the package's `general` missingness
#' procedure and an observed-only `complete-case` visualization mode.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param analysis A string that specifies which inferential path to visualize:
#'   `"general"` uses the package's general-missingness procedure, while
#'   `"complete-case"` restricts the analysis to observed outcomes only.
#' @param missing A string that specifies the missing mechanism. The current
#'   visualization method uses `"general"` when `analysis = "general"`. It is
#'   ignored when `analysis = "complete-case"`.
#' @param class A string that specifies the class of test statistic for the
#'   sharp-null procedure.
#' @param method.list A list that specifies the choice of test statistic.
#' @param block An optional block identifier vector. If supplied, permutations
#'   are generated within blocks while preserving the treated count in each
#'   block.
#' @param cluster An optional cluster identifier vector. If supplied, the
#'   analysis treats clusters as the unit of analysis.
#' @param cluster_outcome A string that specifies how individual outcomes are
#'   aggregated to the cluster level. Currently only `"mean_observed"` is
#'   supported.
#' @param cluster_missing A string that specifies when the cluster-level outcome
#'   is missing.
#' @param nperm A positive integer that specifies the number of permutations.
#' @param style Plot style. `"histogram"` is the default because it directly
#'   represents the empirical permutation distribution. `"density"` is available
#'   as a smoother alternative.
#' @param bins Histogram break rule passed to [graphics::hist()] when
#'   `style = "histogram"`.
#' @param show_density A logical value indicating whether a density curve should
#'   be overlaid when `style = "histogram"`.
#' @param main Optional main title. If `NULL`, the function supplies a default.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param null_col Fill color for the main null distribution.
#' @param tail_col Fill color for the right-tail region used in the p-value.
#' @param density_col Line color for the optional density overlay.
#' @param stat_col Color of the observed-statistic reference line.
#' @param border Border color for histogram bars.
#' @param lwd Width of the observed-statistic reference line.
#' @param ... Additional graphical arguments passed to the underlying plotting
#'   function.
#'
#' @return Invisibly returns an object of class `riattrition_plot` containing
#'   the null distribution, observed statistic, p-value, and plotting metadata.
#' @export
plot_ri_test_stat <- function(Z, Y, c = 0, missing = "general", class = "RS",
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
                              main = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              null_col = "#D7DDE3",
                              tail_col = "#C9B8B4",
                              density_col = "#3E5F7E",
                              stat_col = "#8F3B3B",
                              border = "white",
                              lwd = 2,
                              ...) {
  analysis <- match.arg(analysis)
  style <- match.arg(style)

  if (identical(analysis, "general") && !identical(missing, "general")) {
    stop("When analysis = \"general\", plot_ri_test_stat() currently supports only missing = \"general\".")
  }

  viz <- .ri_visualization_data(
    Z = Z, Y = Y, c = c, analysis = analysis, missing = missing, class = class,
    method.list = method.list, block = block, cluster = cluster,
    cluster_outcome = cluster_outcome, cluster_missing = cluster_missing,
    nperm = nperm
  )

  .ri_draw_null_distribution(
    x = viz, style = style, bins = bins, show_density = show_density,
    main = main, xlab = xlab, ylab = ylab, null_col = null_col,
    tail_col = tail_col, density_col = density_col, stat_col = stat_col,
    border = border, lwd = lwd, ...
  )

  invisible(viz)
}

#----------- New helper: compute visualization data --------------------------#
#
# This helper computes the same sharp-null ingredients used in the printed
# wrapper and stores them in a small object that is easy to plot or inspect.
.ri_visualization_data <- function(Z, Y, c, analysis, missing, class, method.list,
                                   block, cluster, cluster_outcome,
                                   cluster_missing, nperm) {
  analysis_data <- .ri_prepare_analysis_data(
    Z = Z, Y = Y, block = block, cluster = cluster,
    cluster_outcome = cluster_outcome, cluster_missing = cluster_missing
  )

  Z.analysis <- analysis_data$Z
  Y.analysis <- analysis_data$Y
  block.analysis <- analysis_data$block
  counts <- .ri_missing_counts(Z.analysis, Y.analysis)
  test_name <- .ri_test_label(class = class, method.list = method.list)
  if (identical(analysis, "general")) {
    Z.perm <- .ri_default_assignments(
      Z = Z.analysis, Y = Y.analysis, missing = "general",
      block = block.analysis, nperm = nperm
    )
    stat.null <- null_dist(
      n = counts$n_total, m = counts$n_treat, class = class,
      method.list = method.list, Z.perm = Z.perm, nperm = nperm
    )
    sharp_details <- .ri_sharp_details(
      Z = Z.analysis, Y = Y.analysis, c = c, missing = "general", class = class,
      method.list = method.list, stat.null = stat.null, Z.perm = Z.perm,
      nperm = nperm
    )
    mode_label <- "general missingness"
  } else {
    Z.perm <- .ri_default_assignments(
      Z = Z.analysis, Y = Y.analysis, missing = "random",
      block = block.analysis, nperm = nperm
    )
    stat.null <- null_dist(
      n = counts$n_obs, m = counts$n_obs_treat, class = class,
      method.list = method.list, Z.perm = Z.perm, nperm = nperm
    )
    sharp_details <- .ri_sharp_details(
      Z = Z.analysis, Y = Y.analysis, c = c, missing = "random", class = class,
      method.list = method.list, stat.null = stat.null, Z.perm = Z.perm,
      nperm = nperm
    )
    mode_label <- "complete-case analysis"
  }

  structure(
    list(
      call = match.call(),
      analysis = analysis,
      assumption = if (identical(analysis, "general")) "general" else "complete-case",
      mode_label = mode_label,
      stat_null = stat.null,
      stat_obs = sharp_details$stat_obs,
      p_value = sharp_details$p_value,
      test_name = test_name,
      design = analysis_data$design,
      inference_unit = analysis_data$inference_unit,
      counts = counts,
      nperm = nperm
    ),
    class = "riattrition_plot"
  )
}

#----------- New helper: draw histogram or density view ----------------------#
#
# The plot highlights the right tail because the package's sharp-null p-values
# are computed as Pr(T_null >= T_obs).
.ri_draw_null_distribution <- function(x, style, bins, show_density, main,
                                       xlab, ylab, null_col, tail_col,
                                       density_col, stat_col, border, lwd,
                                       ...) {
  null_range <- range(x$stat_null)
  obs_position <- if (x$stat_obs < null_range[1]) {
    "left"
  } else if (x$stat_obs > null_range[2]) {
    "right"
  } else {
    "inside"
  }

  if (is.null(main)) {
    main <- paste("Null distribution under", x$mode_label)
  }
  if (is.null(xlab)) {
    xlab <- x$test_name
  }
  if (is.null(ylab)) {
    ylab <- if (identical(style, "histogram")) "Density" else "Density"
  }

  subtitle_line_1 <- paste0(
    "Obs = ", .ri_format_num(x$stat_obs),
    " | p = ", .ri_format_num(x$p_value)
  )
  subtitle_line_2 <- paste0(
    x$design,
    " | nperm = ", x$nperm
  )

  tail_note <- if (identical(obs_position, "left")) {
    "Observed statistic is left of null support; the entire histogram is in the right tail."
  } else if (identical(obs_position, "right")) {
    "Observed statistic is right of null support; no null draws fall in the right tail."
  } else {
    "Right-tail region is shaded; the p-value is the share of null draws at least as large as the observed statistic."
  }

  if (identical(style, "histogram")) {
    hist_obj <- graphics::hist(
      x$stat_null,
      breaks = bins,
      plot = FALSE
    )

    fill_cols <- rep(null_col, length(hist_obj$counts))
    fill_cols[hist_obj$mids >= x$stat_obs] <- tail_col

    graphics::plot(
      hist_obj,
      freq = FALSE,
      col = fill_cols,
      border = border,
      main = main,
      sub = "",
      xlab = xlab,
      ylab = ylab,
      ...
    )

    if (isTRUE(show_density)) {
      dens <- stats::density(x$stat_null)
      graphics::lines(dens, col = density_col, lwd = 2.5)
    }
  } else {
    dens <- stats::density(x$stat_null)
    y_max <- max(dens$y)

    graphics::plot(
      dens,
      type = "n",
      main = main,
      sub = "",
      xlab = xlab,
      ylab = ylab,
      ylim = c(0, y_max * 1.05),
      ...
    )

    graphics::polygon(
      x = c(dens$x, rev(dens$x)),
      y = c(dens$y, rep(0, length(dens$y))),
      col = null_col,
      border = NA
    )

    tail_idx <- dens$x >= x$stat_obs
    if (any(tail_idx)) {
      graphics::polygon(
        x = c(x$stat_obs, dens$x[tail_idx], max(dens$x[tail_idx])),
        y = c(0, dens$y[tail_idx], 0),
        col = tail_col,
        border = NA
      )
    }

    graphics::lines(dens, col = density_col, lwd = 2)
  }

  usr <- graphics::par("usr")
  y_top <- usr[4]

  if (identical(obs_position, "inside")) {
    graphics::abline(v = x$stat_obs, col = stat_col, lwd = lwd, lty = 2)
  } else if (identical(obs_position, "left")) {
    x_left <- usr[1]
    x_pad <- 0.04 * diff(usr[1:2])
    graphics::segments(x_left, 0, x_left, 0.93 * y_top, col = stat_col, lwd = lwd, lty = 2, xpd = NA)
    graphics::arrows(x_left, 0.93 * y_top, x_left - x_pad, 0.93 * y_top,
                     length = 0.08, code = 2, col = stat_col, lwd = lwd, xpd = NA)
  } else if (identical(obs_position, "right")) {
    x_right <- usr[2]
    x_pad <- 0.04 * diff(usr[1:2])
    graphics::segments(x_right, 0, x_right, 0.93 * y_top, col = stat_col, lwd = lwd, lty = 2, xpd = NA)
    graphics::arrows(x_right, 0.93 * y_top, x_right + x_pad, 0.93 * y_top,
                     length = 0.08, code = 2, col = stat_col, lwd = lwd, xpd = NA)
  }

  graphics::mtext(subtitle_line_1, side = 1, line = 3.2, cex = 0.9)
  graphics::mtext(subtitle_line_2, side = 1, line = 4.3, cex = 0.9)
  graphics::mtext(tail_note, side = 3, line = 0.2, cex = 0.78)
}
