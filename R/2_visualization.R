#----------- Visualization helper for null distributions ---------------------#
#
# This file adds a lightweight plotting layer for applied users.  
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
#' @param off_support A string that specifies how to display the observed test
#'   statistic when it falls outside the support of the permutation null
#'   distribution. `"auto"` extends the axis for modest gaps and uses a split
#'   axis for large gaps. `"arrow"` keeps the current support and points toward
#'   the off-support value. `"extend"` enlarges the x-axis to include the
#'   observed statistic. `"split"` uses a split-axis display so that the
#'   observed statistic remains visible without compressing the null
#'   distribution.
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
                              ...) {
  analysis <- match.arg(analysis)
  style <- match.arg(style)
  off_support <- match.arg(off_support)

  if (identical(analysis, "general") && !identical(missing, "general")) {
    stop("When analysis = \"general\", plot_ri_test_stat() currently supports only missing = \"general\".")
  }

  viz <- .ri_visualization_data(
    Z = Z, Y = Y, c = c, analysis = analysis, missing = missing, class = class,
    method.list = method.list, block = block, cluster = cluster,
    cluster_outcome = cluster_outcome, cluster_missing = cluster_missing,
    nperm = nperm
  )

  plot_meta <- .ri_draw_null_distribution(
    x = viz, style = style, bins = bins, show_density = show_density,
    off_support = off_support,
    main = main, xlab = xlab, ylab = ylab, null_col = null_col,
    tail_col = tail_col, density_col = density_col, stat_col = stat_col,
    border = border, lwd = lwd, ...
  )

  viz$off_support_mode <- plot_meta$off_support_mode
  viz$obs_position <- plot_meta$obs_position

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
                                       off_support,
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
  off_support_mode <- .ri_resolve_off_support_mode(
    stat_obs = x$stat_obs,
    stat_null = x$stat_null,
    off_support = off_support
  )

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

  tail_note <- if (identical(off_support_mode, "split")) {
    "The x-axis is split so the observed statistic and the null support remain simultaneously visible."
  } else if (identical(off_support_mode, "extend")) {
    "The x-axis is extended to include the observed statistic outside the null support."
  } else if (identical(obs_position, "left")) {
    "Observed statistic is left of null support; the entire histogram is in the right tail."
  } else if (identical(obs_position, "right")) {
    "Observed statistic is right of null support; no null draws fall in the right tail."
  } else {
    "Right-tail region is shaded; the p-value is the share of null draws at least as large as the observed statistic."
  }

  plot_data <- .ri_build_plot_data(
    stat_null = x$stat_null,
    stat_obs = x$stat_obs,
    style = style,
    bins = bins,
    show_density = show_density
  )

  if (identical(off_support_mode, "split")) {
    .ri_draw_split_distribution(
      plot_data = plot_data,
      x = x,
      main = main,
      xlab = xlab,
      ylab = ylab,
      null_col = null_col,
      tail_col = tail_col,
      density_col = density_col,
      stat_col = stat_col,
      border = border,
      lwd = lwd,
      subtitle_line_1 = subtitle_line_1,
      subtitle_line_2 = subtitle_line_2,
      tail_note = tail_note,
      ...
    )
  } else {
    xlim <- .ri_main_xlim(
      stat_obs = x$stat_obs,
      stat_null = x$stat_null,
      off_support_mode = off_support_mode
    )
    .ri_draw_single_distribution(
      plot_data = plot_data,
      x = x,
      xlim = xlim,
      main = main,
      xlab = xlab,
      ylab = ylab,
      null_col = null_col,
      tail_col = tail_col,
      density_col = density_col,
      stat_col = stat_col,
      border = border,
      lwd = lwd,
      subtitle_line_1 = subtitle_line_1,
      subtitle_line_2 = subtitle_line_2,
      tail_note = tail_note,
      obs_position = obs_position,
      ...
    )
  }

  invisible(list(
    off_support_mode = off_support_mode,
    obs_position = obs_position
  ))
}

.ri_resolve_off_support_mode <- function(stat_obs, stat_null, off_support) {
  null_range <- range(stat_null)
  if (stat_obs >= null_range[1] && stat_obs <= null_range[2]) {
    return("inside")
  }

  if (!identical(off_support, "auto")) {
    return(off_support)
  }

  support_width <- diff(null_range)
  support_width <- ifelse(support_width <= 0, 1, support_width)
  gap <- min(abs(stat_obs - null_range[1]), abs(stat_obs - null_range[2]))

  if (gap <= 0.5 * support_width) "extend" else "split"
}

.ri_build_plot_data <- function(stat_null, stat_obs, style, bins, show_density) {
  out <- list(
    style = style,
    stat_null = stat_null,
    stat_obs = stat_obs
  )

  if (identical(style, "histogram")) {
    hist_obj <- graphics::hist(stat_null, breaks = bins, plot = FALSE)
    out$hist_obj <- hist_obj
    out$dens <- if (isTRUE(show_density)) stats::density(stat_null) else NULL
    hist_y <- if (length(hist_obj$density)) hist_obj$density else 0
    dens_y <- if (is.null(out$dens)) 0 else out$dens$y
    out$y_max <- max(c(hist_y, dens_y))
  } else {
    dens <- stats::density(stat_null)
    out$dens <- dens
    out$y_max <- max(dens$y)
  }

  out
}

.ri_main_xlim <- function(stat_obs, stat_null, off_support_mode) {
  null_range <- range(stat_null)
  support_width <- diff(null_range)
  support_width <- ifelse(support_width <= 0, 1, support_width)
  null_pad <- 0.04 * support_width

  if (!identical(off_support_mode, "extend")) {
    return(c(null_range[1] - null_pad, null_range[2] + null_pad))
  }

  full_range <- range(c(stat_obs, stat_null))
  full_width <- diff(full_range)
  full_width <- ifelse(full_width <= 0, 1, full_width)
  full_pad <- 0.04 * full_width
  c(full_range[1] - full_pad, full_range[2] + full_pad)
}

.ri_draw_single_distribution <- function(plot_data, x, xlim, main, xlab, ylab,
                                         null_col, tail_col, density_col,
                                         stat_col, border, lwd,
                                         subtitle_line_1, subtitle_line_2,
                                         tail_note, obs_position, ...) {
  if (identical(plot_data$style, "histogram")) {
    fill_cols <- rep(null_col, length(plot_data$hist_obj$counts))
    fill_cols[plot_data$hist_obj$mids >= x$stat_obs] <- tail_col

    graphics::plot(
      plot_data$hist_obj,
      freq = FALSE,
      col = fill_cols,
      border = border,
      main = main,
      sub = "",
      xlab = "",
      ylab = ylab,
      xlim = xlim,
      ylim = c(0, plot_data$y_max * 1.05),
      ...
    )

    if (!is.null(plot_data$dens)) {
      graphics::lines(plot_data$dens, col = density_col, lwd = 2.5)
    }
  } else {
    graphics::plot(
      plot_data$dens,
      type = "n",
      main = main,
      sub = "",
      xlab = "",
      ylab = ylab,
      xlim = xlim,
      ylim = c(0, plot_data$y_max * 1.05),
      ...
    )

    graphics::polygon(
      x = c(plot_data$dens$x, rev(plot_data$dens$x)),
      y = c(plot_data$dens$y, rep(0, length(plot_data$dens$y))),
      col = null_col,
      border = NA
    )

    tail_idx <- plot_data$dens$x >= x$stat_obs
    if (any(tail_idx)) {
      graphics::polygon(
        x = c(x$stat_obs, plot_data$dens$x[tail_idx], max(plot_data$dens$x[tail_idx])),
        y = c(0, plot_data$dens$y[tail_idx], 0),
        col = tail_col,
        border = NA
      )
    }

    graphics::lines(plot_data$dens, col = density_col, lwd = 2)
  }

  usr <- graphics::par("usr")
  y_top <- usr[4]

  if (identical(obs_position, "inside") || identical(x$stat_obs >= usr[1] && x$stat_obs <= usr[2], TRUE)) {
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

  graphics::mtext(xlab, side = 1, line = 1.8, cex = 1)
  graphics::mtext(subtitle_line_1, side = 1, line = 3.4, cex = 0.9)
  graphics::mtext(subtitle_line_2, side = 1, line = 4.6, cex = 0.9)
  graphics::mtext(tail_note, side = 3, line = 0.2, cex = 0.78)
}

.ri_draw_split_distribution <- function(plot_data, x, main, xlab, ylab,
                                        null_col, tail_col, density_col,
                                        stat_col, border, lwd,
                                        subtitle_line_1, subtitle_line_2,
                                        tail_note, ...) {
  null_range <- range(x$stat_null)
  support_width <- diff(null_range)
  support_width <- ifelse(support_width <= 0, 1, support_width)
  obs_position <- if (x$stat_obs < null_range[1]) "left" else "right"
  null_pad <- 0.04 * support_width
  obs_pad <- 0.08 * support_width
  null_xlim <- c(null_range[1] - null_pad, null_range[2] + null_pad)
  obs_xlim <- c(x$stat_obs - obs_pad, x$stat_obs + obs_pad)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (identical(obs_position, "left")) {
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.4, 4.6))
  } else {
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4.6, 1.4))
  }
  graphics::par(oma = c(5.2, 4.2, 4.2, 1.2), mar = c(3.5, 2.5, 1.5, 0.3))

  draw_obs_panel <- function() {
    graphics::plot(
      NA,
      xlim = obs_xlim,
      ylim = c(0, plot_data$y_max * 1.05),
      type = "n",
      axes = FALSE,
      xlab = "",
      ylab = "",
      ...
    )
    graphics::axis(1, at = x$stat_obs, labels = .ri_format_num(x$stat_obs))
    graphics::box()
    graphics::segments(x$stat_obs, 0, x$stat_obs, 0.93 * plot_data$y_max,
                       col = stat_col, lwd = lwd, lty = 2)
    graphics::text(x$stat_obs, 0.98 * plot_data$y_max, labels = "Observed",
                   col = stat_col, cex = 0.85, pos = 3, xpd = NA)
    .ri_add_break_marks(side = if (identical(obs_position, "left")) "right" else "left")
  }

  draw_null_panel <- function() {
    if (identical(plot_data$style, "histogram")) {
      fill_cols <- rep(null_col, length(plot_data$hist_obj$counts))
      fill_cols[plot_data$hist_obj$mids >= x$stat_obs] <- tail_col
      graphics::plot(
        plot_data$hist_obj,
        freq = FALSE,
        col = fill_cols,
        border = border,
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        xlim = null_xlim,
        ylim = c(0, plot_data$y_max * 1.05),
        axes = FALSE,
        ...
      )
      graphics::axis(1)
      graphics::axis(2)
      graphics::box()
      if (!is.null(plot_data$dens)) {
        graphics::lines(plot_data$dens, col = density_col, lwd = 2.5)
      }
    } else {
      graphics::plot(
        plot_data$dens,
        type = "n",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        xlim = null_xlim,
        ylim = c(0, plot_data$y_max * 1.05),
        axes = FALSE,
        ...
      )
      graphics::polygon(
        x = c(plot_data$dens$x, rev(plot_data$dens$x)),
        y = c(plot_data$dens$y, rep(0, length(plot_data$dens$y))),
        col = null_col,
        border = NA
      )
      tail_idx <- plot_data$dens$x >= x$stat_obs
      if (any(tail_idx)) {
        graphics::polygon(
          x = c(x$stat_obs, plot_data$dens$x[tail_idx], max(plot_data$dens$x[tail_idx])),
          y = c(0, plot_data$dens$y[tail_idx], 0),
          col = tail_col,
          border = NA
        )
      }
      graphics::lines(plot_data$dens, col = density_col, lwd = 2)
      graphics::axis(1)
      graphics::axis(2)
      graphics::box()
    }
    .ri_add_break_marks(side = if (identical(obs_position, "left")) "left" else "right")
  }

  if (identical(obs_position, "left")) {
    draw_obs_panel()
    graphics::par(mar = c(3.5, 0.8, 1.5, 0.8))
    draw_null_panel()
  } else {
    draw_null_panel()
    graphics::par(mar = c(3.5, 0.8, 1.5, 0.3))
    draw_obs_panel()
  }

  graphics::mtext(main, side = 3, outer = TRUE, line = 1.1, cex = 1.1)
  graphics::mtext(xlab, side = 1, outer = TRUE, line = 1.6)
  graphics::mtext(ylab, side = 2, outer = TRUE, line = 2.5)
  graphics::mtext(subtitle_line_1, side = 1, outer = TRUE, line = 3.2, cex = 0.9)
  graphics::mtext(subtitle_line_2, side = 1, outer = TRUE, line = 4.3, cex = 0.9)
  graphics::mtext(tail_note, side = 3, outer = TRUE, line = -0.2, cex = 0.78)
}

.ri_add_break_marks <- function(side = c("left", "right")) {
  side <- match.arg(side)
  usr <- graphics::par("usr")
  x_span <- diff(usr[1:2])
  y_span <- diff(usr[3:4])
  x_base <- if (identical(side, "left")) usr[1] else usr[2]
  offset <- if (identical(side, "left")) 1 else -1
  x1 <- x_base + offset * 0.012 * x_span
  x2 <- x_base + offset * 0.040 * x_span
  x3 <- x_base + offset * 0.055 * x_span
  x4 <- x_base + offset * 0.083 * x_span
  y_mid <- usr[3] + 0.10 * y_span
  y_delta <- 0.035 * y_span
  graphics::segments(x1, y_mid - y_delta, x2, y_mid + y_delta, xpd = NA, lwd = 1.2)
  graphics::segments(x3, y_mid - y_delta, x4, y_mid + y_delta, xpd = NA, lwd = 1.2)
}
