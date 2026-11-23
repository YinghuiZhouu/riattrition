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
#' sharp-null p-value. It supports the package's `general`, `mp`, and `mn`
#' missingness procedures as well as an observed-only `complete-case`
#' visualization mode.
#'
#' The figure shows, from top to bottom: the title, the missingness
#' assumption, a legend (right-tail area and kernel density), the null
#' distribution with the observed statistic as a dashed line labeled
#' "Observed statistic = ...", the name of the test statistic on the x-axis, a
#' line with the observed statistic and the p-value, and a line with the design
#' and the number of MC replications. A footnote is added only when the
#' observed statistic lies outside the null distribution.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param analysis A string that specifies which inferential path to visualize:
#'   `"general"`, `"mp"`, and `"mn"` use the corresponding package procedure,
#'   while `"complete-case"` restricts the analysis to observed outcomes only.
#' @param missing A string that specifies the missing mechanism. When supplied,
#'   this should match the requested attrition-aware `analysis`. It is ignored
#'   when `analysis = "complete-case"`.
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
#' @param nperm A positive integer that specifies the number of Monte Carlo
#'   (MC) replications, i.e. random re-assignments drawn to approximate the
#'   randomization distribution (default 10^4).
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
#' @param main Optional main title. If `NULL`, the title is "Upper-tail
#'   randomization test". The missingness assumption is always shown below it.
#' @param xlab Optional x-axis label. By default the name of the test statistic,
#'   e.g. "Wilcoxon rank-sum statistic".
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
                              analysis = c("general", "mp", "mn", "complete-case"),
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
                              null_col = "#DFDED8",
                              tail_col = "#86B6EF",
                              density_col = "#3D3C38",
                              stat_col = "#D95926",
                              border = "white",
                              lwd = 2,
                              ...) {
  analysis <- match.arg(analysis)
  style <- match.arg(style)
  off_support <- match.arg(off_support)

  if (analysis %in% c("general", "mp", "mn") && !identical(missing, analysis)) {
    stop("When analysis is attrition-aware, `missing` should match the requested analysis.")
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
  viz$subtitle <- plot_meta$subtitle

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

  # Random order of the units before ranking, as in ri_sharp_attrition().
  new_order <- .ri_random_order(length(Z.analysis))
  Z.analysis <- Z.analysis[new_order]
  Y.analysis <- Y.analysis[new_order]
  if (!is.null(block.analysis)) {
    block.analysis <- block.analysis[new_order]
  }
  counts <- .ri_missing_counts(Z.analysis, Y.analysis)
  test_name <- .ri_test_label(class = class, method.list = method.list)
  if (analysis %in% c("general", "mp", "mn")) {
    Z.perm <- .ri_default_assignments(
      Z = Z.analysis, Y = Y.analysis, missing = analysis,
      block = block.analysis, nperm = nperm
    )
    stat.null <- null_dist(
      n = counts$n_total, m = counts$n_treat, class = class,
      method.list = method.list, Z.perm = Z.perm, nperm = nperm
    )
    sharp_details <- .ri_sharp_details(
      Z = Z.analysis, Y = Y.analysis, c = c, missing = analysis, class = class,
      method.list = method.list, stat.null = stat.null, Z.perm = Z.perm,
      nperm = nperm
    )
    mode_label <- switch(
      analysis,
      general = "general missingness",
      mp = "monotone positive missingness",
      mn = "monotone negative missingness"
    )
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
      assumption = if (identical(analysis, "complete-case")) "complete-case" else analysis,
      mode_label = mode_label,
      stat_null = stat.null,
      stat_obs = sharp_details$stat_obs,
      p_value = sharp_details$p_value,
      test_name = test_name,
      design = analysis_data$design,
      inference_unit = analysis_data$inference_unit,
      counts = counts,
      mc_replications = length(stat.null)
    ),
    class = "riattrition_plot"
  )
}

#----------- New helper: draw histogram or density view ----------------------#
#
# The plot highlights the right tail because the package's sharp-null p-values
# are computed as Pr(T_null >= T_obs).
#
# Layout, from top to bottom:
#   title     "Upper-tail randomization test" (or `main`)
#   subtitle  the missingness assumption, e.g. "Monotone positive missingness"
#   legend    right-tail area and kernel density
#   plot      the null distribution; the dashed line is labeled
#             "Observed statistic = ..."
#   x-axis    e.g. "Wilcoxon rank-sum statistic"
#   line 1    "Observed statistic = ... | p-value = ..."
#   line 2    "Complete randomization | MC replications = 10,000"
#   footnote  only when the observed statistic lies outside the null
#             distribution (how the x-axis was extended or split)

# Text colors: primary, secondary and muted ink.
.ri_ink <- c(primary = "#0B0B0B", secondary = "#52514E", muted = "#898781")

.ri_capitalize <- function(s) {
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

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

  # Text around the plot.
  if (is.null(main)) {
    main <- "Upper-tail randomization test"
  }
  if (is.null(xlab)) {
    xlab <- paste(x$test_name, "statistic")
  }
  if (is.null(ylab)) {
    ylab <- "Density"
  }
  text <- list(
    main = main,
    assumption = .ri_capitalize(x$mode_label),
    xlab = xlab,
    ylab = ylab,
    obs = paste0("Observed statistic = ", .ri_format_num(x$stat_obs)),
    line_1 = paste0(
      "Observed statistic = ", .ri_format_num(x$stat_obs), "  |  ",
      if (x$p_value < 0.001) "p-value < 0.001" else paste0("p-value = ", .ri_pvalue(x$p_value))
    ),
    line_2 = paste0(
      .ri_capitalize(x$design),
      "  |  MC replications = ", .ri_format_count(x$mc_replications)
    )
  )

  off_support_note <- if (identical(off_support_mode, "split")) {
    "The x-axis is split so that the observed statistic and the null distribution are both visible."
  } else if (identical(off_support_mode, "extend")) {
    "The x-axis is extended to include the observed statistic, which lies outside the null distribution."
  } else if (identical(obs_position, "left")) {
    "The observed statistic is left of the null distribution; all MC replications are in the right tail."
  } else if (identical(obs_position, "right")) {
    "The observed statistic is right of the null distribution; no MC replication is in the right tail."
  } else {
    NULL
  }
  text$footnote <- off_support_note

  plot_data <- .ri_build_plot_data(
    stat_null = x$stat_null,
    stat_obs = x$stat_obs,
    style = style,
    bins = bins,
    show_density = show_density
  )
  colors <- list(null = null_col, tail = tail_col, density = density_col,
                 stat = stat_col, border = border)

  if (identical(off_support_mode, "split")) {
    .ri_draw_split_distribution(plot_data = plot_data, x = x, text = text,
                                colors = colors, lwd = lwd, ...)
  } else {
    xlim <- .ri_main_xlim(
      stat_obs = x$stat_obs,
      stat_null = x$stat_null,
      off_support_mode = off_support_mode
    )
    .ri_draw_single_distribution(plot_data = plot_data, x = x, xlim = xlim,
                                 text = text, colors = colors, lwd = lwd,
                                 obs_position = obs_position, ...)
  }

  invisible(list(
    off_support_mode = off_support_mode,
    obs_position = obs_position,
    subtitle = c(text$line_1, text$line_2)
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

# Draws the null distribution (histogram or density) with the right tail
# shaded, plus the x- and y-axes. Used by both the single and the split layout.
.ri_draw_null_shape <- function(plot_data, stat_obs, xlim, ylim, colors, y_axis = TRUE, ...) {
  if (identical(plot_data$style, "histogram")) {
    # Bars at or above the observed statistic form the right tail.
    fill_cols <- ifelse(plot_data$hist_obj$mids >= stat_obs, colors$tail, colors$null)
    graphics::plot(
      plot_data$hist_obj, freq = FALSE, col = fill_cols, border = colors$border,
      main = "", xlab = "", ylab = "", xlim = xlim, ylim = ylim, axes = FALSE, ...
    )
    if (!is.null(plot_data$dens)) {
      graphics::lines(plot_data$dens, col = colors$density, lwd = 2)
    }
  } else {
    dens <- plot_data$dens
    graphics::plot(
      dens, type = "n", main = "", xlab = "", ylab = "", xlim = xlim, ylim = ylim,
      axes = FALSE, ...
    )
    graphics::polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))),
                      col = colors$null, border = NA)
    tail <- dens$x >= stat_obs
    if (any(tail)) {
      graphics::polygon(c(stat_obs, dens$x[tail], max(dens$x[tail])), c(0, dens$y[tail], 0),
                        col = colors$tail, border = NA)
    }
    graphics::lines(dens, col = colors$density, lwd = 2)
  }
  graphics::axis(1)
  if (y_axis) {
    graphics::axis(2)
  }
}

# Legend for the shaded right tail and the density curve. `x`, `y` give its
# top-left corner, or `x` alone gives a position such as "topleft".
.ri_draw_legend <- function(x, y = NULL, colors, with_density, yjust = 1) {
  labels <- "Right-tail area     "
  pch <- 15
  lty <- NA
  cols <- colors$tail
  if (with_density) {
    labels <- c(labels, "Kernel density")
    pch <- c(pch, NA)
    lty <- c(lty, 1)
    cols <- c(cols, colors$density)
  }
  graphics::legend(
    x = x, y = y, yjust = yjust, legend = labels, pch = pch, pt.cex = 2.2,
    lty = lty, lwd = 2, col = cols, horiz = TRUE, bty = "n", xpd = NA,
    cex = 0.85, text.col = .ri_ink[["secondary"]], seg.len = 1.6, x.intersp = 0.7
  )
}

# Dashed line at the observed statistic with the label "Observed statistic =
# ...". When the statistic lies outside the plotted range ("arrow" mode), the
# line is drawn at the plot edge with an arrow pointing toward it.
.ri_draw_observed <- function(stat_obs, label, colors, lwd) {
  usr <- graphics::par("usr")
  label_y <- usr[3] + 0.93 * (usr[4] - usr[3])
  inside <- stat_obs >= usr[1] && stat_obs <= usr[2]
  # Outside the plotted range: draw the line just inside the plot edge.
  edge_gap <- 0.03 * diff(usr[1:2])
  x_line <- if (inside) stat_obs else if (stat_obs < usr[1]) usr[1] + edge_gap else usr[2] - edge_gap

  graphics::segments(x_line, 0, x_line, label_y, col = colors$stat, lwd = lwd, lty = 2)
  if (!inside) {
    x_edge <- if (stat_obs < usr[1]) usr[1] else usr[2]
    graphics::arrows(x_line, label_y, x_edge, label_y, length = 0.08, col = colors$stat, lwd = lwd)
  }

  # Put the label on the side of the line with more room.
  side <- if (x_line > mean(usr[1:2])) 2 else 4
  graphics::text(x_line, label_y, labels = label, pos = side, cex = 0.85,
                 col = .ri_ink[["primary"]], xpd = NA)
}

# Title, subtitle, axis labels, the two summary lines and the footnote.
# `outer = TRUE` writes them in the outer margins (split layout); `shift`
# moves the lines below the plot up by that many lines.
.ri_draw_text <- function(text, outer = FALSE, shift = 0) {
  ink <- .ri_ink
  graphics::mtext(text$main, side = 3, line = 3.6, adj = 0, cex = 1.2, font = 2,
                  col = ink[["primary"]], outer = outer)
  graphics::mtext(text$assumption, side = 3, line = 2.3, adj = 0, cex = 0.95,
                  col = ink[["secondary"]], outer = outer)
  graphics::mtext(text$xlab, side = 1, line = 2.4 - shift, col = ink[["primary"]], outer = outer)
  graphics::mtext(text$ylab, side = 2, line = 3.3, las = 0, col = ink[["primary"]], outer = outer)
  graphics::mtext(text$line_1, side = 1, line = 4.0 - shift, adj = 0, cex = 0.95,
                  col = ink[["primary"]], outer = outer)
  graphics::mtext(text$line_2, side = 1, line = 5.1 - shift, adj = 0, cex = 0.85,
                  col = ink[["secondary"]], outer = outer)
  for (i in seq_along(text$footnote)) {
    graphics::mtext(text$footnote[[i]], side = 1, line = 5.2 + 0.95 * i - shift, adj = 0, cex = 0.72,
                    col = ink[["muted"]], outer = outer)
  }
}

.ri_draw_single_distribution <- function(plot_data, x, xlim, text, colors, lwd,
                                         obs_position, ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  # Margins (in lines): bottom, left, top, right. Axis lines and tick labels
  # use muted and secondary ink; y-axis tick labels are horizontal.
  graphics::par(
    mar = c(6.6 + length(text$footnote), 4.8, 5.4, 1.2),
    fg = .ri_ink[["muted"]], col.axis = .ri_ink[["secondary"]], las = 1
  )

  # Leave room above the tallest bar for the observed-statistic label.
  ylim <- c(0, plot_data$y_max * 1.18)

  .ri_draw_null_shape(plot_data, x$stat_obs, xlim, ylim, colors, ...)
  .ri_draw_observed(x$stat_obs, text$obs, colors, lwd)

  usr <- graphics::par("usr")
  .ri_draw_legend(x = usr[1], y = usr[4], yjust = 0, colors = colors,
                  with_density = !is.null(plot_data$dens))
  .ri_draw_text(text)
}

.ri_draw_split_distribution <- function(plot_data, x, text, colors, lwd, ...) {
  null_range <- range(x$stat_null)
  support_width <- diff(null_range)
  support_width <- ifelse(support_width <= 0, 1, support_width)
  obs_position <- if (x$stat_obs < null_range[1]) "left" else "right"
  null_xlim <- c(null_range[1] - 0.04 * support_width, null_range[2] + 0.04 * support_width)
  obs_xlim <- c(x$stat_obs - 0.08 * support_width, x$stat_obs + 0.08 * support_width)
  ylim <- c(0, plot_data$y_max * 1.18)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  # Two panels: the null distribution (wide) and the observed statistic
  # (narrow), side by side in the order of the x-axis.
  if (identical(obs_position, "left")) {
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.4, 4.6))
  } else {
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4.6, 1.4))
  }
  graphics::par(
    oma = c(4.4 + length(text$footnote), 4.8, 5.4, 1.2), mar = c(2.2, 0.8, 0.5, 0.8),
    fg = .ri_ink[["muted"]], col.axis = .ri_ink[["secondary"]], las = 1
  )

  draw_obs_panel <- function() {
    graphics::plot(NA, xlim = obs_xlim, ylim = ylim, type = "n", axes = FALSE,
                   xlab = "", ylab = "", ...)
    graphics::axis(1, at = x$stat_obs, labels = .ri_format_num(x$stat_obs))
    if (identical(obs_position, "left")) {
      graphics::axis(2)
    }
    graphics::box()
    usr <- graphics::par("usr")
    graphics::segments(x$stat_obs, 0, x$stat_obs, usr[3] + 0.82 * (usr[4] - usr[3]),
                       col = colors$stat, lwd = lwd, lty = 2)
    graphics::text(x$stat_obs, usr[3] + 0.92 * (usr[4] - usr[3]), labels = "Observed\nstatistic",
                   cex = 0.8, col = .ri_ink[["primary"]])
    .ri_add_break_marks(side = if (identical(obs_position, "left")) "right" else "left")
  }

  draw_null_panel <- function() {
    .ri_draw_null_shape(plot_data, x$stat_obs, null_xlim, ylim, colors,
                        y_axis = identical(obs_position, "right"), ...)
    graphics::box()
    .ri_draw_legend(x = "topleft", colors = colors, with_density = !is.null(plot_data$dens))
    .ri_add_break_marks(side = if (identical(obs_position, "left")) "left" else "right")
  }

  if (identical(obs_position, "left")) {
    draw_obs_panel()
    draw_null_panel()
  } else {
    draw_null_panel()
    draw_obs_panel()
  }

  .ri_draw_text(text, outer = TRUE, shift = 2.2)
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
