#----------- Applied-user Result Wrapper and Printing -------------------------#
#
# This file is intentionally a thin layer on top of the original inference
# engine in:
#   * R/0_function_sharp_null.R
#   * R/0_function_sharp_null_twostep.R
#
# Conceptually, the division of labor is:
#   * Old files: core identification/inference logic
#       - test_stat()
#       - null_dist()
#       - ci_sharp()
#       - ci_sharp_twostep()
#       - psi(), uci(), etc.
#   * This file: applied-user interface and presentation logic
#       - richer result object
#       - blocked-randomization / cluster-randomization helpers
#       - missingness/block diagnostics
#       - print()/summary() methods
#
# In other words, this file mostly asks:
#   "How should we package and display the underlying RI calculations?"
# rather than:
#   "What is the underlying sharp-null / two-step procedure?"
#
#' Attrition inference result for applied users
#'
#' `ri_sharp_attrition()` wraps the existing sharp-null attrition inference
#' functions and returns a richer result object with a default printed summary.
#' The printed output emphasizes general missingness by default, reports the
#' observed test statistic, and omits regression-style columns such as standard
#' errors and z-statistics.
#'
#' @param Z Treatment assignment (\eqn{n \times 1} vector).
#' @param Y Observed outcome (\eqn{n \times 1} vector, including `NA`s).
#' @param c A scalar that specifies the sharp null hypothesis.
#' @param missing A string that specifies the missing mechanism.
#' @param class A string that specifies the class of test statistic for the
#'   sharp-null procedure.
#' @param method.list A list that specifies the choice of test statistic.
#' @param block An optional block identifier vector. If supplied, permutations
#'   are generated within blocks while preserving the treated count in each
#'   block.
#' @param block_missing_sum A logical value indicating whether the printed block
#'   missingness section should include additional summaries beyond the default
#'   two-line display.
#' @param cluster An optional cluster identifier vector. If supplied, the
#'   analysis treats clusters as the unit of analysis.
#' @param cluster_outcome A string that specifies how individual outcomes are
#'   aggregated to the cluster level. Currently only `"mean_observed"` is
#'   supported.
#' @param cluster_missing A string that specifies when the cluster-level outcome
#'   is missing. `"all_missing"` treats a cluster as missing only if all
#'   individual outcomes are missing. `"any_missing"` treats a cluster as
#'   missing if any individual outcome is missing.
#' @param alternative A string that specifies the direction of the confidence
#'   interval.
#' @param stat.null An optional null distribution for the sharp-null procedure.
#' @param Z.perm An optional assignment matrix for the sharp-null procedure.
#' @param nperm A positive integer that specifies the number of permutations.
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A numerical object that specifies the precision of confidence
#'   intervals.
#' @param include_ci A logical value indicating whether confidence intervals
#'   should be computed for the printed result.
#' @param include_twostep A logical value indicating whether the two-step
#'   procedure should be included. If `NULL`, it is included for `general`,
#'   `mp`, and `mn`, and omitted otherwise.
#' @param beta A real number that belongs to \eqn{[0, \alpha]} for the two-step
#'   procedure.
#'
#' @return An object of class `riattrition_result`.
#' @export
ri_sharp_attrition <- function(Z, Y, c = 0, missing = "general", class = "RS",
                               method.list = list(name = "Wilcoxon"),
                               block = NULL,
                               block_missing_sum = TRUE,
                               cluster = NULL,
                               cluster_outcome = "mean_observed",
                               cluster_missing = "all_missing",
                               alternative = "two.sided",
                               stat.null = NULL, Z.perm = NULL, nperm = 10^4,
                               alpha = 0.05, tol = 10^(-3),
                               include_ci = TRUE, include_twostep = NULL,
                               beta = 0.1 * alpha) {
  analysis_data <- .ri_prepare_analysis_data(
    Z = Z, Y = Y, block = block, cluster = cluster,
    cluster_outcome = cluster_outcome, cluster_missing = cluster_missing
  )
  Z.analysis <- analysis_data$Z
  Y.analysis <- analysis_data$Y
  block.analysis <- analysis_data$block

  assumptions_reported <- unique(c("general", missing))
  counts <- .ri_missing_counts(Z.analysis, Y.analysis)
  block_diagnostics <- .ri_block_diagnostics(Z = Z.analysis, Y = Y.analysis, block = block.analysis)
  test_name <- .ri_test_label(class = class, method.list = method.list)
  notes <- character(0)

  # By default, include the two-step procedure only for the missingness
  # assumptions for which the package currently has dedicated support.
  if (is.null(include_twostep)) {
    include_twostep <- missing %in% c("general", "mp", "mn")
  }

  if (length(assumptions_reported) > 1 && (!is.null(Z.perm) || !is.null(stat.null))) {
    notes <- c(
      notes,
      "Custom stat.null/Z.perm inputs were not reused across multiple assumption rows; each row was recomputed with its own design-compatible permutation inputs."
    )
  }

  results_list <- vector("list", length = 0)

  for (assumption_i in assumptions_reported) {
    # The original package expects a permutation matrix Z.perm.  If the user
    # does not supply one, we construct it here using either:
    #   * complete randomization (assign_CRE() from 0_function_sharp_null.R), or
    #   * blocked randomization (.ri_assign_blocked(), added in this file).
    Z.perm_i <- if (length(assumptions_reported) == 1 && !is.null(Z.perm)) {
      Z.perm
    } else {
      .ri_default_assignments(
        Z = Z.analysis, Y = Y.analysis, missing = assumption_i,
        block = block.analysis, nperm = nperm
      )
    }

    stat.null_i <- if (length(assumptions_reported) == 1 && !is.null(stat.null)) {
      stat.null
    } else {
      NULL
    }

    # Sharp-null details are calculated by a wrapper defined below, but that
    # wrapper delegates the core rank-based calculations to:
    #   * null_dist()   from 0_function_sharp_null.R
    #   * test_stat()   from 0_function_sharp_null.R
    sharp_details <- .ri_sharp_details(
      Z = Z.analysis, Y = Y.analysis, c = c, missing = assumption_i, class = class,
      method.list = method.list, stat.null = stat.null_i, Z.perm = Z.perm_i,
      nperm = nperm
    )

    # Confidence intervals are still computed by the original package function
    # ci_sharp(); .ri_safe_ci() is only a small convenience wrapper that turns
    # failures into NA rather than interrupting printing.
    sharp_ci <- .ri_safe_ci(
      fn = ci_sharp,
      Z = Z.analysis, Y = Y.analysis, alternative = alternative, missing = assumption_i,
      class = class, method.list = method.list, stat.null = stat.null_i,
      Z.perm = Z.perm_i, nperm = nperm, alpha = alpha, tol = tol,
      include_ci = include_ci
    )

    results_list[[length(results_list) + 1]] <- data.frame(
      `Missingness assumption` = assumption_i,
      Method = "Sharp",
      `Test stat.` = test_name,
      `Test stat. value` = sharp_details$stat_obs,
      `P-value` = sharp_details$p_value,
      `CI lower` = sharp_ci$lower,
      `CI upper` = sharp_ci$upper,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (include_twostep && assumption_i %in% c("general", "mp", "mn")) {
      # The "two-step" row is optional.  The complicated identification
      # argument still lives in the original / legacy codebase; this file only
      # standardizes the resulting output so it looks parallel to the
      # sharp-null row.
      twostep_details <- tryCatch(
        .ri_twostep_details(
          Z = Z.analysis, Y = Y.analysis, c = c, missing = assumption_i,
          method.list = method.list,
          Z.perm = Z.perm_i, nperm = nperm, beta = beta
        ),
        error = function(e) e
      )

      if (inherits(twostep_details, "error")) {
        notes <- c(notes, paste("Two-step output not added for", assumption_i, ":", conditionMessage(twostep_details)))
      } else {
        # ci_sharp_twostep() is defined in 0_function_sharp_null_twostep.R.
        twostep_ci <- .ri_safe_ci(
          fn = ci_sharp_twostep,
          Z = Z.analysis, Y = Y.analysis, alternative = alternative, missing = assumption_i,
          method.list = method.list, Z.perm = Z.perm_i, nperm = nperm,
          alpha = alpha, beta = beta, tol = tol, include_ci = include_ci
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          `Missingness assumption` = assumption_i,
          Method = "Two-step",
          `Test stat.` = .ri_twostep_test_label(missing = assumption_i, method.list = method.list),
          `Test stat. value` = twostep_details$stat_obs,
          `P-value` = twostep_details$p_value,
          `CI lower` = twostep_ci$lower,
          `CI upper` = twostep_ci$upper,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    }
  }

  results <- do.call(rbind, results_list)

  # Return a richer S3 object than the original scalar p-value interface.  The
  # actual inferential content is unchanged; we are just storing enough metadata
  # to print something that applied users can immediately interpret.
  structure(
    list(
      call = match.call(),
      assumption = missing,
      assumptions_reported = assumptions_reported,
      test_class = class,
      test_name = test_name,
      null_hypothesis = c,
      counts = counts,
      design = analysis_data$design,
      inference_unit = analysis_data$inference_unit,
      cluster_outcome = analysis_data$cluster_outcome,
      cluster_missing = analysis_data$cluster_missing,
      cluster_diagnostics = analysis_data$cluster_diagnostics,
      block_diagnostics = block_diagnostics,
      block_missing_sum = block_missing_sum,
      assumption_hint = .ri_assumption_hint(counts),
      results = results,
      confidence_level = 1 - alpha,
      include_ci = include_ci,
      alternative = alternative,
      notes = notes
    ),
    class = "riattrition_result"
  )
}

# Backward-compatible alias.  The preferred user-facing function name is now
# ri_sharp_attrition().
ri_test <- ri_sharp_attrition

#----------- New helper: choose individual- vs cluster-level analysis ---------#
#
# If cluster is NULL, the package behaves exactly as before.  If cluster is
# supplied, this helper converts individual-level inputs into cluster-level
# vectors and all downstream inference treats clusters as the units.
.ri_prepare_analysis_data <- function(Z, Y, block, cluster,
                                      cluster_outcome = "mean_observed",
                                      cluster_missing = "all_missing") {
  if (is.null(cluster)) {
    return(list(
      Z = Z,
      Y = Y,
      block = block,
      design = if (is.null(block)) "complete randomization" else "blocked randomization",
      inference_unit = "individual",
      cluster_outcome = NULL,
      cluster_missing = NULL,
      cluster_diagnostics = NULL
    ))
  }

  cluster_data <- .ri_cluster_data(
    Z = Z, Y = Y, block = block, cluster = cluster,
    cluster_outcome = cluster_outcome, cluster_missing = cluster_missing
  )

  list(
    Z = cluster_data$Z,
    Y = cluster_data$Y,
    block = cluster_data$block,
    design = if (is.null(cluster_data$block)) "cluster randomization" else "blocked cluster randomization",
    inference_unit = "cluster",
    cluster_outcome = cluster_outcome,
    cluster_missing = cluster_missing,
    cluster_diagnostics = .ri_cluster_diagnostics(cluster_data)
  )
}

#----------- New helper: aggregate individual data to cluster level -----------#
#
# Default behavior follows common cluster-trial practice:
#   * cluster_outcome = "mean_observed": cluster outcome is the mean among
#     observed individual outcomes in that cluster.
#   * cluster_missing = "all_missing": cluster outcome is missing only if the
#     cluster has no observed individual outcomes.
#
# The sensitivity option cluster_missing = "any_missing" treats a cluster as
# missing whenever any individual outcome in that cluster is missing.
.ri_cluster_data <- function(Z, Y, block, cluster,
                             cluster_outcome = "mean_observed",
                             cluster_missing = "all_missing") {
  if (!cluster_outcome %in% "mean_observed") {
    stop("Currently, cluster_outcome must be 'mean_observed'.")
  }
  if (!cluster_missing %in% c("all_missing", "any_missing")) {
    stop("cluster_missing must be 'all_missing' or 'any_missing'.")
  }
  if (length(cluster) != length(Z) || length(Y) != length(Z)) {
    stop("Z, Y, and cluster must have the same length.")
  }
  if (!is.null(block) && length(block) != length(Z)) {
    stop("block must have the same length as Z when supplied.")
  }

  cluster_chr <- as.character(cluster)
  clusters <- unique(cluster_chr)

  cluster_table <- data.frame(
    cluster = clusters,
    n_individuals = NA_integer_,
    n_observed = NA_integer_,
    n_missing = NA_integer_,
    response_rate = NA_real_,
    Z = NA_real_,
    Y = NA_real_,
    block = if (is.null(block)) NA_character_ else NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(clusters)) {
    cl <- clusters[[i]]
    idx <- which(cluster_chr == cl)
    Zi <- unique(stats::na.omit(Z[idx]))
    if (length(Zi) != 1) {
      stop("Treatment assignment must be constant within each cluster.")
    }

    if (!is.null(block)) {
      bi <- unique(stats::na.omit(as.character(block[idx])))
      if (length(bi) != 1) {
        stop("Block assignment must be constant within each cluster when block is supplied.")
      }
      cluster_table$block[i] <- bi[[1]]
    }

    Yi <- Y[idx]
    n_individuals <- length(idx)
    n_observed <- sum(!is.na(Yi))
    n_missing <- n_individuals - n_observed
    response_rate <- n_observed / n_individuals

    cluster_y <- if (n_observed > 0) mean(Yi, na.rm = TRUE) else NA_real_
    if (cluster_missing == "any_missing" && n_missing > 0) {
      cluster_y <- NA_real_
    }

    cluster_table$n_individuals[i] <- n_individuals
    cluster_table$n_observed[i] <- n_observed
    cluster_table$n_missing[i] <- n_missing
    cluster_table$response_rate[i] <- response_rate
    cluster_table$Z[i] <- Zi[[1]]
    cluster_table$Y[i] <- cluster_y
  }

  list(
    Z = cluster_table$Z,
    Y = cluster_table$Y,
    block = if (is.null(block)) NULL else cluster_table$block,
    cluster_table = cluster_table,
    n_individuals = length(Z),
    n_clusters = length(clusters),
    cluster_outcome = cluster_outcome,
    cluster_missing = cluster_missing
  )
}

#----------- New helper: cluster-level diagnostics ---------------------------#
#
# These diagnostics describe individual-level missingness inside randomized
# clusters.  They do not alter the cluster-level inference; they help users see
# how the cluster-level outcome was constructed.
.ri_cluster_diagnostics <- function(cluster_data) {
  tab <- cluster_data$cluster_table

  list(
    n_individuals = cluster_data$n_individuals,
    n_clusters = cluster_data$n_clusters,
    mean_cluster_size = mean(tab$n_individuals),
    median_cluster_size = stats::median(tab$n_individuals),
    mean_cluster_response_rate = mean(tab$response_rate),
    median_cluster_response_rate = stats::median(tab$response_rate),
    clusters_with_no_observed_outcome = sum(tab$n_observed == 0),
    clusters_with_any_missing = sum(tab$n_missing > 0),
    cluster_table = tab
  )
}

#----------- New helper: missingness summary ---------------------------------#
#
# This helper is new in the applied-user wrapper.  It does not perform any
# inference.  Its only job is to summarize how much missingness/attrition exists
# overall and by treatment status.
.ri_missing_counts <- function(Z, Y) {
  M <- as.numeric(!is.na(Y))
  n_total <- length(Z)
  n_treat <- sum(Z)
  n_control <- n_total - n_treat
  n_obs <- sum(M)
  n_miss <- n_total - n_obs
  n_obs_treat <- sum(Z[M == 1])
  n_obs_control <- n_obs - n_obs_treat
  n_miss_treat <- sum(Z[M == 0])
  n_miss_control <- n_miss - n_miss_treat
  attrition_rate_total <- n_miss / n_total
  attrition_rate_treat <- n_miss_treat / n_treat
  attrition_rate_control <- n_miss_control / n_control

  list(
    n_total = n_total,
    n_treat = n_treat,
    n_control = n_control,
    n_obs = n_obs,
    n_miss = n_miss,
    n_obs_treat = n_obs_treat,
    n_obs_control = n_obs_control,
    n_miss_treat = n_miss_treat,
    n_miss_control = n_miss_control,
    attrition_rate_total = attrition_rate_total,
    attrition_rate_treat = attrition_rate_treat,
    attrition_rate_control = attrition_rate_control
  )
}

#----------- New helper: interpretation note for mp vs mn --------------------#
#
# This is also presentation-only logic.  It gives the user a gentle reminder
# about whether the observed attrition pattern is more compatible with mp or mn.
# Importantly, it does NOT claim to identify the true missingness mechanism.
.ri_assumption_hint <- function(counts, tol = 1e-8) {
  diff <- counts$attrition_rate_treat - counts$attrition_rate_control

  if (abs(diff) <= tol) {
    return("Observed treated and control attrition rates are very similar; rates alone do not favor mp or mn.")
  }

  if (diff < 0) {
    return("Treated attrition is lower than control attrition; this is more consistent with mp than mn, although rates alone do not identify the missingness mechanism.")
  }

  "Treated attrition is higher than control attrition; this is more consistent with mn than mp, although rates alone do not identify the missingness mechanism."
}

#----------- New helper: block-level attrition diagnostics --------------------#
#
# This helper exists because the original package assumed complete randomization.
# For blocked designs, users usually want to know not just pooled attrition, but
# also whether attrition differences are systematically positive/negative within
# block.
.ri_block_diagnostics <- function(Z, Y, block) {
  if (is.null(block)) {
    return(NULL)
  }

  block <- as.character(block)
  M <- as.numeric(!is.na(Y))
  blocks <- unique(block)

  block_table <- data.frame(
    block = blocks,
    n = NA_integer_,
    treated_n = NA_integer_,
    control_n = NA_integer_,
    treated_attrition_rate = NA_real_,
    control_attrition_rate = NA_real_,
    attrition_diff = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    idx <- which(block == b)
    Zb <- Z[idx]
    Mb <- M[idx]

    treated_n <- sum(Zb == 1)
    control_n <- sum(Zb == 0)
    treated_missing <- sum(Zb == 1 & Mb == 0)
    control_missing <- sum(Zb == 0 & Mb == 0)

    treated_rate <- if (treated_n > 0) treated_missing / treated_n else NA_real_
    control_rate <- if (control_n > 0) control_missing / control_n else NA_real_

    block_table$n[i] <- length(idx)
    block_table$treated_n[i] <- treated_n
    block_table$control_n[i] <- control_n
    block_table$treated_attrition_rate[i] <- treated_rate
    block_table$control_attrition_rate[i] <- control_rate
    block_table$attrition_diff[i] <- treated_rate - control_rate
  }

  valid_diff <- block_table$attrition_diff[!is.na(block_table$attrition_diff)]
  q <- if (length(valid_diff) > 0) {
    stats::quantile(valid_diff, probs = c(0.25, 0.5, 0.75), names = FALSE)
  } else {
    c(NA_real_, NA_real_, NA_real_)
  }

  list(
    n_blocks = length(blocks),
    # These summaries are deliberately simple: they provide a quick diagnostic
    # of the direction of within-block attrition imbalance.
    mean_within_block_attrition_diff = if (length(valid_diff) > 0) mean(valid_diff) else NA_real_,
    median_within_block_attrition_diff = if (length(valid_diff) > 0) stats::median(valid_diff) else NA_real_,
    min_within_block_attrition_diff = if (length(valid_diff) > 0) min(valid_diff) else NA_real_,
    q1_within_block_attrition_diff = q[[1]],
    q3_within_block_attrition_diff = q[[3]],
    max_within_block_attrition_diff = if (length(valid_diff) > 0) max(valid_diff) else NA_real_,
    block_table = block_table
  )
}

#----------- New helper: choose the permutation scheme ------------------------#
#
# This is the main design-aware addition in this file.
#
# Old behavior:
#   * complete randomization only, via assign_CRE()
#
# New behavior:
#   * if block is supplied, generate permutations within block while preserving
#     the treated count in each block.
#
# For missing = general/mp/mn, inference is defined on all randomized units.
# For other cases, inference is restricted to observed units, matching the logic
# in the legacy functions.
.ri_default_assignments <- function(Z, Y, missing, block, nperm) {
  M <- as.numeric(!is.na(Y))

  if (missing %in% c("general", "mp", "mn")) {
    if (is.null(block)) {
      # assign_CRE() comes from 0_function_sharp_null.R.
      return(assign_CRE(length(Z), sum(Z), nperm))
    }
    return(.ri_assign_blocked(Z = Z, block = block, nperm = nperm))
  }

  Z.obs <- Z[M == 1]
  if (is.null(block)) {
    return(assign_CRE(length(Z.obs), sum(Z.obs), nperm))
  }
  .ri_assign_blocked(Z = Z.obs, block = block[M == 1], nperm = nperm)
}

#----------- New helper: blocked permutation matrix --------------------------#
#
# This function has no analogue in the original complete-randomization code.
# It implements the natural randomization distribution for a blocked experiment:
#   * blocks stay fixed,
#   * the treated count within each block stays fixed,
#   * treatment labels are shuffled only inside each block.
.ri_assign_blocked <- function(Z, block, nperm) {
  block <- as.character(block)
  blocks <- unique(block)
  n <- length(Z)
  Z.perm <- matrix(0, nrow = n, ncol = nperm)

  block_index <- lapply(blocks, function(b) which(block == b))
  block_treat <- vapply(block_index, function(idx) sum(Z[idx]), numeric(1))

  for (iter in seq_len(nperm)) {
    draw <- integer(n)
    for (j in seq_along(block_index)) {
      idx <- block_index[[j]]
      m <- block_treat[[j]]
      if (m > 0) {
        draw[sample(idx, m, replace = FALSE)] <- 1L
      }
    }
    Z.perm[, iter] <- draw
  }

  Z.perm
}

#----------- New helper: user-facing test statistic label ---------------------#
#
# The old package exposes class/method internals.  This helper turns those into
# a readable label for printing.
.ri_test_label <- function(class, method.list) {
  if (class == "RS" && method.list$name == "Wilcoxon") {
    return("Wilcoxon rank-sum")
  }
  if (class == "RS" && method.list$name == "Stephenson") {
    return("Stephenson rank-sum")
  }
  if (class == "MWU+" && method.list$name == "Wilcoxon") {
    return("Wilcoxon MWU+")
  }
  if (class == "MWU-" && method.list$name == "Wilcoxon") {
    return("Wilcoxon MWU-")
  }
  if (class == "MWU+" && method.list$name == "Polynomial") {
    return("Polynomial MWU+")
  }
  if (class == "MWU-" && method.list$name == "Polynomial") {
    return("Polynomial MWU-")
  }

  paste(class, method.list$name)
}

.ri_twostep_test_label <- function(missing, method.list) {
  if (missing %in% c("general", "mp")) {
    return(paste(method.list$name, "MWU+"))
  }
  if (missing == "mn") {
    return(paste(method.list$name, "MWU-"))
  }

  paste(method.list$name, "Two-step")
}

#----------- Bridge to legacy sharp-null engine -------------------------------#
#
# This helper is where the new wrapper meets the old inference engine.
#
# Inputs:
#   * Z, Y, c, missing, class, method.list, Z.perm
# Output:
#   * observed test statistic
#   * p-value
#   * missingness counts
#
# What is inherited from 0_function_sharp_null.R?
#   * null_dist()
#   * test_stat()
#
# What is new here?
#   * user-facing bookkeeping
#   * explicit construction of the completed outcome under each missingness
#     assumption so it can be summarized and routed into the original functions
.ri_sharp_details <- function(Z, Y, c, missing, class, method.list, stat.null, Z.perm, nperm) {
  counts <- .ri_missing_counts(Z, Y)
  M <- as.numeric(!is.na(Y))
  n <- counts$n_total
  n1 <- counts$n_treat

  # Under the sharp null tau = c, these are the imputed potential outcomes:
  #   Y1 = Y0 + c
  #   Y0 = Y1 - c
  # Given observed Y and assignment Z, these formulas recover the missing
  # potential outcome for observed units.
  Y1.imp <- Y + (1 - Z) * c
  Y0.imp <- Y - Z * c

  if (missing %in% c("general", "mp", "mn")) {
    # These default bounds encode the "worst-case" completion rule under each
    # assumption.  They are not new theory; they are a compact way of expressing
    # the completion step before handing the completed outcome to test_stat().
    #
    # b00: units missing under both treatment states
    # b01: observed if treated, missing if control
    # b10: missing if treated, observed if control
    defaults <- switch(
      missing,
      general = list(b00 = 0, b01 = Inf, b10 = -Inf),
      mp = list(b00 = Inf, b01 = Inf, b10 = 0),
      mn = list(b00 = -Inf, b01 = 0, b10 = -Inf)
    )

    Y0.com.imp <- rep(NA_real_, n)

    if (missing == "general") {
      # general:
      #   * observed treated outcomes contribute their imputed control value
      #   * missing treated outcomes are assigned the least favorable value -Inf
      #   * observed control outcomes contribute their observed/imputed control value
      #   * missing control outcomes are assigned the most favorable value Inf
      Y0.com.imp[Z == 1 & M == 1] <- pmin(Y0.imp[Z == 1 & M == 1], defaults$b01)
      Y0.com.imp[Z == 1 & M == 0] <- min(defaults$b00, defaults$b10)
      Y0.com.imp[Z == 0 & M == 1] <- pmax(Y0.imp[Z == 0 & M == 1], defaults$b10)
      Y0.com.imp[Z == 0 & M == 0] <- max(defaults$b00, defaults$b01)
    }

    if (missing == "mp") {
      # mp ("monotone positive"): treatment weakly increases response.
      # Relative to general missingness, the completion reflects the direction
      # restriction M1 >= M0.
      Y0.com.imp[Z == 1 & M == 1] <- pmin(Y0.imp[Z == 1 & M == 1], defaults$b01)
      Y0.com.imp[Z == 1 & M == 0] <- defaults$b00
      Y0.com.imp[Z == 0 & M == 1] <- Y0.imp[Z == 0 & M == 1]
      Y0.com.imp[Z == 0 & M == 0] <- max(defaults$b00, defaults$b01)
    }

    if (missing == "mn") {
      # mn ("monotone negative"): treatment weakly decreases response.
      # This is the mirror-image directional restriction M1 <= M0.
      Y0.com.imp[Z == 1 & M == 1] <- Y0.imp[Z == 1 & M == 1]
      Y0.com.imp[Z == 1 & M == 0] <- min(defaults$b00, defaults$b10)
      Y0.com.imp[Z == 0 & M == 1] <- pmax(Y0.imp[Z == 0 & M == 1], defaults$b10)
      Y0.com.imp[Z == 0 & M == 0] <- defaults$b00
    }

    if (is.null(stat.null)) {
      # null_dist() is imported from 0_function_sharp_null.R.  At this point all
      # design information has already been folded into Z.perm.
      stat.null <- null_dist(n, n1, class = class, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
    }

    # test_stat() is also imported from 0_function_sharp_null.R.
    stat.obs <- test_stat(Z = Z, Y = Y0.com.imp, class = class, method.list = method.list)
    p.value <- mean(stat.null >= stat.obs)

    return(list(
      stat_obs = stat.obs,
      p_value = p.value,
      counts = counts
    ))
  }

  Y0.imp.obs <- Y0.imp[M == 1]
  Z.obs <- Z[M == 1]
  n.obs <- length(Z.obs)
  n1.obs <- sum(Z.obs)

  if (is.null(stat.null)) {
    # In the non-general branch the original procedure works only on units with
    # observed outcomes, so we pass the reduced observed sample to null_dist().
    stat.null <- null_dist(n.obs, n1.obs, class = class, method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  }

  stat.obs <- test_stat(Z = Z.obs, Y = Y0.imp.obs, class = class, method.list = method.list)
  p.value <- mean(stat.null >= stat.obs)

  list(
    stat_obs = stat.obs,
    p_value = p.value,
    counts = counts
  )
}

#----------- Bridge to legacy two-step engine --------------------------------#
#
# This dispatcher is new, but each branch below still relies on identification
# ingredients from the original two-step implementation:
#   * psi()
#   * uci()
#   * null_dist()
#   * test_stat()
.ri_twostep_details <- function(Z, Y, c, missing, method.list, Z.perm, nperm, beta) {
  if (missing == "general") {
    return(.ri_twostep_general_details(Z, Y, c, method.list, Z.perm, nperm, beta))
  }
  if (missing == "mp") {
    return(.ri_twostep_mp_details(Z, Y, c, method.list, Z.perm, nperm, beta))
  }
  if (missing == "mn") {
    return(.ri_twostep_mn_details(Z, Y, c, method.list, Z.perm, nperm, beta))
  }

  stop("Two-step output is currently supported only for general, mp, and mn missingness.")
}

#----------- General-missingness two-step row --------------------------------#
#
# This is the most algebraically involved part of the file.  The formulas are
# not invented here; they are the applied-user wrapper's implementation of the
# existing two-step procedure, with the final output standardized into a simple
# list(stat_obs, p_value).
.ri_twostep_general_details <- function(Z, Y, c, method.list, Z.perm, nperm, beta) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z)
  n.obs <- length(Z[M == 1])
  n1 <- sum(Z)
  n0 <- n - n1
  n11 <- sum(Z[M == 1])
  n10 <- sum(Z[M == 0])
  n01 <- n.obs - n11

  beta1 <- beta / 2
  beta2 <- beta / 2

  if (beta > 0) {
    # WhyperCI_M() gives a confidence region for the number of always-observed /
    # partially observed units.  This is one of the steps that accounts for
    # uncertainty about missingness patterns before computing the worst-case test
    # statistic.
    M.mat <- ExactCIone::WhyperCI_M(x = n01, n = n0, N = n, conf.level = 1 - beta1)
    M1hat <- M.mat$CI[1, 2]
    M2hat <- M.mat$CI[1, 3]
  } else {
    M1hat <- n01
    M2hat <- n1 + n01
  }

  mlbar <- M1hat - n01
  mubar <- M2hat - n01

  m11.grid <- 0:n.obs
  q.hg.vec <- stats::qhyper(p = 1 - beta2, m = m11.grid, n = n - m11.grid, k = n0)
  dubar <- max(n / (n1 * n0) * q.hg.vec - m11.grid / n1)

  ind.treat <- which(Z == 1)
  ind.control <- which(Z == 0)
  ind.control.obs <- which(Z == 0 & M == 1)
  ind.control.miss <- which(Z == 0 & M == 0)
  ind.treat.obs <- which(Z == 1 & M == 1)
  ind.treat.miss <- which(Z == 1 & M == 0)

  r <- rank(Y, ties.method = "first")
  ind.sort <- sort.int(r, index.return = TRUE)$ix
  ind.sort.control.obs <- ind.sort[Z[ind.sort] == 0 & M[ind.sort] == 1]

  phi <- .ri_phi(method.list)

  A <- rep(NA_real_, n)
  A[ind.treat] <- n01 + findInterval(ind.treat, ind.control.miss)

  B.J.precompute.10 <- findInterval(ind.treat.miss, ind.control.obs)

  if (n11 > 0 && n01 > 0) {
    psi.mat <- matrix(0, nrow = n11, ncol = n01)
    for (a in 1:n11) {
      i <- ind.treat.obs[a]
      for (b in 1:n01) {
        j <- ind.sort.control.obs[b]
        psi.mat[a, b] <- psi(i, j, Y[i] - c, Y[j])
      }
    }
    B.J.precompute.11 <- matrix(0, nrow = n11, ncol = n01 + 1)
    psi.rev <- psi.mat[, n01:1, drop = FALSE]
    B.J.precompute.11[, 2:(n01 + 1)] <- t(apply(psi.rev, 1, cumsum))
  } else {
    B.J.precompute.11 <- matrix(0, nrow = n11, ncol = n01 + 1)
  }

  Klbar <- max(0, mlbar - n10)
  Kubar <- min(n11, mubar)

  TK.vec <- rep(NA_real_, Kubar - Klbar + 1)
  for (K in Klbar:Kubar) {
    # Each K indexes a feasible latent configuration of missing potential
    # outcomes.  We evaluate the statistic over that feasible set and take the
    # worst case below.
    J <- min(floor(n0 * (dubar + K / n1)), n01)
    L <- min(mubar - K, n10)

    B.J <- rep(NA_real_, n)
    if (n11 > 0) {
      B.J[ind.treat.obs] <- B.J.precompute.11[, J + 1] + (n01 - J)
    }
    if (n10 > 0) {
      B.J[ind.treat.miss] <- pmax(0, B.J.precompute.10 - J)
    }

    C.J <- phi(A) - phi(B.J)
    r.C <- rank(C.J, ties.method = "first")
    ind.sort.C <- sort.int(r.C, index.return = TRUE)$ix
    ind.sort.C.treat.obs <- ind.sort.C[Z[ind.sort.C] == 1 & M[ind.sort.C] == 1]
    ind.sort.C.treat.miss <- ind.sort.C[Z[ind.sort.C] == 1 & M[ind.sort.C] == 0]

    TK.vec[K - Klbar + 1] <-
      sum(phi(B.J)[Z == 1]) +
      sum(C.J[ind.sort.C.treat.obs[seq_len(n11 - K)]]) +
      sum(C.J[ind.sort.C.treat.miss[seq_len(n10 - L)]])
  }

  # The two-step statistic is the most conservative value across feasible latent
  # configurations.  The null distribution itself is still obtained from the
  # original null_dist() function.
  stat.obs <- min(TK.vec)
  stat.null <- null_dist(n, n1, class = "MWU+", method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  p.value <- min(mean(stat.null >= stat.obs) + beta, 1)

  list(stat_obs = stat.obs, p_value = p.value)
}

#----------- Monotone-positive two-step row ----------------------------------#
.ri_twostep_mp_details <- function(Z, Y, c, method.list, Z.perm, nperm, beta) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z)
  n.obs <- length(Z[M == 1])
  n1 <- sum(Z)
  n0 <- n - n1
  n11 <- sum(Z[M == 1])
  n01 <- n.obs - n11

  Mhat <- uci(N = n, x = n01, n = n0, alpha = beta)
  mlbar <- max(n11 + n01 - Mhat, 0)

  # Under mp, missing outcomes are completed in the direction implied by
  # M1 >= M0, and the ranking-based statistic is then computed using MWU+.
  Y.comp <- ifelse(M == 1, Y, Inf)
  ind.treat.obs <- which(Z == 1 & M == 1)
  ind.control <- which(Z == 0)

  A <- rep(NA_real_, n)
  B <- rep(NA_real_, n)
  for (i in ind.treat.obs) {
    A[i] <- sum(psi(i = i, j = ind.control, x = Inf, y = Y.comp[ind.control]))
    B[i] <- sum(psi(i = i, j = ind.control, x = Y[i] - c, y = Y.comp[ind.control]))
  }

  phi <- .ri_phi(method.list)
  C <- phi(A) - phi(B)
  r <- rank(C, ties.method = "first")
  ind.sort <- sort.int(r, index.return = TRUE)$ix
  ind.sort.treat.obs <- ind.sort[Z[ind.sort] == 1 & M[ind.sort] == 1]

  Y0.com.imp <- rep(NA_real_, n)
  if (mlbar > 0) {
    Y0.com.imp[ind.sort.treat.obs[seq_len(mlbar)]] <- Inf
  }
  if (mlbar < n11) {
    Y0.com.imp[ind.sort.treat.obs[(mlbar + 1):n11]] <- Y[ind.sort.treat.obs[(mlbar + 1):n11]] - c
  }
  Y0.com.imp[Z == 1 & M == 0] <- Inf
  Y0.com.imp[Z == 0 & M == 1] <- Y[Z == 0 & M == 1]
  Y0.com.imp[Z == 0 & M == 0] <- Inf

  stat.null <- null_dist(n, n1, class = "MWU+", method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  stat.obs <- test_stat(Z = Z, Y = Y0.com.imp, class = "MWU+", method.list = method.list)
  p.value <- min(mean(stat.null >= stat.obs) + beta, 1)

  list(stat_obs = stat.obs, p_value = p.value)
}

#----------- Monotone-negative two-step row ----------------------------------#
.ri_twostep_mn_details <- function(Z, Y, c, method.list, Z.perm, nperm, beta) {
  M <- as.numeric(!is.na(Y))
  n <- length(Z)
  n.obs <- length(Z[M == 1])
  n1 <- sum(Z)
  n11 <- sum(Z[M == 1])
  n01 <- n.obs - n11

  Mhat <- uci(N = n, x = n11, n = n1, alpha = beta)
  mlbar <- max(n11 + n01 - Mhat, 0)

  # Under mn, the directional restriction is reversed and the code uses MWU-.
  Y.comp <- ifelse(M == 1, Y - c, -Inf)
  ind.control.obs <- which(Z == 0 & M == 1)
  ind.treat <- which(Z == 1)

  A <- rep(NA_real_, n)
  B <- rep(NA_real_, n)
  for (i in ind.control.obs) {
    A[i] <- sum(psi(i = i, j = ind.treat, x = Y[i], y = Y.comp[ind.treat]))
    B[i] <- sum(psi(i = i, j = ind.treat, x = -Inf, y = Y.comp[ind.treat]))
  }

  phi <- .ri_phi(method.list)
  C <- phi(A) - phi(B)
  r <- rank(C, ties.method = "first")
  ind.sort <- sort.int(r, index.return = TRUE)$ix
  ind.sort.control.obs <- ind.sort[Z[ind.sort] == 0 & M[ind.sort] == 1]

  Y0.com.imp <- rep(NA_real_, n)
  Y0.com.imp[Z == 1 & M == 1] <- Y[Z == 1 & M == 1] - c
  Y0.com.imp[Z == 1 & M == 0] <- -Inf
  if (mlbar < n01) {
    Y0.com.imp[ind.sort.control.obs[(mlbar + 1):n01]] <- Y[ind.sort.control.obs[(mlbar + 1):n01]]
  }
  if (mlbar > 0) {
    Y0.com.imp[ind.sort.control.obs[seq_len(mlbar)]] <- -Inf
  }
  Y0.com.imp[Z == 0 & M == 0] <- -Inf

  stat.null <- null_dist(n, n1, class = "MWU-", method.list = method.list, Z.perm = Z.perm, nperm = nperm)
  stat.obs <- test_stat(Z = Z, Y = Y0.com.imp, class = "MWU-", method.list = method.list)
  p.value <- min(mean(stat.null >= stat.obs) + beta, 1)

  list(stat_obs = stat.obs, p_value = p.value)
}

#----------- New helper: transform for MWU-type two-step statistics ----------#
#
# This helper is local to the applied-user wrapper.  It converts method.list
# into a scalar transformation phi() used inside the MWU-style two-step code.
.ri_phi <- function(method.list) {
  if (method.list$name == "Wilcoxon") {
    return(function(x) x)
  }
  if (method.list$name == "Polynomial") {
    return(function(x) x^(method.list$s - 1))
  }
  stop("Unsupported method.list for two-step output.")
}

#----------- New helper: fail-soft CI handling -------------------------------#
#
# The original CI functions can error in edge cases.  For printing, it is often
# preferable to display NA rather than abort the whole call.
.ri_safe_ci <- function(fn, ..., include_ci) {
  if (!include_ci) {
    return(list(lower = NA_real_, upper = NA_real_))
  }

  ci <- tryCatch(fn(...), error = function(e) e)
  if (inherits(ci, "error")) {
    return(list(lower = NA_real_, upper = NA_real_))
  }

  if (length(ci) == 1) {
    return(list(lower = ci[[1]], upper = NA_real_))
  }

  list(lower = ci[[1]], upper = ci[[2]])
}

#----------- New formatting helpers ------------------------------------------#
#
# These helpers are purely cosmetic.  They never affect inference.
.ri_format_num <- function(x, digits = 3) {
  if (is.na(x)) {
    return("NA")
  }
  formatC(x, digits = digits, format = "f")
}

.ri_format_ci <- function(lower, upper, digits = 3) {
  if (is.na(lower) && is.na(upper)) {
    return("NA")
  }
  if (is.na(upper)) {
    return(paste0("[", .ri_format_num(lower, digits), ", NA]"))
  }
  paste0("[", .ri_format_num(lower, digits), ", ", .ri_format_num(upper, digits), "]")
}

.ri_cluster_outcome_label <- function(cluster_outcome) {
  if (identical(cluster_outcome, "mean_observed")) {
    return("Mean of observed individual outcomes")
  }
  cluster_outcome
}

.ri_cluster_missing_label <- function(cluster_missing) {
  if (identical(cluster_missing, "all_missing")) {
    return("Missing only if all individual outcomes are missing")
  }
  if (identical(cluster_missing, "any_missing")) {
    return("Missing if any individual outcome is missing")
  }
  cluster_missing
}

.ri_summary_table <- function(x) {
  out <- data.frame(
    Assumption = x$results$`Missingness assumption`,
    Null = x$results$Method,
    `Test stat.` = x$results$`Test stat.`,
    Value = vapply(x$results$`Test stat. value`, .ri_format_num, character(1)),
    `P-value` = vapply(x$results$`P-value`, .ri_format_num, character(1)),
    `CI` = vapply(
      seq_len(nrow(x$results)),
      function(i) .ri_format_ci(x$results$`CI lower`[i], x$results$`CI upper`[i]),
      character(1)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!isTRUE(x$include_ci)) {
    out$CI <- NULL
  }

  out
}

#----------- S3 print/summary methods ----------------------------------------#
#
# These are entirely new relative to the legacy package interface.  They turn
# the stored result object into an "rdrobust-style" printed summary emphasizing:
#   * missingness summary
#   * test statistic
#   * p-value
#   * confidence interval
# while omitting regression-style columns such as standard errors and z-stats.
#' @export
print.riattrition_result <- function(x, ...) {
  counts <- x$counts
  summary.table <- .ri_summary_table(x)
  label_width <- 24
  small_col_width <- 10
  wide_label_width <- 38
  if (isTRUE(x$include_ci)) {
    ci.label <- paste0("[", round(100 * x$confidence_level), "% C.I.]")
    names(summary.table)[6] <- ci.label
  }

  cat("Call: ", deparse(x$call), "\n\n", sep = "")
  cat("Inference setup\n", sep = "")
  cat(sprintf("%-*s %s\n", label_width, "Inference unit", x$inference_unit))
  cat(sprintf("%-*s %s\n", label_width, "Design", x$design))
  cat(sprintf("%-*s %s\n", label_width, "Test statistic", x$test_name))
  cat(sprintf("%-*s %s\n\n", label_width, "Null hypothesis", paste0("tau = ", x$null_hypothesis)))

  cat("Overall missingness\n", sep = "")
  cat(sprintf("%-*s %-*s %-*s %-*s\n", label_width, "", small_col_width, "Total", small_col_width, "T", small_col_width, "C"))
  obs_label <- if (identical(x$inference_unit, "cluster")) "Number of clusters" else "Number of Obs."
  miss_label <- if (identical(x$inference_unit, "cluster")) "Missing clusters" else "Missing units"
  cat(sprintf("%-*s %-*s %-*s %-*s\n",
              label_width,
              obs_label,
              small_col_width,
              counts$n_total,
              small_col_width,
              counts$n_treat,
              small_col_width,
              counts$n_control))
  cat(sprintf("%-*s %-*s %-*s %-*s\n",
              label_width,
              miss_label,
              small_col_width,
              counts$n_miss,
              small_col_width,
              counts$n_miss_treat,
              small_col_width,
              counts$n_miss_control))
  cat(sprintf("%-*s %-*s %-*s %-*s\n\n",
              label_width,
              "Attrition rate",
              small_col_width,
              .ri_format_num(counts$attrition_rate_total),
              small_col_width,
              .ri_format_num(counts$attrition_rate_treat),
              small_col_width,
              .ri_format_num(counts$attrition_rate_control)))

  if (!is.null(x$cluster_diagnostics)) {
    cat("Cluster summary\n", sep = "")
    cat(sprintf("%-*s %s\n", wide_label_width, "Number of individuals", x$cluster_diagnostics$n_individuals))
    cat(sprintf("%-*s %s\n", wide_label_width, "Number of clusters", x$cluster_diagnostics$n_clusters))
    cat(sprintf("%-*s %s\n", wide_label_width, "Cluster outcome", .ri_cluster_outcome_label(x$cluster_outcome)))
    cat(sprintf("%-*s %s\n", wide_label_width, "Cluster missingness", .ri_cluster_missing_label(x$cluster_missing)))
    cat(sprintf("%-*s %s\n", wide_label_width, "Mean cluster size", .ri_format_num(x$cluster_diagnostics$mean_cluster_size)))
    cat(sprintf("%-*s %s\n", wide_label_width, "Mean cluster response rate", .ri_format_num(x$cluster_diagnostics$mean_cluster_response_rate)))
    cat(sprintf("%-*s %s\n", wide_label_width, "Clusters with no observed outcome", x$cluster_diagnostics$clusters_with_no_observed_outcome))
    cat(sprintf("%-*s %s\n\n", wide_label_width, "Clusters with any missing", x$cluster_diagnostics$clusters_with_any_missing))
  }

  if (isTRUE(x$block_missing_sum) && !is.null(x$block_diagnostics)) {
    cat("Block missingness\n", sep = "")
    cat(sprintf("%-*s %s\n", label_width, "Number of blocks", x$block_diagnostics$n_blocks))
    if (isTRUE(x$block_missing_sum)) {
      cat(sprintf("%-*s\n", label_width, "Within-block attrition diff (T - C)"))
      cat(sprintf("%-*s %-*s %-*s %-*s %-*s %-*s %-*s\n",
                  label_width, "",
                  small_col_width, "Mean",
                  small_col_width, "Min",
                  small_col_width, "Q1",
                  small_col_width, "Median",
                  small_col_width, "Q3",
                  small_col_width, "Max"))
      cat(sprintf("%-*s %-*s %-*s %-*s %-*s %-*s %-*s\n",
                  label_width, "",
                  small_col_width, .ri_format_num(x$block_diagnostics$mean_within_block_attrition_diff),
                  small_col_width, .ri_format_num(x$block_diagnostics$min_within_block_attrition_diff),
                  small_col_width, .ri_format_num(x$block_diagnostics$q1_within_block_attrition_diff),
                  small_col_width, .ri_format_num(x$block_diagnostics$median_within_block_attrition_diff),
                  small_col_width, .ri_format_num(x$block_diagnostics$q3_within_block_attrition_diff),
                  small_col_width, .ri_format_num(x$block_diagnostics$max_within_block_attrition_diff)))
    }
    cat("\n", sep = "")
  }

  if (isTRUE(x$include_ci)) {
    cat(strrep("=", 120), "\n", sep = "")
    cat(sprintf("%-22s %-12s %-22s %-16s %-12s %-28s\n",
                names(summary.table)[1], names(summary.table)[2],
                names(summary.table)[3], names(summary.table)[4],
                names(summary.table)[5], names(summary.table)[6]))
    cat(strrep("=", 120), "\n", sep = "")
    for (i in seq_len(nrow(summary.table))) {
      cat(sprintf("%-22s %-12s %-22s %-16s %-12s %-28s\n",
                  summary.table$Assumption[i],
                  summary.table$Null[i],
                  summary.table$`Test stat.`[i],
                  summary.table$Value[i],
                  summary.table$`P-value`[i],
                  summary.table[[6]][i]))
    }
    cat(strrep("=", 120), "\n", sep = "")
  } else {
    cat(strrep("=", 90), "\n", sep = "")
    cat(sprintf("%-22s %-12s %-22s %-16s %-12s\n",
                names(summary.table)[1], names(summary.table)[2],
                names(summary.table)[3], names(summary.table)[4],
                names(summary.table)[5]))
    cat(strrep("=", 90), "\n", sep = "")
    for (i in seq_len(nrow(summary.table))) {
      cat(sprintf("%-22s %-12s %-22s %-16s %-12s\n",
                  summary.table$Assumption[i],
                  summary.table$Null[i],
                  summary.table$`Test stat.`[i],
                  summary.table$Value[i],
                  summary.table$`P-value`[i]))
    }
    cat(strrep("=", 90), "\n", sep = "")
  }

  if (length(x$notes) > 0) {
    cat("\nNotes\n", sep = "")
    for (note in x$notes) {
      cat("- ", note, "\n", sep = "")
    }
  }

  cat("\nAssumption note\n", sep = "")
  cat(x$assumption_hint, "\n", sep = "")

  invisible(x)
}

#' @export
summary.riattrition_result <- function(object, ...) {
  class(object) <- c("summary.riattrition_result", class(object))
  object
}

#' @export
print.summary.riattrition_result <- function(x, ...) {
  print.riattrition_result(x, ...)
}
