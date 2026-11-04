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
#' The printed output is a compact plain-text report: a header (outcome,
#' design, test, sharp null), outcome missingness overall and by arm, one row of
#' Monte Carlo p-value and confidence interval per missingness assumption, and
#' notes explaining the assumptions shown.
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
#' @param block_missing_sum A logical value indicating whether the printed
#'   output should also show the distribution (median, IQR, range) of the
#'   within-block treatment - control missingness differences.
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
#' @param nperm A positive integer that specifies the number of Monte Carlo
#'   (MC) replications, i.e. random re-assignments drawn to approximate the
#'   randomization distribution (default 10^4).
#' @param alpha A scalar, where \eqn{1-\alpha} indicates the confidence level.
#' @param tol A numerical object that specifies the precision of confidence
#'   intervals.
#' @param include_ci A logical value indicating whether confidence intervals
#'   should be computed for the printed result.
#' @param include_twostep A logical value indicating whether two-step rows
#'   should be added for `general`, `mp`, and `mn` (default `FALSE`).
#' @param beta A real number that belongs to \eqn{[0, \alpha]} for the two-step
#'   procedure.
#' @param outcome_name Name of the outcome shown in the printed output. By
#'   default it is taken from the `Y` argument, e.g. `data$y` gives `y`.
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
                               include_ci = TRUE, include_twostep = FALSE,
                               beta = 0.1 * alpha, outcome_name = NULL) {
  # Name for the "Outcome:" line of the printed output: by default the part of
  # the Y argument after the last "$", e.g. data$discernment -> discernment.
  if (is.null(outcome_name)) {
    outcome_name <- paste(deparse(substitute(Y)), collapse = "")
    outcome_name <- sub(".*[$]", "", outcome_name)
  }

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

  if (length(assumptions_reported) > 1 && (!is.null(Z.perm) || !is.null(stat.null))) {
    notes <- c(
      notes,
      "Custom stat.null/Z.perm inputs were not reused across multiple assumption rows; each row was recomputed with its own design-compatible permutation inputs."
    )
  }

  results_list <- vector("list", length = 0)
  mc_replications <- integer(0)

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

    # Number of Monte Carlo replications actually behind this row: nperm unless
    # the user supplied their own stat.null / Z.perm.
    mc_replications <- c(
      mc_replications,
      if (!is.null(stat.null_i)) length(stat.null_i) else ncol(Z.perm_i)
    )

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
      outcome_name = outcome_name,
      assumption = missing,
      assumptions_reported = assumptions_reported,
      test_class = class,
      test_name = test_name,
      null_hypothesis = c,
      mc_replications = unique(mc_replications),
      counts = counts,
      design = analysis_data$design,
      inference_unit = analysis_data$inference_unit,
      cluster_outcome = analysis_data$cluster_outcome,
      cluster_missing = analysis_data$cluster_missing,
      cluster_diagnostics = analysis_data$cluster_diagnostics,
      block_diagnostics = block_diagnostics,
      block_missing_sum = block_missing_sum,
      assumption_hint = .ri_assumption_hint(counts, block_diagnostics),
      results = results,
      confidence_level = 1 - alpha,
      beta = beta,
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
#
# In blocked designs the arms are compared within blocks (the block-size
# weighted difference from .ri_block_diagnostics()).  When the treated share
# differs across blocks, the pooled difference mixes block and arm and can even
# have the opposite sign.
.ri_attrition_comparison <- function(counts, block_diagnostics = NULL, tol = 1e-8) {
  pooled_diff <- counts$attrition_rate_treat - counts$attrition_rate_control
  within_diff <- if (is.null(block_diagnostics)) {
    NA_real_
  } else {
    block_diagnostics$weighted_mean_within_block_attrition_diff
  }
  use_blocks <- !is.na(within_diff)
  diff <- if (use_blocks) within_diff else pooled_diff
  list(
    pooled_diff = pooled_diff,
    within_diff = within_diff,
    use_blocks = use_blocks,
    diff = diff,
    reversed = use_blocks && abs(pooled_diff) > tol && abs(diff) > tol &&
      sign(pooled_diff) != sign(diff)
  )
}

.ri_assumption_hint <- function(counts, block_diagnostics = NULL, tol = 1e-8) {
  cmp <- .ri_attrition_comparison(counts, block_diagnostics, tol)
  use_blocks <- cmp$use_blocks
  diff <- cmp$diff
  subject <- if (use_blocks) "Within blocks (block-size weighted), treated" else "Treated"

  if (abs(diff) <= tol) {
    hint <- if (use_blocks) {
      "Within blocks (block-size weighted), treated and control attrition rates are very similar; rates alone do not favor mp or mn."
    } else {
      "Observed treated and control attrition rates are very similar; rates alone do not favor mp or mn."
    }
  } else if (diff < 0) {
    hint <- paste(subject, "attrition is lower than control attrition; this is more consistent with mp than mn, although rates alone do not identify the missingness mechanism.")
  } else {
    hint <- paste(subject, "attrition is higher than control attrition; this is more consistent with mn than mp, although rates alone do not identify the missingness mechanism.")
  }

  if (cmp$reversed) {
    hint <- paste(hint, .ri_reversal_sentence)
  }

  hint
}

.ri_reversal_sentence <- "The pooled attrition rates point the other way because they also reflect differences between blocks (for example, different treated shares); the within-block comparison is the relevant one."

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

  # Blocks without both arms have no within-block difference and are dropped.
  valid <- !is.na(block_table$attrition_diff)
  valid_diff <- block_table$attrition_diff[valid]
  valid_n <- block_table$n[valid]
  q <- if (length(valid_diff) > 0) {
    stats::quantile(valid_diff, probs = c(0.25, 0.5, 0.75), names = FALSE)
  } else {
    c(NA_real_, NA_real_, NA_real_)
  }

  list(
    n_blocks = length(blocks),
    n_blocks_used = sum(valid),
    # The mean weights each block by its share of units (n_b / N, over the blocks
    # used), so it estimates the average effect of treatment on attrition for the
    # whole sample rather than for the average block. The quantiles stay
    # unweighted: they describe how the block-level differences are spread out.
    weighted_mean_within_block_attrition_diff = if (length(valid_diff) > 0) stats::weighted.mean(valid_diff, valid_n) else NA_real_,
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

  # Under the sharp null tau*_i = c for all i, these are the imputed potential outcomes:
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

# Whole numbers such as MC replication counts: "10,000" rather than "1e+04".
.ri_format_count <- function(x) {
  paste(format(x, big.mark = ",", scientific = FALSE, trim = TRUE), collapse = " / ")
}

# Rates stored as proportions are shown as percentages: 0.25 -> "25.0%".
.ri_format_pct <- function(x, digits = 1) {
  if (is.na(x)) {
    return("NA")
  }
  paste0(formatC(100 * x, digits = digits, format = "f"), "%")
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

#----------- Printed output ---------------------------------------------------#
#
# print.riattrition_result() writes a compact plain-text report in four parts:
#   1. a header: outcome, design, test, and sharp null;
#   2. outcome missingness: overall and by arm, the treatment - control
#      difference, and for blocked designs the within-block differences;
#   3. the inference table: one row per missingness assumption;
#   4. notes explaining what is in the table.
# The notes for "general" and "mn" and the note on p-values and confidence
# intervals follow the wording reviewed by the coauthors. The notes for the
# other cases (mp, sharp, random, two-step rows, c != 0, one-sided intervals)
# are drafts written in the same style.

# Names of the missingness assumptions in the inference table.
.ri_assumption_labels <- c(
  general = "Unrestricted (worst case)",
  mp = "Monotone positive, M(1) \u2265 M(0)",
  mn = "Monotone negative, M(1) \u2264 M(0)",
  sharp = "No effect on missingness, M(1) = M(0)",
  random = "Missing at random"
)

# A section title with a line of "-" (or "=") of the same length underneath.
.ri_print_title <- function(title, underline = "-") {
  cat(title, "\n", strrep(underline, nchar(title)), "\n", sep = "")
}

# Add spaces so that `text` takes up `width` characters.
.ri_pad_right <- function(text, width) {
  paste0(text, strrep(" ", max(0, width - nchar(text))))
}

.ri_pad_left <- function(text, width) {
  paste0(strrep(" ", max(0, width - nchar(text))), text)
}

# A number with a real minus sign (Unicode U+2212) instead of a hyphen: -0.5.
.ri_signed_num <- function(x, digits = 3) {
  text <- formatC(abs(x), digits = digits, format = "f")
  if (x < 0) {
    text <- paste0("\u2212", text)
  }
  text
}

# A proportion as a percentage with one decimal: -0.127 -> minus sign, "12.7%".
# With plus = TRUE a positive value gets a "+" sign: 0.011 -> "+1.1%".
.ri_pct <- function(x, plus = FALSE) {
  if (is.na(x)) {
    return("NA")
  }
  value <- 100 * x
  text <- formatC(abs(value), digits = 1, format = "f")
  if (value <= -0.05) {
    text <- paste0("\u2212", text)
  }
  if (plus && value >= 0.05) {
    text <- paste0("+", text)
  }
  paste0(text, "%")
}

# Monte Carlo p-value with three decimals; very small values become "<0.001".
.ri_pvalue <- function(p) {
  if (is.na(p)) {
    return("NA")
  }
  if (p < 0.001) {
    return("<0.001")
  }
  formatC(p, digits = 3, format = "f")
}

# Confidence interval: square brackets for finite ends, round brackets and the
# infinity sign for infinite ends, e.g. "[0.011, 0.066]" or "(-inf, inf)".
.ri_ci_text <- function(lower, upper) {
  if (is.na(lower) || is.na(upper)) {
    return("NA")
  }
  if (is.infinite(lower)) {
    left <- if (lower < 0) "(\u2212\u221e" else "(\u221e"
  } else {
    left <- paste0("[", .ri_signed_num(lower))
  }
  if (is.infinite(upper)) {
    right <- if (upper > 0) "\u221e)" else "\u2212\u221e)"
  } else {
    right <- paste0(.ri_signed_num(upper), "]")
  }
  paste0(left, ", ", right)
}

# The value of the sharp null, e.g. "0" or minus sign + "2".
.ri_null_value <- function(c) {
  if (c < 0) {
    return(paste0("\u2212", format(abs(c))))
  }
  format(c)
}

# Wrap a note to lines of at most `width` characters without breaking a short
# formula such as "M(1) <= M(0)": the spaces around =, <=, >= and the minus sign
# are made non-breaking while wrapping and turned back into spaces afterwards.
.ri_wrap <- function(text, width = 90) {
  keep_together <- gsub(" (=|\u2264|\u2265|\u2212) ", "\u00a0\\1\u00a0", text)
  lines <- strwrap(keep_together, width = width)
  paste(gsub("\u00a0", " ", lines), collapse = "\n")
}

# Note for one non-general missingness assumption shown in the table.
.ri_assumption_note <- function(assumption, x, tol = 1e-8) {
  if (assumption == "random") {
    return(paste(
      "Missing at random means that whether an outcome is missing is unrelated to the",
      "potential outcomes; the test then uses the observed outcomes only."
    ))
  }
  if (assumption == "sharp") {
    return(paste(
      "M(z) = 1 indicates an observed outcome. The assumption M(1) = M(0) means that",
      "treatment does not change whether any outcome is observed; the test then uses the",
      "observed outcomes only."
    ))
  }

  # mp or mn: describe the assumption, then say whether the observed rates
  # point the same way. With blocks the within-block (block-weighted)
  # difference is used, otherwise the unadjusted one.
  if (assumption == "mn") {
    text <- paste(
      "M(z) = 1 indicates an observed outcome. The assumption M(1) \u2264 M(0) means that",
      "every outcome observed under treatment would also be observed under control."
    )
  } else {
    text <- paste(
      "M(z) = 1 indicates an observed outcome. The assumption M(1) \u2265 M(0) means that",
      "every outcome observed under control would also be observed under treatment."
    )
  }

  difference <- .ri_attrition_comparison(x$counts, x$block_diagnostics, tol)$diff
  if (assumption == "mn") {
    rates_agree <- difference >= -tol
  } else {
    rates_agree <- difference <= tol
  }

  if (rates_agree) {
    text <- paste(text, "The observed rates are compatible with this assumption but cannot verify this assumption.")
  } else {
    text <- paste(text, "The observed rates point the other way, which casts doubt on this assumption.")
  }
  text
}

# All notes printed under the inference table, as a vector of paragraphs.
.ri_notes_text <- function(x) {
  assumptions <- unique(x$results$`Missingness assumption`)
  notes <- character(0)

  # 1. The unrestricted (worst-case) row.
  if ("general" %in% assumptions) {
    notes <- c(notes, paste(
      "Unrestricted missingness uses least-favorable imputation for outcomes missing.",
      "It remains valid under arbitrary attrition, but may be highly conservative and",
      "produce an uninformative confidence interval."
    ))
  }

  # 2. Every other assumption in the table.
  for (assumption in assumptions) {
    if (assumption != "general") {
      notes <- c(notes, .ri_assumption_note(assumption, x))
    }
  }

  # 3. Two-step rows, if any.
  if (any(x$results$Method == "Two-step")) {
    notes <- c(notes, paste0(
      "Two-step rows first use the data to bound how many units can be of each ",
      "missingness type, at error level \u03b2 = ", format(x$beta), ", and then use ",
      "least-favorable imputation only within that bound; \u03b2 is added to the p-value."
    ))
  }

  # 4. What the p-values and confidence intervals refer to.
  if (isTRUE(all.equal(x$null_hypothesis, 0))) {
    text <- "P-values test the sharp null that every unit has zero treatment effect."
  } else {
    text <- paste0(
      "P-values test the sharp null that every unit has a treatment effect of exactly ",
      .ri_null_value(x$null_hypothesis), "."
    )
  }
  if (isTRUE(x$include_ci)) {
    level <- round(100 * x$confidence_level)
    if (identical(x$alternative, "two.sided")) {
      interval <- "confidence intervals"
    } else {
      interval <- "one-sided confidence intervals"
    }
    text <- paste0(
      text, " The ", level, "% ", interval, " are obtained by inverting tests of ",
      "constant additive effects, Y_i(1) \u2212 Y_i(0) = \u03c4."
    )
  }
  notes <- c(notes, text)

  # 5. Messages from the computation itself (for example about custom inputs).
  c(notes, x$notes)
}

#----------- S3 print/summary methods ----------------------------------------#
#' @export
print.riattrition_result <- function(x, ...) {
  counts <- x$counts
  label_width <- 24

  # 1. Header -----------------------------------------------------------------
  .ri_print_title("Randomization inference with attrition", underline = "=")
  cat("Outcome: ", x$outcome_name, "\n", sep = "")
  # e.g. "individual-level blocked randomization"; cluster designs already
  # say "cluster randomization".
  if (identical(x$inference_unit, "individual")) {
    cat("Design: individual-level ", x$design, "\n", sep = "")
  } else {
    cat("Design: ", x$design, "\n", sep = "")
  }
  cat("Test: ", x$test_name, " with ", .ri_format_count(x$mc_replications),
      " Monte Carlo permutations\n", sep = "")
  cat("Sharp null: Y_i(1) \u2212 Y_i(0) = ", .ri_null_value(x$null_hypothesis),
      " for every unit\n\n", sep = "")

  # 2. Outcome missingness ------------------------------------------------------
  .ri_print_title("Outcome missingness")
  cat(.ri_pad_right("", label_width), .ri_pad_left("Missing / N", 15), .ri_pad_left("Rate", 9), "\n", sep = "")
  print_missing_row <- function(label, n_missing, n, rate) {
    missing_of_n <- paste0(.ri_format_count(n_missing), " / ", .ri_format_count(n))
    cat(.ri_pad_right(label, label_width), .ri_pad_left(missing_of_n, 15),
        .ri_pad_left(.ri_pct(rate), 9), "\n", sep = "")
  }
  print_missing_row("Overall", counts$n_miss, counts$n_total, counts$attrition_rate_total)
  print_missing_row("Treatment", counts$n_miss_treat, counts$n_treat, counts$attrition_rate_treat)
  print_missing_row("Control", counts$n_miss_control, counts$n_control, counts$attrition_rate_control)

  cat("\nTreatment \u2212 control missingness difference:\n")
  unadjusted <- counts$attrition_rate_treat - counts$attrition_rate_control
  cat(.ri_pad_right("Unadjusted:", label_width), .ri_pad_left(.ri_pct(unadjusted, plus = TRUE), 24), "\n", sep = "")
  blocks <- x$block_diagnostics
  if (!is.null(blocks)) {
    cat(.ri_pad_right("Block-weighted:", label_width),
        .ri_pad_left(.ri_pct(blocks$weighted_mean_within_block_attrition_diff, plus = TRUE), 24), "\n", sep = "")
    cat("\nBlocks: ", blocks$n_blocks, " total; ", blocks$n_blocks_used,
        " contain both treatment arms.\n", sep = "")
    if (isTRUE(x$block_missing_sum)) {
      cat("\nWithin-block treatment \u2212 control differences:\n")
      iqr <- paste0("[", .ri_pct(blocks$q1_within_block_attrition_diff), ", ",
                    .ri_pct(blocks$q3_within_block_attrition_diff), "]")
      range <- paste0("[", .ri_pct(blocks$min_within_block_attrition_diff), ", ",
                      .ri_pct(blocks$max_within_block_attrition_diff), "]")
      cat(.ri_pad_right("Median:", label_width), .ri_pad_left(.ri_pct(blocks$median_within_block_attrition_diff), 24), "\n", sep = "")
      cat(.ri_pad_right("IQR:", label_width), .ri_pad_left(iqr, 24), "\n", sep = "")
      cat(.ri_pad_right("Range:", label_width), .ri_pad_left(range, 24), "\n", sep = "")
    }
  }

  # Cluster designs are not redesigned yet: their summary is printed as before.
  if (!is.null(x$cluster_diagnostics)) {
    wide <- 38
    cat("\nCluster summary\n")
    cat(.ri_pad_right("Number of individuals", wide), x$cluster_diagnostics$n_individuals, "\n", sep = "")
    cat(.ri_pad_right("Number of clusters", wide), x$cluster_diagnostics$n_clusters, "\n", sep = "")
    cat(.ri_pad_right("Cluster outcome", wide), .ri_cluster_outcome_label(x$cluster_outcome), "\n", sep = "")
    cat(.ri_pad_right("Cluster missingness", wide), .ri_cluster_missing_label(x$cluster_missing), "\n", sep = "")
    cat(.ri_pad_right("Mean cluster size", wide), .ri_format_num(x$cluster_diagnostics$mean_cluster_size), "\n", sep = "")
    cat(.ri_pad_right("Mean cluster response rate", wide), .ri_pct(x$cluster_diagnostics$mean_cluster_response_rate), "\n", sep = "")
    cat(.ri_pad_right("Clusters with no observed outcome", wide), x$cluster_diagnostics$clusters_with_no_observed_outcome, "\n", sep = "")
    cat(.ri_pad_right("Clusters with any missing", wide), x$cluster_diagnostics$clusters_with_any_missing, "\n", sep = "")
  }

  # 3. Inference table ------------------------------------------------------------
  name_width <- 40
  cat("\n")
  .ri_print_title("Inference")
  header <- paste0(.ri_pad_right("Missingness assumption", name_width), .ri_pad_left("MC p-value", 10))
  if (isTRUE(x$include_ci)) {
    level <- round(100 * x$confidence_level)
    ci_label <- paste0(level, "% CI for constant effect")
    if (!identical(x$alternative, "two.sided")) {
      ci_label <- paste0(level, "% one-sided CI for constant effect")
    }
    header <- paste0(header, "    ", ci_label)
  }
  cat(header, "\n", sep = "")

  for (i in seq_len(nrow(x$results))) {
    row <- x$results[i, ]
    name <- .ri_assumption_labels[[row$`Missingness assumption`]]
    if (row$Method == "Two-step") {
      name <- paste0(name, ", two-step")
    }
    line <- paste0(.ri_pad_right(name, name_width), .ri_pad_left(.ri_pvalue(row$`P-value`), 10))
    if (isTRUE(x$include_ci)) {
      line <- paste0(line, "    ", .ri_pad_left(.ri_ci_text(row$`CI lower`, row$`CI upper`), nchar(ci_label)))
    }
    cat(line, "\n", sep = "")
  }

  # 4. Notes ------------------------------------------------------------------------
  cat("\n")
  .ri_print_title("Notes")
  notes <- .ri_notes_text(x)
  for (i in seq_along(notes)) {
    if (i > 1) {
      cat("\n")
    }
    cat(.ri_wrap(notes[i]), "\n", sep = "")
  }

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
