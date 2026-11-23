#----------- Test of monotone missingness ------------------------------------#
#
# Monotone missingness (Li, Sheng & Yu 2025, Assumptions 2 and 3) is an
# assumption about how the treatment affects whether an outcome is observed:
#   * mp: M_i(1) >= M_i(0) for every unit (the treatment never makes an
#         outcome missing that would otherwise be observed);
#   * mn: M_i(1) <= M_i(0) for every unit.
#
# M_i is observed for every unit, so the effect of Z on M can be tested like
# the effect on any other outcome. mp says that every unit-level effect
# M_i(1) - M_i(0) is >= 0: a bounded null hypothesis. A randomization test of
# "no effect on M" that uses an effect-increasing statistic (here the
# difference in response rates) is valid for this bounded null (Caughey,
# Dafoe, Li & Miratrix 2023). A small p-value is evidence that the treatment
# made some outcomes missing (against mp) or observed (against mn).
#
# A large p-value does not show that the assumption holds: the test only sees
# the net difference in response rates, so units that violate the assumption
# can be hidden by units with an effect in the other direction.


# Block-size weighted difference in response rates (treatment - control),
# using only blocks that contain both arms. `treated_observed` is the number of
# treated units with an observed outcome in each block; it can be a vector (one
# assignment) or a matrix with one row per Monte Carlo draw.
.ri_weighted_response_diff <- function(treated_observed, blocks) {
  w <- blocks$n / sum(blocks$n)
  slope <- w * (1 / blocks$n_treat + 1 / blocks$n_control)
  # Treated rate - control rate in block b equals
  #   T_b / n1_b - (O_b - T_b) / n0_b = T_b * (1 / n1_b + 1 / n0_b) - O_b / n0_b.
  as.numeric(treated_observed %*% slope) - sum(w * blocks$n_observed / blocks$n_control)
}

#' Test of monotone missingness
#'
#' Tests the monotone missingness assumption used by `missing = "mp"` or
#' `missing = "mn"`: H0: M_i(1) >= M_i(0) for every unit (`"mp"`) or
#' H0: M_i(1) <= M_i(0) for every unit (`"mn"`), where M_i = 1 if unit i's
#' outcome is observed.
#'
#' @param Z Treatment assignment (0/1 vector).
#' @param Y Outcome vector; `NA` marks a missing outcome.
#' @param assumption `"mp"` or `"mn"`.
#' @param block Optional block identifier vector. If supplied, treatment is
#'   reshuffled within blocks.
#' @param nperm Number of Monte Carlo replications (only used with `block`).
#'
#' @return An object of class `monotone_missingness_test`.
#' @export
test_monotone_missingness <- function(Z, Y, assumption = c("mp", "mn"), block = NULL,
                                      nperm = 10^4) {
  assumption <- match.arg(assumption)
  if (length(Z) != length(Y)) {
    stop("`Z` and `Y` must have the same length.")
  }
  if (!all(Z %in% c(0, 1))) {
    stop("`Z` must be a 0/1 vector.")
  }
  M <- as.numeric(!is.na(Y))
  rate_treat <- mean(M[Z == 1])
  rate_control <- mean(M[Z == 0])

  if (is.null(block)) {
    # Complete randomization. Under "no effect on M", the number of treated
    # units with an observed outcome is hypergeometric, so the p-value is exact
    # (a one-sided Fisher exact test).
    treated_observed <- sum(M[Z == 1])
    n_observed <- sum(M)
    n_missing <- sum(M == 0)
    n_treat <- sum(Z)
    if (assumption == "mp") {
      # Evidence against mp: few treated units with an observed outcome.
      p_value <- stats::phyper(treated_observed, n_observed, n_missing, n_treat)
    } else {
      # Evidence against mn: many treated units with an observed outcome.
      p_value <- stats::phyper(treated_observed - 1, n_observed, n_missing, n_treat,
                               lower.tail = FALSE)
    }
    design <- "complete randomization"
    method <- "exact"
    statistic <- rate_treat - rate_control
    mc_replications <- NA
  } else {
    # Blocked randomization. Summarize each block; blocks with only one arm
    # carry no information about the effect of Z and are left out.
    block <- as.character(block)
    blocks <- data.frame(
      n = as.numeric(table(block)),
      n_treat = as.numeric(tapply(Z, block, sum)),
      n_observed = as.numeric(tapply(M, block, sum)),
      treated_observed = as.numeric(tapply(M * Z, block, sum))
    )
    blocks$n_control <- blocks$n - blocks$n_treat
    blocks <- blocks[blocks$n_treat > 0 & blocks$n_control > 0, , drop = FALSE]

    statistic <- .ri_weighted_response_diff(blocks$treated_observed, blocks)

    # Under "no effect on M", the number of treated units with an observed
    # outcome in each block is hypergeometric, independently across blocks.
    # Drawing these counts is the same as reshuffling treatment within blocks.
    draws <- sapply(seq_len(nrow(blocks)), function(b) {
      stats::rhyper(nperm, blocks$n_observed[b], blocks$n[b] - blocks$n_observed[b],
                    blocks$n_treat[b])
    })
    draws <- matrix(draws, nrow = nperm)
    stat_null <- .ri_weighted_response_diff(draws, blocks)

    # A small tolerance keeps ties (equal statistics up to rounding) as ties.
    tol <- 1e-12
    p_value <- if (assumption == "mp") {
      mean(stat_null <= statistic + tol)
    } else {
      mean(stat_null >= statistic - tol)
    }
    design <- "blocked randomization"
    method <- "Monte Carlo"
    mc_replications <- nperm
  }

  structure(
    list(
      assumption = assumption,
      hypothesis = if (assumption == "mp") {
        "M(1) \u2265 M(0) for every unit"
      } else {
        "M(1) \u2264 M(0) for every unit"
      },
      design = design,
      method = method,
      mc_replications = mc_replications,
      response_rate_treat = rate_treat,
      response_rate_control = rate_control,
      statistic = statistic,
      p_value = p_value
    ),
    class = "monotone_missingness_test"
  )
}

#' @export
print.monotone_missingness_test <- function(x, ...) {
  diff_label <- if (identical(x$design, "blocked randomization")) "block-weighted difference" else "difference"
  p_text <- if (x$p_value < 0.001) "p < 0.001" else paste0("p = ", .ri_pvalue(x$p_value))
  method_text <- if (identical(x$method, "exact")) {
    "exact test, complete randomization"
  } else {
    paste0("Monte Carlo, ", .ri_format_count(x$mc_replications), " replications within blocks")
  }

  cat("Test of monotone missingness (", x$assumption, ")\n", sep = "")
  cat("H0: ", x$hypothesis, ", where M = 1 means the outcome is observed\n", sep = "")
  cat("Response rate: treatment ", .ri_pct(x$response_rate_treat), ", control ",
      .ri_pct(x$response_rate_control), " (", diff_label, " ",
      .ri_pct(x$statistic, plus = TRUE), ")\n", sep = "")
  cat(p_text, " (", method_text, ")\n\n", sep = "")
  cat("A small p-value is evidence against the assumption. A large p-value does not show\n")
  cat("that the assumption holds: units that violate it can be offset by units with an\n")
  cat("effect in the other direction.\n")
  invisible(x)
}

# Runs `expr` without changing the random-number state, so that adding the
# monotonicity test to ri_sharp_attrition() leaves all other results unchanged.
.ri_keep_seed <- function(expr) {
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (has_seed) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
  }
  expr
}
