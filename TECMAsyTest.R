TECMAsyTest <- function(model, robust = FALSE) {
  # Get all coefficient names
  all_coefs <- names(coef(model))

  # Initialize results storage
  adjustment_table <- NULL
  short_term_table <- NULL

  # 1. Adjustment Path Asymmetry (Long-run) Test
  if ("ect_pos" %in% all_coefs && "ect_neg" %in% all_coefs) {
    hypothesis_adj <- "ect_pos = ect_neg"

    if (robust) {
      vcov_mat <- sandwich::vcovHC(model)
    } else {
      vcov_mat <- vcov(model)
    }

    test_result_adj <- car::linearHypothesis(
      model = model,
      hypothesis.matrix = hypothesis_adj,
      vcov. = vcov_mat,
      test = "Chisq"
    )

    adjustment_table <- tibble::tibble(
      Test = "Adjustment Path",
      Hypothesis = hypothesis_adj,
      `Chi-Square` = test_result_adj$Chisq[2],
      df = test_result_adj$Df[2],
      `P-value` = test_result_adj$"Pr(>Chisq)"[2],
      Robust = robust
    )
  }

  # 2. Short-term Asymmetry Tests
  lag_coefs <- all_coefs[grepl("^L\\(", all_coefs)]
  base_vars <- unique(
    gsub("^L\\((.*?)_(pos|neg),\\s*\\d+\\)$", "\\1", lag_coefs)
  )
  base_vars <- base_vars[!base_vars %in% c("", "ect")]

  short_term_rows <- list()

  for (var in base_vars) {
    var_coefs <- all_coefs[grepl(paste0("^L\\(", var, "_"), all_coefs)]
    pos_coefs <- var_coefs[grepl("_pos", var_coefs)]
    neg_coefs <- var_coefs[grepl("_neg", var_coefs)]

    if (length(pos_coefs) == 0 || length(neg_coefs) == 0) next

    # Individual lag tests
    lag_map <- unique(
      as.numeric(gsub(".*,\\s*(\\d+)\\)$", "\\1", var_coefs))
    )

    for (lg in lag_map) {
      lg_coefs <- var_coefs[grepl(paste0(",\\s*", lg, "\\)$"), var_coefs)]
      pos_name <- lg_coefs[grepl("_pos", lg_coefs)]
      neg_name <- lg_coefs[grepl("_neg", lg_coefs)]

      if (length(pos_name) == 1 && length(neg_name) == 1) {
        hypothesis <- paste(pos_name, "=", neg_name)

        if (robust) {
          vcov_mat <- sandwich::vcovHC(model)
        } else {
          vcov_mat <- vcov(model)
        }

        test_result <- car::linearHypothesis(
          model = model,
          hypothesis.matrix = hypothesis,
          vcov. = vcov_mat,
          test = "Chisq"
        )

        short_term_rows[[paste(var, lg)]] <- tibble::tibble(
          Test = "Short-term",
          Variable = var,
          Lag = as.character(lg),  # Convert to character
          Hypothesis = hypothesis,
          `Chi-Square` = test_result$Chisq[2],
          df = test_result$Df[2],
          `P-value` = test_result$"Pr(>Chisq)"[2],
          Robust = robust
        )
      }
    }

    # Joint test
    if (length(pos_coefs) == length(neg_coefs)) {
      joint_hyp <- paste0(var, ":Cumulative asymetry effect")

      if (robust) {
        vcov_mat <- sandwich::vcovHC(model)
      } else {
        vcov_mat <- vcov(model)
      }

      test_result_joint <- car::linearHypothesis(
        model = model,
        hypothesis.matrix = paste(pos_coefs, "=", neg_coefs),
        vcov. = vcov_mat,
        test = "Chisq"
      )

      short_term_rows[[paste(var, "joint")]] <- tibble::tibble(
        Test = "Short-term",
        Variable = var,
        Lag = "Joint",
        Hypothesis = joint_hyp,
        `Chi-Square` = test_result_joint$Chisq[2],
        df = test_result_joint$Df[2],
        `P-value` = test_result_joint$"Pr(>Chisq)"[2],
        Robust = robust
      )
    }
  }

  # Combine all results
  if (length(short_term_rows) > 0) {
    short_term_table <- dplyr::bind_rows(short_term_rows)
  }

  # Create final output
  result_list <- list(
    adjustment_test = adjustment_table,
    short_term_tests = short_term_table,
    robust = robust,
    tested_at = Sys.time(),
    model_formula = formula(model)
  )

  # Custom print method
  class(result_list) <- "TECMAsyTest"
  return(result_list)
}

#Final output
print.TECMAsyTest <- function(x, digits = 3, ...) {
  cat("\nECM Asymmetry Test Results\n")
  cat("Model:", deparse(x$model_formula), "\n")
  cat("Tested at:", format(x$tested_at), "\n")
  cat("Robust standard errors:", x$robust, "\n\n")

  if (!is.null(x$adjustment_test)) {
    cat("--- Adjustment Path Asymmetry Test ---\n")
    print(knitr::kable(x$adjustment_test, digits = digits))
    cat("\n")
  }

  if (!is.null(x$short_term_tests)) {
    cat("--- Short-term Asymmetry Tests ---\n")
    print(knitr::kable(x$short_term_tests, digits = digits))
    cat("\n")
  }

  invisible(x)
}

