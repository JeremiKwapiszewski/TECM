TECMLag <- function(y, X, lag.max = 5,
                    model = c("linear", "tar", "mtar"),
                    thresh = 0,
                    split = TRUE, break.start = NULL, break.end = NULL) {

  model <- match.arg(model)
  lag_range <- 0:lag.max
  aic_vals <- numeric(length(lag_range))
  bic_vals <- numeric(length(lag_range))
  models_list <- vector("list", length(lag_range))

  for (i in seq_along(lag_range)) {
    lag_i <- lag_range[i]

    fit <- tryCatch({
      TECMFit(y, X, lag = lag_i, model = model, thresh = thresh,
              split = split, break.start = break.start, break.end = break.end)
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      models_list[[i]] <- fit
      aic_vals[i] <- AIC(fit$ECM_model)
      bic_vals[i] <- BIC(fit$ECM_model)
    } else {
      aic_vals[i] <- NA
      bic_vals[i] <- NA
    }
  }

  info_criteria <- tibble::tibble(
    Lag = lag_range,
    AIC = round(aic_vals, 3),
    BIC = round(bic_vals, 3)
  )

  best_aic_lag <- info_criteria$Lag[which.min(info_criteria$AIC)]
  best_bic_lag <- info_criteria$Lag[which.min(info_criteria$BIC)]

  result <- list(
    table = info_criteria,
    best_lag_aic = best_aic_lag,
    best_lag_bic = best_bic_lag,
    best_model_aic = models_list[[which.min(aic_vals)]],
    best_model_bic = models_list[[which.min(bic_vals)]],
    model_type = model,
    threshold = thresh,
    tested_at = Sys.time()
  )

  class(result) <- "TECMLag"
  return(result)
}

print.TECMLag <- function(x, digits = 3, ...) {
  cat("\nOptimal Lag Selection for TECM Model\n")
  cat("Model Type:", x$model_type, "\n")
  cat("Threshold Value:", x$threshold, "\n")
  cat("Tested at:", format(x$tested_at), "\n\n")

  cat("--- Information Criteria Table ---\n")
  print(knitr::kable(x$table, digits = digits))
  cat("\n")

  cat("Selected Lag (AIC):", x$best_lag_aic, "\n")
  cat("Selected Lag (BIC):", x$best_lag_bic, "\n")

  invisible(x)
}
