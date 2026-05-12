TECMFit <- function (y, X, lag = 0, model = c("linear", "tar", "mtar"),
          thresh = 0, split = TRUE, break.start = NULL, break.end = NULL)
{
  model <- match.arg(model)
  if (!is.matrix(X))
    X <- as.matrix(X)
  if (is.null(colnames(X)))
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  if (!is.ts(y))
    y <- ts(y)
  
  time_index <- time(y)
  start_ts   <- start(y)
  freq       <- frequency(y)
  
  to_decimal <- function(date) date[1] + (date[2] - 1) / freq
  
  dummy_list <- list()
  if (!is.null(break.start)) {
    if (!is.list(break.start)) break.start <- list(break.start)
    if (is.null(break.end)) {
      break.end <- vector("list", length(break.start))
    } else if (!is.list(break.end)) {
      break.end <- list(break.end)
    }
    for (i in seq_along(break.start)) {
      start_i <- to_decimal(break.start[[i]])
      end_i   <- if (!is.null(break.end[[i]])) to_decimal(break.end[[i]]) else max(time_index)
      dummy_list[[paste0("dummy", i)]] <- ifelse(
        time_index >= start_i & time_index <= end_i, 1, 0)
    }
  }
  
  #Long-run
  data_list <- list(y = y)
  for (j in seq_len(ncol(X))) {
    data_list[[colnames(X)[j]]] <- ts(X[, j], start = start_ts, frequency = freq)
  }
  rhs_terms <- colnames(X)
  if (length(dummy_list) > 0) {
    for (i in seq_along(dummy_list)) {
      name <- names(dummy_list)[i]
      data_list[[name]] <- dummy_list[[i]]
      rhs_terms <- c(rhs_terms, name)
    }
  }
  formula_LR <- as.formula(paste("y ~", paste(rhs_terms, collapse = " + ")))
  data_ts_LR <- do.call(cbind, data_list)
  LR.model   <- lm(formula_LR, data = data_ts_LR)
  ect        <- resid(LR.model)
  
  #ECT
  if (model == "tar") {
    ect_term1 <- ts(ect * (ect >  thresh), start = start_ts, frequency = freq)
    ect_term2 <- ts(ect * (ect <= thresh), start = start_ts, frequency = freq)
    ect_names <- c("ect_pos", "ect_neg")
  } else if (model == "mtar") {
    d_ect     <- c(NA, diff(ect))
    ect_term1 <- ts(ect * (d_ect >  thresh), start = start_ts, frequency = freq)
    ect_term2 <- ts(ect * (d_ect <= thresh), start = start_ts, frequency = freq)
    ect_names <- c("ect_pos", "ect_neg")
  } else {
    ect_term1 <- ts(ect, start = start_ts, frequency = freq)
    ect_term2 <- NULL
    ect_names <- "ect"
  }
  
  #Short-run dynamics
  dy <- ts(diff(y), start = start_ts, frequency = freq)
  dX <- diff(X)
  
  data_list <- list(dy = dy)
  data_list[[ect_names[1]]] <- ect_term1
  if (!is.null(ect_term2))
    data_list[[ect_names[2]]] <- ect_term2
  
  rhs_terms <- c(-1, ect_names)
  
  if (length(dummy_list) > 0) {
    for (i in seq_along(dummy_list)) {
      dname <- paste0("d", names(dummy_list)[i])
      data_list[[dname]] <- ts(diff(dummy_list[[i]]), start = start_ts, frequency = freq)
      rhs_terms <- c(rhs_terms, dname)
    }
  }
  
  for (j in seq_len(ncol(X))) {
    xname <- colnames(X)[j]
    if (split) {
      data_list[[paste0(xname, "_pos")]] <- ts(ifelse(dX[, j] >  0, dX[, j], 0),
                                               start = start_ts, frequency = freq)
      data_list[[paste0(xname, "_neg")]] <- ts(ifelse(dX[, j] <= 0, dX[, j], 0),
                                               start = start_ts, frequency = freq)
      for (l in 0:lag) {
        rhs_terms <- c(rhs_terms,
                       paste0("L(", xname, "_pos, ", l, ")"),
                       paste0("L(", xname, "_neg, ", l, ")"))
      }
    } else {
      data_list[[xname]] <- ts(dX[, j], start = start_ts, frequency = freq)
      for (l in 0:lag) {
        rhs_terms <- c(rhs_terms, paste0("L(", xname, ", ", l, ")"))
      }
    }
  }
  
  formula_ECM <- as.formula(paste("dy ~", paste(rhs_terms, collapse = " + ")))
  data_ts_ECM <- do.call(cbind, data_list)
  ECM.model   <- dynlm(formula_ECM, data = data_ts_ECM)
  
  #Output
  out <- list(
    formula       = formula_ECM,
    ECM_model     = ECM.model,
    LR_model      = LR.model,
    ECM           = summary(ECM.model),
    LR            = summary(LR.model),
    LR_residuals  = ect,
    ECM_residuals = resid(ECM.model),
    LR_fitted     = fitted(LR.model),   
    ECM_fitted    = fitted(ECM.model)   
  )
  if (model != "linear") {
    out$ect_pos <- ect_term1
    out$ect_neg <- ect_term2
  } else {
    out$ect <- ect_term1
  }
  return(out)
}
