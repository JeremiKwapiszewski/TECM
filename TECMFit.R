# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

TECMFit <- function(y, X, lag = 0, model = c("linear", "tar", "mtar"),
                    thresh = 0, split = TRUE, break.start = NULL, break.end = NULL) {

  model <- match.arg(model)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  if (!is.ts(y)) y <- ts(y)
  time_index <- time(y)
  start_ts <- start(y)[1]
  freq <- frequency(y)

  # Dummy variable
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
      if (length(break.end) >= i && !is.null(break.end[[i]])) {
        end_i <- to_decimal(break.end[[i]])
      } else {
        end_i <- max(time_index)
      }
      dummy_i <- ifelse(time_index >= start_i & time_index <= end_i, 1, 0)
      dummy_name <- paste0("dummy", i)
      dummy_list[[dummy_name]] <- dummy_i
    }
  }



  #LR
  data_list <- list(y = y)
  for (j in seq_len(ncol(X))) {
    X_ts <- ts(X[, j], start = start_ts, frequency = freq)
    data_list[[colnames(X)[j]]] <- X_ts
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
  LR.model <- lm(formula_LR, data = data_ts_LR)
  ect <- resid(LR.model)

  #1st Differences
  dy <- diff(y)
  dX <- diff(X)

  #Model selection
  if (model == "tar") {
    ect_term1 <- ect * (ect > thresh)
    ect_term2 <- ect * (ect < thresh)
    ect_names <- c("ect_pos", "ect_neg")
  } else if (model == "mtar") {
    ect_term1 <- ect[-1] * (diff(ect) > thresh)
    ect_term2 <- ect[-1] * (diff(ect) < thresh)
    ect_names <- c("ect_pos", "ect_neg")
  } else { # linear
    ect_term1 <- ect
    ect_term2 <- NULL
    ect_names <- "ect"
  }

  #ECM Data preparation
  dy <- ts(dy, start = start_ts, frequency = freq)
  data_list <- list(dy = dy)
  data_list[[ect_names[1]]] <- ts(ect_term1, start = start_ts, frequency = freq)
  if (!is.null(ect_term2)) {
    data_list[[ect_names[2]]] <- ts(ect_term2, start = start_ts, frequency = freq)
  }

  rhs_terms <- c("-1", ect_names)

  for (j in seq_len(ncol(X))) {
    xname <- colnames(X)[j]
    if (split) {
      d_pos <- ts(ifelse(dX[, j] > 0, dX[, j], 0), start = start_ts, frequency = freq)
      d_neg <- ts(ifelse(dX[, j] <= 0, dX[, j], 0), start = start_ts, frequency = freq)
      data_list[[paste0(xname, "_pos")]] <- d_pos
      data_list[[paste0(xname, "_neg")]] <- d_neg

      for (l in 0:lag) {
        rhs_terms <- c(rhs_terms,
                       paste0("L(", xname, "_pos, ", l, ")"),
                       paste0("L(", xname, "_neg, ", l, ")"))
      }
    } else {
      dx_ts <- ts(dX[, j], start = start_ts, frequency = freq)
      data_list[[paste0(xname, ", ")]] <- dx_ts

      for (l in 0:lag) {
        rhs_terms <- c(rhs_terms,
                       paste0("L(", xname, ", ", l, ")"))
      }
    }
  }

  #ECM
  formula_ECM <- as.formula(paste("dy ~", paste(rhs_terms, collapse = " + ")))
  data_ts_ECM <- do.call(cbind, data_list)
  ECM.model <- dynlm(formula_ECM, data = data_ts_ECM)

  return(list(
    formula = formula_ECM,
    ECM_model = ECM.model,
    LR_model = LR.model,
    ECM = summary(ECM.model),
    LR = summary(LR.model),
    LR_residuals = ect,
    ect_pos = ect_term1,
    ect_neg = ect_term2,
    ECM_residuals = resid(ECM.model)
  ))
}
}
