TECMThresh <- function(y, X, model = c("tar", "mtar"),
                       thresh.grid = NULL,
                       grid.n = 100,
                       trim = 0.15,
                       lag = 0,
                       split = TRUE,
                       break.start = NULL,
                       break.end = NULL) {

  model <- match.arg(model)
  if (!is.ts(y)) y <- ts(y)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))

  #LR
  LR.model <- lm(y ~ X)
  ect <- resid(LR.model)

  #Grid
  ect_sorted <- sort(ect)
  n <- length(ect_sorted)
  lower <- floor(n * trim)
  upper <- ceiling(n * (1 - trim))

  if (is.null(thresh.grid)) {
    thresh.grid <- ect_sorted[seq(lower, upper, length.out = grid.n)]
  }

  #RSS
  rss_vec <- numeric(length(thresh.grid))

  for (i in seq_along(thresh.grid)) {
    threshold <- thresh.grid[i]

    model_fit <- try(TECMFit(y = y, X = X,
                             lag = lag,
                             model = model,
                             thresh = threshold,
                             split = split,
                             break.start = break.start,
                             break.end = break.end), silent = TRUE)

    if (!inherits(model_fit, "try-error")) {
      rss_vec[i] <- sum(resid(model_fit$ECM_model)^2, na.rm = TRUE)
    } else {
      rss_vec[i] <- NA
    }
  }

  best_idx <- which.min(rss_vec)
  best_thresh <- thresh.grid[best_idx]

  return(list(
    threshold = best_thresh,
    rss = rss_vec,
    grid = thresh.grid
  ))
}
