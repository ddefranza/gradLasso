#' Generate Lambda Sequence
#'
#' @keywords internal
get_lambda_seq <- function(X, y, family, n_lambda = 20, ratio = 0.001) {
  n <- nrow(X)
  p <- ncol(X)

  # Initialize dummy coefficients (all zeros)
  if (family$name == "zinb") {
    dummy <- c(rep(0, p), 0)
  } else if (family$name == "negbin") {
    dummy <- c(rep(0, p), 0)
  } else {
    dummy <- rep(0, p)
  }

  # Calculate gradient at zero to determine max lambda
  g <- family$gradient(X, y, dummy)
  lambda_max <- max(abs(g))
  lambda_min <- lambda_max * ratio

  return(exp(seq(log(lambda_max), log(lambda_min), length.out = n_lambda)))
}

#' Cross-Validation for gradLasso
#'
#' @param object Matrix X (predictors).
#' @param data Vector y (response).
#' @param family Family object (e.g., grad_gaussian, grad_zinb).
#' @param lambdas Vector of lambda values to test. If NULL, a sequence is generated.
#' @param nfolds Integer. Number of CV folds (default 5).
#' @param batch_size Integer. Mini-batch size for SGD.
#' @param subsample Integer. Number of rows to use for CV (if NULL, uses all data).
#' @param parallel Logical. If TRUE, runs folds in parallel.
#' @param verbose Logical. Print progress to console?
#' @return A list containing CV results (mean error, SD, optimal lambdas).
#' @export
cv.gradLasso <- function(object, data = NULL, family, lambdas = NULL, nfolds = 5,
                         batch_size = NULL, subsample = NULL, parallel = FALSE, verbose = FALSE) {
  if (inherits(object, "formula")) stop("cv.gradLasso expects X matrix and y vector.")

  X <- object
  y <- data
  n_total <- nrow(X)

  # Subsampling for speed ----
  if (!is.null(subsample) && subsample < n_total) {
    if (verbose) cat(sprintf("Subsampling CV: Using %d rows out of %d...\n", subsample, n_total))
    idx_sub <- sample(n_total, subsample)
    X <- X[idx_sub, , drop = FALSE]
    y <- y[idx_sub]
    n <- subsample
  } else {
    n <- n_total
  }

  if (n < nfolds) nfolds <- n
  if (is.null(lambdas)) lambdas <- get_lambda_seq(X, y, family, n_lambda = 20)

  folds <- sample(rep(1:nfolds, length.out = n))
  cv_err <- matrix(NA, length(lambdas), nfolds)

  if (verbose) {
    cat(sprintf("Running %d-fold CV over %d lambdas", nfolds, length(lambdas)))
    if (parallel) cat(" (Parallel)...\n") else cat("...\n")
  }

  # Cross-validation loop ----
  for (i in seq_along(lambdas)) {
    lam <- lambdas[i]

    if (parallel) {
      # Parallel execution: distributing folds across workers
      # Note: Relies on cluster registered in parent gradLasso() call
      err_vec <- foreach::foreach(k = 1:nfolds, .combine = c, .packages = "gradLasso") %dopar% {
        test_idx <- which(folds == k)
        X_train <- X[-test_idx, , drop = FALSE]
        y_train <- y[-test_idx]
        X_test <- X[test_idx, , drop = FALSE]
        y_test <- y[test_idx]

        coefs <- fit_gradlasso(X_train, y_train, family, lam, batch_size = batch_size)
        family$deviance(X_test, y_test, coefs)
      }
      cv_err[i, ] <- err_vec
    } else {
      # Sequential execution
      for (k in 1:nfolds) {
        test_idx <- which(folds == k)
        X_train <- X[-test_idx, , drop = FALSE]
        y_train <- y[-test_idx]
        X_test <- X[test_idx, , drop = FALSE]
        y_test <- y[test_idx]

        coefs <- fit_gradlasso(X_train, y_train, family, lam, batch_size = batch_size)
        cv_err[i, k] <- family$deviance(X_test, y_test, coefs)
      }
    }
  }

  # Calculate stats and optimal lambda ----
  mean_cv <- rowMeans(cv_err, na.rm = TRUE)
  sd_cv <- apply(cv_err, 1, sd, na.rm = TRUE) / sqrt(nfolds)

  # Lambda min: minimizes CV error
  best_idx <- which.min(mean_cv)

  # Lambda 1se: largest lambda within 1 SE of the minimum
  cutoff <- mean_cv[best_idx] + sd_cv[best_idx]
  valid_idx <- which(mean_cv <= cutoff)
  se_idx <- valid_idx[which.max(lambdas[valid_idx])]

  structure(list(
    lambda = lambdas,
    cvm = mean_cv,
    cvsd = sd_cv,
    lambda.min = lambdas[best_idx],
    lambda.1se = lambdas[se_idx]
  ), class = "cv.gradLasso")
}
