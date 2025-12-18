#' Gradient Descent LASSO with Stability Selection
#'
#' @param formula Formula object. Supports pipes for ZINB (e.g., y ~ x1 + x2 | z1).
#' @param data Data frame.
#' @param family Family object.
#' @param lambda Optional fixed lambda.
#' @param lambda_cv Configuration for CV.
#' @param standardize Logical. Standardize predictors?
#' @param cv_subsample Integer. Speedup for CV.
#' @param parallel Logical. Enable parallel processing?
#' @param n_cores Integer. Number of cores.
#' @param boot Logical. Run stability selection?
#' @param n_boot Number of bootstraps.
#' @param boot_ci Vector of two probabilities for CIs.
#' @param batch_size Integer. Mini-batch SGD.
#' @param warm_start Logical. Warm start bootstraps.
#' @importFrom foreach %dopar%
#' @export
gradLasso <- function(formula, data = NULL, family = grad_gaussian(),
                      lambda = NULL,
                      lambda_cv = TRUE,
                      standardize = TRUE,
                      cv_subsample = NULL,
                      parallel = FALSE,
                      n_cores = NULL,
                      boot = TRUE, n_boot = 50,
                      boot_ci = c(0.025, 0.975),
                      batch_size = NULL,
                      warm_start = TRUE) {
  # Parallel setup ----
  if (parallel) {
    if (!requireNamespace("foreach", quietly = TRUE) || !requireNamespace("doParallel", quietly = TRUE)) {
      stop("Packages 'foreach' and 'doParallel' are required for parallel=TRUE.")
    }
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1

    # Register cluster and ensure it closes on exit
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  # Formula parsing ----
  if (is.null(data)) data <- environment(formula)

  f_text <- paste(deparse(formula), collapse = " ")

  # Handle ZINB pipe syntax (y ~ count | zero)
  if (grepl("\\|", f_text)) {
    if (family$name != "zinb") stop("Pipe syntax '|' is only supported for family='zinb'.")

    parts <- strsplit(f_text, "\\|")[[1]]
    f_count <- as.formula(parts[1], env = environment(formula))

    # Extract response variable to construct the zero-inflation formula
    tt <- terms(f_count, data = data)
    resp_idx <- attr(tt, "response")
    resp_name <- if (resp_idx > 0) deparse(attr(tt, "variables")[[resp_idx + 1]]) else NULL

    mf <- model.frame(f_count, data)
    y <- model.response(mf)

    # Do not remove intercept column; it is needed for correct un-standardization later
    X_mu <- model.matrix(f_count, data)

    f_zero_str <- if (!is.null(resp_name)) paste("~", parts[2], "-", resp_name) else paste("~", parts[2])
    f_zero <- as.formula(f_zero_str, env = environment(formula))
    X_pi <- model.matrix(f_zero, data)
  } else {
    # Standard formula handling
    mf <- model.frame(formula, data)
    y <- model.response(mf)
    X_raw <- model.matrix(formula, data)

    if (family$name == "zinb") {
      X_mu <- X_raw
      X_pi <- X_raw
    } else {
      X <- X_raw
    }
  }

  # Combine ZINB matrices if applicable
  if (family$name == "zinb") {
    if (nrow(X_mu) != nrow(X_pi)) stop("Row mismatch between count and zero matrices.")
    colnames(X_mu) <- paste0("count_", colnames(X_mu))
    colnames(X_pi) <- paste0("zero_", colnames(X_pi))
    X <- cbind(X_mu, X_pi)
  }

  if (nrow(X) == 0) stop("Data has 0 rows.")
  feature_names <- colnames(X)

  # Standardization ----
  X_orig <- X
  scaling_info <- NULL

  if (standardize) {
    X_mean <- colMeans(X)
    X_sd <- apply(X, 2, sd)

    # Identify constant columns to prevent division by zero
    const_cols <- which(X_sd == 0)

    if (length(const_cols) > 0) {
      is_intercept <- apply(X[, const_cols, drop = FALSE], 2, function(v) all(v == 1))
      intercept_indices <- const_cols[is_intercept]

      # Force intercepts to mean=0, sd=1 so scale() leaves them untouched
      if (length(intercept_indices) > 0) {
        X_mean[intercept_indices] <- 0
        X_sd[intercept_indices] <- 1
      }

      other_consts <- const_cols[!is_intercept]
      if (length(other_consts) > 0) {
        X_sd[other_consts] <- 1
        warning(sprintf("Columns %s are constant; scaling skipped.", paste(names(other_consts), collapse = ", ")))
      }
    }

    X <- scale(X, center = X_mean, scale = X_sd)
    scaling_info <- list(mean = X_mean, sd = X_sd)
  }

  # Lambda selection ----
  run_cv <- FALSE
  cv_lambdas_to_test <- NULL

  if (is.logical(lambda_cv) && isTRUE(lambda_cv)) {
    run_cv <- TRUE
  } else if (is.numeric(lambda_cv)) {
    run_cv <- TRUE
    if (length(lambda_cv) == 1 && lambda_cv >= 1 && lambda_cv %% 1 == 0) {
      # Generate sequence if a single integer is provided
      cv_lambdas_to_test <- get_lambda_seq(X, y, family, n_lambda = lambda_cv)
    } else {
      # Use provided vector directly
      cv_lambdas_to_test <- lambda_cv
    }
  } else if (is.null(lambda) && isFALSE(lambda_cv)) {
    run_cv <- TRUE
  }

  final_lambda <- lambda
  cv_res <- NULL

  if (run_cv) {
    cat("Selecting lambda via Cross-Validation...\n")
    cv_res <- cv.gradLasso(X, y, family,
      lambdas = cv_lambdas_to_test,
      batch_size = batch_size, subsample = cv_subsample,
      parallel = parallel
    )
    final_lambda <- cv_res$lambda.min
    cat(sprintf("Optimal Lambda: %f\n", final_lambda))
  }

  if (is.null(final_lambda)) stop("No lambda selected.")

  # Final model fit ----
  cat(sprintf("Fitting Final Model (Lambda=%f)...\n", final_lambda))
  main_beta <- fit_gradlasso(X, y, family, final_lambda, batch_size = batch_size)

  # Assign names to coefficients if missing
  if (is.null(names(main_beta))) {
    n_beta <- length(main_beta)
    n_cols <- ncol(X)

    if (n_beta == n_cols) {
      names(main_beta) <- colnames(X)
    } else if (n_beta == n_cols + 1) {
      # ZINB/NegBin models include an extra dispersion parameter (Theta)
      names(main_beta) <- c(colnames(X), "(Theta)")
    } else {
      warning(sprintf("Coef length (%d) != cols (%d). Names not assigned.", n_beta, n_cols))
    }
  }

  # Stability selection (Bootstrap) ----
  boot_mat <- NULL

  if (boot) {
    cat(sprintf("Running %d Bootstraps...\n", n_boot))
    n_params <- length(main_beta)

    if (parallel) {
      # Foreach handles the parallel backend automatically
      boot_mat <- foreach::foreach(
        b = 1:n_boot, .combine = rbind,
        .packages = "gradLasso"
      ) %dopar% {
        idx <- sample(nrow(X), replace = TRUE)
        fit_gradlasso(X[idx, ], y[idx], family, final_lambda,
          batch_size = batch_size, init_beta = NULL, tol = 1e-5
        )
      }
    } else {
      # Sequential execution allows for warm starts
      boot_mat <- matrix(0, nrow = n_boot, ncol = n_params)
      current_beta <- main_beta

      for (b in 1:n_boot) {
        if (b %% 10 == 0) cat(sprintf("  Iter %d\n", b))

        idx <- sample(nrow(X), replace = TRUE)
        init_b <- if (warm_start) current_beta else NULL
        tol_b <- if (warm_start) 1e-4 else 1e-5

        boot_mat[b, ] <- fit_gradlasso(X[idx, ], y[idx], family, final_lambda,
          batch_size = batch_size, init_beta = init_b, tol = tol_b
        )

        if (warm_start) current_beta <- boot_mat[b, ]
      }
    }

    if (!is.null(names(main_beta)) && !is.null(boot_mat)) {
      colnames(boot_mat) <- names(main_beta)
    }
  }

  # Un-standardization ----
  # Helper function to restore coefficients to original scale
  unscale_coefs <- function(beta_vec, info) {
    if (is.null(info) || is.null(names(beta_vec))) {
      return(beta_vec)
    }
    new_beta <- beta_vec

    # 1. Slope Adjustment: beta_orig = beta_scaled / sd
    common_names <- intersect(names(info$sd), names(beta_vec))
    if (length(common_names) > 0) {
      sd_ordered <- info$sd[common_names]
      new_beta[common_names] <- new_beta[common_names] / sd_ordered
    }

    # 2. Intercept Adjustment: beta_0 = beta_0_scaled - sum(beta_i * mean_i)
    int_indices <- grep("Intercept", names(beta_vec), ignore.case = TRUE)
    if (length(int_indices) == 0) {
      return(new_beta)
    }

    # Handle split intercepts for ZINB
    count_int <- grep("count_.*Intercept", names(beta_vec))
    zero_int <- grep("zero_.*Intercept", names(beta_vec))

    if (length(count_int) > 0 || length(zero_int) > 0) {
      if (length(count_int) > 0) {
        c_preds <- grep("^count_", names(info$mean), value = TRUE)
        c_preds <- c_preds[!grepl("Intercept", c_preds)]

        if (length(c_preds) > 0) {
          adj <- sum(new_beta[c_preds] * info$mean[c_preds])
          new_beta[count_int] <- new_beta[count_int] - adj
        }
      }

      if (length(zero_int) > 0) {
        z_preds <- grep("^zero_", names(info$mean), value = TRUE)
        z_preds <- z_preds[!grepl("Intercept", z_preds)]

        if (length(z_preds) > 0) {
          adj <- sum(new_beta[z_preds] * info$mean[z_preds])
          new_beta[zero_int] <- new_beta[zero_int] - adj
        }
      }
    } else {
      # Standard single intercept
      pred_names <- names(info$mean)
      pred_names <- pred_names[!grepl("Intercept", pred_names, ignore.case = TRUE)]

      if (length(pred_names) > 0) {
        adj <- sum(new_beta[pred_names] * info$mean[pred_names])
        new_beta[int_indices] <- new_beta[int_indices] - adj
      }
    }

    return(new_beta)
  }

  if (standardize && !is.null(scaling_info)) {
    main_beta <- unscale_coefs(main_beta, scaling_info)

    # Apply unscaling to bootstrap matrix slopes
    if (!is.null(boot_mat)) {
      common_cols <- intersect(names(scaling_info$sd), colnames(boot_mat))
      if (length(common_cols) > 0) {
        boot_mat[, common_cols] <- sweep(
          boot_mat[, common_cols, drop = FALSE],
          2, scaling_info$sd[common_cols], `/`
        )
      }
    }
  }

  # Return object ----
  # Strip theta parameter from beta vector if necessary for prediction
  beta_for_pred <- main_beta
  if (length(main_beta) == ncol(X_orig) + 1) {
    beta_for_pred <- main_beta[1:ncol(X_orig)]
  }

  fitted <- family$predict(X_orig, beta_for_pred, type = "response")
  final_deviance <- family$deviance(X_orig, y, main_beta)
  residuals <- y - fitted
  n_obs <- nrow(X_orig)

  structure(list(
    coefficients = main_beta,
    fitted.values = fitted,
    residuals = residuals,
    boot_matrix = boot_mat,
    lambda = final_lambda,
    family = family,
    feature_names = feature_names,
    formula = formula,
    cv_results = cv_res,
    boot_ci = boot_ci,
    deviance = final_deviance,
    nobs = n_obs,
    standardized = standardize
  ), class = "gradLasso")
}
