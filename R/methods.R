#' Predict method for gradLasso
#'
#' @param object A gradLasso fitted object.
#' @param newdata Optional new data frame for prediction. If missing, returns fitted values.
#' @param type Type of prediction: "response" (default), "link", "count" (mu), or "zero" (pi).
#' @param ... Additional arguments passed to methods.
#' @return A vector or matrix of predictions.
#' @export
predict.gradLasso <- function(object, newdata, type = c("response", "link", "count", "zero"), ...) {
  type <- match.arg(type)

  if (missing(newdata)) {
    if (type == "response") {
      return(object$fitted.values)
    }
    stop("For prediction types other than 'response', please provide 'newdata'.")
  }

  # Construct model matrix from new data
  mf <- model.frame(object$terms, newdata, na.action = na.pass)
  X <- model.matrix(object$terms, mf)

  # Ensure column consistency with training data
  if (!is.null(object$feature_names)) {
    valid_cols <- object$feature_names
    if (all(valid_cols %in% colnames(X))) {
      X <- X[, valid_cols, drop = FALSE]
    } else {
      stop("Mismatch in columns between model and new data.")
    }
  }

  pred <- object$family$predict(X, object$coefficients, type = type)
  return(pred)
}

#' Summary method for gradLasso
#'
#' @param object A gradLasso fitted object.
#' @param ... Additional arguments passed to methods.
#' @return A list containing the coefficient table, fit statistics (AIC/BIC), and stability selection results.
#' @export
summary.gradLasso <- function(object, ...) {
  beta <- object$coefficients
  names <- object$feature_names
  n_total <- length(beta)

  # Calculate fit statistics
  k <- sum(abs(beta) > 1e-5)
  n <- object$nobs
  dev <- object$deviance

  if (!is.null(n) && !is.null(dev)) {
    aic <- dev + 2 * k
    bic <- dev + k * log(n)
  } else {
    aic <- NULL
    bic <- NULL
  }

  # Organize components for printing
  components <- rep("Main", n_total)
  predictors <- if (length(names) == n_total) names else c(names, rep("Extra", n_total - length(names)))

  # Handle special families (ZINB, NegBin)
  if (object$family$name %in% c("zinb", "negbin")) {
    components[n_total] <- "Dispersion"
    predictors[n_total] <- "(Theta)"

    if (object$family$name == "zinb") {
      if (any(grepl("^zero_", names))) {
        components[grep("^count_", names)] <- "Count"
        components[grep("^zero_", names)] <- "Zero-Infl"
        predictors <- gsub("^count_|^zero_", "", names)
        predictors[n_total] <- "(Theta)"
      } else if (length(names) < (n_total - 1)) {
        # Fallback for when names are not prefixed but split exists
        k_len <- length(names)
        predictors <- c(names, names, "(Theta)")
        components[1:k_len] <- "Count"
        components[(k_len + 1):(2 * k_len)] <- "Zero-Infl"
      }
    }
  }

  # Build statistics table
  if (!is.null(object$boot_matrix)) {
    sel_prob <- colMeans(abs(object$boot_matrix) > 1e-5)
    boot_mean <- colMeans(object$boot_matrix)
    ci_probs <- if (!is.null(object$boot_ci)) object$boot_ci else c(0.025, 0.975)

    ci_lower <- apply(object$boot_matrix, 2, quantile, probs = ci_probs[1])
    ci_upper <- apply(object$boot_matrix, 2, quantile, probs = ci_probs[2])

    res_df <- data.frame(
      Component = components,
      Predictor = predictors,
      Estimate = beta,
      Selection_Prob = sel_prob,
      Boot_Mean = boot_mean,
      CI_Low = ci_lower,
      CI_High = ci_upper
    )
    ci_names <- paste0("CI_", format(ci_probs * 100, trim = TRUE, digits = 3))
    colnames(res_df)[6:7] <- ci_names
    res_df <- res_df[order(-res_df$Selection_Prob), ]
  } else {
    res_df <- data.frame(Component = components, Predictor = predictors, Estimate = beta)
  }

  ans <- list(
    coefficients = res_df, family = object$family$name,
    lambda = object$lambda, deviance = dev, aic = aic, bic = bic, df = k,
    is_cv = !is.null(object$cv_results),
    boot = !is.null(object$boot_matrix),
    n_boot = if (!is.null(object$boot_matrix)) nrow(object$boot_matrix) else 0,
    boot_ci = if (!is.null(object$boot_matrix)) object$boot_ci else NULL
  )
  class(ans) <- "summary.gradLasso"
  return(ans)
}

#' Print method for summary
#'
#' @param x A summary.gradLasso object.
#' @param ... Additional arguments passed to print.
#' @return Invisibly returns the input object.
#' @export
print.summary.gradLasso <- function(x, ...) {
  cat("\n------------------------------------------------\n")
  cat("gradLasso Model Summary\n")
  cat("------------------------------------------------\n")
  cat(sprintf("Family:   %s\n", x$family))

  if (!is.null(x$deviance)) {
    cat(sprintf("Deviance: %.2f\n", x$deviance))
    cat(sprintf("AIC:      %.2f\n", x$aic))
    cat(sprintf("BIC:      %.2f\n", x$bic))
    cat(sprintf("DF:       %d\n", x$df))
  } else {
    cat("[Fit statistics missing. Re-run model to view AIC/BIC]\n")
  }

  cat("------------------------------------------------\n")
  lam_str <- format(x$lambda, digits = 4)
  cv_str <- if (x$is_cv) "(Selected via CV)" else "(User-Defined)"
  cat(sprintf("Lambda:   %s %s\n", lam_str, cv_str))

  if (x$boot) {
    cat(sprintf("Method:   Stability Selection (%d bootstraps)\n", x$n_boot))
    if (!is.null(x$boot_ci)) {
      cat(sprintf("Interval: %s%% - %s%%\n", x$boot_ci[1] * 100, x$boot_ci[2] * 100))
    }
  } else {
    cat("Method:   Single LASSO Fit\n")
  }

  # Filter for printing: keep selected vars or those with high selection probability
  df <- x$coefficients
  if ("Selection_Prob" %in% colnames(df)) {
    df <- df[df$Selection_Prob > 0 | abs(df$Estimate) > 1e-5, ]
  } else {
    df <- df[abs(df$Estimate) > 1e-5, ]
  }

  print_block <- function(sub_df, title) {
    cat(sprintf("\n--- %s ---\n", title))
    if (nrow(sub_df) > 0) {
      sub_df$Component <- NULL
      row.names(sub_df) <- NULL
      nums <- sapply(sub_df, is.numeric)
      sub_df[nums] <- round(sub_df[nums], 4)
      print(sub_df, row.names = FALSE)
    } else {
      cat("No coefficients selected.\n")
    }
  }

  if (x$family == "zinb") {
    print_block(df[df$Component == "Count", ], "Count Model Coefficients")
    print_block(df[df$Component == "Zero-Infl", ], "Zero-Inflation Model Coefficients")
    if (any(df$Component == "Dispersion")) cat(sprintf("\n--- Dispersion ---\nTheta: %.4f\n", df[df$Component == "Dispersion", "Estimate"]))
  } else if (x$family == "negbin") {
    print_block(df[df$Component %in% c("Main", "Count"), ], "Negative Binomial Coefficients")
    if (any(df$Component == "Dispersion")) cat(sprintf("\n--- Dispersion ---\nTheta: %.4f\n", df[df$Component == "Dispersion", "Estimate"]))
  } else {
    print_block(df, "Selected Coefficients")
  }
  cat("------------------------------------------------\n")
  invisible(x)
}

#' Print method for gradLasso object
#'
#' @param x A gradLasso fitted object.
#' @param ... Additional arguments passed to print.
#' @return Invisibly returns the input object.
#' @export
print.gradLasso <- function(x, ...) {
  cat("\ngradLasso Fitted Object\n")
  cat("Family:", x$family$name, "\n")
  cat("Lambda:", format(x$lambda, digits = 4), "\n")
  if (!is.null(x$deviance)) cat("Deviance:", format(x$deviance, digits = 2), "\n")
  cat("Use plot() to view diagnostics or summary() for coefficients.\n")
}

#' Master Plot Method
#'
#' Diagnostic plots for gradLasso objects (Stability, CV, Residuals).
#'
#' @param x A gradLasso fitted object.
#' @param which Integer vector specifying which plots to draw (1:5).
#' @param ... Additional arguments passed to plotting functions.
#' @return Invisibly returns NULL.
#' @export
plot.gradLasso <- function(x, which = c(1L:5L), ...) {
  plots <- list()
  plot_names <- character(5)

  # 1. Stability Selection Plot
  if (!is.null(x$boot_matrix)) {
    plot_names[1] <- "Stability Selection"
    plots[[1]] <- function() {
      # Save and restore graphical parameters
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      sel <- colMeans(abs(x$boot_matrix) > 1e-5)
      sel <- sort(sel, decreasing = TRUE)
      if (length(sel) > 30) sel <- sel[1:30]
      par(mar = c(7, 4, 4, 2))
      barplot(sel,
        main = "Stability Selection (Top 30)", ylab = "Selection Probability",
        col = ifelse(sel >= 0.8, "steelblue", "lightgrey"), border = NA, las = 2, cex.names = 0.8
      )
      abline(h = 0.8, col = "red", lty = 2)
      legend("topright", legend = c("Stable (>0.8)", "Unstable"), fill = c("steelblue", "lightgrey"), bty = "n")
    }
  }

  # 2. CV Deviance Plot
  if (!is.null(x$cv_results)) {
    plot_names[2] <- "CV Deviance"
    plots[[2]] <- function() {
      # Save and restore graphical parameters
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      par(mar = c(5, 4, 4, 2) + 0.1)
      cv <- x$cv_results
      plot(log(cv$lambda), cv$cvm, type = "b", pch = 19, col = "darkgrey", xlab = "Log(Lambda)", ylab = "CV Deviance", main = "Cross-Validation Error")
      abline(v = log(cv$lambda.min), col = "red", lty = 2, lwd = 2)
      if (!is.null(cv$lambda.1se)) abline(v = log(cv$lambda.1se), col = "blue", lty = 2, lwd = 2)
      legend("topright", legend = c("Min", "1-SE"), col = c("red", "blue"), lty = 2)
    }
  }

  # 3-5. Residual Plots
  if (!is.null(x$residuals)) {
    y_hat <- x$fitted.values
    residuals <- x$residuals

    plot_names[3] <- "Residuals vs Fitted"
    plots[[3]] <- function() {
      # Save and restore graphical parameters
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      par(mar = c(5, 4, 4, 2) + 0.1)
      plot(y_hat, residuals, main = "Residuals vs Fitted", xlab = "Fitted Values", ylab = "Residuals", pch = 19, col = rgb(0, 0, 0, 0.3))
      abline(h = 0, col = "red", lty = 2)
      try(
        {
          lines(lowess(y_hat, residuals), col = "blue")
        },
        silent = TRUE
      )
    }

    plot_names[4] <- "Normal Q-Q"
    plots[[4]] <- function() {
      # Save and restore graphical parameters
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      par(mar = c(5, 4, 4, 2) + 0.1)
      qqnorm(residuals, main = "Normal Q-Q Plot", pch = 19, col = rgb(0, 0, 0, 0.5))
      qqline(residuals, col = "red", lwd = 2)
    }

    plot_names[5] <- "Residual Histogram"
    plots[[5]] <- function() {
      # Save and restore graphical parameters
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      par(mar = c(5, 4, 4, 2) + 0.1)
      hist(residuals, breaks = 30, col = "steelblue", border = "white", main = "Histogram of Residuals", xlab = "Residuals")
    }
  }

  available_indices <- intersect(which(!sapply(plots, is.null)), which)
  if (length(available_indices) == 0) {
    message("No valid plots available to display.")
    return(invisible())
  }

  # Direct output for file devices, non-interactive sessions, OR CI environments
  is_file_device <- any(grepl("pdf|png|jpeg|tiff|bmp|postscript", names(dev.cur()), ignore.case = TRUE))

  # FIX: Added `|| !interactive()` to ensure CI environments never hit menu()
  if (is_file_device || !interactive()) {
    for (i in available_indices) plots[[i]]()
    return(invisible())
  }

  cat("\nDiagnostic Plots Available:\n")
  choices <- plot_names[available_indices]
  repeat {
    if (exists("flush.console")) flush.console()
    selection <- menu(choices, title = "Select a plot (0 to exit):")
    if (selection == 0) break
    plots[[available_indices[selection]]]()
  }
}

#' Print CV results
#'
#' @param x A cv.gradLasso fitted object.
#' @param ... Additional arguments passed to print.
#' @return Invisibly returns the input object.
#' @export
print.cv.gradLasso <- function(x, ...) {
  cat("Cross-Validation Results\n")
  cat("Min Lambda:", x$lambda.min, "\n")
  cat("Min Deviance:", min(x$cvm), "\n")
}

#' Plot CV results (Standalone)
#'
#' @param x A cv.gradLasso fitted object.
#' @param ... Additional arguments passed to plot.
#' @return Invisibly returns NULL.
#' @export
plot.cv.gradLasso <- function(x, ...) {
  plot(log(x$lambda), x$cvm, type = "b", pch = 19, main = "CV Deviance")
  abline(v = log(x$lambda.min), col = "red", lty = 2)
}

# New S3 methods ----

#' Extract Model Coefficients
#'
#' @param object A gradLasso fitted object.
#' @param ... Additional arguments.
#' @return A numeric vector of coefficients.
#' @export
coef.gradLasso <- function(object, ...) {
  return(object$coefficients)
}

#' Extract Fitted Values
#'
#' @param object A gradLasso fitted object.
#' @param ... Additional arguments.
#' @return A numeric vector of fitted values.
#' @export
fitted.gradLasso <- function(object, ...) {
  return(object$fitted.values)
}

#' Extract Residuals
#'
#' @param object A gradLasso fitted object.
#' @param ... Additional arguments.
#' @return A numeric vector of residuals.
#' @export
residuals.gradLasso <- function(object, ...) {
  if (!is.null(object$residuals)) {
    return(object$residuals)
  }
  warning("Residuals not stored in object.")
  return(NULL)
}

#' Extract Log-Likelihood
#'
#' @param object A gradLasso fitted object.
#' @param ... Additional arguments.
#' @return An object of class \code{logLik}.
#' @export
logLik.gradLasso <- function(object, ...) {
  if (is.null(object$deviance)) {
    return(NA)
  }

  # Deviance = -2 * logLik
  val <- -object$deviance / 2

  # Calculate degrees of freedom (non-zero coefficients)
  k <- sum(abs(object$coefficients) > 1e-5)

  attr(val, "nobs") <- object$nobs
  attr(val, "df") <- k
  class(val) <- "logLik"
  return(val)
}
