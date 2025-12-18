#' Internal FISTA Solver
#' 
#' @keywords internal
fit_gradlasso <- function(X, y, family, lambda, alpha = 1, 
                          max_iter = 2000, tol = 1e-5, step_size = 0.001,
                          batch_size = NULL, init_beta = NULL) {
  
  n_samples <- nrow(X)
  p <- ncol(X)
  
  # Initialize coefficients ----
  if (!is.null(init_beta)) {
    beta <- init_beta
  } else {
    if (family$name == "zinb" || family$name == "negbin") {
      # ZINB and NegBin include coefficients matching columns plus one theta parameter
      beta <- c(rep(0, p), 0)
    } else {
      # Gaussian and Binomial coefficients match columns exactly
      beta <- rep(0, p)
    }
  }
  
  beta_prev <- beta
  t_k <- 1
  
  # Optimization loop ----
  for (i in 1:max_iter) {
    
    # Handle mini-batch sampling if batch_size is provided
    if (is.null(batch_size)) {
      X_b <- X
      y_b <- y
    } else {
      actual_size <- min(batch_size, n_samples)
      idx <- sample(n_samples, actual_size)
      X_b <- X[idx, , drop=FALSE]
      y_b <- y[idx]
    }
    
    # Compute gradient and apply descent step
    grad <- family$gradient(X_b, y_b, beta)
    
    beta_smooth <- beta - step_size * grad
    thresh <- step_size * lambda * alpha
    
    beta_next <- beta_smooth
    
    # Apply soft-thresholding penalization
    # Note: Theta (dispersion parameter) is never penalized
    if (family$name %in% c("zinb", "negbin")) {
       # Exclude the last element (Theta) from penalization
       pen_idx <- 1:(length(beta) - 1)
       beta_next[pen_idx] <- sign(beta_smooth[pen_idx]) * pmax(0, abs(beta_smooth[pen_idx]) - thresh)
       beta_next[length(beta)] <- beta_smooth[length(beta)]
    } else {
       # Penalize all coefficients for standard families
       beta_next <- sign(beta_smooth) * pmax(0, abs(beta_smooth) - thresh)
    }
    
    # Apply Nesterov momentum (FISTA) if not using mini-batches
    if (is.null(batch_size)) {
      t_next <- (1 + sqrt(1 + 4*t_k^2)) / 2
      beta <- beta_next + ((t_k - 1) / t_next) * (beta_next - beta_prev)
      
      # Check convergence
      if (max(abs(beta - beta_prev)) < tol) break
      
      beta_prev <- beta_next
      t_k <- t_next
    } else {
      # No momentum for stochastic updates
      beta <- beta_next
    }
  }
  return(beta)
}
