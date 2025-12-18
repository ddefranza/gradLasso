#' Gaussian Family (Least Squares)
#'
#' @return A list containing gradient, deviance, and prediction functions for Gaussian regression.
#' @export
grad_gaussian <- function() {
  list(
    name = "gaussian",
    gradient = function(X, y, beta) {
      resid <- (X %*% beta) - y
      crossprod(X, resid) / length(y)
    },
    deviance = function(X, y, beta) {
      resid <- (X %*% beta) - y
      sum(resid^2)
    },
    predict = function(X, beta, type = "response") {
      eta <- drop(X %*% beta)
      if (type == "link") {
        return(eta)
      }
      eta
    }
  )
}

#' Binomial Family (Logistic Regression)
#'
#' @return A list containing gradient, deviance, and prediction functions for logistic regression.
#' @export
grad_binomial <- function() {
  list(
    name = "binomial",
    gradient = function(X, y, beta) {
      probs <- 1 / (1 + exp(-X %*% beta))
      crossprod(X, probs - y) / length(y)
    },
    deviance = function(X, y, beta) {
      probs <- 1 / (1 + exp(-X %*% beta))

      # Clamp probabilities to avoid log(0)
      probs <- pmax(1e-15, pmin(1 - 1e-15, probs))

      ll <- sum(y * log(probs) + (1 - y) * log(1 - probs))
      -2 * ll
    },
    predict = function(X, beta, type = "response") {
      eta <- drop(X %*% beta)
      if (type == "link") {
        return(eta)
      }
      1 / (1 + exp(-eta))
    }
  )
}

#' Negative Binomial Family
#'
#' @return A list containing gradient, deviance, and prediction functions for Negative Binomial regression.
#' @export
grad_negbin <- function() {
  list(
    name = "negbin",
    deviance = function(X, y, beta) {
      # Extract parameters: last element is log(theta)
      theta <- exp(tail(beta, 1))
      beta_mu <- head(beta, -1)

      mu <- exp(X %*% beta_mu)

      ll <- lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
        theta * log(theta) + y * log(mu) -
        (y + theta) * log(theta + mu)

      -2 * sum(ll)
    },
    gradient = function(X, y, beta) {
      n <- length(y)

      # Extract parameters
      theta <- exp(tail(beta, 1))
      beta_mu <- head(beta, -1)

      mu <- exp(drop(X %*% beta_mu))

      # Gradient wrt linear predictor (eta)
      grad_eta <- y - mu * (y + theta) / (mu + theta)

      # Gradient wrt theta
      dL_dtheta <- digamma(y + theta) - digamma(theta) +
        log(theta) - log(theta + mu) +
        1 - (y + theta) / (mu + theta)

      # Chain rule for log(theta) parameterization: dL/d(log_theta) = dL/dtheta * theta
      grad_zeta <- dL_dtheta * theta

      g_mu <- -crossprod(X, grad_eta) / n
      g_theta <- -sum(grad_zeta) / n

      c(g_mu, g_theta)
    },
    predict = function(X, beta, type = "response") {
      beta_mu <- head(beta, -1)
      eta <- drop(X %*% beta_mu)

      if (type == "link") {
        return(eta)
      }
      exp(eta)
    }
  )
}

#' Zero-Inflated Negative Binomial Family
#'
#' @return A list containing gradient, deviance, and prediction functions for ZINB regression.
#' @export
grad_zinb <- function() {
  list(
    name = "zinb",

    # Helper to separate count and zero components from the unified design matrix
    split_components = function(X, beta) {
      cols <- colnames(X)

      # Identify column indices by prefix
      idx_mu <- grep("^count_", cols)
      idx_pi <- grep("^zero_", cols)

      # Extract parameters (last element is always theta)
      theta <- exp(tail(beta, 1))
      beta_mu <- beta[idx_mu]
      beta_pi <- beta[idx_pi]

      # Compute linear predictors
      eta_mu <- drop(X[, idx_mu, drop = FALSE] %*% beta_mu)
      eta_pi <- drop(X[, idx_pi, drop = FALSE] %*% beta_pi)

      list(
        eta_mu = eta_mu, eta_pi = eta_pi, theta = theta,
        idx_mu = idx_mu, idx_pi = idx_pi
      )
    },
    deviance = function(X, y, beta) {
      comps <- grad_zinb()$split_components(X, beta)

      mu <- exp(comps$eta_mu)
      pi <- 1 / (1 + exp(-comps$eta_pi))
      theta <- comps$theta

      # Probability of zero from the Negative Binomial part
      prob_nb_0 <- (theta / (theta + mu))^theta

      # Total probability of zero (Structural + NB)
      prob_zero <- pi + (1 - pi) * prob_nb_0

      # Log-likelihood calculation
      ll_zero <- log(pmax(1e-15, prob_zero))
      ll_nonzero <- log(pmax(1e-15, 1 - pi)) +
        lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
        theta * log(theta / (theta + mu)) +
        y * log(mu / (theta + mu))

      log_lik <- ifelse(y == 0, ll_zero, ll_nonzero)
      -2 * sum(log_lik)
    },
    gradient = function(X, y, beta) {
      n <- length(y)

      comps <- grad_zinb()$split_components(X, beta)
      theta <- comps$theta

      mu <- exp(comps$eta_mu)
      pi <- 1 / (1 + exp(-comps$eta_pi))

      prob_nb_0 <- (theta / (theta + mu))^theta
      prob_zero <- pi + (1 - pi) * prob_nb_0

      # Gradient wrt Count Model (mu) ----
      grad_eta_mu_nz <- y - mu * (y + theta) / (mu + theta)
      grad_eta_mu_z <- -((1 - pi) * prob_nb_0 / prob_zero) * (theta * mu) / (theta + mu)
      grad_eta_mu <- ifelse(y == 0, grad_eta_mu_z, grad_eta_mu_nz)

      # Gradient wrt Zero Model (pi) ----
      grad_eta_pi_nz <- -pi
      grad_eta_pi_z <- (1 - prob_nb_0) / prob_zero * pi * (1 - pi)
      grad_eta_pi <- ifelse(y == 0, grad_eta_pi_z, grad_eta_pi_nz)

      # Gradient wrt Theta ----
      dL_dtheta_nz <- digamma(y + theta) - digamma(theta) +
        log(theta) - log(theta + mu) +
        1 - (y + theta) / (mu + theta)

      d_pnb0_dtheta <- prob_nb_0 * (log(theta / (theta + mu)) + mu / (theta + mu))
      dL_dtheta_z <- ((1 - pi) / prob_zero) * d_pnb0_dtheta

      dL_dtheta <- ifelse(y == 0, dL_dtheta_z, dL_dtheta_nz)
      grad_zeta <- dL_dtheta * theta

      # Map gradients back to unified vector ----
      g_full <- numeric(length(beta))

      g_full[comps$idx_mu] <- -crossprod(X[, comps$idx_mu, drop = FALSE], grad_eta_mu) / n
      g_full[comps$idx_pi] <- -crossprod(X[, comps$idx_pi, drop = FALSE], grad_eta_pi) / n
      g_full[length(beta)] <- -sum(grad_zeta) / n

      g_full
    },
    predict = function(X, beta, type = "response") {
      comps <- grad_zinb()$split_components(X, beta)

      if (type == "link") {
        return(cbind(count = comps$eta_mu, zero = comps$eta_pi))
      }

      mu <- exp(comps$eta_mu)
      pi <- 1 / (1 + exp(-comps$eta_pi))

      if (type == "count") {
        return(mu)
      }
      if (type == "zero") {
        return(pi)
      }
      if (type == "response") {
        return((1 - pi) * mu)
      }
    }
  )
}
