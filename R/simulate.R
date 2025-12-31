#' Simulate Data for gradLasso
#'
#' Generates synthetic data for Gaussian, Binomial, Negative Binomial, or ZINB models
#' with correlated predictors.
#'
#' @param n Number of observations.
#' @param p Number of predictors.
#' @param family Model family: "gaussian", "binomial", "negbin", or "zinb".
#' @param rho Correlation coefficient between predictors (Toeplitz structure).
#' @param k Number of non-zero coefficients (sparsity) for single-part models.
#' @param k_mu Number of non-zero coefficients for Count part (ZINB only).
#' @param k_pi Number of non-zero coefficients for Zero part (ZINB only).
#' @param theta Dispersion parameter for NegBin and ZINB.
#' @param intercept_mu Intercept for main model (or count part).
#' @param intercept_pi Intercept for zero-inflation part.
#' @param snr Signal-to-noise ratio (Gaussian only).
#' @return A list containing the following components:
#' \item{X}{A matrix of predictor variables with induced correlation.}
#' \item{y}{A vector of the simulated response variable.}
#' \item{family}{The family string used for simulation.}
#' \item{truth}{A list containing the true parameters used to generate the data (e.g., \code{beta}, \code{theta}, \code{sigma}).}
#' @export
simulate_data <- function(n = 1000, p = 20, family = "gaussian", rho = 0.2,
                          k = 5, k_mu = 5, k_pi = 5, theta = 1.0,
                          intercept_mu = 0, intercept_pi = -1.0, snr = 3) {
  family <- match.arg(family, c("gaussian", "binomial", "negbin", "zinb"))

  # Generate correlated predictors (X) ----
  # Construct Toeplitz covariance matrix
  Sigma <- matrix(0, p, p)
  for (i in 1:p) for (j in 1:p) Sigma[i, j] <- rho^abs(i - j)

  # Induce correlation via Cholesky decomposition
  L <- chol(Sigma)
  X <- matrix(rnorm(n * p), n, p) %*% L
  colnames(X) <- paste0("Var", 1:p)

  # Helper to generate sparse coefficient vectors ----
  create_beta <- function(k_active, shift = 0) {
    b <- rep(0, p)
    if (k_active > 0) {
      # Use alternating signs and varying magnitudes for realism
      vals <- rep_len(c(0.5, -0.5, 0.7, -0.3, 0.4, -0.4), k_active)

      # Shift indices to ensure different models select different variables
      idx <- ((0:(k_active - 1) + shift) %% p) + 1
      b[idx] <- vals
    }
    return(b)
  }

  y <- NULL
  truth <- list()

  # Generate response variable (y) ----
  if (family == "gaussian") {
    # Gaussian family
    beta <- create_beta(k)
    signal <- X %*% beta

    # Scale noise to achieve target Signal-to-Noise Ratio (SNR)
    sigma_err <- sqrt(var(signal) / snr)
    y <- intercept_mu + signal + rnorm(n, sd = sigma_err)
    truth <- list(beta = beta, sigma = sigma_err)
  } else if (family == "binomial") {
    # Binomial family
    beta <- create_beta(k)
    eta <- intercept_mu + X %*% beta
    probs <- 1 / (1 + exp(-eta))
    y <- rbinom(n, 1, probs)
    truth <- list(beta = beta)
  } else if (family == "negbin") {
    # Negative Binomial family
    beta <- create_beta(k)
    eta <- intercept_mu + X %*% beta

    # Clip linear predictor for numerical stability
    mu <- exp(pmin(eta, 20))
    y <- rnbinom(n, size = theta, mu = mu)
    truth <- list(beta = beta, theta = theta)
  } else if (family == "zinb") {
    # Zero-Inflated Negative Binomial (ZINB) family
    beta_mu <- create_beta(k_mu, shift = 0)

    # Shift zero model variables so they differ from count model variables
    beta_pi <- create_beta(k_pi, shift = 2)

    eta_mu <- intercept_mu + X %*% beta_mu
    eta_pi <- intercept_pi + X %*% beta_pi

    mu <- exp(pmin(eta_mu, 20))
    pi <- 1 / (1 + exp(-eta_pi))

    # Mixture process: structural zeros vs negative binomial counts
    y <- ifelse(rbinom(n, 1, pi) == 1, 0, rnbinom(n, size = theta, mu = mu))
    truth <- list(beta_mu = beta_mu, beta_pi = beta_pi, theta = theta)
  }

  return(list(X = X, y = drop(y), family = family, truth = truth))
}
