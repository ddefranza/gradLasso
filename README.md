# gradLasso

[![R-CMD-check](https://github.com/ddefranza/gradLasso/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ddefranza/gradLasso/actions)
**gradLasso** implements an efficient gradient descent solver for LASSO-penalized regression models. It is designed to handle large $N$ data and complex Generalized Linear Models (GLMs) including:

* **Gaussian** (Standard Least Squares)
* **Binomial** (Logistic Regression)
* **Negative Binomial** (Overdispersed counts)
* **Zero-Inflated Negative Binomial (ZINB)** (Excess zeros + overdispersion)

Key features include built-in **Stability Selection** for robust variable identification and **Cross-Validation** for automatic lambda tuning.

## Installation

You can install the development version of gradLasso from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("ddefranza/gradLasso")
```
## Quick Start

## 1. Gaussian Regression

Simulate data and fit a standard LASSO model with stability selection (bootstrapping).

```r
library(gradLasso)

# Simulate Gaussian data
set.seed(123)
sim <- simulate_data(n = 200, p = 20, family = "gaussian", k = 5, snr = 3.0)
df <- data.frame(y = sim$y, sim$X)

# Fit model with CV and 50 bootstraps
fit <- gradLasso(y ~ ., data = df, lambda_cv = TRUE, boot = TRUE, n_boot = 50)

# Print summary
summary(fit)
```

### Diagnostic Plots:

```r
# Plot Stability Selection and CV Deviance
plot(fit, which = c(1, 2))
```

## 2. Zero-Inflated Negative Binomial (ZINB)

`gradLasso` supports a pipe syntax (`|`) to specify separate predictors for the Count model (Negative Binomial) and the Zero-Inflation model (Logit).

```r
# Simulate ZINB data
set.seed(456)
sim_zinb <- simulate_data(n = 500, p = 20, family = "zinb", 
                          k_mu = 5, k_pi = 5, theta = 2.0)
df_zinb <- data.frame(y = sim_zinb$y, sim_zinb$X)

# Fit ZINB model
# Syntax: y ~ count_predictors | zero_predictors
fit_zinb <- gradLasso(y ~ . | ., data = df_zinb, 
                      family = grad_zinb(), 
                      n_boot = 10,
                      lambda = 0.05)

summary(fit_zinb)
```

## Features

* Parallel Processing: Automatically parallelize Cross-Validation and Bootstrapping using `parallel = TRUE` and `n_cores = X`.
* Sparse Solutions: Uses Proximal Gradient Descent (FISTA) to ensure true sparsity in coefficients.
* Tidy Summaries: `summary()` returns a clean data frame of coefficients, selection probabilities, and confidence intervals.

## Citation

If you use `gradLasso` in your research, please cite it as:

> DeFranza, D. (2025). gradLasso: Gradient Descent LASSO with Stability Selection. R package version 0.1.0. https://github.com/ddefranza/gradLasso

## License

MIT

