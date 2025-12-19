## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gradLasso)

## -----------------------------------------------------------------------------
set.seed(42)
# Simulate 200 obs, 20 predictors, 5 active
sim <- simulate_data(n = 200, p = 20, family = "gaussian", k = 5, snr = 3.0)
df <- data.frame(y = sim$y, sim$X)

# Check the first few rows
head(df[, 1:6])

## -----------------------------------------------------------------------------
fit <- gradLasso(y ~ ., data = df, lambda_cv = TRUE, boot = TRUE, n_boot = 50)

print(fit)

## -----------------------------------------------------------------------------
summary(fit)

