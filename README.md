# gradLasso

[![R-CMD-check](https://github.com/ddefranza/gradLasso/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ddefranza/gradLasso/actions)
**gradLasso** implements an efficient gradient descent solver for LASSO-penalized regression models. It is designed to handle complex Generalized Linear Models (GLMs) and large $N$ data including:

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
