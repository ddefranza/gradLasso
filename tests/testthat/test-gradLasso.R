test_that("gradLasso fits a simple gaussian model", {
  set.seed(123)
  sim <- simulate_data(n = 100, p = 5, family = "gaussian")

  df <- data.frame(y = sim$y, sim$X)
  colnames(df) <- c("y", paste0("Var", 1:5))

  expect_equal(ncol(df), 6) # y + 5 Vars

  fit <- gradLasso(y ~ ., data = df, boot = FALSE)

  expect_s3_class(fit, "gradLasso")

  expect_equal(length(coef(fit)), 6)

  expect_true("(Intercept)" %in% names(coef(fit)))
})

test_that("gradLasso fits a ZINB model with pipe", {
  set.seed(123)
  sim <- simulate_data(n = 200, p = 5, family = "zinb")

  df <- data.frame(y = sim$y, sim$X)
  colnames(df) <- c("y", paste0("Var", 1:5))

  fit <- gradLasso(y ~ . | ., data = df, family = grad_zinb(), lambda = 0.1, boot = FALSE)

  expect_s3_class(fit, "gradLasso")

  c_names <- names(coef(fit))

  expect_true(any(grepl("count_Var1", c_names)))
  expect_true(any(grepl("zero_Var1", c_names)))
})

test_that("S3 methods work (print, summary, logLik)", {
  sim <- simulate_data(n = 50, p = 3, family = "gaussian")
  df <- data.frame(y = sim$y, sim$X)
  colnames(df) <- c("y", paste0("Var", 1:3))

  fit <- gradLasso(y ~ ., data = df, boot = FALSE)

  expect_output(print(fit), "gradLasso Fitted Object")

  s <- summary(fit)
  expect_s3_class(s, "summary.gradLasso")
  expect_output(print(s), "gradLasso Model Summary")

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_type(as.numeric(ll), "double")
})
