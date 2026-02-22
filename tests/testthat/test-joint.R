test_that("joint credible bands are wider than marginal", {
  set.seed(42)
  draws <- matrix(rnorm(3000), ncol = 3)
  
  res_joint <- joint(draws, alpha = 0.05)
  
  # Marginal 95%
  res_marginal <- t(apply(draws, 2, quantile, probs = c(0.025, 0.975)))
  
  # Joint must be wider (lower is lower, upper is higher)
  expect_true(all(res_joint[, "lower"] <= res_marginal[, 1]))
  expect_true(all(res_joint[, "upper"] >= res_marginal[, 2]))
})

test_that("joint handles the resolution warning", {
  small_draws <- matrix(rnorm(20), ncol = 2) 
  
  # This regex will catch the warning regardless of case or punctuation
  expect_warning(
    joint(small_draws, alpha = 0.001), 
    regexp = "smaller than.*resolution"
  )
})