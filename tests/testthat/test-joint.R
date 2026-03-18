test_that("joint credible bands are wider than marginal", {
  set.seed(42)
  draws <- matrix(rnorm(3000), ncol = 3)
  
  res_joint    <- joint(draws, prob = 0.95)
  res_marginal <- apply(draws, 2, quantile, probs = c(0.025, 0.975))
  
  # joint lower bounds must be at or below marginal lower bounds
  expect_true(all(res_joint$lower <= res_marginal[1, ]))
  # joint upper bounds must be at or above marginal upper bounds
  expect_true(all(res_joint$upper >= res_marginal[2, ]))
})

test_that("joint returns correct structure", {
  set.seed(42)
  draws <- matrix(rnorm(1000), ncol = 4)
  colnames(draws) <- paste0("theta", 1:4)
  res <- joint(draws, prob = 0.95)
  
  expect_named(res, c("lower", "upper", "prob", "cq", "est", "est.FUN"))
  expect_length(res$lower, 4)
  expect_length(res$upper, 4)
  expect_true(res$cq < 0.025)          #
  expect_true(all(res$lower < res$upper))  # bounds are ordered
})
