test_that("prior_adj adjusts pd values correctly", {
  result <- prior_adj(c(0.9, 0.95), q = 0.5)
  expect_length(result, 2)
  expect_true(all(result >= 0.5 & result <= 1))
})

test_that("prior_adj validates inputs", {
  expect_error(prior_adj(c(0.3, 0.9), q = 0.5))
  expect_error(prior_adj(c(0.9), q = 1.5))
})

