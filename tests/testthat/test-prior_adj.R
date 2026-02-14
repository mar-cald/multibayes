test_that("pd.adjust adjusts pd values correctly", {
  result <- pd.adjust(c(0.9, 0.95), q = 0.5)
  expect_length(result, 2)
  expect_true(all(result >= 0.5 & result <= 1))
})

test_that("pd.adjust validates inputs", {
  expect_error(pd.adjust(c(0.3, 0.9), q = 0.5))
  expect_error(pd.adjust(c(0.9), q = 1.5))
})

