test_that("pd.adjust handles basic logic and math", {
  # Use q=0.7, m=1 so H0 = 0.7 (which is > 0.5, allowing adjustment)
  pd_val <- 0.8
  q_val <- 0.7
  m_val <- 1
  
  res <- pd.adjust(pd = pd_val, q = q_val, m = m_val)
  
  # Manual calculation:
  # H0 = 0.7, H1 = 0.3
  # pd_adj = (0.8 * 0.3) / (0.8 * 0.3 + 0.2 * 0.7) 
  # pd_adj = 0.24 / (0.24 + 0.14) = 0.24 / 0.38 ≈ 0.63157
  expect_equal(res$pd_adj, 0.24 / 0.38, tolerance = 1e-5)
})

test_that("pd.adjust Cheverud correction works", {
  pd_vals <- c(0.9, 0.9, 0.9)
  
  # Identity matrix: m_eff should equal length(pd)
  res_ident <- pd.adjust(pd = pd_vals, R = diag(3), q = 0.8)
  expect_equal(res_ident$m[1], 3)
  
  # High correlation: m_eff should be lower than 3
  R_high <- matrix(0.95, 3, 3); diag(R_high) <- 1
  res_corr <- pd.adjust(pd = pd_vals, R = R_high, q = 0.8)
  expect_lt(res_corr$m[1], 3)
})

test_that("pd.adjust catches invalid inputs and triggers warnings", {
  # Should error: pd must be >= 0.5
  expect_error(pd.adjust(pd = 0.4)) 
  
  # Should warn: q=0.4, m=1 means H0=0.4 (which is <= 0.5)
  expect_warning(pd.adjust(pd = 0.9, q = 0.4, m = 1), "Non-conservative prior")
})