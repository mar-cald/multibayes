test_that("pd.adjust() errors cleanly for invalid p", {
  expect_error(
    pd.adjust(pd = 0.9, p = 0),
    regexp = "`p`: must be a single number in \\(0, 1\\)"
  )
  expect_error(
    pd.adjust(pd = 0.9, p = 1),
    regexp = "`p`: must be a single number in \\(0, 1\\)"
  )
  expect_error(
    pd.adjust(pd = 0.9, p = c(0.4, 0.5)),
    regexp = "`p`: must be a single number in \\(0, 1\\)"
  )
})

test_that("pd.adjust() errors for invalid pd range", {
  expect_error(
    pd.adjust(pd = c(0.4, 0.6), p = 0.4),
    "`pd` must be numeric in \\[0.5, 1\\]."
  )
  expect_error(
    pd.adjust(pd = c(0.6, 1.1), p = 0.4),
    "`pd` must be numeric in \\[0.5, 1\\]."
  )
})

test_that("pd.adjust() ignores direction when pd is supplied", {
  expect_warning(
    res <- pd.adjust(pd = c(0.9, 0.8), p = 0.4, direction = "greater"),
    "`direction` cannot be specified, fixed to `two.sided`"
  )
  expect_equal(res$direction, c("two.sided", "two.sided"))
})

test_that("pd.adjust() returns unadjusted pd when prior_H0 <= 0.5", {
  # p small enough that p^(1/m) <= 0.5
  pd  <- c(0.9, 0.8, 0.7, 0.6)
  p   <- 0.01
  out <- expect_warning(
    pd.adjust(pd = pd, p = p),
    "Pr\\(H0_i\\) <= 0.5 \\(non-conservative prior\\); returning unadjusted pd\\."
  )
  expect_equal(out$pd, out$pd.adj)
})

test_that("pd.adjust() shrinks pd toward 0.5 for two.sided tests", {
  pd  <- c(0.6, 0.8, 0.95)
  p   <- 0.4
  out <- pd.adjust(pd = pd, p = p)
  
  # adjusted pd should be strictly between 0.5 and original pd
  expect_true(all(out$pd.adj >= 0.5))
  expect_true(all(out$pd.adj <= out$pd))
  # at least one should actually shrink
  expect_true(any(out$pd.adj < out$pd))
})

test_that("pd.adjust() with draws reproduces pd logic for directions", {
  set.seed(1)
  Sigma <- diag(3)
  mu    <- c(1, -1, 0)  # clearly positive, clearly negative, centered
  draws <- MASS::mvrnorm(n = 5000, mu = mu, Sigma = Sigma)
  colnames(draws) <- c("H_pos", "H_neg", "H_twosided")
  
  out <- pd.adjust(
    draws     = draws,
    p         = 0.4,
    null.value= 0,
    direction = c("greater", "less", "two.sided")
  )
  
  # H_pos: pd ~ P(theta>0) should be close to mean(draws[,1]>0)
  pd_manual_pos <- mean(draws[, 1] > 0)
  expect_equal(out$pd[ out$direction == "greater" ],
               pd_manual_pos,
               tolerance = 0.02)
  
  # H_neg: pd ~ P(theta<0) should be close to mean(draws[,2]<0)
  pd_manual_neg <- mean(draws[, 2] < 0)
  expect_equal(out$pd[ out$direction == "less" ],
               pd_manual_neg,
               tolerance = 0.02)
  
  # H_twosided: pd ~ max(P(theta>0), P(theta<0)) for 3rd column
  p_gt0 <- mean(draws[, 3] > 0)
  p_lt0 <- mean(draws[, 3] < 0)
  pd_manual_two <- max(p_gt0, p_lt0)
  expect_equal(out$pd[ out$direction == "two.sided" ],
               pd_manual_two,
               tolerance = 0.02)
})

test_that("pd.adjust() floors two.sided pd.adj at 0.5", {
  set.seed(2)
  draws <- matrix(rnorm(4000 * 2, mean = 0, sd = 1), ncol = 2)
  colnames(draws) <- c("H1", "H2")
  
  out <- expect_warning(
    pd.adjust(draws = draws, p = 0.9, null.value = 0,
              direction = c("two.sided", "two.sided")),
    "some pd.adj have been floored to 0.5\\."
  )
  
  expect_true(all(out$pd.adj >= 0.5))
})

test_that("pd.adjust() handles mixed directions and parameter-specific nulls", {
  set.seed(3)
  Sigma <- diag(4)
  mu    <- c(1, -0.5, 0, 0.5)
  draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
  colnames(draws) <- paste0("H", 1:4)
  
  out <- pd.adjust(
    draws      = draws,
    p          = 0.3,
    null.value = c(0, 0, 0.2, 0.2),
    direction  = c("greater", "less", "two.sided", "greater")
  )
  
  expect_equal(nrow(out), 4L)
  expect_equal(out$direction,
               c("greater", "less", "two.sided", "greater"))
  expect_equal(out$null.value, c(0, 0, 0.2, 0.2))
  expect_true(all(out$pd.adj >= 0))        # valid range
  expect_true(all(out$pd.adj <= 1))
})

test_that("pd.adjust() reports prior pieces correctly", {
  pd <- c(0.9, 0.8, 0.7)
  p  <- 0.4
  out <- pd.adjust(pd = pd, p = p)
  
  # m and pm should reflect number of tests and per-test prior
  expect_true(all(out$m == length(pd)))
  expect_true(all(out$pm == round(p^(1/length(pd)),4)))
  expect_true(all(out$p  == p))
})