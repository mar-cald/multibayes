test_that("pd.adjust() errors cleanly for invalid p0", {
  expect_error(
    pd.adjust(pd = 0.9, p0 = 0),
    regexp = "`p0`: must be a single number in \\(0, 1\\)"
  )
  expect_error(
    pd.adjust(pd = 0.9, p0 = 1),
    regexp = "`p0`: must be a single number in \\(0, 1\\)"
  )
  expect_error(
    pd.adjust(pd = 0.9, p0 = c(0.4, 0.5)),
    regexp = "`p0`: must be a single number in \\(0, 1\\)"
  )
})

test_that("pd.adjust() errors for invalid pd range", {
  expect_error(
    pd.adjust(pd = c(0.4, 0.6), p0 = 0.4),
    "`pd` must be numeric in \\[0.5, 1\\]."
  )
  expect_error(
    pd.adjust(pd = c(0.6, 1.1), p0 = 0.4),
    "`pd` must be numeric in \\[0.5, 1\\]."
  )
})

test_that("pd.adjust() ignores direction when pd is supplied", {
  expect_warning(
    res <- pd.adjust(pd = c(0.9, 0.8), p0 = 0.4, direction = "greater"),
    "`direction` cannot be specified, fixed to `two.sided`"
  )
  expect_equal(res$direction, c("two.sided", "two.sided"))
})

test_that("pd.adjust() returns unadjusted pd when prior_H0 <= 0.5", {
  # p0 small enough that p0^(1/m) <= 0.5
  pd  <- c(0.9, 0.8, 0.7, 0.6)
  p0  <- 0.01
  out <- expect_warning(
    pd.adjust(pd = pd, p0 = p0),
    "non-conservative prior"
  )
  expect_equal(out$pd, out$pd.adj)
})

test_that("pd.adjust() shrinks pd toward 0.5 for two.sided tests", {
  pd  <- c(0.6, 0.8, 0.95)
  p0  <- 0.4
  out <- pd.adjust(pd = pd, p0 = p0)

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
    p0        = 0.4,
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
    pd.adjust(draws = draws, p0 = 0.9, null.value = 0,
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
    p0         = 0.3,
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
  p0 <- 0.4
  out <- pd.adjust(pd = pd, p0 = p0)

  # m and pi0 should reflect number of tests and per-test prior
  expect_true(all(out$m == length(pd)))
  expect_true(all(out$pi0 == round(p0^(1/length(pd)), 4)))
  # the global null p0 is intentionally not returned
  expect_false("p0" %in% names(out))
})

test_that("pd.adjust() accepts pi0 directly and matches the p0 interface", {
  pd  <- c(0.99, 0.9, 0.8, 0.7)
  m   <- length(pd)
  p0  <- 0.4
  pi0 <- p0^(1 / m)

  from_p0  <- suppressWarnings(pd.adjust(pd = pd, p0 = p0))
  from_pi0 <- suppressWarnings(pd.adjust(pd = pd, pi0 = pi0))

  # identical adjustment regardless of which input is used
  expect_equal(from_pi0$pd.adj, from_p0$pd.adj)
  # per-test prior is reported consistently (global null p0 is not returned)
  expect_equal(from_pi0$pi0, round(rep(pi0, m), 4))
})

test_that("pd.adjust() requires exactly one of p0 or pi0", {
  pd <- c(0.9, 0.8)
  expect_error(pd.adjust(pd = pd), "exactly one")
  expect_error(pd.adjust(pd = pd, p0 = 0.4, pi0 = 0.8), "exactly one")
  expect_error(pd.adjust(pd = pd, pi0 = 1.2), "must be a single number")
})

test_that("pd.adjust() returns a classed object and keeps $ access", {
  out <- suppressWarnings(pd.adjust(pd = c(0.9, 0.8, 0.7), p0 = 0.4))
  expect_s3_class(out, "pd_adjust")
  expect_s3_class(out, "data.frame")
  expect_true(is.numeric(out$pd.adj))      # column access still works
  expect_true(all(c("pd", "pd.adj", "pi0", "m") %in% names(out)))
  expect_false("p0" %in% names(out))        # global null not returned
})

test_that("pd.adjust() adds joint_cum only when joint = TRUE", {
  pd <- c(0.99, 0.9, 0.8)
  expect_null(suppressWarnings(pd.adjust(pd = pd, p0 = 0.4))$joint_cum)
  out <- suppressWarnings(pd.adjust(pd = pd, p0 = 0.4, joint = TRUE))
  expect_true("joint_cum" %in% names(out))
  # pd-only: joint_cum is the running product of pd (input order)
  expect_equal(out$joint_cum, round(cumprod(pd), 4))
})

test_that("pd.adjust() joint from draws is the exact intersection (dependence-aware)", {
  set.seed(10)
  Sigma <- matrix(0.4, 4, 4); diag(Sigma) <- 1
  draws <- MASS::mvrnorm(20000, mu = c(1, 1, 1, 1), Sigma = Sigma)
  colnames(draws) <- paste0("H", 1:4)

  out <- suppressWarnings(pd.adjust(draws = draws, p0 = 0.4, joint = TRUE))
  side <- draws > 0                                  # all directional claims are theta > 0 here
  manual <- mean(rowSums(side) == 4)
  expect_equal(unname(out$joint_cum[4]), round(manual, 4))
  # under positive dependence the joint exceeds the independence product
  expect_gt(out$joint_cum[4], prod(out$pd) - 1e-8)
})

test_that("print.pd_adjust returns invisibly and prints a header", {
  out <- suppressWarnings(pd.adjust(pd = c(0.9, 0.8), p0 = 0.4, joint = TRUE))
  expect_output(print(out), "Prior-odds adjusted probability of direction")
  expect_output(print(out), "joint_cum")
  expect_identical(suppressWarnings(print(out)), out)
})
