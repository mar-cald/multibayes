test_that("joint() validates its inputs", {
  d <- matrix(rnorm(2000 * 3), ncol = 3)
  expect_error(joint(d[1:10, ]), "at least 1000 rows")
  expect_error(joint(matrix(rnorm(2000), ncol = 1)), "at least 2 columns")
  expect_error(joint(d, interval = "yes"), "must be TRUE or FALSE")
  expect_error(joint(d, interval = TRUE, prob = 0), "must be a single number in \\(0, 1\\)")
  expect_error(joint(d, interval = TRUE, prob = 1), "must be a single number in \\(0, 1\\)")
  expect_error(joint(d, interval = TRUE, est.FUN = 42), "must be a function")
})

# ---- Cumulative joint statement (default) ----------------------------------

test_that("joint() default returns the cumulative joint statement", {
  set.seed(1)
  d <- MASS::mvrnorm(3000, mu = c(2.5, 1, -1.5, 0.2), Sigma = diag(4))
  colnames(d) <- c("a", "b", "c", "e")
  out <- joint(d, null.value = 0)

  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("median.est", "null.value", "pd", "direction", "joint_cum"))
  expect_equal(nrow(out), 4L)
  # ordered by decreasing pd, and joint_cum is non-increasing down the rows
  expect_false(is.unsorted(rev(out$pd)))
  expect_false(is.unsorted(rev(out$joint_cum)))
  # first joint_cum equals the strongest pd
  expect_equal(out$joint_cum[1], out$pd[1])
})

test_that("joint() cumulative is the exact intersection (dependence-aware)", {
  set.seed(10)
  Sigma <- matrix(0.4, 4, 4); diag(Sigma) <- 1
  draws <- MASS::mvrnorm(20000, mu = c(1, 1, 1, 1), Sigma = Sigma)
  colnames(draws) <- paste0("H", 1:4)

  out    <- joint(draws, null.value = 0)
  manual <- mean(rowSums(draws > 0) == 4)            # all four on the theta>0 side
  expect_equal(unname(out$joint_cum[4]), round(manual, 4))
  # positive dependence: the joint exceeds the independence product of the pd
  expect_gt(out$joint_cum[4], prod(out$pd) - 1e-8)
})

test_that("joint() cumulative honours direction and validates lengths", {
  set.seed(2)
  draws <- MASS::mvrnorm(3000, mu = c(1, -1, 0.5), Sigma = diag(3))
  colnames(draws) <- c("p1", "p2", "p3")
  out <- joint(draws, direction = c("greater", "less", "two.sided"))
  expect_equal(nrow(out), 3L)
  expect_error(joint(draws, direction = c("greater", "less")), "scalar or a vector")
  expect_error(joint(draws, null.value = c(0, 0)), "scalar or a vector")
})

# ---- Simultaneous credible intervals (interval = TRUE) ---------------------

test_that("joint(interval = TRUE) returns simultaneous intervals", {
  set.seed(1)
  d <- MASS::mvrnorm(3000, mu = c(2, 0, -1), Sigma = diag(3))
  colnames(d) <- c("a", "b", "c")
  out <- joint(d, interval = TRUE, prob = 0.95)

  expect_equal(names(out), c("lower", "est", "upper", "prob", "cq"))
  expect_equal(rownames(out), c("a", "b", "c"))
  expect_true(all(out$lower < out$upper))
  expect_true(all(out$prob == 0.95))
})

test_that("joint(interval = TRUE) bands are wider than marginal and calibrated", {
  set.seed(2)
  p <- 5
  Sigma <- 0.4^abs(outer(seq_len(p), seq_len(p), "-")); diag(Sigma) <- 1
  mu    <- c(3, 1, 0, -1, 0.5)
  cal   <- MASS::mvrnorm(20000, mu, Sigma)
  colnames(cal) <- paste0("t", seq_len(p))

  band <- joint(cal, interval = TRUE, prob = 0.95)
  marg <- t(apply(cal, 2, stats::quantile, probs = c(0.025, 0.975)))
  expect_true(all((band$upper - band$lower) >= (marg[, 2] - marg[, 1]) - 1e-8))

  fresh  <- MASS::mvrnorm(50000, mu, Sigma)
  inside <- rowSums(t(t(fresh) >= band$lower & t(fresh) <= band$upper)) == p
  expect_equal(mean(inside), 0.95, tolerance = 0.02)
})

test_that("joint(interval = TRUE) accepts a data frame and a custom estimator", {
  set.seed(3)
  d   <- as.data.frame(MASS::mvrnorm(2000, mu = c(1, -1), Sigma = diag(2)))
  out <- joint(d, interval = TRUE, prob = 0.90, est.FUN = mean)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$prob == 0.90))
})
