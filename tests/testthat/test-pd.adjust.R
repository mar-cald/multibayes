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
    "`direction` is ignored when `pd` is supplied"
  )
  expect_equal(res$direction, c("two.sided", "two.sided"))
})

test_that("pd.adjust() leaves the threshold unchanged when pi0 <= 0.5", {
  # p0 small enough that p0^(1/m) <= 0.5
  pd  <- c(0.9, 0.8, 0.7, 0.6)
  out <- expect_warning(
    pd.adjust(pd = pd, p0 = 0.01),
    "no penalty is applied"
  )
  expect_equal(attr(out, "threshold.adj"), attr(out, "threshold"))
})

test_that("pd.adjust() reports the threshold c*", {
  pd  <- c(0.6, 0.8, 0.95)
  out <- pd.adjust(pd = pd, p0 = 0.4)

  expect_equal(out$pd, pd)
  expect_gt(attr(out, "threshold.adj"), attr(out, "threshold"))
})

test_that("pd.adj is the posterior under the most skeptical prior and matches c*", {
  pi0 <- 0.7; c0 <- 0.975
  pd  <- c(0.6, 0.97, 0.98, 0.985, 0.99, 0.999)
  out <- pd.adjust(pd = pd, pi0 = pi0, threshold = c0)
  expected <- pd * (1 - pi0) / (pd * (1 - pi0) + (1 - pd) * pi0)
  expect_equal(out$pd.adj, round(expected, 4))
  # posterior odds under the skeptical prior = Bayes factor x prior odds
  expect_equal(expected / (1 - expected), pd / (1 - pd) * (1 - pi0) / pi0)
  # pd.adj >= c if and only if pd >= c*
  expect_identical(expected >= c0, pd >= attr(out, "threshold.adj"))
  # no floor: values below .5 are kept
  expect_lt(out$pd.adj[1], 0.5)
  # pi0 = .5 leaves pd unchanged
  expect_equal(pd.adjust(pd = pd, pi0 = 0.5)$pd.adj, round(pd, 4))
  txt <- capture.output(print(out))
  expect_true(any(grepl("pd.adj: Pr(claimed direction | y)", txt, fixed = TRUE)))
  expect_true(any(grepl("Claim a direction when pd.adj >= 0.975", txt, fixed = TRUE)))
})

test_that("the threshold reproduces the decision under the most skeptical prior of the class", {
  # Penalising the odds by (1 - pi0)/pi0 equals reweighting the draws so the
  # prior mass either side of the null moves to (pi0, 1 - pi0). The reported
  # threshold must reproduce the decision taken under that shifted prior.
  set.seed(7)
  pi0 <- 0.8
  mus <- c(0.5, 1.5, 2.0, 2.4, 2.8, 3.5)
  dr  <- sapply(mus, function(mu) rnorm(2e4, mu, 1))
  colnames(dr) <- paste0("H", seq_along(mus))

  out <- pd.adjust(draws = dr, pi0 = pi0, direction = "greater")
  shifted <- apply(dr, 2, function(x) {
    w <- ifelse(x > 0, 2 * (1 - pi0), 2 * pi0)
    sum(w * (x > 0)) / sum(w)
  })

  expect_equal(unname(shifted > attr(out, "threshold")),
               unname(colMeans(dr > 0) > attr(out, "threshold.adj")))
  # the grid must contain both outcomes for the check to mean anything
  expect_true(any(shifted > attr(out, "threshold")))
  expect_true(any(shifted <= attr(out, "threshold")))
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

test_that("directional pd is the one-sided mass and may fall below 0.5", {
  set.seed(2)
  draws <- matrix(rnorm(4000 * 2, mean = -0.3, sd = 1), ncol = 2)
  colnames(draws) <- c("H1", "H2")

  out <- pd.adjust(draws = draws, p0 = 0.9, null.value = 0,
                   direction = c("greater", "less"))

  expect_equal(out$pd[1], round(mean(draws[, 1] > 0), 4))
  expect_equal(out$pd[2], round(mean(draws[, 2] < 0), 4))
  expect_lt(out$pd[1], 0.5)        # posterior against the predicted side
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
  expect_true(all(out$pd >= 0 & out$pd <= 1))   # valid range
})

test_that("pd.adjust() reports prior pieces correctly", {
  pd <- c(0.9, 0.8, 0.7)
  p0 <- 0.4
  out <- pd.adjust(pd = pd, p0 = p0)

  # m and pi0 should reflect the number of tests and the per-test value pi0
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

  # identical threshold regardless of which input is used
  expect_equal(attr(from_pi0, "threshold.adj"), attr(from_p0, "threshold.adj"))
  # pi0 is reported consistently (p0 is not returned)
  expect_equal(from_pi0$pi0, round(rep(pi0, m), 4))
})

test_that("pd.adjust() requires exactly one of p0 or pi0", {
  pd <- c(0.9, 0.8)
  expect_error(pd.adjust(pd = pd), "exactly one")
  expect_error(pd.adjust(pd = pd, p0 = 0.4, pi0 = 0.8), "exactly one")
  expect_error(pd.adjust(pd = pd, pi0 = 1.2), "must be a single number")
})

test_that("pd.adjust() reports the adjusted threshold and effective alpha", {
  pd  <- c(0.999, 0.95, 0.80)
  pi0 <- 0.7
  out <- pd.adjust(pd = pd, pi0 = pi0)

  cc  <- 0.975
  cst <- (cc * pi0) / ((1 - pi0) * (1 - cc) + cc * pi0)   # closed-form c_*
  expect_equal(attr(out, "threshold"),     0.975)
  expect_equal(attr(out, "threshold.adj"), cst)
  expect_equal(attr(out, "alpha.nominal"), 2 * (1 - cc))
  expect_equal(attr(out, "alpha.eff"),     2 * (1 - cst))

  # equivalence: penalized odds >= c/(1 - c)  <=>  pd >= c*
  shifted <- out$pd * (1 - pi0) / (out$pd * (1 - pi0) + (1 - out$pd) * pi0)
  expect_equal(shifted >= 0.975, out$pd >= cst)
})

test_that("pd.adjust() honours a custom threshold and validates it", {
  pd  <- c(0.999, 0.95, 0.80)
  out <- pd.adjust(pd = pd, pi0 = 0.7, threshold = 0.95)
  cst <- (0.95 * 0.7) / ((1 - 0.7) * (1 - 0.95) + 0.95 * 0.7)
  expect_equal(attr(out, "threshold.adj"), cst)
  shifted <- out$pd * (1 - 0.7) / (out$pd * (1 - 0.7) + (1 - out$pd) * 0.7)
  expect_equal(shifted >= 0.95, out$pd >= cst)

  expect_error(pd.adjust(pd = pd, pi0 = 0.7, threshold = 0.4),
               "`threshold`: must be a single number")
  expect_error(pd.adjust(pd = pd, pi0 = 0.7, threshold = c(0.9, 0.95)),
               "`threshold`: must be a single number")
})

test_that("pd.adjust() leaves the threshold unchanged for pi0 <= 0.5", {
  pd  <- c(0.9, 0.8, 0.7, 0.6)
  out <- suppressWarnings(pd.adjust(pd = pd, pi0 = 0.4))   # pi0 <= 0.5: no adjustment
  expect_equal(attr(out, "threshold.adj"), 0.975)
  expect_equal(attr(out, "alpha.eff"), attr(out, "alpha.nominal"))
})

test_that("pd.adjust() accepts alpha as an alternative to threshold", {
  pd <- c(0.999, 0.95, 0.80)

  # alpha = 0.05 is equivalent to the default threshold = 0.975
  by_alpha <- pd.adjust(pd = pd, pi0 = 0.7, alpha = 0.05)
  by_thr   <- pd.adjust(pd = pd, pi0 = 0.7, threshold = 0.975)
  expect_equal(attr(by_alpha, "threshold"),     attr(by_thr, "threshold"))
  expect_equal(attr(by_alpha, "threshold.adj"), attr(by_thr, "threshold.adj"))
  expect_equal(attr(by_alpha, "alpha.nominal"), 0.05)

  # alpha overrides threshold:  c = 1 - alpha/2
  ov <- pd.adjust(pd = pd, pi0 = 0.7, threshold = 0.99, alpha = 0.10)
  expect_equal(attr(ov, "threshold"), 0.95)

  expect_error(pd.adjust(pd = pd, pi0 = 0.7, alpha = 0),
               "`alpha`: must be a single number")
  expect_error(pd.adjust(pd = pd, pi0 = 0.7, alpha = 1),
               "`alpha`: must be a single number")
})

test_that("pd.adjust() returns a classed object and keeps $ access", {
  out <- suppressWarnings(pd.adjust(pd = c(0.9, 0.8, 0.7), p0 = 0.4))
  expect_s3_class(out, "pd_adjust")
  expect_s3_class(out, "data.frame")
  expect_true(is.numeric(out$pd))          # column access still works
  expect_true(all(c("pd", "pi0", "m", "direction") %in% names(out)))
  expect_false("p0" %in% names(out))        # global null not returned
})

test_that("pd.adjust() no longer accepts the removed joint argument", {
  expect_error(suppressWarnings(pd.adjust(pd = c(0.9, 0.8), p0 = 0.4, joint = TRUE)),
               "unused argument")
})

test_that("print.pd_adjust returns invisibly and prints a header", {
  out <- suppressWarnings(pd.adjust(pd = c(0.9, 0.8), p0 = 0.4))
  expect_output(print(out), "Prior-odds threshold for directional claims")
  expect_identical(suppressWarnings(print(out)), out)
})

test_that("pi0 below 0.5 is set to 0.5, and pi0 = 0.5 is accepted silently", {
  pd <- c(0.9, 0.8, 0.7, 0.6)
  out <- expect_warning(pd.adjust(pd = pd, pi0 = 0.4), "set to 0.5")
  expect_true(all(out$pi0 == 0.5))
  expect_equal(attr(out, "threshold.adj"), attr(out, "threshold"))
  expect_no_warning(out5 <- pd.adjust(pd = pd, pi0 = 0.5))
  expect_equal(attr(out5, "threshold.adj"), attr(out5, "threshold"))
})

test_that("print.pd_adjust reports c, c* and the constant columns in the header", {
  set.seed(1)
  dr <- matrix(rnorm(4000 * 3, c(0, 1, 3), 1), ncol = 3, byrow = TRUE)
  colnames(dr) <- c("a", "b", "c")
  out <- pd.adjust(draws = dr, pi0 = 0.7)
  txt <- capture.output(print(out))
  expect_true(any(grepl("m = 3, pi0 = 0.7", txt, fixed = TRUE)))
  expect_true(any(grepl(sprintf("c = 0.975, c* = %.4g", attr(out, "threshold.adj")),
                        txt, fixed = TRUE)))
  expect_true(any(grepl("null value = 0, direction = two.sided", txt, fixed = TRUE)))
  expect_false(any(grepl("Pr(H0)|alpha", txt)))
  # constant columns are not repeated in the table
  colhdr <- txt[grepl("median.est", txt)]
  expect_length(colhdr, 1)
  expect_false(grepl("direction|null.value", colhdr))
})
