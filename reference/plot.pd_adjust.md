# Plot method for `pd_adjust` objects

For each hypothesis, draws the range of the posterior probability of the
favored (or predicted) direction across the class of priors behind
[`pd.adjust`](https://mar-cald.github.io/multibayes/reference/pd.adjust.md).
The left end of each bar is `pd.adj`, the posterior probability under
the prior most skeptical of that direction; the right end is the
posterior probability under the prior most favorable to it,
\\pd\\\bar{\pi} / (pd\\\bar{\pi} + (1 - pd)(1 - \bar{\pi}))\\. The point
marks `pd`, the value under the estimation prior, and the two ends of
each range are also printed on the right. A direction is claimed when
the whole bar lies at or above the threshold \\c\\ (dashed line), that
is, when `pd.adj` \\\ge c\\, which is the same as \\pd \ge c^\*\\. A bar
that crosses 0.5 means that the priors in the class disagree on which
direction is more probable.

## Usage

``` r
# S3 method for class 'pd_adjust'
plot(x, labels = rownames(x), xlim = NULL, digits = 3, ...)
```

## Arguments

- x:

  A `pd_adjust` object returned by
  [`pd.adjust`](https://mar-cald.github.io/multibayes/reference/pd.adjust.md).

- labels:

  Character vector of labels for the hypotheses, in the order of the
  rows of `x`. Defaults to the row names.

- xlim:

  Range of the x-axis. Defaults to the range of the bars, extended to
  include \\c\\.

- digits:

  Number of digits for the printed ranges. Defaults to 3.

- ...:

  Further graphical arguments passed to `plot.default`.

## Value

Invisibly, a data frame with one row per hypothesis and columns `lower`
(equal to `pd.adj`), `pd`, `upper`, and `claim`.

## Examples

``` r
set.seed(1)
draws <- MASS::mvrnorm(4000, mu = c(0.4, 0.15, 0.05),
                       Sigma = diag(c(0.12, 0.12, 0.12)^2))
colnames(draws) <- c("b1", "b2", "b3")
res <- pd.adjust(draws = draws, pi.max = 0.8)
plot(res)

```
