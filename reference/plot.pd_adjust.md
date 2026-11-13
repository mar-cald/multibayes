# Plot method for `pd_adjust` objects

Shows, for each test, how the posterior probability of the favored (or
predicted) direction changes with the prior. Each bar runs from
`pd.adj`, the value when that direction has prior probability \\1 -
\bar{\pi}\\, to the value when it has prior probability \\\bar{\pi}\\,
which is \\pd\\\bar{\pi} / (pd\\\bar{\pi} + (1 - pd)(1 - \bar{\pi}))\\.
The point marks `pd`, and the two ends of each bar are printed on the
right.

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

  Labels of the tests, in the order of the rows of `x`. Defaults to the
  row names.

- xlim:

  Range of the x-axis. By default it covers all bars and \\c\\.

- digits:

  Number of digits of the values printed on the right. Defaults to 3.

- ...:

  Further arguments passed to `plot.default`.

## Value

Invisibly, a data frame with one row per test and the columns `lower`
(equal to `pd.adj`), `pd`, `upper`, and `claim`.

## Details

A direction is claimed when the whole bar lies at or above the threshold
\\c\\ (dashed line), that is, when `pd.adj` is at least \\c\\. A bar
that crosses 0.5 means that the more probable direction depends on the
prior.

## Examples

``` r
set.seed(1)
draws <- MASS::mvrnorm(4000, mu = c(0.4, 0.15, 0.05),
                       Sigma = diag(c(0.12, 0.12, 0.12)^2))
colnames(draws) <- c("b1", "b2", "b3")
res <- pd.adjust(draws = draws, pi.max = 0.8)
plot(res)

```
