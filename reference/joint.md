# Simultaneous credible intervals from joint posterior draws Computes simultaneous equitailed credible intervals for a set of parameters, guaranteeing that **all** parameters stay within their bounds simultaneously with a user specified probability.

Simultaneous credible intervals from joint posterior draws Computes
simultaneous equitailed credible intervals for a set of parameters,
guaranteeing that **all** parameters stay within their bounds
simultaneously with a user specified probability.

## Usage

``` r
joint(draws, prob = 0.95, est.FUN = median, null.value = NULL)
```

## Arguments

- draws:

  Numeric matrix or data frame of posterior draws (rows = draws, columns
  = parameters).

- prob:

  Numeric scalar in \\(0, 1)\\. Joint coverage probability (default
  `0.95`).

- est.FUN:

  Function for point estimates (default `median`).

- null.value:

  Optional numeric vector of null values to test against the
  simultaneous credible intervals, one per parameter. If `NULL`
  (default), no test is performed. When provided, a logical column
  `test` is appended: `TRUE` means the null value is **excluded** from
  the interval.

## Value

A dataframe with:

- lower:

  Lower bounds, one per parameter.

- est:

  Point estimates via `est.FUN`.

- upper:

  Upper bounds, one per parameter.

- prob:

  Requested joint coverage probability.

- cq:

  Critical value.

- sig:

  Logical; `TRUE` if `null.value` is excluded from the interval. Only
  present when `null.value` is supplied.

## Details

Simultaneous coverage is calibrated by examining how extreme each draw
is across all parameters jointly. For each draw, the minimum tail
probability across all parameters is computed. The simultaneous
threshold is the \\\alpha\\-quantile of these minima: only \\\alpha\\ of
draws are more extreme than this value in at least one parameter
simultaneously.

A closely related implementation is `sim.cred.band` in the credsubs
package (Schnell et al., 2020).

## Note

If you use this function in published research, please cite:

- The package: <https://github.com/mar-cald/multibayes>.

## Examples

``` r
mu <- c(4, 0, -2)
Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.3, 0.5, 0.3, 1), 3, 3)
draws <- MASS::mvrnorm(2000, mu, Sigma)
colnames(draws) <- c("theta1", "theta2", "theta3")

joint(draws, prob = 0.95)
#>            lower          est     upper prob     cq
#> theta1  1.704391  4.020594092 6.2265003 0.95 0.0095
#> theta2 -2.413244  0.002074391 2.2378859 0.95 0.0095
#> theta3 -4.213334 -1.959627191 0.3455705 0.95 0.0095
joint(draws, prob = 0.95, null.value = c(0, 0, 0))
#>            lower          est     upper prob     cq  test
#> theta1  1.704391  4.020594092 6.2265003 0.95 0.0095  TRUE
#> theta2 -2.413244  0.002074391 2.2378859 0.95 0.0095 FALSE
#> theta3 -4.213334 -1.959627191 0.3455705 0.95 0.0095 FALSE
```
