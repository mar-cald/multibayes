# Simultaneous credible intervals and joint directional statements

Two ways of making directional claims about several parameters from
their joint posterior. By default, `joint()` returns simultaneous
credible intervals. With `interval = FALSE` it returns, for each \\k\\,
the joint posterior probability that the \\k\\ strongest directional
claims all hold.

## Usage

``` r
joint(
  draws,
  interval = TRUE,
  prob = 0.95,
  null.value = 0,
  direction = NULL,
  est.FUN = median
)
```

## Arguments

- draws:

  Matrix or data frame of posterior draws, one column per parameter,
  with at least 1000 rows and 2 columns.

- interval:

  Logical. If `TRUE` (the default), simultaneous credible intervals are
  returned. If `FALSE`, the joint directional statement.

- prob:

  Number in \\(0, 1)\\. Joint posterior probability of the intervals.
  Defaults to `0.95`. Used only when `interval = TRUE`.

- null.value:

  Number or vector, the null value of each parameter. A single number
  applies to all parameters. Defaults to `0`. Used only when
  `interval = FALSE`.

- direction:

  Character vector with elements `"two.sided"`, `"greater"`, or
  `"less"`, one per parameter, or a single value for all. `NULL` (the
  default) means `"two.sided"` for all, that is, the side favored by the
  posterior. Used only when `interval = FALSE`.

- est.FUN:

  Function that computes a point estimate from the draws of one
  parameter. Defaults to
  [`median`](https://rdrr.io/r/stats/median.html). Used only when
  `interval = TRUE`.

## Value

A data frame with one row per parameter and a print method. With
`interval = TRUE`, an object of class `joint_interval` with the columns
`lower`, `est`, `upper`, and `prob`. With `interval = FALSE`, an object
of class `joint_directional`, sorted by decreasing *pd*, with the
columns `median.est`, `null.value`, `pd`, `direction`, and `joint_cum`.

## Details

**Simultaneous credible intervals** (`interval = TRUE`). Equal-tailed
intervals such that all parameters lie inside their interval at the same
time with posterior probability `prob`. For every draw, the function
takes the smallest tail probability across parameters, and the intervals
are set at the \\(1 - \texttt{prob})\\-quantile of these values. A
direction can be claimed for a parameter when its interval excludes the
null value.

The ends of the intervals lie further in the tails as the number of
parameters grows, so more draws are needed. Each end rests on about \\S
(1 - \texttt{prob}^{1/m}) / 2\\ draws, where \\S\\ is the number of
draws and \\m\\ the number of parameters. With fewer than 10, the
intervals are too narrow on average and `joint()` gives a warning. With
`prob = 0.95`, this means at least about 1,600 draws for 4 parameters,
7,800 for 20, and 15,600 for 40. `sim.cred.band` in the credsubs package
(Schnell et al., 2020) is a closely related implementation.

**Joint directional statement** (`interval = FALSE`). The direction of
each parameter is the side of `null.value` favored by its posterior, or
the side given in `direction`. Parameters are sorted by decreasing *pd*,
and `joint_cum` gives, for each \\k\\, the proportion of draws in which
the first \\k\\ directions all hold. It therefore accounts for the
posterior dependence among parameters. This is a probability about the
signs of the parameters under the prior of the fitted model, with no
adjustment of the prior odds. In families with many null effects it
leads to many false claims.

## References

Box, G. E. P., & Tiao, G. C. (1992). Bayesian inference in statistical
analysis (1st ed.). *Wiley*.

Besag, J., Green, P., Higdon, D., & Mengersen, K. (1995). Bayesian
computation and stochastic systems. *Statistical Science*, 10(1), 3–41.

Held, L. (2004). Simultaneous posterior probability statements from
Monte Carlo output. *Journal of Computational and Graphical Statistics*,
13(1), 20–35.

Gelman, A., Tuerlinckx, F. (2000). Type S error rates for classical and
Bayesian single and multiple comparison procedures. *Computational
Statistics* 15, 373–390 (2000). https://doi.org/10.1007/s001800000040

Schnell, P. et al. (2020). credsubs: Credible Subsets. R package.
<https://CRAN.R-project.org/package=credsubs>.

## See also

[`pd.adjust`](https://mar-cald.github.io/multibayes/reference/pd.adjust.md)
for the per-test prior-odds threshold.

## Examples

``` r
mu <- c(1.5, 0, -2)
Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.3, 0.5, 0.3, 1), 3, 3)
draws <- MASS::mvrnorm(2000, mu, Sigma)
colnames(draws) <- c("theta1", "theta2", "theta3")

joint(draws)                     # simultaneous credible intervals (default)
#> Simultaneous credible intervals (joint probability 0.95)
#> 
#>          lower     est  upper
#> theta1 -0.7416  1.4721 3.8524
#> theta2 -2.3636 -0.0241 2.1471
#> theta3 -4.3930 -1.9897 0.3948
joint(draws, interval = FALSE)   # joint directional statement
#> Joint directional statement
#>   null value = 0, direction = two.sided
#>   joint_cum: Pr(the top k claims all hold | y)
#> 
#>        median.est     pd joint_cum
#> theta3    -1.9897 0.9750    0.9750
#> theta1     1.4721 0.9360    0.9115
#> theta2    -0.0241 0.5095    0.4415
```
