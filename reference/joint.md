# Summaries from joint posterior draws

Computes either **simultaneous credible intervals** (the default), or,
with `interval = FALSE`, the **joint directional statement**: for each
\\k\\, the joint posterior probability that the \\k\\ strongest
directional claims hold simultaneously.

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

  Numeric matrix or data frame of posterior draws (rows = draws, columns
  = parameters). At least 1000 rows and 2 columns.

- interval:

  Logical (default `TRUE`). If `TRUE`, return simultaneous credible
  intervals; if `FALSE`, return the joint directional statement.

- prob:

  Numeric scalar in \\(0, 1)\\. Joint posterior probability that all
  parameters lie within their intervals (default `0.95`). Used only when
  `interval = TRUE`.

- null.value:

  Numeric scalar or vector giving the null reference value(s) used to
  define the claimed direction (a scalar is recycled). Used only when
  `interval = FALSE`. Defaults to `0`.

- direction:

  Character vector of `"greater"`, `"less"`, or `"two.sided"` (or
  `NULL`), the testing mode per parameter. A scalar is recycled; `NULL`
  (default) uses `"two.sided"` (dominant side) for all. Used only when
  `interval = FALSE`.

- est.FUN:

  Function returning a point estimate from a numeric vector (default
  [`median`](https://rdrr.io/r/stats/median.html)). Used only when
  `interval = TRUE`.

## Value

A `data.frame`. For `interval = FALSE`, an object of class
`joint_directional` (a `data.frame` with a print method), with one row
per parameter ordered by decreasing *pd*, with columns `median.est`,
`null.value`, `pd`, `direction`, and `joint_cum` (the joint posterior
probability that the \\k\\ strongest claims all hold). For
`interval = TRUE`, an object of class `joint_interval` (a `data.frame`
with a print method), with one row per parameter and columns `lower`,
`est`, `upper`, and `prob`.

## Details

**Joint directional statement** (`interval = FALSE`). Each parameter's
claimed direction is its dominant posterior side relative to
`null.value` (or the side set by `direction`). Parameters are ordered by
decreasing *pd*, and the `joint_cum` column reports, for each \\k\\, the
joint posterior probability that the \\k\\ strongest claims hold
simultaneously. It is computed directly from the draws by intersecting
the per-draw directional events, so it reflects the posterior dependence
among parameters.

**Simultaneous credible intervals** (`interval = TRUE`, the default).
Equi-tailed intervals calibrated so that **all** parameters lie within
their bounds simultaneously with probability `prob`, via the
quantile-based simultaneous credible band. The intervals are set by
examining how extreme each draw is across all parameters jointly: for
every draw the minimum (nearer) tail probability across parameters is
taken, and the threshold is the \\(1 - \texttt{prob})\\-quantile of
these minima. A direction can be claimed for a parameter when its band
excludes the null value. A closely related implementation is
`sim.cred.band` in the credsubs package (Schnell et al., 2020).

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
