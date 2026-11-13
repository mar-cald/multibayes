# Prior-odds adjustment for the probability of direction

Decides which directional claims to make in a family of tests, from the
probability of direction (*pd*) of each parameter. Two separate choices
enter the decision. The threshold \\c\\ (`threshold`) is the posterior
probability required to make a claim. The bound \\\bar{\pi}\\ (`pi.max`)
is the largest prior probability allowed for either direction of a
parameter, and sets how skeptical the prior may be about a claim. For
each test the function returns `pd` and `pd.adj`, the posterior
probability of the favored direction when that direction has prior
probability \\1 - \bar{\pi}\\. A direction is claimed when `pd.adj` is
at least \\c\\.

## Usage

``` r
pd.adjust(
  pd = NULL,
  draws = NULL,
  p0 = NULL,
  pi.max = NULL,
  null.value = 0,
  direction = NULL,
  order = FALSE,
  threshold = 0.975,
  alpha = NULL
)
```

## Arguments

- pd:

  Numeric vector of *pd* values in \\\[0.5, 1\]\\. Every test is then
  treated as `"two.sided"` and `direction` is ignored. Ignored if
  `draws` is supplied.

- draws:

  Matrix or data frame of posterior draws, one column per parameter.
  *pd* is computed from the draws according to `direction` and
  `null.value`.

- p0:

  Number in \\(0, 1)\\. The largest prior probability that no claim in
  the family is correct, \\\Pr(H_0)\\, from which \\\bar{\pi} =
  \Pr(H_0)^{1/m}\\. Supply either `p0` or `pi.max`, not both.

- pi.max:

  Number in \\(0, 1)\\. The largest prior probability allowed for either
  direction of a parameter, chosen before the data. See Details for how
  to choose it. Values below 0.5 are set to 0.5, with a warning, so that
  no adjustment is made. Supply either `pi.max` or `p0`, not both.

- null.value:

  Number or vector, the null value of each parameter on the scale of the
  draws. A single number applies to all parameters. Ignored when `pd` is
  supplied. Defaults to `0`.

- direction:

  Character vector with elements `"two.sided"`, `"greater"`, or
  `"less"`, one per parameter, or a single value for all. `NULL` (the
  default) means `"two.sided"` for all.

- order:

  Logical. If `TRUE`, rows are sorted by decreasing `pd`. If `FALSE`
  (the default), rows keep the order of the input.

- threshold:

  Number in \\(0.5, 1)\\. The threshold \\c\\, the posterior probability
  required to make a claim. Defaults to `0.975`.

- alpha:

  Optional number in \\(0, 1)\\. If supplied, the threshold is set to
  \\c = 1 - \alpha/2\\ and `threshold` is ignored.

## Value

An object of class `pd_adjust`, a data frame with one row per test and
the columns `pd`, `pd.adj`, `pi.max`, `m` (the number of tests), and
`direction`. With `draws`, it also has `median.est` (the posterior
median) and `null.value`.

The thresholds are stored as attributes, `threshold` (\\c\\) and
`threshold.adj` (\\c^\*\\), together with `alpha.nominal` (\\2(1-c)\\)
and `alpha.eff` (\\2(1-c^\*)\\). The last two are per-test error rates
only if *pd* is uniform on \\\[0.5, 1\]\\ for a null effect, which holds
approximately under diffuse priors.

See
[`print.pd_adjust`](https://mar-cald.github.io/multibayes/reference/print.pd_adjust.md)
and
[`plot.pd_adjust`](https://mar-cald.github.io/multibayes/reference/plot.pd_adjust.md).

## Details

**What is computed.** *pd* comes from the fitted model, whose prior (the
estimation prior) must be continuous and symmetric about `null.value`.
Flat priors and symmetric priors centred on `null.value`, such as normal
or Student-t priors, qualify. Each direction then has prior probability
0.5, and \\pd/(1 - pd)\\ is the Bayes factor of the favored direction
against the other.

With several tests, the prior odds of each claim are lowered (Jeffreys,
1938; Westfall et al., 1997). The direction of a claim is not known
before the data, so the function considers every prior obtained by
reweighting the two sides of the estimation prior so that each direction
has prior probability between \\1 - \bar{\pi}\\ and \\\bar{\pi}\\. A
claim must hold under all of these priors, and it is enough to check the
one most skeptical of the favored direction, under which \$\$pd\_{adj} =
\frac{pd\\(1-\bar{\pi})}{pd\\(1-\bar{\pi}) + (1-pd)\\\bar{\pi}}.\$\$ The
model is not refitted, and the posterior used for estimation is not
changed.

**The threshold on pd.** The threshold \\c\\ is not changed by the
adjustment. Because \\pd\_{adj}\\ increases with *pd*, \\pd\_{adj} \ge
c\\ holds exactly when *pd* reaches \$\$c^\* =
\frac{c\\\bar{\pi}}{(1-c)(1-\bar{\pi}) + c\\\bar{\pi}},\$\$ the value of
*pd* required by \\c\\ and \\\bar{\pi}\\ together. The adjustment is
thus equivalent to raising the threshold on *pd* from \\c\\ to \\c^\*\\,
reported as `threshold.adj`.

**Choosing pi.max.** `pi.max` says how skeptical you are of a
directional claim before the data. With `pi.max = 0.5` no adjustment is
made. With `pi.max = 0.8`, a claim needs a Bayes factor four times
larger. The same `pi.max` is used for every test, so the parameters
should be exchangeable in direction, that is, no parameter is judged
more likely than another to be positive or negative. They can be on
different scales, because the sign of a parameter does not depend on its
unit. The proportion of effects in the family that could be null is a
reference value. If up to 70\\ direction, so each direction should be
allowed a prior probability as low as 0.3, which gives `pi.max = 0.7`. A
larger value gives fewer false claims at a cost in power. The priors
considered are continuous, so a null effect has probability zero under
all of them, and `pd.adj` is the posterior probability of a sign. To
test whether an effect is exactly null, use a prior that puts
probability on the null value, such as a spike-and-slab prior. A
narrower prior in the fitted model is no substitute, because a symmetric
prior gives each direction prior probability 0.5 however narrow it is.
Alternatively, `p0` gives the largest prior probability that no claim in
the family is correct. For parameters independent a priori, \\\bar{\pi}
= \Pr(H_0)^{1/m}\\, which grows with the number of tests \\m\\.

**Direction.** With `"two.sided"`, *pd* is the larger of \\\Pr(\theta \>
\theta\_\text{null} \mid y)\\ and \\\Pr(\theta \< \theta\_\text{null}
\mid y)\\, between 0.5 and 1. With `"greater"` or `"less"`, *pd* is the
posterior probability of the predicted side, between 0 and 1, and
`pd.adj` is computed under the prior most skeptical of that side. The
mode can differ from one test to another.

## References

Berger, J. O. (1994). An overview of robust Bayesian analysis. *Test*,
3(1), 5–124. https://doi.org/10.1007/BF02562676

Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A Bayesian
perspective on the Bonferroni adjustment. *Biometrika*, 84(2), 419–427.
https://doi.org/10.1093/biomet/84.2.419

Jeffreys, H. (1938). Significance tests when several degrees of freedom
arise simultaneously. *Proceedings of the Royal Society of London.
Series A. Mathematical and Physical Sciences*, 165(921), 161–198.
https://doi.org/10.1098/rspa.1938.0052

Storey, J. D. (2003). The positive false discovery rate: A Bayesian
interpretation and the q-value. *The Annals of Statistics*, 31(6),
2013–2035. https://doi.org/10.1214/aos/1074290335

Scott, J. G., & Berger, J. O. (2006). An exploration of aspects of
Bayesian multiple testing. *Journal of Statistical Planning and
Inference*, 136(7), 2144–2162.
https://doi.org/10.1016/j.jspi.2005.08.031

Scott, J. G., & Berger, J. O. (2010). Bayes and empirical-Bayes
multiplicity adjustment in the variable- selection problem. *The Annals
of Statistics*, 38(5), 2587–2619. https://doi.org/10.1214/10-AOS792

## Examples

``` r
# From a vector of pd values
pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763,
H5 = 0.891, H6 = 0.987)
pd.adjust(pd = pd_values, p0 = 0.4)
#> Prior-odds adjustment for directional claims
#>   m = 6, pi.max = 0.8584
#>   c = 0.975, c* = 0.9958
#>   direction = two.sided
#> 
#>       pd pd.adj
#> H1 0.999 0.9940
#> H2 0.946 0.7430
#> H3 0.813 0.4177
#> H4 0.763 0.3469
#> H5 0.891 0.5742
#> H6 0.987 0.9261
#> 
#> pd.adj: pd when the favored direction has prior probability 1 - pi.max.

# The same, giving pi.max directly
pd.adjust(pd = pd_values, pi.max = 0.4^(1/6))
#> Prior-odds adjustment for directional claims
#>   m = 6, pi.max = 0.8584
#>   c = 0.975, c* = 0.9958
#>   direction = two.sided
#> 
#>       pd pd.adj
#> H1 0.999 0.9940
#> H2 0.946 0.7430
#> H3 0.813 0.4177
#> H4 0.763 0.3469
#> H5 0.891 0.5742
#> H6 0.987 0.9261
#> 
#> pd.adj: pd when the favored direction has prior probability 1 - pi.max.

# Simulate draws
Sigma <- matrix(0, nrow = 6, ncol = 6); diag(Sigma) <- 1
mu    <- c(1, -0.1, 0.8, 0, 2, 3)
draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")

# From posterior draws
pd.adjust(draws = draws, pi.max = 0.7, null.value = 0)
#> Prior-odds adjustment for directional claims
#>   m = 6, pi.max = 0.7
#>   c = 0.975, c* = 0.9891
#>   null value = 0, direction = two.sided
#> 
#>    median.est     pd pd.adj
#> H1     0.9860 0.8358 0.6856
#> H2    -0.1175 0.5350 0.3302
#> H3     0.7761 0.7822 0.6062
#> H4    -0.0089 0.5030 0.3025
#> H5     1.9500 0.9762 0.9463
#> H6     2.9866 0.9995 0.9988
#> 
#> pd.adj: pd when the favored direction has prior probability 1 - pi.max.

# Rows sorted by decreasing pd
pd.adjust(draws = draws, pi.max = 0.7, null.value = 0, order = TRUE)
#> Prior-odds adjustment for directional claims
#>   m = 6, pi.max = 0.7
#>   c = 0.975, c* = 0.9891
#>   null value = 0, direction = two.sided
#> 
#>    median.est     pd pd.adj
#> H6     2.9866 0.9995 0.9988
#> H5     1.9500 0.9762 0.9463
#> H1     0.9860 0.8358 0.6856
#> H3     0.7761 0.7822 0.6062
#> H2    -0.1175 0.5350 0.3302
#> H4    -0.0089 0.5030 0.3025
#> 
#> pd.adj: pd when the favored direction has prior probability 1 - pi.max.

# Different directions and null values across parameters
pd.adjust(draws = draws, pi.max = 0.6, null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
          direction = c("greater", "two.sided", "greater",
          "two.sided", "greater", "greater"))
#> Prior-odds adjustment for directional claims
#>   m = 6, pi.max = 0.6
#>   c = 0.975, c* = 0.9832
#> 
#>    median.est null.value     pd pd.adj direction
#> H1     0.9860        0.2 0.7820 0.7051   greater
#> H2    -0.1175        0.0 0.5350 0.4341 two.sided
#> H3     0.7761        0.2 0.7230 0.6350   greater
#> H4    -0.0089        0.0 0.5030 0.4029 two.sided
#> H5     1.9500        0.5 0.9285 0.8965   greater
#> H6     2.9866        0.5 0.9948 0.9921   greater
#> 
#> pd.adj: pd when the favored direction has prior probability 1 - pi.max.
```
