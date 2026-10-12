# Prior-odds threshold for the Probability of Direction *pd*

Converts an upper bound, set by the researcher, on the prior probability
of the sign opposite to the claimed one into the threshold that *pd*
must exceed before a direction is claimed. The function accepts either a
vector of pre-computed *pd* values or a matrix of posterior draws, from
which *pd* values are computed internally. Both direction-agnostic and
directional tests are supported: the `direction` argument controls how
*pd* is computed per hypothesis. The threshold is governed by \\\pi_0\\
(argument `pi0`), the upper bound on the prior probability that a
parameter lies on the side of `null.value` opposite to the claimed one
(or equals `null.value`). Equivalently, one may supply \\\Pr(H_0)\\
(argument `p0`), the largest prior probability that no claim in the
family is correct, from which \\\pi_0 = \Pr(H_0)^{1/m}\\ is recovered if
the parameters are independent a priori. For each hypothesis the
function also reports `pd.adj`, the posterior probability of the claimed
direction (favoured or predicted) under the prior in the class most
skeptical of it; under the same prior the other side has probability
`1 - pd.adj`, so the two sum to one. A direction is claimed when
`pd.adj` exceeds \\c\\, which is the same as \\pd \> c^\*\\.

## Usage

``` r
pd.adjust(
  pd = NULL,
  draws = NULL,
  p0 = NULL,
  pi0 = NULL,
  null.value = 0,
  direction = NULL,
  order = FALSE,
  threshold = 0.975,
  alpha = NULL
)
```

## Arguments

- pd:

  Numeric vector of direction-agnostic *pd* values in \\\[0.5, 1\]\\.
  With `pd` input every test is treated as direction-agnostic and
  `direction` is ignored; for directional tests, supply `draws`. Ignored
  if `draws` is supplied.

- draws:

  Optional matrix or data frame of posterior draws (columns =
  parameters). If provided, *pd* values are computed automatically from
  the draws according to `direction` and `null.value`.

- p0:

  Numeric scalar in \\(0, 1)\\. The largest prior probability that
  **no** directional claim in the family is correct (\\\Pr(H_0)\\).
  Supply **either** `pi0` **or** `p0`, not both. When `p0` is given,
  \\\pi_0\\ is derived, for parameters independent a priori, as \\\pi_0
  = \Pr(H_0)^{1/m}\\.

- pi0:

  Numeric scalar in \\(0, 1)\\. An upper bound, set by the researcher,
  on the prior probability that a parameter lies on the side of
  `null.value` opposite to the claimed one (or equals `null.value`).
  Because the same `pi0` is applied to every test, it is also an upper
  bound on the expected proportion of parameters in the family that are
  not on the claimed side, null effects included, so a judgment about
  how many null effects the family could contain can guide its choice.
  Values below 0.5 are set to 0.5, with a warning, so that no penalty is
  applied. Supply **either** `pi0` **or** `p0`, not both.

- null.value:

  Numeric scalar or vector. The null (reference) value against which the
  posterior is evaluated, specified on the scale of the posterior. A
  single scalar applies the same null to all parameters; a vector of
  length `ncol(draws)` assigns a distinct null to each parameter.
  Ignored when `pd` is supplied directly. Defaults to `0`. The class of
  priors behind the threshold is defined for a within-model prior
  symmetric about the null value (see Details), which holds for flat
  priors and for priors centred on it.

- direction:

  Character vector of `"greater"`, `"less"`, or `"two.sided"` (or
  `NULL`). Specifies the testing mode for each hypothesis: `"greater"`
  for a positive directional test (\\\Pr(\theta \>
  \theta\_\text{null})\\), `"less"` for a negative directional test
  (\\\Pr(\theta \< \theta\_\text{null})\\), and `"two.sided"` for
  direction-agnostic testing (maximum over both sides). A scalar is
  recycled to match the number of parameters; a mixed vector allows
  different modes across hypotheses. Defaults to `NULL`
  (direction-agnostic for all parameters).

- order:

  Logical. If `TRUE`, the returned rows are ordered by decreasing `pd`
  (strongest directional evidence first). If `FALSE` (the default), rows
  are returned in input order, so the output stays aligned with the
  columns of `draws` (or the elements of `pd`).

- threshold:

  Numeric scalar in \\(0.5, 1)\\. The threshold \\c\\ on *pd* without
  penalty (default `0.975`, a two-sided per-test \\\alpha = 0.05\\). The
  function reports the threshold implied by \\c\\ and `pi0`, \\c^\* =
  c\\\pi_0 / ((1-\pi_0)(1-c) + c\\\pi_0)\\.

- alpha:

  Optional numeric scalar in \\(0, 1)\\. A convenience parameterization
  of `threshold` on the error scale: when supplied, the two-sided cutoff
  is set to \\c = 1 - \alpha/2\\ (so `alpha = 0.05` gives
  `threshold = 0.975`), **overriding** any `threshold` value.

## Value

An object of class `pd_adjust` (a `data.frame` with a dedicated
[`print.pd_adjust`](https://mar-cald.github.io/multibayes/reference/print.pd_adjust.md)
method), with one row per hypothesis. Columns: `pd` (the probability of
direction, to be compared with `threshold.adj`), `pd.adj` (the posterior
probability of the claimed direction under the prior in the class most
skeptical of it, to be compared with `threshold`), `pi0`, `m`, and
`direction`. When `order = TRUE` the rows are sorted by decreasing `pd`.
When `draws` are supplied, the output additionally includes `median.est`
(posterior median per parameter) and `null.value`. The `print` method
shows \\m\\, \\\pi_0\\, the thresholds \\c\\ and \\c^\*\\, and the
columns that are constant across rows (`null.value`, `direction`) in a
header, followed by the per-test table; all columns remain available for
programmatic access. The object also carries, as attributes, the
threshold `threshold` (\\c\\), the adjusted threshold `threshold.adj`
(\\c^\*\\), the corresponding two-sided per-test Type I error rates
`alpha.nominal` (\\2(1-c)\\) and `alpha.eff` (\\2(1-c^\*)\\); retrieve
them with, e.g., `attr(x, "threshold.adj")`. The threshold
`threshold.adj` is an exact function of \\c\\ and \\\pi_0\\, but the
reported per-test *rate* `alpha.eff` uses the identity \\\alpha = 2(1 -
c)\\, which holds **only when the null *pd* is approximately uniform on
\\\[0.5, 1\]\\**.

## Details

Following Jeffreys (1938) and Westfall et al. (1997), the prior odds of
each directional claim in the family are penalized. With `pi0` the
penalty is the same for any number of tests; with `p0`, \\\pi_0 =
\Pr(H_0)^{1/m}\\ and the penalty grows with the number of tests. Because
the direction of a claim is chosen by the data, the penalty is applied
over a class of priors: the within-model prior is reweighted on the two
sides of `null.value`, so that the prior probability of each direction
lies between \\1-\pi_0\\ and \\\pi_0\\, and a claim is made only if its
posterior probability exceeds \\c\\ under every prior in the class. For
parameters with independent posteriors, this is equivalent to \$\$pd \>
c^\*, \qquad c^\* = \frac{c\\\pi_0}{(1-c)(1-\pi_0) + c\\\pi_0},\$\$ and
the function reports \\c^\*\\ (`threshold.adj`). The smallest posterior
probability is reached under the prior that gives the claimed direction
probability \\1 - \pi_0\\; under that prior the posterior probability of
the claimed direction is \$\$pd\_{adj} =
\frac{pd\\(1-\pi_0)}{pd\\(1-\pi_0) + (1-pd)\\\pi_0},\$\$ reported as
`pd.adj`, and \\pd\_{adj} \> c\\ if and only if \\pd \> c^\*\\. With
correlated posteriors, the bound over the class can be stricter than
this rule. The posterior used for estimation is not modified. The class
is defined for a within-model prior that is symmetric about
`null.value`, as are priors centred on it and, as a limit, flat priors.

**Direction-agnostic tests** (`"two.sided"`): *pd* is
\\\max\\\big(\Pr(\theta \> \theta\_\text{null} \mid y),\\ \Pr(\theta \<
\theta\_\text{null} \mid y)\big)\\ on \\\[0.5, 1\]\\, and a claim is
made in the favoured direction when \\pd \> c^\*\\.

**Directional tests** (`"greater"` or `"less"`): *pd* is the posterior
mass on the predicted side, on \\\[0, 1\]\\, and the predicted direction
is claimed when \\pd \> c^\*\\.

Mixed use of directional and direction-agnostic tests within the same
call is supported: each element of `direction` is handled independently,
and the same threshold \\c^\*\\ applies to all hypotheses regardless of
their directionality.

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
#> Prior-odds threshold for directional claims
#>   m = 6, pi0 = 0.8584
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
#> pd.adj: Pr(claimed direction | y) under the prior in the class most skeptical
#>         of it; 1 - pd.adj is the probability of the other side under that prior.
#> Claim a direction when pd.adj > 0.975, or equivalently pd > 0.9958.

# Equivalent call specifying pi0 directly
pd.adjust(pd = pd_values, pi0 = 0.4^(1/6))
#> Prior-odds threshold for directional claims
#>   m = 6, pi0 = 0.8584
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
#> pd.adj: Pr(claimed direction | y) under the prior in the class most skeptical
#>         of it; 1 - pd.adj is the probability of the other side under that prior.
#> Claim a direction when pd.adj > 0.975, or equivalently pd > 0.9958.

# Simulate draws
Sigma <- matrix(0, nrow = 6, ncol = 6); diag(Sigma) <- 1
mu    <- c(1, -0.1, 0.8, 0, 2, 3)
draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")

# From posterior draws
pd.adjust(draws = draws, pi0 = 0.7, null.value = 0)
#> Prior-odds threshold for directional claims
#>   m = 6, pi0 = 0.7
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
#> pd.adj: Pr(claimed direction | y) under the prior in the class most skeptical
#>         of it; 1 - pd.adj is the probability of the other side under that prior.
#> Claim a direction when pd.adj > 0.975, or equivalently pd > 0.9891.

# Order the rows by strength of evidence
pd.adjust(draws = draws, pi0 = 0.7, null.value = 0, order = TRUE)
#> Prior-odds threshold for directional claims
#>   m = 6, pi0 = 0.7
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
#> pd.adj: Pr(claimed direction | y) under the prior in the class most skeptical
#>         of it; 1 - pd.adj is the probability of the other side under that prior.
#> Claim a direction when pd.adj > 0.975, or equivalently pd > 0.9891.

# Mix of directional and agnostic tests with parameter-specific nulls
pd.adjust(draws = draws, pi0 = 0.6, null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
          direction = c("greater", "two.sided", "greater",
          "two.sided", "greater", "greater"))
#> Prior-odds threshold for directional claims
#>   m = 6, pi0 = 0.6
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
#> pd.adj: Pr(claimed direction | y) under the prior in the class most skeptical
#>         of it; 1 - pd.adj is the probability of the other side under that prior.
#> Claim a direction when pd.adj > 0.975, or equivalently pd > 0.9832.
```
