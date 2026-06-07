# Prior-odds adjustment for Probability of Direction *pd*

The function accepts either a vector of pre-computed *pd* values or a
matrix of posterior draws, from which *pd* values are computed
internally. Both direction-agnostic and directional tests are supported:
the `direction` argument controls which formulation is applied per
hypothesis. The adjustment is governed by the per-test prior probability
that an individual effect is null, \\\pi_0\\ (argument `pi0`).
Equivalently, one may supply the global null \\\Pi_0\\ (argument `p0`),
the prior that all \\m\\ hypotheses are null, from which \\\pi_0 =
\Pi_0^{1/m}\\ is recovered under independence.

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
  joint = FALSE
)
```

## Arguments

- pd:

  Numeric vector of *pd* values. For direction-agnostic tests, values
  must be in \\\[0.5, 1\]\\. For directional tests, values are raw
  one-sided probabilities in \\\[0, 1\]\\. Ignored if `draws` is
  supplied.

- draws:

  Optional matrix or data frame of posterior draws (columns =
  parameters). If provided, *pd* values are computed automatically from
  the draws according to `direction` and `null.value`.

- p0:

  Numeric scalar in \\(0, 1)\\. The prior probability that **all**
  hypotheses are null simultaneously (the global null, \\\Pi_0\\).
  Supply **either** `pi0` **or** `p0`, not both. When `p0` is given, the
  per-test prior is derived under independence as \\\pi_0 =
  \Pi_0^{1/m}\\.

- pi0:

  Numeric scalar in \\(0, 1)\\. The per-test prior probability that an
  individual effect is null, applied to every test. Supply **either**
  `pi0` **or** `p0`, not both. When `pi0` is given, the implied global
  null is \\\Pi_0 = \pi_0^{m}\\ (valid under independence) and is
  returned in `p0`.

- null.value:

  Numeric scalar or vector. The null (reference) value against which the
  posterior is evaluated, specified on the scale of the posterior. A
  single scalar applies the same null to all parameters; a vector of
  length `ncol(draws)` assigns a distinct null to each parameter.
  Ignored when `pd` is supplied directly. Defaults to `0`.

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
  columns of `draws` (or the elements of `pd`); use this for
  programmatic work that pairs `pd.adj` with external vectors. Ordering
  is enabled automatically when `joint = TRUE`.

- joint:

  Logical. If `TRUE`, the output *additionally* includes the cumulative
  joint probability `joint_cum` (see Details). Because that cumulative
  joint is only interpretable when accumulated from the strongest claim
  down, `joint = TRUE` orders the rows by decreasing `pd` (i.e. it
  implies `order = TRUE`). Defaults to `FALSE`, in which case only the
  per-test correction is returned.

## Value

An object of class `pd_adjust` (a `data.frame` with a dedicated
[`print.pd_adjust`](https://mar-cald.github.io/multibayes/reference/print.pd_adjust.md)
method), with one row per hypothesis. Columns: `pd` (values used in the
adjustment), `pd.adj` (adjusted values), `pi0` (per-test null prior
\\\pi_0\\), and `m` (number of tests). The global null \\\Pi_0 =
\pi_0^m\\ is **not** returned: it is exact only under independence, and
reporting it alongside per-test results can mislead when the tests are
dependent; the per-test prior \\\pi_0\\ is the robust quantity (the
expected proportion of nulls under exchangeability). When `joint = TRUE`
a `joint_cum` column is added; when `order = TRUE` the rows are sorted
by decreasing `pd`. For direction-agnostic tests, both `pd` and `pd.adj`
are bounded in \\\[0.5, 1\]\\; for directional tests, both are on \\\[0,
1\]\\, with values below \\0.5\\ indicating that the data (and the
adjustment) favoured the opposite direction. When `draws` are supplied,
the output additionally includes `median.est` (posterior median per
parameter), `null.value`, and `direction`. The `print` method summarises
the constant quantities (`pi0`, `m`) in a header and displays the
per-test table; the columns themselves remain available for programmatic
access.

## Details

The adjustment follows from Bayes' theorem. Given a per-hypothesis prior
\\\pi_0 = \Pi_0^{1/m}\\ and its complement \\\pi_1 = 1 - \pi_0\\, the
adjusted *pd* is: \$\$ pd\_{adj} = \frac{pd \pi_1}{pd \pi_1 +
(1-pd)\pi_0} \$\$

Because the prior is conservative (\\\pi_0 \> \pi_1\\), the adjustment
always shrinks *pd* toward its lower bound.

**Direction-agnostic tests** (`"two.sided"`): *pd* is defined as
\\\max\\\big(\Pr(\hat\theta \> \theta\_\text{null}),\\ \Pr(\hat\theta \<
\theta\_\text{null})\big)\\ and is bounded in \\\[0.5, 1\]\\ by
construction. \\pd\_{adj}\\ is also floored at \\0.5\\, so the
adjustment produces shrinkage toward \\0.5\\.

**Directional tests** (`"greater"` or `"less"`): *pd* is the raw
posterior probability mass on the predicted side, \\\Pr(\hat\theta \>
\theta\_\text{null})\\ or \\\Pr(\hat\theta \< \theta\_\text{null})\\,
and is defined on \\\[0, 1\]\\. Values of *pd* below \\0.5\\ indicate
that the posterior is concentrated opposite to the predicted direction;
the adjustment will further shrink such values toward \\0\\, reflecting
the combined weight of the data and the conservative prior against the
hypothesis.

Mixed use of directional and direction-agnostic tests within the same
call is supported: each element of `direction` is handled independently,
and the same prior-odds adjustment is applied uniformly across all
hypotheses regardless of their directionality.

**Per-test correction and joint statement.** By default the function
returns the per-test correction only (the marginal `pd` and its adjusted
value `pd.adj`). Setting `joint = TRUE` *additionally* returns
`joint_cum` and orders the rows by decreasing `pd`, so that `joint_cum`
is the cumulative joint probability that the \\k\\ strongest claims all
hold in their claimed direction. When `draws` are supplied this joint is
computed directly from the draws by intersecting the per-draw
directional events, so it reflects the posterior dependence among
parameters; when only `pd` is supplied it is the running product of the
*pd* values, valid only under independence. The claimed direction for a
`"two.sided"` test is its dominant posterior side, and the last entry of
`joint_cum` is the joint probability that all \\m\\ claims hold
simultaneously.

## Examples

``` r
# From a vector of pd values
pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763,
H5 = 0.891, H6 = 0.987)
pd.adjust(pd = pd_values, p0 = 0.4)
#> Warning: some pd.adj have been floored to 0.5.
#> Prior-odds adjusted probability of direction
#> Tests: 6  | pi0 = 0.8584
#> 
#>       pd pd.adj direction
#> H1 0.999 0.9940 two.sided
#> H2 0.946 0.7430 two.sided
#> H3 0.813 0.5000 two.sided
#> H4 0.763 0.5000 two.sided
#> H5 0.891 0.5742 two.sided
#> H6 0.987 0.9261 two.sided

# Equivalent call specifying the per-test prior pi0 directly
pd.adjust(pd = pd_values, pi0 = 0.4^(1/6))
#> Warning: some pd.adj have been floored to 0.5.
#> Prior-odds adjusted probability of direction
#> Tests: 6  | pi0 = 0.8584
#> 
#>       pd pd.adj direction
#> H1 0.999 0.9940 two.sided
#> H2 0.946 0.7430 two.sided
#> H3 0.813 0.5000 two.sided
#> H4 0.763 0.5000 two.sided
#> H5 0.891 0.5742 two.sided
#> H6 0.987 0.9261 two.sided

# Simulate draws
Sigma <- matrix(0, nrow = 6, ncol = 6); diag(Sigma) <- 1
mu    <- c(1, -0.1, 0.8, 0, 2, 3)
draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")

# Per-test correction only (default)
pd.adjust(draws = draws, p0 = 0.4, null.value = 0)
#> Warning: some pd.adj have been floored to 0.5.
#> Prior-odds adjusted probability of direction
#> Tests: 6  | pi0 = 0.8584
#> 
#>    median.est null.value     pd pd.adj direction
#> H1     0.9602          0 0.8410 0.5000 two.sided
#> H2    -0.0793          0 0.5330 0.5000 two.sided
#> H3     0.7481          0 0.7708 0.5000 two.sided
#> H4    -0.0349          0 0.5118 0.5000 two.sided
#> H5     2.0256          0 0.9755 0.8679 two.sided
#> H6     2.9828          0 0.9998 0.9985 two.sided

# Also return the cumulative joint statement, ordered by evidence
pd.adjust(draws = draws, p0 = 0.4, null.value = 0, order = TRUE, joint = TRUE)
#> Warning: some pd.adj have been floored to 0.5.
#> Prior-odds adjusted probability of direction
#> Tests: 6  | pi0 = 0.8584
#> 
#>    median.est null.value     pd pd.adj direction joint_cum
#> H6     2.9828          0 0.9998 0.9985 two.sided    0.9998
#> H5     2.0256          0 0.9755 0.8679 two.sided    0.9752
#> H1     0.9602          0 0.8410 0.5000 two.sided    0.8202
#> H3     0.7481          0 0.7708 0.5000 two.sided    0.6328
#> H2    -0.0793          0 0.5330 0.5000 two.sided    0.3405
#> H4    -0.0349          0 0.5118 0.5000 two.sided    0.1663
#> 
#> joint_cum: cumulative joint Pr(first k claims hold), in the order shown.

# Mix of directional and agnostic tests with parameter-specific nulls
pd.adjust(draws = draws, p0 = 0.2, null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
          direction = c("greater", "two.sided", "greater",
          "two.sided", "greater", "greater"))
#> Warning: some pd.adj have been floored to 0.5.
#> Prior-odds adjusted probability of direction
#> Tests: 6  | pi0 = 0.7647
#> 
#>    median.est null.value     pd pd.adj direction
#> H1     0.9602        0.2 0.7832 0.5265   greater
#> H2    -0.0793        0.0 0.5330 0.5000 two.sided
#> H3     0.7481        0.2 0.7035 0.4220   greater
#> H4    -0.0349        0.0 0.5118 0.5000 two.sided
#> H5     2.0256        0.5 0.9368 0.8200   greater
#> H6     2.9828        0.5 0.9958 0.9863   greater
```
