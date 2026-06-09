# Prior-odds adjustment for Probability of Direction *pd*

The function accepts either a vector of pre-computed *pd* values or a
matrix of posterior draws, from which *pd* values are computed
internally. Both direction-agnostic and directional tests are supported:
the `direction` argument controls which formulation is applied per
hypothesis. The adjustment is governed by the per-test prior probability
that an individual effect is null, \\\pi_0\\ (argument `pi0`).
Equivalently, one may supply the global null \\\Pr(H_0)\\ (argument
`p0`), the prior that all \\m\\ hypotheses are null, from which \\\pi_0
= \Pr(H_0)^{1/m}\\ is recovered under independence.

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
  alpha = NULL,
  fwer = FALSE
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
  hypotheses are null simultaneously (the global null, \\\Pr(H_0)\\).
  Supply **either** `pi0` **or** `p0`, not both. When `p0` is given, the
  per-test prior is derived under independence as \\\pi_0 =
  \Pr(H_0)^{1/m}\\.

- pi0:

  Numeric scalar in \\(0, 1)\\. The per-test prior probability that an
  individual effect is null, applied to every test. Supply **either**
  `pi0` **or** `p0`, not both.

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
  columns of `draws` (or the elements of `pd`).

- threshold:

  Numeric scalar in \\(0.5, 1)\\. The decision cutoff applied to *pd*
  (default `0.975`, a two-sided per-test \\\alpha = 0.05\\). It does
  **not** change `pd.adj`; it is used only to report the equivalent
  **adjusted threshold on the raw** *pd*, since rejecting \\pd\_{adj} \>
  c\\ is identical to rejecting raw \\pd \> c\_\*\\ with \\c\_\* =
  c\\\pi_0 / ((1-\pi_0)(1-c) + c\\\pi_0)\\.

- alpha:

  Optional numeric scalar in \\(0, 1)\\. A convenience parameterization
  of `threshold` on the error scale: when supplied, the two-sided cutoff
  is set to \\c = 1 - \alpha/2\\ (so `alpha = 0.05` gives
  `threshold = 0.975`), **overriding** any `threshold` value.

- fwer:

  Logical (default `FALSE`). If `TRUE`, the output also reports a
  **family-wise** cutoff on the raw *pd* that controls the
  prior-averaged FWER at the nominal level over the \\m\\ tests, given
  \\\pi_0\\. Standard corrections (e.g. Bonferroni) implicitly assume
  every test is null (\\\pi_0 = 1\\); here the affordable per-test rate
  is \\\[1 - (1 - \alpha)^{1/m}\] / \pi_0\\, relaxing the standard
  correction by \\1/\pi_0\\ and reducing to it as \\\pi_0 \to 1\\. This
  is a prior-averaged (Bayesian) FWER under the assumed \\\pi_0\\,
  **not** strong frequentist control.

## Value

An object of class `pd_adjust` (a `data.frame` with a dedicated
[`print.pd_adjust`](https://mar-cald.github.io/multibayes/reference/print.pd_adjust.md)
method), with one row per hypothesis. Columns: `pd` (values used in the
adjustment), `pd.adj` (adjusted values), `pi0` (per-test null prior
\\\pi_0\\), and `m` (number of tests). When `order = TRUE` the rows are
sorted by decreasing `pd`. For direction-agnostic tests, both `pd` and
`pd.adj` are bounded in \\\[0.5, 1\]\\; for directional tests, both are
on \\\[0, 1\]\\, with values below \\0.5\\ indicating that the data (and
the adjustment) favoured the opposite direction. When `draws` are
supplied, the output additionally includes `median.est` (posterior
median per parameter), `null.value`, and `direction`. The `print` method
summarises the constant quantities (`pi0`, `m`) in a header and displays
the per-test table; the columns themselves remain available for
programmatic access. The object also carries, as attributes, the
decision cutoff `threshold` (\\c\\), the equivalent adjusted cutoff on
the raw *pd* `threshold.adj` (\\c\_\*\\), and the corresponding
two-sided per-test Type I rates `alpha.nominal` (\\2(1-c)\\) and
`alpha.eff` (\\2(1-c\_\*)\\); retrieve them with, e.g.,
`attr(x, "threshold.adj")`. The cutoffs `threshold.adj` and
`threshold.fwer` are exact transformations of the decision rule, but the
reported per-test *rates* (`alpha.eff`, and the `fwer` calibration) use
the identity \\\alpha = 2(1 - c)\\, which holds **only for a diffuse
within-model prior** or **large \\n\\** (when the null *pd* is
approximately uniform on \\\[0.5, 1\]\\). Under an informative prior the
true per-test rate is smaller, so these rates are conservative and
`threshold.fwer` controls the FWER at or below the stated level.

## Details

The adjustment follows from Bayes' theorem. Given a per-hypothesis prior
\\\pi_0\\ and its complement \\\pi_1 = 1 - \pi_0\\, the adjusted *pd*
is: \$\$ pd\_{adj} = \frac{pd \pi_1}{pd \pi_1 + (1-pd)\pi_0} \$\$

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

`pd.adjust` returns the per-test correction only (the marginal `pd` and
its adjusted value `pd.adj`). Family-wise summaries across hypotheses,
the cumulative joint statement and simultaneous credible intervals, are
provided separately by
[`joint`](https://mar-cald.github.io/multibayes/reference/joint.md).

## References

Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A bayesian
perspective on the bonferroni adjustment. *Biometrika*, 84(2), 419–427.

Jeffreys, H. (1938). Significance tests when several degrees of freedom
arise simultaneously. *Proceedings of the Royal Society of London.
Series A. Mathematical and Physical Sciences*, 165(921), 161–198.
https://doi.org/10.1098/rspa.1938.0052

Storey, J. D. (2003). The positive false discovery rate: A bayesian
interpretation and the q-value. *The Annals of Statistics*, 31(6),
2013–2035. https://doi.org/10.1214/aos/1074290335

Scott, J. G., & Berger, J. O. (2006). An exploration of aspects of
bayesian multiple testing. *Journal of Statistical Planning and
Inference*, 136(7), 2144–2162.

Scott, J. G., & Berger, J. O. (2010). Bayes and empirical-Bayes
multiplicity adjustment in the variable- selection problem. *The Annals
of Statistics*, 38(5), 2587–2619. https://doi.org/10.1214/10-AOS792

## Examples

``` r
# From a vector of pd values
pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763,
H5 = 0.891, H6 = 0.987)
pd.adjust(pd = pd_values, p0 = 0.4)
#> Warning: some pd.adj have been floored to 0.5.
#> Prior-odds adjusted probability of direction
#> Tests: 6  |  pi0 = 0.8584
#> Decision: raw pd > 0.9958 (adjusted from 0.975); effective per-test alpha
#> 0.00843 (from 0.05)
#> 
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
#> Tests: 6  |  pi0 = 0.8584
#> Decision: raw pd > 0.9958 (adjusted from 0.975); effective per-test alpha
#> 0.00843 (from 0.05)
#> 
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
#> Tests: 6  |  pi0 = 0.8584
#> Decision: raw pd > 0.9958 (adjusted from 0.975); effective per-test alpha
#> 0.00843 (from 0.05)
#> 
#> 
#>    median.est null.value     pd pd.adj direction
#> H1     0.9860          0 0.8358 0.5000 two.sided
#> H2    -0.1175          0 0.5350 0.5000 two.sided
#> H3     0.7761          0 0.7822 0.5000 two.sided
#> H4    -0.0089          0 0.5030 0.5000 two.sided
#> H5     1.9500          0 0.9762 0.8715 two.sided
#> H6     2.9866          0 0.9995 0.9970 two.sided

# Order the rows by strength of evidence
pd.adjust(draws = draws, p0 = 0.4, null.value = 0, order = TRUE)
#> Warning: some pd.adj have been floored to 0.5.
#> Prior-odds adjusted probability of direction
#> Tests: 6  |  pi0 = 0.8584
#> Decision: raw pd > 0.9958 (adjusted from 0.975); effective per-test alpha
#> 0.00843 (from 0.05)
#> 
#> 
#>    median.est null.value     pd pd.adj direction
#> H6     2.9866          0 0.9995 0.9970 two.sided
#> H5     1.9500          0 0.9762 0.8715 two.sided
#> H1     0.9860          0 0.8358 0.5000 two.sided
#> H3     0.7761          0 0.7822 0.5000 two.sided
#> H2    -0.1175          0 0.5350 0.5000 two.sided
#> H4    -0.0089          0 0.5030 0.5000 two.sided

# Family-wise summaries are obtained separately with joint()
joint(draws, null.value = 0)               # cumulative joint statement
#>    median.est null.value     pd direction joint_cum
#> H6     2.9866          0 0.9995 two.sided    0.9995
#> H5     1.9500          0 0.9762 two.sided    0.9758
#> H1     0.9860          0 0.8358 two.sided    0.8155
#> H3     0.7761          0 0.7822 two.sided    0.6398
#> H2    -0.1175          0 0.5350 two.sided    0.3420
#> H4    -0.0089          0 0.5030 two.sided    0.1735
joint(draws, interval = TRUE)              # simultaneous credible intervals
#>         lower          est    upper prob
#> H1 -1.5133941  0.985996121 3.595471 0.95
#> H2 -2.8220596 -0.117524024 2.471335 0.95
#> H3 -1.8565114  0.776129478 3.336798 0.95
#> H4 -2.5959143 -0.008943627 2.611432 0.95
#> H5 -0.7100979  1.950044477 4.707048 0.95
#> H6  0.4556620  2.986561857 5.535382 0.95

# Mix of directional and agnostic tests with parameter-specific nulls
pd.adjust(draws = draws, p0 = 0.2, null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
          direction = c("greater", "two.sided", "greater",
          "two.sided", "greater", "greater"))
#> Warning: some pd.adj have been floored to 0.5.
#> Prior-odds adjusted probability of direction
#> Tests: 6  |  pi0 = 0.7647
#> Decision: raw pd > 0.9922 (adjusted from 0.975); effective per-test alpha
#> 0.0157 (from 0.05)
#> 
#> 
#>    median.est null.value     pd pd.adj direction
#> H1     0.9860        0.2 0.7820 0.5246   greater
#> H2    -0.1175        0.0 0.5350 0.5000 two.sided
#> H3     0.7761        0.2 0.7230 0.4454   greater
#> H4    -0.0089        0.0 0.5030 0.5000 two.sided
#> H5     1.9500        0.5 0.9285 0.7998   greater
#> H6     2.9866        0.5 0.9948 0.9831   greater
```
