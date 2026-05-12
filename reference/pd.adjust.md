# Prior-odds adjustment for Probability of Direction (pd)

The function accepts either a vector of pre-computed *pd* values or a
matrix of posterior draws, from which *pd* values are computed
internally. Both direction-agnostic and directional tests are supported:
the `direction` argument controls which formulation is applied per
hypothesis. The global prior probability that all tested hypotheses are
null, \\q\\, is decomposed into a per-hypothesis prior \\P(H_0) =
q^{1/m}\\, where \\m\\ is the number of hypotheses.

## Usage

``` r
pd.adjust(pd = NULL, draws = NULL, q = NULL, null.value = 0, direction = NULL)
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

- q:

  Numeric scalar in \\(0, 1)\\. The prior probability that **all**
  hypotheses are null simultaneously.

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

## Value

A `data.frame` with one row per hypothesis, containing: `pd` (values
used in the adjustment), `pd.adj` (adjusted values), `q` (prior
probability of the global null), and `m` (number of tests), `qm` (prior
probability for the null for each test). For direction-agnostic tests,
both `pd` and `pd.adj` are bounded in \\\[0.5, 1\]\\; for directional
tests, both are on \\\[0, 1\]\\, with values below \\0.5\\ indicating
that the data (and the adjustment) favoured the opposite direction. When
`draws` are supplied, the output additionally includes `mean.est`
(posterior mean per parameter), `null.value` (null reference values),
and `direction`.

## Details

The adjustment follows from Bayes' theorem. Given a per-hypothesis prior
\\P(H_0) = q^{1/m}\\ and its complement \\P(H_1) = 1 - P(H_0)\\, the
adjusted *pd* is: \$\$ pd\_{adj} = \frac{pd P(H_1)}{pd P(H_1) + (1 - pd)
P(H_0)} \$\$

Because the prior is conservative (\\P(H_0) \> P(H_1)\\), the adjustment
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

## References

Jeffreys, H. (1938). Significance tests when several degrees of freedom
arise simultaneously. *Proceedings of the Royal Society of London.
Series A. Mathematical and Physical Sciences, 165*(921), 161–198.
<https://doi.org/10.1098/rspa.1938.0052>

Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A Bayesian
perspective on the Bonferroni adjustment. *Biometrika, 84*(2), 419–427.
<https://doi.org/10.2307/2337467>

## Examples

``` r
if (FALSE) { # \dontrun{
# From a vector of pd values 
pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763, 
H5 = 0.891, H6 = 0.987)
pd.adjust(pd = pd_values, q = 0.4)

# Simulate draws
Sigma <- matrix(0, nrow = 6, ncol = 6); diag(Sigma) <- 1
mu    <- c(1, -0.1, 0.8, 0, 2, 3)
draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")

# From posterior draws
pd.adjust(draws = draws, q = 0.4, null.value = 0)

# Mix of directional and agnostic tests with parameter-specific nulls
pd.adjust(draws = draws, q = 0.2, null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
          direction = c("greater", "two.sided", "greater", 
          "two.sided", "greater", "greater"), R = TRUE)
} # }
```
