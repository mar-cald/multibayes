# Prior-odds adjustment for Probability of Direction (pd)

The function accepts either a vector of pre-computed *pd* values or a
matrix of posterior draws, from which *pd* values are computed
internally. Both direction-agnostic and directional tests are supported:
the `direction` argument controls which formulation is applied per
hypothesis. The global prior probability that all tested hypotheses are
null, \\q\\, is decomposed into a per-hypothesis prior \\P(H_0) =
q^{1/m}\\, where \\m\\ is the number of hypotheses or, when the
correlation structure among parameters is taken into account, the
effective number of tests \\m\_{\text{eff}}\\.

## Usage

``` r
pd.adjust(
  pd = NULL,
  draws = NULL,
  q = 0.4,
  null.value = 0,
  direction = NULL,
  R = NULL
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

- q:

  Numeric scalar in \\(0, 1)\\. The prior probability that **all**
  hypotheses are null simultaneously. Defaults to `0.4`.

- null.value:

  Numeric scalar or vector. The null (reference) value against which the
  posterior is evaluated. A scalar applies the same null to all
  parameters; a vector of length equal to `ncol(draws)` allows a
  different null per parameter. Ignored if `pd` is supplied directly.
  Defaults to `0`.

- direction:

  Integer vector of `1`, `-1`, or `0` (or `NULL`). Specifies the testing
  mode for each hypothesis: `1` for a positive directional test
  (\\\Pr(\theta \> \theta\_\text{null})\\), `-1` for a negative
  directional test (\\\Pr(\theta \< \theta\_\text{null})\\), and `0` for
  direction-agnostic testing (maximum over both sides). A scalar is
  recycled to match the number of parameters; mixed vectors allow
  different modes across hypotheses. Defaults to `NULL`
  (direction-agnostic for all parameters).

- R:

  Optional correlation information for computing \\m\_{\text{eff}}\\.
  Accepts `TRUE` (correlation estimated from `draws`), a numeric scalar
  (assumed uniform correlation applied to all parameter pairs), or a
  full correlation matrix. When provided, \\m\_{\text{eff}}\\ replaces
  the nominal \\m\\.

## Value

A `data.frame` with one row per hypothesis, containing: `pd` (values
used in the adjustment), `pd.adj` (adjusted values), `q` (prior
probability of the global null), and `m` (nominal or effective number of
tests). For direction-agnostic tests, both `pd` and `pd.adj` are bounded
in \\\[0.5, 1\]\\; for directional tests, both are on \\\[0, 1\]\\, with
values below \\0.5\\ indicating that the data (and the adjustment)
favoured the opposite direction. When `draws` are supplied, the output
additionally includes `mean.est` (posterior mean per parameter),
`null.value` (null reference values), and `direction`.

## Details

The adjustment follows from Bayes' theorem. Given a per-hypothesis prior
\\P(H_0) = q^{1/m}\\ and its complement \\P(H_1) = 1 - P(H_0)\\, the
adjusted *pd* is: \$\$ pd\_{adj} = \frac{pd \cdot P(H_1)}{pd \cdot
P(H_1) + (1 - pd) \cdot P(H_0)} \$\$

Because the prior is conservative (\\P(H_0) \> P(H_1)\\), the adjustment
always shrinks *pd* toward its lower bound.

**Direction-agnostic tests** (`direction = 0`): *pd* is defined as
\\\max\\\big(\Pr(\hat\theta \> \theta\_\text{null}),\\ \Pr(\hat\theta \<
\theta\_\text{null})\big)\\ and is bounded in \\\[0.5, 1\]\\ by
construction. \\pd\_{adj}\\ is also floored at \\0.5\\, so the
adjustment produces shrinkage toward \\0.5\\.

**Directional tests** (`direction = 1` or `-1`): *pd* is the raw
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

When `R` is supplied, the effective number of tests \\m\_{\text{eff}}\\
is estimated from the eigenvalues \\\lambda\\ of the correlation matrix
(Cheverud, 2001): \$\$ m\_{\text{eff}} = K \left( 1 -
\frac{(K-1)\\\text{Var}(\lambda)}{K^2} \right) \$\$ where \\K\\ is the
number of hypotheses.

## References

Jeffreys, H. (1938). Significance tests when several degrees of freedom
arise simultaneously. *Proceedings of the Royal Society of London.
Series A. Mathematical and Physical Sciences, 165*(921), 161–198.
<https://doi.org/10.1098/rspa.1938.0052>

Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A Bayesian
perspective on the Bonferroni adjustment. *Biometrika, 84*(2), 419–427.
<https://doi.org/10.2307/2337467>

Cheverud, J. (2001). A simple correction for multiple comparisons in
interval mapping genome scans. *Heredity, 87*, 52–58.
<https://doi.org/10.1046/j.1365-2540.2001.00901.x>
