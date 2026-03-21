# Prior-odds adjustment for Probability of Direction (pd)

The function accepts either a vector of pre-computed *pd* values or a
matrix of posterior draws, from which *pd* values are computed
internally. The global prior probability that all tested hypotheses are
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
  mu0 = 0,
  direction = NULL,
  R = NULL
)
```

## Arguments

- pd:

  Numeric vector of *pd* values in \\\[0.5, 1\]\\.

- draws:

  Optional matrix or data frame of posterior draws (columns =
  parameters). If provided, `pd` is calculated automatically.

- q:

  Numeric scalar in \\(0, 1)\\. The prior probability that **all**
  hypotheses are null. Defaults to `0.4`.

- mu0:

  Numeric scalar or vector. The null (reference) value against which the
  posterior is evaluated. A scalar applies the same null to all
  parameters; a vector of length equal to `ncol(draws)` allows a
  different null per parameter. Ignored if `pd` is supplied directly.
  Defaults to `0`.

- direction:

  Integer vector of `1`, `-1`, or `0` (or `NULL`). Specifies the
  expected direction of each effect: `1` for positive (tests
  `draws > mu0`), `-1` for negative (tests `draws < mu0`), and `0` for
  direction-agnostic (takes the maximum over both sides). Should be a
  vector of length `ncol(draws)` to specify a different direction per
  parameter; a scalar is recycled across all parameters. Should be
  specified when `mu0 != 0`. Defaults to `NULL` (direction-agnostic for
  all parameters).

- R:

  Optional correlation information for computing \\m\_{\text{eff}}\\.
  Accepts `TRUE` (correlation estimated from `draws`), a numeric scalar
  (assumed uniform correlation applied to all parameter pairs), or a
  full correlation matrix. When provided, \\m\_{\text{eff}}\\ replaces
  the nominal \\m\\.

## Value

A `data.frame` with one row per hypothesis, containing the following
columns: `pd` (original values), `pd_adj` (adjusted values), `q` (prior
probability of the global null), and `m` (number of hypotheses or
effective number of tests). When `draws` are supplied, `mean_est`
(posterior mean per parameter), `mu0` (null reference values), and
`direction` (declared direction per hypothesis; `0` indicates
direction-agnostic) are also returned.

## Details

The adjustment follows from Bayes' theorem. Given a per-hypothesis prior
\\P(H_0) = q^{1/m}\\ and its complement \\P(H_1) = 1 - P(H_0)\\, the
adjusted *pd* is: \$\$ pd\_{adj} = \frac{pd \cdot P(H_1)}{pd \cdot
P(H_1) + (1 - pd) \cdot P(H_0)} \$\$

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
