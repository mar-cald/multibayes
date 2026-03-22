# multibayes

**multibayes** provides tools for Bayesian multiple testing, currently
centred on prior-odds adjustment of the Probability of Direction (*pd*).

## Installation

``` r
# install.packages("remotes")
remotes::install_github("mar-cald/multibayes")
```

## Overview

| Function                                                                      | Input                          | Correction type                |
|-------------------------------------------------------------------------------|--------------------------------|--------------------------------|
| [`pd.adjust()`](https://mar-cald.github.io/multibayes/reference/pd.adjust.md) | Posterior draws or *pd* vector | Prior-odds adjustment for *pd* |

------------------------------------------------------------------------

## `pd.adjust()`: Prior-Odds Adjustment for *pd*

The Probability of Direction (*pd*) is often used as a decision rule: a
hypothesis is accepted when *pd* exceeds a threshold. When this decision
is made repeatedly across several hypotheses, the family-wise error
probability can be much higher than intended.
[`pd.adjust()`](https://mar-cald.github.io/multibayes/reference/pd.adjust.md)
addresses this decision-level problem by making the implicit prior
probability of the global null explicit, following the prior-odds
framework of Jeffreys (1938) and Westfall et al. (1997).

The global prior probability that **all** hypotheses are null, *q*, is
decomposed into a per-hypothesis prior:

$$P\left( H_{0} \right) = q^{1/m}$$

Each *pd* is then reweighted by Bayes’ theorem:

$$pd_{\text{adj}} = \frac{pd \cdot P\left( H_{1} \right)}{pd \cdot P\left( H_{1} \right) + (1 - pd) \cdot P\left( H_{0} \right)}$$

Because the prior is conservative
($P\left( H_{0} \right) > P\left( H_{1} \right)$), the adjustment always
moves *pd* toward 0; a floor at 0.5 is applied so the effective result
is shrinkage toward 0.5. When parameters are correlated, the effective
number of tests $m_{\text{eff}}$ (Cheverud, 2001) is used in place of
*m*, producing a less conservative adjustment.

### Usage

``` r
library(multibayes)

# From a vector of pd values (independence assumed)
pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763, H5 = 0.891, H6 = 0.987)
pd.adjust(pd = pd_values, q = 0.4)

# Simulate correlated posterior draws
Sigma <- matrix(0.4, nrow = 6, ncol = 6); diag(Sigma) <- 1
mu    <- c(1, -0.1, 0.8, 0, 2, 3)
draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")

# From posterior draws: pd and correlation estimated automatically
pd.adjust(draws = draws, q = 0.4, mu0 = 0, R = TRUE)

# Different null value per parameter (e.g., minimum effect of practical interest)
pd.adjust(draws = draws, q = 0.4, mu0 = c(0.2, 0, 0.2, 0, 0.5, 0.5),
          direction = c(1, 0, 1, 0, 1, 1), R = TRUE)

# When draws are unavailable, supply an assumed mean correlation
pd.adjust(pd = pd_values, q = 0.4, R = 0.4)
```

### Output

When `draws` are supplied, the output is a `data.frame` with one row per
hypothesis:

| Column      | Description                                                                                                                                                                                    |
|-------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `mean_est`  | Posterior mean per parameter                                                                                                                                                                   |
| `mu0`       | Null reference value used                                                                                                                                                                      |
| `direction` | Declared direction (`1`, `-1`, or `0` for agnostic)                                                                                                                                            |
| `pd`        | *pd* used in the adjustment, floored at 0.5                                                                                                                                                    |
| `pd_adj`    | Adjusted *pd* after prior-odds correction, floored at 0.5                                                                                                                                      |
| `pd_raw`    | Raw directional probability before flooring (only for directional tests; `NA` for agnostic tests); values below 0.5 indicate the posterior is concentrated opposite to the predicted direction |
| `q`         | Global null probability used                                                                                                                                                                   |
| `m`         | Family size used (*m* or $m_{\text{eff}}$)                                                                                                                                                     |

When a `pd` vector is supplied directly, only `pd`, `pd_adj`, `q`, and
`m` are returned.

### Choosing *q*, *mu0*, and *direction*

- **`q`** encodes your prior belief that all tested hypotheses are
  simultaneously null. A value of `0.4` is a skeptical default.
- **`mu0`** sets the null reference value for each parameter. A scalar
  applies the same null to all parameters; a vector allows a different
  null per parameter, for instance when testing against a minimum effect
  of practical interest.
- **`direction`** specifies the expected sign of each effect (`1` for
  positive, `-1` for negative, `0` for agnostic). When specified, *pd*
  is the probability mass on the predicted side and is floored at 0.5;
  the raw probability before flooring is returned in `pd_raw`. A value
  of `pd = 0.5` indicates absence of support for the predicted
  direction; `pd_raw` allows the researcher to assess whether the floor
  was triggered and how strongly the data contradicted the prediction.

------------------------------------------------------------------------

## Citation

If you use **multibayes** in published research, please cite:

> Calderan, M., Gambarota, F., Toffalini, E., & Altoè, G. (2026).
> Adjusting the Probability of Direction for Multiple Testing.

------------------------------------------------------------------------

## References

- Jeffreys, H. (1938). Significance tests when several degrees of
  freedom arise simultaneously. *Proceedings of the Royal Society of
  London. Series A. Mathematical and Physical Sciences, 165*(921),
  161–198. <https://doi.org/10.1098/rspa.1938.0052>
- Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A Bayesian
  perspective on the Bonferroni adjustment. *Biometrika, 84*(2),
  419–427. <https://doi.org/10.2307/2337467>
- Cheverud, J. (2001). A simple correction for multiple comparisons in
  interval mapping genome scans. *Heredity, 87*, 52–58.
  <https://doi.org/10.1046/j.1365-2540.2001.00901.x>

------------------------------------------------------------------------

## License

GNU GENERAL PUBLIC LICENSE, Version 3
