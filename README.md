# multibayes

**multibayes** provides tools for Bayesian multiple testing, currently centred on prior-odds adjustment of the Probability of Direction (*pd*). 
Both direction-agnostic and directional tests are supported, and can be mixed freely within a single call.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mar-cald/multibayes")
```

## Overview

| Function | Input | Correction type |
|---|---|---|
| `pd.adjust()` | Posterior draws or *pd* vector | Prior-odds adjustment for *pd* (agnostic or directional) |

***

## `pd.adjust()`: Prior-Odds Adjustment for *pd*

The Probability of Direction (*pd*) is often used as a decision rule: a hypothesis is accepted when *pd* exceeds a threshold. When this decision is made repeatedly across several hypotheses, the family-wise error probability can be much higher than intended. `pd.adjust()` addresses this problem by making the implicit prior probability of the global null explicit, following the prior-odds framework of Jeffreys (1938) and Westfall et al. (1997).

The global prior probability that **all** hypotheses are null, *q*, is decomposed into a per-hypothesis prior:

$$P(H_0) = q^{1/m}$$

Each *pd* is then reweighted by Bayes' theorem:

$$pd_{\text{adj}} = \frac{pd P(H_1)}{pd P(H_1) + (1 - pd) P(H_0)}$$

Because the prior is conservative ( $P(H_0) > P(H_1)$ ), the adjustment always shrinks *pd* toward its lower bound. When parameters are correlated, the effective number of tests $m_\text{eff}$ (Cheverud, 2001) is used in place of *m*, producing a less conservative adjustment.

The function supports two testing modes, which can be mixed across hypotheses within the same call:

- **Direction-agnostic** (`direction = "two.sided"`): *pd* = $\max\big(\Pr(\hat\theta > \theta_\text{null}),\, \Pr(\hat\theta < \theta_\text{null})\big)$, bounded in $[0.5, 1]$ by construction; $pd_\text{adj}$ is also floored at $0.5$.
- **Directional** (`direction = "greater"` or `"less"`): *pd* is the raw one-sided posterior probability on the predicted side, on $[0, 1]$. Values below $0.5$ indicate that the posterior is concentrated opposite to the predicted direction; the adjustment will further shrink such values, reflecting the combined weight of the data and the conservative prior against the hypothesis.

### Usage

```r
library(multibayes)

# From a vector of pd values (independence assumed, direction-agnostic)
pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763, H5 = 0.891, H6 = 0.987)
pd.adjust(pd = pd_values, q = 0.4)

# Simulate correlated posterior draws
Sigma <- matrix(0.4, nrow = 6, ncol = 6); diag(Sigma) <- 1
mu    <- c(1, -0.1, 0.8, 0, 2, 3)
draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")

# From posterior draws: pd and correlation estimated automatically
pd.adjust(draws = draws, q = 0.4, null.value = 0, R = TRUE)

# Mix of directional and agnostic tests with parameter-specific nulls
pd.adjust(draws = draws, q = 0.4, null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
           direction = c("greater", "two.sided", "greater", "two.sided", 
           "greater", "greater"), R = TRUE)

# When draws are unavailable, supply an assumed mean correlation
pd.adjust(pd = pd_values, q = 0.4, R = 0.4)
```

### Output

When `draws` are supplied, the output is a `data.frame` with one row per hypothesis:

| Column | Description |
|---|---|
| `mean.est` | Posterior mean per parameter |
| `null.value` | Null reference value used |
| `direction` | Testing mode: `greater`, `less`, `two.sided` |
| `pd` | *pd* used in the adjustment; in $[0.5, 1]$ for agnostic tests, $[0, 1]$ for directional tests |
| `pd.adj` | Adjusted *pd* after prior-odds correction; same bounds as `pd` |
| `q` | Global null probability used |
| `m` | Family size used (*m* or $m_\text{eff}$) |

When a `pd` vector is supplied directly, only `pd`, `pd.adj`, `q`, and `m` are returned.

### Choosing *q*, *null.value*, and *direction*

- **`q`** encodes your prior belief that all tested hypotheses are simultaneously null. A value of `0.4` is a skeptical default.
- **`null.value`** sets the null reference value for each parameter. A scalar applies the same null to all parameters; a vector allows a different null per parameter, for instance when testing against a minimum effect of practical interest.
- **`direction`** selects the testing mode per hypothesis (`greater`, `less`, `two.sided`). A scalar is recycled; a mixed vector applies different modes across hypotheses within the same call. For directional tests, values of `pd` and `pd.adj` below `0.5` indicate that the data contradicted the predicted direction; the adjustment amplifies this evidence under the conservative prior.

---

## Citation

If you use **multibayes** in published research, please cite:

> Calderan, M., Gambarota, F., Toffalini, E., & Altoè, G. (2026). Adjusting the Probability of Direction for Multiple Testing.

---

## References

- Jeffreys, H. (1938). Significance tests when several degrees of freedom arise simultaneously. *Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences, 165*(921), 161–198. <https://doi.org/10.1098/rspa.1938.0052>
- Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A Bayesian perspective on the Bonferroni adjustment. *Biometrika, 84*(2), 419–427. <https://doi.org/10.2307/2337467>
- Cheverud, J. (2001). A simple correction for multiple comparisons in interval mapping genome scans. *Heredity, 87*, 52–58. <https://doi.org/10.1046/j.1365-2540.2001.00901.x>

---

## License

GNU GENERAL PUBLIC LICENSE, Version 3