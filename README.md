# multibayes

**multibayes** provides tools for Bayesian inference and decision-making across multiple hypotheses simultaneously. 
The package offers two functions rooted in different but related principles: `pd.adjust()` addresses 
the decision-level multiplicity problem that arises when the Probability of 
Direction (*pd*) is used as a repeated decision rule, while `joint()` provides 
simultaneous credible intervals that are faithful to the joint posterior distribution, 
enabling inference about several parameters at once.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mar-cald/multibayes")
```

## Overview

| Function | Input | Correction type |
|---|---|---|
| `pd.adjust()` | Posterior draws or *pd* vector | Prior-odds adjustment for *pd* |
| `joint()` | Posterior draws | Simultaneous credible intervals |

Both functions operate directly on posterior draws from different modeling framework
(e.g., `brms`, `rstan`, `rstanarm`).

---

## `pd.adjust()`: Prior-Odds Adjustment for *pd*

The Probability of Direction (*pd*) is often used as a decision rule: a hypothesis 
is accepted when *pd* exceeds a threshold. When this decision is made repeatedly 
across several hypotheses, the family-wise error probability can be much higher 
than intended. `pd.adjust()` addresses this decision-level problem by making the 
implicit prior probability of the global null explicit, following the prior-odds 
framework of Jeffreys (1938) and Westfall et al. (1997).

The global prior probability that **all** hypotheses are null, *q*, is
decomposed into a per-hypothesis prior:

$$P(H_0) = q^{1/m}$$

Each *pd* is then reweighted by Bayes' theorem:

$$pd_{\text{adj}} = \frac{pd \cdot P(H_1)}{pd \cdot P(H_1) + (1 - pd) \cdot P(H_0)}$$

When parameters are correlated, the effective number of tests $m_\text{eff}$
(Cheverud, 2001) is used in place of *m*, producing a less conservative
adjustment.

### Usage

```r
library(multibayes)

# From a vector of pd values (independence assumed)
pd.adjust(pd = c(0.99, 0.87, 0.76, 0.97, 0.83, 0.98), q = 0.4)

# From posterior draws, pd computed internally, correlation estimated automatically
draws <- matrix(rnorm(6000), ncol = 6)
pd.adjust(draws = draws, q = 0.4, R = TRUE)

# Supplying a scalar correlation (e.g., assumed mean r = 0.3)
pd.adjust(pd = c(0.92, 0.87, 0.76, 0.95, 0.83, 0.91), q = 0.4, R = 0.3)
```

### Output

A `data.frame` with one row per hypothesis:

| Column | Description |
|---|---|
| `pd` | Original unadjusted *pd* |
| `pd_adj` | Adjusted *pd* after prior-odds correction |
| `q` | Global null probability used |
| `m` | Family size used (*m* or $m_\text{eff}$) |

### Choosing *q* and *m*

- **`q`** encodes your prior belief that all tested hypotheses are
  simultaneously null. A value of `0.4` is a moderately skeptical default.
  Lower values are more lenient; higher values are more conservative.
- **`m`** should reflect only the hypotheses for which directional claims
  are made. Nuisance parameters (random effects, control covariates) should
  not be counted.

---

## `joint()`: Simultaneous Credible Intervals

`joint()` computes equitailed credible intervals that hold **jointly** across
all parameters: the posterior probability that every parameter lies within its
interval *simultaneously* equals the requested level. This is stricter than
the standard marginal credible interval, which guarantees coverage only
parameter-by-parameter.

Simultaneous coverage is calibrated by finding the critical quantile *cq*
such that only a fraction $\alpha$ of posterior draws are more extreme than
*cq* in at least one parameter simultaneously.

### Usage

```r
library(MASS)

mu    <- c(4, 0, -2)
Sigma <- matrix(c(1, 0.8, 0.5,
                  0.8, 1, 0.3,
                  0.5, 0.3, 1), 3, 3)
draws <- mvrnorm(2000, mu, Sigma)
colnames(draws) <- c("theta1", "theta2", "theta3")

joint(draws, prob = 0.95)
```

### Output

A named list:

| Element | Description |
|---|---|
| `lower` | Lower bounds, one per parameter |
| `upper` | Upper bounds, one per parameter |
| `est` | Point estimates (default: median) |
| `prob` | Requested joint posterior probability |
| `cq` | Critical quantile used for calibration |

---

## Citation

If you use **multibayes** in published research, please cite:

For *pd.adjust()*:

- Calderan, M., Gambarota, F., Toffalini, E., & Altoè, G. (2026). Adjusting the Probability of Direction for Multiple Testing.

For *joint()*:

- Goeman, J., Calderan, M., Solari, A. (2026). Bonferroni for Bayesian: Multiplicity Correction Based on Joint Credibility.

---

## References

- Jeffreys, H. (1938). Significance tests when several degrees of freedom arise simultaneously.Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences,165(921), 161–198.
  <https://doi.org/10.1098/rspa.1938.0052>
- Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A Bayesian Perspective on the Bonferroni Adjustment.
  Biometrika, 84(2), 419–427. <http://www.jstor.org/stable/2337467>.
- Cheverud, J. A simple correction for multiple comparisons in interval mapping genome scans. Heredity 87, 52–58 (2001). <https://doi.org/10.1046/j.1365-2540.2001.00901.x>.


## License

