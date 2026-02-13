# pdAdjust

Bayesian multiplicity adjustment using probability of direction (*pd*) values.

## Overview

he `pdAdjust` package provides a principled approach to Bayesian multiplicity adjustment based on prior-odds reweighting of probability of direction values (*pd*). 
The main function adjusts *pd* values using a global prior probability that all tested hypotheses are null, converting it into per-hypothesis prior odds that account for family size. It also supports adjusting for dependent hypotheses using a posterior correlation matrix.

## Installation

You can install the development version of pdAdjust from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Mar-Cald/pdAdjust")
```
## Usage
```r
library(pdAdjust)

# 1. Standard adjustment
pd <- c(0.985, 0.999, 0.975)
prior_adj(pd, q = 0.4)
#> [1] 0.9591114 0.9972055 0.9330259

# 2. Adjustment with a costum q
prior_adj(pd, q = 0.2)
#> [1] 0.9790012 0.9985921 0.9651435

# 3. Adjustment accounting for dependence
corr_mat <- matrix(c(1.0, 0.5, 0.3,
                     0.5, 1.0, 0.4,
                     0.3, 0.4, 1.0), nrow = 3)
prior_adj(pd, q = 0.4, post_corr = corr_mat)
#> [1] 0.9641900 0.9975646 0.9411455
```

## Function Details

```r
prior_adj(pd, q = 0.4, m = length(pd), post_corr = NULL)
```

**Arguments**:

-  **pd**: Numeric vector of pd values (in [0.5, 1]). Each entry represents the posterior mass on the favored direction.

-  **q**: Numeric scalar in (0, 1). Prior probability that all hypotheses in the family are null (default: 0.4).

-  **m**: Family size (default: length of pd vector).

- **post_corr**: Optional correlation matrix (NULL by default). If provided, it must be a square matrix matching the length of *pd*. It is used to calculate an effective family size, accounting for dependence between tests.

**Returns**:

Numeric vector of adjusted pd values.

**Details**:

The function computes per-hypothesis prior odds as $\Pr(H_{0_i}) = q^{1/m}$, where $m$ is the family size. 
If $\Pr(H_{0_i}) \geq 0.5$, each *pd* value is adjusted as:

$pd_{adj_i} = (pd_i \cdot \Pr(H_{1_i})) / (pd_i \cdot \Pr(H_{1_i}) + (1 - pd_i) \cdot \Pr(H_{0_i}))$

where $\Pr(H_{1_i}) = 1 - \Pr(H_{0_i})$.

Otherwise original *pd* values are reported.

## Citation
If you use this package in your research, please cite:

Calderan, M., Gambarota, F., Toffalini, E., & Altoè, G. (2026). Multiple Probabilities of Direction: A Rationale for Prior-odds Adjustment,
https://github.com/Mar-Cald/pdAdjust

## License
MIT License

## Issues and Contributions
Please report issues at: https://github.com/Mar-Cald/pdAdjust/issues


