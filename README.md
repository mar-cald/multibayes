# pdAdjust

Bayesian multiplicity adjustment using probability of direction (*pd*) values.

## Overview

The `pdAdjust` package provides a principled approach to Bayesian multiplicity adjustment based on prior-odds reweighting of probability of direction values (*pd*). 
The main function adjusts *pd* values using a global prior probability that all tested hypotheses are null, converting it into per-hypothesis prior odds that account for family size.


## Installation

You can install the development version of pdAdjust from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Mar-Cald/pdAdjust")
```
## Usage
```r
library(pdAdjust)
# Adjust a vector of pd values
pd <- c(0.55, 0.80, 0.97)
prior_adj(pd, q = 0.4)
#> [0.3039065 0.5882800 0.9203171]
```
The adjustment becomes more conservative as family size increases:
```r
# With 20 hypotheses
pd_large <- rep(0.97, 10)
prior_adj(pd_large, q = 0.4)

# Custom family size and q
prior_adj(c(0.90, 0.95, 0.97), q = 0.5, m = 2)
```

## Function Details

```r
prior_adj(pd, q = 0.4, m = length(pd))
```

**Arguments**:

-  **pd**: Numeric vector of pd values (in [0.5, 1]). Each entry represents the posterior mass on the favored direction.

-  **q**: Numeric scalar in (0, 1). Prior probability that all hypotheses in the family are null (default: 0.4).

-  **m**: Family size (default: length of pd vector).

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
R package version 0.1.0. https://github.com/Mar-Cald/pdAdjust

## License
MIT License

## Issues and Contributions
Please report issues at: https://github.com/Mar-Cald/pdAdjust/issues


