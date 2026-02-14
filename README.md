# multibayes

Bayesian Multiplicity Adjustment via Prior Odds.

## Overview

`multibayes` implements a Bayesian multiplicity adjustment based on calibrating
hypothesis-specific prior odds from a global prior probability for the complete null.
In practice, the main functions adjusts `posterio odds` values by translating a global prior
probability that *all* tested hypotheses are null into a per-hypothesis prior null
probability that depends on the family size.

Optionally, the adjustment can account for dependence across tests via a posterior
correlation matrix, which is used to compute an effective family size.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("mar-cald/multibayes")
```

