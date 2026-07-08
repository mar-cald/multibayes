# multibayes

**multibayes** provides tools for **Bayesian multiplicity adjustments**. 
Direction-agnostic and directional tests are supported, and can be mixed freely within a single call.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mar-cald/multibayes")
```

## Overview

| Function | Input | Correction type |
|---|---|---|
| `pd.adjust()` | Posterior draws or *pd* vector | Prior-odds adjustment for *pd* (per-test) |
| `joint()` | Posterior draws | Joint probability (directional) statements (default) or simultaneous credible intervals (`interval = TRUE`) |



## Citation

> Calderan, M., Gambarota, F., Finos, L., & Altoè, G. (n.d.). A Prior-Odds Adjustment for the Probability of Direction in Multiple Testing. Retrieved from osf.io/preprints/psyarxiv/8zwx2_v1

Run `citation("multibayes")` in R for the full reference (article and package entries).

---

## License

GNU General Public License, Version 3
