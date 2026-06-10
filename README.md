# multibayes

**multibayes** provides tools for **Bayesian multiplicity adjustments** — correcting and summarising decisions taken across many parameters at once. It currently offers a prior-odds adjustment of the Probability of Direction (*pd*) via `pd.adjust()`, and two family-wise summaries via `joint()`: a cumulative **joint statement** (default) and **simultaneous credible intervals** (`interval = TRUE`).
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
| `joint()` | Posterior draws | Family-wise: cumulative joint statement (default) or simultaneous credible intervals (`interval = TRUE`) |

---

## `pd.adjust()`: Prior-Odds Adjustment for *pd*

The Probability of Direction (*pd*) is often used as a decision rule: a hypothesis 
is accepted when *pd* exceeds a threshold. 

The primary input is the per-test prior $\pi_0$, the probability that an individual 
effect is null. Under exchangeability this is the expected proportion of null 
effects in the family. Equivalently, you may specify the global null, *p0* ($\Pr(H_0)$), 
the prior probability that **all** $m$ hypotheses are null, from which the per-test prior is recovered 
under independence:

$$\Pr(H_0) = \pi_0^{m}, \qquad \pi_0 = \Pr(H_0)^{1/m}.$$

Each *pd* is then reweighted by Bayes' theorem:

$$pd_{\text{adj}} = \frac{pd \cdot \pi_1}{pd \cdot \pi_1 + (1 - pd) \cdot \pi_0}$$

Because the prior is conservative ($\pi_0 > \pi_1$), the adjustment always 
shrinks *pd* toward its lower bound. 

The function supports two testing modes, which can be mixed across hypotheses within the same call:

- **Direction-agnostic** (`direction = "two.sided"`): *pd* = $\max\big(\Pr(\hat\theta > \theta_\text{null}),\, \Pr(\hat\theta < \theta_\text{null})\big)$, bounded in $[0.5, 1]$ by construction; $pd_\text{adj}$ is also floored at $0.5$.
- **Directional** (`direction = "greater"` or `"less"`): *pd* is the raw one-sided posterior probability on the predicted side, on $[0, 1]$. Values below $0.5$ indicate that the posterior is concentrated opposite to the predicted direction; the adjustment will further shrink such values, reflecting the combined weight of the data and the conservative prior against the hypothesis.

### Usage

```r
library(multibayes)

# From a vector of pd values (direction-agnostic)
pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763, H5 = 0.891, H6 = 0.987)

# Specify the per-test null prior pi0 directly (primary interface)
pd.adjust(pd = pd_values, pi0 = 0.65)

# ...or specify the global null p = Pr(all null); pi0 is then p^(1/m)
pd.adjust(pd = pd_values, p0 = 0.4)

# Simulate posterior draws
Sigma <- matrix(0, nrow = 6, ncol = 6); diag(Sigma) <- 1
mu <- c(1, -0.1, 0.8, 0, 2, 3)
draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")

# From posterior draws (per-test correction only)
pd.adjust(draws = draws, p0 = 0.4, null.value = 0)

# Mix of directional and agnostic tests with parameter-specific nulls
pd.adjust(
  draws = draws,
  pi0 = 0.4,
  null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
  direction  = c("greater", "two.sided", "greater", "two.sided", "greater", "greater")
)
```

### Output

`pd.adjust()` returns a `pd_adjust` object — a `data.frame` with a tidy `print` method that summarises the constant quantities (`pi0`, `m`) in a header and shows the per-test table. With `draws`, it has one row per hypothesis:

| Column | Description |
|---|------|
| `median.est` | Posterior median per parameter |
| `null.value` | Null reference value used |
| `pd` | *pd* used in the adjustment; in $[0.5, 1]$ for agnostic tests, $[0, 1]$ for directional tests |
| `pd.adj` | Adjusted *pd* after prior-odds correction; same bounds as `pd` |
| `pi0` | Per-test null prior $\pi_0$ |
| `m` | Family size used |
| `direction` | Testing mode: `"greater"`, `"less"`, or `"two.sided"` |


`pd.adjust()` returns the **per-test correction only**. **`order = TRUE`** sorts the rows by decreasing `pd`; when a `pd` vector is supplied directly, `median.est`, `null.value`, and `direction` are not returned. 

### The adjusted threshold (what your `pi0` implies)

Rejecting `pd.adj > c` is identical to rejecting the raw `pd > c*`, where

$$c^* = \frac{c\,\pi_0}{(1-\pi_0)(1-c) + c\,\pi_0}.$$

The print header reports this: e.g. with `pi0 = 0.7` and the default `c = 0.975`,

```
Decision: raw pd > 0.9891 (adjusted from 0.975);  effective per-test alpha 0.0217 (from 0.05)
```

> **Note.** The *cutoff* $c^*$ is an exact transformation of the rule, but the per-test *rate* (the `0.022`) uses $\alpha = 2(1-c)$, which holds for a diffuse prior or large $n$. Under an informative prior the true per-test rate is smaller, so the reported rate is conservative.

```r
pd.adjust(pd = pd_values, pi0 = 0.7, alpha = 0.05)   # same as threshold = 0.975
```

**Family-wise mode.** Set `fwer = TRUE` to additionally report a **family-wise** conservative cutoff on the raw *pd* that controls the (prior-averaged $\pi_0$) FWER over the *m* tests. 
Standard corrections assume every test is null ($\pi_0 = 1$); using your $\pi_0$ relaxes that by $1/\pi_0$. It is a prior-averaged (Bayesian) FWER, *not* strong frequentist control.

```r
pd.adjust(pd = pd_values, pi0 = 0.7, fwer = TRUE)
#> Decision:    raw pd > 0.9891 (adjusted from 0.975); effective per-test alpha 0.0217 (from 0.05)
#> Family-wise: raw pd > 0.9959 (FWER 0.05 over m = ... under the assumed pi0; per-test alpha 0.0081)
```

### Choosing `pi0` (or `p0`), `null.value`, and `direction`

- **`pi0`** (primary) is your prior probability that any single effect is null — 
equivalently, the expected proportion of null effects in the family. 
Alternatively, **`p0`** encodes the prior belief that **all** tested hypotheses 
are simultaneously null, with $\Pr(H_0) = \pi_0^m$. Supply one of the two.
- **`null.value`** sets the null reference value for each parameter. 
A scalar applies the same null to all parameters; a vector allows a different 
null per parameter.
- **`direction`** selects the testing mode per hypothesis 
(`"greater"`, `"less"`, `"two.sided"`). A scalar is recycled across 
all hypotheses; a mixed vector applies different modes within the same call. 
For directional tests, values of `pd` and `pd.adj` below `0.5` 
indicate that the data contradicted the predicted direction; 
the adjustment amplifies this evidence under the conservative prior.

---

## `joint()`: Family-wise statements

**Cumulative joint statement (default).** The posterior probability that the *k* strongest directional claims all hold *at once*. 
Rows are ordered strongest-first; `joint_cum` accumulates down the column, and `1 − joint_cum` is the posterior probability that *any* of those *k* claims is mis-signed (a family-wise, sign-error / Type S, statement). 
It is computed from the draws, so it respects their dependence.

```r
joint(draws, null.value = 0)
#>        median.est null.value     pd direction joint_cum
#> theta1     ...            0    ...  two.sided      ...
#> ...
```

**Simultaneous credible intervals (`interval = TRUE`).** Equitailed intervals such that **all** parameters lie within their bounds *simultaneously* with probability `prob` (the quantile-based simultaneous credible band; closely related to `credsubs::sim.cred.band()`, Schnell et al., 2020). 
Declare an effect when its band excludes the null value.

```r
joint(draws, interval = TRUE, prob = 0.95)
#>         lower    est  upper prob
#> theta1  ...      ...  ...   0.95
#> ...
```

Both are pure posterior quantities — they do **not** use `pi0`. 

---

## Citation

If you use **multibayes** in published research, please cite the **article** (preferred over citing the software itself):

> Calderan, M., Gambarota, F., Finos, L., & Altoè, G. (2026). 
A Prior-Odds Adjustment for the Probability of Direction in Multiple Testing

Run `citation("multibayes")` in R for the full reference (article and package entries).

---

## License

GNU General Public License, Version 3
