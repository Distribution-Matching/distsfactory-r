# distsfactory

An R package for constructing probability distributions from partial specifications — moments, quantiles, mode, and more.

Part of the DistributionsFactories family (alongside [DistributionsFactories.jl](https://github.com/Ron-Ash/DistributionsFactories.jl) for Julia and [distsfactory-python](https://github.com/yoninazarathy/distsfactory-python)).

## Design

- **Zero heavy dependencies** — works with base R's built-in distribution functions (`dgamma`, `pnorm`, `qbeta`, etc.)
- Specify what you know (mean, variance, quantiles, mode) and get back a distribution object
- The returned object wraps base R's `d/p/q/r` functions with parameters baked in
- Raw parameters are always accessible via `$params`

## Quick start

```r
library(distsfactory)

# Construct from moments
d <- make_dist("gamma", mean = 5, var = 3)

# Use the built-in d/p/q/r methods (mirrors base R convention)
d$d(2)          # density:  dgamma(2, shape=..., rate=...)
d$p(5)          # CDF:      pgamma(5, shape=..., rate=...)
d$q(0.5)        # quantile: qgamma(0.5, shape=..., rate=...)
d$r(100)        # random:   rgamma(100, shape=..., rate=...)

# Access raw parameters for use with base R directly
d$params        # list(shape = 8.33, rate = 1.67)
```

## Supported distributions

| Distribution | Supported specifications |
|---|---|
| **Gamma** | mean+var, mean+mode, mode+var, mode+quantile, mode+iqr, two quantiles, mean+quantile |
| **Exponential** | mean, mean+var, var, single quantile |
| **Logistic** | mean+var, two quantiles, mode+iqr, mean+quantile |
| **Beta** | mean+var, mean+mode, two quantiles, mean+quantile |

## Specification styles

```r
# Moment-based
make_dist("gamma", mean = 5, var = 3)
make_dist("gamma", mean = 5, std = 2)          # std -> var
make_dist("gamma", mean = 5, cv = 0.5)         # coefficient of variation
make_dist("gamma", mean = 4, scv = 0.5)        # squared CV
make_dist("exponential", mean = 3)              # 1-parameter family

# Quantile-based
make_dist("exponential", median = 2.0)
make_dist("logistic", q1 = 2, q3 = 8)
make_dist("gamma", quantiles = list(c(0.1, 1.0), c(0.9, 10.0)))
make_dist("beta", mean = 0.4, median = 0.38)

# Mode-based
make_dist("gamma", mean = 5, mode = 3)
make_dist("beta", mean = 0.4, mode = 0.35)
make_dist("gamma", mode = 3, iqr = 4)
make_dist("logistic", mode = 5, iqr = 4)
```

## Feasibility checks

```r
dist_exists("beta", mean = 0.5, var = 0.1)        # TRUE
dist_exists("beta", mean = 0.5, var = 0.3)        # FALSE (var too large)
dist_exists("exponential", mean = 2.5, var = 6.25) # TRUE (var == mean^2)
dist_exists("exponential", mean = 2.5, var = 1.5)  # FALSE
```

## Discovery

```r
available_distributions(mean = 5, var = 3)
# [1] "gamma"    "logistic"

available_distributions(mean = 5, var = 25)
# [1] "gamma"       "exponential" "logistic"

available_distributions(mean = 0.5, var = 0.05)
# [1] "gamma"    "logistic" "beta"
```

## Installation

Not yet published. For development:

```r
# install.packages("devtools")
devtools::install_github("yoninazarathy/distsfactory-r")
```

## Authors

Ron Ashri, Sarat Moka, Yoni Nazarathy
