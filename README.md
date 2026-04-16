# distsfactory

An R package for constructing probability distributions from partial specifications — moments, quantiles, mode, and more.

Part of the DistributionsFactories family (alongside [DistributionsFactories.jl](https://github.com/Ron-Ash/DistributionsFactories.jl) for Julia and [distsfactory-python](https://github.com/yoninazarathy/distsfactory-python)).

## Design

- **Zero heavy dependencies** — works with base R's built-in distribution functions (`dnorm`, `pgamma`, `qbeta`, etc.)
- Specify what you know (mean, variance, quantiles, mode) and get back a distribution object
- The returned object wraps base R's `d/p/q/r` functions with parameters baked in
- Raw parameters are always accessible via `$params`

## Planned interface

```r
library(distsfactory)

# Construct from moments
d <- make_dist("gamma", mean = 5, var = 3)

# Use the built-in d/p/q/r methods (mirrors base R convention)
d$d(2)          # density:  dgamma(2, shape=..., rate=...)
d$p(0.95)       # CDF:      pgamma(0.95, shape=..., rate=...)
d$q(0.5)        # quantile: qgamma(0.5, shape=..., rate=...)
d$r(100)        # random:   rgamma(100, shape=..., rate=...)

# Access raw parameters for use with base R directly
d$params        # list(shape = 8.33, rate = 1.67)

# Various specification styles
make_dist("normal", mean = 10, std = 2)
make_dist("normal", q1 = 10, q3 = 30)
make_dist("beta", mean = 0.4, median = 0.35)
make_dist("exponential", mean = 3)
make_dist("gamma", mean = 5, cv = 0.5)

# Check if a distribution can match the given constraints
dist_exists("beta", mean = 0.5, var = 0.1)       # TRUE
dist_exists("exponential", mean = 2.5, var = 1.5) # FALSE (variance must equal mean^2)

# Discovery: which distributions fit these constraints?
available_distributions(mean = 5, var = 3)
```

## Installation

Not yet published. Development in progress.

## Authors

Ron Ashri, Sarat Moka, Yoni Nazarathy
