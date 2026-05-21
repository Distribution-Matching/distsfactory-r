# distsfactory

An R package for constructing probability distributions from partial specifications — moments, quantiles, mode, and support.

Part of the DistributionsFactories family (alongside [DistributionsFactories.jl](https://github.com/Distribution-Matching/DistributionsFactories.jl) for Julia and [distsfactory-python](https://github.com/Distribution-Matching/distsfactory-python) for Python). The Julia package is the parameterization master; this R port mirrors its behaviour and is cross-validated against it.

## Design

- **Zero heavy dependencies** — works with base R's built-in `stats` distribution functions (`dgamma`, `pnorm`, `qbeta`, …) and ships inline `d/p/q/r` for the families base R lacks (Laplace, Gumbel, Pareto, Frechet, Rayleigh, Inverse-Gamma, Chi, Folded-Normal, Sym-Triangular, Triangular, and the two discrete-triangular families).
- Specify what you know (mean, variance, quantiles, mode, support) and get back a ready-to-use distribution object.
- The returned object exposes `$d`, `$p`, `$q`, `$r` (mirroring base R's d/p/q/r convention) and `$params` for the raw parameters.

## Quick start

```r
library(distsfactory)

# Construct from moments
d <- make_dist("gamma", mean = 5, var = 3)
class(d)         # "distsfactory_dist"

# Use the distribution: d/p/q/r mirror base R's convention
d$d(2)           # density
d$p(0.95)        # CDF
d$q(0.5)         # quantile
d$r(100)         # random samples

# Summary statistics — same surface as Julia (mean(d), var(d), ...)
# and Python (d.mean(), d.var(), ...)
d$mean()         # 5
d$var()          # 3
d$std()          # 1.732
d$median()       # 4.671
mean(d); median(d); quantile(d, c(0.25, 0.75))   # S3 generics also dispatch

# Inspect the canonical R parameters (what dgamma/pgamma see)
d$params         # list(shape = 8.33, rate = 1.67)
```

### Truncated Normal on `[-1, 4]` with mean 1 and standard deviation 0.8

```r
d <- make_dist("normal", mean = 1.0, std = 0.8, support = c(-1, 4))
d$name           # "truncated_normal"

# The distribution's own moments — what you asked for.
d$mean()         # 1.0
d$std()          # 0.8

# The support and the parent are exposed as sibling slots. A wrapped
# dist has no $params of its own (the parent's parameters are not the
# wrapped distribution's moments, so flattening them would be confusing).
d$support        # c(-1, 4)
d$parent         # the un-truncated Normal (also a distsfactory_dist)
d$parent$params  # list(mean = 0.9822, sd = 0.8232)
                 # — parent Normal(mu, sigma) solved so the *truncated*
                 #   moments hit the requested (mean = 1, std = 0.8)
```

The same `$parent` accessor works on every wrapped distribution — affine shift (Gamma on `[a, Inf)`), affine flip (Gamma on `(-Inf, b]`), affine scale (Beta on `[a, b]`), and truncation. Sibling slots describe the transform: `$shift`, `$flip_point`, `$scale_loc` + `$scale_width`, or `$support` for truncation.

For real-line families (Normal / Laplace / Logistic) on a bounded interval, the constructor solves a 2-D Newton system on the parent `(mu, sigma)` so the truncated moments match, gating feasibility on the Langevin envelope. For positive-support families on a half-infinite interval the transform is an affine shift; on a bounded interval it's a generic 2-D truncation solver. Beta on an arbitrary `[a, b]` is an affine scale.

## Supported distributions

30 families across continuous (real, positive, unit) and discrete supports. Each family supports a subset of specification types — `mean+var` is universal; quantile- and mode-based forms are implemented where they exist in the Julia package.

### Continuous on `(-Inf, Inf)`

| Distribution | Free params | Methods |
|---|---|---|
| Normal | 2 | mean+var, q1+q3, mode+var |
| Student's T | 1 | mean+var (mu=0); arbitrary location-scale via `partial_dist("tdist", df=...)` (with support) |
| Cauchy | 2 | two quantiles, median+iqr (moments undefined) |
| Laplace | 2 | mean+var, two quantiles, mode+var |
| Logistic | 2 | mean+var, two quantiles, mode+iqr, mean+quantile |
| Gumbel | 2 | mean+var, two quantiles, mean+quantile, mean+mode, mode+var |
| Uniform | 2 | mean+var, two quantiles, mean+quantile |
| Symmetric Triangular | 2 | mean+var, mode+var |
| Triangular (asymmetric) | 3 | mean+var+mode |

### Continuous on `[0, Inf)`

| Distribution | Free params | Methods |
|---|---|---|
| Gamma | 2 | mean+var, mean+mode, mode+var, mode+iqr, mode+quantile, two quantiles, mean+quantile |
| Erlang | 2 | mean+var (integer shape required) |
| Exponential | 1 | mean, var, single quantile, mean+var |
| Chi-squared | 1 | mean, mean+var |
| Chi | 1 | mean, mean+var |
| Rayleigh | 1 | mean, mode, median, mean+var |
| Log-normal | 2 | mean+var, two quantiles, mean+quantile, mean+mode |
| Weibull | 2 | mean+var, two quantiles, mean+quantile |
| Frechet | 2 | mean+var, two quantiles |
| F | 2 | mean+var |
| Inverse Gamma | 2 | mean+var |
| Pareto | 2 | mean+var, two quantiles |
| Folded Normal | 2 | mean+var (2-D Newton) |

### Continuous on `[0, 1]`

| Distribution | Free params | Methods |
|---|---|---|
| Beta | 2 | mean+var, mean+mode, two quantiles, mean+quantile |

### Discrete

| Distribution | Support | Methods |
|---|---|---|
| Binomial | `{0, ..., n}` | mean+var |
| Discrete Uniform | `{a, ..., b}` | mean+var |
| Discrete Symmetric Triangular | `{mu-n, ..., mu+n}` | mean+var |
| Discrete Triangular | `{a, ..., b}` (mode at c) | mean+var+mode (approximate) |
| Poisson | `{0, 1, 2, ...}` | mean, var, mean+var |
| Negative Binomial | `{0, 1, 2, ...}` | mean+var |
| Geometric | `{0, 1, 2, ...}` | mean, var, single quantile, mean+var |

Aliases follow the conventional R / scipy naming where relevant (`norm`, `gauss`, `lnorm`, `invgamma`, `nbinom`, `binom`, `geom`, `pois`, `dunif`, `student`, `t`, `f`, ...).

## Specification styles

```r
# Moment-based
make_dist("gamma", mean = 5, var = 3)
make_dist("gamma", mean = 5, std = 2)              # std -> var
make_dist("gamma", mean = 5, cv = 0.5)             # coefficient of variation
make_dist("gamma", mean = 4, scv = 0.5)            # squared CV
make_dist("gamma", mean = 5, second_moment = 28)   # E[X^2] -> var
make_dist("exponential", mean = 3)                  # 1-parameter family

# Quantile-based
make_dist("exponential", median = 2.0)
make_dist("logistic", q1 = 2, q3 = 8)
make_dist("normal", q1 = -1, q3 = 1)
make_dist("gamma", quantiles = list(c(0.1, 1.0), c(0.9, 10.0)))
make_dist("beta", mean = 0.4, median = 0.38)
make_dist("normal", median = 5, iqr = 2)

# Mode-based
make_dist("rayleigh", mode = 2)
make_dist("gamma", mean = 5, mode = 3)
make_dist("beta", mean = 0.4, mode = 0.35)
make_dist("gamma", mode = 3, iqr = 4)
make_dist("normal", mode = 3, var = 4)
make_dist("logistic", mode = 5, iqr = 4)

# 3-parameter triangular (mean + var + mode)
make_dist("triangular", mean = 5, var = 2, mode = 4)
make_dist("discrete_triangular", mean = 5, var = 2, mode = 5)
```

## Support — affine transforms and truncation

The `support=` keyword places a distribution on an arbitrary support. The package chooses between an affine transform (when the requested support has the same shape as the natural one) and truncation (when it's strictly contained).

```r
# Affine shift — Gamma on [3, Inf)
make_dist("gamma", mean = 8, var = 3, support = c(3, Inf))

# Affine flip — Gamma on (-Inf, 10]
make_dist("gamma", mean = 5, var = 3, support = c(-Inf, 10))

# Affine scale — Beta on [2, 7]
make_dist("beta", mean = 3.5, var = 0.5, support = c(2, 7))

# Truncation — Normal on [-0.5, 0.5] with moment matching (2-D Newton)
make_dist("normal", mean = 0.1, var = 0.05, support = c(-0.5, 0.5))

# Truncation — generic positive-family on a bounded interval
make_dist("gamma", mean = 3, var = 1, support = c(0, 10))
make_dist("beta",  mean = 0.5, var = 0.02, support = c(0.2, 0.8))
```

Discrete-support truncation (e.g. truncated Poisson on `{2, ..., 10}`) is not yet implemented — see [#3](https://github.com/Distribution-Matching/distsfactory-r/issues/3). The `partial_dist + support=` combination is also a tracked gap ([#2](https://github.com/Distribution-Matching/distsfactory-r/issues/2)).

## Partial specifications — `partial_dist`

The R analog of Julia's `@dist` macro. Pin some canonical parameters and leave the rest to be solved from moment constraints.

```r
# Pin Gamma's shape, solve rate from mean
spec <- partial_dist("gamma", shape = 3.0)
d <- make_dist(spec, mean = 5.0)
d$params       # $shape = 3.0,  $rate = 0.6

# Pin Gamma's shape, solve rate from variance
d <- make_dist(spec, var = 3.0)

# Pin Normal's mean, solve sd from variance
spec <- partial_dist("normal", mean = 2.0)
make_dist(spec, var = 4.0)$params   # $mean = 2.0, $sd = 2.0

# No pins — solve both from (mean, var)
make_dist(partial_dist("gamma"), mean = 5, var = 3)
```

`partial_dist` uses the canonical (Julia-compatible) parameter set for each family. See `canonical_params(name)` for the list of names you can pin.

## Feasibility checks

```r
dist_exists("beta", mean = 0.5, var = 0.1)          # TRUE
dist_exists("beta", mean = 0.5, var = 0.3)          # FALSE (var too large)
dist_exists("exponential", mean = 2.5, var = 6.25)  # TRUE (var == mean^2)
dist_exists("exponential", mean = 2.5, var = 1.5)   # FALSE
dist_exists("tdist", mean = 1, var = 2)             # FALSE (TDist requires mean = 0)
```

When `make_dist` fails for the same input it raises an error carrying the same reason.

## Discovery

```r
available_distributions(mean = 5, var = 3)
# [1] "logistic" "gumbel"   "gamma"    "lognormal" "weibull" "frechet" ...

available_distributions(mean = 5, var = 25)
# Includes "exponential" (var = mean^2)

available_distributions(mean = 0.5, var = 0.05)
# Includes "beta"

# With support: filter to families that fit on the requested interval
available_distributions(mean = 3.5, var = 0.5, support = c(2, 7))
# [1] "beta"
```

## Testing

The test suite is self-contained — `devtools::test()` covers the package end-to-end with no Julia install required.

One file, `tests/testthat/test-cross_julia.R`, reads a checked-in JSON oracle (`tests/testthat/data/cross_oracle.json`) of reference values generated by the Julia package. It catches any cross-language numerical drift. The oracle is produced by [`scripts/build_cross_oracle.jl`](https://github.com/Distribution-Matching/DistributionsFactories.jl/blob/main/scripts/build_cross_oracle.jl) in the Julia repo, which writes the oracle into both the Python and R sibling-package test directories in a single run. Cases for families this package does not yet implement are filtered at load time. Regenerate after material changes with:

```
cd ../DistributionsFactories.jl
julia --project=. scripts/build_cross_oracle.jl
```

## Installation

Not yet on CRAN. For development:

```r
# install.packages("devtools")
devtools::install_github("Distribution-Matching/distsfactory-r")
```

## Authors

Ron Ashri, Sarat Moka, Yoni Nazarathy
