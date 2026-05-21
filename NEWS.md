# distsfactory 0.2.0

## New: summary-statistic methods

Every `distsfactory_dist` now exposes `$mean()`, `$var()`, `$std()`, and `$median()`. The base-R S3 generics `mean(d)`, `median(d)`, and `quantile(d, probs)` also dispatch. This mirrors Julia's `mean(d)`, `var(d)`, `std(d)`, `median(d)` and Python's `d.mean()`, `d.var()`, `d.std()`, `d.median()`. Use `d$var()` / `d$std()` for variance and standard deviation because base R's `var()` is not an S3 generic. For Cauchy (no finite moments) these return `NaN`.

## Breaking: wrapped distributions expose `$parent` and `$support` instead of a flat `$params`

The 0.1.0 design flattened a wrapped distribution's parent parameters into its `$params` slot. For a truncated Normal that meant `$params$mean` was the *parent's* `mu`, not the truncated distribution's mean — confusing and easy to misuse.

In 0.2.0, every wrapped distribution exposes:

- `$parent`: the un-transformed `distsfactory_dist`, with its own `$params`, `$mean()`, `$var()`, etc.
- `$support`: a length-2 numeric vector giving the wrapped distribution's support.
- A transform-specific sibling slot: `$shift` (affine shift), `$flip_point` (affine flip), `$scale_loc` + `$scale_width` (affine scale), or nothing extra for truncation.

Wrapped distributions no longer have a `$params` slot — reading it returns `NULL`, which fails visibly rather than the prior silent mislabelling.

Before:

```r
d <- make_dist("normal", mean = 1.0, std = 0.8, support = c(-1, 4))
d$params   # list(mean = 0.9822, sd = 0.8232, .lo = -1, .hi = 4)  -- ambiguous
```

After:

```r
d$mean()              # 1.0           -- the truncated distribution's mean
d$support             # c(-1, 4)
d$parent$params       # list(mean = 0.9822, sd = 0.8232)  -- the parent's params
d$parent$mean()       # 0.9822        -- the parent's mean
```

To migrate: replace `d$params$<key>` with `d$parent$params$<key>` on any wrapped distribution. The `$params` slot is unchanged on plain (un-wrapped) family distributions.

## Other

- New family: `cauchy` (quantile-based construction only; no finite moments). Closes [#1](https://github.com/Distribution-Matching/distsfactory-r/issues/1). Brings the registered family count to 30.
- `make_truncated()` now takes a parent `distsfactory_dist` rather than raw `name + params + d/p/q/r`.
- `print()` on a wrapped distribution now shows the transform (e.g. `truncated to [-1, 4]`) and the parent's parameters on a second line.
- 1017 testthat assertions pass (up from 925 in 0.1.0). 56 of those are new regression tests in `tests/testthat/test-wrapper_accessors.R` covering the new slots and the print method.


# distsfactory 0.1.0

First release. R port of [DistributionsFactories.jl](https://github.com/Distribution-Matching/DistributionsFactories.jl) (Julia is the parameterization master). Cross-validated against the same Julia oracle as the Python sibling package.

## Distribution coverage

29 families across continuous (real, positive, unit) and discrete supports. Each family supports a subset of specification types — `mean+var` is universal; quantile- and mode-based forms are implemented where they exist in the Julia package.

- **Real-line continuous:** `normal`, `laplace`, `logistic`, `gumbel`, `tdist`, `uniform`, `sym_triangular`, `triangular`
- **Positive continuous:** `gamma`, `erlang`, `exponential`, `lognormal`, `weibull`, `frechet`, `chi`, `chisq`, `rayleigh`, `fdist`, `inverse_gamma`, `pareto`, `folded_normal`
- **Unit-interval continuous:** `beta`
- **Discrete:** `binomial`, `poisson`, `negative_binomial`, `geometric`, `discrete_uniform`, `discrete_sym_triangular`, `discrete_triangular`

Aliases follow the conventional R / scipy naming (`norm`, `gauss`, `lnorm`, `invgamma`, `nbinom`, `binom`, `geom`, `pois`, `dunif`, `student`, `t`, `f`, ...).

## Specification styles

`make_dist(name, ...)` accepts:

- Moments and their alternatives: `mean`, `var`, `std`, `cv`, `scv`, `second_moment`
- Quantiles: `median`, `q1`, `q3`, `iqr`, `quantiles = list(c(p1, q1), c(p2, q2))`
- Mode (alone or combined with `mean`, `var`, `iqr`, or a quantile)
- 3-parameter triangular: `mean + var + mode`

`parse_spec` validates that conflicting dispersion measures agree and rejects non-finite values up front.

## Arbitrary supports

`support = c(lo, hi)` places a distribution on a non-natural interval (either endpoint may be `Inf`). The right transform is selected automatically:

- **Affine shift / flip** when the requested shape matches the natural one (e.g. Gamma on `[3, Inf)`, Gamma on `(-Inf, 10]`)
- **Affine scale** for unit-interval families on a bounded interval (e.g. Beta on `[2, 7]`)
- **Truncation** otherwise — for Normal/Laplace/Logistic the constructor solves a 2-D Newton system on the parent `(mu, sigma)` so the truncated moments match, gating feasibility on the Langevin envelope. For positive families on bounded intervals a generic 2-D Newton with numeric quadrature is used.

The `d`, `p`, `q` methods on the resulting distribution honour `log=`, `lower.tail=`, and `log.p=`.

## Partial distributions

`partial_dist(family, ...)` pins some canonical parameters; passing the spec to `make_dist` with moment constraints solves the remaining parameters:

- 1 free parameter: 1-D bracketed brentq with family-aware bracket candidates (excluding divergent regions for `tdist`, `pareto`, `frechet`, `inverse_gamma`)
- 2 free parameters with `mean + var`: 2-D damped Newton seeded from the family's own `from_mean_var` constructor (and from `from_mean_var_mode` when the user pins the triangular mode)

## Discovery and feasibility

`dist_exists(name, ..., support=)` and `available_distributions(..., support=)` answer "is this feasible" and "which families fit" respectively. `available_distributions()` with no constraints returns the full registered set.

## Testing

902 testthat assertions in `tests/testthat/`. The suite is self-contained — `devtools::test()` does not need a Julia install. The cross-package oracle in `tests/testthat/data/cross_oracle.json` (generated by `scripts/build_cross_oracle.jl` in the Julia repo) provides numerical-parity coverage against Julia for all 65 constructor cases and all 23 feasibility cases that involve supported families.

GitHub Actions runs R CMD check on Ubuntu (release + devel) and macOS on every push.

## Known gaps

Tracked openly as GitHub issues:

- [#1](https://github.com/Distribution-Matching/distsfactory-r/issues/1) Cauchy family
- [#2](https://github.com/Distribution-Matching/distsfactory-r/issues/2) `partial_dist + support=` combination
- [#3](https://github.com/Distribution-Matching/distsfactory-r/issues/3) Discrete-support truncation (truncated Poisson, bounded Binomial subset)
- [#4](https://github.com/Distribution-Matching/distsfactory-r/issues/4) Per-family spec breadth audit vs the Python sibling
