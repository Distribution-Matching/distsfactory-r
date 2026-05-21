## Submission summary

This is a new submission. `distsfactory` constructs probability distributions
from partial specifications (moments, quantiles, mode) by wrapping the base R
`d*`/`p*`/`q*`/`r*` functions and returning convenient distribution objects.

No package dependencies beyond base R + `stats`. `testthat` and `jsonlite` are
in `Suggests` for the test suite only.

## Test environments

* Local: macOS aarch64-apple-darwin20, R 4.5.1 (2025-06-13)
* GitHub Actions: ubuntu-latest (R release + devel), macos-latest (R release)
* win-builder: R-devel and R-release
* R-hub: linux (R-devel), macos (R-release), windows (R-devel)

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE

The single NOTE is the standard `New submission` notice.

The local R CMD check additionally reports the environmental NOTEs
`unable to verify current time`, "pandoc not installed" and "HTML Tidy
not recent enough" — these are local-system limitations that do not
apply on CRAN's check farm and are absent from the GitHub Actions / R-hub
/ win-builder runs.

## Downstream dependencies

There are no downstream dependencies on CRAN — this is a new submission.
