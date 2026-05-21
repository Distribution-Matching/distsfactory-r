#' distsfactory: Construct Probability Distributions from Partial Specifications
#'
#' @description
#' Construct probability distributions from partial specifications such as
#' moments (mean, variance), quantiles, mode, or arbitrary supports via
#' truncation and affine transforms. Wraps base R's \code{stats} distribution
#' functions and ships inline implementations for families not provided by
#' base R (Laplace, Gumbel, Pareto, Frechet, Rayleigh, InverseGamma, Chi,
#' FoldedNormal, SymTriangular, Triangular, and the two discrete-triangular
#' families).
#'
#' Part of the DistributionsFactories family — the
#' \href{https://github.com/Distribution-Matching/DistributionsFactories.jl}{Julia
#' package} is the parameterization master; this R package and the sibling
#' Python package (\code{distsfactory}) cross-check their numerical output
#' against a shared oracle.
#'
#' @keywords internal
#' @name distsfactory-package
#' @aliases distsfactory
#' @importFrom stats dbeta dbinom dcauchy dchisq dexp df dgamma dgeom dlnorm
#'   dlogis dnbinom dnorm dpois dt dunif dweibull integrate median pbeta
#'   pbinom pcauchy pchisq pexp pf pgamma pgeom plnorm plogis pnbinom pnorm
#'   ppois pt punif pweibull qbeta qbinom qcauchy qchisq qexp qf qgamma qgeom
#'   qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom
#'   rcauchy rchisq rexp rf rgamma rgeom rlnorm rlogis rnbinom rnorm rpois rt
#'   runif rweibull uniroot
"_PACKAGE"
