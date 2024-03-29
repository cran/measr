% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loo-methods.R
\name{loo.measrfit}
\alias{loo.measrfit}
\title{Efficient approximate leave-one-out cross-validation (LOO)}
\usage{
\method{loo}{measrfit}(x, ..., r_eff = NA, force = FALSE)
}
\arguments{
\item{x}{A \link{measrfit} object.}

\item{...}{Additional arguments passed to \code{\link[loo:loo]{loo::loo.array()}}.}

\item{r_eff}{Vector of relative effective sample size estimates for the
likelihood (\code{exp(log_lik)}) of each observation. This is related to
the relative efficiency of estimating the normalizing term in
self-normalizing importance sampling when using posterior draws obtained
with MCMC. If MCMC draws are used and \code{r_eff} is not provided then
the reported PSIS effective sample sizes and Monte Carlo error estimates
will be over-optimistic. If the posterior draws are independent then
\code{r_eff=1} and can be omitted. The warning message thrown when \code{r_eff} is
not specified can be disabled by setting \code{r_eff} to \code{NA}. See the
\code{\link[loo:relative_eff]{relative_eff()}} helper functions for computing \code{r_eff}.}

\item{force}{If the LOO criterion has already been added to the model object
with \code{\link[=add_criterion]{add_criterion()}}, should it be recalculated. Default is \code{FALSE}.}
}
\value{
The object returned by \code{\link[loo:loo]{loo::loo.array()}}.
}
\description{
A \code{\link[loo:loo]{loo::loo()}} method that is customized for \code{measrfit} objects. This is a
simple wrapper around \code{\link[loo:loo]{loo::loo.array()}}. See the \strong{loo} package
\href{https://mc-stan.org/loo/articles/}{vignettes} for details.
}
