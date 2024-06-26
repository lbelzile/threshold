% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{simulate_parametric}
\alias{simulate_parametric}
\title{Simulate observations from parametric models}
\usage{
simulate_parametric(n, model, rounding = 0, seed = NULL, ar = 0)
}
\arguments{
\item{n}{[integer] sample size}

\item{model}{[list] list with arguments \code{family}, a string indicating the parametric family and a list \code{args} for the distribution parameters}

\item{rounding}{[character] string indicating the level of rounding of observations. Default to \code{none}}

\item{seed}{[integer] seed for the simulation}

\item{ar}{[numeric] scalar indicating the strength of the first-order autocorrelation. If missing, default to \code{0} (uncorrelated samples)}
}
\description{
Simulate observations from parametric models
}
\examples{
simulate_parametric(
  n = 100,
  model = list(family = "weibull", args = list(shape = 1.25)),
  rounding = 0.5
)
}
