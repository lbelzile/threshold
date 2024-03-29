% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{rrjbmixt}
\alias{rrjbmixt}
\title{Random number generation from mixture model}
\usage{
rrjbmixt(n, kappa, beta, sigma, xi, u, eps)
}
\arguments{
\item{n}{sample size}

\item{kappa}{shape of Weibull}

\item{beta}{scale of Weibull}

\item{sigma}{scale of generalized Pareto}

\item{xi}{shape of generalized Pareto}

\item{u}{threshold above which mixture}

\item{eps}{amount of overlap}
}
\description{
The samples are generated using the quantile transform, which is evaluated numerically using
root-finding methods based on the distribution function.
}
\references{
Roth, Jongbloef, Buishand (2016), \emph{Threshold selection for regional peaks-over-threshold data}, Journal of Applied Statistics
}
\keyword{internal}
