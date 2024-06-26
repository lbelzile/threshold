% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WanDavis.R
\name{WD.diag}
\alias{WD.diag}
\title{Energy statistic tests for multivariate thresholds}
\usage{
WD.diag(
  dat,
  thresh,
  norm = c("lp", "proj", "max", "min"),
  int = 1L,
  B = 999L,
  margtransfo = c("Pareto", "Frechet"),
  ties.method = c("average", "first", "last", "random", "max", "min"),
  plot = TRUE
)
}
\arguments{
\item{dat}{(matrix) an \code{n} by \code{p} matrix of observations}

\item{thresh}{(vector) probability levels for quantiles}

\item{norm}{(string) one of componentwise-projection (\code{proj}), \code{min}, \code{max} or \code{lp} norm}

\item{int}{(integer) the index of the component if \code{norm='proj'}, or the exponent of the \code{lp} norm. \code{Inf} returns the maximum component.}

\item{B}{number of bootstrap replications for each thresholds}

\item{margtransfo}{(string) marginal distribution, either generalized Pareto above zero (\code{Pareto}) or unit Frechet (\code{Frechet})}

\item{ties.method}{(string) method for handling ties. See \link[base]{rank}}

\item{plot}{(logical) if \code{TRUE} (default), produce a plot of the p-value path}
}
\description{
Energy statistic tests for multivariate thresholds
}
\references{
Wan, P. and R.A. Davis (2018), \emph{Threshold selection for multivariate heavy-tailed data}, Extremes \bold{22}(1), pp. 131-166.
}
