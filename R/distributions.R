# NOT USED BECAUSE BIMODAL.

#' Gumbel-generalized Pareto mixture with contamination
#'
#' The model follows a Weibull distribution below \code{u} and a generalized Pareto above \code{u}, but blends continuously on the interval [\code{u}, \code{u}+\code{eps}]
#'
#' @param kappa shape of Weibull
#' @param beta scale of Weibull
#' @param sigma scale of generalized Pareto
#' @param xi shape of generalized Pareto
#' @param u threshold above which mixture
#' @param eps amount of overlap
#' @param lower.tail logical; if \code{TRUE}, return distribution function and survival function otherwise
#' @keywords internal
#' @references Roth, Jongbloef, Buishand (2016), \emph{Threshold selection for regional peaks-over-threshold data}, Journal of Applied Statistics
#' @examples
#' q <- seq(0, 16, length.out = 200)
#' plot(x = q,
#'      y = prjbmixt(q = q,
#'           kappa = 0.7,
#'           beta = 2,
#'           sigma = 2,
#'           xi = (1-0.7)/0.7*((9/2)^(-0.7)),
#'           u = 9,
#'           eps = 1),
#'     xlab = "quantile",
#'     ylab = "distribution function",
#'     main = "Weibull-GP mixture",
#'     type = "l")
prjbmixt <- function(q, kappa, beta,
                     sigma, xi, u, eps,
                     lower.tail = TRUE){
  # Compute cumulative hazard and return survival function
  H1 <- function(x){
    beta^(-kappa)*x^kappa
  }
  H2 <- function(x){
    1/xi*log1p(pmax(0,xi*(x-u)/sigma))
  }
  # Numerical integration using Sage
  Ht <- function(x){
    kappa*(beta*(x/beta)^kappa/kappa - 3*beta*u^2*(x/beta)^kappa/(eps^2*kappa) -
             2*beta*u^3*(x/beta)^kappa/(eps^3*kappa) +
             6*beta^(-kappa + 1)*u*x^(kappa + 1)/(eps^2*(kappa + 1)) +
             6*beta^(-kappa + 1)*u^2*x^(kappa + 1)/(eps^3*(kappa + 1)) -
             3*beta^(-kappa + 1)*x^(kappa + 2)/(eps^2*(kappa + 2)) -
             6*beta^(-kappa + 1)*u*x^(kappa + 2)/(eps^3*(kappa + 2)) +
             2*beta^(-kappa + 1)*x^(kappa + 3)/(eps^3*(kappa + 3)))/beta +
      1/6*(3*(3*eps + 4*u)*x^2 - 4*x^3- 6*(3*eps*u + 2*u^2)*x)/(eps^3*sigma*xi)
  }
  H1_u <- H1(u)
  Ht_u <- Ht(u)
  Ht_upe <- Ht(u + eps)
  survvec <- ifelse(q <= u, exp(-H1(q)),
                    ifelse(u < q & q < u+eps,
                           exp(-(H1_u + Ht(q) - Ht_u)),
                           exp(-(H1_u + Ht_upe - Ht_u + H2(q)))))
  if(lower.tail){
    return(1-survvec)
  } else{
    return(survvec)
  }
}

#' Random number generation from mixture model
#'
#' The samples are generated using the quantile transform, which is evaluated numerically using
#' root-finding methods based on the distribution function.
#' @inheritParams prjbmixt
#' @param n sample size
#' @references Roth, Jongbloef, Buishand (2016), \emph{Threshold selection for regional peaks-over-threshold data}, Journal of Applied Statistics
#' @keywords internal
rrjbmixt <- function(n, kappa, beta,
                     sigma, xi, u, eps){
  sapply(runif(n), function(ui){
    uniroot(f = function(x){
      prjbmixt(x, kappa = kappa,
               beta = beta,
               sigma = sigma,
               xi = xi,
               u = u,
               eps = eps) - ui},
      lower = 0,
      upper = 1e10)$root})
}
