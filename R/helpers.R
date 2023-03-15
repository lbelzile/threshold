

#' Quantile function of folded Student-t
#'
#' @param p vector of probabilities
#' @param df degrees of freedom
#' @param lower.tail logical; if \code{TRUE}, returns the lower tail (distribution function)
#' @param log.p logical; if \code{TRUE}, argument \code{p} contains log-probabilities.
#' @keywords internal
qft <- function(p, df = 1, lower.tail = TRUE, log.p = FALSE){
  qt(0.5+p/2, df = df, lower.tail = lower.tail,  log.p =log.p)
}

#' Random generation for folded Student-t
#'
#' Generate random sample from half-Student above zero via folding
#' @param n sample size
#' @keywords internal
rft <- function(n, df){
  abs(rt(n = n, df = df))
}

# Functions for the threshold selection simulation study

#' Simulate observations from parametric models
#'
#' @param n [integer] sample size
#' @param model [list] list with arguments \code{family}, a string indicating the parametric family and a list \code{args} for the distribution parameters
#' @param rounding [character] string indicating the level of rounding of observations. Default to \code{none}
#' @param seed [integer] seed for the simulation
#' @param ar [numeric] scalar indicating the strength of the first-order autocorrelation. If missing, default to \code{0.5}
#' @examples
#' simulate_parametric(
#'   n = 100,
#'   model = list(family = "weibull", args = list(shape = 1.25)),
#'   rounding = 0.5
#' )
simulate_parametric <- function(
    n,
    model,
    rounding = 0,
    seed = NULL,
    # nobs = 1000L,
    # p = 0.5,
    ar = 0
){
  stopifnot(is.numeric(rounding),
            is.finite(rounding),
            length(rounding) == 1L,
            rounding < 1,
            rounding >= 0)
  stopifnot(is.numeric(ar),
            is.finite(ar),
            isTRUE(all(ar < 1)),
            isTRUE(length(ar) >= 1))
  serialdep <- !isTRUE(all.equal(ar, 0,
                                 check.attributes = FALSE))
  if(!is.null(seed)){
     stopifnot(length(seed) == 1L,
               isTRUE(is.finite(seed)))
    set.seed(as.integer(seed))
  }
  # Users provide a list with arguments
  # 'family': the name of the model,
  #  with quantile/rng/df functions
  stopifnot(is.list(model))
  stopifnot(!is.null(model$family))
  family <- model$family
  args <- model$args
  # Get formal arguments
  argsnames <- formalArgs(paste0("r", family))[-1]
  if(!is.null(args)){
    # Check arguments names
    stopifnot(isTRUE(all(
      names(args) %in% argsnames)))
  }
  if(!serialdep){
    dat <- do.call(paste0("r", family),
                   args = c(args, n = n))
  } else {
    # Simulate from an AR(1) model
    # Then use df / quantile transform
    # to change the marginal distribution
    nar <- length(ar)
    # This uses probability integral transform
    # to change the margins from normal
    # a poorman's attempt to create clustering
    args$p <-
     pnorm(
       arima.sim(model = list(
         order = c(nar, 0, 0),
         ar = ar),
         rand.gen = rnorm,
         n = n))
   dat <- do.call(paste0("q", family),
                  args = args)
  }
  if(rounding > 0){
  # Rounding
  # It is tempting to round quantiles on the uniform scale,
  #  but then all of the largest observations are equal
  #  and this in turn gives shape estimates of -1...
  # Instead, map the 95% to 100, then round to 0.1 or 0.5
  # and back-transform
  args$p <- 0.95
  hquant <- 100/do.call(paste0("q", family), args = args)
  round_fn <- function(x, c){
    floor(x/c)*c + 0.5*c
  }
  dat <- round_fn(hquant*dat, c = 0.1)/hquant
  }
 return(dat)
}

#' Risk measures for selected parametric models
#'
#' Given a parametric distribution, compute selected risk
#' measures. These are used to benchmark the threshold methods
#' by comparing the bias, variance, mean absolute error,
#' root mean squared error relative to an oracle that picks the
#' threshold with the best performance, or the correct model if the tail
#' is exactly generalized Pareto above a threshold
#'
#' @param model [list] a list with argument \code{family} and additional arguments within list \code{args}
#' @param p [double] probability level for quantile and expectiles
#' @param risk [string] choice of risk function, either `quantile` or `maxquant` for `penultimate`
#' @param mlen [integer] length of block maximum
#' @examples
#' # Compute return level of 1000 obs
#' risk_measures(
#'  model = list(family = "weibull",
#'               args = list(shape = 1.25)),
#'  p = 0.366,
#'  m = 1000L,
#'  risk = "maxquant")
#' risk_measures(
#'  model = list(family = "weibull",
#'               args = list(shape = 1.25)),
#'  p = 1-1/1000,
#'  risk = "quantile")
#' risk_measures(
#'  model = list(family = "weibull",
#'               args = list(shape = 1.25)),
#'  p = 0.366,
#'  mlen = 1000L,
#'  risk = "penultimate")
risk_measures <- function(
      model,
      p = 0.999,
      risk = c("quantile",
               "maxquant",
               "penultimate"),
      mlen = 1L){
  risk <- match.arg(risk)
  if(risk != "quantile"){
   mlen <- as.integer(mlen)
  stopifnot(length(mlen) == 1L,
            isTRUE(mlen >= 1L),
            is.finite(mlen))
  }
  stopifnot(is.numeric(p),
            is.finite(p),
            isTRUE(all(p < 1 & p > 0)))
  stopifnot(is.list(model))
  stopifnot(!is.null(model$family))
  family <- model$family
  args <- model$args
  # Get formal arguments
  argsnames <- formalArgs(paste0("r", family))[-1]
  if(!is.null(args)){
    # Check arguments names
    stopifnot(isTRUE(all(
      names(args) %in% argsnames)))
  }
  if(risk == "quantile"){
    do.call(what = paste0("q", family),
            args = c(args, p = p))
  } else if(risk == "maxquant"){
    do.call(what = paste0("q", family),
            args = c(args, p = p^(1/mlen)))
  } else {
    penult <-
      do.call(mev::smith.penult,
              args = c(list(family = family,
                            method = "bm",
                            m = mlen),
                       model$args))
    return(revdbayes::qgev(p = p,
              loc = penult$loc,
              scale = penult$scale,
              shape = penult$shape))
  }
}


#' Quantiles from generalized Pareto
#'
#' This function fits the generalized Pareto distribution
#' at each threshold and returns the estimated quantile of the \code{mlen}
#' \code{p}th quantile (default is median) as measure of risk.
#' The function is vectorized over \code{thresh}.
#'
#' @param data [numeric] vector of observations
#' @param thresh [numeric] vector of thresholds
#' @param p [numeric] scalar probability level for quantile
#' @param mlen [integer] scalar number of observations for the extrapolation.
compute_summary <-
  function(data,
           thresh,
           p = 0.5,
           mlen = 1000L){
    stopifnot(is.finite(thresh),
              is.numeric(p),
              is.finite(p),
              is.integer(mlen),
              length(mlen) == 1L)
  thresh <- sort(thresh)
  risk <- rep(NA, length(thresh))
  # Fit generalized Pareto
  for(i in seq_along(thresh)){
    fit <- try(mev::fit.gpd(data,
                            threshold = thresh[i]))
    if(!inherits(fit, "try-error")){
    nexc <- sum(data > thresh[i])
    # Get exponent
    pow <- 1/(mlen*nexc/length(data))
    # Compute risk measure
    risk[i] <- revdbayes::qgp(p = p^pow,
              scale = coef(fit)['scale'],
              shape = coef(fit)['shape'])
    }
  }
  return(risk + thresh)
}


