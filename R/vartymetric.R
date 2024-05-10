
#' Metric-based threshold selection
#'
#' Poor man's adaptation of Varty et al. for the
#' independent and identically distributed case
#' with no rounding.
#'
#' The algorithm proceeds by first computing the maximum
#' likelihood algorithm and then simulating datasets from
#' replication with parameters drawn from a bivariate normal
#' approximation to the maximum likelihood estimator distribution.
#'
#' For each bootstrap sample, we refit the
#'  model and convert the quantiles to
#' exponential or uniform variates.
#' The mean absolute or mean squared distance
#' is calculated on these. The threshold
#' returned is the one with the lowest value
#' of the metric.
#'
#' @export
#' @param xdat [numeric] vector of observations
#' @param thresh [numeric] vector of thresholds
#' @param B [integer] number of bootstrap replications
#' @param type [string] type of graph, either \code{pp} for probability-probability plots, \code{qq} for quantile-quantile plots or \code{exp} for unit exponential Q-Q plots
#' @param neval [integer] number of points at which to estimate the metric. Default to 1000
#' @param dist [string] the distance used, either absolute distance (\code{l1}) or Euclidean distance (\code{l2})
#' @param paruq [logical] if \code{TRUE}, generate bootstrap samples accounting for the sampling distribution of parameters
vmetric.diag <- function(
      xdat,
      thresh,
      B = 199L,
      type = c("exp", "qq", "pp"),
      dist = c("l1", "l2"),
      paruq = FALSE,
      neval = 1000L){
  dist <- match.arg(dist)
  type = match.arg(type)
  B <- as.integer(B)
  stopifnot(is.finite(B),
            isTRUE(B > 1),
            is.numeric(xdat))
  xdat <- xdat[is.finite(xdat)]
  thresh <- sort(thresh[is.finite(thresh)])
  nt <- length(thresh)


  # Compute metric - xdat should be exponential/gp/uniform
  compute_metric <- function(xdat,
                             ppoints,
                             type = c("exp", "qq", "pp"),
                             dist = c("l1", "l2"),
                             pars = c(1, 0)){
    dist <- match.arg(dist)
    type <- match.arg(type)
    if(type == "exp" & dist == "l1"){
      m <- mean(abs(-log(1-ppoints) - quantile(xdat, ppoints)))
    } else if(type == "exp" & dist == "l2"){
      m <- mean((-log(1-ppoints) - quantile(xdat, ppoints))^2)
    } else if(type == "qq" & dist == "l1"){
      m <- mean(abs(mev::qgp(ppoints, scale = pars[1], shape = pars[2]) - quantile(xdat, ppoints)))
    } else if(type == "qq" & dist == "l2"){
      m <- mean((mev::qgp(ppoints, scale = pars[1], shape = pars[2]) - quantile(xdat, ppoints))^2)
    } else if(type == "pp" & dist == "l1"){
      m <- mean((ppoints*(1-ppoints)/sqrt(length(xdat)))^(-1/2)*abs(ppoints - ecdf(xdat)(ppoints)))
    } else if(type == "pp" & dist == "l2"){
      m <- mean((ppoints*(1-ppoints)/sqrt(length(xdat)))^(-1/2)*(ppoints - ecdf(xdat)(ppoints))^2)
    }
    return(m)
  }
  # Container for results at each threshold
  metric <- numeric(nt)
  ps <- ppoints(n = neval)
  for(i in seq_along(thresh)){
   # 1) Fit model
   mle <- mev::fit.gpd(
     xdat = xdat,
     threshold = thresh[i])
   if(isTRUE(paruq)){
   boot_par <- mev::gpd.boot(
     object = mle,
     B = B,
     method = "post")
   }
   spars <- as.numeric(coef(mle))
   # Create GP sample, refit model
   stat <- rep(NA, B)
   for(j in seq_len(B)){
     if(isTRUE(paruq)){
      spars <-  boot_par[j,]
     }
     boot_samp <- mev::rgp(
       n = nobs(mle),
       scale = spars[1],
       shape = spars[2])
     boot_mle <- try(mev::fit.gpd(boot_samp,
                              threshold = 0))
     if(!inherits(boot_mle, "try-error")){
       if(type %in% c("exp","pp")){
         boot_samp <- revdbayes::pgp(q = boot_samp,
                 scale = coef(boot_mle)[1],
                 shape = coef(boot_mle)[2])
         if(type == "exp"){
          boot_samp <- qexp(boot_samp)
         }
       }
    stat[j] <- compute_metric(xdat = boot_samp,
                              ppoints = ps,
                              type = type,
                              dist = dist,
                              pars = coef(boot_mle))
     } else { # if the fit failed, retry again!
       j = j - 1
     }
   }
   metric[i] <- mean(stat[is.finite(stat)],
                     na.rm = TRUE)
  }
 invisible(list(thresh = thresh[which.min(metric)],
      cthresh = thresh,
      metric = metric,
      type = type,
      dist = dist))
}

#' Bootstrap approximation for generalized Pareto parameters
#'
#' Given an object of class \code{mev_gpd},
#' returns a matrix of parameter values to mimic
#' the estimation uncertainty.
#'
#' Two options are available: a normal approximation to
#' the scale and shape based on the maximum likelihood
#' estimates and the observed information matrix.
#' This method uses forward sampling to simulate
#' from a bivariate normal distribution that satisfies
#' the support and positivity constraints
#'
#' The second approximation uses the ratio-of-uniforms
#' method to obtain samples from the posterior
#' distribution with uninformative priors, thus
#' mimicking the joint distribution of maximum likelihood.
#' The benefit of the latter is that it is closer to the
#' sampling distribution in small samples and when the shape is negative.
#'
#' @param object object of class \code{mev_gpd}
#' @param B number of pairs to sample
#' @param method string; one of \code{'norm'} for the
#' normal approximation or \code{'post'} (default) for posterior sampling
#' @return a matrix of size B by 2 whose columns contain scale and shape parameters
#' @export
gpd.boot <- function(object,
                     B = 1000L,
                     method = c("post","norm")){
  method <- match.arg(method)
  B <- as.integer(B)
  stopifnot(B > 1L,
            inherits(object, "mev_gpd"))
  if(is.null(object$exceedances)){
    stop("Exported object does not contain exceedances.")
  }
  if(method == "post"){
    if (!requireNamespace("revdbayes", quietly = TRUE)) {
    stop(
      "Package \"revdbayes\" must be installed to use this function.",
      call. = FALSE
    )
    }
    rpostsamp <- suppressWarnings(
      try(revdbayes::rpost(n = B,
                     model = "gp",
                     prior = revdbayes::set_prior(
                       prior = "flat",
                       model = "gp",
                       min_xi = -1),
                     thresh = 0,
                     data = object$exceedances,
                     init_ests = coef(object),
                     trans = "BC")$sim_vals))
    if(inherits(rpostsamp, "try-error")){
      stop("Ratio-of-uniform method failed.")
    }
    boot_par <-rpostsamp$sim_vals
  } else if (method == "norm"){
    if(isTRUE(coef(object)[2] < -0.5)){
      stop("Observed information undefined:\ncannot use normal approximation")
    }
    boot_par <- matrix(NA, ncol = 2, nrow = B)
    stopifnot(isTRUE(all(eigen(object$vcov,
              only.values = TRUE)$values > 0)))
    boot_par[,2] <- rnorm(n = B,
                          mean = coef(object)[2],
                          sd = object$std.err[2])
    vmat <- vcov(object)
    cmean <- coef(object)[1] +
      vmat[1,2]/vmat[2,2]*
      (boot_par[,2] - coef(object)[2])
    csd <- sqrt(vmat[1,1] - vmat[1,2]^2/vmat[2,2])
    # This breaks down if the mean is too small,
    # below -8.3 lower bound, once standardised
    #  but the cases we consider here will have
    #  positive mean
    maxexc <- max(object$exceedances)
    # Sample one-sided truncated normal
    rltnorm <- function(n, mean, sd, lb){
      stopifnot(isTRUE(length(lb) %in% c(1L, n)),
                isTRUE(length(mean) %in% c(1L, n)),
                isTRUE(length(sd) %in% c(1L, n)))
      lbs <- (lb - mean)/sd
      mean + sd*qnorm(
        pnorm(lbs) +
          pnorm(lbs, lower.tail = FALSE)*runif(n))
    }
    if (requireNamespace("TruncatedNormal", quietly = TRUE)) {
      boot_par[,1] <-
        TruncatedNormal::rtnorm(
          n = 1,
          mu = cmean,
          sd = csd,
          lb = ifelse(boot_par[,2]< 0,
                            -boot_par[,2]*maxexc,
                            0),
          ub = rep(Inf, B))
    } else{
    # This works most of the time, but try-catch?
    boot_par[,1] <-
      rltnorm(n = B,
              mean = cmean,
              sd = csd,
              lb = ifelse(boot_par[,2]< 0,
                          -boot_par[,2]*maxexc,
                          0))
}
  colnames(boot_par) <- c("scale", "shape")
  return(boot_par)
 }
}
