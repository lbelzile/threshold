#' Pickands' order statistics selection method
#'
#' Restricting to the largest fourth of the data, consider
#' the minimum Kolmogorov-Smirnov statistic, i.e., the maximum
#' absolute difference between the estimated generalized Pareto
#' and the empirical distribution of exceedances.
#'
#' @references  James Pickands III (1975). \emph{Statistical inference using extreme order statistics}, Annals of Statistics, 3(\bold{1}) 119-131, \url{https://doi.org/10.1214/aos/1176343003}
#' @param xdat [numeric] vector of observations
#' @param method [string] estimation method, either the quartiles of Pickands (1975), maximum likelihood, probability weighted moments or L-moments
#' @return a vector with elements \code{k}, the number of order statistics to keep, and \code{thresh}, the numerical value of the latter
#' @export
pickands.diag <- function(xdat, method = c("quartiles","mle","pwm","lmom")){
  method <- match.arg(method)
  xdat <- as.vector(xdat)
  xdat <- sort(xdat[is.finite(xdat)], decreasing = TRUE)
  n <- length(xdat)
  mmax <- floor(n/4)
  mmin <- 10L
  if(mmin > mmax){
    stop("Not enough observations for threshold selection method.")
  }
  m_candidate <- mmin:mmax
  shape <- scale <- dist <- numeric(length(m_candidate))
  for(i in seq_along(m_candidate)){
    m <- m_candidate[i]
    samp <- xdat[seq_len(m-1)] - xdat[m]
    if(method == "quartiles"){
      quants <- as.numeric(quantile(samp, probs = c(0.5,0.75)))
      shape[i] <- (log(diff(quants)) - log(quants[1]))/log(2)
      scale[i] <- quants[1]*shape[i]/(2^shape[i]-1)
    } else if(method == "mle"){
     coefs <- coef(mev::fit.gpd(xdat = samp, threshold = 0))
     shape[i] <- coefs['shape']
     scale[i] <- coefs['scale']
    }  else if(method == "pwm"){
      pars <- gp.pwm(samp, sorted = TRUE)
      scale[i] <- pars['scale']
      shape[i] <- pars['shape']
    } else if(method == "lmom"){
      pars <- gp.lmom(samp, sorted = TRUE)
      scale[i] <- pars['scale']
      shape[i] <- pars['shape']
    }
    dist[i] <- max(abs(rank(samp)/n - mev::pgp(samp, scale = scale[i], shape = shape[i])))
  }
  k = as.integer(m_candidate[which.min(dist)])
  return(c(k = k, thresh = xdat[k]))
}

#' Automatic L-moment ratio selection method
#'
#' Given a sample of observations, calculate the L-skewness and L-kurtosis
#' over a set of candidate thresholds. For each threshold candidate, we
#' find the L-skewness that minimizes the sum of squared distance between the
#' theoretical L-skewness and L-kurtosis of the generalized Pareto distribution,
#' \deqn{\min_{\tau_3} (t_3-\tau_3)^2 + [t_4 - \tau_3(1+5\tau_3)/(5+\tau_3)]^2.}
#' The function returns the threshold with the minimum distance.
#'
#' @param xdat [numeric] vector of observations
#' @param thresh [numeric] vector of candidate thresholds. If missing, 20 sample quantiles starting at the 0.25 quantile in increments of 3.75 percent.
#' @param plot [logical] if \code{TRUE}, return a plot of the sample L-kurtosis against the L-skewness, along with the theoretical generalized Pareto curve.
#' @return a numeric threshold
#' Silva Lomba, J., Fraga Alves, M.I. (2020). \emph{L-moments for automatic threshold selection in extreme value analysis}. Stoch Environ Res Risk Assess, 34, 465–491. \url{https://doi.org/10.1007/s00477-020-01789-x}
#' @export
alrsm.diag <- function(xdat, thresh, plot = FALSE){
  xdat <- as.numeric(xdat[is.finite(xdat)])
  xdat <- sort(xdat)
 if(missing(thresh)){
  thresh <- quantile(xdat, probs = seq(0.25, length.out = 20, by = 0.0375))
 }
  lkurtosis_gpd <- function(tau3) {
    tau3 * (1 + 5 * tau3)/(5 + tau3)
  }
  dist <- numeric(length(thresh))
  lratios <- matrix(nrow = length(thresh), ncol = 2)
  for(i in seq_along(dist)){
    lmom <- lmoments(xdat[xdat > thresh[i]] - thresh[i], sorted = TRUE)
    t3 <- lmom[4]/lmom[3]
    t4 <- lmom[4]/lmom[2]
    lratios[i,] <- c(t3, t4)
    obj <- function(x){ (t3 - x)^2 + (t4 - lkurtosis_gpd(x))^2 }
    optim <- nlm(f = obj, p = t3)
    dist[i] <- optim$minimum
  }
  if(isTRUE(plot)){
    # Theoretical L-kurtosis for the GP based on L-skewness

    plot(lkurtosis_gpd,
         xlab = expression(tau[3]),
         ylab = expression(tau[4]),
         col = "grey", lwd = 2, bty = "l",
         xlim = c(0, 1))
    mtext(text = "L-moments", adj = 0)
    lines(x = lratios[,1],
          y = lratios[,2],
          type = "b",
          col = ifelse(dist == min(dist), "black","grey"))
  }
  return(as.numeric(thresh[which.min(dist)]))
}

# Probability weighted moments

pwm <- function(xdat, sorted = FALSE){
  if(!isTRUE(sorted)){
    xdat <- sort(xdat)
  }
  n <- length(xdat)
  lg <- lgamma(seq_len(n + 1L))
  # w <- lchoose(n-seq_len(n),r)
  as <- numeric(4)
  for(r in 0:3){
    as[r+1L] <- sum(xdat * exp(c(lg[n+1-seq_len(n-r)] - lg[n+1-seq_len(n-r)-r] - log(n), rep(-Inf, r)))) * exp(lg[n-r] - lg[n])
  }
  return(as)
}

# Compute sample l-moments (unbiased estimator)

#lmom::samlmu(xdat, ratios = FALSE) # faster because in Fortran
lmoments <- function(xdat, sorted = FALSE){
  as <- pwm(xdat, sorted = sorted)
  c(as[1], as[1]-2*as[2], as[1]-6*as[2]+6*as[3], as[1]-12*as[2]+30*as[3]-20*as[4])
}

#' Probability weighted moment estimator
#'
#' Given a vector of sample exceedances, returns a parameter estimates for the generalized Pareto distribution
#' based on the two first empirical probability weighted moments.
#' @inheritParams gp.lmom
#' @export
#' @return a vector of length two with the scale and shape estimates
gp.pwm <- function(xdat, sorted = FALSE){
  as <- pwm(xdat, sorted = sorted)
  c(scale = 2*as[1]*as[2]/(as[1]-2*as[2]), shape = as[1]/(2*as[2]-as[1]) + 2)
}

#' Estimation of generalized Pareto parameters via L-moments
#'
#' Given a sample of exceedances, compute the first four L-moments and
#' use L-skewness and L-scale to compute the scale and shape of the generalized Pareto distribution
#'
#' @param xdat [numeric] vector of observations
#' @param sorted [logical] if \code{TRUE}, observations are sorted in increasing order
#' @return a vector of length two with the scale and shape estimates
gp.lmom <- function(xdat, sorted = FALSE){
 lmom <- lmoments(xdat, sorted = sorted)
 t3 <- lmom[3]/lmom[2]
  xi <- (1-3*t3)/(1+t3)
  sigma <- (1 - xi)*(2 -xi)*lmom[2]
 return(c(scale = sigma, shape = xi))
}


#' Mahalanobis distance-based methodology
#'
#' Compute the Mahalanobis distance-based threshold method over a grid of thresholds by
#' transforming data from generalized Pareto to unit exponential based on probability weighted moment estimates,
#' then computing the first L-moment and the L-skewness. The latter are compared to the
#' theoretical counterparts from a unit exponential sample of the same size, which is used to compute
#' the Mahalanobis distance. The threshold returned is the one which minimizes the distance.
#'
#' @inheritParams alrsm.diag
#' @param approx [string] method to use to obtain moments of first L-moment
#' @param B [integer] number of replications for Monte Carlo approximation
#' @return a vector of length 2 with the selected threshold and a corresponding \emph{p}-value based on a chi-square approximation to the test statistic
#' @references Kiran, K. G. and Srivinas, V.V. (2021). \emph{A Mahalanobis distance-based automatic threshold selection method for peaks over threshold model.} Water Resources Research 57. \url{https://doi.org/10.1029/2020WR027534}
#' @export
mahadist.diag <- function(xdat, thresh, approx = c("asymptotic","mc"),  B = 1e3L){
  approx <- match.arg(approx)
  xdat <- as.numeric(xdat[is.finite(xdat)])
  xdat <- sort(xdat)
  if(missing(thresh)){
    thresh <- quantile(xdat, probs = seq(0.7, 1-20/length(xdat), length.out = 20))
  }
  dist <- numeric(length = length(thresh))
  for(i in seq_along(thresh)){
    samp <- xdat[xdat > thresh[i]] - thresh[i]
    nexc <- length(samp)
    pars <- gp.pwm(xdat = samp, sorted = TRUE)
    zl <- lmoments(-log(mev::pgp(q = samp, scale = pars['scale'], shape = pars['shape'], lower.tail = TRUE)))
    xp <- c(zl[1], zl[3]/zl[2])

    if(approx == "mc"){
     results <- matrix(nrow = B, ncol = 2L)
     for(j in seq_len(B)){
      lmom <- lmoments(rexp(nexc))
      results[j,] <- c(lmom[1], lmom[3]/lmom[2])
     }
     mus <- colMeans(results)
     covs <- cov(results)
    } else if(approx == "asymptotic"){
     mus <- c(1, 1/3-0.307*nexc^(-1.076))
     covs <- diag(c(1/sqrt(nexc), (166.507*nexc + 102.384)/(347.652*nexc^(1.498)-43.041*nexc-575.098)))
    }
    dist[i] <- mahalanobis(x = xp, center = mus, cov = covs)
  }
  return(c(thresh = as.numeric(thresh[which.min(dist)]), pval = qchisq(min(dist), df = 2)))
}
