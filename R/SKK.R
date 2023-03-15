
#' Smooth inverse Hill stat
#' Threshold selection via smooth inverse Hill statistic using \code{eBsc} nonparametric smoothing. As the latter is computationally intensive, it is recommended to limit the number of order statistics considered
#'
#' @references Schneider, L.F., Krajina, A. & Krivobokova, T. (2021). \emph{Threshold selection in univariate extreme value analysis}, Extremes, \bold{24}, 881–913 \url{https://doi.org/10.1007/s10687-021-00405-7}
#' @param data vector of positive exceedances
#' @param maxk maximum number of exceedances for the method, default to 500.
#' @return a list with components
#' \describe{
#' \item{\code{k0}}{order statistic corresponding to threshold (number of exceedances)}
#' \item{\code{tail.index}}{Hill's estimator of the tail index based on k0 exceedances}
#' \item{\code{threshold}}{numerical value of the threshold, the n-k0+1 order statistic of the original sample}
#' }
sihs <- function(data, maxk = 500L){
  # Keep log exceedances and sort
  data <- data[is.finite(data) & data > 0]
  n <- length(data)
  logdat <- sort(log(data), decreasing = TRUE)
  cumlogdat <- cumsum(logdat)
  k <- 1:(n-1)
  # Hill estimator for k=2, ..., n
  gamma_hill <- cumlogdat[k]/k - logdat[k+1]
  ihs <-  (4-k)/(2*gamma_hill*k)
  sihs <- eBsc::eBsc(ihs[1:maxk], method = "N")
  k0 <- which.min(sihs$f.hat)
  return(list(k0 = k0,
              threshold = data[k0],
              tail.index = gamma_hill[k0]))
              # sihs = sihs$f.hat,
              # ihs = ihs[1:maxk],
              # gamma_hill = gamma_hill[1:maxk])
}


#' Threshold selection via SAMSEE
#'
#' Smooth asymptotic mean squared error estimator
#' of Schneider et al. (2021) for threshold selection.
#' The implementation uses a second-order regular variation index of -1.
#'
#' @references Schneider, L.F., Krajina, A. & Krivobokova, T. (2021). \emph{Threshold selection in univariate extreme value analysis}, Extremes, \bold{24}, 881–913 \url{https://doi.org/10.1007/s10687-021-00405-7}
#' @param data vector of positive exceedances
#' @return a list with elements
#' \describe{
#' \item{\code{k}}{optimal number of exceedances}
#' #' \item{\code{gamma_hill}}{Hill's estimator of the tail index}
#' \item{\code{gamma_v}}{de Vries estimator of the tail index}
#' \item{\code{gamma_gj}}{generalized jackknife estimator of the tail index}
#' }
samsee <- function(data){
  # Keep log exceedances and sort
  data <- data[is.finite(data) & data > 0]
  n <- length(data)
  logdat <- sort(log(data), decreasing = TRUE)
  cumlogdat <- cumsum(logdat)
  k <- 1:(n-1)
  # Hill estimator for k=2, ..., n
  gamma_hill <- cumlogdat[k]/k - logdat[k+1]
  # Compute square log spacing
  Mn <- vector(mode = "numeric",
               length = length(gamma_hill))
  for(i in seq_along(Mn)){
    Mn[i] <- mean((logdat[1:k[i]]-logdat[k[i]+1])^2)
  }
  # Compute de Vries and generalized jackknife estimators
  gamma_v <- 0.5 * Mn/gamma_hill
  gamma_gj <- 2 * gamma_v - gamma_hill
  # Bias via averaging of Hill's estimator
  mean_gamma_hill <- cumsum(gamma_hill)/k

  # Minimum value for K > k
  Kmin <- 5L
  # Compute MSE of gamma_v using bias term
  AD <- vector(mode = "numeric",
               length = length(gamma_hill) - Kmin)
  # Compute bar_gamma(1, K) to bar_gamma(K-1,K)
  gamma_kK <- function(K, hill = gamma_hill){
    rev(cumsum(hill[K:1]))/(K-1:K+1)
  }
  # Compute mean squared error approx
  for(i in seq_along(AD)){
    K <- i + Kmin - 1
   AD[i] <- mean((gamma_v[1:K] + gamma_kK(K) - mean_gamma_hill[K] - gamma_hill[1:K])^2)
  }
  # Compute average 'derivative' of AD around K
  #  via finite differences (vectorized)
  i <- 3:(length(AD)-2)
  D_AD <- abs(AD[i+2]-AD[i])/2 + abs(AD[i+1]-AD[i]) + abs(AD[i-2]-AD[i])/2 + abs(AD[i-1]-AD[i])
  # Find optimal K
  Kstar <- which.min(D_AD) + Kmin
  # Compute SAMSEE
  SAMSEE <- gamma_gj[Kstar]^2/(1:Kstar) + 4 * (gamma_kK(Kstar) - mean_gamma_hill[Kstar])^2
  kstar <- which.min(SAMSEE[-c(1, length(SAMSEE))]) + 1L
  return(list(
    k0 = kstar,
    gamma_hill = gamma_hill,
    gamma_v = gamma_v,
    gamma_gj = gamma_gj
    ))
}


#' Minimum distance selection procedure
#'
#'
#' @references Clauset, A., Shalizi, C.R. and Newman, M.E.J. (2009). \emph{Power-Law Distributions in Empirical Data}. SIAM Review. Society for Industrial and Applied Mathematics, \bold{51}, 661-703, \url{https://doi.org/10.1137/070710111}
#' @param data vector of positive exceedances
#' @return a list with components
#' \describe{
#' \item{\code{k0}}{order statistic corresponding to threshold (number of exceedances)}
#' \item{\code{tail.index}}{Hill's estimator of the tail index based on k0 exceedances}
#' \item{\code{threshold}}{numerical value of the threshold, the n-k0+1 order statistic of the original sample}
#' }
mdps <- function(data){
  data <- data[is.finite(data) & data > 0]
  n <- length(data)
  data <- sort(data, decreasing = TRUE)
  logdat <- log(data)
  cumlogdat <- cumsum(logdat)
  kseq <- 2:n
  # Hill estimator for k=2, ..., n
  hill <- cumlogdat[kseq-1]/(kseq-1) - logdat[kseq]
  Dk <- numeric(length = n - 1L)
  for(k in seq_along(Dk)){
    j <- 1:(k-1)
    powS <- (data[j]/data[k])^(-1/hill[k])
   Dk[k] <- max(powS -(n-j)/(k-1), (n-j+1)/(k-1) - powS)
  }
  k0 <- which.min(Dk)
  list(k0 = k0,
       tail.index = hill[k0],
       threshold = data[k0])
}
