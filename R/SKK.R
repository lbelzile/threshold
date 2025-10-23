#' Smooth inverse Hill threshold selection
#'
#' Threshold selection via smooth inverse Hill statistic using \code{eBsc} nonparametric smoothing. As the latter is computationally intensive, it is recommended to limit the number of order statistics considered
#'
#' @references Schneider, L.F., Krajina, A. & Krivobokova, T. (2021). \emph{Threshold selection in univariate extreme value analysis}, Extremes, \bold{24}, 881–913 \url{https://doi.org/10.1007/s10687-021-00405-7}
#' @param xdat vector of positive exceedances
#' @param kmax maximum number of exceedances for the method, default to 500.
#' @return a list with components
#' \describe{
#' \item{\code{k0}}{order statistic corresponding to threshold (number of exceedances)}
#' \item{\code{shape}}{Hill's estimator of the tail index based on k0 exceedances}
#' \item{\code{thresh0}}{numerical value of the threshold, the n-k0+1 order statistic of the original sample}
#' }
#' @export
thselect.sihs <- function(
  xdat,
  kmax = 500L,
  method = c("mgcv", "ebsc")
) {
  method <- match.arg(method)
  # Keep log exceedances and sort
  xdat <- xdat[is.finite(xdat) & xdat > 0]
  n <- length(xdat)
  logdat <- sort(log(xdat), decreasing = TRUE)
  cumlogdat <- cumsum(logdat)
  k <- 1:(n - 1)
  kmax <- min(kmax, n)
  # Hill estimator for k=2, ..., n
  gamma_hill <- cumlogdat[k] / k - logdat[k + 1]
  ihs <- (4 - k) / (2 * gamma_hill * k)
  if (method == "eBsc") {
    sihs <- eBsc::eBsc(
      y = ihs[1:kmax],
      tol.lambda = 1e-6,
      method = "N"
    )
    fhat <- sihs$f.hat
  } else if (method == "mgcv") {
    library(mgcv)
    sihs <- mgcv::gamm(
      formula = y ~ s(x, bs = "cr", k = 100),
      method = "REML",
      # knots = seq(1, kmax, length.out = 50),
      data = data.frame(y = ihs[1:kmax], x = 1:kmax),
      correlation = corAR1()
    )
    plot(y = ihs[1:kmax], x = 1:kmax)
    fhat <- predict(sihs$gam)
    lines(fhat)
  }
  k0 <- which.min(as.numeric(fhat))
  return(list(k0 = k0, threshold = xdat[k0], tail.index = gamma_hill[k0]))
  # sihs = sihs$f.hat,
  # ihs = ihs[1:kmax],
  # gamma_hill = gamma_hill[1:kmax])
}
