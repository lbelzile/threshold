# Threshold selection algorithms from
# Wan and Davis (2018), Extremes
#
# Procedure:
# 1) transform observations to regularly varying scale
# 2) form pseudo-observations (radius and angles)
# 3) select sequence of candidate thresholds
# 4) compute energy test of independance
# 5) use WBS to fit median splines and detect changepoint


#' Energy statistic tests for multivariate thresholds
#'
#' @export
#' @references Wan, P. and R.A. Davis (2018), \emph{Threshold selection for multivariate heavy-tailed data}, Extremes \bold{22}(1), pp. 131-166.
#' @param dat (matrix) an \code{n} by \code{p} matrix of observations
#' @param thresh (vector) probability levels for quantiles
#' @param norm (string) one of componentwise-projection (\code{proj}), \code{min}, \code{max} or \code{lp} norm
#' @param int (integer) the index of the component if \code{norm='proj'}, or the exponent of the \code{lp} norm. \code{Inf} returns the maximum component.
#' @param B number of bootstrap replications for each thresholds
#' @param margtransfo (string) marginal distribution, either generalized Pareto above zero (\code{Pareto}) or unit Frechet (\code{Frechet})
#' @param ties.method (string) method for handling ties. See \link[base]{rank}
#' @param plot (logical) if \code{TRUE} (default), produce a plot of the p-value path
WD.diag <- function(dat,
         thresh,
         norm = c("lp","proj", "max", "min"),
         int = 1L,
         B = 999L,
         margtransfo = c("Pareto","Frechet"),
         ties.method = c("average", "first", "last", "random", "max", "min"),
         plot = TRUE
         ){

          if (!requireNamespace("breakfast", quietly = TRUE)) {
        warning("\"breakfast\" package is not installed.")
        return(invisible(NULL) )
      }
        if (!requireNamespace("energy", quietly = TRUE)) {
        warning("\"energy\" package is not installed.")
        return(invisible(NULL) )
      }
  margtransfo <- match.arg(margtransfo)
  norm <- match.arg(norm)
  ties.method = match.arg(ties.method)
  stopifnot(is.matrix(dat),
            is.integer(int) | is.infinite(int),
            int > 0,
            isTRUE(all(thresh > 0)),
            isTRUE(all(thresh < 1)),
            length(thresh) > 1L)

  # 1) transform observations to regularly varying scale
  # Create heavy-tailed pseudo observations through
  # probability integral transform
  if(margtransfo == "Pareto"){
  pseudo <- 1/(1-apply(dat, 2, rank, ties.method = ties.method)/(nrow(dat) + 1L)) - 1
  } else if(margtransfo == "Frechet"){
    pseudo <- -1/log(apply(dat, 2, rank, ties.method = ties.method)/(nrow(dat) + 1L))
  }
  if(norm == "lp" && is.infinite(int)){
    norm == "max"
  }
  # 2) form pseudo-observations (radius and angles)
  if(norm == "proj"){
    stopifnot(int <= ncol(pseudo),
              int > 0)
    radius <- pseudo[,int]
    angles <- pseudo[,-int, drop = FALSE]/radius
  } else if(norm == "lp"){
    if(int == 1L){
     radius <- rowSums(pseudo)
     angles <- pseudo[,-1, drop = FALSE]/radius
    } else{
      radius <- rowSums(pseudo^int)^(1/int)
      angles <- pseudo[,-1, drop = FALSE]/radius
    }
  } else if(norm == "min"){
    radius <- apply(pseudo, 1, min)
    angles <- as.matrix(apply(pseudo, 1, function(x){x[-which.min(x)]})/radius)
  } else if(norm == "max"){
    radius <- apply(pseudo, 1, max)
    angles <- as.matrix(apply(pseudo, 1, function(x){x[-which.max(x)]})/radius)
  }
  # Sort entries
  # od <- order(radius, decreasing = TRUE)
  # angles <- angles[od,]
  log_radius <- log(radius) #log(radius[od])
  thresh <- sort(thresh, decreasing = TRUE)
  u <- quantile(log_radius, thresh)
  pval <- numeric(length = length(u))
  # 3) compute energy test of independance energy::indep.test
  for(i in seq_along(u)){
    pval[i] <-  energy::indep.test(x = log_radius[log_radius > u[i]],
                       angles[log_radius > u[i],],
                       method = "dcov",
                       R = B)$p.value
  }
  #4) use WBS to fit median splines and detect changepoint
  pval_chpt <- breakfast::model.thresh(breakfast::sol.wbs(pval))
  # Median splines fit: pval_chpt$est
  retval <-
    data.frame(
      thresh = thresh,
      pval = pval,
      # changepoints = pval_chpt$cpts,
      median_pval = pval_chpt$est)
  if(plot){
  plot(x = 1 - thresh,
       y = rev(pval),
       pch = 20,
       bty = "l",
       ylim = c(0,1),
       yaxs = "i",
       xlab = "tail probability",
       ylab = "p-value")
  lines(x = 1 - thresh,
        y = rev(pval_chpt$est),
        col = 2)
  abline(v = rev(1 - thresh)[pval_chpt$cpts[1]])
  }
  return(invisible(retval))
}


