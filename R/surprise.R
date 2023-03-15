#' Bayesian measure of surprise
#'
#' Threshold stability plots proposed by Lee, Fan and Sisson (2015) using measures of goodness-of-fit for threshold selection.
#' Values around 1/2 imply that the generalized Pareto model is a good approximation to the threshold exceedances.
#' @details The statistics are computed in three steps for each threshold:
#' \itemize{
#' \item Sample from the posterior distribution of the scale and shape parameters
#' \item For each parameter pair drawn of the posterior, simulate a sample from the posterior predictive distribution with as many observations as the number of exceedances above the threshold
#' \item Evaluate the statistic for both the exceedances and the simulated sample
#' \item Repeat a large number of time the previous steps.
#' \item Compute Bayesian p-value as proportion of times the value exceeds that for observed data
#' }
#'
#' @param dat [vector] data
#' @param u [vector] threshold
#' @param prior [string] name of prior for the \code{revdbayes} package (default to maximal data information)
#' @param B [integer] number of posterior replications (default to 10 000)
#' @param nrep [integer] number of replications for the p-value calculation
#' @param stat [string] the statistic, either reciprocal likelihood (\code{reciplik}) or jth empirical quantile(\code{quantile})
#' @param os [integer] order statistic for the \code{quantile} statistic. Default to maximum (\code{os=1})
#' @param plot [logical] produce a plot of the Bayesian p-value against threshold
#' @param type [string] type of likelihood, either \code{partial} or \code{posterior}. See the article for references
#' @return an invisible list with elements
#' \describe{
#' \item{\code{threshold}}{\code{k} vector of candidate thresholds}
#' \item{\code{pval}}{\code{k} by \code{nrep} matrix of p-values}
#' \item{\code{mean_pval}}{average p-value for each threshold},
#' \item{\code{stat}}{string indicating the test statistic}
#' \item{\code{type}}{string indicating whether the partial or full posterior are used}
#' \item{\code{nrep}}{integer number of replication}
#' }
surprise.diag <-
  function(dat,
           u,
           type = c("posterior", "partial"),
           prior = "mdi",
           B = 1e3L,
           nrep = 100,
           stat = c("reciplik", "quantile", "rlargest"),
           os = 1,
           plot = TRUE) {
    stopifnot(
      "Data is not a vector" = is.vector(dat),
      "Threshold is not a vector" = is.vector(u),
      "Threshold is negative" = isTRUE(all(u >= 0)),
      "Threshold exceeds maximum observation" = max(dat) > max(u)
    )
    u <- sort(u)
    dat <- sort(dat)
    stat <- match.arg(stat)
    type <- match.arg(type)
    os <- as.integer(os)

    dallosgp <- function(y, scale, shape, k, log = TRUE){
      stopifnot(length(scale) == 1L,
                length(shape) == 1L,
                k > 0,
                k < length(y))
      if(is.unsorted(y)){
        y <- sort(y)
      }
      n <- length(y)
      if (y[1] < 0 || scale < 0 || (shape < 0 && shape * y[n] + scale < 0)) {
        return(-1e15)
      }
      # lgamma(n + 1L) -
      #   lgamma(n - k + 1) +
      sum(revdbayes::dgp(x= y[(n - k + 1):n],
                         loc = 0,
                         scale = scale,
                         shape = shape,
                         log = TRUE)) +
        (n - k) * revdbayes::pgp(q = y[n - k + 1L],
                                 loc = 0,
                                 scale = scale,
                                 shape = shape,
                                 log.p = TRUE)
    }

    dosgp <- function(yj, scale, shape, j, n, log = TRUE) {
      stopifnot(
        "User should only pass order statistic" = length(yj) == 1,
        length(scale) == 1L,
        length(shape) == 1L
      )
      # Out of bound constraints
      if (yj < 0 || scale < 0 || (shape < 0 && shape * yj + scale < 0)) {
        return(-1e15)
      }
      lgamma(n + 1) - lgamma(j) - lgamma(n - j + 1) + revdbayes::dgp(
        yj,
        loc = 0,
        scale = scale,
        shape = shape,
        log = TRUE
      ) + (j - 1) * log(revdbayes::pgp(
        yj,
        loc = 0,
        scale = scale,
        shape = shape
      )) + (n - j) * log(revdbayes::pgp(
        yj,
        loc = 0,
        scale = scale,
        shape = shape,
        lower.tail = FALSE
      ))
    }

    log_post_partial <- function(x, dat, yj, os) {
      if(x[1] < 0){
        return(1e-15)
      }
      sum(revdbayes::dgp(
        x = dat,
        loc = 0,
        scale = x[1],
        shape = x[2],
        log = TRUE
      )) - dosgp(
        yj = yj,
        scale = x[1],
        shape = x[2],
        j = length(dat) - os + 1,
        n = length(dat),
        log = TRUE
      ) + revdbayes::gp_mdi(x)
    }

    log_postos_partial <- function(x, dat, os) {
      # This assumes observations are ordered
      if(x[2] < 0 && x[2] * dat[length(dat)] + x[1] < 0){
        return(-1e15)
      }
      sum(revdbayes::dgp(
        x = dat,
        loc = 0,
        scale = x[1],
        shape = x[2],
        log = TRUE
      )) - dallosgp(
        y = dat,
        scale = x[1],
        shape = x[2],
        k = 10,
        log = TRUE
      ) + revdbayes::gp_mdi(x)
    }

    fp <- revdbayes::set_prior(prior = prior,
                               model = "gp")
    stat_fn1 <- function(x, scale, shape) {
      sum(revdbayes::dgp(
        x = x,
        loc = 0,
        scale = scale,
        shape = shape,
        log = TRUE
      ))
    }
    stat_fn2 <- function(x, scale, shape, os) {
      rk <- length(x) - os + 1
      yj <- sort(x)[rk]
      # lgamma(length(x) + 1) -
      #   lgamma(rk) -
      #   lgamma(length(x) - rk + 1) +
        revdbayes::dgp(
          yj,
          loc = 0,
          scale = scale,
          shape = shape,
          log = TRUE
        ) +
        (rk - 1) * log(revdbayes::pgp(
          yj,
          loc = 0,
          scale = scale,
          shape = shape
        )) +
        (length(x) - rk) *
        log(revdbayes::pgp(
          yj,
          loc = 0,
          scale = scale,
          shape = shape,
          lower.tail = FALSE
        ))
    }
    stat_fn3 <- function(x, scale, shape, os) {
      dallosgp(y = x,
               scale = scale,
               shape = shape,
               k = os,
               log = TRUE)
    }
    post_pval <- matrix(0, nrow = length(u), ncol = nrep)
    for (i in seq_along(u)) {
      exc <- dat[dat > u[i]] - u[i]
      if (stat %in% c("quantile", "rlargest") && os > length(exc)) {
        warning("Order statistic is smaller than the number of observations.")
        post_pval[i, ] <- rep(NA, nrep)
      } else {
        if (type == "posterior") {
          gpg <- revdbayes::rpost(
            n = B * nrep,
            model = "gp",
            prior = fp,
            thresh = 0,
            data = exc
          )
        } else if (type == "partial") {
          if (stat == "reciplik") {
            stop("The distribution of the test statistics is intractable.")
          }
          if(stat == "quantile"){
            ss <- rust::gpd_sum_stats(exc)
            temp <- do.call(rust::gpd_init, ss)
            min_phi <- pmax(0, temp$init_phi - 2 * temp$se_phi)
            max_phi <- pmax(0, temp$init_phi + 2 * temp$se_phi)
            phi_to_theta <- function(phi) c(phi[1], phi[2] - phi[1] / ss$xm)
            log_j <- function(x){0}
            lambda <- rust::find_lambda(logf = rust::gpd_logpost,
                                        ss = ss,
                                        d = 2,
                                        min_phi = min_phi,
                                        max_phi = max_phi,
                                        phi_to_theta = phi_to_theta,
                                        log_j = log_j)
            gpg <- rust::ru(
              logf = log_post_partial,
              n = B * nrep,
              d = 2L,
              trans = "BC",
              lambda = lambda,
              dat = exc,
              os = os
            )
          } else if(stat == "rlargest"){
            # These functions from the rust vignette
            # show how to transform the raw generalized Pareto
            # to improve sampling with a Box-Cox transformation
            ss <- rust::gpd_sum_stats(exc)
            temp <- do.call(rust::gpd_init, ss)
            min_phi <- pmax(0, temp$init_phi - 2 * temp$se_phi)
            max_phi <- pmax(0, temp$init_phi + 2 * temp$se_phi)
            phi_to_theta <- function(phi) c(phi[1], phi[2] - phi[1] / ss$xm)
            log_j <- function(x){0}
            lambda <- rust::find_lambda(logf = rust::gpd_logpost,
                                        ss = ss,
                                        d = 2,
                                        min_phi = min_phi,
                                        max_phi = max_phi,
                                        phi_to_theta = phi_to_theta,
                                        log_j = log_j)
            gpg <- rust::ru(
              logf = log_postos_partial,
              n = B * nrep,
              d = 2L,
              trans = "BC",
              lambda = lambda,
              dat = exc,
              os = os
            )
          }
        }
        pval_eval <- function(pars, stat) {
          if (stat == "reciplik") {
            stat_fn1(x = exc,
                     scale = pars[1],
                     shape = pars[2]) <=
              stat_fn1(
                x = revdbayes::rgp(
                  n = length(exc),
                  loc = 0,
                  scale = pars[1],
                  shape = pars[2]
                ),
                scale = pars[1],
                shape = pars[2]
              )
          } else if (stat == "quantile") {
            stat_fn2(
              x = exc,
              scale = pars[1],
              shape = pars[2],
              os = os
            ) <=
              stat_fn2(
                x = revdbayes::rgp(
                  n = length(exc),
                  loc = 0,
                  scale = pars[1],
                  shape = pars[2]
                ),
                scale = pars[1],
                shape = pars[2],
                os = os
              )
          } else if (stat == "rlargest") {
            stat_fn3(
              x = exc,
              scale = pars[1],
              shape = pars[2],
              os = os
            ) <=
              stat_fn3(
                x = revdbayes::rgp(
                  n = length(exc),
                  loc = 0,
                  scale = pars[1],
                  shape = pars[2]
                ),
                scale = pars[1],
                shape = pars[2],
                os = os
              )
          }
        }
        for (j in seq_len(nrep)) {
          post_pval[i, j] <-
            mean(apply(gpg$sim_vals[(j - 1L) * B + 1:B, ], 1, function(x) {
              pval_eval(pars = x, stat = stat)
            }))
        }
      }
    }
    ret <- list(
      threshold = u,
      pval = post_pval,
      mean_pval = colMeans(post_pval),
      stat = stat,
      type = type,
      nrep = nrep
    )
    class(ret) <- "mev_surprise"
    if (plot) {
      if (nrep > 10L) {
        with(
          ret,
          boxplot(
            pval ~ threshold,
            panel.first = {
              abline(a = 0.5, b = 0, col = "grey")
            },
            ylim = c(0, 1),
            yaxs = 'i',
            bty = "l",
            ylab = paste("Bayesian p-value for", switch(
              stat,
              reciplik = "reciprocal likelihood",
              quantile = "order statistic"
            ))
          )
        )
      } else{
        with(
          ret,
          plot(
            x = threshold,
            y = mean_pval,
            panel.first = {
              abline(a = 0.5, b = 0, col = "grey")
            },
            ylim = c(0, 1),
            outline = FALSE,
            yaxs = 'i',
            bty = "l",
            ylab = paste("Bayesian p-value for\n", switch(
              stat,
              reciplik = "reciprocal likelihood",
              quantile = "order statistic"
            ))
          )
        )
      }
    }
    invisible(ret)
  }
