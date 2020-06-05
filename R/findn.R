#' findn
#'
#' findn estimates the sample size for a certain target function based on repeated simulations using a model
#' based approach.
#'
#' @param fun The target function that estimates the power of a trial for a sample size of 
#' \code{n} based on \code{k} simulations.
#' @param targ The target power.
#' @param start An initial guess for the sample size.
#' @param k How often the target function is evaluated at each design point. The default is 25.
#' @param initevals How many evaluations the first model is based on. The default is 100.
#' @param r A multiplicator for the range of the initial design points. The default is 4.
#' @param stop The stopping criterion. One of \code{"evals"}, \code{"power.ci"}, 
#' \code{"abs.unc"}, \code{"rel.unc"}.
#' @param maxevals Maximum number of function evaluations. The default is 2000.
#' @param level Significance level for the confidence intervals if \code{stop} is something other
#' than \code{"evals"}. Also used to determine the levels for the confidence intervals that are printed 
#' if \code{verbose = TRUE}.
#' @param power.ci.tol Tolerance parameter if \code{stop = "power.ci"}. The default is 0.02.
#' @param abs.unc.tol Tolerance parameter if \code{stop = "abs.unc"}. The default is 10.
#' @param rel.unc.tol Tolerance parameter if \code{stop is "rel.unc"}. The default is 0.1.
#' @param var_alpha Variance of the prior distribution for the intercept. The default is 0.05.
#' @param var_beta Variance of the prior distribution for the slope. The default is 1.
#' @param alpha The significance level of the underlying test. This is used to compute the mean of the prior
#' distribution of the intercept. The default is 0.05.
#' @param alternative Either "two.sided" or "one.sided". The default is "two.sided".
#' @param minx The minimum sample size that \code{fun} can be evaluated for.
#' @param verbose If \code{TRUE}, the current sample size estimate, the predicted power and its 
#' \code{level} percent confidence is returned.
#' interval is printed after every iteration.
#' @param ... Further arguments to be passed to \code{findn}.
#'
#' @details \code{findn} estimates the sample size for a target function that simulates a statistical
#' test or a trial. The function must have at least two arguments, \code{n}, the sample size for which 
#' the trial is simulated, and \code{k}, that states how often the trial is simulated. The function 
#' must return an estimate for the power of the trial for the sample size \code{n} based on \code{k}
#' Monte Carlo simulations.
#' 
#' \code{findn} uses an algorithm that assumes a probit model and computes Bayesian parameter estimates.
#' The mean for the prior distribution of the intercept is computed from the significance level 
#' \code{alpha} of the underlying test and the \code{alternative}. The mean for the prior distribution
#' of the slope is computed from the initial guess for the sample size - \code{start}. The variances
#' can be adjusted using the arguments \code{var_alpha} and \code{var_beta}.
#' 
#' There are four different stopping criterions. When \code{stop = "evals"} the algorithm stops
#' when the target function was evaluated \code{maxevals} times. When \code{stop = "power.ci"} 
#' the algorithm stops when the \code{level} percent confidence interval of the predicted power
#' at the current sample size estimate is within the interval \code{targ} plus and minus 
#' \code{power.ci.tol}. When \code{stop = "abs.unc"} the algorithm stops when the number of sample
#' sizes in the uncertainty set smaller than \code{abs.unc.tol}. The uncertainty set is defined
#' as the set that contains all sample sizes for which the \code{level} percent confidence interval 
#' for the predicted power contains \code{targ}. When \code{stop = "rel.unc"} the algorithm stops when
#' the relative uncertainty range is smaller than \code{rel.unc.tol}. The relative uncertainty range
#' is defined as the greatest integer in the uncertainty set minus the smallest integer in the 
#' uncertainty set, divided by the smallest number in the uncertainty set. The algorithm also
#' stops when \code{stop} is either \code{"power.ci"}, \code{"abs.unc"} or \code{"rel.unc"} and the
#' stopping criterion couldn't be satisfied within \code{maxevals} evaluations.
#' 
#' @return findn returns an object of class "\code{findn}". See also \code{\link{print.findn}}.
#' @export
#'
#' @examples
#' # Function that simulates the outcomes of a two-sample t-test
#' ttest <- function(mu1 = 0, mu2 = 1, sd, n, k) {
#'   sample1 <- matrix(rnorm(n = ceiling(n) * k, mean = mu1, sd = sd),
#'     ncol = k)
#'   mean1 <- apply(sample1, 2, mean)
#'   sd1_hat <- apply(sample1, 2, sd)
#'   sample2 <- matrix(rnorm(n = ceiling(n) * k, mean = mu2, sd = sd),
#'     ncol = k)
#'   mean2 <- apply(sample2, 2, mean)
#'   sd2_hat <- apply(sample2, 2, sd)
#'   sd_hat <- sqrt((sd1_hat^2 + sd2_hat^2) / 2)
#'   teststatistic <- (mean1 - mean2) / (sd_hat * sqrt(2 / n))
#'   crit <- qt(1 - 0.025, 2*n - 2)
#'   return(mean(teststatistic < -crit))
#' }
#' 
#' res.ttest <- findn(fun = ttest, targ = 0.8, k = 25, start = 100, 
#'   initevals = 100, r = 4, stop = "evals", maxevals = 2000, 
#'   level = 0.05, var_alpha = 0.05, var_beta = 1, alpha = 0.025, 
#'   alternative = "one.sided", sd = 2, verbose = FALSE)

findn <- function (fun, targ, start, k = 25, initevals = 100, r = 4,
                    stop = c("evals", "power.ci", "abs.unc", "rel.unc"),
                    maxevals = 2000, level = 0.05,
                    power.ci.tol = 0.02, abs.unc.tol = 10, rel.unc.tol = 0.1,
                    var_alpha = 0.05, var_beta = 1,
                    alpha = 0.05, alternative = c("two.sided", "one.sided"),
                    minx = 2, verbose = FALSE, ...) {
  x <- y <- xest <- numeric()
  stop <- match.arg(stop)
  alternative <- match.arg(alternative)
  alpha <- ifelse(alternative == "two.sided", alpha / 2, alpha)
  startno <- round((initevals - k) / k)
  maxiter <- floor((maxevals / k - (startno + 1)) / 2)

  if (maxiter < 1) stop("maxevals is too low")

  func <- function(x) fun(n = x, k = k, ...)
  startvals <- pmax(round(start * seq(from = 1 / r, to = r,
    length.out = startno)), minx)
  x <- pmax(startvals, minx)
  y <- sapply(x[1:startno], function(x) func(x))
  ttarg <- stats::qnorm(targ)
  ycur <- y[1:startno]
  xcur <- x[1:startno]
  init_par <- get_init_par(x = xcur, y = ycur, k = k, alpha = alpha)
  par_logbeta <- get_par_logbeta(
    n = start, alpha = alpha, targ = targ, var_beta = var_beta
    )

  fit <- fitMod(
    x = xcur, y = ycur, k = k, weights = rep(1, length(xcur)),
    start_par = init_par, mu_alpha = stats::qnorm(alpha), sd_alpha = sqrt(var_alpha),
    mu_logbeta = par_logbeta[1], sd_logbeta = par_logbeta[2]
    )

  xest[1] <- ceiling(get_est(fit = fit, ttarg = ttarg))
  x[startno + 1] <- pmax(xest[1], minx)
  y[startno + 1] <- func(x[startno + 1])
  count <- startno + 1

  for (i in 1:maxiter) {
    ind <- 1:count
    xcur <- x[ind]
    ycur <- y[ind]

    if (sum(ycur) == length(ycur)) {
      new_x <- pmax(round(c(min(xcur) / 2, min(xcur) / 4)), 2)
    } else if (sum(ycur) == 0) {
      new_x <- pmax(round(c(max(xcur) * 2, max(xcur) * 4)), 2)
    } else {
      tycurpred <- predict_fit(fit, xcur, se = FALSE)
      weights <- wgts(typred = tycurpred, ttarg = ttarg)
      fit <- fitMod(
        x = xcur, y = ycur, k = k, weights = weights, start_par = fit$cf,
        mu_alpha = stats::qnorm(alpha), sd_alpha = sqrt(var_alpha),
        mu_logbeta = par_logbeta[1], sd_logbeta = par_logbeta[2]
        )
      xest[i] <- get_est(fit, ttarg)
      
      if (verbose == TRUE) print_verbose(fit, xest[i], level)

      if (stop == "power.ci") {
        stp <- stop_power(
          tol = power.ci.tol, fit = fit, xest = xest[i], targ = targ,
          level = level
        )
        if (stp$stop) break
      }

      if (stop == "abs.unc") {
        stp <- stop_uncertainty(
          tol = abs.unc.tol, fit = fit, xest = xest[i], targ = targ,
          level = level, type = "absolute"
          )
        if (stp$stop) break
      }

      if (stop == "rel.unc") {
        stp <- stop_uncertainty(
          tol = rel.unc.tol, fit = fit, xest = xest[i], targ = targ,
          level = level, type = "relative"
        )
        if (stp$stop) break
      }

      new_x <- pmax(get_new_points(fit, xest[i], targ), minx)
    }
    new_y <- sapply(new_x, function(x) func(x))
    ind <- (count + 1):(count + 2)
    x[ind] <- new_x
    y[ind] <- c(new_y[1], new_y[2])
    count <- count + 2

    if (i == maxiter) {
      ind <- 1:count
      xcur <- x[ind]
      ycur <- y[ind]
      tycurpred <- predict_fit(fit, xcur, se = FALSE)
      weights <- wgts(typred = tycurpred, ttarg = ttarg)
      fit <- fitMod(
        x = xcur, y = ycur, k = k, weights = weights, start_par = fit$cf,
        mu_alpha = stats::qnorm(alpha), sd_alpha = sqrt(var_alpha),
        mu_logbeta = par_logbeta[1], sd_logbeta = par_logbeta[2]
      )
      xest[i+1] <- get_est(fit, ttarg)
    }
  }

  root <- ceiling(xest[length(xest)])

  if (stop != "evals" && !stp$stop ) {
    exit.mes <- "Maximum number of evaluations reached without sufficient accuracy"
  } else if (stop != "evals" && stp$stop) {
    exit.mes <- "Stopping criterion fulfilled"
  } else {
    exit.mes <- "Normal completion"
  }

  out <- list(
    sample_size = root,
    fit = fit,
    targ = targ,
    level = level,
    exit.mes = exit.mes
  )

  class(out) <- "findn"
  out
}
