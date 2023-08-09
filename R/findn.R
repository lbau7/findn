#' Find the Sample Size for a trial based repeated simulation using a
#' model based approach
#'
#' \code{findn} estimates the sample size to achieve a pre-defined
#' power, when the power can only be evaluated using
#' simulations. \code{findn} uses a model-based approach for this
#' purpose.
#'
#' @template fun
#' @template targ
#' @template start
#' @template k
#' @template init_evals
#' @param r A multiplicator for the range of the initial design points.
#' @param stop The stopping criterion. One of \code{"evals"}, \code{"power_ci"}, 
#'   \code{"abs_unc"}, \code{"rel_unc"}.
#' @template max_evals
#' @param level Significance level for the confidence intervals if \code{stop} 
#'   is something other than \code{"evals"}. Also used to determine the levels 
#'   for the confidence intervals that are printed if \code{verbose = TRUE}.
#' @param power_ci_tol Tolerance parameter if \code{stop = "power_ci"}.
#' @param abs_unc_tol Tolerance parameter if \code{stop = "abs_unc"}.
#' @param rel_unc_tol Tolerance parameter if \code{stop is "rel_unc"}.
#' @param var_alpha Variance of the prior distribution for the intercept.
#' @param var_beta Variance of the prior distribution for the slope.
#' @param alpha The significance level of the underlying test. This is used to 
#'   compute the mean of the prior distribution of the intercept.
#' @param alternative Either "two.sided" or "one.sided". This is only used to
#'   determine the mean of the intercept prior.
#' @param min_x The minimum sample size that \code{fun} can be evaluated for.
#' @param verbose If \code{TRUE}, the current sample size estimate, the 
#'   predicted power and its \code{level} percent confidence is returned after
#'   every iteration.
#' @template dotdotdot
#'
#' @details \code{findn} estimates the sample size for a target function that 
#' returns a simulated power value for a test or a trial. The target function 
#' must have at least two arguments, \code{n}, the sample size for which the 
#' trial is simulated, and \code{k}, that specifies how often the trial is 
#' simulated. Note that depending on how \code{fun} is written, \code{n} can 
#' either be the sample size per group or the total sample size.
#' The function has to return an estimate for the power of the trial 
#' for the sample size \code{n} based on \code{k} Monte Carlo simulations.
#' 
#' \code{findn} uses an algorithm that assumes a probit model and computes
#' Bayesian parameter estimates. The mean of the prior distribution of the
#' intercept is computed from the significance level \code{alpha} of the
#' underlying test and the \code{alternative}. The mean of the prior
#' distribution of the slope is computed from the initial guess for the sample
#' size - \code{start}. The variances of the prior distributions can be
#' adjusted using the arguments \code{var_alpha} and \code{var_beta}.
#' 
#' There are four different stopping criteria. When \code{stop = "evals"} 
#' the algorithm stops when the target function was evaluated \code{max_evals}
#' times. When \code{stop = "power_ci"}the algorithm stops when the \code{level}
#' percent confidence interval of the predicted power at the current sample 
#' size estimate is within the interval \code{targ} plus and minus
#' \code{power_ci_tol}. When \code{stop = "abs_unc"} the algorithm stops when 
#' the number of sample sizes in the uncertainty set smaller than
#' \code{abs_unc_tol}. The uncertainty set is defined as the set that contains 
#' all sample sizes for which the \code{level} percent confidence interval for
#' the predicted power contains \code{targ}. When \code{stop = "rel_unc"} the
#' algorithm stops when the relative uncertainty range is smaller than
#' \code{rel_unc_tol}. The relative uncertainty range is defined as the greatest
#' integer in the uncertainty set minus the smallest integer in the uncertainty
#' set, divided by the smallest number in the uncertainty set. The algorithm also 
#' stops when \code{stop} is either \code{"power_ci"}, \code{"abs_unc"} or
#' \code{"rel_unc"} and the stopping criterion couldn't be satisfied within
#' \code{max_evals} evaluations.
#' 
#' @return \code{findn} returns an object of class "\code{findn}". 
#'   See also \code{\link{print.findn}}.
#' @export
#'
#' @examples
#' # Function that simulates the outcomes of a two-sample t-test
#' ttest <- function(n, k, mu1 = 0, mu2 = 1, sd = 2) {
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
#'   crit <- qt(1 - 0.025, 2 * n - 2)
#'   return(mean(teststatistic < -crit))
#' }
#' 
#' res_ttest <- findn(fun = ttest, targ = 0.8, k = 25, start = 100, 
#'   init_evals = 100, r = 4, stop = "evals", max_evals = 2000, 
#'   level = 0.05, var_alpha = 1, var_beta = 0.1, alpha = 0.025, 
#'   alternative = "one.sided", verbose = FALSE)
findn <- function(fun, targ, start, k = 25, init_evals = 100, r = 4,
                  stop = c("evals", "power_ci", "abs_unc", "rel_unc"),
                  max_evals = 2000, level = 0.05,
                  power_ci_tol = 0.02, abs_unc_tol = 10, rel_unc_tol = 0.1,
                  var_alpha = 0.05, var_beta = 1,
                  alpha = 0.05, alternative = c("two.sided", "one.sided"),
                  min_x = 2, verbose = FALSE, ...) {
  x <- y <- xest <- numeric()
  stop <- match.arg(stop)
  alternative <- match.arg(alternative)
  alpha <- ifelse(alternative == "two.sided", alpha / 2, alpha)
  start_no <- round((init_evals - k) / k)
  if (start_no == 0) {
    stop("init_evals is too small, either decrease k or increase init_evals")
  }
  maxiter <- floor((max_evals / k - (start_no + 1)) / 2)
  if (maxiter < 1) {
    stop("max_evals is too small")
  }
  func <- function(x) fun(n = x, k = k, ...)
  start_vals <- pmax(round(start * seq(from = 1 / r, to = r,
    length.out = start_no)), min_x)
  x <- pmax(start_vals, min_x)
  y <- sapply(x[1:start_no], function(x) func(x))
  ttarg <- stats::qnorm(targ)
  ycur <- y[1:start_no]
  xcur <- x[1:start_no]
  init_par <- get_init_par(x = xcur, y = ycur, k = k, alpha = alpha)
  par_logbeta <- get_par_logbeta(
    n = start, alpha = alpha, targ = targ, var_beta = var_beta
    )

  fit <- fit_mod(
    x = xcur, y = ycur, k = k, weights = rep(1, length(xcur)),
    start_par = init_par, mu_alpha = stats::qnorm(alpha), 
    sd_alpha = sqrt(var_alpha), mu_logbeta = par_logbeta[1], 
    sd_logbeta = par_logbeta[2]
    )

  xest[1] <- ceiling(get_est(fit = fit, ttarg = ttarg))
  x[start_no + 1] <- pmax(xest[1], min_x)
  y[start_no + 1] <- func(x[start_no + 1])
  count <- start_no + 1

  for (i in 1:maxiter) {
    ind <- 1:count
    xcur <- x[ind]
    ycur <- y[ind]

    if (sum(ycur) == length(ycur)) {
      new_x <- pmax(round(c(min(xcur) / 2, min(xcur) / 4)), 2)
    } else if (sum(ycur) == 0) {
      # new_x <- pmax(round(c(max(xcur) * 2, max(xcur) * 4)), 2)
    } else {
      tycurpred <- predict_fit(fit, xcur, se = FALSE)
      weights <- wgts(typred = tycurpred, ttarg = ttarg)
      fit <- fit_mod(
        x = xcur, y = ycur, k = k, weights = weights, start_par = fit$cf,
        mu_alpha = stats::qnorm(alpha), sd_alpha = sqrt(var_alpha),
        mu_logbeta = par_logbeta[1], sd_logbeta = par_logbeta[2]
        )
      xest[i] <- get_est(fit, ttarg)
      
      if (verbose == TRUE) print_verbose(fit, xest[i], level)

      if (stop == "power_ci") {
        stp <- stop_power(
          tol = power_ci_tol, fit = fit, xest = xest[i], targ = targ,
          level = level
        )
        if (stp$stop) break
      }

      if (stop == "abs_unc") {
        stp <- stop_uncertainty(
          tol = abs_unc_tol, fit = fit, xest = xest[i], targ = targ,
          level = level, type = "absolute"
          )
        if (stp$stop) break
      }

      if (stop == "rel_unc") {
        stp <- stop_uncertainty(
          tol = rel_unc_tol, fit = fit, xest = xest[i], targ = targ,
          level = level, type = "relative"
        )
        if (stp$stop) break
      }

      new_x <- pmax(get_new_points(fit, xest[i], targ), min_x)
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
      fit <- fit_mod(
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
    all_evals = x,
    targ = targ,
    level = level,
    exit.mes = exit.mes
  )

  class(out) <- "findn"
  out
}
