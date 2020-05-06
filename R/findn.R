logPost <- function(par, suc, fail, n, w,
                    mu_alpha, sd_alpha, mu_logbeta, sd_logbeta) {
  alpha <- par[1]
  beta <- par[2]
  p_est <- pnorm(alpha + beta * sqrt(n))
  loglik <- sum(w * dbinom(suc, suc + fail, p_est, log = TRUE))
  logprior <- dnorm(alpha, mu_alpha, sd_alpha, log = TRUE)
  logprior <- logprior + dlnorm(beta, mu_logbeta, sd_logbeta, log = TRUE)
  - loglik - logprior
}


fitMod <- function(x, y, k, weights, start_par,
                   mu_alpha, sd_alpha, mu_logbeta, sd_logbeta) {
  suc <- y * k
  fail <- (1 - y) * k

  opt <- optim(
    start_par, logPost, suc = suc, fail = fail, n = x, w = weights,
    mu_alpha = mu_alpha, sd_alpha = sd_alpha,
    mu_logbeta = mu_logbeta, sd_logbeta = sd_logbeta
    )

  cf <- opt$par
  vcov <- try(solve(optimHess(
    cf, logPost, suc=suc, fail=fail, n=x, w=weights, mu_alpha=mu_alpha,
    sd_alpha=sd_alpha, mu_logbeta=mu_logbeta, sd_logbeta=sd_logbeta
    )))

  list(cf = cf, vcov = vcov)
}

predict_fit <- function(fit, new_n, se) {
  X <- cbind(1, sqrt(new_n))
  pred <- X %*% fit$cf

  if(se == TRUE) {
    if(is.matrix(fit$vcov)) {
      se <- sqrt(diag(X %*% fit$vcov %*% t(X)))
    } else {
      se <- 0.5
    }
    data.frame(pred = pred, se = se)
  } else {
    pred
  }
}

wgts <- function (typred, ttarg) {
  dists <- abs(typred - ttarg)
  if (sum(dists < 1 & dists > 1e-04) < 3) {
    return(rep(1, length(typred)))
  }
  (1 - dists^3)^3 * (dists < 1)
}

stop_power <- function (tol, fit, xest, targ, level) {
  if (!is.matrix(fit$vcov)) {
   return(list(stop = FALSE))
  } else {
    pred <- predict_fit(fit, xest, se = TRUE)
    crit <- qnorm(1 - level / 2)
    cond1 <- pred$pred[1] + crit * pred$se[1] < qnorm(targ + tol)
    cond2 <- pred$pred[1] - crit * pred$se[1] > qnorm(targ - tol)
    return(list(stop = cond1 & cond2))
  }
}

stop_uncertainty <- function(tol, fit, xest, targ, level, type) {
  if (!is.matrix(fit$vcov)) {
    return(list(stop = FALSE))
  } else {
    x <- 1:(3*xest)
    pred <- predict_fit(fit, x, se = TRUE)
    crit <- qnorm(1 - level / 2)
    pred.lowercl <- pnorm(pred$pred - crit * pred$se)
    pred.uppercl <- pnorm(pred$pred + crit * pred$se)
    x.unc <- x[which(pred.uppercl > 0.8 & pred.lowercl < 0.8)]

    if (type == "absolute") {
      cond <- length(x.unc) < tol
    } else if (type == "relative") {
      rel.unc <- (x.unc[length(x.unc)] - x.unc[1]) / x.unc[1]
      cond <- rel.unc < tol
    }
    return(list(stop = cond))
  }
}

print_verbose <- function(fit, xest, level) {
  if(is.matrix(fit$vcov)) {
    pred <- predict_fit(fit, ceiling(xest), se = TRUE)
    crit <- qnorm(1 - level / 2)
    cat(
      "n_Estimate: ", ceiling(xest), " ",
      "Predicted Power: ", round(pnorm(pred$pred), 3), " ",
      "[", round(pnorm(pred$pred - crit * pred$se), 3), 
      "; ",
      round(pnorm(pred$pred + crit * pred$se), 3), "]",
      "\n", sep =""
    )
  } else {
    pred <- predit_fit(fit, ceiling(xest), se = FALSE)
    cat(
      "n_Estimate", ceiling(xest),
      "Predicted Power:", round(pnorm(pred$pred), 3)
    )
  }
}
        
get_init_par <- function(x, y, k, alpha) {
  ymat <- cbind(y * k, (1 - y) * k)
  off.par <- rep(qnorm(alpha), length(x))
  init.mod <- glm(
    ymat ~ -1 + sqrt(x) + offset(off.par),
    family = binomial("probit")
    )
  beta <- ifelse(coef(init.mod) < 0, 0.1, coef(init.mod))
  c(qnorm(alpha), beta)
}

get_par_logbeta <- function(n, alpha, targ, var_beta) {
  mu_beta <- (qnorm(targ) - qnorm(alpha)) / sqrt(n)
  mu_logbeta <- log(mu_beta^2 / sqrt(var_beta + mu_beta^2))
  sd_logbeta <- sqrt(log(var_beta / mu_beta^2 + 1))
  par_logbeta <- c(mu_logbeta, sd_logbeta)
  par_logbeta
}

get_est <- function(fit, ttarg) {
  cf <- fit$cf
  (as.numeric(ttarg - cf[1]) / cf[2])^2
}

get_new_points <- function(fit, xest, targ, interval) {
  pred <- predict_fit(fit, xest, se=T)
  cf <- fit$cf

  xLB <- floor(((pred$pred - pred$se - cf[1]) / cf[2])^2)
  xUB <- ceiling(((pred$pred + pred$se - cf[1]) / cf[2])^2)

  c(xLB, xUB)
}

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
  ttarg <- qnorm(targ)
  ycur <- y[1:startno]
  xcur <- x[1:startno]
  init_par <- get_init_par(x = xcur, y = ycur, k = k, alpha = alpha)
  par_logbeta <- get_par_logbeta(
    n = start, alpha = alpha, targ = targ, var_beta = var_beta
    )

  fit <- fitMod(
    x = xcur, y = ycur, k = k, weights = rep(1, length(xcur)),
    start_par = init_par, mu_alpha = qnorm(alpha), sd_alpha = sqrt(var_alpha),
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
        mu_alpha = qnorm(alpha), sd_alpha = sqrt(var_alpha),
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

      new_x <- pmax(get_new_points(fit, xest[i], targ, interval), minx)
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
        mu_alpha = qnorm(alpha), sd_alpha = sqrt(var_alpha),
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
