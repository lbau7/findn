#' Find the Sample Size Using the Algorithm by Maruo et al.
#' 
#' \code{findn_maruo} estimates the sample size for a certain target function 
#' based on repeated simulations using a model based approach proposed by Maruo 
#' et al. (2018).
#'
#' @template fun 
#' @param targ The target power. Must be either 0.8 or 0.9.
#' @param start Starting value for the algorithm. Maruo et al. suggest to
#'   use 10. 
#' @template k
#' @template dotdotdot
#'
#' @references Maruo, K., Tada, K., Ishil, R. and Gosho M. (2018) An
#'   Efficient Procedure for Calculating Sample Size Through
#'   Statistical Simulations, Statistics in Biopharmaceutical Research
#'   10, 1-8.

#' @return \code{findn_maruo} returns a list containing the point estimate
#' for the sample size and a list of all sample sizes that for which the trial
#' function was evaluated.
#' @export
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
#' findn_maruo(fun = ttest, targ = 0.8)
findn_maruo <- function(fun, targ, start = 10, k = 100, ...) {
  if(targ == 0.8) {
    boundary <- 0.9
  } else if (targ == 0.9) {
    boundary <- 0.95
  } else {
    stop ("targ must be 0.8 or 0.9")
  }
  
  func <- function(x) fun(n = x, k = k, ...)
  ttarg <- stats::qnorm(targ)
  n <- start
  power <- func(n)
  
  if (power > 0.5) {
    n <- c(n, 2)
    power <- func(2)
    step <- 1
  } else if (power > 0.3) {
    step <- 2
  } else {
    step <- 5
  }

  counter <- 0
  while (counter < 2) {
    n_loop <- n[length(n)] + step
    pow_loop <- func(n_loop)
    power <- c(power, pow_loop)
    n <- c(n, n_loop)
    if (pow_loop > boundary) {
      counter <- counter + 1
    }
    if ((n_loop >= 100) & (pow_loop < 0.3)) {
      step <- 10
    }
  }
  
  n_fit <- n[which(power > 0.6)]
  pow_fit <- power[which(power > 0.6)]
  cf <- get_par_maruo(x = n_fit, y = pow_fit, k = k)
  samp_size <- as.numeric(ceiling((ttarg - cf[1])^2 / cf[2]^2))
  return(list(
    Point_Estimate_n = samp_size,
    All_n = n
  ))
}
