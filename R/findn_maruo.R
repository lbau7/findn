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
#'   10, 1-8.  \doi{https://doi.org/10.1080/19466315.2017.1349689}
#'
#' @export
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
    if ((n_loop => 100) & (pow_loop < 0.3)) {
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
