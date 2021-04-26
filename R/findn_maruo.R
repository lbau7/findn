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
    if ((n_loop > 100) & (pow_loop < 0.3)) {
      step <- 10
    }
  }
  
  sel_index <- min(which(power > 0.6))
  n_fit <- n[sel_index:length(n)]
  pow_fit <- power[sel_index:length(power)]
  cf <- get_par_maruo(x = n_fit, y = pow_fit, k = k)
  samp_size <- as.numeric(ceiling((ttarg - cf[1])^2 / cf[2]^2))
  return(list(
    Point_Estimate_n = samp_size,
    All_n = n
  ))
}