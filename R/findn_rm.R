#' Find the Sample Size Using the Robbing-Monro Algorithm
#' 
#' \code{findn_rm} estimates the sample size for a certain target function
#' using a modified version of the Robbins-Monro algorithm.
#'
#' @param fun A target function that estimates the power of a trial. 
#' @param targ The target power. Must bei either 0.8 or 0.9. 
#' @param start An initial guess for the sample size.
#' @param k How often the target function is evaluated at each design point.
#' @param maxevals Maximum number of function evaluations.
#' @param j A tuning parameter for choosing the step size. In the original
#'   Robbins-Monro algorithm j is 1.
#' @param avg A tuning parameter for the averaged veresion of the Robbins-Monro
#'   algorithm. The sample size estimate is the average of the last avg% of
#'   design points in the algorithm.
#' @param ... Further arguments to be passed to \code{findn_rm}. 
#'
#' @export
findn_rm <- function(fun, targ, start, k = 50, maxevals = 1000, j = 1,
                     avg = 0.5, ...) {
  if (j <= 0.5 | j > 1) stop("j has to be in (0.5, 1]")
  if (avg <= 0 | avg > 1) stop("avg has to be in (0, 1]")
  if (maxevals < 3 * k) stop("too little evals")
  if (maxevals %% k != 0) stop("maxevals has to be multiple of k")
  
  func <- function(x) fun(n = x, k = k, ...)
  ttarg <- stats::qnorm(targ)
  maxiter <- maxevals / k - 2
  
  x <- ty <- numeric(maxiter)
  x[2] <- start
  ty[2] <- transf_y(func(x[2]))
  x[1] <- ifelse(ty[2] < ttarg, x[2] * 4, x[2]  / 4)
  ty[1] <- transf_y(func(x[1]))
  
  b <- get_b(x = x[1:2], y = ty[1:2])
  
  for (i in 3:(maxiter + 2)) {
    x[i] <- max(x[i - 1] + (1 / (b * (i - 2)^j)) * (ttarg - ty[i - 1]), 2)
    ty.real <- transf_y(func(round(x[i])))
    ty[i] <- ty.real + b * (x[i] - round(x[i]))
    b <- get_b(x = x[1:i], y = ty[1:i])
  }
  
  list(
    Point_Estimate_n = ceiling(x[length(x)]), 
    Point_Estimate_n_avg = ceiling(mean(x[round(length(x) * avg):length(x)], 
      na.rm = TRUE)),
    All_n = x
  )
}