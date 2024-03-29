#' Find the Sample Size Using the Robbins-Monro Algorithm
#' 
#' \code{findn_rm} estimates the sample size to achieve a specified
#' power using a modified version of the Robbins-Monro algorithm.
#'
#' @template fun
#' @template targ
#' @template start
#' @template k
#' @template max_evals
#' @param j A tuning parameter for choosing the step size. In the original
#'   Robbins-Monro algorithm j is 1.
#' @param avg A tuning parameter for the averaged version of the Robbins-Monro
#'   algorithm. The sample size estimate is the average of the last avg% of
#'   design points in the algorithm.
#' @template dotdotdot
#'
#' @export
findn_rm <- function(fun, targ, start, k = 50, max_evals = 1000, j = 1,
                     avg = 0.5, ...) {
  if (j <= 0.5 | j > 1) stop("j has to be in (0.5, 1]")
  if (avg <= 0 | avg > 1) stop("avg has to be in (0, 1]")
  if (max_evals < 3 * k) stop("max_evals has to be at least 3*k")
  if (max_evals %% k != 0) stop("max_evals has to be multiple of k")
  
  func <- function(x) fun(n = x, k = k, ...)
  ttarg <- stats::qnorm(targ)
  maxiter <- max_evals / k - 2
  
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
