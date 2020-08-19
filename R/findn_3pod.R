#' Find the Sample Size Using the 3pod Algorithm
#' 
#' \code{findn_3pod} estimates the sample size for a certain target function 
#' based on a modified version of the 3pod algorithm (Wu & Tian, 2014).
#'
#' @param fun A target function that estimates the power of a trial. 
#' @param targ The target power. Must bei either 0.8 or 0.9. 
#' @param start An initial guess for the sample size.
#' @param k How often the target function is evaluated at each design point.
#' @param startevals How many evaluations are used for the first phase.
#' @param maxevals Maximum number of function evaluations.
#' @param ... Further arguments to be passed to \code{findn_3pod}. 
#'
#' @export
findn_3pod <- function(fun, targ, start, k = 50, startevals = 100, 
                       maxevals = 1000, ...) {
  func <- function(x) fun(n = x, k = k, ...)
  ttarg <- stats::qnorm(targ)
  n <- floor(maxevals / k)
  n1 <- startevals / k
  n2 <- ifelse(maxevals == 250, (maxevals * 0.48) / (2 * k), 
    (maxevals * 0.4) / (2 * k))
  n3 <- n - n1 - 2 * n2
  
  # Phase 1
  x <- round(seq(from = 0.25 * start, to = 4 * start, length.out = n1))
  y <- sapply(x, func)
  
  # Phase 2
  for (i in 1:n2) {
    ph2.mod <- fit_mod_3pod(x, y, k)
    new_x <- pmax(get_new_points_3pod(ph2.mod), 2)
    new_y <- sapply(new_x, func)
    x <- c(x, new_x)
    y <- c(y, new_y)
  }
  
  # Phase 3
  ph2.final <- fit_mod_final_3pod(x, y, k)
  b <- stats::coef(ph2.final)[2]
  b <- ifelse(b < 0.0001, 0.0001, b)
  
  x <- ty <- numeric(n3)
  x[1] <- pmin(pmax(get_final_point_3pod(ph2.final, ttarg), 2), 1000)
  ty[1] <- transf_y(func(x[1]))
  
  for(i in 2:n3) {
    x[i] <- max(x[i - 1] + (1 / (b * (i - 1))) * (ttarg - ty[i - 1]), 2)
    ty.real <- transf_y(func(round(x[i])))
    ty[i] <- ty.real + b * (x[i] - round(x[i]))
  }
  
  list(
    Point_Estimate_n = ceiling(x[length(x)]),
    All_n = x
  )
}
