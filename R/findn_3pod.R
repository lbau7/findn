#' Find the Sample Size Using the 3pod Algorithm
#' 
#' \code{findn_3pod} estimates the sample size for a certain target function 
#' based on a modified version of the 3pod algorithm (Wu & Tian, 2014).
#'
#' @template fun
#' @template targ
#' @template start
#' @template k
#' @template init_evals
#' @template max_evals
#' @template dotdotdot
#'
#' @export
findn_3pod <- function(fun, targ, start, k = 50, init_evals = 100, 
                       max_evals = 1000, 
                       alpha = 0.05, alternative = c("two.sided", "one.sided"),
                       min_x = 2, max_x = 1000, ...) {
  alternative <- match.arg(alternative)
  alpha <- ifelse(alternative == "two.sided", alpha / 2, alpha)
  
  func <- function(x) fun(n = x, k = k, ...)
  ttarg <- stats::qnorm(targ)
  n <- floor(max_evals / k)
  n1 <- init_evals / k
  n2 <- ifelse(max_evals == 250, (max_evals * 0.48) / (2 * k), 
    (max_evals * 0.4) / (2 * k))
  n3 <- n - n1 - 2 * n2
  
  # Phase 1
  x <- pmax(round(seq(from = 0.25 * start, to = 4 * start, length.out = n1)),
    min_x)
  y <- sapply(x, func)
  
  # Phase 2
  for (i in 1:n2) {
    ph2.mod <- fit_mod_3pod(x = x, y = y, k = k, alpha = alpha)
    new_x <- pmax(get_new_points_3pod(ph2.mod), min_x)
    new_y <- sapply(new_x, func)
    x <- c(x, new_x)
    y <- c(y, new_y)
  }
  
  # Phase 3
  ph2.final <- fit_mod_final_3pod(x, y, k)
  b <- stats::coef(ph2.final)[2]
  b <- ifelse(b < 0.0001, 0.0001, b)
  
  x <- ty <- numeric(n3)
  x[1] <- pmin(pmax(get_final_point_3pod(ph2.final, ttarg), min_x), max_x)
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
