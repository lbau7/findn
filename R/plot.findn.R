#' Plot of a findn Object
#'
#' @param x object of class \code{findn}.
#' @param min_n lower limit of the x-axis.
#' @param max_n upper limit of the x-axis. The default is \code{NULL}.
#' @param power_lim if \code{max_n} is \code{NULL} then the upper limit of the x-axis is the 
#' smallest sample size for which the lower limit of the \code{level} percent confidence interval
#' for the predicted power exceeds the value of \code{power_lim}. The default is 0.95.
#' @param ... Further arguments.
#'
#' @return None.
#' @importFrom rlang .data
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
#' # Create a findn object
#' res.ttest <- findn(fun = ttest, targ = 0.8, k = 25, start = 100, 
#'   init_evals = 100, r = 4, stop = "evals", max_evals = 2000, 
#'   level = 0.05, var_alpha = 0.05, var_beta = 1, alpha = 0.025, 
#'   alternative = "one.sided", sd = 2, verbose = FALSE)
#'
#' # plot with default settings
#' plot(res.ttest, power_lim = 0.95)
plot.findn <- function(x, min_n = 1, max_n = NULL, power_lim = 0.95, ...) {
  if(is.null(max_n)) {
    max_n <- 5 * get_est(x$fit, power_lim)
    data <- get_details(x, max_n = max_n)
    max_rows <- ifelse(length(which(data$Lower.CL > power_lim)) > 0,
      min(which(data$Lower.CL > power_lim)), max_n)
    data <- data[min_n:max_rows, ]
  } else {
    data <- get_details(x, max_n = max_n)
  }
  
  data$Rating <- factor(data$Rating, levels = c("Too Low", "Uncertain", "Sufficient"))
  n_unc <- data$n[which(data$Rating == "Uncertain")]
  
  if (length(n_unc) > 0) {
    low_bound <- min(n_unc)
    up_bound <- max(n_unc)
    
    ggplot2::ggplot(data, ggplot2::aes(x = .data$n)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$Est.Power)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Lower.CL, 
        ymax = .data$Upper.CL, fill = .data$Rating), 
        alpha = 0.3) +
      ggplot2::geom_hline(yintercept = x$targ, color = "darkred") +
      ggplot2::geom_vline(xintercept = low_bound, col = "darkgreen", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = up_bound, col = "darkgreen", alpha = 0.5) + 
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
      ggplot2::scale_fill_discrete(name = "Rating") +
      ggplot2::xlab("n") +
      ggplot2::ylab("Estimated Power") +
      ggplot2::theme_bw()
  } else {
    ggplot2::ggplot(data, ggplot2::aes(x = .data$n)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$Est.Power)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Lower.CL, 
        ymax = .data$Upper.CL, fill = .data$Rating), 
        alpha = 0.3) +
      ggplot2::geom_hline(yintercept = x$targ, color = "darkred") +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
      ggplot2::scale_fill_discrete(name = "Rating") +
      ggplot2::xlab("n") +
      ggplot2::ylab("Estimated Power") +
      ggplot2::theme_bw()
  }
    
}
