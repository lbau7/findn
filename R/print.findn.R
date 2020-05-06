get_details <- function(out, details = c("high", "low"), max_n = NULL) {
  details <- match.arg(details)
  if (!is.matrix(out$fit$vcov)) {
    if(is.null(max_n)) {
      x <- 1:(5 * out$sample_size)
    } else {
      x <- 1:max_n
    }
    pred <- pnorm(predict_fit(out$fit, x, se = FALSE))
    print.df <- data.frame("n" = x, "Est.Power" = pred)
    if(is.null(max_n)) {
      print.df <- with(print.df,
        print.df[which(Est.Power > (out$targ - 0.02) & Est.Power < (out$targ + 0.02)), ]
      )
    }
  } else {
    if(is.null(max_n)) {
      x <- 1:(5 * out$sample_size)
    } else {
      x <- 1:max_n
    }
    pred <- predict_fit(out$fit, x, se = TRUE)
    crit <- qnorm(1 - out$level / 2)
    pred.lowercl <- pnorm(pred$pred - crit * pred$se)
    pred.uppercl <- pnorm(pred$pred + crit * pred$se)

    rating <- ifelse(pred.uppercl < out$targ, "Too Low",
      ifelse(pred.lowercl > out$targ, "Sufficient", "Uncertain"))

    if(details == "low") {
      print.df <- x[min(which(rating == "Sufficient"))]
    } else {
      print.df <- data.frame(
        "n" = x,
        "Est.Power" = pnorm(pred$pred),
        "Lower CL" = pred.lowercl,
        "Upper CL" = pred.uppercl,
        "Rating" = rating
      )

      if(is.null(max_n)) {
        if (length(print.df$Rating == "Uncertain") > 0) {
          lower.unc <- min(which(print.df$Rating == "Uncertain"))
          upper.unc <- max(which(print.df$Rating == "Uncertain"))
          print.df <- print.df[(lower.unc - 3):(upper.unc + 3), ]
        } else {
          lower.suf <- min(which(print.df$Rating == "Sufficient"))
          print.df <- print.df[(lower.suf - 3):(lower.suf + 2), ]
        }
      }
    }
  }
  print.df
}



#' Printing a findn Object
#' 
#' Displays details about a sample size estimation from a \code{findn} object.
#'
#' @param x Object of class \code{findn}.
#' @param details Either \code{"low"} (default) or \code{"high"}. See also 'Details'.
#' @param max_n If \code{details = "high"} the predicted power values and confidence intervals 
#' are shown for all sample sizes from 1 to \code{max_n} if \code{max_n} is non-\code{NULL}. 
#' See also 'Details'.
#' @param digits Number of decimal places to be shown.
#'
#' @details When \code{details = "low"}, only the point estimate (i.e., the smallest sample 
#' for which the predicted power exceeds the target power), the "minimum sufficient sample 
#' size" (i.e., the smallest sample size for which the lower limit of the \code{level}% confidence
#' interval for the predicted power exceeds the target power) and an exit message. The exit 
#' message shows whether the chosen stopping rule was satisfied. If \code{details = "high"} 
#' then the default behaviour (i.e. when \code{max_n = NULL}) is to display all sample sizes, 
#' their predicted power values and the alpha% confidence intervals, for which it is uncertain
#' whether their power exceeds the target power, and the three largest sample sizes that 
#' are smaller than the smallest sample size that is rated uncertain and the three smallest
#' sample sizes which are greater than the smallest sample size that is rated uncertain.
#' If \code{details = "high"} and \code{max_n} is non-\code{NULL}, then the sample sizes,
#' their predicted power values and the confidence intervals for the predicted power values 
#' from 1 to \code{max_n} are displayed.
#'
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
#'   initevals = 100, r = 4, stop = "evals", maxevals = 2000, 
#'   level = 0.05, var_alpha = 0.05, var_beta = 1, alpha = 0.025, 
#'   alternative = "one.sided", sd = 2, verbose = FALSE)
#'
#' # print with default settings
#' print(res.ttest, details = "low", digits = 3)
print.findn <- function(x, details = c("low", "high"), max_n = NULL, digits = 3) {
  details <- match.arg(details)
  detail.df <- get_details(x, details, max_n)

  if (details == "low") {
    x_list <- list(
      Point_Estimate_n = x$sample_size,
      Minimum_Sufficient_n = format(detail.df, digits = digits),
      Message = x$exit.mes
    )
  } else {
    x_list <- list(
      Point_Estimate_n = x$sample_size,
      Details = format(detail.df, digits = digits),
      Message = x$exit.mes
    )
  }
  print(x_list)
}
