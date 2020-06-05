logPost <- function(par, suc, fail, n, w,
  mu_alpha, sd_alpha, mu_logbeta, sd_logbeta) {
  alpha <- par[1]
  beta <- par[2]
  p_est <- stats::pnorm(alpha + beta * sqrt(n))
  loglik <- sum(w * stats::dbinom(suc, suc + fail, p_est, log = TRUE))
  logprior <- stats::dnorm(alpha, mu_alpha, sd_alpha, log = TRUE)
  logprior <- logprior + stats::dlnorm(beta, mu_logbeta, sd_logbeta, log = TRUE)
  - loglik - logprior
}


fitMod <- function(x, y, k, weights, start_par,
  mu_alpha, sd_alpha, mu_logbeta, sd_logbeta) {
  suc <- y * k
  fail <- (1 - y) * k
  
  opt <- stats::optim(
    start_par, logPost, suc = suc, fail = fail, n = x, w = weights,
    mu_alpha = mu_alpha, sd_alpha = sd_alpha,
    mu_logbeta = mu_logbeta, sd_logbeta = sd_logbeta
  )
  
  cf <- opt$par
  vcov <- try(solve(stats::optimHess(
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
    crit <- stats::qnorm(1 - level / 2)
    cond1 <- pred$pred[1] + crit * pred$se[1] < stats::qnorm(targ + tol)
    cond2 <- pred$pred[1] - crit * pred$se[1] > stats::qnorm(targ - tol)
    return(list(stop = cond1 & cond2))
  }
}

stop_uncertainty <- function(tol, fit, xest, targ, level, type) {
  if (!is.matrix(fit$vcov)) {
    return(list(stop = FALSE))
  } else {
    x <- 1:(3*xest)
    pred <- predict_fit(fit, x, se = TRUE)
    crit <- stats::qnorm(1 - level / 2)
    pred.lowercl <- stats::pnorm(pred$pred - crit * pred$se)
    pred.uppercl <- stats::pnorm(pred$pred + crit * pred$se)
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
    crit <- stats::qnorm(1 - level / 2)
    cat(
      "n_Estimate: ", ceiling(xest), " ",
      "Predicted Power: ", round(stats::pnorm(pred$pred), 3), " ",
      "[", round(stats::pnorm(pred$pred - crit * pred$se), 3), 
      "; ",
      round(stats::pnorm(pred$pred + crit * pred$se), 3), "]",
      "\n", sep =""
    )
  } else {
    pred <- predict_fit(fit, ceiling(xest), se = FALSE)
    cat(
      "n_Estimate", ceiling(xest),
      "Predicted Power:", round(stats::pnorm(pred$pred), 3)
    )
  }
}

get_init_par <- function(x, y, k, alpha) {
  ymat <- cbind(y * k, (1 - y) * k)
  off.par <- rep(stats::qnorm(alpha), length(x))
  init.mod <- stats::glm(
    ymat ~ -1 + sqrt(x) + offset(off.par),
    family = stats::binomial("probit")
  )
  beta <- ifelse(stats::coef(init.mod) < 0, 0.1, stats::coef(init.mod))
  c(stats::qnorm(alpha), beta)
}

get_par_logbeta <- function(n, alpha, targ, var_beta) {
  mu_beta <- (stats::qnorm(targ) - stats::qnorm(alpha)) / sqrt(n)
  mu_logbeta <- log(mu_beta^2 / sqrt(var_beta + mu_beta^2))
  sd_logbeta <- sqrt(log(var_beta / mu_beta^2 + 1))
  par_logbeta <- c(mu_logbeta, sd_logbeta)
  par_logbeta
}

get_est <- function(fit, ttarg) {
  cf <- fit$cf
  (as.numeric(ttarg - cf[1]) / cf[2])^2
}

get_new_points <- function(fit, xest, targ) {
  pred <- predict_fit(fit, xest, se=T)
  cf <- fit$cf
  
  xLB <- floor(((pred$pred - pred$se - cf[1]) / cf[2])^2)
  xUB <- ceiling(((pred$pred + pred$se - cf[1]) / cf[2])^2)
  
  c(xLB, xUB)
}

get_details <- function(out, details = c("high", "low"), max_n = NULL) {
  details <- match.arg(details)
  if (!is.matrix(out$fit$vcov)) {
    if(is.null(max_n)) {
      x <- 1:(5 * out$sample_size)
    } else {
      x <- 1:max_n
    }
    pred <- stats::pnorm(predict_fit(out$fit, x, se = FALSE))
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
    crit <- stats::qnorm(1 - out$level / 2)
    pred.lowercl <- stats::pnorm(pred$pred - crit * pred$se)
    pred.uppercl <- stats::pnorm(pred$pred + crit * pred$se)
    
    rating <- ifelse(pred.uppercl < out$targ, "Too Low",
      ifelse(pred.lowercl > out$targ, "Sufficient", "Uncertain"))
    
    if(details == "low") {
      print.df <- x[min(which(rating == "Sufficient"))]
    } else {
      print.df <- data.frame(
        "n" = x,
        "Est.Power" = stats::pnorm(pred$pred),
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