# Helper function for findn_maruo
get_par_maruo <- function(x, y, k) {
  y_mat <- cbind(y, 1 - y) * k
  mod <- stats::glm(y_mat ~ sqrt(x), family = stats::binomial(link = "probit"))
  cf <- unname(stats::coef(mod))
  return(cf)
}

# Helper functions for findn_3pod
transf_y <- function(x) {
  if (x == 0) {
    stats::qnorm(0.001) 
  } else if (x == 1) {
    stats::qnorm(0.999) 
  } else {
    stats::qnorm(x)
  }
}

fit_mod_3pod <- function (x, y, k, alpha) {
  ymat <- cbind(y, 1 - y) * k
  off.par <- rep(stats::qnorm(alpha), length(x))
  fit <- stats::glm(ymat ~ x, family = stats::binomial("probit"))
  if(stats::coef(fit)[2] < 0.0001) {
    fit <- stats::glm(ymat ~ x - 1 + offset(off.par),
      family = stats::binomial("probit"))
  }
  return(fit)
}

fit_mod_final_3pod <- function (x, y, k) {
  ymat <- cbind(y, 1 - y) * k
  fit <- stats::glm(ymat ~ x, family = stats::binomial("probit"))
  return(fit)
}

get_new_points_3pod <- function(fit, alpha) {
  cf <- stats::coef(fit)
  if (length(cf) == 1) cf <- c(stats::qnorm(alpha), cf)
  xLB <- floor((stats::qnorm(0.128) - cf[1]) / cf[2])
  xUB <- ceiling((stats::qnorm(0.872) - cf[1]) / cf[2])
  c(xLB, xUB)
}

get_final_point_3pod <- function(fit, ttarg) {
  cf <- stats::coef(fit)
  round((ttarg - cf[1]) / cf[2])
}

# Helper function for findn_rm
get_b <- function(x, y) {
  rm.fit <- stats::lm(y ~ x)
  b <- stats::coef(rm.fit)[2]
  unname(ifelse(b < 0.0001, 0.0001, b))
}

# Helper functions for findn
log_post <- function(par, suc, fail, n, w,
  mu_alpha, sd_alpha, mu_logbeta, sd_logbeta) {
  alpha <- par[1]
  beta <- par[2]
  p_est <- stats::pnorm(alpha + beta * sqrt(n))
  loglik <- sum(w * stats::dbinom(suc, suc + fail, p_est, log = TRUE))
  logprior <- stats::dnorm(alpha, mu_alpha, sd_alpha, log = TRUE)
  logprior <- logprior + stats::dlnorm(beta, mu_logbeta, sd_logbeta, log = TRUE)
  - loglik - logprior
}

fit_mod <- function(x, y, k, weights, start_par,
  mu_alpha, sd_alpha, mu_logbeta, sd_logbeta) {
  suc <- y * k
  fail <- (1 - y) * k
  
  opt <- stats::optim(
    start_par, log_post, suc = suc, fail = fail, n = x, w = weights,
    mu_alpha = mu_alpha, sd_alpha = sd_alpha,
    mu_logbeta = mu_logbeta, sd_logbeta = sd_logbeta
  )
  
  cf <- opt$par
  vcov <- try(solve(stats::optimHess(
    cf, log_post, suc=suc, fail=fail, n=x, w=weights, mu_alpha=mu_alpha,
    sd_alpha=sd_alpha, mu_logbeta=mu_logbeta, sd_logbeta=sd_logbeta
  )), silent = TRUE)
  
  if (class(vcov) == "try-error") {
    stop("vc-matrix cannot be computed.
    Try a different value for start or smaller prior variances.")
  }
  
  list(cf = cf, vcov = vcov)
}

predict_fit <- function(fit, new_n, se) {
  X <- cbind(1, sqrt(new_n))
  pred <- X %*% fit$cf
  
  if(se == TRUE) {
    se <- sqrt(diag(X %*% fit$vcov %*% t(X)))
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
  pred <- predict_fit(fit, xest, se = TRUE)
  crit <- stats::qnorm(1 - level / 2)
  cond1 <- pred$pred[1] + crit * pred$se[1] < stats::qnorm(targ + tol)
  cond2 <- pred$pred[1] - crit * pred$se[1] > stats::qnorm(targ - tol)
  return(list(stop = cond1 & cond2))
}

stop_uncertainty <- function(tol, fit, xest, targ, level, type) {
  x <- 1:(3*xest)
  pred <- predict_fit(fit, x, se = TRUE)
  crit <- stats::qnorm(1 - level / 2)
  pred.lowercl <- stats::pnorm(pred$pred - crit * pred$se)
  pred.uppercl <- stats::pnorm(pred$pred + crit * pred$se)
  x.unc <- x[which(pred.uppercl > 0.8 & pred.lowercl < 0.8)]
  
  if (type == "absolute") {
    cond <- length(x.unc) <= tol
  } else if (type == "relative") {
    rel.unc <- (x.unc[length(x.unc)] - x.unc[1]) / x.unc[1]
    cond <- rel.unc <= tol
  }
  return(list(stop = cond))
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