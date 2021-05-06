test_that("get_b works", {
  set.seed(20210330)
  x <- stats::rnorm(100, sd = 10)
  y1 <- 10 + 3*x + stats::rnorm(100, sd = 0.1)
  y2 <- 10 + stats::rnorm(100, sd = 0.01)
  mod1 <- stats::lm(y1 ~ x)
  
  b1 <- get_b(x = x, y = y1)
  b2 <- get_b(x = x, y = y2)
  
  expect_equal(b1, unname(stats::coef(mod1)[2]))
  expect_equal(b2, 1e-4)
})

test_that("get_par_maruo works", {
  set.seed(20210330)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  x <- 5:30
  y <- sapply(x, function(x) fun_ttest(n = x, k = 50))
  y_mat <- cbind(y * 50, (1 - y) * 50)
  mod <- glm(y_mat ~ sqrt(x), family = stats::binomial(link = "probit"))
  cf <- unname(stats::coef(mod))
  par <- get_par_maruo(x = x, y = y, k = 50)
  
  expect_equal(cf, par)
  expect_length(par, 2)
})
  
test_that("transform_y works", {
  y_transf1 <- transf_y(0)
  y_transf2 <- transf_y(1)
  y_transf3 <- transf_y(0.5)
  
  expect_equal(y_transf1, stats::qnorm(0.001))
  expect_equal(y_transf2, stats::qnorm(0.999))
  expect_equal(y_transf3, stats::qnorm(0.5))
})

test_that("fit_mod_3pod works", {
  set.seed(20210330)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }

  x1 <- 5:30
  y1 <- sapply(x1, function(x) fun_ttest(n = x, k = 50))
  y_mat1 <- cbind(y1 * 50, (1 - y1) * 50)

  x2 <- 200:250
  y2 <- sapply(x2, function(x) fun_ttest(n = x, k = 50))
  y_mat2 <- cbind(y2 * 50, (1 - y2) * 50)

  mod1 <- glm(y_mat1 ~ x1, family = stats::binomial(link = "probit"))
  cf1 <- unname(stats::coef(mod1))
  par1 <- unname(stats::coef(fit_mod_3pod(x = x1, 
    y = y1, k = 50, alpha = 0.025)))

  off.par <- rep(stats::qnorm(0.025), length(x2))
  mod2 <- suppressWarnings(glm(y_mat2 ~ x2 - 1 + offset(off.par), 
    family = stats::binomial(link = "probit")))
  cf2 <- unname(stats::coef(mod2))
  par2 <- suppressWarnings(unname(stats::coef(fit_mod_3pod(x = x2, 
    y = y2, k = 50, alpha = 0.025))))
  
  expect_equal(cf1, par1)
  expect_equal(cf2, par2)
  expect_length(par1, 2)
  expect_length(par2, 1)
})
  
  