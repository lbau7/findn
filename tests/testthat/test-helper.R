# Helper function for findn_maruo
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
  
  
  