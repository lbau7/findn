test_that("findn_rm finds approximately the correct sample size", {
  set.seed(20210330)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_rm1 <- findn_rm(fun = fun_ttest, targ = 0.8, start = 100, k = 50, 
    max_evals = 1000, j = 1, avg = 0.5)
  res_rm2 <- findn_rm(fun = fun_ttest, targ = 0.9, start = 100, k = 50, 
    max_evals = 1000, j = 1, avg = 0.5)
  
  n_true1 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.8)$n)
  n_true2 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.9)$n)
  
  expect_equal(res_rm1$Point_Estimate_n, n_true1, tol = 2)
  expect_equal(res_rm1$Point_Estimate_n_avg, n_true1, tol = 2)
  expect_equal(res_rm2$Point_Estimate_n, n_true2, tol = 2)
  expect_equal(res_rm2$Point_Estimate_n_avg, n_true2, tol = 2)
})

test_that("findn_rm stops when it's supposed to", {
  expect_error(findn_rm(fun = function(n, k) n, targ = 0.8, start = 100, k = 50,
    max_evals = 1000, j = 0.4, avg = 0.5))
  expect_error(findn_rm(fun = function(n, k) x, targ = 0.8, start = 100, k = 50,
    max_evals = 1000, j = 1, avg = 2))
  expect_error(findn_rm(fun = function(x) x, targ = 0.8, start = 100, k = 350,
    max_evals = 1000, j = 1, avg = 1))
  expect_error(findn_rm(fun = function(x) x, targ = 0.8, start = 100, k = 99,
    max_evals = 1000, j = 1, avg = 1))
})