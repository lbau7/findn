test_that("findn_3pod finds approximately the correct sample size", {
  set.seed(20210426)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_3pod1 <- suppressWarnings(findn_3pod(fun = fun_ttest, targ = 0.8, 
    start = 100, k = 50))
  res_3pod2 <- suppressWarnings(findn_3pod(fun = fun_ttest, targ = 0.9, 
    start = 100, k = 50))
  
  n_true1 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.8)$n)
  n_true2 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.9)$n)
  
  expect_equal(res_3pod1$Point_Estimate_n, n_true1, tol = 2)
  expect_equal(res_3pod2$Point_Estimate_n, n_true2, tol = 2)
})