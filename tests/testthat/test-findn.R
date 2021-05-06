test_that("findn finds approximately the correct sample size", {
  set.seed(20210426)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_bll1 <- suppressWarnings(findn(fun = fun_ttest, targ = 0.8, start = 100))
  res_bll2 <- suppressWarnings(findn(fun = fun_ttest, targ = 0.9, start = 100)) 
  
  n_true1 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.8)$n)
  n_true2 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.9)$n)
  
  expect_equal(res_bll1$sample_size, n_true1, tol = 2)
  expect_equal(res_bll2$sample_size, n_true2, tol = 2)
})