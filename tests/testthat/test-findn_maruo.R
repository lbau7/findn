test_that("findn_maruo finds approximately the correct sample size", {
  set.seed(20210330)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_maruo1 <- findn_maruo(fun = fun_ttest, targ = 0.8, start = 10, k = 100)
  res_maruo2 <- findn_maruo(fun = fun_ttest, targ = 0.9, start = 10, k = 100)
  
  n_true1 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.8)$n)
  n_true2 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.9)$n)
  
  expect_equal(res_maruo1$Point_Estimate_n, n_true1, tolerance = 2)
  expect_equal(res_maruo2$Point_Estimate_n, n_true2, tolerance = 2)
})

test_that("findn_maruo stop when power is neither 0.8 nor 0.9", {
  expect_error(findn_maruo(fun = function(n, k) n, targ = 0.85, start = 10,
    k = 100))
})

test_that("findn_maruo increases step-size correctly", {
  set.seed(20210506)
  fun_ttest1 <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 1.5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  fun_ttest2 <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 10, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  fun_ttest3 <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 6, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  pow1 <- fun_ttest1(10, k = 100)
  pow2 <- fun_ttest2(10, k = 100)
  pow3 <- fun_ttest3(10, k = 100)
  
  res_maruo1 <- findn_maruo(fun = fun_ttest1, targ = 0.8, start = 10, k = 100)
  res_maruo2 <- findn_maruo(fun = fun_ttest2, targ = 0.8, start = 10, k = 100)
  res_maruo3 <- findn_maruo(fun = fun_ttest3, targ = 0.8, start = 10, k = 100)
  
  expect_equal(res_maruo1$All_n[2] + 5, res_maruo1$All_n[3])
  expect_equal(res_maruo1$All_n[length(res_maruo1$All_n) - 1] + 10, 
    res_maruo1$All_n[length(res_maruo1$All_n)])
  expect_equal(res_maruo2$All_n[2] + 1, res_maruo2$All_n[3])
  expect_equal(res_maruo3$All_n[2] + 2, res_maruo3$All_n[3])
})
  