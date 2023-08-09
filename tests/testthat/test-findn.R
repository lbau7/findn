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
  
  expect_equal(length(res_bll1$all_evals) * 25, 2000)
  expect_equal(res_bll1$sample_size, n_true1, tol = 2)
  expect_equal(res_bll2$sample_size, n_true2, tol = 2)
})
  
test_that("findn works when the starting value is way off", {
  set.seed(20210514)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  res_bll1 <- suppressWarnings(findn(fun = fun_ttest, targ = 0.8, 
    start = 5000))
  n_true1 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.8)$n)
  
  expect_equal(res_bll1$sample_size, n_true1, tol = 2)
})

test_that("stopping rules of findn work", {
  set.seed(20210514)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_bll1 <- suppressWarnings(findn(fun = fun_ttest, targ = 0.8, 
    start = 100, stop = "power_ci", power_ci_tol = 0.02))
  det_bll1 <- print(res_bll1, details = "high", invisible = TRUE)$Details
  len_bll1 <- diff(as.numeric(det_bll1[which(det_bll1$n == 
      res_bll1$sample_size), c(2, 4)]))
  
  res_bll2 <- suppressWarnings(findn(fun = fun_ttest, targ = 0.8, 
    start = 100, stop = "abs_unc", abs_unc_tol = 10))
  det_bll2 <- print(res_bll2, details = "high", invisible = TRUE)$Details
  len_bll2 <- nrow(det_bll2[which(det_bll2$Ratin == "Uncertain"), ])
  
  res_bll3 <- suppressWarnings(findn(fun = fun_ttest, targ = 0.8, 
    start = 100, stop = "rel_unc", rel_unc_tol = 0.1))
  det_bll3 <- print(res_bll3, details = "high", invisible = TRUE)$Details
  min_bll3 <- as.numeric(det_bll3[min(which(det_bll3$Rating == 
      "Uncertain")), 1])
  max_bll3 <- as.numeric(det_bll3[max(which(det_bll3$Rating == 
      "Uncertain")), 1])
  len_bll3 <- (max_bll3 - min_bll3) / min_bll3
  
  expect_lt(len_bll1, 0.02)
  expect_lte(len_bll2, 10)
  expect_lte(len_bll3, 0.1)
})

test_that("findn returns estimate after every iteration when verbose = TRUE", {
  set.seed(20210514)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_bll <- capture.output(suppressWarnings(findn(fun = fun_ttest, targ = 0.8, 
    start = 100, k = 25, init_evals = 100, max_evals = 2000, verbose = TRUE)))

  expect_equal(length(res_bll) - 9, (2000 - 100) / 25 / 2)
})
  
test_that("findn stops when it's supposed to", {
  set.seed(20210514)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  expect_error(findn(fun = function(n, k) n, targ = 0.8, start = 100,
    k = 50, init_evals = 25))
  expect_error(findn(fun = function(n, k) n, targ = 0.8, start = 100,
    k = 50, max_evals = 100))
  expect_error(suppressWarnings(findn(fun = fun_ttest, targ = 0.8,
    start = 5000, var_beta = 1000, var_alpha = 1000)))
})

test_that("findn shows correct exit message", {
  set.seed(20210514)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_bll1 <- suppressWarnings(findn(fun = fun_ttest, targ = 0.8, start = 100,
    stop = "power_ci", power_ci_tol = 0.05))
  res_bll2 <- suppressWarnings(findn(fun = fun_ttest, targ = 0.8, start = 100,
    stop = "power_ci", power_ci_tol = 0.001))
  
  
  expect_equal(res_bll1$exit.mes, "Stopping criterion fulfilled")
  expect_equal(res_bll2$exit.mes,
    "Maximum number of evaluations reached without sufficient accuracy")
})
  