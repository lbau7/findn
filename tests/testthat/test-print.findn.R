test_that("plot.findn works", {
  set.seed(20210520)
  fun_ttest1 <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  fun_ttest2 <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 2.5), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_bll1 <- suppressWarnings(findn(fun = fun_ttest1, targ = 0.8, start = 100))
  det_bll11 <- print(res_bll1, details = "low", invisible = TRUE)
  pred_bll <- predict_fit(res_bll1$fit, det_bll11$Minimum_Sufficient_n, 
    se = TRUE)
  low_cl <- stats::pnorm(pred_bll$pred - stats::qnorm(1 - (res_bll1$level / 
      2)) * pred_bll$se)
  det_bll12 <- print(res_bll1, details = "high", max_n = 10, invisible = TRUE)
  det_bll13 <- print(res_bll1, details = "high", invisible = TRUE)
  res_bll2 <- suppressWarnings(findn(fun = fun_ttest2, targ = 0.8, start = 100,
    max_evals = 3000))
  det_bll2 <- print(res_bll2, details = "high", invisible = TRUE)
    
  expect_gt(low_cl, 0.8)
  expect_equal(nrow(det_bll12$Details), 10)
  expect_equal(sum(det_bll13$Details$Rating == "Too Low"), 3)
  expect_equal(sum(det_bll13$Details$Rating == "Sufficient"), 3)
  expect_lt(max(as.numeric(det_bll13$Details$Upper.CL[which(
    det_bll13$Details$Rating == "Too Low")])), 0.8)
  expect_gt(max(as.numeric(det_bll13$Details$Lower.CL[which(
    det_bll13$Details$Rating == "Sufficient")])), 0.8)
  expect_equal(sum(det_bll2$Details$Rating == "Uncertain"), 0)
})

  