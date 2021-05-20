test_that("plot.findn works", {
  set.seed(20210520)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_bll <- suppressWarnings(findn(fun = fun_ttest, targ = 0.8, start = 100))
  det_bll1 <- print(res_bll, details = "low", invisible = TRUE)
  pred_bll <- predict_fit(res_bll$fit, det_bll1$Minimum_Sufficient_n, se = TRUE)
  low_cl <- stats::pnorm(pred_bll$pred - stats::qnorm(1 - (res_bll$level / 2)) * 
      pred_bll$se)
  det_bll2 <- print(res_bll, details = "high", max_n = 10, invisible = TRUE)
  det_bll3 <- print(res_bll, details = "high", invisible = TRUE)
    
  expect_gt(low_cl, 0.8)
  expect_equal(nrow(det_bll2$Details), 10)
  expect_equal(sum(det_bll3$Details$Rating == "Too Low"), 3)
  expect_equal(sum(det_bll3$Details$Rating == "Sufficient"), 3)
  expect_lt(max(as.numeric(det_bll3$Details$Upper.CL[which(
    det_bll3$Details$Rating == "Too Low")])), 0.8)
  expect_gt(max(as.numeric(det_bll3$Details$Lower.CL[which(
    det_bll3$Details$Rating == "Sufficient")])), 0.8)
})

  