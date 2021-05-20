test_that("plot.findn works", {
  set.seed(20210520)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_bll <- suppressWarnings(findn(fun = fun_ttest, targ = 0.8, start = 100))
  plot_bll1 <- plot(res_bll, min_n = 1, power_lim = 0.95)
  plot_bll2 <- plot(res_bll, min_n = 1, max_n = 100)
  
  expect_is(plot_bll1, "ggplot")
  expect_gt(max(plot_bll1$data$Lower.CL), 0.95)
  expect_true(all(plot_bll1$data$Upper.CL[which(
    plot_bll1$data$Rating == "Too Low")] < 0.8))
  expect_true(all(plot_bll1$data$Lower.CL[which(
    plot_bll1$data$Rating == "Sufficient")] > 0.8))
  expect_equal(max(plot_bll2$data$n), 100)
  expect_equal(min(plot_bll2$data$n), 1)
})

  
  