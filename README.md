
<!-- README.md is generated from README.Rmd. -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/lbau7/findn/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lbau7/findn/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/lbau7/findn/branch/master/graph/badge.svg)](https://codecov.io/gh/lbau7/findn?branch=master)
<!-- badges: end -->

# Introduction

findn estimates the sample size for a certain target functions that
simulates a statistical test or a trial. The main function of the
package is findn, which uses the Bayesian Local Linear Regression (BLL)
algorithm. The BLL algorithm uses repeated simulations to estimate the
power for various sample sizes and fits a weighted probit regression
model with Bayesian parameter estimates to all simulation outcomes to
estimate the sample size that corresponds to the target power.

# Installation

To install findn use:

``` r
# install.packages("devtools")
devtools::install_github("lbau7/textools")
```

# Usage

findn can be used to estimate the sample size for a statistical test or
a trial when no closed formula to calculate the power or the sample size
is available. To use findn we have to provide a target function that
simulates the power of the test or the trial for a certain sample size.
The target function has to take at least two arguments: n, the sample
size and k, the number of simulations that the power estimate should be
based on. The target function must return an estimate for the power at
sample size n, i.e. the proportion of trials for which the goal of the
trial was achieved.

Consider for example a trial where three independent hypotheses are
tested with a two-sample t-test. We want to control the family wise
error rate by using the Bonferroni-Holm procedure. We are interested to
estimate the sample size that corresponds to a target power of 80%,
where power is defined as the probability to reject at least one of the
three null hypotheses.

``` r
sim_bh <- function(n, k) {
  success <- numeric(k)
  for (i in 1:k) {
    obs.exp1 <- rnorm(n = n / 6, mean = 20, sd = 5)
    obs.cont1 <- rnorm(n = n / 6, mean = 22, sd = 5)
    p1 <- t.test(obs.exp1, obs.cont1)$p.value
    
    obs.exp2 <- rnorm(n = n / 6, mean = 30, sd = 8)
    obs.cont2 <- rnorm(n = n / 6, mean = 34, sd = 8)
    p2 <- t.test(obs.exp2, obs.cont2)$p.value
    
    obs.exp3 <- rnorm(n = n / 6, mean = 50, sd = 10)
    obs.cont3 <- rnorm(n = n / 6, mean = 56, sd = 10)
    p3 <- t.test(obs.exp3, obs.cont3)$p.value
    
    success[i] <- sum(p.adjust(p = c(p1, p2, p3), method = "holm") <= 0.05) > 0
  }
  return(mean(success))
}
```

Note that n corresponds to the total sample size in this example, so
it’s easy to adapt the target function for instance for a design with
unbalanced sample size allocation. But it’s also possible to set up the
target function such that n corresponds to the sample size per group.

Now we can use findn to estimate the sample size:

``` r
res_bh <- findn(fun = sim_bh, targ = 0.8, start = 100)
res_bh
# $Point_Estimate_n
# [1] 224
# 
# $Minimum_Sufficient_n
# [1] "235"
# 
# $Message
# [1] "Normal completion"
```

Point_Estimate_n corresponds to the smallest sample size for which the
estimated power exceeds the target power, while Minimum_Sufficient_n
corresponds to the smallest sample size for which the lower limit of the
95% confidence interval for the estimated power exceeds the target
power. If a different level for the confidence interval is desired, the
parameter level can be changed. For more details about the estimated
power values and their confidence intervals for various sample sizes we
can use the print function and set details to “high”.

``` r
print(res_bh, details = "high")
# $Point_Estimate_n
# [1] 219
# 
# $Details
#       n Est.Power Lower.CL Upper.CL     Rating
# 206 206     0.776    0.757    0.795    Too Low
# 207 207     0.778    0.759    0.796    Too Low
# 208 208     0.780    0.761    0.798    Too Low
# 209 209     0.782    0.763    0.800  Uncertain
# 210 210     0.784    0.765    0.802  Uncertain
# 211 211     0.786    0.767    0.804  Uncertain
# 212 212     0.788    0.769    0.805  Uncertain
# 213 213     0.789    0.771    0.807  Uncertain
# 214 214     0.791    0.773    0.809  Uncertain
# 215 215     0.793    0.774    0.811  Uncertain
# 216 216     0.795    0.776    0.812  Uncertain
# 217 217     0.797    0.778    0.814  Uncertain
# 218 218     0.798    0.780    0.816  Uncertain
# 219 219     0.800    0.782    0.818  Uncertain
# 220 220     0.802    0.784    0.819  Uncertain
# 221 221     0.804    0.786    0.821  Uncertain
# 222 222     0.805    0.787    0.823  Uncertain
# 223 223     0.807    0.789    0.824  Uncertain
# 224 224     0.809    0.791    0.826  Uncertain
# 225 225     0.811    0.793    0.828  Uncertain
# 226 226     0.812    0.794    0.829  Uncertain
# 227 227     0.814    0.796    0.831  Uncertain
# 228 228     0.816    0.798    0.832  Uncertain
# 229 229     0.817    0.799    0.834  Uncertain
# 230 230     0.819    0.801    0.836 Sufficient
# 231 231     0.820    0.803    0.837 Sufficient
# 232 232     0.822    0.804    0.839 Sufficient
# 
# $Message
# [1] "Normal completion"
```

The rating “Too Low” means, that the upper limit of the confidence
interval doesn’t exceed the target power. “Uncertain” means, that the
confidence interval contains the target power. “Sufficient” means that
the lower limit of the confidence interval exceeds the target power.
