
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sfaR: Stochastic Frontier Analysis Using R

<!-- badges: start -->

[![CodeFactor](https://www.codefactor.io/repository/github/hdakpo/sfaR/badge)](https://www.codefactor.io/repository/github/hdakpo/sfaR)
[![R-CMD-check](https://github.com/hdakpo/sfaR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hdakpo/sfaR/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/sfaR)](https://CRAN.R-project.org/package=sfaR)
[![](https://img.shields.io/badge/devel%20version-1.0.0-darkred.svg)](https://github.com/hdakpo/sfaR)
[![](https://img.shields.io/badge/license-GPL-blue)](https://github.com/hdakpo/sfaR)
[![Downloads](https://cranlogs.r-pkg.org/badges/sfaR)](https://CRAN.R-project.org/package=sfaR)
[![](https://img.shields.io/github/languages/code-size/hdakpo/sfaR.svg)](https://github.com/hdakpo/sfaR)
<!-- badges: end -->

*sfaR* provides a set of tools (maximum likelihood and maximum simulated
likelihood) for various specifications of stochastic frontier analysis.

Three categories of models are available in *sfaR*:

1.  **Classic Stochastic Frontier Model**

This model allows the estimation of the frontier for a cross-sectional
or pooled data. Basically we have

$$y_i = \mathbf{x_i'}\boldsymbol{\beta} + v_i - Su_i$$

where $S = 1$ for production function and $S = -1$ for cost function.
$v$ follows a normal distribution $\mathcal{N}(0, \sigma_v^2)$. For $u$
ten different distributions are available. These distributions include:

- Half-Normal
- Truncated Normal
- Exponential
- Rayleigh
- Gamma
- Generalized Exponential
- Lognormal
- Truncated Skewed Laplace
- Uniform
- Weibull

In the case of the Gamma, lognormal and Weibull distributions, maximum
simulated likelihood is used with the possibility of four possibilities
to construct the draws: Halton, Generalized Halton, Sobol and uniform.

Heteroscedasticity in both error terms can be implemented, in addition
to heterogeneity in the truncated mean parameter in the case of the
truncated normal and lognormal distributions. In addition, in the case
of the truncated normal distribution, the scaling property can be
estimated. The main function for this class of model is `sfacross`.

2.  **Latent Class Stochastic Frontier Model**

This model accounts for technological heterogeneity by splitting the
observations into a maximum number of five classes. The classification
operates based on a logit functional form that can be specified using
some covariates (namely, the separating variables allowing the
separation of observations in several classes). Only the half normal
distribution is available for the one-sided error term.
Heteroscedasticity in both error terms is possible. The choice of the
number of classes can be guided by several information criteria
(i.e. AIC, BIC or HQIC). The main function for this class of model is
`sfalcmcross`.

3.  **Sample Selection Correction Stochastic Frontier Model**

This model solves the selection bias due to the correlation between the
two-sided errors terms in both the selection and the frontier equations,
in the case of cross-sectional or pooled data. The main function for
this class of model is `sfaselectioncross`.

An important features of *sfaR* is to provide eleven different
optimization algorithms. For complex problem, several algorithms can be
combined especially non-gradient based in a first step.

## Installation

You can install the development version of sfaR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hdakpo/sfaR")
```

Install the current version on CRAN with

``` r
# install.packages("sfaR")
```

## Example

This subsection provides set of examples introducing some important
features of *sfaR*.

``` r
library(sfaR)
#> * Please cite the 'sfaR' package as:
#>   Dakpo KH., Desjeux Y., Henningsen A., and Latruffe L. (2023). sfaR: Stochastic Frontier Analysis Routines. R package version 1.0.0.
#> 
#> See also: citation("sfaR")
#> 
#> * For any questions, suggestions, or comments on the 'sfaR' package, please make use of Tracker facilities at:
#>   https://github.com/hdakpo/sfaR/issues
## basic examples code

## let's estimate the classic frontier for different distributions using 
## the utility dataset, which contains data on fossil fuel fired 
## steam electric power generation plants in the United States
hlf <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
 log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
 I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
 udist = 'hnormal', uhet = ~ regu, data = utility, S = -1, method = 'bfgs')

trnorm <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
 log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
 I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
 udist = 'tnormal', muhet = ~ regu, data = utility, S = -1, method = 'bfgs')

tscal <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
 log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
 I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
 udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, 
 S = -1, method = 'bfgs', scaling = TRUE)

expo <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
 log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
 I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
 udist = 'exponential', uhet = ~ regu, data = utility, S = -1, method = 'bfgs')
```

Outputs of estimation can be exported using the *texreg* package. For
instance, using the command `screenreg(list(hlf, trnorm, tscal, expo))`
yields the following output

![sfacross](https://user-images.githubusercontent.com/29732089/235988357-90a74e12-7695-47ae-8b29-3591ca221bcd.png)

``` r
## For the latent class stochastic frontier we have:
lcm2c1 <- sfalcmcross(formula = ly ~ lk + ll + yr, thet = ~initStat, 
 data = worldprod)
#> Initialization: SFA + halfnormal - normal distributions...
#> LCM 2 Classes Estimation...
lcm2c2 <- sfalcmcross(formula = ly ~ lk + ll + yr, uhet = ~initStat, 
 data = worldprod)
#> Initialization: SFA + halfnormal - normal distributions...
#> LCM 2 Classes Estimation...
```

The command `screenreg(list(lcm2c1, lcm2c2))` generates the following

![sfalcmcross](https://user-images.githubusercontent.com/29732089/236163537-dc12e886-84c2-49d4-a943-9d61cfb82000.png)

``` r
## The following simulation is used for the sample selection
 N <- 2000  # sample size
 set.seed(12345)
 z1 <- rnorm(N)
 z2 <- rnorm(N)
 v1 <- rnorm(N)
 v2 <- rnorm(N)
 e1 <- v1
 e2 <- 0.7071 * (v1 + v2)
 ds <- z1 + z2 + e1
 d <- ifelse(ds > 0, 1, 0)
 u <- abs(rnorm(N))
 x1 <- rnorm(N)
 x2 <- rnorm(N)
 y <- x1 + x2 + e2 - u
 data <- cbind(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2, d = d)
 
 ## Estimation using quadrature (Gauss-Kronrod)
 
 selecRes1 <- sfaselectioncross(selectionF = d ~ z1 + z2, frontierF = y ~ x1 + x2, 
 modelType = 'greene10', method = 'bfgs',
 logDepVar = TRUE, data = as.data.frame(data),
 S = 1L, udist = 'hnormal', lType = 'kronrod', Nsub = 100, uBound = Inf,
 simType = 'halton', Nsim = 300, prime = 2L, burn = 10, antithetics = FALSE,
 seed = 12345, itermax = 2000, printInfo = FALSE)
#> First step probit model...
#> Second step Frontier model...
 
 ## Estimation using quadrature (Gauss-Hermite)
 
 selecRes2 <- sfaselectioncross(selectionF = d ~ z1 + z2, frontierF = y ~ x1 + x2, 
 modelType = 'greene10', method = 'bfgs',
 logDepVar = TRUE, data = as.data.frame(data),
 S = 1L, udist = 'hnormal', lType = 'ghermite', Nsub = 100, uBound = Inf,
 simType = 'halton', Nsim = 300, prime = 2L, burn = 10, antithetics = FALSE,
 seed = 12345, itermax = 2000, printInfo = FALSE)
#> First step probit model...
#> Second step Frontier model...
```

The command `screenreg(list(selecRes1, selecRes2))`

![selectioncross](https://user-images.githubusercontent.com/29732089/236043261-af06f359-bbfc-46e0-8bde-2951cc80813c.png)
