
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sfaR

<!-- badges: start -->
<!-- badges: end -->

*sfaR* provides a set of tools (maximum likelihood and maximum simulated
likelihood) for various specifications of stochastic frontier analysis.

Three categories of models are available in *sfaR*:

-   **Classic Stochastic Frontier Model**: This model allows the
    estimation of the frontier for a cross-sectional or pooled data.
    Basically we have

-   **Latent Class Stochastic Frontier Model** (LCM): This model
    accounts for technological heterogeneity by splitting the
    observations into a maximum number of five classes. The
    classification operates based on a logit functional form that can be
    specified using some covariates (namely, the separating variables
    allowing the separation of observations in several classes). Only
    the half normal distribution is available for the one-sided error
    term. Heteroscedasticity in both error terms is possible. The choice
    of the number of classes can be guided by several information
    criteria (i.e. AIC, BIC or HQIC).

-   **Sample Selection Correction Stochastic Frontier Model**: This
    model solves the selection bias due to the correlation between the
    two-sided errors terms in both the selection and the frontier
    equations, in the case of cross-sectional or pooled data.

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
#>   Dakpo KH., Desjeux Y., and Latruffe L. (2022). sfaR: Stochastic Frontier Analysis Routines. R package version 1.0.0.
#> See also: citation("sfaR")
#> 
#> * For any questions, suggestions, or comments on the 'sfaR' package, please make use of Tracker facilities at:
#>   https://github.com/hdakpo/sfaR/issues
## basic example code
```
