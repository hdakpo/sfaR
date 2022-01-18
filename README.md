
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sfaR

<!-- badges: start -->
<!-- badges: end -->

*sfaR* provides a set of tools (maximum likelihood and maximum simulated
likelihood) for various specifications of stochastic frontier analysis.

Three categories of models are available in *sfaR*:

-   **Classic Stochastic Frontier Model**: This model allows the
    estimation of the frontier for a cross-sectional or pooled data.

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
install.packages("sfaR")
#> Installation du package dans 'C:/Users/Dakpo/AppData/Local/Temp/RtmpoNSUnI/temp_libpath9d010ca2ffe'
#> (car 'lib' n'est pas spécifié)
#> package 'sfaR' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\Dakpo\AppData\Local\Temp\RtmpGiJWtA\downloaded_packages
```

## Example

This subsection provides set of examples introducing some important
features of *sfaR*.

``` r
library(sfaR)
#> * Please cite the 'sfaR' package as:
#>   Dakpo KH., Desjeux Y. and Latruffe L. (2021). sfaR: Stochastic Frontier Analysis using R. R package version 0.1.0.
#> See also: citation("sfaR")
#> 
#> * For any questions, suggestions, or comments on the 'sfaR' package, please make use of Tracker facilities at:
#>   https://r-forge.r-project.org/projects/sfar/
## basic example code
```
