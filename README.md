
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

``` r
katex::katex_html(katex::example_math())
```

<span class="katex-display"><span class="katex"><span class="katex-mathml"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><semantics><mrow><mi>f</mi><mo stretchy="false">(</mo><mi>x</mi><mo stretchy="false">)</mo><mo>=</mo><mfrac><mn>1</mn><mrow><mi>σ</mi><msqrt><mrow><mn>2</mn><mi>π</mi></mrow></msqrt></mrow></mfrac><msup><mi>e</mi><mrow><mo>−</mo><mfrac><mn>1</mn><mn>2</mn></mfrac><mo stretchy="false">(</mo><mfrac><mrow><mi>x</mi><mo>−</mo><mi>μ</mi></mrow><mi>σ</mi></mfrac><msup><mo stretchy="false">)</mo><mn>2</mn></msup></mrow></msup></mrow><annotation encoding="application/x-tex">f(x)= {\frac{1}{\sigma\sqrt{2\pi}}}e^{- {\frac {1}{2}} (\frac {x-\mu}{\sigma})^2}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height:1em;vertical-align:-0.25em;"></span><span class="mord mathnormal" style="margin-right:0.10764em;">f</span><span class="mopen">(</span><span class="mord mathnormal">x</span><span class="mclose">)</span><span class="mspace" style="margin-right:0.2777777777777778em;"></span><span class="mrel">=</span><span class="mspace" style="margin-right:0.2777777777777778em;"></span></span><span class="base"><span class="strut" style="height:2.25144em;vertical-align:-0.93em;"></span><span class="mord"><span class="mord"><span class="mopen nulldelimiter"></span><span class="mfrac"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:1.32144em;"><span style="top:-2.2027799999999997em;"><span class="pstrut" style="height:3em;"></span><span class="mord"><span class="mord mathnormal" style="margin-right:0.03588em;">σ</span><span class="mord sqrt"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.90722em;"><span class="svg-align" style="top:-3em;"><span class="pstrut" style="height:3em;"></span><span class="mord" style="padding-left:0.833em;"><span class="mord">2</span><span class="mord mathnormal" style="margin-right:0.03588em;">π</span></span></span><span style="top:-2.86722em;"><span class="pstrut" style="height:3em;"></span><span class="hide-tail" style="min-width:0.853em;height:1.08em;"><svg width='400em' height='1.08em' viewBox='0 0 400000 1080' preserveAspectRatio='xMinYMin slice'><path d='M95,702
c-2.7,0,-7.17,-2.7,-13.5,-8c-5.8,-5.3,-9.5,-10,-9.5,-14
c0,-2,0.3,-3.3,1,-4c1.3,-2.7,23.83,-20.7,67.5,-54
c44.2,-33.3,65.8,-50.3,66.5,-51c1.3,-1.3,3,-2,5,-2c4.7,0,8.7,3.3,12,10
s173,378,173,378c0.7,0,35.3,-71,104,-213c68.7,-142,137.5,-285,206.5,-429
c69,-144,104.5,-217.7,106.5,-221
l0 -0
c5.3,-9.3,12,-14,20,-14
H400000v40H845.2724
s-225.272,467,-225.272,467s-235,486,-235,486c-2.7,4.7,-9,7,-19,7
c-6,0,-10,-1,-12,-3s-194,-422,-194,-422s-65,47,-65,47z
M834 80h400000v40h-400000z'/></svg></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.13278em;"><span></span></span></span></span></span></span></span><span style="top:-3.23em;"><span class="pstrut" style="height:3em;"></span><span class="frac-line" style="border-bottom-width:0.04em;"></span></span><span style="top:-3.677em;"><span class="pstrut" style="height:3em;"></span><span class="mord"><span class="mord">1</span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.93em;"><span></span></span></span></span></span><span class="mclose nulldelimiter"></span></span></span><span class="mord"><span class="mord mathnormal">e</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:1.0369199999999998em;"><span style="top:-3.4130000000000003em;margin-right:0.05em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mtight">−</span><span class="mord mtight"><span class="mord mtight"><span class="mopen nulldelimiter sizing reset-size3 size6"></span><span class="mfrac"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.8443142857142858em;"><span style="top:-2.656em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size3 size1 mtight"><span class="mord mtight"><span class="mord mtight">2</span></span></span></span><span style="top:-3.2255000000000003em;"><span class="pstrut" style="height:3em;"></span><span class="frac-line mtight" style="border-bottom-width:0.049em;"></span></span><span style="top:-3.384em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size3 size1 mtight"><span class="mord mtight"><span class="mord mtight">1</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.344em;"><span></span></span></span></span></span><span class="mclose nulldelimiter sizing reset-size3 size6"></span></span></span><span class="mopen mtight">(</span><span class="mord mtight"><span class="mopen nulldelimiter sizing reset-size3 size6"></span><span class="mfrac"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height:0.87905em;"><span style="top:-2.656em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size3 size1 mtight"><span class="mord mtight"><span class="mord mathnormal mtight" style="margin-right:0.03588em;">σ</span></span></span></span><span style="top:-3.2255000000000003em;"><span class="pstrut" style="height:3em;"></span><span class="frac-line mtight" style="border-bottom-width:0.049em;"></span></span><span style="top:-3.4623857142857144em;"><span class="pstrut" style="height:3em;"></span><span class="sizing reset-size3 size1 mtight"><span class="mord mtight"><span class="mord mathnormal mtight">x</span><span class="mbin mtight">−</span><span class="mord mathnormal mtight">μ</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height:0.344em;"><span></span></span></span></span></span><span class="mclose nulldelimiter sizing reset-size3 size6"></span></span><span class="mclose mtight"><span class="mclose mtight">)</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height:0.8913142857142857em;"><span style="top:-2.931em;margin-right:0.07142857142857144em;"><span class="pstrut" style="height:2.5em;"></span><span class="sizing reset-size3 size1 mtight"><span class="mord mtight">2</span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span></span>
<!-- $$y_i = \mathbf{x_i'}\boldsymbol{\beta} + v_i - Su_i$$ -->

where *S* = 1 for production function and *S* =  − 1 for cost function.
*v* follows a normal distribution 𝒩(0,*σ*<sub>*v*</sub><sup>2</sup>).
For *u* ten different distributions are available. These distributions
include:

-   Half-Normal
-   Truncated Normal
-   Exponential
-   Rayleigh
-   Gamma
-   Generalized Exponential
-   Lognormal,
-   Truncated Skewed Laplace
-   Uniform
-   Weibull

In the case of the Gamma, lognormal and Weibull distributions, maximum
simulated likelihood is used with the possibility of four specific
distributions to construct the draws: Halton, Generalized Halton, Sobol
and uniform.

Heteroscedasticity in both error terms can be implemented, in addition
to heterogeneity in the truncated mean parameter in the case of the
truncated normal and lognormal distributions. In addition, in the case
of the truncated normal distribution, the scaling property can be
estimated.

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
