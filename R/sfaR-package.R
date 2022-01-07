################################################################################
#                                                                              #
# sfaR package doc                                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# sfaR package overview                                                        #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' sfaR: A package for estimating stochastic frontier models
#'
#' The \pkg{sfaR} package provides a set of tools (maximum likelihood - ML and
#' maximum simulated likelihood - MSL) for various specifications of stochastic
#' frontier analysis (SFA).
#'
#' Two categories of important functions are available: \code{\link{sfacross}}
#' and \code{\link{lcmcross}}, which estimate different types of frontiers and
#' offer eleven alternative optimization algorithms (i.e. "bfgs", "bhhh", "nr",
#' "nm", "cg", "sann", "ucminf", "mla", "sr1", "sparse", "nlminb").
#'
#' @name sfaR-package
#'
#' @aliases sfaR-package sfaR
#'
#' @docType package
#'
#' @section lcmcross: \code{\link{lcmcross}} estimates latent class stochastic
#' frontier models (LCM), which accounts for technological heterogeneity by
#' splitting the observations into a maximum number of five classes. The
#' classification operates based on a logit functional form that can be
#' specified using some covariates (namely, the separating variables allowing
#' the separation of observations in several classes). Only the half normal
#' distribution is available for the one-sided error term. Heteroscedasticity
#' in both error terms is possible. The choice of the number of classes can be
#' guided by several information criteria (i.e. AIC, BIC or HQIC).
#'
#' @section sfacross: \code{\link{sfacross}} estimates the frontier for cross-sectional data
#'   and allows for ten different distributions for the one-sided error term.
#'   These distributions include the exponential, the Gamma, the generalized
#'   exponential, the half normal, the lognormal, the truncated normal, the
#'   truncated skewed Laplace, the Rayleigh, the uniform and the Weibull distributions.
#'   In the case of the Gamma, lognormal and Weibull distributions, maximum simulated
#'   likelihood (MSL) is used with the possibility of four specific distributions to
#'   construct the draws: Halton, Generalized Halton, Sobol and uniform.
#'   Heteroscedasticity in both error terms can be implemented, in addition to
#'   heterogeneity in the truncated mean parameter in the case of the truncated
#'   normal and lognormal distributions. In addition, in the case of the truncated normal
#'   distribution, the scaling property can be estimated.
#'
#' @section Bugreport: Any bug or suggestion can be reported using the \code{sfaR}'\code{s}
#' tracker facilities at: \url{https://github.com/hdakpo/sfaR/issues}
#'
#' @author K Hervé Dakpo, Yann Desjeux and Laure Latruffe
#'
#' @importFrom stats coefficients dnorm lm model.frame
#' @importFrom stats model.matrix model.response nlminb
#' @importFrom stats pnorm qnorm delete.response fitted
#' @importFrom stats logLik residuals terms vcov formula
#' @importFrom stats integrate runif model.weights
#' @importFrom stats na.pass printCoefmat pt qt
#' @importFrom methods as
#' @importFrom Formula as.Formula Formula
#' @importFrom primes generate_primes is_prime
#' @importFrom randtoolbox sobol
#' @importFrom qrng ghalton
#' @importFrom MASS ginv
#' @importFrom gsl erf erfc
#' @import maxLik
#' @importFrom ucminf ucminf
#' @importFrom trustOptim trust.optim
#' @importFrom marqLevAlg mla
#' @importFrom dplyr bind_cols as_tibble mutate select
#' @importFrom numDeriv jacobian
#' @importFrom nleqslv nleqslv
#' @importFrom emdbook qchibarsq
#' @importFrom fBasics dagoTest
#' @importFrom cubature hcubature
#' @importFrom katex math_to_rd
NULL
