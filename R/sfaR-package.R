################################################################################
#                                                                              #
# sfaR package doc                                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# sfaR package overview                                                        #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Sample selection correction for Stochastic Frontier Model         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' sfaR: A package for estimating stochastic frontier models
#'
#' The \pkg{sfaR} package provides a set of tools (maximum likelihood - ML and
#' maximum simulated likelihood - MSL) for various specifications of stochastic
#' frontier analysis (SFA).
#'
#' Three categories of functions are available: \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfaselectioncross}},
#' which estimate different types of frontiers and offer eleven alternative 
#' optimization algorithms (i.e., "bfgs", "bhhh", "nr", "nm", "cg", "sann", 
#' "ucminf", "mla", "sr1", "sparse", "nlminb").
#'
#' @name sfaR-package
#'
#' @aliases sfaR-package sfaR
#'
#' @docType package
#' 
#' @section sfacross: \code{\link{sfacross}} estimates the basic stochastic 
#' frontier analysis (SFA) for cross-sectional or pooled data and allows for 
#' ten different distributions for the one-sided error term. These distributions 
#' include the exponential, the gamma, the generalized exponential,
#' the half normal, the lognormal, the truncated normal, the truncated skewed 
#' Laplace, the Rayleigh, the uniform, and the Weibull distributions. 
#' In the case of the gamma, lognormal, and Weibull distributions, maximum 
#' simulated likelihood (MSL) is used with the possibility of four specific 
#' distributions to construct the draws: halton, generalized halton, sobol and 
#' uniform. Heteroscedasticity in both error terms can be implemented, in 
#' addition to heterogeneity in the truncated mean parameter in the case of the 
#' truncated normal and lognormal distributions. In addition, in the case of the
#'  truncated normal distribution, the scaling property can be estimated.
#'  
#' @section sfalcmcross: \code{\link{sfalcmcross}} estimates latent class 
#' stochastic frontier models (LCM) for cross-sectional or pooled data. 
#' It accounts for technological heterogeneity by splitting the observations 
#' into a maximum number of five classes. The classification operates based on 
#' a logit functional form that can be specified using some covariates (namely, 
#' the separating variables allowing the separation of observations in several 
#' classes). Only the half normal distribution is available for the one-sided 
#' error term. Heteroscedasticity in both error terms is possible. The choice of 
#' the number of classes can be guided by several information criteria (i.e., 
#' AIC, BIC, or HQIC).
#'  
#' @section sfaselectioncross: \code{\link{sfaselectioncross}} estimates the 
#' frontier for cross-sectional or pooled data in the presence of sample 
#' selection. The model solves the selection bias due to the correlation 
#' between the two-sided error terms in both the selection and the frontier 
#' equations. The likelihood can be estimated using five different
#' possibilities: gauss-kronrod quadrature, adaptive integration over hypercubes 
#' (hcubature and pcubature), gauss-hermite quadrature, and
#'  maximum simulated likelihood. Only the half normal
#' distribution is available for the one-sided error term. Heteroscedasticity
#' in both error terms is possible.
#' 
#' @section Bugreport: Any bug or suggestion can be reported using the 
#' \code{sfaR} tracker facilities at: 
#' \url{https://github.com/hdakpo/sfaR/issues}
#'
#' @author K Hervé Dakpo, Yann Desjeux, Arne Henningsen and Laure Latruffe
#'
# @importFrom base standardGeneric
#' @importFrom stats coefficients dnorm lm model.frame
#' @importFrom stats model.matrix model.response nlminb
#' @importFrom stats pnorm qnorm delete.response fitted
#' @importFrom stats logLik residuals terms vcov formula
#' @importFrom stats integrate runif model.weights nobs
#' @importFrom stats na.pass printCoefmat pt qt na.omit
#' @importFrom stats pchisq qchisq uniroot complete.cases
#' @importFrom methods as new
#' @importFrom Formula as.Formula Formula
#' @importFrom randtoolbox get.primes sobol
#' @importFrom qrng ghalton
#' @import maxLik
#' @importFrom ucminf ucminf
#' @importFrom trustOptim trust.optim
#' @importFrom marqLevAlg mla
#' @importFrom nleqslv nleqslv
#' @importFrom fastGHQuad gaussHermiteData
#' @importFrom cubature hcubature pcubature
#' @importFrom sandwich bread estfun
# @importFrom calculus jacobian
#' @importFrom plm pdata.frame index
#' @importFrom texreg extract createTexreg
#' @importFrom mnorm pmnorm
NULL
