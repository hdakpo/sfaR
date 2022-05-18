################################################################################
#                                                                              #
# sfaR package doc                                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# sfaR package overview                                                        #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
#         -Sample selection correction                                         #
#         -Zero inefficiency stochastic frontier                               #
#         -Contaminated noise stochastic frontier                              #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' sfaR: A package for estimating stochastic frontier models
#'
#' The \pkg{sfaR} package provides a set of tools (maximum likelihood - ML and
#' maximum simulated likelihood - MSL) for various specifications of stochastic
#' frontier analysis (SFA).
#'
#' Five categories of important functions are available: \code{\link{cnsfcross}}, 
#' \code{\link{lcmcross}}, \code{\link{sfacross}}, 
#' \code{\link{sfaselectioncross}}, and \code{\link{zisfcross}}, which estimate 
#' different types of frontiers and offer eleven alternative optimization 
#' algorithms (i.e. "bfgs", "bhhh", "nr", "nm", "cg", "sann", "ucminf", "mla", 
#' "sr1", "sparse", "nlminb").
#'
#' @name sfaR-package
#'
#' @aliases sfaR-package sfaR
#'
#' @docType package
#' 
#' @section cnsfcross: \code{\link{cnsfcross}} the contaminated noise stochastic 
#' frontier (CNSF) allows for the possibility of outliers in the two-sided error
#' term. The model splits the observations into two groups with different 
#' variance for the two-sided error term. The probability of belonging to one 
#' of these classes is based on either a logit, probit, cauchit or cloglog 
#' functional form that can be specified using some covariates. As in the case 
#' of \code{\link{sfacross}} ten different distributions for the one-sided error 
#' term are implemented. These distributions include the exponential, the Gamma, 
#' the generalized exponential, the half normal, the lognormal, the truncated 
#' normal, the truncated skewed Laplace, the Rayleigh, the uniform and the 
#' Weibull distributions. In the case of the Gamma, lognormal and Weibull 
#' distributions, maximum simulated likelihood (MSL) is used with the 
#' possibility of four specific distributions to construct the draws: Halton, 
#' Generalized Halton, Sobol and uniform. Heteroscedasticity in both error 
#' terms can be implemented, in addition to heterogeneity in the truncated mean 
#' parameter in the case of the truncated normal and lognormal distributions. 
#' Two variants of the CNSF are available: one that allows common one-sided 
#' error term variance for both classes and the other one that allows different 
#' one-sided error term variances between the two aforementionned classes.
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
#' @section sfacross: \code{\link{sfacross}} estimates the basic stochastic 
#' frontier analysis (SFA) for cross-sectional data and allows for ten 
#' different distributions for the one-sided error term. These distributions 
#' include the exponential, the Gamma, the generalized exponential,
#' the half normal, the lognormal, the truncated normal, the truncated skewed 
#' Laplace, the Rayleigh, the uniform and the Weibull distributions. 
#' In the case of the Gamma, lognormal and Weibull distributions, maximum 
#' simulated likelihood (MSL) is used with the possibility of four specific 
#' distributions to construct the draws: Halton, Generalized Halton, Sobol and 
#' uniform. Heteroscedasticity in both error terms can be implemented, in 
#' addition to heterogeneity in the truncated mean parameter in the case of the 
#' truncated normal and lognormal distributions. In addition, in the case of the
#'  truncated normal distribution, the scaling property can be estimated.
#'   
#' @section sfaselectioncross: \code{\link{sfaselectioncross}} estimates the 
#' frontier for cross-sectional data in the presence of sample selection. 
#' The model solves the selection bias due to the correlation between the
#' two-sided errors terms in both the selection and the frontier equations. 
#' The likelihood can be estimated using five different
#' possibilities: Gauss-Kronrod quadrature, adaptive integration over hypercubes 
#' (hcubature and pcubature), Gauss-Hermite quadrature, and
#'  maximum simulated likelihood. Only the half normal
#' distribution is available for the one-sided error term. Heteroscedasticity
#' in both error terms is possible.
#' 
#' @section zisfcross: \code{\link{zisfcross}} relaxes the assumption that all 
#' firms are inefficient (zero inefficiency stochastic frontier - ZISF). 
#' The model splits the observations into two groups:
#' one efficient and the other inefficient. The probability of belonging to one 
#' of these classes is based on either a logit, probit, cauchit or cloglog 
#' functional form that can be specified using some covariates. As in the case 
#' of \code{\link{sfacross}} ten different distributions for the one-sided error 
#' term are implemented. These distributions include the exponential, the Gamma, 
#' the generalized exponential, the half normal, the lognormal, the truncated 
#' normal, the truncated skewed Laplace, the Rayleigh, the uniform and the 
#' Weibull distributions. In the case of the Gamma, lognormal and Weibull 
#' distributions, maximum simulated likelihood (MSL) is used with the 
#' possibility of four specific distributions to construct the draws: Halton, 
#' Generalized Halton, Sobol and uniform. Heteroscedasticity in both error 
#' terms can be implemented, in addition to heterogeneity in the truncated mean 
#' parameter in the case of the truncated normal and lognormal distributions. 
#' Two variants of the ZISF are available: one that allows common two-sided 
#' error term variance for both classes (inefficient vs efficient) and the other
#' one that allows different two-sided error term variances between the two
#' aforementionned classes.
#'
#' @section Bugreport: Any bug or suggestion can be reported using the 
#' \code{sfaR} tracker facilities at: \url{https://github.com/hdakpo/sfaR/issues}
#'
#' @author K Hervé Dakpo, Yann Desjeux, Laure Latruffe, and Arne Henningsen
#'
#' @importFrom stats coefficients dnorm lm model.frame
#' @importFrom stats model.matrix model.response nlminb
#' @importFrom stats pnorm qnorm delete.response fitted
#' @importFrom stats logLik residuals terms vcov formula
#' @importFrom stats integrate runif model.weights nobs
#' @importFrom stats na.pass printCoefmat pt qt na.omit
#' @importFrom stats pchisq qchisq uniroot
#' @importFrom methods as new
#' @importFrom Formula as.Formula Formula
#' @importFrom randtoolbox get.primes sobol
#' @importFrom qrng ghalton
#' @import maxLik
#' @importFrom ucminf ucminf
#' @importFrom trustOptim trust.optim
#' @importFrom marqLevAlg mla
#' @importFrom dplyr bind_cols as_tibble mutate select pull
#' @importFrom magrittr %>%
#' @importFrom nleqslv nleqslv
#' @importFrom katex math_to_rd
#' @importFrom fastGHQuad gaussHermiteData
#' @importFrom cubature hcubature pcubature
#' @importFrom sandwich bread estfun
#' @importFrom calculus jacobian
#' @importFrom plm pdata.frame index
NULL
