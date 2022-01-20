################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Model: Sample Selection Stochastic Frontier Analysis                         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Sample Selection in Stochastic frontier estimation using cross-section data
#'
#' @description
#' \code{\link{selectioncross}} is a symbolic formula based function for the
#' estimation of the stochastic frontier model in the presence of sample
#' selection. The model accommodates cross-section or pooled cross-section data. 
#' The model can be estimated using different quadrature approaches or 
#' maximum simulated likelihood (MSL). See Greene (2010).
#'
#' Only the half-normal distribution is possible for the one-sided error term.
#' Eleven optimization algorithms are available.
#'
#' The function also accounts for heteroscedasticity in both one-sided and
#' two-sided error terms, as in Reifschneider and Stevenson (1991), Caudill and
#' Ford (1993), Caudill \emph{et al.} (1995) and Hadri (1999).
#'
#' @aliases selectioncross print.selectioncross
#'
#' @param selectionF A symbolic (formula) description of the selection equation.
#' @param frontierF A symbolic (formula) description of the outcome (frontier) equation.
#' @param uhet A one-part formula to consider heteroscedasticity in the
#' one-sided error variance (see section \sQuote{Details}).
#' @param vhet A one-part formula to consider heteroscedasticity in the
#' two-sided error variance (see section \sQuote{Details}).
#' @param modelType Character string. Model used to solve the selection bias. Only the 
#' model discussed in Greene (2010) is available.
#' @param logDepVar Logical. Informs whether the dependent variable is logged
#' (\code{TRUE}) or not (\code{FALSE}). Default = \code{TRUE}.
#' @param data The data frame containing the data.
#' @param subset An optional vector specifying a subset of observations to be
#' used in the optimization process.
#' @param weights An optional vector of weights to be used for weighted log-likelihood.
#' Should be \code{NULL} or numeric vector with positive values. When \code{NULL}, 
#' a numeric vector of 1 is used.
#' @param wscale Logical. When \code{weights} is not \code{NULL}, a scaling transformation
#' is used such that the the \code{weights} sums to the sample size. Default \code{TRUE}.
#' When \code{FALSE} no scaling is used.
#' @param S If \code{S = 1} (default), a production (profit) frontier is
#' estimated: \eqn{\epsilon_i = v_i-u_i}. If \code{S = -1}, a cost frontier is
#' estimated: \eqn{\epsilon_i = v_i+u_i}.
#' @param udist Character string. Distribution specification for the one-sided
#' error term. Only the half normal distribution \code{'hnormal'} is currently
#' implemented.
#' @param start Numeric vector. Optional starting values for the maximum
#' likelihood (ML) estimation.
#' @param method Optimization algorithm used for the estimation.  Default =
#' \code{'bfgs'}. 11 algorithms are available: \itemize{ \item \code{'bfgs'},
#' for Broyden-Fletcher-Goldfarb-Shanno (see
#' \code{\link[maxLik:maxBFGS]{maxBFGS}}) \item \code{'bhhh'}, for
#' Berndt-Hall-Hall-Hausman (see \code{\link[maxLik:maxBHHH]{maxBHHH}}) \item
#' \code{'nr'}, for Newton-Raphson (see \code{\link[maxLik:maxNR]{maxNR}})
#' \item \code{'nm'}, for Nelder-Mead (see \code{\link[maxLik:maxNM]{maxNM}}) 
#' \item \code{'cg'}, for Conjugate Gradient (see \code{\link[maxLik:maxCG]{maxCG}})
#' \item \code{'sann'}, for Simulated Annealing (see \code{\link[maxLik:maxSANN]{maxSANN}})
#' \item \code{'ucminf'}, implements a quasi-Newton type with BFGS updating of the
#' inverse Hessian and soft line search with a trust region type monitoring of
#' the input to the line search algorithm (see \code{\link[ucminf:ucminf]{ucminf}})
#' \item \code{'mla'}, for general-purpose optimization based on
#' Marquardt-Levenberg algorithm (see \code{\link[marqLevAlg:mla]{mla}})
#' \item \code{'sr1'}, for Symmetric Rank 1 (see
#' \code{\link[trustOptim:trust.optim]{trust.optim}}) \item \code{'sparse'}, for trust
#' regions and sparse Hessian (see \code{\link[trustOptim:trust.optim]{trust.optim}}) \item
#' \code{'nlminb'}, for optimization using PORT routines (see
#' \code{\link[stats:nlminb]{nlminb}})}
#' @param hessianType Integer. If \code{1} (Default), analytic Hessian is
#' returned. If \code{2}, bhhh Hessian is estimated (\eqn{g'g}).
#' @param lType Specify the way the likelihood is estimated. Five possibilities are
#' available to the user: \code{kronrod} for Gauss-Kronrod quadrature 
#' (see \code{\link[stats:integrate]{integrate}}), \code{hcubature} and
#' \code{pcubature} for adaptive integration over hypercubes 
#' (see \code{\link[cubature:hcubature]{hcubature}} and 
#' \code{\link[cubature:pcubature]{pcubature}}), \code{ghermite} for Gauss-Hermite
#' quadrature (see \code{\link[fastGHQuad:gaussHermiteData]{gaussHermiteData}}), and
#' \code{msl} for maximum simulated likelihood. Default \code{kronrod}.
#' @param Nsub Integer. Number of subdivisions used for quadrature approaches. 
#' Default \code{Nsub = 100}.
#' @param uBound Numeric. Upper bound for the inefficiency component when solving
#' integrals using quadrature approaches except Gauss-Hermite for which the upper
#' bound is automatically infinite (\code{Inf}). Default \code{uBound = Inf}.
#' @param simType Character string. If \code{simType = 'halton'} (Default),
#' Halton draws are used for maximum simulated likelihood (MSL). If
#' \code{simType = 'ghalton'}, Generalized-Halton draws are used for MSL. If
#' \code{simType = 'sobol'}, Sobol draws are used for MSL. If \code{simType =
#' 'uniform'}, uniform draws are used for MSL. (see section \sQuote{Details}).
#' @param Nsim Number of draws for MSL.
#' @param prime Prime number considered for Halton and Generalized-Halton
#' draws. Default = \code{2}.
#' @param burn Number of the first observations discarded in the case of Halton
#' draws. Default = \code{10}.
#' @param antithetics Logical. Default = \code{FALSE}. If \code{TRUE},
#' antithetics counterpart of the uniform draws is computed. (see section
#' \sQuote{Details}).
#' @param seed Numeric. Seed for the random draws.
#' @param itermax Maximum number of iterations allowed for optimization.
#' Default = \code{2000}.
#' @param printInfo Logical. Print information during optimization. Default =
#' \code{FALSE}.
#' @param intol Numeric. Integration tolerance for quadrature approaches 
#' (\code{kronrod, hcubature, pcubature}). 
#' @param tol Numeric. Convergence tolerance. Default = \code{1e-12}.
#' @param gradtol Numeric. Convergence tolerance for gradient. Default =
#' \code{1e-06}.
#' @param stepmax Numeric. Step max for \code{ucminf} algorithm. Default =
#' \code{0.1}.
#' @param qac Character. Quadratic Approximation Correction for \code{'bhhh'}
#' and \code{'nr'} algorithms. If \code{'stephalving'}, the step length is
#' decreased but the direction is kept. If \code{'marquardt'} (default), the
#' step length is decreased while also moving closer to the pure gradient
#' direction. See \code{\link[maxLik:maxBHHH]{maxBHHH}} and
#' \code{\link[maxLik:maxNR]{maxNR}}.
#' @param x an object of class selectioncross (returned by the function \code{\link{selectioncross}}).
#' @param ... additional arguments of frontier are passed to selectioncross; 
#' additional arguments of the print, bread, estfun, nobs methods are currently ignored.
#' 
#' @details 
#' The current model is an extension of Heckman(1976, 1979) sample selection model to
#' nonlinear models particularly stochastic frontier model. The model has first been discussed in 
#' Greene (2010), and an application can be found in Dakpo et al. (2021). Practically, we have:
#' 
#' \deqn{
#' y_{1i} =  \left\{ \begin{array}{ll}
#' 1 & \mbox{if} \quad y_{1i}^* > 0  \\
#' 0 & \mbox{if} \quad y_{1i}^* \leq 0 \\
#' \end{array}
#' \right.
#' }
#' 
#' where 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('y_{1i}^*=\\\mathbf{Z_{si}^{\\\prime}} \\\bm{\\\gamma} + w_i, \\\quad w_i \\\sim \\\mathcal{N}(0, 1)')
#' } 
#' 
#' and
#' 
#' \deqn{
#' y_{2i} =  \left\{ \begin{array}{ll}
#' y_{2i}^* & \mbox{if} \quad y_{1i}^* > 0  \\
#' NA & \mbox{if} \quad y_{1i}^* \leq 0 \\
#' \end{array}
#' \right.
#' }
#'
#' where 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('y_{2i}^*=\\\mathbf{x_{i}^{\\\prime}} \\\bm{\\\beta} + v_i - Su_i, \\\quad v_i = \\\sigma_vV_i \\\quad \\\wedge \\\quad V_i \\\sim \\\mathcal{N}(0, 1), \\\quad u_i = \\\sigma_u|U_i| \\\quad \\\wedge \\\quad U_i \\\sim \\\mathcal{N}(0, 1)')
#' }
#' 
#' \eqn{y_{1i}} describes the selection equation while \eqn{y_{2i}} represents the frontier equation. The selection
#' bias emerges from the correlation between the two symmetric random components \eqn{v_i} and \eqn{w_i} with
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('(v_i, w_i) \\\sim \\\mathcal{N}_2\\\\left\\\\lbrack(0,0), (1, \\\rho \\\sigma_v, \\\sigma_v^2) \\\\right\\\\rbrack')
#' }
#' 
#' Conditionaly on \eqn{|U_i|}, the probability associated to each observation is
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('Pr \\\\left\\\\lbrack y_{1i}^* \\\\leq 0 \\\\right\\\\rbrack^{1-y_{1i}} \\\cdot \\\\left\\\\lbrace f(y_{2i}|y_{1i}^* > 0) \\\times Pr\\\\left\\\\lbrack y_{1i}^* > 0\\\\right\\\\rbrack \\\\right\\\\rbrace^{y_{1i}}')
#' }
#' 
#' Using the conditional probability formula
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('P\\\\left(A\\\cap B\\\\right) = P(A) \\\cdot P(B|A) = P(B) \\\cdot P(A|B)')
#' }
#' 
#' Therefore
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('f(y_{2i}|y_{1i}^* \\\geq 0) \\\cdot Pr\\\\left\\\\lbrack y_{1i}^* \\\geq 0\\\\right\\\\rbrack = f(y_{2i}) \\\cdot Pr(y_{1i}^* \\\geq 0|y_{2i})')
#' }
#' 
#' Using the properties of a bivariate normal distribution, we have
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('y_{i1}^* | y_{i2} \\\sim N\\\\left(\\\mathbf{Z_{si}^{\\\prime}} \\\bm{\\\gamma}+\\\frac{\\\rho}{\\\sigma_v}v_i, 1-\\\rho^2\\\\right)')
#' }
#' 
#' Hence conditionally on \eqn{|U_i|}, we have
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('f(y_{2i}|y_{1i}^* \\\geq 0) \\\cdot Pr\\\\left\\\\lbrack y_{1i}^* \\\geq 0\\\\right\\\\rbrack = \\\frac{1}{\\\sigma_v}\\\phi\\\\left(\\\frac{v_i}{\\\sigma_v}\\\\right)\\\Phi\\\\left(\\\frac{\\\mathbf{Z_{si}^{\\\prime}} \\\bm{\\\gamma}+\\\frac{\\\rho}{\\\sigma_v}v_i}{\\\sqrt{1-\\\rho^2}}\\\\right)')
#' }
#' 
#' The conditional likelihood is equal to
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('L_i\\\big||U_i| = \\\Phi(-\\\mathbf{Z_{si}^{\\\prime}} \\\bm{\\\gamma})^{1-y_{1i}} \\\times \\\\left\\\\lbrace \\\frac{1}{\\\sigma_v}\\\phi\\\\left(\\\frac{y_{2i}^*-\\\mathbf{x_{i}^{\\\prime}} \\\bm{\\\beta} + S\\\sigma_u|U_i|}{\\\sigma_v}\\\\right)\\\Phi\\\\left(\\\frac{\\\mathbf{Z_{si}^{\\\prime}} \\\bm{\\\gamma}+\\\frac{\\\rho}{\\\sigma_v}\\\\left(y_{2i}-\\\mathbf{x_{i}^{\\\prime}} \\\bm{\\\beta} + S\\\sigma_u|U_i|\\\\right)}{\\\sqrt{1-\\\rho^2}}\\\\right) \\\\right\\\\rbrace ^{y_{1i}}') 
#' }
#' 
#' Since the non-selected observations do not bring any additional information, thus 
#' the conditional likelihood to be considered is
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('L_i\\\big||U_i| = \\\frac{1}{\\\sigma_v}\\\phi\\\\left(\\\frac{y_{2i}^*-\\\mathbf{x_{i}^{\\\prime}} \\\bm{\\\beta} + S\\\sigma_u|U_i|}{\\\sigma_v}\\\\right) \\\Phi\\\\left(\\\frac{\\\mathbf{Z_{si}^{\\\prime}} \\\bm{\\\gamma}+\\\frac{\\\rho}{\\\sigma_v}\\\\left(y_{2i}-\\\mathbf{x_{i}^{\\\prime}} \\\bm{\\\beta} + S\\\sigma_u|U_i|\\\\right)}{\\\sqrt{1-\\\rho^2}}\\\\right) ') 
#' }
#' 
#' The unconditional likelihood is obtained by integrating \eqn{|U_i|} out of the conditional likelihood. Thus
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('L_i\\ = \\\int_{|U_i|} \\\frac{1}{\\\sigma_v}\\\phi\\\\left(\\\frac{y_{2i}^*-\\\mathbf{x_{i}^{\\\prime}} \\\bm{\\\beta} + S\\\sigma_u|U_i|}{\\\sigma_v}\\\\right) \\\Phi\\\\left(\\\frac{\\\mathbf{Z_{si}^{\\\prime}} \\\bm{\\\gamma}+ \\\frac{\\\rho}{\\\sigma_v}\\\\left(y_{2i}-\\\mathbf{x_{i}^{\\\prime}} \\\bm{\\\beta} + S\\\sigma_u|U_i|\\\\right)}{\\\sqrt{1-\\\rho^2}}\\\\right) ') 
#' }
#' 
#' To simplifiy the estimation, the likelihood can be estimated using a two-step approach.
#' In the first step, the probit model can be run an estimate of \eqn{\gamma} can be obtained.
#' Then, in the second step, the following model is estimated:
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('L_i\\ = \\\int_{|U_i|} \\\frac{1}{\\\sigma_v}\\\phi\\\\left(\\\frac{y_{2i}^*-\\\mathbf{x_{i}^{\\\prime}} \\\bm{\\\beta} + S\\\sigma_u|U_i|}{\\\sigma_v}\\\\right) \\\Phi\\\\left(\\\frac{a_i \\\bm{\\\gamma}+ \\\frac{\\\rho}{\\\sigma_v}\\\\left(y_{2i}-\\\mathbf{x_{i}^{\\\prime}} \\\bm{\\\beta} + S\\\sigma_u|U_i|\\\\right)}{\\\sqrt{1-\\\rho^2}}\\\\right) ') 
#' }
#' 
#' where \eqn{a_i = \mathbf{Z_{si}^{\prime}} \hat{\bm{\gamma}}}. This likelihood can be estimated using 
#' five different approaches: Gauss-Kronrod quadrature, adaptive integration over hypercubes 
#' (hcubature and pcubature), Gauss-Hermite quadrature, and
#'  maximum simulated likelihood. We also use the BHHH estimator to obtain 
#'  the asymptotic standard errors for the parameter estimators.
#'  
#'  \code{selectioncross} allows for the maximization of weighted log-likelihood.
#' When option \code{weights} is specified and \code{wscale = TRUE}, the weights
#' is scaled as 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('new_{weights} = sample_{size} \\\times \\\frac{old_{weights}}{\\\sum(old_{weights})}')
#' }
#' 
#' For difficult problems, non-gradient methods (e.g. \code{nm} or \code{sann}) can be 
#' used to warm start the optimization and zoom in the neighborhood of the 
#' solution. Then a gradient-based methods is recommanded in the second step. In the case
#' of \code{sann}, we recommand to significantly increase the iteration limit 
#' (e.g. \code{itermax = 20000}). The Conjugate Gradient (\code{cg}) can also be used
#' in the first stage.
#' 
#' A set of extractor functions for fitted model objects is available for objects of class
#' \code{'selectioncross'} including methods to the generic functions \code{\link[=print.selectioncross]{print}},
#' \code{\link[=summary.selectioncross]{summary}}, \code{\link[=coef.selectioncross]{coef}}, 
#' \code{\link[=fitted.selectioncross]{fitted}}, \code{\link[=logLik.selectioncross]{logLik}}, 
#' \code{\link[=residuals.selectioncross]{residuals}}, \code{\link[=vcov.selectioncross]{vcov}}, 
#' \code{\link[=efficiencies.selectioncross]{efficiencies}}, \code{\link[=ic.selectioncross]{ic}}, 
#' \code{\link[=marginal.selectioncross]{marginal}}, 
#' \code{\link[=estfun.selectioncross]{estfun}} and 
#' \code{\link[=bread.selectioncross]{bread}} (from the \CRANpkg{sandwich} package), 
#' [lmtest::coeftest()] (from the \CRANpkg{lmtest} package).
#' 
#' @return \code{\link{selectioncross}} returns a list of class \code{'selectioncross'}
#' containing the following elements:
#'
#' \item{call}{The matched call.}
#'
#' \item{selectionF}{The selection equation formula.}
#' 
#' \item{frontierF}{The frontier equation formula.}
#'
#' \item{S}{The argument \code{'S'}. See the section \sQuote{Arguments}.}
#'
#' \item{typeSfa}{Character string. 'Stochastic Production/Profit Frontier, e =
#' v - u' when \code{S = 1} and 'Stochastic Cost Frontier, e = v + u' when
#' \code{S = -1}.}
#'
#' \item{Nobs}{Number of observations used for optimization.}
#'
#' \item{nXvar}{Number of explanatory variables in the production or cost
#' frontier.}
#'
#' \item{logDepVar}{The argument \code{'logDepVar'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{nuZUvar}{Number of variables explaining heteroscedasticity in the
#' one-sided error term.}
#'
#' \item{nvZVvar}{Number of variables explaining heteroscedasticity in the
#' two-sided error term.}
#'
#' \item{nParm}{Total number of parameters estimated.}
#'
#' \item{udist}{The argument \code{'udist'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{startVal}{Numeric vector. Starting value for M(S)L estimation.}
#'
#' \item{dataTable}{A data frame (tibble format) containing information on data
#' used for optimization along with residuals and fitted values of the OLS and
#' M(S)L estimations, and the individual observation log-likelihood. When \code{weights}
#' is specified an additional variable is also provided in \code{dataTable}.}
#' 
#' \item{lpmObj}{Linear probability model used for initializing the first step
#' probit model}
#' 
#' \item{probitObj}{Probit model. Object of class \code{'maxLik'} and \code{'maxim'}.}
#'
#' \item{ols2stepParam}{Numeric vector. OLS second step estimates for 
#' selection correction. Inverse Mills Ratio is introduced as an additional 
#' explanatory variable.}
#'
#' \item{ols2stepStder}{Numeric vector. Standard errors of OLS second step estimates.}
#'
#' \item{ols2stepSigmasq}{Numeric. Estimated variance of OLS second step random error.}
#'
#' \item{ols2stepLoglik}{Numeric. Log-likelihood value of OLS second step estimation.}
#'
#' \item{ols2stepSkew}{Numeric. Skewness of the residuals of the OLS second step estimation.}
#'
#' \item{ols2stepM3Okay}{Logical. Indicating whether the residuals of the OLS
#' second step estimation have the expected skewness.}
#'
#' \item{CoelliM3Test}{Coelli's test for OLS residuals skewness. (See Coelli,
#' 1995).}
#'
#' \item{AgostinoTest}{D'Agostino's test for OLS residuals skewness. (See
#' D'Agostino and Pearson, 1973).}
#' 
#' \item{isWeights}{Logical. If \code{TRUE} weighted log-likelihood is
#' maximized.}
#' 
#' \item{lType}{Type of likelihood estimated. See the section \sQuote{Arguments}.}
#'
#' \item{optType}{Optimization algorithm used.}
#'
#' \item{nIter}{Number of iterations of the ML estimation.}
#'
#' \item{optStatus}{Optimization algorithm termination message.}
#'
#' \item{startLoglik}{Log-likelihood at the starting values.}
#'
#' \item{mlLoglik}{Log-likelihood value of the M(S)L estimation.}
#'
#' \item{mlParam}{Parameters obtained from M(S)L estimation.}
#'
#' \item{gradient}{Each variable gradient of the M(S)L estimation.}
#'
#' \item{gradL_OBS}{Matrix. Each variable individual observation gradient of
#' the M(S)L estimation.}
#'
#' \item{gradientNorm}{Gradient norm of the M(S)L estimation.}
#'
#' \item{invHessian}{Covariance matrix of the parameters obtained from the
#' M(S)L estimation.}
#'
#' \item{hessianType}{The argument \code{'hessianType'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{mlDate}{Date and time of the estimated model.}
#'
#' \item{simDist}{The argument \code{'simDist'}, only if \code{lType =
#' 'msl'}. See the section \sQuote{Arguments}.}
#'
#' \item{Nsim}{The argument \code{'Nsim'}, only if \code{lType = 'msl'}.
#' See the section \sQuote{Arguments}.}
#'
#' \item{FiMat}{Matrix of random draws used for MSL, only if \code{lType =
#' 'msl'}.}
#' 
#' \item{gHermiteData}{List. Gauss-Hermite quadrature rule as provided by 
#' \code{\link[fastGHQuad:gaussHermiteData]{gaussHermiteData}}. Only if \code{lType = 
#' 'ghermite'}.}
#' 
#' \item{Nsub}{Number of subdivisions used for quadrature approaches.}
#' 
#' \item{uBound}{Upper bound for the inefficiency component when solving
#' integrals using quadrature approaches except Gauss-Hermite for which the upper
#' bound is automatically infinite (\code{Inf}).}
#' 
#' \item{intol}{Integration tolerance for quadrature approaches except Gauss-Hermite.}
#'
#' @note For the Halton draws, the code is adapted from the \pkg{mlogit}
#' package.
#'
# @author K Hervé Dakpo, Yann Desjeux, and Laure Latruffe
#'
#' @seealso \code{\link[=print.selectioncross]{print}} for printing \code{selectioncross} object.
#' 
#' \code{\link[=summary.selectioncross]{summary}} for creating and printing
#' summary results.
#'
#' \code{\link[=coef.selectioncross]{coef}} for extracting coefficients of the
#' estimation.
#'
#' \code{\link[=efficiencies.selectioncross]{efficiencies}} for computing
#' (in-)efficiency estimates.
#'
#' \code{\link[=fitted.selectioncross]{fitted}} for extracting the fitted frontier
#' values.
#'
#' \code{\link[=ic.selectioncross]{ic}} for extracting information criteria.
#'
#' \code{\link[=logLik.selectioncross]{logLik}} for extracting log-likelihood
#' value(s) of the estimation.
#'
#' \code{\link[=marginal.selectioncross]{marginal}} for computing marginal effects of
#' inefficiency drivers.
#'
#' \code{\link[=residuals.selectioncross]{residuals}} for extracting residuals of the
#' estimation.
#'
#' \code{\link[=vcov.selectioncross]{vcov}} for computing the variance-covariance
#' matrix of the coefficients.
#' 
#' \code{\link[=bread.selectioncross]{bread}} for bread for sandwich estimator.
#' 
#' \code{\link[=estfun.selectioncross]{estfun}} for gradient extraction for each 
#' observation.
#'
#' @references Caudill, S. B., and Ford, J. M. 1993. Biases in frontier estimation due to
#' heteroscedasticity. \emph{Economics Letters}, \bold{41}(1), 17--20.
#'
#' Caudill, S. B., Ford, J. M., and Gropper, D. M. 1995. Frontier estimation
#' and firm-specific inefficiency measures in the presence of
#' heteroscedasticity. \emph{Journal of Business & Economic Statistics},
#' \bold{13}(1), 105--111.
#' 
#' Coelli, T. 1995. Estimators and hypothesis tests for a stochastic frontier
#' function - a Monte-Carlo analysis. \emph{Journal of Productivity Analysis},
#' \bold{6}:247--268.
#'
#' D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
#' Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
#' \emph{Biometrika}, \bold{60}:613--622.
#' 
#' Dakpo, K. H., Latruffe, L., Desjeux, Y., Jeanneaux, P., 2021. 
#' Modeling heterogeneous technologies in the presence of sample selection: 
#' The case of dairy farms and the adoption of agri-environmental schemes in France.
#' \emph{Agricultural Economics}. n/a.
#' 
#' Greene, W., 2010. A stochastic frontier model with correction 
#' for sample selection. \emph{Journal of Productivity Analysis}. \bold{34}, 15--24.
#' 
#' Hadri, K. 1999. Estimation of a doubly heteroscedastic stochastic frontier
#' cost function. \emph{Journal of Business & Economic Statistics},
#' \bold{17}(3), 359--363.
#' 
#' Heckman, J., 1976. Discrete, qualitative and limited dependent variables. 
#' \emph{Ann Econ Soc Meas.} \bold{4}, 475--492.
#' 
#' Heckman, J., 1979. Sample Selection Bias as a Specification Error.
#' \emph{Econometrica}. \bold{47}, 153--161.
#' 
#' Reifschneider, D., and Stevenson, R. 1991. Systematic departures from the
#' frontier: A framework for the analysis of firm inefficiency.
#' \emph{International Economic Review}, \bold{32}(3), 715--723.
#' 
#' @keywords models optimize cross-section likelihood
#'
#' @examples
#' 
#' ## Simulated example
#' 
#' N <- 2000  # sample size
#' set.seed(12345)
#' z1 <- rnorm(N)
#' z2 <- rnorm(N)
#' v1 <- rnorm(N)
#' v2 <- rnorm(N)
#' e1 <- v1
#' e2 <- 0.7071 * (v1 + v2)
#' ds <- z1 + z2 + e1
#' d <- ifelse(ds > 0, 1, 0)
#' u <- abs(rnorm(N))
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' y <- x1 + x2 + e2 - u
#' data <- cbind(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2, d = d)
#' 
#' ## Estimation using quadrature (Gauss-Kronrod)
#' 
#' selecRes1 <- selectioncross(selectionF = d ~ z1 + z2, frontierF = y ~ x1 + x2, 
#' modelType = 'greene10', method = 'bfgs',
#' logDepVar = TRUE, data = as.data.frame(data),
#' S = 1L, udist = 'hnormal', lType = 'kronrod', Nsub = 100, uBound = Inf,
#' simType = 'halton', Nsim = 300, prime = 2L, burn = 10, antithetics = FALSE,
#' seed = 12345, itermax = 2000, printInfo = TRUE)
#' 
#' summary(selecRes1)
#' 
#' ## Estimation using maximum simulated likelihood
#' 
#' selecRes2 <- selectioncross(selectionF = d ~ z1 + z2, frontierF = y ~ x1 + x2, 
#' modelType = 'greene10', method = 'bfgs',
#' logDepVar = TRUE, data = as.data.frame(data),
#' S = 1L, udist = 'hnormal', lType = 'msl', Nsub = 100, uBound = Inf,
#' simType = 'halton', Nsim = 300, prime = 2L, burn = 10, antithetics = FALSE,
#' seed = 12345, itermax = 2000, printInfo = TRUE)
#' 
#' summary(selecRes2)
#' 
#' @export
selectioncross <- function(selectionF, frontierF, uhet, vhet,
  modelType = "greene10", logDepVar = TRUE, data, subset, weights,
  wscale = TRUE, S = 1L, udist = "hnormal", start = NULL, method = "bfgs",
  hessianType = 2L, lType = "kronrod", Nsub = 100, uBound = Inf,
  simType = "halton", Nsim = 300, prime = 2L, burn = 10, antithetics = FALSE,
  seed = 12345, itermax = 2000, printInfo = FALSE, intol = 1e-06,
  tol = 1e-12, gradtol = 1e-06, stepmax = 0.1, qac = "marquardt") {
  # u distribution check -------
  udist <- tolower(udist)
  if (udist != "hnormal") {
    stop("Currently selection model only handles half-normal distribution ... ",
      call. = FALSE)
  }
  # selection model check -------
  modelType <- tolower(modelType)
  if (modelType != "greene10") {
    stop("Currently only Greene (2010) selection model is available ... ",
      call. = FALSE)
  }
  # method to solve selection model -------
  lType <- tolower(lType)
  if (!(lType %in% c("kronrod", "hcubature", "pcubature", "ghermite",
    "msl"))) {
    stop("Unknown inefficiency distribution: ", paste(lType),
      call. = FALSE)
  }
  # Formula manipulation -------
  if (length(Formula(selectionF))[2] != 1) {
    stop("argument 'selectionF' must have one RHS part",
      call. = FALSE)
  }
  if (length(Formula(frontierF))[2] != 1) {
    stop("argument 'frontierF' must have one RHS part", call. = FALSE)
  }
  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  ## Selection formula -------
  m_S <- match(c("selectionF", "data", "subset", "weights"),
    names(mc), nomatch = 0L)
  mc_S <- mc[c(1L, m_S)]
  mc_S$drop.unused.levels <- TRUE
  selectionF <- interCheckSelection(formula = selectionF)
  ## Generate required datasets -------
  names(mc_S)[2] <- "formula"
  mc_S$na.action <- na.omit
  mc_S[[1L]] <- quote(stats::model.frame)
  mc_S <- eval(mc_S, parent.frame())
  Yvar_S <- model.response(mc_S)
  selecLevels <- levels(as.factor(Yvar_S))
  if (length(selecLevels) != 2) {
    stop("the dependent variable of the 'selection' equation has to contain",
      "exactly two levels (e.g. FALSE and TRUE)")
  }
  Yvar_S <- as.integer(Yvar_S == selecLevels[2])
  mtX_S <- terms(selectionF, data = data, rhs = 1)
  Xvar_S <- model.matrix(mtX_S, mc_S)
  nXvar_S <- ncol(Xvar_S)
  N_S <- nrow(Xvar_S)
  if (N_S == 0L) {
    stop("0 (non-NA) cases in selection model", call. = FALSE)
  }
  wProbit <- as.vector(model.weights(mc_S))
  if (length(wscale) != 1 || !is.logical(wscale[1])) {
    stop("argument 'wscale' must be a single logical value",
      call. = FALSE)
  }
  if (!is.null(wProbit)) {
    if (!is.numeric(wProbit)) {
      stop("'weights' must be a numeric vector", call. = FALSE)
    } else {
      if (any(wProbit < 0 | is.na(wProbit)))
        stop("missing or negative weights not allowed",
          call. = FALSE)
    }
    if (wscale) {
      wProbit <- wProbit/sum(wProbit) * N_S
    }
  } else {
    wProbit <- rep(1, N_S)
  }
  ## frontier formula -------
  m_F <- match(c("frontierF", "data", "subset", "weights"),
    names(mc), nomatch = 0L)
  mc_F <- mc[c(1L, m_F)]
  mc_F$drop.unused.levels <- TRUE
  frontierF <- interCheckMain(formula = frontierF)
  ## Generate required datasets -------
  names(mc_F)[2] <- "formula"
  mc_F$na.action <- na.omit
  mc_F[[1L]] <- quote(stats::model.frame)
  mc_F <- eval(mc_F, parent.frame())
  Yvar_F <- model.response(mc_F, "numeric")[Yvar_S == 1]
  mtX_F <- terms(frontierF, data = data, rhs = 1)
  Xvar_F <- model.matrix(mtX_F, mc_F)[Yvar_S == 1, , drop = FALSE]
  nXvar_F <- ncol(Xvar_F)
  N_F <- nrow(Xvar_F)
  if (N_F == 0L) {
    stop("0 (non-NA) cases in frontier model", call. = FALSE)
  }
  wHvar <- wProbit[Yvar_S == 1]
  if (!missing(uhet)) {
    uhet <- clhsCheck_u(formula = uhet, scaling = FALSE)
  } else {
    uhet <- ~1
  }
  if (!missing(vhet)) {
    vhet <- clhsCheck_v(formula = vhet)
  } else {
    vhet <- ~1
  }
  formula <- formDist_selectioncross(udist = udist, formula = frontierF,
    uhet = uhet, vhet = vhet)
  mtuH <- delete.response(terms(formula, data = data, rhs = 2))
  uHvar <- model.matrix(mtuH, mc_F)[Yvar_S == 1, , drop = FALSE]
  nuZUvar <- ncol(uHvar)
  mtvH <- delete.response(terms(formula, data = data, rhs = 3))
  vHvar <- model.matrix(mtvH, mc_F)[Yvar_S == 1, , drop = FALSE]
  nvZVvar <- ncol(vHvar)
  # Check other supplied options -------
  if (length(S) != 1 || !(S %in% c(-1L, 1L))) {
    stop("argument 'S' must equal either 1 or -1: 1 for production or profit frontier
   and -1 for cost frontier",
      call. = FALSE)
  }
  typeSfa <- if (S == 1L) {
    "Stochastic Production/Profit Frontier, e = v - u"
  } else {
    "Stochastic Cost Frontier, e = v + u"
  }
  if (length(logDepVar) != 1 || !is.logical(logDepVar[1])) {
    stop("argument 'logDepVar' must be a single logical value",
      call. = FALSE)
  }
  # Number of parameters -------
  nParm <- nXvar_F + nuZUvar + nvZVvar + 1
  # Checking starting values when provided -------
  if (!is.null(start)) {
    if (length(start) != nParm) {
      stop("Wrong number of initial values: model has ",
        nParm, " parameters", call. = FALSE)
    }
  }
  if (nParm > N_F) {
    stop("Model has more parameters than observations", call. = FALSE)
  }
  # Check algorithms -------
  method <- tolower(method)
  if (!(method %in% c("ucminf", "bfgs", "bhhh", "nr", "nm",
    "sr1", "mla", "sparse", "nlminb"))) {
    stop("Unknown or non-available optimization algorithm: ",
      paste(method), call. = FALSE)
  }
  # Check hessian type
  if (length(hessianType) != 1 || !(hessianType %in% c(1L,
    2L))) {
    stop("argument 'hessianType' must equal either 1 or 2",
      call. = FALSE)
  }
  # Draws for MSL -------
  if (lType == "msl") {
    if (!(simType %in% c("halton", "ghalton", "sobol", "uniform"))) {
      stop("Unknown or non-available random draws method",
        call. = FALSE)
    }
    if (!is.numeric(Nsim) || length(Nsim) != 1) {
      stop("argument 'Nsim' must be a single numeric scalar",
        call. = FALSE)
    }
    if (!is.numeric(burn) || length(burn) != 1) {
      stop("argument 'burn' must be a single numeric scalar",
        call. = FALSE)
    }
    if (!is_prime(prime)) {
      stop("argument 'prime' must be a single prime number",
        call. = FALSE)
    }
    if (length(antithetics) != 1 || !is.logical(antithetics[1])) {
      stop("argument 'antithetics' must be a single logical value",
        call. = FALSE)
    }
    if (antithetics && (Nsim%%2) != 0) {
      Nsim <- Nsim + 1
    }
    simDist <- if (simType == "halton") {
      "Halton"
    } else {
      if (simType == "ghalton") {
        "Generalized Halton"
      } else {
        if (simType == "sobol") {
          "Sobol"
        } else {
          if (simType == "uniform") {
          "Uniform"
          }
        }
      }
    }
    cat("Initialization of", Nsim, simDist, "draws per observation ...\n")
    FiMat <- drawMat(N = N_F, Nsim = Nsim, simType = simType,
      prime = prime, burn = burn + 1, antithetics = antithetics,
      seed = seed)
  }
  # Other optimization options -------
  if (!is.numeric(Nsub) || length(Nsub) != 1) {
    stop("argument 'Nsub' must be a single numeric scalar",
      call. = FALSE)
  }
  if (Nsub != round(Nsub)) {
    stop("argument 'Nsub' must be an integer", call. = FALSE)
  }
  if (Nsub <= 0) {
    stop("argument 'Nsub' must be positive", call. = FALSE)
  }
  Nsub <- as.integer(Nsub)
  uBound <- abs(uBound)
  if (!is.numeric(itermax) || length(itermax) != 1) {
    stop("argument 'itermax' must be a single numeric scalar",
      call. = FALSE)
  }
  if (itermax != round(itermax)) {
    stop("argument 'itermax' must be an integer", call. = FALSE)
  }
  if (itermax <= 0) {
    stop("argument 'itermax' must be positive", call. = FALSE)
  }
  itermax <- as.integer(itermax)
  if (length(printInfo) != 1 || !is.logical(printInfo[1])) {
    stop("argument 'printInfo' must be a single logical value",
      call. = FALSE)
  }
  if (!is.numeric(intol) || length(intol) != 1) {
    stop("argument 'intol' must be numeric", call. = FALSE)
  }
  if (intol < 0) {
    stop("argument 'intol' must be non-negative", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1) {
    stop("argument 'tol' must be numeric", call. = FALSE)
  }
  if (tol < 0) {
    stop("argument 'tol' must be non-negative", call. = FALSE)
  }
  if (!is.numeric(gradtol) || length(gradtol) != 1) {
    stop("argument 'gradtol' must be numeric", call. = FALSE)
  }
  if (gradtol < 0) {
    stop("argument 'gradtol' must be non-negative", call. = FALSE)
  }
  if (!is.numeric(stepmax) || length(stepmax) != 1) {
    stop("argument 'stepmax' must be numeric", call. = FALSE)
  }
  if (stepmax < 0) {
    stop("argument 'stepmax' must be non-negative", call. = FALSE)
  }
  if (!(qac %in% c("marquardt", "stephalving"))) {
    stop("argument 'qac' must be either 'marquardt' or 'stephalving'",
      call. = FALSE)
  }
  # Nodes for Gauss Hermite quadrature -------
  if (lType == "ghermite") {
    gH <- gaussHermiteData(n = Nsub * 2)
  }
  # Step 1: 2step Heckman correction -------

  ## OLS with linear probability model -------
  lpm <- if (colnames(Xvar_S)[1] == "(Intercept)") {
    if (dim(Xvar_S)[2] == 1) {
      lm(Yvar_S ~ 1)
    } else {
      lm(Yvar_S ~ ., data = as.data.frame(Xvar_S[, -1,
        drop = FALSE]), weights = wProbit)
    }
  } else {
    lm(Yvar_S ~ -1 + ., data = as.data.frame(Xvar_S), weights = wProbit)
  }
  if (any(is.na(lpm$coefficients))) {
    stop("at least one of the 'selectionF' coefficients is NA: ",
      paste(colnames(Xvar_S)[is.na(lpm$coefficients)],
        collapse = ", "), "This may be due to a singular matrix
   due to potential perfect multicollinearity",
      call. = FALSE)
  }
  ## Solution for probit model -------
  cat("First step probit model...\n")
  probRes <- maxLik(logLik = probit_likelihood, grad = probit_gradient,
    start = coefficients(lpm), method = "bfgs", control = list(iterlim = itermax,
      printLevel = printInfo, reltol = tol, tol = tol,
      qac = qac), Yvar = Yvar_S, Xvar = Xvar_S, wHvar = wProbit)
  probitParam <- probRes$estimate
  probitLoglik <- logLik(probRes)[1]
  ## heckit model using IMR and OLS -------
  if (inherits(data, "plm.dim")) {
    dataTable <- data[Yvar_S == 1, 1:2]
  } else {
    dataTable <- data.frame(IdObs = c(1:N_F))
  }
  dataTable <- as_tibble(cbind(dataTable, data[Yvar_S == 1,
    c(all.vars(terms(selectionF)), all.vars(terms(frontierF)))],
    weights = wHvar))
  dataTable$PROBIT_PREDICTIONS <- as.numeric(crossprod(matrix(probitParam),
    t(Xvar_S[Yvar_S == 1, , drop = FALSE])))
  dataTable$IMR <- dnorm(dataTable[["PROBIT_PREDICTIONS"]])/pnorm(dataTable[["PROBIT_PREDICTIONS"]])
  ols2step <- if (colnames(Xvar_F)[1] == "(Intercept)") {
    if (dim(Xvar_F)[2] == 1) {
      lm(Yvar_F ~ 1)
    } else {
      lm(Yvar_F ~ ., data = as.data.frame(cbind(Xvar_F[,
        -1, drop = FALSE], IMR = dataTable$IMR)), weights = wHvar)
    }
  } else {
    lm(Yvar_F ~ -1 + ., data = as.data.frame(cbind(Xvar_F,
      IMR = dataTable$IMR)), weights = wHvar)
  }
  if (any(is.na(ols2step$coefficients))) {
    stop("at least one of the 'frontierF' coefficients is NA: ",
      paste(colnames(Xvar_F)[is.na(ols2step$coefficients)],
        collapse = ", "), "This may be due to a singular matrix
   due to potential perfect multicollinearity",
      call. = FALSE)
  }
  ols2stepParam <- ols2step$coefficients
  ols2stepSigmasq <- summary(ols2step)$sigma^2
  ols2stepStder <- sqrt(diag(vcov(ols2step)))
  ols2stepLoglik <- logLik(ols2step)[1]
  dataTable <- mutate(dataTable, ols2stepResiduals = residuals(ols2step),
    olsFitted = fitted(ols2step))
  ## skewness is run on the 2step model -------
  ols2stepSkew <- skewness(dataTable[["ols2stepResiduals"]])
  ols2stepM3Okay <- if (S * ols2stepSkew < 0) {
    "Residuals have the expected skeweness"
  } else {
    "Residuals do not have the expected skeweness"
  }
  if (S * ols2stepSkew > 0) {
    warning("The residuals of the OLS are ", if (S == 1) {
      " right"
    } else {
      "left"
    }, "-skewed. This may indicate the absence of inefficiency or
  model misspecification or sample 'bad luck'",
      call. = FALSE)
  }
  CoelliM3Test <- c(z = sum(dataTable[["ols2stepResiduals"]]^3)/N_F/sqrt(6 *
    (sum(dataTable[["ols2stepResiduals"]]^2)/N_F)^3/N_F),
    p.value = 2 * pnorm(-abs(sum(dataTable[["ols2stepResiduals"]]^3)/N_F/sqrt(6 *
      (sum(dataTable[["ols2stepResiduals"]]^2)/N_F)^3/N_F))))
  AgostinoTest <- dagoTest(dataTable[["ols2stepResiduals"]])
  class(AgostinoTest) <- "dagoTest"
  # Step 2: MLE arguments -------
  FunArgs <- if (lType %in% c("kronrod", "hcubature", "pcubature")) {
    list(start = start, olsParam = ols2stepParam, dataTable = dataTable,
      nXvar = nXvar_F, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar_F, Xvar = Xvar_F,
      wHvar = wHvar, S = S, PREDICTIONS = dataTable[["PROBIT_PREDICTIONS"]],
      uBound = uBound, subdivisions = Nsub, intol, method = method,
      printInfo = printInfo, itermax = itermax, stepmax = stepmax,
      tol = tol, gradtol = gradtol, hessianType = hessianType,
      qac = qac)
  } else {
    if (lType == "ghermite") {
      list(start = start, olsParam = ols2stepParam, dataTable = dataTable,
        nXvar = nXvar_F, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar_F,
        Xvar = Xvar_F, wHvar = wHvar, S = S, PREDICTIONS = dataTable[["PROBIT_PREDICTIONS"]],
        gH = gH, N = N_F, method = method, printInfo = printInfo,
        itermax = itermax, stepmax = stepmax, tol = tol,
        gradtol = gradtol, hessianType = hessianType,
        qac = qac)
    } else {
      if (lType == "msl") {
        list(start = start, olsParam = ols2stepParam,
          dataTable = dataTable, nXvar = nXvar_F, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar_F, Xvar = Xvar_F, wHvar = wHvar,
          S = S, PREDICTIONS = dataTable[["PROBIT_PREDICTIONS"]],
          FiMat = FiMat, N = N_F, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol,
          gradtol = gradtol, hessianType = hessianType,
          qac = qac)
      }
    }
  }
  ## MLE run -------
  cat("Second step Frontier model...\n")
  mleList <- tryCatch(switch(lType, kronrod = do.call(halfnormAlgOpt_ss_GK,
    FunArgs), hcubature = do.call(halfnormAlgOpt_ss_HCUB,
    FunArgs), pcubature = do.call(halfnormAlgOpt_ss_PCUB,
    FunArgs), ghermite = do.call(halfnormAlgOpt_ss_GH, FunArgs),
    msl = do.call(halfnormAlgOpt_ss_MSL, FunArgs)), error = function(e) e)
  if (inherits(mleList, "error")) {
    stop("The current error occurs during optimization:\n",
      mleList$message, call. = FALSE)
  }
  # Inverse Hessian + other -------
  mleList$invHessian <- vcovObj(mleObj = mleList$mleObj, hessianType = hessianType,
    method = method, nParm = nParm)
  mleList <- c(mleList, if (method == "ucminf") {
    list(type = "ucminf max.", nIter = unname(mleList$mleObj$info["neval"]),
      status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$value,
      gradient = mleList$mleObj$gradient)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
      list(type = substr(mleList$mleObj$type, 1, 27), nIter = mleList$mleObj$iterations,
        status = mleList$mleObj$message, mleLoglik = mleList$mleObj$maximum,
        gradient = mleList$mleObj$gradient)
    } else {
      if (method == "sr1") {
        list(type = "SR1 max.", nIter = mleList$mleObj$iterations,
          status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
          gradient = mleList$mleObj$gradient)
      } else {
        if (method == "mla") {
          list(type = "Lev. Marquardt max.", nIter = mleList$mleObj$ni,
          status = switch(mleList$mleObj$istop, `1` = "convergence criteria were satisfied",
            `2` = "maximum number of iterations was reached",
            `4` = "algorithm encountered a problem in the function computation"),
          mleLoglik = -mleList$mleObj$fn.value, gradient = mleList$mleObj$grad)
        } else {
          if (method == "sparse") {
          list(type = "Sparse Hessian max.", nIter = mleList$mleObj$iterations,
            status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
            gradient = mleList$mleObj$gradient)
          } else {
          if (method == "nlminb") {
            list(type = "nlminb max.", nIter = mleList$mleObj$iterations,
            status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$objective,
            gradient = mleList$mleObj$gradient)
          }
          }
        }
      }
    }
  })
  # quick renaming -------
  names(mleList$startVal) <- fName_uvr_selectioncross(Xvar = Xvar_F,
    uHvar = uHvar, vHvar = vHvar)
  names(mleList$mlParam) <- names(mleList$startVal)
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
  names(mleList$gradient) <- names(mleList$mlParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
  # Return object -------
  mlDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M")
  dataTable$mlResiduals <- Yvar_F - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar_F]),
    t(Xvar_F)))
  dataTable$mlFitted <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar_F]),
    t(Xvar_F)))
  dataTable$logL_OBS <- mleList$mleObj$logL_OBS
  returnObj <- list()
  returnObj$call <- cl
  returnObj$selectionF <- selectionF
  returnObj$frontierF <- frontierF
  returnObj$formula <- formula
  returnObj$S <- S
  returnObj$typeSfa <- typeSfa
  returnObj$Ninit <- N_S
  returnObj$Nobs <- N_F
  returnObj$nXvar <- nXvar_F
  returnObj$logDepVar <- logDepVar
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- mleList$startVal
  returnObj$dataTable <- dataTable
  returnObj$lpmObj <- lpm
  returnObj$probitObj <- probRes
  returnObj$ols2stepParam <- ols2stepParam
  returnObj$ols2stepSigmasq <- ols2stepSigmasq
  returnObj$ols2stepStder <- ols2stepStder
  returnObj$ols2stepLoglik <- ols2stepLoglik
  returnObj$ols2stepSkew <- ols2stepSkew
  returnObj$ols2stepM3Okay <- ols2stepM3Okay
  returnObj$CoelliM3Test <- CoelliM3Test
  returnObj$AgostinoTest <- AgostinoTest
  returnObj$isWeights <- !all.equal(wHvar, rep(1, N_F))
  returnObj$lType <- lType
  returnObj$optType <- mleList$type
  returnObj$nIter <- mleList$nIter
  returnObj$optStatus <- mleList$status
  returnObj$startLoglik <- mleList$startLoglik
  returnObj$mlLoglik <- mleList$mleLoglik
  returnObj$mlParam <- mleList$mlParam
  returnObj$gradient <- mleList$gradient
  returnObj$gradL_OBS <- mleList$mleObj$gradL_OBS
  returnObj$gradientNorm <- sqrt(sum(mleList$gradient^2))
  returnObj$invHessian <- mleList$invHessian
  returnObj$hessianType <- if (hessianType == 1) {
    "Analytic Hessian"
  } else {
    if (hessianType == 2) {
      "BHHH Hessian"
    }
  }
  returnObj$mlDate <- mlDate
  if (lType == "msl") {
    returnObj$simDist <- simDist
    returnObj$Nsim <- Nsim
    returnObj$FiMat <- FiMat
  } else {
    if (lType == "ghermite") {
      returnObj$gHermiteData <- gH
    } else {
      returnObj$Nsub <- Nsub
      returnObj$intol <- intol
      returnObj$uBound <- uBound
    }
  }
  rm(mleList)
  class(returnObj) <- "selectioncross"
  return(returnObj)
}

# print for selectioncross ----------
#' @rdname selectioncross
#' @export
print.selectioncross <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat("Normal-Half Normal Sample Selection Stochastic Frontier Model",
    "\n")
  cat("Data Type: Cross Section \n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}

# Bread for Sandwich Estimator ----------
#' @rdname selectioncross
#' @export
bread.selectioncross <- function(x, ...) {
  if (x$hessianType == "Analytic Hessian") {
    return(x$invHessian * x$Nobs)
  } else {
    cat("Computing Analytical Hessian")
    Yvar <- model.response(model.frame(x$formula, data = x$dataTable))
    Xvar <- model.matrix(x$formula, rhs = 1, data = x$dataTable)
    uHvar <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
    vHvar <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
    if (x$lType == "kronrod") {
      hessAnalytical <- chesshalfnormlike_ss_GK(x$mlParam,
        nXvar = ncol(Xvar), nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = x$dataTable$weights, S = x$S, PREDICTIONS = x$dataTable[["PROBIT_PREDICTIONS"]],
        uBound = x$uBound, subdivisions = x$Nsub, intol = x$intol)
    } else {
      if (x$lType == "hcubature") {
        hessAnalytical <- chesshalfnormlike_ss_HCUB(x$mlParam,
          nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
          nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
          S = x$S, PREDICTIONS = x$dataTable[["PROBIT_PREDICTIONS"]],
          uBound = x$uBound, subdivisions = x$Nsub, intol = x$intol)
      } else {
        if (x$lType == "pcubature") {
          hessAnalytical <- chesshalfnormlike_ss_PCUB(x$mlParam,
          nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
          nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
          S = x$S, PREDICTIONS = x$dataTable[["PROBIT_PREDICTIONS"]],
          uBound = x$uBound, subdivisions = x$Nsub,
          intol = x$intol)
        } else {
          if (x$lType == "ghermite") {
          hessAnalytical <- chesshalfnormlike_ss_GH(x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
            S = x$S, PREDICTIONS = x$dataTable[["PROBIT_PREDICTIONS"]],
            gH = x$gHermiteData, N = x$Nobs)
          } else {
          if (x$lType == "msl") {
            hessAnalytical <- chesshalfnormlike_ss_MSL(x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar,
            vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
            wHvar = x$dataTable$weights, S = x$S,
            PREDICTIONS = x$dataTable[["PROBIT_PREDICTIONS"]],
            FiMat = x$FiMat, N = x$Nobs)
          }
          }
        }
      }
    }
    invHess <- invHess_fun(hess = hessAnalytical)
    colnames(invHess) <- rownames(invHess) <- names(x$mlParam)
    return(invHess * x$Nobs)
  }
}

# Gradients Evaluated at each Observation ----------
#' @rdname selectioncross
#' @export
estfun.selectioncross <- function(x, ...) {
  return(x$gradL_OBS)
}
