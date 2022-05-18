################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Model: Contaminated Noise Stochastic Frontier Analysis                       #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Contaminated noise stochastic frontier estimation using cross-section data
#'
#' @description
#' \code{\link{cnsfcross}} is a symbolic formula-based function for the
#' estimation of stochastic frontier models in the presence of efficient and 
#' inefficient groups of observation, in the case of cross-sectional or
#' pooled cross-section data, using maximum (simulated) likelihood - M(S)L.
#'
#' The function accounts for heteroscedasticity in both one-sided and two-sided
#' error terms as in Reifschneider and Stevenson (1991), Caudill and Ford
#' (1993), Caudill \emph{et al.} (1995) and Hadri (1999), but also
#' heterogeneity in the mean of the pre-truncated distribution as in Kumbhakar
#' \emph{et al.} (1991), Huang and Liu (1994) and Battese and Coelli (1995).
#'
#' Ten distributions are possible for the one-sided error term and nine
#' optimization algorithms are available.
#'
#' @aliases cnsfcross print.cnsfcross
#'
#' @param formula A symbolic description of the model to be estimated based on
#' the generic function \code{formula} (see section \sQuote{Details}).
#' @param muhet A one-part formula to consider heterogeneity in the mean of the
#' pre-truncated distribution (see section \sQuote{Details}).
#' @param uhet A one-part formula to consider heteroscedasticity in the
#' one-sided error variance (see section \sQuote{Details}).
#' @param vhet A one-part formula to consider heteroscedasticity in the
#' two-sided error variance (see section \sQuote{Details}).
#' @param thet A one-part formula to account for inefficient and efficient 
#' groups of observations (two classes).
#' @param logDepVar Logical. Informs whether the dependent variable is logged
#' (\code{TRUE}) or not (\code{FALSE}). Default = \code{TRUE}.
#' @param data The data frame containing the data.
#' @param subset An optional vector specifying a subset of observations to be
#' used in the optimization process.
#' @param weights An optional vector of weights to be used for weighted 
#' log-likelihood. Should be \code{NULL} or numeric vector with positive values.
#' When \code{NULL}, a numeric vector of 1 is used.
#' @param wscale Logical. When \code{weights} is not \code{NULL}, a scaling 
#' transformation is used such that the the \code{weights} sums to the sample 
#' size. Default \code{TRUE}. When \code{FALSE} no scaling is used.
#' @param S If \code{S = 1} (default), a production (profit) frontier is
#' estimated: \eqn{\epsilon_i = v_i-u_i}. If \code{S = -1}, a cost frontier is
#' estimated: \eqn{\epsilon_i = v_i+u_i}.
#' @param udist Character string. Default = \code{'hnormal'}. Distribution
#' specification for the one-sided error term. 10 different distributions are
#' available: \itemize{ \item \code{'hnormal'}, for the half normal
#' distribution (Aigner \emph{et al.} 1977, Meeusen and Vandenbroeck 1977)
#' \item \code{'exponential'}, for the exponential distribution \item
#' \code{'tnormal'} for the truncated normal distribution (Stevenson 1980)
#' \item \code{'rayleigh'}, for the Rayleigh distribution (Hajargasht 2015)
#' \item \code{'uniform'}, for the uniform distribution (Li 1996, Nguyen 2010)
#' \item \code{'gamma'}, for the Gamma distribution (Greene 2003) \item
#' \code{'lognormal'}, for the log normal distribution (Migon and Medici 2001,
#' Wang and Ye 2020) \item \code{'weibull'}, for the Weibull distribution
#' (Tsionas 2007) \item \code{'genexponential'}, for the generalized
#' exponential distribution (Papadopoulos 2020) \item \code{'tslaplace'}, for
#' the truncated skewed Laplace distribution (Wang 2012). }
#' @param sigmauType Character string. Default = \code{'common'}. Nature of the 
#' two-sided error term between classes. Two possibilities: \code{'common'} for 
#' common error term \eqn{v} and \code{'different'} for different error term 
#' \eqn{v}.
#' @param linkF. Character string. Link function used for class membership
#' probability. Four possibilities are available: \code{logit} (default), 
#' \code{probit}, \code{cauchit}, \code{cloglog}. 
#' @param start Numeric vector. Optional starting values for the maximum
#' likelihood (ML) estimation.
#' @param method Optimization algorithm used for the estimation.  Default =
#' \code{'bfgs'}. 11 algorithms are available: \itemize{ \item \code{'bfgs'},
#' for Broyden-Fletcher-Goldfarb-Shanno (see
#' \code{\link[maxLik:maxBFGS]{maxBFGS}}) \item \code{'bhhh'}, for
#' Berndt-Hall-Hall-Hausman (see \code{\link[maxLik:maxBHHH]{maxBHHH}}) \item
#' \code{'nr'}, for Newton-Raphson (see \code{\link[maxLik:maxNR]{maxNR}})
#' \item \code{'nm'}, for Nelder-Mead (see \code{\link[maxLik:maxNM]{maxNM}})
#' \item \code{'cg'}, for Conjugate Gradient (see 
#' \code{\link[maxLik:maxCG]{maxCG}}) \item \code{'sann'}, for Simulated 
#' Annealing (see \code{\link[maxLik:maxSANN]{maxSANN}})
#'  \item \code{'ucminf'}, implements a quasi-Newton type with BFGS updating 
#'  of the inverse Hessian and soft line search with a trust region type 
#'  monitoring of the input to the line search algorithm (see 
#'  \code{\link[ucminf:ucminf]{ucminf}})
#' \item \code{'mla'}, for general-purpose optimization based on
#' Marquardt-Levenberg algorithm (see \code{\link[marqLevAlg:mla]{mla}})
#' \item \code{'sr1'}, for Symmetric Rank 1 (see
#' \code{\link[trustOptim:trust.optim]{trust.optim}}) \item \code{'sparse'}, 
#' for trust regions and sparse Hessian (see 
#' \code{\link[trustOptim:trust.optim]{trust.optim}})
#' \item \code{'nlminb'}, for optimization using PORT routines (see
#' \code{\link[stats:nlminb]{nlminb}})}
#' @param hessianType Integer. If \code{1} (Default), analytic/numeric Hessian
#' is returned for all the distributions. If \code{2},  bhhh Hessian is 
#' estimated (\eqn{g'g}).
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
#' @param x an object of class cnsfcross (returned by the function 
#' \code{\link{cnsfcross}}).
#' @param ... additional arguments of frontier are passed to cnsfcross; 
#' additional arguments of the print, bread, estfun, nobs methods are currently 
#' ignored.
#'
#' @details
#' Following Wheat \emph{et al.} (2019) the model can be defined as: 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('y_i = \\\alpha + \\\mathbf{x_i^{\\\prime}}\\\bm{\\\beta} + 
#' v_{1i} - Su_i \\\quad \\\hbox{with probability} \\\quad p,')
#' }
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('y_i = \\\alpha + \\\mathbf{x_i^{\\\prime}}\\\bm{\\\beta} + 
#' v_{2i} - Su_i \\\quad \\\hbox{with probability} \\\quad 1-p')
#' }
#' 
#' where \eqn{i} is the observation, \eqn{y} is the
#' output (cost, revenue, profit), \eqn{x} is the vector of main explanatory
#' variables (inputs and other control variables), \eqn{u} is the one-sided
#' error term with variance \eqn{\sigma_{u}^2}, and \eqn{v_1} and \eqn{v_2} are
#' the two-sided error terms with variance \eqn{\sigma_{v_1}^2} and 
#' \eqn{\sigma_{v_2}^2} respectively.
#'
#' \code{S = 1} in the case of production (profit) frontier function and
#' \code{S = -1} in the case of cost frontier function.
#' 
#' The prior probability of belonging to class \eqn{1} can depend on some
#' covariates using a logit specification: 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('p_i = \\\frac{\\\exp{(Z_{hi}^{\\\prime}\\\theta)}}{1-
#' \\\exp{(Z_{hi}^{\\\prime}\\\theta)}}')
#' }
#' 
#' with \eqn{Z_h} the covariates, \eqn{\bm{\theta}} the vector of coefficients 
#' estimated for the covariates. 
#' 
#' Other link functions are also available. In the case of the probit we have:
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('p_i = \\\Phi\\\\left(Z_{hi}^{\\\prime}\\\theta\\\\right)')
#' }
#' 
#' when the cauchit link is retained we have:
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('p_i = 1/\\\pi\\\arctan(Z_{hi}^{\\\prime}\\\theta)+1/2')
#' }
#' 
#' and finally in the case of the cloglog link we have:
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('p_i = 1-\\\exp\\\\left(-\\\exp(Z_{hi}^{\\\prime}
#' \\\theta)\\\\right)')
#' }
#' 
#' Let 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('\\\epsilon_i = v_i -Su_i')
#' }
#' 
#' In the case of the truncated normal distribution, the convolution of
#' \eqn{u_i} and \eqn{v_i} is:
#'
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('f(\\\epsilon_i)=\\\frac{p_i}{\\\sqrt{\\\\sigma_u^2 + 
#' \\\sigma_{v_1^2}}}\\\phi\\\\left(\\\frac{S\\\epsilon_i + \\\mu}{\\\sqrt{
#' \\\sigma_u^2 + \\\sigma_{v_1^2}}}\\\\right)\\\Phi\\\\left(\\\frac{
#' \\\mu_{1i\\\ast}}{\\\sigma_{1\\\ast}}\\\\right) \\\Big/\\\Phi\\\\left(
#' \\\frac{\\\mu}{\\\sigma_u}\\\\right) + \\\frac{1-p_i}{\\\sqrt{\\\\sigma_u^2 + 
#' \\\sigma_{v_2^2}}}\\\phi\\\\left(\\\frac{S\\\epsilon_i + \\\mu}{\\\sqrt{
#' \\\sigma_u^2 + \\\sigma_{v_2^2}}}\\\\right)\\\Phi\\\\left(\\\frac{
#' \\\mu_{2i\\\ast}}{\\\sigma_{2\\\ast}}\\\\right) \\\Big/\\\Phi\\\\left(
#' \\\frac{\\\mu}{\\\sigma_u}\\\\right)')
#' }
#'
#' where
#'
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('\\\mu_{1i*}=\\\frac{\\\mu\\\sigma_{v_1^2} - 
#' S\\\epsilon_i\\\sigma_u^2}{\\\sigma_u^2 + \\\sigma_{v_1^2}}')
#' }
#'
#' and
#'
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('\\\sigma_{1*}^2 = \\\frac{\\\sigma_u^2 
#' \\\sigma_{v_1^2}}{\\\sigma_u^2 + \\\sigma_{v_1^2}}')
#' }
#' 
#' \eqn{\mu_{2i*}} and \eqn{\sigma_{2*}^2} are similarly obtained.
#'
#' In the case of the half normal distribution the convolution is obtained by
#' setting \eqn{\mu=0}.
#' 
#' Class assignment is based on the largest posterior probability. This
#' probability is obtained using Bayes' rule, as follows for the inefficient 
#' class
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('p_i^\\\ast = \\\frac{p_i \\\times 
#' \\\pi(i, 1)}{f(\\\epsilon_i)}')
#' }
#' 
#' where
#' 
#'  \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('\\\pi(i, 1)=\\\frac{1}{\\\sqrt{\\\\sigma_u^2 + 
#' \\\sigma_v^2}}\\\phi\\\\left(\\\frac{S\\\epsilon_i + \\\mu}{\\\sqrt{
#' \\\sigma_u^2 + \\\sigma_v^2}}\\\\right)\\\Phi\\\\left(\\\frac{
#' \\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)\\\Big/\\\Phi\\\\left(
#' \\\frac{\\\mu}{\\\sigma_u}\\\\right)')
#' }
#' 
#' The model presented so far assumed common inefficiency term \eqn{u} for the 
#' two classes of observations. By setting argument 
#' \code{sigmauType = 'different'}, different inefficiency terms are imposed for 
#' the two groups. The new estimated model is presented as follows:
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('y_i = \\\alpha + \\\mathbf{x_i^{\\\prime}}\\\bm{\\\beta} + 
#' v_{1i} - Su_{1i} \\\quad \\\hbox{with probability} \\\quad p,')
#' }
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('y_i = \\\alpha + \\\mathbf{x_i^{\\\prime}}\\\bm{\\\beta} + 
#' v_{2i} - Su_{2i} \\\quad \\\hbox{with probability} \\\quad 1-p')
#' }
#' 
#' \code{cnsfcross} allows for the maximization of weighted log-likelihood.
#' When option \code{weights} is specified and \code{wscale = TRUE}, the weights
#' is scaled as 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('new_{weights} = sample_{size} \\\times 
#' \\\frac{old_{weights}}{\\\sum(old_{weights})}')
#' }
#' 
#' For difficult problems, non-gradient methods (e.g. \code{nm} or \code{sann}) 
#' can be used to warm start the optimization and zoom in the neighborhood of 
#' the solution. Then a gradient-based methods is recommanded in the second 
#' step. In the case of \code{sann}, we recommand to significantly increase the 
#' iteration limit (e.g. \code{itermax = 20000}). The Conjugate Gradient 
#' (\code{cg}) can also be used in the first stage.
#' 
#' A set of extractor functions for fitted model objects is available for 
#' objects of class \code{'cnsfcross'} including methods to the generic 
#' functions \code{\link[=print.cnsfcross]{print}}, 
#' \code{\link[=summary.cnsfcross]{summary}}, 
#' \code{\link[=coef.cnsfcross]{coef}}, 
#' \code{\link[=fitted.cnsfcross]{fitted}}, 
#' \code{\link[=logLik.cnsfcross]{logLik}}, 
#' \code{\link[=residuals.cnsfcross]{residuals}}, 
#' \code{\link[=vcov.cnsfcross]{vcov}}, 
#' \code{\link[=efficiencies.cnsfcross]{efficiencies}}, 
#' \code{\link[=ic.cnsfcross]{ic}}, 
#' \code{\link[=marginal.cnsfcross]{marginal}},
#' \code{\link[=estfun.cnsfcross]{estfun}} and 
#' \code{\link[=bread.cnsfcross]{bread}} (from the \CRANpkg{sandwich} package), 
#' [lmtest::coeftest()] (from the \CRANpkg{lmtest} package).
#' 
#' @return \code{\link{cnsfcross}} returns a list of class \code{'cnsfcross'}
#' containing the following elements:
#'
#' \item{call}{The matched call.}
#'
#' \item{formula}{The estimated model.}
#'
#' \item{S}{The argument \code{'S'}. See the section \sQuote{Arguments}.}
#'
#' \item{typeSfa}{Character string. 'Contaminated Noise Production/Profit 
#' Frontier, e = v - u' when \code{S = 1} and 'Contaminated Noise Cost Frontier, 
#' e = v + u' when \code{S = -1}.}
#'
#' \item{Nobs}{Number of observations used for optimization.}
#'
#' \item{nXvar}{Number of explanatory variables in the production or cost
#' frontier.}
#'
#' \item{nmuZUvar}{Number of variables explaining heterogeneity in the
#' truncated mean, only if \code{udist = 'tnormal'} or \code{'lognormal'}.}
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
#' M(S)L estimations, and the individual observation log-likelihood. When 
#' \code{weights} is specified an additional variable is also provided in 
#' \code{dataTable}.}
#' 
#' \item{sigmauType}{Character string. Nature of the 
#' two-sided error term between classes. See **Arguments** section.}
#' 
#' \item{linkF}{Character string. Link function retained to define
#' class membership probability. See **Arguments** section.}
#' 
#' \item{initHalf}{When \code{start = NULL} and \code{udist = 'hnormal'}. 
#' Initial ML estimation with half normal distribution for the one-sided error 
#' term. Model to construct the starting values for the cnsf estimation. 
#' Object of class \code{'maxLik'} and \code{'maxim'} returned. 
#' For \code{udist = 'exponential'} it returns \code{initExpo}, 
#' for \code{udist = 'tnormal'} it returns \code{initTrunc}, 
#' for \code{udist = 'rayleigh'} it returns \code{initRay}, 
#' for \code{udist = 'uniform'} it returns \code{initUni}, 
#' for \code{udist = 'gamma'} it returns \code{initGamma}, 
#' for \code{udist = 'lognormal'} it returns \code{initLog}, 
#' for \code{udist = 'genexponential'} it returns \code{initGenexpo}, 
#' for \code{udist = 'tslaplace'} it returns \code{initTSL}.}
#' 
#' \item{isWeights}{Logical. If \code{TRUE} weighted log-likelihood is
#' maximized.}
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
#' \item{simDist}{The argument \code{'simDist'}, only if \code{udist =
#' 'gamma'}, \code{'lognormal'} or , \code{'weibull'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{Nsim}{The argument \code{'Nsim'}, only if \code{udist = 'gamma'},
#' \code{'lognormal'} or , \code{'weibull'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{FiMat}{Matrix of random draws used for MSL, only if \code{udist =
#' 'gamma'}, \code{'lognormal'} or , \code{'weibull'}.}
#'
#' @note For the Halton draws, the code is adapted from the \pkg{mlogit}
#' package.
#'
# @author K Hervé Dakpo
#'
#' @seealso \code{\link[=print.cnsfcross]{print}} for printing \code{cnsfcross} 
#' object.
#' 
#' \code{\link[=summary.cnsfcross]{summary}} for creating and printing
#' summary results.
#'
#' \code{\link[=coef.cnsfcross]{coef}} for extracting coefficients of the
#' estimation.
#'
#' \code{\link[=efficiencies.cnsfcross]{efficiencies}} for computing
#' (in-)efficiency estimates.
#'
#' \code{\link[=fitted.cnsfcross]{fitted}} for extracting the fitted frontier
#' values.
#'
#' \code{\link[=ic.cnsfcross]{ic}} for extracting information criteria.
#'
#' \code{\link[=logLik.cnsfcross]{logLik}} for extracting log-likelihood
#' value(s) of the estimation.
#'
#' \code{\link[=marginal.cnsfcross]{marginal}} for computing marginal effects of
#' inefficiency drivers.
#'
#' \code{\link[=residuals.cnsfcross]{residuals}} for extracting residuals of the
#' estimation.
#'
#' \code{\link[=vcov.cnsfcross]{vcov}} for computing the variance-covariance
#' matrix of the coefficients.
#' 
#' \code{\link[=bread.cnsfcross]{bread}} for bread for sandwich estimator.
#' 
#' \code{\link[=estfun.cnsfcross]{estfun}} for gradient extraction for each 
#' observation.
#'
#' @references Aigner, D., Lovell, C. A. K., and Schmidt, P. 1977. Formulation
#' and estimation of stochastic frontier production function models.
#' \emph{Journal of Econometrics}, \bold{6}(1), 21--37.
#'
#' Battese, G. E., and Coelli, T. J. 1995. A model for technical inefficiency
#' effects in a stochastic frontier production function for panel data.
#' \emph{Empirical Economics}, \bold{20}(2), 325--332.
#'
#' Caudill, S. B., and Ford, J. M. 1993. Biases in frontier estimation due to
#' heteroscedasticity. \emph{Economics Letters}, \bold{41}(1), 17--20.
#'
#' Caudill, S. B., Ford, J. M., and Gropper, D. M. 1995. Frontier estimation
#' and firm-specific inefficiency measures in the presence of
#' heteroscedasticity. \emph{Journal of Business & Economic Statistics},
#' \bold{13}(1), 105--111.
#'
#' Greene, W. H. 2003. Simulated likelihood estimation of the normal-Gamma
#' stochastic frontier function. \emph{Journal of Productivity Analysis},
#' \bold{19}(2-3), 179--190.
#'
#' Hadri, K. 1999. Estimation of a doubly heteroscedastic stochastic frontier
#' cost function. \emph{Journal of Business & Economic Statistics},
#' \bold{17}(3), 359--363.
#'
#' Hajargasht, G. 2015. Stochastic frontiers with a Rayleigh distribution.
#' \emph{Journal of Productivity Analysis}, \bold{44}(2), 199--208.
#'
#' Huang, C. J., and Liu, J.-T. 1994. Estimation of a non-neutral stochastic
#' frontier production function. \emph{Journal of Productivity Analysis},
#' \bold{5}(2), 171--180.
#' 
#' Kumbhakar, S. C., Parmeter, C. F., and Tsionas, E. G. 2013) 
#' A zero inefficiency stochastic frontier model. 
#' \emph{Journal of Econometrics}, \bold{172}(1), 66--76.
#'
#' Li, Q. 1996. Estimating a stochastic production frontier when the adjusted
#' error is symmetric. \emph{Economics Letters}, \bold{52}(3), 221--228.
#'
#' Meeusen, W., and Vandenbroeck, J. 1977. Efficiency estimation from
#' Cobb-Douglas production functions with composed error. \emph{International
#' Economic Review}, \bold{18}(2), 435--445.
#'
#' Migon, H. S., and Medici, E. V. 2001. Bayesian hierarchical models for
#' stochastic production frontier. Lacea, Montevideo, Uruguay.
#'
#' Nguyen, N. B. 2010. Estimation of technical efficiency in stochastic
#' frontier analysis. PhD dissertation, Bowling Green State University, August.
#'
#' Papadopoulos, A. 2021. Stochastic frontier models using the generalized
#' exponential distribution. \emph{Journal of Productivity Analysis},
#' \bold{55}:15--29.
#'
#' Reifschneider, D., and Stevenson, R. 1991. Systematic departures from the
#' frontier: A framework for the analysis of firm inefficiency.
#' \emph{International Economic Review}, \bold{32}(3), 715--723.
#'
#' Stevenson, R. E. 1980. Likelihood Functions for Generalized Stochastic
#' Frontier Estimation. \emph{Journal of Econometrics}, \bold{13}(1), 57--66.
#'
#' Tsionas, E. G. 2007. Efficiency measurement with the Weibull stochastic
#' frontier. \emph{Oxford Bulletin of Economics and Statistics}, \bold{69}(5),
#' 693--706.
#'
#' Wang, K., and Ye, X. 2020. Development of alternative stochastic frontier
#' models for estimating time-space prism vertices. \emph{Transportation}.
#'
#' Wang, J. 2012. A normal truncated skewed-Laplace model in stochastic
#' frontier analysis. Master thesis, Western Kentucky University, May.
#' 
#' Wheat, P., Stead, A. D., & Greene, W. H. 2019. Controlling for Outliers in 
#' Efficiency Analysis: A Contaminated Normal-Half Normal Stochastic Frontier 
#' Model. Working Paper.
#'
#' @keywords models optimize cross-section likelihood
#'
#' @examples
#' 
#' ## Using data on Spanish dairy farms
#' # Cobb Douglas (production function) half normal distribution
#' 
#' cb_cnsf_h <- cnsfcross(formula = YIT ~ X1 + X2 + X3 + X4, udist = 'hnormal', 
#' data = dairyspain, S = 1, method = 'bfgs')
#' 
#' summary(cb_cnsf_h)
#' 
#' @export 
cnsfcross <- function(formula, muhet, uhet, vhet, thet, logDepVar = TRUE,
  data, subset, weights, wscale = TRUE, S = 1L, udist = "hnormal",
  sigmauType = "common", linkF = "logit", start = NULL, method = "bfgs",
  hessianType = 1, simType = "halton", Nsim = 300, prime = 2L,
  burn = 10, antithetics = FALSE, seed = 12345, itermax = 2000L,
  printInfo = FALSE, tol = 1e-12, gradtol = 1e-06, stepmax = 0.1,
  qac = "marquardt") {
  # sigma_v model check -------
  sigmauType <- tolower(sigmauType)
  if (!(sigmauType %in% c("common", "different"))) {
    stop("Unknown 'sigmauType': ", paste(sigmauType), call. = FALSE)
  }
  # link functions check -------
  linkF <- tolower(linkF)
  if (!(linkF %in% c("logit", "probit", "cauchit", "cloglog"))) {
    stop("Unknown 'linkF': ", paste(linkF), call. = FALSE)
  }
  # u distribution check -------
  udist <- tolower(udist)
  if (!(udist %in% c("hnormal", "exponential", "tnormal", "rayleigh",
    "uniform", "gamma", "lognormal", "weibull", "genexponential",
    "tslaplace"))) {
    stop("Unknown inefficiency distribution: ", paste(udist),
      call. = FALSE)
  }
  # Formula manipulation -------
  if (length(Formula(formula))[2] != 1) {
    stop("argument 'formula' must have one RHS part", call. = FALSE)
  }
  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"), names(mc),
    nomatch = 0L)
  mc <- mc[c(1L, m)]
  mc$drop.unused.levels <- TRUE
  formula <- interCheckMain(formula = formula)
  if (!missing(muhet)) {
    muhet <- clhsCheck_mu(formula = muhet, scaling = FALSE)
  } else {
    muhet <- ~1
  }
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
  if (!missing(thet)) {
    thet <- clhsCheck_q(formula = thet)
  } else {
    thet <- ~1
  }
  formula <- formDist_cnsfcross(udist = udist, formula = formula,
    muhet = muhet, uhet = uhet, vhet = vhet, qhet = thet)
  # Generate required datasets -------
  if (missing(data)) {
    data <- environment(formula)
  }
  mc$formula <- formula
  mc$na.action <- na.pass
  mc[[1L]] <- quote(model.frame)
  mc <- eval(mc, parent.frame())
  validObs <- rowSums(is.na(mc) | is.infinite.data.frame(mc)) ==
    0
  Yvar <- model.response(mc, "numeric")
  Yvar <- Yvar[validObs]
  mtX <- terms(formula, data = data, rhs = 1)
  Xvar <- model.matrix(mtX, mc)
  Xvar <- Xvar[validObs, , drop = FALSE]
  nXvar <- ncol(Xvar)
  N <- nrow(Xvar)
  if (N == 0L) {
    stop("0 (non-NA) cases", call. = FALSE)
  }
  wHvar <- as.vector(model.weights(mc))
  if (length(wscale) != 1 || !is.logical(wscale[1])) {
    stop("argument 'wscale' must be a single logical value",
      call. = FALSE)
  }
  if (!is.null(wHvar)) {
    if (!is.numeric(wHvar)) {
      stop("'weights' must be a numeric vector", call. = FALSE)
    } else {
      if (any(wHvar < 0 | is.na(wHvar)))
        stop("missing or negative weights not allowed",
          call. = FALSE)
    }
    if (wscale) {
      wHvar <- wHvar/sum(wHvar) * N
    }
  } else {
    wHvar <- rep(1, N)
  }
  if (length(Yvar) != nrow(Xvar)) {
    stop(paste("the number of observations of the dependent variable (",
      length(Yvar), ") must be the same to the number of observations of the 
      exogenous variables (",
      nrow(Xvar), ")", sep = ""), call. = FALSE)
  }
  if (udist %in% c("tnormal", "lognormal")) {
    mtmuH <- delete.response(terms(formula, data = data,
      rhs = 2))
    muHvar <- model.matrix(mtmuH, mc)
    muHvar <- muHvar[validObs, , drop = FALSE]
    nmuZUvar <- ncol(muHvar)
    mtuH <- delete.response(terms(formula, data = data, rhs = 3))
    uHvar <- model.matrix(mtuH, mc)
    uHvar <- uHvar[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar)
    mtvH <- delete.response(terms(formula, data = data, rhs = 4))
    vHvar <- model.matrix(mtvH, mc)
    vHvar <- vHvar[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar)
    mtZ <- delete.response(terms(formula, data = data, rhs = 5))
    Zvar <- model.matrix(mtZ, mc)
    Zvar <- Zvar[validObs, , drop = FALSE]
    nZHvar <- ncol(Zvar)
  } else {
    mtuH <- delete.response(terms(formula, data = data, rhs = 2))
    uHvar <- model.matrix(mtuH, mc)
    uHvar <- uHvar[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar)
    mtvH <- delete.response(terms(formula, data = data, rhs = 3))
    vHvar <- model.matrix(mtvH, mc)
    vHvar <- vHvar[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar)
    mtZ <- delete.response(terms(formula, data = data, rhs = 4))
    Zvar <- model.matrix(mtZ, mc)
    Zvar <- Zvar[validObs, , drop = FALSE]
    nZHvar <- ncol(Zvar)
  }
  # Check other supplied options -------
  if (length(S) != 1 || !(S %in% c(-1L, 1L))) {
    stop("argument 'S' must equal either 1 or -1: 1 for production or profit 
    frontier and -1 for cost frontier",
      call. = FALSE)
  }
  typeSfa <- if (S == 1L) {
    "Contaminated Noise Production/Profit Frontier, e = v - u"
  } else {
    "Contaminated Noise Cost Frontier, e = v + u"
  }
  if (length(logDepVar) != 1 || !is.logical(logDepVar[1])) {
    stop("argument 'logDepVar' must be a single logical value",
      call. = FALSE)
  }
  # Number of parameters -------
  if (sigmauType == "common") {
    nParm <- if (udist %in% c("tnormal", "lognormal")) {
      nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar
    } else {
      if (udist %in% c("gamma", "weibull", "tslaplace")) {
        nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1
      } else {
        nXvar + nuZUvar + 2 * nvZVvar + nZHvar
      }
    }
  } else {
    if (sigmauType == "different") {
      nParm <- if (udist %in% c("tnormal", "lognormal")) {
        nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar +
          nZHvar
      } else {
        if (udist %in% c("gamma", "weibull", "tslaplace")) {
          nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar +
          1
        } else {
          nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar
        }
      }
    }
  }
  # Checking starting values when provided -------
  if (!is.null(start)) {
    if (length(start) != nParm) {
      stop("Wrong number of initial values: model has ",
        nParm, " parameters", call. = FALSE)
    }
  }
  if (nParm > N) {
    stop("Model has more parameters than observations", call. = FALSE)
  }
  # Check algorithms -------
  method <- tolower(method)
  if (!(method %in% c("ucminf", "bfgs", "bhhh", "nr", "nm",
    "cg", "sann", "sr1", "mla", "sparse", "nlminb"))) {
    stop("Unknown or non-available optimization algorithm: ",
      paste(method), call. = FALSE)
  }
  # Check hessian type
  if (length(hessianType) != 1 || !(hessianType %in% c(1L,
    2L))) {
    stop("argument 'hessianType' must equal either 1 or 2",
      call. = FALSE)
  }
  # Draws for SML -------
  if (udist %in% c("gamma", "lognormal", "weibull")) {
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
    FiMat <- drawMat(N = N, Nsim = Nsim, simType = simType,
      prime = prime, burn = burn + 1, antithetics = antithetics,
      seed = seed)
  }
  # Other optimization options -------
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
  # Step 1: OLS -------
  olsRes <- if (colnames(Xvar)[1] == "(Intercept)") {
    if (dim(Xvar)[2] == 1) {
      lm(Yvar ~ 1)
    } else {
      lm(Yvar ~ ., data = as.data.frame(Xvar[, -1]), weights = wHvar)
    }
  } else {
    lm(Yvar ~ -1 + ., data = as.data.frame(Xvar), weights = wHvar)
  }
  if (any(is.na(olsRes$coefficients))) {
    stop("at least one of the OLS coefficients is NA: ",
      paste(colnames(Xvar)[is.na(olsRes$coefficients)],
        collapse = ", "), "This may be due to a singular matrix
   due to potential perfect multicollinearity",
      call. = FALSE)
  }
  olsParam <- c(olsRes$coefficients)
  if (inherits(data, "plm.dim")) {
    dataTable <- data[validObs, 1:2]
  } else {
    dataTable <- data.frame(IdObs = c(1:sum(validObs)))
  }
  dataTable <- as_tibble(cbind(dataTable, data[, all.vars(terms(formula))],
    weights = wHvar))
  dataTable <- mutate(dataTable, olsResiduals = residuals(olsRes),
    olsFitted = fitted(olsRes))
  olsSkew <- skewness(dataTable[["olsResiduals"]])
  if (S * olsSkew > 0) {
    warning("The residuals of the OLS are", if (S == 1) {
      " right"
    } else {
      " left"
    }, "-skewed. This may indicate the absence of inefficiency or
  model misspecification or sample 'bad luck'",
      call. = FALSE)
  }
  # Step 2: MLE arguments -------
  FunArgs <- if (udist == "tnormal") {
    list(start = start, olsParam = olsParam, dataTable = dataTable,
      nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Zvar = Zvar, nZHvar = nZHvar,
      Xvar = Xvar, S = S, wHvar = wHvar, method = method,
      printInfo = printInfo, itermax = itermax, stepmax = stepmax,
      tol = tol, gradtol = gradtol, hessianType = hessianType,
      qac = qac)
  } else {
    if (udist == "lognormal") {
      list(start = start, olsParam = olsParam, dataTable = dataTable,
        nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Zvar = Zvar, nZHvar = nZHvar,
        Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
        method = method, printInfo = printInfo, itermax = itermax,
        stepmax = stepmax, tol = tol, gradtol = gradtol,
        hessianType = hessianType, qac = qac)
    } else {
      if (udist %in% c("gamma", "weibull")) {
        list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Zvar = Zvar, nZHvar = nZHvar, Xvar = Xvar,
          S = S, wHvar = wHvar, N = N, FiMat = FiMat,
          method = method, printInfo = printInfo, itermax = itermax,
          stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, qac = qac)
      } else {
        list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Zvar = Zvar, nZHvar = nZHvar, Xvar = Xvar,
          S = S, wHvar = wHvar, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol,
          gradtol = gradtol, hessianType = hessianType,
          qac = qac)
      }
    }
  }
  ## MLE run -------
  if (linkF == "logit") {
    if (sigmauType == "common") {
      mleList <- tryCatch(switch(udist, hnormal = do.call(cnsfhalfnormAlgOpt_logit,
        FunArgs), exponential = do.call(cnsfexponormAlgOpt_logit,
        FunArgs), tnormal = do.call(cnsftruncnormAlgOpt_logit,
        FunArgs), rayleigh = do.call(cnsfraynormAlgOpt_logit,
        FunArgs), gamma = do.call(cnsfgammanormAlgOpt_logit,
        FunArgs), uniform = do.call(cnsfuninormAlgOpt_logit,
        FunArgs), lognormal = do.call(cnsflognormAlgOpt_logit,
        FunArgs), weibull = do.call(cnsfweibullnormAlgOpt_logit,
        FunArgs), genexponential = do.call(cnsfgenexponormAlgOpt_logit,
        FunArgs), tslaplace = do.call(cnsftslnormAlgOpt_logit,
        FunArgs)), error = function(e) e)
    } else {
      if (sigmauType == "different") {
        mleList <- tryCatch(switch(udist, hnormal = do.call(mcesfhalfnormAlgOpt_logit,
          FunArgs), exponential = do.call(mcesfexponormAlgOpt_logit,
          FunArgs), tnormal = do.call(mcesftruncnormAlgOpt_logit,
          FunArgs), rayleigh = do.call(mcesfraynormAlgOpt_logit,
          FunArgs), gamma = do.call(mcesfgammanormAlgOpt_logit,
          FunArgs), uniform = do.call(mcesfuninormAlgOpt_logit,
          FunArgs), lognormal = do.call(mcesflognormAlgOpt_logit,
          FunArgs), weibull = do.call(mcesfweibullnormAlgOpt_logit,
          FunArgs), genexponential = do.call(mcesfgenexponormAlgOpt_logit,
          FunArgs), tslaplace = do.call(mcesftslnormAlgOpt_logit,
          FunArgs)), error = function(e) e)
      }
    }
  } else {
    if (linkF == "probit") {
      if (sigmauType == "common") {
        mleList <- tryCatch(switch(udist, hnormal = do.call(cnsfhalfnormAlgOpt_probit,
          FunArgs), exponential = do.call(cnsfexponormAlgOpt_probit,
          FunArgs), tnormal = do.call(cnsftruncnormAlgOpt_probit,
          FunArgs), rayleigh = do.call(cnsfraynormAlgOpt_probit,
          FunArgs), gamma = do.call(cnsfgammanormAlgOpt_probit,
          FunArgs), uniform = do.call(cnsfuninormAlgOpt_probit,
          FunArgs), lognormal = do.call(cnsflognormAlgOpt_probit,
          FunArgs), weibull = do.call(cnsfweibullnormAlgOpt_probit,
          FunArgs), genexponential = do.call(cnsfgenexponormAlgOpt_probit,
          FunArgs), tslaplace = do.call(cnsftslnormAlgOpt_probit,
          FunArgs)), error = function(e) e)
      } else {
        if (sigmauType == "different") {
          mleList <- tryCatch(switch(udist, hnormal = do.call(mcesfhalfnormAlgOpt_probit,
          FunArgs), exponential = do.call(mcesfexponormAlgOpt_probit,
          FunArgs), tnormal = do.call(mcesftruncnormAlgOpt_probit,
          FunArgs), rayleigh = do.call(mcesfraynormAlgOpt_probit,
          FunArgs), gamma = do.call(mcesfgammanormAlgOpt_probit,
          FunArgs), uniform = do.call(mcesfuninormAlgOpt_probit,
          FunArgs), lognormal = do.call(mcesflognormAlgOpt_probit,
          FunArgs), weibull = do.call(mcesfweibullnormAlgOpt_probit,
          FunArgs), genexponential = do.call(mcesfgenexponormAlgOpt_probit,
          FunArgs), tslaplace = do.call(mcesftslnormAlgOpt_probit,
          FunArgs)), error = function(e) e)
        }
      }
    } else {
      if (linkF == "cauchit") {
        if (sigmauType == "common") {
          mleList <- tryCatch(switch(udist, hnormal = do.call(cnsfhalfnormAlgOpt_cauchit,
          FunArgs), exponential = do.call(cnsfexponormAlgOpt_cauchit,
          FunArgs), tnormal = do.call(cnsftruncnormAlgOpt_cauchit,
          FunArgs), rayleigh = do.call(cnsfraynormAlgOpt_cauchit,
          FunArgs), gamma = do.call(cnsfgammanormAlgOpt_cauchit,
          FunArgs), uniform = do.call(cnsfuninormAlgOpt_cauchit,
          FunArgs), lognormal = do.call(cnsflognormAlgOpt_cauchit,
          FunArgs), weibull = do.call(cnsfweibullnormAlgOpt_cauchit,
          FunArgs), genexponential = do.call(cnsfgenexponormAlgOpt_cauchit,
          FunArgs), tslaplace = do.call(cnsftslnormAlgOpt_cauchit,
          FunArgs)), error = function(e) e)
        } else {
          if (sigmauType == "different") {
          mleList <- tryCatch(switch(udist, hnormal = do.call(mcesfhalfnormAlgOpt_cauchit,
            FunArgs), exponential = do.call(mcesfexponormAlgOpt_cauchit,
            FunArgs), tnormal = do.call(mcesftruncnormAlgOpt_cauchit,
            FunArgs), rayleigh = do.call(mcesfraynormAlgOpt_cauchit,
            FunArgs), gamma = do.call(mcesfgammanormAlgOpt_cauchit,
            FunArgs), uniform = do.call(mcesfuninormAlgOpt_cauchit,
            FunArgs), lognormal = do.call(mcesflognormAlgOpt_cauchit,
            FunArgs), weibull = do.call(mcesfweibullnormAlgOpt_cauchit,
            FunArgs), genexponential = do.call(mcesfgenexponormAlgOpt_cauchit,
            FunArgs), tslaplace = do.call(mcesftslnormAlgOpt_cauchit,
            FunArgs)), error = function(e) e)
          }
        }
      } else {
        if (linkF == "cloglog") {
          if (sigmauType == "common") {
          mleList <- tryCatch(switch(udist, hnormal = do.call(cnsfhalfnormAlgOpt_cloglog,
            FunArgs), exponential = do.call(cnsfexponormAlgOpt_cloglog,
            FunArgs), tnormal = do.call(cnsftruncnormAlgOpt_cloglog,
            FunArgs), rayleigh = do.call(cnsfraynormAlgOpt_cloglog,
            FunArgs), gamma = do.call(cnsfgammanormAlgOpt_cloglog,
            FunArgs), uniform = do.call(cnsfuninormAlgOpt_cloglog,
            FunArgs), lognormal = do.call(cnsflognormAlgOpt_cloglog,
            FunArgs), weibull = do.call(cnsfweibullnormAlgOpt_cloglog,
            FunArgs), genexponential = do.call(cnsfgenexponormAlgOpt_cloglog,
            FunArgs), tslaplace = do.call(cnsftslnormAlgOpt_cloglog,
            FunArgs)), error = function(e) e)
          } else {
          if (sigmauType == "different") {
            mleList <- tryCatch(switch(udist, hnormal = do.call(mcesfhalfnormAlgOpt_cloglog,
            FunArgs), exponential = do.call(mcesfexponormAlgOpt_cloglog,
            FunArgs), tnormal = do.call(mcesftruncnormAlgOpt_cloglog,
            FunArgs), rayleigh = do.call(mcesfraynormAlgOpt_cloglog,
            FunArgs), gamma = do.call(mcesfgammanormAlgOpt_cloglog,
            FunArgs), uniform = do.call(mcesfuninormAlgOpt_cloglog,
            FunArgs), lognormal = do.call(mcesflognormAlgOpt_cloglog,
            FunArgs), weibull = do.call(mcesfweibullnormAlgOpt_cloglog,
            FunArgs), genexponential = do.call(mcesfgenexponormAlgOpt_cloglog,
            FunArgs), tslaplace = do.call(mcesftslnormAlgOpt_cloglog,
            FunArgs)), error = function(e) e)
          }
          }
        }
      }
    }
  }
  if (inherits(mleList, "error")) {
    stop("The current error occurs during optimization:\n",
      mleList$message, call. = FALSE)
  }
  # Inverse Hessian + other -------
  mleList$invHessian <- vcovObj(mleObj = mleList$mleObj, hessianType = hessianType,
    method = method, nParm = nParm)
  mleList <- c(mleList, if (method == "ucminf") {
    list(type = "ucminf maximization", nIter = unname(mleList$mleObj$info["neval"]),
      status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$value,
      gradient = mleList$mleObj$gradient)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
      list(type = mleList$mleObj$type, nIter = mleList$mleObj$iterations,
        status = mleList$mleObj$message, mleLoglik = mleList$mleObj$maximum,
        gradient = mleList$mleObj$gradient)
    } else {
      if (method == "sr1") {
        list(type = "SR1 maximization", nIter = mleList$mleObj$iterations,
          status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
          gradient = mleList$mleObj$gradient)
      } else {
        if (method == "mla") {
          list(type = "Levenberg-Marquardt maximization",
          nIter = mleList$mleObj$ni, status = switch(mleList$mleObj$istop,
            `1` = "convergence criteria were satisfied",
            `2` = "maximum number of iterations was reached",
            `4` = "algorithm encountered a problem in the function computation"),
          mleLoglik = -mleList$mleObj$fn.value, gradient = mleList$mleObj$grad)
        } else {
          if (method == "sparse") {
          list(type = "Sparse Hessian maximization",
            nIter = mleList$mleObj$iterations, status = mleList$mleObj$status,
            mleLoglik = -mleList$mleObj$fval, gradient = mleList$mleObj$gradient)
          } else {
          if (method == "nlminb") {
            list(type = "nlminb maximization", nIter = mleList$mleObj$iterations,
            status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$objective,
            gradient = mleList$mleObj$gradient)
          }
          }
        }
      }
    }
  })
  # quick renaming -------
  if (udist %in% c("tnormal", "lognormal")) {
    names(mleList$startVal) <- fName_mu_cnsfcross(Xvar = Xvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Zvar = Zvar,
      sigmauType)
  } else {
    names(mleList$startVal) <- fName_uv_cnsfcross(Xvar = Xvar,
      udist = udist, uHvar = uHvar, vHvar = vHvar, Zvar = Zvar,
      sigmauType)
  }
  names(mleList$mlParam) <- names(mleList$startVal)
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
  names(mleList$gradient) <- names(mleList$mlParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
  # Return object -------
  mleDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M")
  dataTable$mlResiduals <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$mlFitted <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$logL_OBS <- mleList$mleObj$logL_OBS
  returnObj <- list()
  returnObj$call <- cl
  returnObj$formula <- formula
  returnObj$S <- S
  returnObj$typeSfa <- typeSfa
  returnObj$Nobs <- N
  returnObj$nXvar <- nXvar
  if (udist %in% c("tnormal", "lognormal")) {
    returnObj$nmuZUvar <- nmuZUvar
  }
  returnObj$nZHvar <- nZHvar
  returnObj$logDepVar <- logDepVar
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- mleList$startVal
  returnObj$dataTable <- dataTable
  returnObj$sigmauType <- sigmauType
  returnObj$linkF <- linkF
  if (is.null(start)) {
    if (udist == "hnormal") {
      returnObj$initHalf <- mleList$initHalf
    } else {
      if (udist == "exponential") {
        returnObj$initExpo <- mleList$initExpo
      } else {
        if (udist == "tnormal") {
          returnObj$initTrunc <- mleList$initTrunc
        } else {
          if (udist == "rayleigh") {
          returnObj$initRay <- mleList$initRay
          } else {
          if (udist == "uniform") {
            returnObj$initUni <- mleList$initUni
          } else {
            if (udist == "gamma") {
            returnObj$initGamma <- mleList$initGamma
            } else {
            if (udist == "lognormal") {
              returnObj$initLog <- mleList$initLog
            } else {
              if (udist == "weibull") {
              returnObj$initWeibull <- mleList$initWeibull
              } else {
              if (udist == "genexponential") {
                returnObj$initGenExpo <- mleList$initGenExpo
              } else {
                if (udist == "tslaplace") {
                returnObj$initTSL <- mleList$initTSL
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  returnObj$isWeights <- !all.equal(wHvar, rep(1, N))
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
    "Analytic/Numeric Hessian"
  } else {
    if (hessianType == 2) {
      "BHHH Hessian"
    }
  }
  returnObj$mlDate <- mleDate
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    returnObj$simDist <- simDist
    returnObj$Nsim <- Nsim
    returnObj$FiMat <- FiMat
  }
  rm(mleList)
  class(returnObj) <- "cnsfcross"
  return(returnObj)
}

# print for cnsfcross ----------
#' @rdname cnsfcross
#' @export
print.cnsfcross <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call, width.cutoff = 500))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat(cnsfdist(x$udist), "\n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}

# Bread for Sandwich Estimator ----------
#' @rdname cnsfcross
#' @export
bread.cnsfcross <- function(x, ...) {
  if (x$hessianType == "Analytic/Numeric Hessian") {
    return(x$invHessian * x$Nobs)
  } else {
    cat("Computing Analytical/Numerical Hessian \n")
    Yvar <- model.response(model.frame(x$formula, data = x$dataTable))
    Xvar <- model.matrix(x$formula, rhs = 1, data = x$dataTable)
    if (x$udist %in% c("tnormal", "lognormal")) {
      muHvar <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
      uHvar <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
      vHvar <- model.matrix(x$formula, rhs = 4, data = x$dataTable)
      Zvar <- model.matrix(x$formula, rhs = 5, data = x$dataTable)
    } else {
      uHvar <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
      vHvar <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
      Zvar <- model.matrix(x$formula, rhs = 4, data = x$dataTable)
    }
    if (x$linkF == "logit") {
      if (x$sigmauType == "common") {
        if (x$udist == "hnormal") {
          hessAnalytical <- chesscnsfhalfnormlike_logit(parm = x$mlParam,
          nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
          nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, Zvar = Zvar, nZHvar = ncol(Zvar),
          wHvar = x$dataTable$weights, S = x$S)
        } else {
          if (x$udist == "exponential") {
          hessAnalytical <- chesscnsfexponormlike_logit(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, Zvar = Zvar,
            nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
            S = x$S)
          } else {
          if (x$udist == "tnormal") {
            hessAnalytical <- chesscnsftruncnormlike_logit(parm = x$mlParam,
            nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
            nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, Zvar = Zvar,
            nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
            S = x$S)
          } else {
            if (x$udist == "rayleigh") {
            hessAnalytical <- chesscnsfraynormlike_logit(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
              S = x$S)
            } else {
            if (x$udist == "uniform") {
              hessAnalytical <- chesscnsfuninormlike_logit(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar),
              wHvar = x$dataTable$weights, S = x$S)
            } else {
              if (x$udist == "genexponential") {
              hessAnalytical <- chesscnsfgenexponormlike_logit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights, S = x$S)
              } else {
              if (x$udist == "tslaplace") {
                hessAnalytical <- chesscnsftslnormlike_logit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights,
                S = x$S)
              } else {
                if (x$udist == "gamma") {
                hessAnalytical <- jacobian(function(parm) colSums(cgradcnsfgammanormlike_logit(parm,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar, nZHvar = ncol(Zvar),
                  wHvar = x$dataTable$weights,
                  S = x$S, N = x$Nobs, FiMat = x$FiMat)),
                  unname(x$mlParam))
                } else {
                if (x$udist == "weibull") {
                  hessAnalytical <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_logit(parm,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar,
                  nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
                  S = x$S, N = x$Nobs, FiMat = x$FiMat)),
                  unname(x$mlParam))
                } else {
                  if (x$udist == "lognormal") {
                  hessAnalytical <- jacobian(function(parm) colSums(cgradcnsflognormlike_logit(parm,
                    nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
                    nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    muHvar = muHvar, uHvar = uHvar,
                    vHvar = vHvar, Yvar = Yvar,
                    Xvar = Xvar, Zvar = Zvar,
                    nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs, FiMat = x$FiMat)),
                    unname(x$mlParam))
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        }
      } else {
        if (x$sigmauType == "different") {
          if (x$udist == "hnormal") {
          hessAnalytical <- chessmcesfhalfnormlike_logit(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, Zvar = Zvar,
            nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
            S = x$S)
          } else {
          if (x$udist == "exponential") {
            hessAnalytical <- chessmcesfexponormlike_logit(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar,
            vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
            Zvar = Zvar, nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
            S = x$S)
          } else {
            if (x$udist == "tnormal") {
            hessAnalytical <- chessmcesftruncnormlike_logit(parm = x$mlParam,
              nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
              nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
              muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
              Yvar = Yvar, Xvar = Xvar, Zvar = Zvar,
              nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
              S = x$S)
            } else {
            if (x$udist == "rayleigh") {
              hessAnalytical <- chessmcesfraynormlike_logit(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar),
              wHvar = x$dataTable$weights, S = x$S)
            } else {
              if (x$udist == "uniform") {
              hessAnalytical <- chessmcesfuninormlike_logit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights, S = x$S)
              } else {
              if (x$udist == "genexponential") {
                hessAnalytical <- chessmcesfgenexponormlike_logit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights,
                S = x$S)
              } else {
                if (x$udist == "tslaplace") {
                hessAnalytical <- chessmcesftslnormlike_logit(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar, nZHvar = ncol(Zvar),
                  wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                if (x$udist == "gamma") {
                  hessAnalytical <- jacobian(function(parm) colSums(cgradmcesfgammanormlike_logit(parm,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar,
                  nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
                  S = x$S, N = x$Nobs, FiMat = x$FiMat)),
                  unname(x$mlParam))
                } else {
                  if (x$udist == "weibull") {
                  hessAnalytical <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_logit(parm,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs, FiMat = x$FiMat)),
                    unname(x$mlParam))
                  } else {
                  if (x$udist == "lognormal") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradmcesflognormlike_logit(parm,
                    nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
                    nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    muHvar = muHvar, uHvar = uHvar,
                    vHvar = vHvar, Yvar = Yvar,
                    Xvar = Xvar, Zvar = Zvar,
                    nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs,
                    FiMat = x$FiMat)), unname(x$mlParam))
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        }
      }
    } else {
      if (x$linkF == "probit") {
        if (x$sigmauType == "common") {
          if (x$udist == "hnormal") {
          hessAnalytical <- chesscnsfhalfnormlike_probit(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, Zvar = Zvar,
            nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
            S = x$S)
          } else {
          if (x$udist == "exponential") {
            hessAnalytical <- chesscnsfexponormlike_probit(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar,
            vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
            Zvar = Zvar, nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
            S = x$S)
          } else {
            if (x$udist == "tnormal") {
            hessAnalytical <- chesscnsftruncnormlike_probit(parm = x$mlParam,
              nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
              nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
              muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
              Yvar = Yvar, Xvar = Xvar, Zvar = Zvar,
              nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
              S = x$S)
            } else {
            if (x$udist == "rayleigh") {
              hessAnalytical <- chesscnsfraynormlike_probit(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar),
              wHvar = x$dataTable$weights, S = x$S)
            } else {
              if (x$udist == "uniform") {
              hessAnalytical <- chesscnsfuninormlike_probit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights, S = x$S)
              } else {
              if (x$udist == "genexponential") {
                hessAnalytical <- chesscnsfgenexponormlike_probit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights,
                S = x$S)
              } else {
                if (x$udist == "tslaplace") {
                hessAnalytical <- chesscnsftslnormlike_probit(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar, nZHvar = ncol(Zvar),
                  wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                if (x$udist == "gamma") {
                  hessAnalytical <- jacobian(function(parm) colSums(cgradcnsfgammanormlike_probit(parm,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar,
                  nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
                  S = x$S, N = x$Nobs, FiMat = x$FiMat)),
                  unname(x$mlParam))
                } else {
                  if (x$udist == "weibull") {
                  hessAnalytical <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_probit(parm,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs, FiMat = x$FiMat)),
                    unname(x$mlParam))
                  } else {
                  if (x$udist == "lognormal") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradcnsflognormlike_probit(parm,
                    nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
                    nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    muHvar = muHvar, uHvar = uHvar,
                    vHvar = vHvar, Yvar = Yvar,
                    Xvar = Xvar, Zvar = Zvar,
                    nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs,
                    FiMat = x$FiMat)), unname(x$mlParam))
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        } else {
          if (x$sigmauType == "different") {
          if (x$udist == "hnormal") {
            hessAnalytical <- chessmcesfhalfnormlike_probit(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar,
            vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
            Zvar = Zvar, nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
            S = x$S)
          } else {
            if (x$udist == "exponential") {
            hessAnalytical <- chessmcesfexponormlike_probit(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
              S = x$S)
            } else {
            if (x$udist == "tnormal") {
              hessAnalytical <- chessmcesftruncnormlike_probit(parm = x$mlParam,
              nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
              nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
              muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
              Yvar = Yvar, Xvar = Xvar, Zvar = Zvar,
              nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
              S = x$S)
            } else {
              if (x$udist == "rayleigh") {
              hessAnalytical <- chessmcesfraynormlike_probit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights, S = x$S)
              } else {
              if (x$udist == "uniform") {
                hessAnalytical <- chessmcesfuninormlike_probit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights,
                S = x$S)
              } else {
                if (x$udist == "genexponential") {
                hessAnalytical <- chessmcesfgenexponormlike_probit(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar, nZHvar = ncol(Zvar),
                  wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                if (x$udist == "tslaplace") {
                  hessAnalytical <- chessmcesftslnormlike_probit(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar,
                  nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                  if (x$udist == "gamma") {
                  hessAnalytical <- jacobian(function(parm) colSums(cgradmcesfgammanormlike_probit(parm,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs, FiMat = x$FiMat)),
                    unname(x$mlParam))
                  } else {
                  if (x$udist == "weibull") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_probit(parm,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs,
                    FiMat = x$FiMat)), unname(x$mlParam))
                  } else {
                    if (x$udist == "lognormal") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradmcesflognormlike_probit(parm,
                      nXvar = ncol(Xvar),
                      nmuZUvar = ncol(muHvar),
                      nuZUvar = ncol(uHvar),
                      nvZVvar = ncol(vHvar),
                      muHvar = muHvar, uHvar = uHvar,
                      vHvar = vHvar, Yvar = Yvar,
                      Xvar = Xvar, Zvar = Zvar,
                      nZHvar = ncol(Zvar),
                      wHvar = x$dataTable$weights,
                      S = x$S, N = x$Nobs,
                      FiMat = x$FiMat)),
                      unname(x$mlParam))
                    }
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        }
      } else {
        if (x$linkF == "cauchit") {
          if (x$sigmauType == "common") {
          if (x$udist == "hnormal") {
            hessAnalytical <- chesscnsfhalfnormlike_cauchit(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar,
            vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
            Zvar = Zvar, nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
            S = x$S)
          } else {
            if (x$udist == "exponential") {
            hessAnalytical <- chesscnsfexponormlike_cauchit(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
              S = x$S)
            } else {
            if (x$udist == "tnormal") {
              hessAnalytical <- chesscnsftruncnormlike_cauchit(parm = x$mlParam,
              nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
              nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
              muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
              Yvar = Yvar, Xvar = Xvar, Zvar = Zvar,
              nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
              S = x$S)
            } else {
              if (x$udist == "rayleigh") {
              hessAnalytical <- chesscnsfraynormlike_cauchit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights, S = x$S)
              } else {
              if (x$udist == "uniform") {
                hessAnalytical <- chesscnsfuninormlike_cauchit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights,
                S = x$S)
              } else {
                if (x$udist == "genexponential") {
                hessAnalytical <- chesscnsfgenexponormlike_cauchit(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar, nZHvar = ncol(Zvar),
                  wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                if (x$udist == "tslaplace") {
                  hessAnalytical <- chesscnsftslnormlike_cauchit(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar,
                  nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                  if (x$udist == "gamma") {
                  hessAnalytical <- jacobian(function(parm) colSums(cgradcnsfgammanormlike_cauchit(parm,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs, FiMat = x$FiMat)),
                    unname(x$mlParam))
                  } else {
                  if (x$udist == "weibull") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_cauchit(parm,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs,
                    FiMat = x$FiMat)), unname(x$mlParam))
                  } else {
                    if (x$udist == "lognormal") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradcnsflognormlike_cauchit(parm,
                      nXvar = ncol(Xvar),
                      nmuZUvar = ncol(muHvar),
                      nuZUvar = ncol(uHvar),
                      nvZVvar = ncol(vHvar),
                      muHvar = muHvar, uHvar = uHvar,
                      vHvar = vHvar, Yvar = Yvar,
                      Xvar = Xvar, Zvar = Zvar,
                      nZHvar = ncol(Zvar),
                      wHvar = x$dataTable$weights,
                      S = x$S, N = x$Nobs,
                      FiMat = x$FiMat)),
                      unname(x$mlParam))
                    }
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          } else {
          if (x$sigmauType == "different") {
            if (x$udist == "hnormal") {
            hessAnalytical <- chessmcesfhalfnormlike_cauchit(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
              S = x$S)
            } else {
            if (x$udist == "exponential") {
              hessAnalytical <- chessmcesfexponormlike_cauchit(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar),
              wHvar = x$dataTable$weights, S = x$S)
            } else {
              if (x$udist == "tnormal") {
              hessAnalytical <- chessmcesftruncnormlike_cauchit(parm = x$mlParam,
                nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
                nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
                muHvar = muHvar, uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights, S = x$S)
              } else {
              if (x$udist == "rayleigh") {
                hessAnalytical <- chessmcesfraynormlike_cauchit(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights,
                S = x$S)
              } else {
                if (x$udist == "uniform") {
                hessAnalytical <- chessmcesfuninormlike_cauchit(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar, nZHvar = ncol(Zvar),
                  wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                if (x$udist == "genexponential") {
                  hessAnalytical <- chessmcesfgenexponormlike_cauchit(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar,
                  nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                  if (x$udist == "tslaplace") {
                  hessAnalytical <- chessmcesftslnormlike_cauchit(parm = x$mlParam,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S)
                  } else {
                  if (x$udist == "gamma") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradmcesfgammanormlike_cauchit(parm,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs,
                    FiMat = x$FiMat)), unname(x$mlParam))
                  } else {
                    if (x$udist == "weibull") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_cauchit(parm,
                      nXvar = ncol(Xvar),
                      nuZUvar = ncol(uHvar),
                      nvZVvar = ncol(vHvar),
                      uHvar = uHvar, vHvar = vHvar,
                      Yvar = Yvar, Xvar = Xvar,
                      Zvar = Zvar, nZHvar = ncol(Zvar),
                      wHvar = x$dataTable$weights,
                      S = x$S, N = x$Nobs,
                      FiMat = x$FiMat)),
                      unname(x$mlParam))
                    } else {
                    if (x$udist == "lognormal") {
                      hessAnalytical <- jacobian(function(parm) colSums(cgradmcesflognormlike_cauchit(parm,
                      nXvar = ncol(Xvar),
                      nmuZUvar = ncol(muHvar),
                      nuZUvar = ncol(uHvar),
                      nvZVvar = ncol(vHvar),
                      muHvar = muHvar,
                      uHvar = uHvar, vHvar = vHvar,
                      Yvar = Yvar, Xvar = Xvar,
                      Zvar = Zvar, nZHvar = ncol(Zvar),
                      wHvar = x$dataTable$weights,
                      S = x$S, N = x$Nobs,
                      FiMat = x$FiMat)),
                      unname(x$mlParam))
                    }
                    }
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        } else {
          if (x$linkF == "cloglog") {
          if (x$sigmauType == "common") {
            if (x$udist == "hnormal") {
            hessAnalytical <- chesscnsfhalfnormlike_cloglog(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
              S = x$S)
            } else {
            if (x$udist == "exponential") {
              hessAnalytical <- chesscnsfexponormlike_cloglog(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar),
              wHvar = x$dataTable$weights, S = x$S)
            } else {
              if (x$udist == "tnormal") {
              hessAnalytical <- chesscnsftruncnormlike_cloglog(parm = x$mlParam,
                nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
                nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
                muHvar = muHvar, uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights, S = x$S)
              } else {
              if (x$udist == "rayleigh") {
                hessAnalytical <- chesscnsfraynormlike_cloglog(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights,
                S = x$S)
              } else {
                if (x$udist == "uniform") {
                hessAnalytical <- chesscnsfuninormlike_cloglog(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar, nZHvar = ncol(Zvar),
                  wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                if (x$udist == "genexponential") {
                  hessAnalytical <- chesscnsfgenexponormlike_cloglog(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar,
                  nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                  if (x$udist == "tslaplace") {
                  hessAnalytical <- chesscnsftslnormlike_cloglog(parm = x$mlParam,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S)
                  } else {
                  if (x$udist == "gamma") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradcnsfgammanormlike_cloglog(parm,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S, N = x$Nobs,
                    FiMat = x$FiMat)), unname(x$mlParam))
                  } else {
                    if (x$udist == "weibull") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_cloglog(parm,
                      nXvar = ncol(Xvar),
                      nuZUvar = ncol(uHvar),
                      nvZVvar = ncol(vHvar),
                      uHvar = uHvar, vHvar = vHvar,
                      Yvar = Yvar, Xvar = Xvar,
                      Zvar = Zvar, nZHvar = ncol(Zvar),
                      wHvar = x$dataTable$weights,
                      S = x$S, N = x$Nobs,
                      FiMat = x$FiMat)),
                      unname(x$mlParam))
                    } else {
                    if (x$udist == "lognormal") {
                      hessAnalytical <- jacobian(function(parm) colSums(cgradcnsflognormlike_cloglog(parm,
                      nXvar = ncol(Xvar),
                      nmuZUvar = ncol(muHvar),
                      nuZUvar = ncol(uHvar),
                      nvZVvar = ncol(vHvar),
                      muHvar = muHvar,
                      uHvar = uHvar, vHvar = vHvar,
                      Yvar = Yvar, Xvar = Xvar,
                      Zvar = Zvar, nZHvar = ncol(Zvar),
                      wHvar = x$dataTable$weights,
                      S = x$S, N = x$Nobs,
                      FiMat = x$FiMat)),
                      unname(x$mlParam))
                    }
                    }
                  }
                  }
                }
                }
              }
              }
            }
            }
          } else {
            if (x$sigmauType == "different") {
            if (x$udist == "hnormal") {
              hessAnalytical <- chessmcesfhalfnormlike_cloglog(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              Zvar = Zvar, nZHvar = ncol(Zvar),
              wHvar = x$dataTable$weights, S = x$S)
            } else {
              if (x$udist == "exponential") {
              hessAnalytical <- chessmcesfexponormlike_cloglog(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights, S = x$S)
              } else {
              if (x$udist == "tnormal") {
                hessAnalytical <- chessmcesftruncnormlike_cloglog(parm = x$mlParam,
                nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
                nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
                muHvar = muHvar, uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                Zvar = Zvar, nZHvar = ncol(Zvar),
                wHvar = x$dataTable$weights,
                S = x$S)
              } else {
                if (x$udist == "rayleigh") {
                hessAnalytical <- chessmcesfraynormlike_cloglog(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar, nZHvar = ncol(Zvar),
                  wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                if (x$udist == "uniform") {
                  hessAnalytical <- chessmcesfuninormlike_cloglog(parm = x$mlParam,
                  nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                  nvZVvar = ncol(vHvar), uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, Zvar = Zvar,
                  nZHvar = ncol(Zvar), wHvar = x$dataTable$weights,
                  S = x$S)
                } else {
                  if (x$udist == "genexponential") {
                  hessAnalytical <- chessmcesfgenexponormlike_cloglog(parm = x$mlParam,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S)
                  } else {
                  if (x$udist == "tslaplace") {
                    hessAnalytical <- chessmcesftslnormlike_cloglog(parm = x$mlParam,
                    nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                    nvZVvar = ncol(vHvar),
                    uHvar = uHvar, vHvar = vHvar,
                    Yvar = Yvar, Xvar = Xvar,
                    Zvar = Zvar, nZHvar = ncol(Zvar),
                    wHvar = x$dataTable$weights,
                    S = x$S)
                  } else {
                    if (x$udist == "gamma") {
                    hessAnalytical <- jacobian(function(parm) colSums(cgradmcesfgammanormlike_cloglog(parm,
                      nXvar = ncol(Xvar),
                      nuZUvar = ncol(uHvar),
                      nvZVvar = ncol(vHvar),
                      uHvar = uHvar, vHvar = vHvar,
                      Yvar = Yvar, Xvar = Xvar,
                      Zvar = Zvar, nZHvar = ncol(Zvar),
                      wHvar = x$dataTable$weights,
                      S = x$S, N = x$Nobs,
                      FiMat = x$FiMat)),
                      unname(x$mlParam))
                    } else {
                    if (x$udist == "weibull") {
                      hessAnalytical <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_cloglog(parm,
                      nXvar = ncol(Xvar),
                      nuZUvar = ncol(uHvar),
                      nvZVvar = ncol(vHvar),
                      uHvar = uHvar, vHvar = vHvar,
                      Yvar = Yvar, Xvar = Xvar,
                      Zvar = Zvar, nZHvar = ncol(Zvar),
                      wHvar = x$dataTable$weights,
                      S = x$S, N = x$Nobs,
                      FiMat = x$FiMat)),
                      unname(x$mlParam))
                    } else {
                      if (x$udist == "lognormal") {
                      hessAnalytical <- jacobian(function(parm) colSums(cgradmcesflognormlike_cloglog(parm,
                        nXvar = ncol(Xvar),
                        nmuZUvar = ncol(muHvar),
                        nuZUvar = ncol(uHvar),
                        nvZVvar = ncol(vHvar),
                        muHvar = muHvar,
                        uHvar = uHvar,
                        vHvar = vHvar,
                        Yvar = Yvar, Xvar = Xvar,
                        Zvar = Zvar, nZHvar = ncol(Zvar),
                        wHvar = x$dataTable$weights,
                        S = x$S, N = x$Nobs,
                        FiMat = x$FiMat)),
                        unname(x$mlParam))
                      }
                    }
                    }
                  }
                  }
                }
                }
              }
              }
            }
            }
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
#' @rdname cnsfcross
#' @export
estfun.cnsfcross <- function(x, ...) {
  return(x$gradL_OBS)
}
