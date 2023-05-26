################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Model: Stochastic Frontier Analysis                                          #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Stochastic frontier estimation using cross-sectional data
#'
#' @description
#' \code{\link{sfacross}} is a symbolic formula-based function for the
#' estimation of stochastic frontier models in the case of cross-sectional or
#' pooled cross-sectional data, using maximum (simulated) likelihood - M(S)L.
#'
#' The function accounts for heteroscedasticity in both one-sided and two-sided
#' error terms as in Reifschneider and Stevenson (1991), Caudill and Ford
#' (1993), Caudill \emph{et al.} (1995) and Hadri (1999), but also
#' heterogeneity in the mean of the pre-truncated distribution as in Kumbhakar
#' \emph{et al.} (1991), Huang and Liu (1994) and Battese and Coelli (1995).
#'
#' Ten distributions are possible for the one-sided error term and eleven
#' optimization algorithms are available.
#'
#' The truncated normal - normal distribution with scaling property as in Wang
#' and Schmidt (2002) is also implemented.
#'
#' @aliases sfacross print.sfacross
#'
#' @param formula A symbolic description of the model to be estimated based on
#' the generic function \code{formula} (see section \sQuote{Details}).
#' @param muhet A one-part formula to consider heterogeneity in the mean of the
#' pre-truncated distribution (see section \sQuote{Details}).
#' @param uhet A one-part formula to consider heteroscedasticity in the
#' one-sided error variance (see section \sQuote{Details}).
#' @param vhet A one-part formula to consider heteroscedasticity in the
#' two-sided error variance (see section \sQuote{Details}).
#' @param logDepVar Logical. Informs whether the dependent variable is logged
#' (\code{TRUE}) or not (\code{FALSE}). Default = \code{TRUE}.
#' @param data The data frame containing the data.
#' @param subset An optional vector specifying a subset of observations to be
#' used in the optimization process.
#' @param weights An optional vector of weights to be used for weighted 
#' log-likelihood. Should be \code{NULL} or numeric vector with positive values. 
#' When \code{NULL}, a numeric vector of 1 is used.
#' @param wscale Logical. When \code{weights} is not \code{NULL}, a scaling 
#' transformation is used such that the \code{weights} sum to the sample 
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
#' @param scaling Logical. Only when \code{udist = 'tnormal'} and \code{scaling
#' = TRUE}, the scaling property model (Wang and Schmidt 2002) is estimated.
#' Default = \code{FALSE}. (see section \sQuote{Details}).
#' @param start Numeric vector. Optional starting values for the maximum
#' likelihood (ML) estimation.
#' @param method Optimization algorithm used for the estimation. Default =
#' \code{'bfgs'}. 11 algorithms are available: \itemize{ \item \code{'bfgs'},
#' for Broyden-Fletcher-Goldfarb-Shanno (see
#' \code{\link[maxLik:maxBFGS]{maxBFGS}}) \item \code{'bhhh'}, for
#' Berndt-Hall-Hall-Hausman (see \code{\link[maxLik:maxBHHH]{maxBHHH}}) \item
#' \code{'nr'}, for Newton-Raphson (see \code{\link[maxLik:maxNR]{maxNR}})
#' \item \code{'nm'}, for Nelder-Mead (see \code{\link[maxLik:maxNM]{maxNM}})
#' \item \code{'cg'}, for Conjugate Gradient 
#' (see \code{\link[maxLik:maxCG]{maxCG}}) \item \code{'sann'}, for Simulated 
#' Annealing (see \code{\link[maxLik:maxSANN]{maxSANN}}) \item \code{'ucminf'}, 
#' for a quasi-Newton type optimisation with BFGS updating of the inverse Hessian and 
#' soft line search with a trust region type monitoring of the input to the line 
#' search algorithm (see \code{\link[ucminf:ucminf]{ucminf}})
#' \item \code{'mla'}, for general-purpose optimization based on
#' Marquardt-Levenberg algorithm (see \code{\link[marqLevAlg:mla]{mla}})
#' \item \code{'sr1'}, for Symmetric Rank 1 (see
#' \code{\link[trustOptim:trust.optim]{trust.optim}}) \item \code{'sparse'}, 
#' for trust regions and sparse Hessian 
#' (see \code{\link[trustOptim:trust.optim]{trust.optim}}) \item
#' \code{'nlminb'}, for optimization using PORT routines (see
#' \code{\link[stats:nlminb]{nlminb}})}
#' @param hessianType Integer. If \code{1} (Default), analytic Hessian is
#' returned for all the distributions. If \code{2},
#' bhhh Hessian is estimated (\eqn{g'g}).
#' @param simType Character string. If \code{simType = 'halton'} (Default),
#' Halton draws are used for maximum simulated likelihood (MSL). If
#' \code{simType = 'ghalton'}, Generalized-Halton draws are used for MSL. If
#' \code{simType = 'sobol'}, Sobol draws are used for MSL. If \code{simType =
#' 'uniform'}, uniform draws are used for MSL. (see section \sQuote{Details}).
#' @param Nsim Number of draws for MSL. Default 100.
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
#' @param x an object of class sfacross (returned by the function 
#' \code{\link{sfacross}}).
#' @param ... additional arguments of frontier are passed to sfacross; 
#' additional arguments of the print, bread, estfun, nobs methods are currently
#'  ignored.
#'
#' @details
#' The stochastic frontier model for the cross-sectional data is defined as: 
#' 
#' \deqn{y_i = \alpha + \mathbf{x_i^{\prime}}\bm{\beta} + v_i - Su_i}
#' 
#' with
#' 
#' \deqn{\epsilon_i = v_i -Su_i}
#'
#' where \eqn{i} is the observation, \eqn{y} is the
#' output (cost, revenue, profit), \eqn{\mathbf{x}} is the vector of main explanatory
#' variables (inputs and other control variables), \eqn{u} is the one-sided
#' error term with variance \eqn{\sigma_{u}^2}, and \eqn{v} is the two-sided
#' error term with variance \eqn{\sigma_{v}^2}.
#'
#' \code{S = 1} in the case of production (profit) frontier function and
#' \code{S = -1} in the case of cost frontier function.
#'
#' The model is estimated using maximum likelihood (ML) for most distributions
#' except the Gamma, Weibull and log-normal distributions for which maximum
#' simulated likelihood (MSL) is used. For this latter, several draws can be
#' implemented namely Halton, Generalized Halton, Sobol and uniform. In the
#' case of uniform draws, antithetics can also be computed: first \code{Nsim/2}
#' draws are obtained, then the \code{Nsim/2} other draws are obtained as
#' counterpart of one (\code{1-draw}).
#'
#' To account for heteroscedasticity in the variance parameters of the error
#' terms, a single part (right) formula can also be specified. To impose the
#' positivity to these parameters, the variances are modelled as:
#' \eqn{\sigma^2_u = \exp{(\bm{\delta}'\mathbf{Z}_u)}} or \eqn{\sigma^2_v =
#' \exp{(\bm{\phi}'\mathbf{Z}_v)}}, where \eqn{\mathbf{Z}_u} and \eqn{\mathbf{Z}_v} are the heteroscedasticity
#' variables (inefficiency drivers in the case of \eqn{\mathbf{Z}_u}) and \eqn{\bm{\delta}}
#' and \eqn{\bm{\phi}} the coefficients. In the case of heterogeneity in the
#' truncated mean \eqn{\mu}, it is modelled as \eqn{\mu=\bm{\omega}'\mathbf{Z}_{\mu}}. The
#' scaling property can be applied for the truncated normal distribution:
#' \eqn{u \sim h(\mathbf{Z}_u, \delta)u} where \eqn{u} follows a truncated normal
#' distribution \eqn{N^+(\tau, \exp{(cu)})}.
#'
#' In the case of the truncated normal distribution, the convolution of
#' \eqn{u_i} and \eqn{v_i} is:
#' 
#' \deqn{f(\epsilon_i)=\frac{1}{\sqrt{\sigma_u^2 + 
#' \sigma_v^2}}\phi\left(\frac{S\epsilon_i + \mu}{\sqrt{
#' \sigma_u^2 + \sigma_v^2}}\right)\Phi\left(\frac{
#' \mu_{i*}}{\sigma_*}\right)\Big/\Phi\left(\frac{
#' \mu}{\sigma_u}\right)}
#'
#' where
#' 
#' \deqn{\mu_{i*}=\frac{\mu\\\sigma_v^2 - 
#' S\epsilon_i\sigma_u^2}{\sigma_u^2 + \sigma_v^2}}
#'
#' and
#' 
#' \deqn{\sigma_*^2 = \frac{\sigma_u^2 
#' \sigma_v^2}{\sigma_u^2 + \sigma_v^2}}
#' 
#' In the case of the half normal distribution the convolution is obtained by
#' setting \eqn{\mu=0}.
#' 
#' \code{sfacross} allows for the maximization of weighted log-likelihood.
#' When option \code{weights} is specified and \code{wscale = TRUE}, the weights
#' are scaled as: 
#' 
#' \deqn{new_{weights} = sample_{size} \times 
#' \frac{old_{weights}}{\sum(old_{weights})}}
#' 
#' For complex problems, non-gradient methods (e.g. \code{nm} or \code{sann}) 
#' can be used to warm start the optimization and zoom in the neighborhood of 
#' the solution. Then a gradient-based methods is recommended in the second 
#' step. In the case of \code{sann}, we recommend to significantly increase the
#'  iteration limit (e.g. \code{itermax = 20000}). The Conjugate Gradient 
#'  (\code{cg}) can also be used in the first stage.
#' 
#' A set of extractor functions for fitted model objects is available for 
#' objects of class \code{'sfacross'} including methods to the generic functions 
#' \code{\link[=print.sfacross]{print}},
#' \code{\link[=summary.sfacross]{summary}}, \code{\link[=coef.sfacross]{coef}}, 
#' \code{\link[=fitted.sfacross]{fitted}},
#'  \code{\link[=logLik.sfacross]{logLik}}, 
#' \code{\link[=residuals.sfacross]{residuals}}, 
#' \code{\link[=vcov.sfacross]{vcov}}, 
#' \code{\link[=efficiencies.sfacross]{efficiencies}}, 
#' \code{\link[=ic.sfacross]{ic}}, 
#' \code{\link[=marginal.sfacross]{marginal}}, 
#' \code{\link[=skewnessTest]{skewnessTest}}, 
#' \code{\link[=estfun.sfacross]{estfun}} and 
#' \code{\link[=bread.sfacross]{bread}} (from the \CRANpkg{sandwich} package), 
#' [lmtest::coeftest()] (from the \CRANpkg{lmtest} package).
#'
#' @return \code{\link{sfacross}} returns a list of class \code{'sfacross'}
#' containing the following elements:
#'
#' \item{call}{The matched call.}
#'
#' \item{formula}{The estimated model.}
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
#' \item{nmuZUvar}{Number of variables explaining heterogeneity in the
#' truncated mean, only if \code{udist = 'tnormal'} or \code{'lognormal'}.}
#'
#' \item{scaling}{The argument \code{'scaling'}. See the section
#' \sQuote{Arguments}.}
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
#' \item{olsParam}{Numeric vector. OLS estimates.}
#'
#' \item{olsStder}{Numeric vector. Standard errors of OLS estimates.}
#'
#' \item{olsSigmasq}{Numeric. Estimated variance of OLS random error.}
#'
#' \item{olsLoglik}{Numeric. Log-likelihood value of OLS estimation.}
#'
#' \item{olsSkew}{Numeric. Skewness of the residuals of the OLS estimation.}
#'
#' \item{olsM3Okay}{Logical. Indicating whether the residuals of the OLS
#' estimation have the expected skewness.}
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
#' @seealso \code{\link[=print.sfacross]{print}} for printing \code{sfacross} 
#' object.
#' 
#' \code{\link[=summary.sfacross]{summary}} for creating and printing
#' summary results.
#'
#' \code{\link[=coef.sfacross]{coef}} for extracting coefficients of the
#' estimation.
#'
#' \code{\link[=efficiencies.sfacross]{efficiencies}} for computing
#' (in-)efficiency estimates.
#'
#' \code{\link[=fitted.sfacross]{fitted}} for extracting the fitted frontier
#' values.
#'
#' \code{\link[=ic.sfacross]{ic}} for extracting information criteria.
#'
#' \code{\link[=logLik.sfacross]{logLik}} for extracting log-likelihood
#' value(s) of the estimation.
#'
#' \code{\link[=marginal.sfacross]{marginal}} for computing marginal effects of
#' inefficiency drivers.
#'
#' \code{\link[=residuals.sfacross]{residuals}} for extracting residuals of the
#' estimation.
#' 
#' \code{\link[=skewnessTest]{skewnessTest}} for conducting residuals
#' skewness test.
#'
#' \code{\link[=vcov.sfacross]{vcov}} for computing the variance-covariance
#' matrix of the coefficients.
#' 
#' \code{\link[=bread.sfacross]{bread}} for bread for sandwich estimator.
#' 
#' \code{\link[=estfun.sfacross]{estfun}} for gradient extraction for each 
#' observation.
#'
#' \code{\link{skewnessTest}} for implementing skewness test.
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
#' Coelli, T. 1995. Estimators and hypothesis tests for a stochastic frontier
#' function - a Monte-Carlo analysis. \emph{Journal of Productivity Analysis},
#' \bold{6}:247--268.
#'
#' D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
#' Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
#' \emph{Biometrika}, \bold{60}:613--622.
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
#' Kumbhakar, S. C., Ghosh, S., and McGuckin, J. T. 1991) A generalized
#' production frontier approach for estimating determinants of inefficiency in
#' U.S. dairy farms. \emph{Journal of Business & Economic Statistics},
#' \bold{9}(3), 279--286.
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
#' Wang, H.J., and Schmidt, P. 2002. One-step and two-step estimation of the
#' effects of exogenous variables on technical efficiency levels. \emph{Journal
#' of Productivity Analysis}, \bold{18}:129--144.
#'
#' Wang, J. 2012. A normal truncated skewed-Laplace model in stochastic
#' frontier analysis. Master thesis, Western Kentucky University, May.
#'
#' @keywords models optimize cross-section likelihood
#'
#' @examples
#'
#' ## Using data on fossil fuel fired steam electric power generation plants in the U.S.
#' # Translog (cost function) half normal with heteroscedasticity
#' tl_u_h <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'hnormal', uhet = ~ regu, data = utility, S = -1, method = 'bfgs')
#' summary(tl_u_h)
#'
#' # Translog (cost function) truncated normal with heteroscedasticity
#' tl_u_t <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, data = utility, S = -1, method = 'bhhh')
#' summary(tl_u_t)
#'
#' # Translog (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = 'mla')
#' summary(tl_u_ts)
#'
#' ## Using data on Philippine rice producers
#' # Cobb Douglas (production function) generalized exponential, and Weibull 
#' # distributions
#'
#' cb_p_ge <- sfacross(formula = log(PROD) ~ log(AREA) + log(LABOR) + log(NPK) +
#' log(OTHER), udist = 'genexponential', data = ricephil, S = 1, method = 'bfgs')
#' summary(cb_p_ge)
#'
#' ## Using data on U.S. electric utility industry
#' # Cost frontier Gamma distribution
#' tl_u_g <- sfacross(formula = log(cost/fprice) ~ log(output) + I(log(output)^2) +
#' I(log(lprice/fprice)) + I(log(cprice/fprice)), udist = 'gamma', uhet = ~ 1,
#' data = electricity, S = -1, method = 'bfgs', simType = 'halton', Nsim = 200,
#' hessianType = 2)
#' summary(tl_u_g)
#'
#' @export
sfacross <- function(formula, muhet, uhet, vhet, logDepVar = TRUE,
  data, subset, weights, wscale = TRUE, S = 1L, udist = "hnormal",
  scaling = FALSE, start = NULL, method = "bfgs", hessianType = 1L,
  simType = "halton", Nsim = 100, prime = 2L, burn = 10, antithetics = FALSE,
  seed = 12345, itermax = 2000, printInfo = FALSE, tol = 1e-12,
  gradtol = 1e-06, stepmax = 0.1, qac = "marquardt") {
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
  formula <- interCheckMain(formula = formula, data = data)
  if (!missing(muhet)) {
    muhet <- clhsCheck_mu(formula = muhet)
  } else {
    muhet <- ~1
  }
  if (!missing(uhet)) {
    uhet <- clhsCheck_u(formula = uhet)
  } else {
    uhet <- ~1
  }
  if (!missing(vhet)) {
    vhet <- clhsCheck_v(formula = vhet)
  } else {
    vhet <- ~1
  }
  formula <- formDist_sfacross(udist = udist, formula = formula,
    muhet = muhet, uhet = uhet, vhet = vhet)
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
  # if subset is non-missing and there are NA, force data to
  # change
  data <- data[row.names(data) %in% attr(mc, "row.names"),
    ]
  data <- data[validObs, ]
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
  } else {
    mtuH <- delete.response(terms(formula, data = data, rhs = 2))
    uHvar <- model.matrix(mtuH, mc)
    uHvar <- uHvar[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar)
    mtvH <- delete.response(terms(formula, data = data, rhs = 3))
    vHvar <- model.matrix(mtvH, mc)
    vHvar <- vHvar[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar)
  }
  # Check other supplied options -------
  if (length(S) != 1 || !(S %in% c(-1L, 1L))) {
    stop("argument 'S' must equal either 1 or -1: 1 for production or profit 
    frontier
   and -1 for cost frontier",
      call. = FALSE)
  }
  typeSfa <- if (S == 1L) {
    "Stochastic Production/Profit Frontier, e = v - u"
  } else {
    "Stochastic Cost Frontier, e = v + u"
  }
  if (length(scaling) != 1 || !is.logical(scaling[1])) {
    stop("argument 'scaling' must be a single logical value",
      call. = FALSE)
  }
  if (scaling) {
    if (udist != "tnormal") {
      stop("argument 'udist' must be 'tnormal' when scaling option is TRUE",
        call. = FALSE)
    }
    if (nuZUvar != nmuZUvar) {
      stop("argument 'muhet' and 'uhet' must have the same length",
        call. = FALSE)
    }
    if (!all(colnames(uHvar) == colnames(muHvar))) {
      stop("argument 'muhet' and 'uhet' must contain the same variables",
        call. = FALSE)
    }
    if (nuZUvar == 1 || nmuZUvar == 1) {
      if (attr(terms(muhet), "intercept") == 1 || attr(terms(uhet),
        "intercept") == 1) {
        stop("at least one exogeneous variable must be provided for the scaling 
             option",
          call. = FALSE)
      }
    }
  }
  if (length(logDepVar) != 1 || !is.logical(logDepVar[1])) {
    stop("argument 'logDepVar' must be a single logical value",
      call. = FALSE)
  }
  # Number of parameters -------
  nParm <- if (udist == "tnormal") {
    if (scaling) {
      if (attr(terms(muhet), "intercept") == 1 || attr(terms(uhet),
        "intercept") == 1) {
        nXvar + (nmuZUvar - 1) + 2 + nvZVvar
      } else {
        nXvar + nmuZUvar + 2 + nvZVvar
      }
    } else {
      nXvar + nmuZUvar + nuZUvar + nvZVvar
    }
  } else {
    if (udist == "lognormal") {
      nXvar + nmuZUvar + nuZUvar + nvZVvar
    } else {
      if (udist %in% c("gamma", "weibull", "tslaplace")) {
        nXvar + nuZUvar + nvZVvar + 1
      } else {
        nXvar + nuZUvar + nvZVvar
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
    # burn + 1 cause halton fn always starts with 0
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
      lm(Yvar ~ ., data = as.data.frame(Xvar[, -1, drop = FALSE]),
        weights = wHvar)
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
  names(olsRes$coefficients) <- colnames(Xvar)
  olsParam <- c(olsRes$coefficients)
  olsSigmasq <- summary(olsRes)$sigma^2
  olsStder <- sqrt(diag(vcov(olsRes)))
  olsLoglik <- logLik(olsRes)[1]
  if (inherits(data, "pdata.frame")) {
    dataTable <- data[, names(index(data))]
  } else {
    dataTable <- data.frame(IdObs = 1:sum(validObs))
  }
  dataTable <- cbind(dataTable, data[, all.vars(terms(formula))],
    weights = wHvar)
  dataTable <- cbind(dataTable, olsResiduals = residuals(olsRes),
    olsFitted = fitted(olsRes))
  # possibility to have duplicated columns if ID or TIME
  # appears in ols in the case of panel data
  dataTable <- dataTable[!duplicated(as.list(dataTable))]
  olsSkew <- skewness(dataTable[["olsResiduals"]])
  olsM3Okay <- if (S * olsSkew < 0) {
    "Residuals have the expected skeweness"
  } else {
    "Residuals do not have the expected skeweness"
  }
  if (S * olsSkew > 0) {
    warning("The residuals of the OLS are", if (S == 1) {
      " right"
    } else {
      " left"
    }, "-skewed. This may indicate the absence of inefficiency or
  model misspecification or sample 'bad luck'",
      call. = FALSE)
  }
  m2 <- mean((dataTable[["olsResiduals"]] - mean(dataTable[["olsResiduals"]]))^2)
  m3 <- mean((dataTable[["olsResiduals"]] - mean(dataTable[["olsResiduals"]]))^3)
  CoelliM3Test <- c(z = m3/sqrt(6 * m2^3/N), p.value = 2 *
    pnorm(-abs(m3/sqrt(6 * m2^3/N))))
  AgostinoTest <- dagoTest(dataTable[["olsResiduals"]])
  class(AgostinoTest) <- "dagoTest"
  # Step 2: MLE arguments -------
  FunArgs <- if (udist == "tnormal") {
    if (scaling) {
      list(start = start, olsParam = olsParam, dataTable = dataTable,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, method = method, printInfo = printInfo,
        itermax = itermax, stepmax = stepmax, tol = tol,
        gradtol = gradtol, hessianType = hessianType,
        qac = qac)
    } else {
      list(start = start, olsParam = olsParam, dataTable = dataTable,
        nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, method = method, printInfo = printInfo,
        itermax = itermax, stepmax = stepmax, tol = tol,
        gradtol = gradtol, hessianType = hessianType,
        qac = qac)
    }
  } else {
    if (udist == "lognormal") {
      list(start = start, olsParam = olsParam, dataTable = dataTable,
        nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, method = method,
        printInfo = printInfo, itermax = itermax, stepmax = stepmax,
        tol = tol, gradtol = gradtol, hessianType = hessianType,
        qac = qac)
    } else {
      if (udist %in% c("gamma", "weibull")) {
        list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
          method = method, printInfo = printInfo, itermax = itermax,
          stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, qac = qac)
      } else {
        list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, S = S, wHvar = wHvar, method = method,
          printInfo = printInfo, itermax = itermax, stepmax = stepmax,
          tol = tol, gradtol = gradtol, hessianType = hessianType,
          qac = qac)
      }
    }
  }
  ## MLE run -------
  mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt,
    FunArgs), exponential = do.call(exponormAlgOpt, FunArgs),
    tnormal = if (scaling) {
      do.call(truncnormscalAlgOpt, FunArgs)
    } else {
      do.call(truncnormAlgOpt, FunArgs)
    }, rayleigh = do.call(raynormAlgOpt, FunArgs), gamma = do.call(gammanormAlgOpt,
      FunArgs), uniform = do.call(uninormAlgOpt, FunArgs),
    lognormal = do.call(lognormAlgOpt, FunArgs), weibull = do.call(weibullnormAlgOpt,
      FunArgs), genexponential = do.call(genexponormAlgOpt,
      FunArgs), tslaplace = do.call(tslnormAlgOpt, FunArgs)),
    error = function(e) print(e))
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
          gradient = -mleList$mleObj$gradient)
      } else {
        if (method == "mla") {
          list(type = "Lev. Marquardt max.", nIter = mleList$mleObj$ni,
          status = switch(mleList$mleObj$istop, `1` = "convergence criteria were satisfied",
            `2` = "maximum number of iterations was reached",
            `4` = "algorithm encountered a problem in the function computation"),
          mleLoglik = -mleList$mleObj$fn.value, gradient = -mleList$mleObj$grad)
        } else {
          if (method == "sparse") {
          list(type = "Sparse Hessian max.", nIter = mleList$mleObj$iterations,
            status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
            gradient = -mleList$mleObj$gradient)
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
  # quick renaming (only when start is provided) -------
  if (!is.null(start)) {
    if (udist %in% c("tnormal", "lognormal")) {
      names(mleList$startVal) <- fName_mu_sfacross(Xvar = Xvar,
        udist = udist, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, scaling = scaling)
    } else {
      names(mleList$startVal) <- fName_uv_sfacross(Xvar = Xvar,
        udist = udist, uHvar = uHvar, vHvar = vHvar)
    }
    names(mleList$mlParam) <- names(mleList$startVal)
  }
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
  names(mleList$gradient) <- names(mleList$mlParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
  # Return object -------
  mlDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M")
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
  returnObj$scaling <- scaling
  returnObj$logDepVar <- logDepVar
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- mleList$startVal
  returnObj$dataTable <- dataTable
  returnObj$olsParam <- olsParam
  returnObj$olsStder <- olsStder
  returnObj$olsSigmasq <- olsSigmasq
  returnObj$olsLoglik <- olsLoglik
  returnObj$olsSkew <- olsSkew
  returnObj$olsM3Okay <- olsM3Okay
  returnObj$CoelliM3Test <- CoelliM3Test
  returnObj$AgostinoTest <- AgostinoTest
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
    "Analytic Hessian"
  } else {
    if (hessianType == 2) {
      "BHHH Hessian"
    }
  }
  returnObj$mlDate <- mlDate
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    returnObj$simDist <- simDist
    returnObj$Nsim <- Nsim
    returnObj$FiMat <- FiMat
  }
  rm(mleList)
  class(returnObj) <- "sfacross"
  return(returnObj)
}

# print for sfacross ----------
#' @rdname sfacross
#' @export
print.sfacross <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat(sfacrossdist(x$udist), "\n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}

# Bread for Sandwich Estimator ----------
#' @rdname sfacross
#' @export
bread.sfacross <- function(x, ...) {
  if (x$hessianType == "Analytic Hessian") {
    return(x$invHessian * x$Nobs)
  } else {
    cat("Computing Analytical Hessian \n")
    Yvar <- model.response(model.frame(x$formula, data = x$dataTable))
    Xvar <- model.matrix(x$formula, rhs = 1, data = x$dataTable)
    if (x$udist %in% c("tnormal", "lognormal")) {
      muHvar <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
      uHvar <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
      vHvar <- model.matrix(x$formula, rhs = 4, data = x$dataTable)
    } else {
      uHvar <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
      vHvar <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
    }
    if (x$udist == "hnormal") {
      hessAnalytical <- chesshalfnormlike(parm = x$mlParam,
        nXvar = ncol(Xvar), nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = x$dataTable$weights, S = x$S)
    } else {
      if (x$udist == "exponential") {
        hessAnalytical <- chessexponormlike(parm = x$mlParam,
          nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
          nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
          S = x$S)
      } else {
        if (x$udist == "tnormal") {
          if (x$scaling == TRUE) {
          hessAnalytical <- chesstruncnormscalike(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
            S = x$S)
          } else {
          hessAnalytical <- chesstruncnormlike(parm = x$mlParam,
            nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
            nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
            S = x$S)
          }
        } else {
          if (x$udist == "rayleigh") {
          hessAnalytical <- chessraynormlike(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
            S = x$S)
          } else {
          if (x$udist == "uniform") {
            hessAnalytical <- chessuninormlike(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar,
            vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
            wHvar = x$dataTable$weights, S = x$S)
          } else {
            if (x$udist == "genexponential") {
            hessAnalytical <- chessgenexponormlike(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              wHvar = x$dataTable$weights, S = x$S)
            } else {
            if (x$udist == "tslaplace") {
              hessAnalytical <- chesstslnormlike(parm = x$mlParam,
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
              nvZVvar = ncol(vHvar), uHvar = uHvar,
              vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              wHvar = x$dataTable$weights, S = x$S)
            } else {
              if (x$udist == "gamma") {
              hessAnalytical <- chessgammanormlike(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                wHvar = x$dataTable$weights, S = x$S,
                N = x$Nobs, FiMat = x$FiMat)
              } else {
              if (x$udist == "weibull") {
                hessAnalytical <- chessweibullnormlike(parm = x$mlParam,
                nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar,
                vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                wHvar = x$dataTable$weights,
                S = x$S, N = x$Nobs, FiMat = x$FiMat)
              } else {
                if (x$udist == "lognormal") {
                hessAnalytical <- chesslognormlike(parm = x$mlParam,
                  nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
                  nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
                  muHvar = muHvar, uHvar = uHvar,
                  vHvar = vHvar, Yvar = Yvar,
                  Xvar = Xvar, wHvar = x$dataTable$weights,
                  S = x$S, N = x$Nobs, FiMat = x$FiMat)
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
#' @rdname sfacross
#' @export
estfun.sfacross <- function(x, ...) {
  return(x$gradL_OBS)
}
