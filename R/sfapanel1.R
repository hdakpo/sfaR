################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Model: Stochastic Frontier Analysis                                          #
# First generation panel models                                                #
# Data: Panel data                                                             #
#------------------------------------------------------------------------------#

#' Stochastic frontier estimation using panel data
#' 
#' @description
#' \code{\link{sfapanel1}} is a symbolic formula-based function for the
#' estimation of stochastic frontier models in the case of panel data, using 
#' maximum (simulated) likelihood - M(S)L. Several panel models are implemented:
#' the time invariant inefficiency model discussed in Pitt and Lee (1981), the 
#' time varying efficiency models suggested in Kumbhakar (1990), 
#' Battese and Coelli (1992), Cuesta (2000), Cuesta and Orea (2002), 
#' Kumbhakar and Wang (2005), Alvarez \emph{et al.} (2006), 
#' Feng and Serletis (2009). We also suggested a modified version of the model
#' by Lee and Schmidt (1993) - see \sQuote{Details} section.
#'
#' Depending on the specification, the function accounts for heteroscedasticity 
#' in both one-sided and two-sided error terms as in 
#' Reifschneider and Stevenson (1991), Caudill and Ford (1993), 
#' Caudill \emph{et al.} (1995) and Hadri (1999), but also
#' heterogeneity in the mean of the pre-truncated distribution as in Kumbhakar
#' \emph{et al.} (1991), Huang and Liu (1994) and Battese and Coelli (1995).
#' Alvarez \emph{et al.} (2006) implements a version of the time varying 
#' inefficiency using the scaling property as in Wang and Schmidt (2002).
#'
#' Ten distributions are possible for the one-sided error term and eleven
#' optimization algorithms are available.
#'
#' @aliases sfapanel1 bread.sfapanel1 estfun.sfapanel1 print.sfapanel1
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
#' @param data The data frame containing the data. The data can be an object of
#' class \code{'pdata.frame'}.
#' @param idVar Character. When \code{'data'} is not an object of class 
#' \code{'pdata.frame'}, \code{'idVar'} defines the cross sections 
#' identification. Default = \code{NULL}.
#' @param timeVar Character. When \code{'data'} is not an object of class 
#' \code{'pdata.frame'}, \code{'timeVar'} defines the periods identification.
#' Default = \code{NULL}.
#' @param subset An optional vector specifying a subset of observations to be
#' used in the optimization process.
#' @param weights An optional vector of weights to be used for weighted 
#' log-likelihood. Should be \code{NULL} or numeric vector with positive values. 
#' When \code{NULL}, a numeric vector of 1 is used.
#' @param wscale Logical. When \code{weights} is not \code{NULL}, a scaling 
#' transformation is used such that the the \code{weights} sum to the sample 
#' size. Default = \code{TRUE}. When \code{FALSE} no scaling is used.
#' @param S If \code{S = 1} (default), a production (profit) frontier is
#' estimated. If \code{S = -1}, a cost frontier is estimated.
#' @param modelType Character string specifying the panel model to be estimated.
#' See \sQuote{Details} section.
#' \itemize{
#' \item Default = \code{'bc92a'} for the model discussed in Battese and Coelli 
#' (1992).
#' \item A modified version (\code{'mbc92'}) of the previous model was also 
#' presented in Battese and Coelli (1992).
#' \item \code{'bc92b'} implements the model presented in Cuesta and Orea (2002), 
#' and Feng and Serletis (2009).
#' \item \code{'bc92c'} for the model in Alvarez \emph{et al.} (2006).
#' \item \code{'c00'} for the model in Cuesta (2000).
#' \item \code{'kw05'} for the model in Kumbhakar and Wang (2005).
#' \item \code{'mols93'} for the modified version of the model in 
#' Lee and Schmidt (1993).
#' \item \code{'pl81'} for the time invariant inefficiency in
#' Pitt and Lee (1981).
#' }
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
#' @param start Numeric vector. Optional starting values for the maximum
#' likelihood (ML) estimation.
#' @param randStart Logical. Define if random starting values should be used for
#' M(S)L estimation. New starting values are obtained as old ones + draws from
#' normal distribution with std. deviation of 0.01. \code{'seed'} is not 
#' considered here, then each run will provide different starting values 
#' (unless a seed is set by the user before the run).
#' @param whichStart Integer. If \code{'whichStart = 1'}, the starting values 
#' are obtained from the method of moments applied to the pooled data. 
#' When \code{'whichStart = 2'} (Default), the model is initialized by solving 
#' the classic SFA model for the homoscedastic pooled cross section data. 
#' \code{'whichStart = 1'} can be fast especially in the case of maximum 
#' simulated likelihood.
#' @param initAlg Character string specifying the algorithm used for 
#' initializing and obtaining the starting values (when \code{'whichStart = 2'}).
#' Only \pkg{maxLik} package algorithms are available: 
#' \itemize{ \item \code{'bfgs'}, for Broyden-Fletcher-Goldfarb-Shanno 
#' (see \code{\link[maxLik:maxBFGS]{maxBFGS}})
#'  \item \code{'bhhh'}, for Berndt-Hall-Hall-Hausman 
#'  (see \code{\link[maxLik:maxBHHH]{maxBHHH}}) 
#'  \item \code{'nr'}, for Newton-Raphson (see \code{\link[maxLik:maxNR]{maxNR}})
#' \item \code{'nm'}, for Nelder-Mead - Default - 
#'  (see \code{\link[maxLik:maxNM]{maxNM}})
#' \item \code{'cg'}, for Conjugate Gradient 
#' (see \code{\link[maxLik:maxCG]{maxCG}}) \item \code{'sann'}, for Simulated 
#' Annealing (see \code{\link[maxLik:maxSANN]{maxSANN}})
#' }
#' @param initIter Maximum number of iterations for initializing algorithm.
#' Default \code{500}.
#' @param invariance Integer. Except in the case of \code{'modelType = 'bc92c''},
#' in the presence of heteroscedasticity, variables are supposed to be constant
#' over time for each cross section. When this is not the case, three options 
#' are suggested to create constant variables: when \code{'invariance = 1L'}, 
#' the first period observation is used; when \code{'invariance = 2L'}, the last
#' period observation is used; and when \code{'invariance = 3L'}, the period
#' mean is used. The \code{'invariance'} option also applies to the 
#' \code{'weights'} variable.
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
#' for a quasi-Newton type optimization with BFGS updating of the inverse 
#' Hessian and soft line search with a trust region type monitoring of the input 
#' to the line search algorithm (see \code{\link[ucminf:ucminf]{ucminf}})
#' \item \code{'mla'}, for general-purpose optimization based on
#' Marquardt-Levenberg algorithm (see \code{\link[marqLevAlg:mla]{mla}})
#' \item \code{'sr1'}, for Symmetric Rank 1 (see
#' \code{\link[trustOptim:trust.optim]{trust.optim}}) \item \code{'sparse'}, 
#' for trust regions and sparse Hessian 
#' (see \code{\link[trustOptim:trust.optim]{trust.optim}}) \item
#' \code{'nlminb'}, for optimization using PORT routines (see
#' \code{\link[stats:nlminb]{nlminb}})}
#' @param hessianType Integer. If \code{1} (Default), analytic/numeric Hessian 
#' is returned for all the distributions. If \code{2}, bhhh Hessian is 
#' estimated (\eqn{g'g}).
#' @param simType Character string. If \code{simType = 'halton'} (Default),
#' Halton draws are used for maximum simulated likelihood (MSL). If
#' \code{simType = 'ghalton'}, Generalized-Halton  draws are used for MSL. If
#' \code{simType = 'sobol'} or \code{simType = 'rsobol'}, Sobol draws or 
#' randomized Sobol draws are used for MSL, respectively. If 
#' \code{simType = 'richtmyer'}, or \code{simType = 'rrichtmyer'}, 
#' Richtmyer sequence or randomized Richtmyer sequence is used for MSL 
#' estimation, respectively. If \code{simType = 'mlhs'}, modified latin 
#' hypercube sampling is performed. If \code{simType = 'uniform'}, uniform draws 
#' are used for MSL. (see section \sQuote{Details}).
#' @param Nsim Number of draws for MSL. Default 100.
#' @param prime Prime number considered for quasi-random sequences. The first
#' 10,000 primes are available. Default = \code{2}.
#' @param burn Number of the first observations discarded in the case of Halton,
#' Richtmyer, and Sobol draws. For the randomized versions, the burn is not
#' considered. Default = \code{10}.
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
#' @param x an object of class sfapanel1 (returned by the function 
#' \code{\link{sfapanel1}}).
#' @param ... additional arguments of frontier are passed to sfapanel; 
#' additional arguments of the print, bread, estfun, nobs methods are currently
#'  ignored.
#'
#' @details
#' The stochastic frontier model for the panel data in the case of the time 
#' invariant inefficiency model (Pitt and Lee (1981): \code{'pl81'}) is defined 
#' as: 
#' 
#' \deqn{y_{it} = \alpha + \mathbf{x_{it}^{\prime}}\bm{\beta} +  v_{it} - Su_i}
#' 
#' with
#' 
#' \deqn{\epsilon_{it} = v_{it} -Su_i}
#'
#' where \eqn{i} is the cross section and \eqn{t} the period the cross section 
#' is observed, \eqn{y} is the output (cost, revenue, profit), \eqn{\mathbf{x}} 
#' is the vector of main explanatory variables (inputs and other control 
#' variables), \eqn{u} is the one-sided error term with variance 
#' \eqn{\sigma_{u}^2}, and \eqn{v} is the two-sided error term with variance 
#' \eqn{\sigma_{v}^2}.
#'
#' \code{S = 1} in the case of production (profit) frontier function and
#' \code{S = -1} in the case of cost frontier function.
#' 
#' All models are estimated using maximum likelihood (ML) for most distributions
#' except the Gamma, Weibull and log-normal distributions for which maximum
#' simulated likelihood (MSL) is used. For this latter, several draws can be
#' implemented e.g., Halton, Generalized Halton, Sobol, or uniform. In the
#' case of uniform draws, antithetics can also be computed: first \code{Nsim/2}
#' draws are obtained, then the \code{Nsim/2} other draws are obtained as
#' counterpart of one (\code{1-draw}).
#' 
#' In the case of the truncated normal distribution, the density of the 
#' convolution for each cross section is obtained as:
#' 
#' \deqn{
#' f(\bm{\epsilon}_i)=\frac{\sigma_*}{\sigma_u\sigma_v^T(2\pi)^{T/2}
#' \Phi\left(\frac{\mu}{\sigma_u}\right)}\exp{\left[-\frac{1}{2}\left(
#' -\frac{\mu_{i*}^2}{\sigma_*^2} + \frac{\sum_t\epsilon_{it}^2}{\sigma_v^2} + 
#' \frac{\mu^2}{\sigma_u^2}\right)\right]}\Phi\left(\frac{\mu_{i*}}{\sigma_*}
#' \right)
#' }
#' 
#' where
#' 
#' \deqn{\mu_{i*}=\frac{\mu\sigma_v^2-S\sigma_u^2\sum_t\epsilon_{it}}{
#' T\sigma_u^2+\sigma_v^2}}
#' 
#' \deqn{\sigma_*^2=\frac{\sigma_u^2\sigma_v^2}{T\sigma_u^2+\sigma_v^2}}
#' 
#' and \eqn{T} is the last period of observation of cross section \eqn{i}.
#' 
#' To account for heteroscedasticity in the variance parameters of the error
#' terms, a single part (right) formula can also be specified. To impose the
#' positivity to these parameters, the variances are modelled as:
#' \eqn{\sigma^2_{ui} = \exp{(\bm{\delta}'\mathbf{Z}_{ui})}} or 
#' \eqn{\sigma^2_{vi} = \exp{(\bm{\phi}'\mathbf{Z}_{vi})}}, where 
#' \eqn{\mathbf{Z}_{ui}} and \eqn{\mathbf{Z}_{vi}} are the heteroscedasticity
#' variables (inefficiency drivers in the case of \eqn{\mathbf{Z}_{ui}}) and 
#' \eqn{\bm{\delta}} and \eqn{\bm{\phi}} the coefficients. In the case of 
#' heterogeneity in the truncated mean \eqn{\mu_i}, it is modelled as 
#' \eqn{\mu_i=\bm{\omega}'\mathbf{Z}_{\mu_i}}. The \eqn{\mathbf{Z}} variables 
#' are such that they are constant over time for each cross section \eqn{i}. 
#' When this is not the case, one of the \code{'invariance'} option is 
#' activated.
#' 
#' In the case of the time varying inefficiency, the model writes as follows:
#' 
#' \deqn{y_{it} = \alpha + \mathbf{x_{it}^{\prime}}\bm{\beta} +
#'  v_{it} - SG(t)u_i}
#' 
#' where \eqn{G(t)} is defined following different specifications:
#' 
#' \deqn{\begin{array}{lll}
#' \text{References} & \text{Options in function} & \text{Formula} \\
#' \text{Kumbhakar (1990)} & \text{\code{'modelType = 'k90''}} & 
#' G(t)=\left[1+\exp{\left(\eta_1t + \eta_2t^2\right)}\right]^{-1}\\ 
#' \text{Battese and Coelli (1992)} & \text{\code{'modelType = 'bc92a''}} & 
#' G(t)=\exp{\left[-\eta(t-T)\right]} \\
#' \text{Cuesta and Orea (2002),  Feng and Serletis (2009)} & 
#' \text{\code{'modelType = 'bc92b''}} & 
#' G(t)=\exp{\left[-\eta_1(t-T)-\eta_2(t-T)^2\right]} \\
#' \text{Alvarez et al. (2006)} & \text{\code{'modelType = 'bc92c''}} & 
#' G(t)=\exp{\left(\bm{\eta}'\mathbf{Z}_{git}\right)} \\
#' \text{Battese and Coelli (1992)} & \text{\code{'modelType = 'mbc92''}} & 
#' G(t)= 1 + \eta_1(t-T) + \eta_2(t-T)^2\\
#' \text{Cuesta (2000)} & \text{\code{'modelType = 'c00''}} & 
#' G(t)=\exp{\left[-\eta_i(t-T)\right]} \\
#' \text{Kumbhakar and Wang (2005)} & \text{\code{'modelType = 'kw05''}} & 
#' G(t)=\exp{\left[-\eta(t-t_1)\right]} \\
#' \text{Modified Lee and Schmidt (1993)} & 
#' \text{\code{'modelType = 'mols93''}} & G(t)=\exp{\left[-\eta_t(t-T)\right]}
#' \end{array}
#' }
#' 
#' For the Modified Lee and Schmidt (1993) model, the last period parameter is not
#' identifiable so we set \eqn{\eta_T=1}.
#' 
#' In the case of truncated normal distribution, the density is similar to the 
#' previous equation except that:
#' 
#' \deqn{\mu_{i*}=\frac{\mu\sigma_v^2-S\sigma_u^2\sum_tG(t)\epsilon_{it}}{\sigma_u^2\sum_tG(t)^2+\sigma_v^2}}
#' 
#' \deqn{\sigma_*^2=\frac{\sigma_u^2\sigma_v^2}{\sigma_u^2\sum_tG(t)^2+\sigma_v^2}}
#' 
#' In the case, \code{'modelType = 'k90''}, or \code{'modelType = 'bc92a''},or
#' \code{'modelType = 'bc92b''}, or \code{'modelType = 'mbc92''}, or 
#' \code{'modelType = 'kw05''}, or \code{'modelType = 'mols93''}, specifying 
#' \code{'muhet'} or \code{'uhet'} options does not have any impact in the 
#' estimation because for these models, only the heteroscedasticity in the 
#' two-sided error term is considered.
#' 
#' \code{sfapanel1} allows for the maximization of weighted log-likelihood.
#' When option \code{weights} is specified and \code{wscale = TRUE}, the weights
#' are scaled as 
#' 
#' \deqn{new_{weights} = sample_{size} \times \frac{old_{weights}}{\sum(old_{weights})}}
#' 
#' For the weight variable, the \code{'invariance'} option is fist activated
#' before the scaling.
#' 
#' For complex problems, non-gradient methods (e.g. \code{nm} or \code{sann}) 
#' can be used to warm start the optimization and zoom in the neighborhood of 
#' the solution. Then a gradient-based methods is recommended in the second 
#' step. In the case of \code{sann}, we recommend to significantly increase the
#'  iteration limit (e.g. \code{itermax = 20000}). The Conjugate Gradient 
#'  (\code{cg}) can also be used in the first stage.
#' 
#' A set of extractor functions for fitted model objects is available for 
#' objects of class \code{'sfapanel1'} including methods to the generic functions 
#' \code{\link[=print.sfapanel1]{print}}
#'
#' @return \code{\link{sfapanel1}} returns a list of class \code{'sfapanel1'}
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
#' \item{Nobs}{Total number of observations.}
#' 
#' \item{Nid}{Number of cross-sections in the panel data.}
#' 
#' \item{Vtime}{Vector of total number of periods of presence of each 
#' cross-section.}
#' 
#' \item{Ntime}{Average periods of presence of each cross-section.}
#' 
#' \item{nXvar}{Number of explanatory variables in the production or cost
#' frontier.}
#'
#' \item{nmuZUvar}{Number of variables explaining heterogeneity in the
#' truncated mean, only if \code{udist = 'tnormal'} or \code{'lognormal'}.}
#' 
#' \item{nuZUvar}{Number of variables explaining heteroscedasticity in the
#' one-sided error term.}
#'
#' \item{nvZVvar}{Number of variables explaining heteroscedasticity in the
#' two-sided error term.}
#' 
#' \item{logDepVar}{The argument \code{'logDepVar'}. See the section
#' \sQuote{Arguments}.}
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
#' \item{modelType}{Spefication used for the panel model.}
#' 
#' \item{gHvar}{Matrix of variables used in the case of time varying 
#' inefficiency when \code{'modelType %in% c('bc92a', 'bc92b', 'bc92c', 'kw05',
#' 'c00', 'mols93')'} (for internal use).}
#' 
#' \item{ngZGvar}{Number of columns in gHvar when available.}
#' 
#' \item{invariance}{Methodology used to obtain time invariant exogenous 
#' variables (in the case of heteroscedasticity).}
#' 
#' \item{initHalf}{When \code{'whichStart = 2L'} and \code{'udist = 'hnormal''}. 
#' Initial ML estimation with half normal distribution for the one-sided error 
#' term. Model to construct the starting values. In the case of other 
#' distributions, we have: \code{'initExpo'} for \code{'udist = 'exponential''}; 
#' \code{'initTrunc'} for \code{'udist = 'tnormal''}; 
#' \code{'initRay'} for \code{'udist = 'rayleigh''}
#' \code{'initUni'} for \code{'udist = 'uniform''}
#' \code{'initGamma'} for \code{'udist = 'gamma''}
#' \code{'initLog'} for \code{'udist = 'lognormal''}
#' \code{'initWeibull'} for \code{'udist = 'weibull''}
#' \code{'initGenExpo'} for \code{'udist = 'genexponential''}
#' \code{'initTSL'} for \code{'udist = 'tslaplace''}. Object of class 'maxLik' 
#' and 'maxim' returned.}
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
#' \item{conditionNums}{Matrix. Condition number adding columns one by one.}
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
#'@seealso \code{\link[=print.sfapanel1]{print}} for printing \code{sfapanel1} 
#' object.
#' 
#' @export
#' 
#' @references Aigner, D., Lovell, C. A. K., and Schmidt, P. 1977. Formulation
#' and estimation of stochastic frontier production function models.
#' \emph{Journal of Econometrics}, \bold{6}(1), 21--37.
#' 
#' Alvarez, A., Amsler, C., Orea, L., & Schmidt, P. (2006). 
#' Interpreting and Testing the Scaling Property in Models where Inefficiency 
#' Depends on Firm Characteristics. \emph{Journal of Productivity Analysis}, 
#' \bold{25}(3), 201-212. doi: 10.1007/s11123-006-7639-3.
#' 
#' Battese, G. E., & Coelli, T. J. (1992). Frontier production 
#' functions, technical efficiency and panel data: With application to paddy 
#' farmers in India. \emph{Journal of Productivity Analysis}, \bold{3}(1-2), 
#' 153-169. doi: 10.1007/bf00158774.
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
#' Cuesta, R. A. (2000). A production model with firm-specific temporal 
#' variation in technical inefficiency: With application to Spanish dairy farms. 
#' \emph{Journal of Productivity Analysis}, \bold{13}(2), 139-158. 
#' doi: Doi 10.1023/A:1017297831646.
#' 
#' Cuesta, R. A., & Orea, L. (2002). Mergers and technical efficiency in Spanish 
#' savings banks: A stochastic distance function approach. 
#' \emph{Journal of Banking & Finance}, \bold{26}(12), 2231-2247. 
#' Doi 10.1016/S0378-4266(01)00184-4.
#' 
#' Feng, G., & Serletis, A. (2009). Efficiency and productivity of the US 
#' banking industry, 1998–2005: evidence from the Fourier cost function 
#' satisfying global regularity conditions.
#' \emph{Journal of Applied Econometrics}, \bold{24}(1), 105-138. 
#' doi: https://doi.org/10.1002/jae.1021.
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
#' Kumbhakar, S. C. (1990). Production Frontiers, Panel Data, and Time-Varying 
#' Technical Inefficiency. \emph{Journal of Econometrics}, \bold{46}(1-2), 
#' 201-211. doi: Doi 10.1016/0304-4076(90)90055-X.
#' 
#' Kumbhakar, S. C., Ghosh, S., and McGuckin, J. T. 1991) A generalized
#' production frontier approach for estimating determinants of inefficiency in
#' U.S. dairy farms. \emph{Journal of Business & Economic Statistics},
#' \bold{9}(3), 279--286.
#' 
#' Kumbhakar, S. C., & Wang, H.-J. (2005). Estimation of growth convergence 
#' using a stochastic production frontier approach. 
#' \emph{Economics Letters}, \bold{88}(3), 300-305. 
#' doi: https://doi.org/10.1016/j.econlet.2005.01.023.
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
#' Pitt, M. M., & Lee, L. F. (1981). The Measurement and Sources of Technical 
#' Inefficiency in the Indonesian Weaving Industry. 
#' \emph{Journal of Development Economics}, \bold{9}(1), 43-64. 
#' doi: Doi 10.1016/0304-3878(81)90004-3.
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
#' Wang, H.J., and Schmidt, P. 2002. One-step and two-step estimation of the
#' effects of exogenous variables on technical efficiency levels. \emph{Journal
#' of Productivity Analysis}, \bold{18}:129--144.
#' 
#' Wang, K., and Ye, X. 2020. Development of alternative stochastic frontier
#' models for estimating time-space prism vertices. \emph{Transportation}.
#'
#' @examples
#' ## Swiss railways data
#' # Pitt and Lee (1981) time invariant inefficiency
#' data <- swissrailways
#' # remove observations only present one year
#' data <- data[data$ID != names(which(table(data$ID) == 1)), ]
#' data <- plm::pdata.frame(data, index = c('ID', 'YEAR'))
#' res_pl81 <- sfapanel1(formula = LNCT ~ LNQ2 + LNQ3 + LNNET + LNPK + LNPL, 
#' modelType = 'pl81', S = -1, data = data)
sfapanel1 <- function(formula, muhet, uhet, vhet, logDepVar = TRUE, data, idVar = NULL,
  timeVar = NULL, subset, weights, wscale = TRUE, S = 1L, modelType = "bc92a",
  udist = "hnormal", start = NULL, randStart = FALSE, whichStart = 2L, initAlg = "nm",
  initIter = 500, invariance = 2L, method = "bfgs", hessianType = 1L, simType = "halton",
  Nsim = 100, prime = 2L, burn = 10, antithetics = FALSE, seed = 12345, itermax = 2000,
  printInfo = FALSE, tol = 1e-12, gradtol = 1e-06, stepmax = 0.1, qac = "marquardt") {
  # panel model check -------
  modelType <- tolower(modelType)
  if (!(modelType %in% c("pl81", "bc92a", "bc92b", "bc92c", "mbc92", "k90", "kw05",
    "c00", "mols93"))) {
    stop("Unknown SFA panel model: ", paste(modelType), call. = FALSE)
  }
  # u distribution check -------
  udist <- tolower(udist)
  if (!(udist %in% c("hnormal", "exponential", "tnormal", "rayleigh", "uniform",
    "gamma", "lognormal", "weibull", "genexponential", "tslaplace")))
    stop("Unknown inefficiency distribution: ", paste(udist), call. = FALSE)
  # check data for panel dimension -------
  if (missing(data))
    data <- environment(formula)
  if (!inherits(data, "pdata.frame")) {
    if (is.null(idVar) & is.null(timeVar)) {
      stop("'data' must be of class 'pdata.frame' or arguments 'idVar' & 'timeVar' must be provided",
        call. = FALSE)
    } else {
      if (is.null(idVar) & !is.null(timeVar)) {
        stop("Argument 'idVar' must be provided", call. = FALSE)
      } else {
        if (!is.null(idVar) & is.null(timeVar)) {
          stop("Arguments 'timeVar' must be provided", call. = FALSE)
        } else {
          if (!is.null(idVar) & !is.null(timeVar)) {
          data$origin_names <- rownames(data)
          data <- pdata.frame(data, index = c(idVar, timeVar))
          }
        }
      }
    }
  }
  # Formula manipulation -------
  if (length(Formula(formula))[2] != 1) {
    stop("argument 'formula' must have one RHS part", call. = FALSE)
  }
  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"), names(mc), nomatch = 0L)
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
  if (modelType == "bc92c") {
    if (udist %in% c("tnormal", "lognormal")) {
      if (!identical(muhet, uhet)) {
        stop("argument 'muhet' and 'uhet' must contain the same variables for modelType = 'bc92c'",
          call. = FALSE)
      } else {
        if (length(attr(terms(muhet), "term.labels")) == 0 || length(attr(terms(uhet),
          "term.labels")) == 0) {
          stop("at least one exogeneous variable must be provided for the scaling option for modelType = 'bc92c'",
          call. = FALSE)
        } else {
          muhet <- ~1
          ghet <- plhsCheck_u_bc92c(formula = uhet)
          uhet <- ~1
        }
      }
    } else {
      if (length(attr(terms(uhet), "term.labels")) == 0) {
        stop("at least one exogeneous variable must be provided for the scaling option for modelType = 'bc92c'",
          call. = FALSE)
      } else {
        ghet <- plhsCheck_u_bc92c(formula = uhet)
        uhet <- ~1
      }
    }
    if (!missing(vhet)) {
      vhet <- clhsCheck_v(formula = vhet)
    } else {
      vhet <- ~1
    }
  }
  if (modelType == "bc92c") {
    formula <- formDist_sfapanel1_bc92c(udist = udist, formula = formula, muhet = muhet,
      uhet = uhet, vhet = vhet, ghet = ghet)
  } else {
    formula <- formDist_sfacross(udist = udist, formula = formula, muhet = muhet,
      uhet = uhet, vhet = vhet)
  }
  # Generate required datasets -------
  mc$formula <- formula
  mc$na.action <- na.omit
  mc[[1L]] <- quote(model.frame)
  mc <- eval(mc, parent.frame())
  validObs <- rowSums(is.na(mc) | is.infinite.data.frame(mc)) == 0
  Yvar <- model.response(mc, "numeric")
  Yvar <- Yvar[validObs]
  mtX <- terms(formula, data = data, rhs = 1)
  Xvar <- model.matrix(mtX, mc)
  Xvar <- Xvar[validObs, , drop = FALSE]
  nXvar <- ncol(Xvar)
  # if subset is non-missing and there NA, force data to change
  pdimNames <- names(index(data))
  if (!is.null(idVar) & !is.null(timeVar)) {
    data <- data[data$origin_names %in% attr(mc, "row.names"), ]
  } else {
    data <- data[row.names(data) %in% attr(mc, "row.names"), ]
  }
  data <- data[validObs, ]
  # since data is pdata.frame I can drop unused levels
  data[, pdimNames[1]] <- droplevels(data[, pdimNames[1]])
  data[, pdimNames[2]] <- droplevels(data[, pdimNames[2]])
  pindex <- index(data)
  NT <- nrow(Xvar)
  TT <- as.numeric(table(pindex[, 1]))
  N <- length(TT)
  if (NT == 0L) {
    stop("0 (non-NA) cases", call. = FALSE)
  }
  if (!(invariance %in% 1:3))
    stop("'invariance' can either be 1, 2 or 3", call. = FALSE)
  wHvar <- as.vector(model.weights(mc))
  if (length(wscale) != 1 || !is.logical(wscale[1])) {
    stop("argument 'wscale' must be a single logical value", call. = FALSE)
  }
  if (!is.null(wHvar)) {
    if (!is.numeric(wHvar)) {
      stop("'weights' must be a numeric vector", call. = FALSE)
    } else {
      if (any(wHvar < 0 | is.na(wHvar)))
        stop("missing or negative weights not allowed", call. = FALSE)
    }
    if (wscale) {
      wHvar_c <- wHvar/sum(wHvar) * NT
      varw_p <- if (invariance == 1) {
        as.vector(tapply(wHvar, pindex[, 1], function(u) u[1]))
      } else {
        if (invariance == 2) {
          as.vector(tapply(wHvar, pindex[, 1], function(u) u[length(u)]))
        } else {
          if (invariance == 3) {
          as.vector(tapply(wHvar, pindex[, 1], mean))
          }
        }
      }
      wHvar_p <- N * varw_p/sum(varw_p)
    }
  } else {
    wHvar_c <- rep(1, NT)
    wHvar_p <- rep(1, N)
  }
  if (length(Yvar) != nrow(Xvar)) {
    stop(paste("the number of observations of the dependent variable (", length(Yvar),
      ") must be the same to the number of observations of the exogenous variables (",
      nrow(Xvar), ")", sep = ""), call. = FALSE)
  }
  # check for hetero in U (first, last, mean) -------
  if (udist %in% c("tnormal", "lognormal")) {
    mtmuH <- delete.response(terms(formula, data = data, rhs = 2))
    muHvar_c <- model.matrix(mtmuH, mc)
    muHvar_c <- muHvar_c[validObs, , drop = FALSE]
    nmuZUvar <- ncol(muHvar_c)
    mtuH <- delete.response(terms(formula, data = data, rhs = 3))
    uHvar_c <- model.matrix(mtuH, mc)
    uHvar_c <- uHvar_c[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar_c)
    mtvH <- delete.response(terms(formula, data = data, rhs = 4))
    vHvar_c <- model.matrix(mtvH, mc)
    vHvar_c <- vHvar_c[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar_c)
    if (modelType == "bc92c") {
      mtgH <- delete.response(terms(formula, data = data, rhs = 5))
    }
    if (invariance == 1) {
      muHvar_p <- apply(muHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
    } else {
      if (invariance == 2) {
        muHvar_p <- apply(muHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
      } else {
        if (invariance == 3) {
          muHvar_p <- apply(muHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
          uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
        }
      }
    }
  } else {
    mtuH <- delete.response(terms(formula, data = data, rhs = 2))
    uHvar_c <- model.matrix(mtuH, mc)
    uHvar_c <- uHvar_c[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar_c)
    mtvH <- delete.response(terms(formula, data = data, rhs = 3))
    vHvar_c <- model.matrix(mtvH, mc)
    vHvar_c <- vHvar_c[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar_c)
    if (modelType == "bc92c") {
      mtgH <- delete.response(terms(formula, data = data, rhs = 4))
    }
    if (invariance == 1) {
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
    } else {
      if (invariance == 2) {
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
      } else {
        if (invariance == 3) {
          uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) mean(u))
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) mean(u))
          })
        }
      }
    }
  }
  if (modelType == "bc92a") {
    gHvar <- cbind(eta = unlist(lapply(TT, FUN = function(x) -rev(-(0:(x - 1))))))
    ngZGvar <- dim(gHvar)[2]
  } else {
    if (modelType == "bc92b") {
      gHvar <- cbind(eta1 = unlist(lapply(TT, FUN = function(x) -rev(-(0:(x -
        1))))), eta2 = unlist(lapply(TT, FUN = function(x) -rev(-(0:(x -
        1)))^2)))
      ngZGvar <- dim(gHvar)[2]
    } else {
      if (modelType == "bc92c") {
        gHvar <- model.matrix(mtgH, mc)
        gHvar <- gHvar[validObs, , drop = FALSE]
        ngZGvar <- dim(gHvar)[2]
      } else {
        if (modelType == "kw05") {
          gHvar <- cbind(eta = unlist(lapply(TT, FUN = function(x) (seq(1,
          x) - 1))))
          ngZGvar <- dim(gHvar)[2]
        } else {
          if (modelType == "c00") {
          iHvar <- model.matrix(~-1 + as.factor(data[, names(pindex)[1]]),
            rhs = 1)
          nameID <- levels(data[, names(pindex)[1]])
          colnames(iHvar) <- nameID
          gHvar <- sweep(iHvar, MARGIN = 1, STATS = unlist(lapply(TT, FUN = function(x) -rev(-(0:(x -
            1))))), FUN = "*")
          ngZGvar <- dim(gHvar)[2]
          } else {
          if (modelType == "mols93") {
            tHvar <- model.matrix(~-1 + as.factor(data[, names(pindex)[2]]),
            rhs = 1)
            nameIT <- levels(data[, names(pindex)[2]])
            colnames(tHvar) <- nameIT
            gHvar <- sweep(tHvar, MARGIN = 1, STATS = unlist(lapply(TT,
            FUN = function(x) -rev(-(0:(x - 1))))), FUN = "*")
            ngZGvar <- dim(gHvar)[2]
          }
          }
        }
      }
    }
  }
  # check other supplied options -------
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
    stop("argument 'logDepVar' must be a single logical value", call. = FALSE)
  }
  # Number of parameters -------
  if (modelType == "pl81") {
    nParm <- if (udist %in% c("tnormal", "lognormal")) {
      nXvar + nmuZUvar + nuZUvar + nvZVvar
    } else {
      if (udist %in% c("gamma", "weibull", "tslaplace")) {
        nXvar + nuZUvar + nvZVvar + 1
      } else {
        nXvar + nuZUvar + nvZVvar
      }
    }
  } else {
    if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00")) {
      nParm <- if (udist %in% c("tnormal", "lognormal")) {
        nXvar + nmuZUvar + nuZUvar + ngZGvar + nvZVvar
      } else {
        if (udist %in% c("gamma", "weibull", "tslaplace")) {
          nXvar + ngZGvar + nuZUvar + nvZVvar + 1
        } else {
          nXvar + ngZGvar + nuZUvar + nvZVvar
        }
      }
    } else {
      if (modelType %in% c("k90", "mbc92")) {
        nParm <- if (udist %in% c("tnormal", "lognormal")) {
          nXvar + nmuZUvar + nuZUvar + nvZVvar + 2
        } else {
          if (udist %in% c("gamma", "weibull", "tslaplace")) {
          nXvar + nuZUvar + nvZVvar + 1 + 2
          } else {
          nXvar + nuZUvar + nvZVvar + 2
          }
        }
      } else {
        if (modelType == "mols93") {
          nParm <- if (udist %in% c("tnormal", "lognormal")) {
          nXvar + nmuZUvar + nuZUvar + (ngZGvar - 1) + nvZVvar
          } else {
          if (udist %in% c("gamma", "weibull", "tslaplace")) {
            nXvar + nuZUvar + (ngZGvar - 1) + nvZVvar + 1
          } else {
            nXvar + nuZUvar + (ngZGvar - 1) + nvZVvar
          }
          }
        }
      }
    }
  }
  # checking starting values when provided -------
  if (!is.null(start)) {
    if (length(start) != nParm) {
      stop("Wrong number of initial values: model has ", nParm, " parameters",
        call. = FALSE)
    }
  }
  if (nParm > NT) {
    stop("Model has more parameters than observations", call. = FALSE)
  }
  # Set std. Error when random start is allowed -------
  sdStart <- if (randStart)
    0.01 else NULL
  # check whichStart -------
  if (length(whichStart) != 1 || !(whichStart %in% c(1L, 2L))) {
    stop("argument 'whichStart' must equal either 1 or 2", call. = FALSE)
  }
  # check algorithms -------
  method <- tolower(method)
  if (!(method %in% c("ucminf", "bfgs", "bhhh", "nr", "nm", "cg", "sann", "sr1",
    "mla", "sparse", "nlminb"))) {
    stop("Unknown or non-available optimization algorithm: ", paste(method),
      call. = FALSE)
  }
  # check hessian type
  if (length(hessianType) != 1 || !(hessianType %in% c(1L, 2L))) {
    stop("argument 'hessianType' must equal either 1 or 2", call. = FALSE)
  }
  # Draws for SML -------
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    if (!(simType %in% c("halton", "ghalton", "rhalton", "sobol", "rsobol", "richtmyer",
      "rrichtmyer", "uniform", "mlhs"))) {
      stop("Unknown or non-available random draws method", call. = FALSE)
    }
    if (!is.numeric(Nsim) || length(Nsim) != 1) {
      stop("argument 'Nsim' must be a single numeric scalar", call. = FALSE)
    }
    if (!is.numeric(burn) || length(burn) != 1) {
      stop("argument 'burn' must be a single numeric scalar", call. = FALSE)
    }
    if (!is_prime(prime)) {
      stop("argument 'prime' must be a single prime number", call. = FALSE)
    }
    if (length(antithetics) != 1 || !is.logical(antithetics[1])) {
      stop("argument 'antithetics' must be a single logical value", call. = FALSE)
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
          if (simType == "rsobol") {
          "Randomized Sobol"
          } else {
          if (simType == "richtmyer") {
            "Richtmyer"
          } else {
            if (simType == "rrichtmyer") {
            "Randomized Richtmyer"
            } else {
            if (simType == "uniform") {
              "Uniform"
            } else {
              if (simType == "mlhs") {
              "Modified Latin Hypercube Sampling"
              }
            }
            }
          }
          }
        }
      }
    }
    cat("Initialization of", Nsim, simDist, "draws per observation/cross-section ...\n")
    FiMat_N <- drawMatUniDim(N = N, Nsim = Nsim, simType = simType, prime = prime,
      burn = burn + 1, antithetics = antithetics, seed = seed)
    FiMat_NT <- drawMatUniDim(N = NT, Nsim = Nsim, simType = simType, prime = prime,
      burn = burn + 1, antithetics = antithetics, seed = seed)
  }
  # Other optimization options -------
  if (!is.numeric(initIter) || length(initIter) != 1) {
    stop("argument 'initIter' must be a single numeric scalar", call. = FALSE)
  }
  if (initIter != round(initIter)) {
    stop("argument 'initIter' must be an integer", call. = FALSE)
  }
  if (initIter <= 0) {
    stop("argument 'initIter' must be positive", call. = FALSE)
  }
  initIter <- as.integer(initIter)
  if (!is.numeric(itermax) || length(itermax) != 1) {
    stop("argument 'itermax' must be a single numeric scalar", call. = FALSE)
  }
  initAlg <- tolower(initAlg)
  if (!(initAlg %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann"))) {
    stop("Unknown or non-available optimization algorithm: ", paste(initAlg),
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
    stop("argument 'printInfo' must be a single logical value", call. = FALSE)
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
    stop("argument 'qac' must be either 'marquardt' or 'stephalving'", call. = FALSE)
  }
  # Step 1: OLS -------
  olsRes <- if (colnames(Xvar)[1] == "(Intercept)") {
    # shall I initialize with lm or plm with RE???
    if (dim(Xvar)[2] == 1) {
      lm(Yvar ~ 1)
    } else {
      lm(Yvar ~ ., data = as.data.frame(Xvar[, -1]), weights = wHvar_c)
    }
  } else {
    lm(Yvar ~ -1 + ., data = as.data.frame(Xvar), weights = wHvar_c)
  }
  if (any(is.na(olsRes$coefficients))) {
    stop("at least one of the OLS coefficients is NA: ", paste(colnames(Xvar)[is.na(olsRes$coefficients)],
      collapse = ", "), "This may be due to a singular matrix
   due to potential perfect multicollinearity",
      call. = FALSE)
  }
  names(olsRes$coefficients) <- colnames(Xvar)
  olsParam <- c(olsRes$coefficients)
  olsSigmasq <- summary(olsRes)$sigma^2
  olsStder <- sqrt(diag(vcov(olsRes)))
  olsLoglik <- logLik(olsRes)[1]
  dataTable <- data[, names(pindex)]
  dataTable <- cbind(dataTable, data[, all.vars(terms(formula))], weights = wHvar_c)
  dataTable <- cbind(dataTable, olsResiduals = residuals(olsRes), olsFitted = fitted(olsRes))
  # possibility to have duplicated columns if ID or TIME appears in ols
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
  CoelliM3Test <- c(z = m3/sqrt(6 * m2^3/NT), p.value = 2 * pnorm(-abs(m3/sqrt(6 *
    m2^3/NT))))
  AgostinoTest <- dagoTest(dataTable[["olsResiduals"]])
  class(AgostinoTest) <- "dagoTest"
  # Step 2: MLE arguments -------
  if (modelType %in% c("pl81", "k90", "mbc92")) {
    FunArgs <- if (udist == "tnormal") {
      list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam,
        dataTable = dataTable, nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, muHvar_c = muHvar_c, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
        muHvar_p = muHvar_p, uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, S = S, wHvar_c = wHvar_c, wHvar_p = wHvar_p, pindex = pindex,
        TT = TT, method = method, printInfo = printInfo, itermax = itermax,
        stepmax = stepmax, tol = tol, gradtol = gradtol, hessianType = hessianType,
        initAlg = initAlg, initIter = initIter, whichStart = whichStart,
        qac = qac)
    } else {
      if (udist == "lognormal") {
        list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam,
          dataTable = dataTable, nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, muHvar_c = muHvar_c, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
          muHvar_p = muHvar_p, uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
          Xvar = Xvar, wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S, pindex = pindex,
          TT = TT, N = N, NT = NT, FiMat_N = FiMat_N, FiMat_NT = FiMat_NT,
          method = method, printInfo = printInfo, itermax = itermax, stepmax = stepmax,
          tol = tol, gradtol = gradtol, hessianType = hessianType, initAlg = initAlg,
          initIter = initIter, whichStart = whichStart, qac = qac)
      } else {
        if (udist %in% c("gamma", "weibull")) {
          list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam,
          dataTable = dataTable, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar_c = uHvar_c, vHvar_c = vHvar_c, uHvar_p = uHvar_p, vHvar_p = vHvar_p,
          Yvar = Yvar, Xvar = Xvar, wHvar_c = wHvar_c, wHvar_p = wHvar_p,
          S = S, pindex = pindex, TT = TT, N = N, NT = NT, FiMat_N = FiMat_N,
          FiMat_NT = FiMat_NT, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, initAlg = initAlg, initIter = initIter,
          whichStart = whichStart, qac = qac)
        } else {
          list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam,
          dataTable = dataTable, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar_c = uHvar_c, vHvar_c = vHvar_c, uHvar_p = uHvar_p, vHvar_p = vHvar_p,
          Yvar = Yvar, Xvar = Xvar, wHvar_c = wHvar_c, wHvar_p = wHvar_p,
          S = S, pindex = pindex, TT = TT, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, initAlg = initAlg, initIter = initIter,
          whichStart = whichStart, qac = qac)
        }
      }
    }
  } else {
    if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00")) {
      FunArgs <- if (udist == "tnormal") {
        list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam,
          dataTable = dataTable, nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, muHvar_c = muHvar_c, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
          muHvar_p = muHvar_p, uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
          Xvar = Xvar, modelType = modelType, ngZGvar = ngZGvar, gHvar = gHvar,
          S = S, wHvar_c = wHvar_c, wHvar_p = wHvar_p, pindex = pindex, TT = TT,
          method = method, printInfo = printInfo, itermax = itermax, stepmax = stepmax,
          tol = tol, gradtol = gradtol, hessianType = hessianType, initAlg = initAlg,
          initIter = initIter, whichStart = whichStart, qac = qac)
      } else {
        if (udist == "lognormal") {
          list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam,
          dataTable = dataTable, nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, muHvar_c = muHvar_c, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
          muHvar_p = muHvar_p, uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
          Xvar = Xvar, wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S, modelType = modelType,
          ngZGvar = ngZGvar, gHvar = gHvar, pindex = pindex, TT = TT, N = N,
          NT = NT, FiMat_N = FiMat_N, FiMat_NT = FiMat_NT, method = method,
          printInfo = printInfo, itermax = itermax, stepmax = stepmax,
          tol = tol, gradtol = gradtol, hessianType = hessianType, initAlg = initAlg,
          initIter = initIter, whichStart = whichStart, qac = qac)
        } else {
          if (udist %in% c("gamma", "weibull")) {
          list(start = start, randStart = randStart, sdStart = sdStart,
            olsParam = olsParam, dataTable = dataTable, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
            uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar,
            wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S, modelType = modelType,
            ngZGvar = ngZGvar, gHvar = gHvar, pindex = pindex, TT = TT,
            N = N, NT = NT, FiMat_N = FiMat_N, FiMat_NT = FiMat_NT, method = method,
            printInfo = printInfo, itermax = itermax, stepmax = stepmax,
            tol = tol, gradtol = gradtol, hessianType = hessianType, initAlg = initAlg,
            initIter = initIter, whichStart = whichStart, qac = qac)
          } else {
          list(start = start, randStart = randStart, sdStart = sdStart,
            olsParam = olsParam, dataTable = dataTable, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
            uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar,
            wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S, modelType = modelType,
            ngZGvar = ngZGvar, gHvar = gHvar, pindex = pindex, TT = TT,
            method = method, printInfo = printInfo, itermax = itermax,
            stepmax = stepmax, tol = tol, gradtol = gradtol, hessianType = hessianType,
            initAlg = initAlg, initIter = initIter, whichStart = whichStart,
            qac = qac)
          }
        }
      }
    } else {
      if (modelType == "mols93") {
        FunArgs <- if (udist == "tnormal") {
          list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam,
          dataTable = dataTable, nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, muHvar_c = muHvar_c, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
          muHvar_p = muHvar_p, uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
          Xvar = Xvar, ngZGvar = ngZGvar, gHvar = gHvar, S = S, wHvar_c = wHvar_c,
          wHvar_p = wHvar_p, pindex = pindex, TT = TT, method = method,
          printInfo = printInfo, itermax = itermax, stepmax = stepmax,
          tol = tol, gradtol = gradtol, hessianType = hessianType, initAlg = initAlg,
          initIter = initIter, whichStart = whichStart, qac = qac)
        } else {
          if (udist == "lognormal") {
          list(start = start, randStart = randStart, sdStart = sdStart,
            olsParam = olsParam, dataTable = dataTable, nXvar = nXvar,
            nmuZUvar = nmuZUvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
            muHvar_c = muHvar_c, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
            muHvar_p = muHvar_p, uHvar_p = uHvar_p, vHvar_p = vHvar_p,
            Yvar = Yvar, Xvar = Xvar, wHvar_c = wHvar_c, wHvar_p = wHvar_p,
            S = S, ngZGvar = ngZGvar, gHvar = gHvar, pindex = pindex, TT = TT,
            N = N, NT = NT, FiMat_N = FiMat_N, FiMat_NT = FiMat_NT, method = method,
            printInfo = printInfo, itermax = itermax, stepmax = stepmax,
            tol = tol, gradtol = gradtol, hessianType = hessianType, initAlg = initAlg,
            initIter = initIter, whichStart = whichStart, qac = qac)
          } else {
          if (udist %in% c("gamma", "weibull")) {
            list(start = start, randStart = randStart, sdStart = sdStart,
            olsParam = olsParam, dataTable = dataTable, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar_c = uHvar_c,
            vHvar_c = vHvar_c, uHvar_p = uHvar_p, vHvar_p = vHvar_p,
            Yvar = Yvar, Xvar = Xvar, wHvar_c = wHvar_c, wHvar_p = wHvar_p,
            S = S, ngZGvar = ngZGvar, gHvar = gHvar, pindex = pindex,
            TT = TT, N = N, NT = NT, FiMat_N = FiMat_N, FiMat_NT = FiMat_NT,
            method = method, printInfo = printInfo, itermax = itermax,
            stepmax = stepmax, tol = tol, gradtol = gradtol, hessianType = hessianType,
            initAlg = initAlg, initIter = initIter, whichStart = whichStart,
            qac = qac)
          } else {
            list(start = start, randStart = randStart, sdStart = sdStart,
            olsParam = olsParam, dataTable = dataTable, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar_c = uHvar_c,
            vHvar_c = vHvar_c, uHvar_p = uHvar_p, vHvar_p = vHvar_p,
            Yvar = Yvar, Xvar = Xvar, wHvar_c = wHvar_c, wHvar_p = wHvar_p,
            S = S, ngZGvar = ngZGvar, gHvar = gHvar, pindex = pindex,
            TT = TT, method = method, printInfo = printInfo, itermax = itermax,
            stepmax = stepmax, tol = tol, gradtol = gradtol, hessianType = hessianType,
            initAlg = initAlg, initIter = initIter, whichStart = whichStart,
            qac = qac)
          }
          }
        }
      }
    }
  }
  ## MLE run -------
  if (modelType == "pl81") {
    mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt_pl81,
      FunArgs), exponential = do.call(exponormAlgOpt_pl81, FunArgs), tnormal = do.call(truncnormAlgOpt_pl81,
      FunArgs), rayleigh = do.call(raynormAlgOpt_pl81, FunArgs), gamma = do.call(gammanormAlgOpt_pl81,
      FunArgs), uniform = do.call(uninormAlgOpt_pl81, FunArgs), lognormal = do.call(lognormAlgOpt_pl81,
      FunArgs), weibull = do.call(weibullnormAlgOpt_pl81, FunArgs), genexponential = do.call(genexponormAlgOpt_pl81,
      FunArgs), tslaplace = do.call(tslnormAlgOpt_pl81, FunArgs)), error = function(e) print(e))
  } else {
    if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00")) {
      mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt_gzit,
        FunArgs), exponential = do.call(exponormAlgOpt_gzit, FunArgs), tnormal = do.call(truncnormAlgOpt_gzit,
        FunArgs), rayleigh = do.call(raynormAlgOpt_gzit, FunArgs), gamma = do.call(gammanormAlgOpt_gzit,
        FunArgs), uniform = do.call(uninormAlgOpt_gzit, FunArgs), lognormal = do.call(lognormAlgOpt_gzit,
        FunArgs), weibull = do.call(weibullnormAlgOpt_gzit, FunArgs), genexponential = do.call(genexponormAlgOpt_gzit,
        FunArgs), tslaplace = do.call(tslnormAlgOpt_gzit, FunArgs)), error = function(e) print(e))
    } else {
      if (modelType == "k90") {
        mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt_k90,
          FunArgs), exponential = do.call(exponormAlgOpt_k90, FunArgs), tnormal = do.call(truncnormAlgOpt_k90,
          FunArgs), rayleigh = do.call(raynormAlgOpt_k90, FunArgs), gamma = do.call(gammanormAlgOpt_k90,
          FunArgs), uniform = do.call(uninormAlgOpt_k90, FunArgs), lognormal = do.call(lognormAlgOpt_k90,
          FunArgs), weibull = do.call(weibullnormAlgOpt_k90, FunArgs), genexponential = do.call(genexponormAlgOpt_k90,
          FunArgs), tslaplace = do.call(tslnormAlgOpt_k90, FunArgs)), error = function(e) print(e))
      } else {
        if (modelType == "mbc92") {
          mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt_mbc92,
          FunArgs), exponential = do.call(exponormAlgOpt_mbc92, FunArgs),
          tnormal = do.call(truncnormAlgOpt_mbc92, FunArgs), rayleigh = do.call(raynormAlgOpt_mbc92,
            FunArgs), gamma = do.call(gammanormAlgOpt_mbc92, FunArgs),
          uniform = do.call(uninormAlgOpt_mbc92, FunArgs), lognormal = do.call(lognormAlgOpt_mbc92,
            FunArgs), weibull = do.call(weibullnormAlgOpt_mbc92, FunArgs),
          genexponential = do.call(genexponormAlgOpt_mbc92, FunArgs), tslaplace = do.call(tslnormAlgOpt_mbc92,
            FunArgs)), error = function(e) print(e))
        } else {
          if (modelType == "mols93") {
          mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt_mols93,
            FunArgs), exponential = do.call(exponormAlgOpt_mols93, FunArgs),
            tnormal = do.call(truncnormAlgOpt_mols93, FunArgs), rayleigh = do.call(raynormAlgOpt_mols93,
            FunArgs), gamma = do.call(gammanormAlgOpt_mols93, FunArgs),
            uniform = do.call(uninormAlgOpt_mols93, FunArgs), lognormal = do.call(lognormAlgOpt_mols93,
            FunArgs), weibull = do.call(weibullnormAlgOpt_mols93, FunArgs),
            genexponential = do.call(genexponormAlgOpt_mols93, FunArgs),
            tslaplace = do.call(tslnormAlgOpt_mols93, FunArgs)), error = function(e) print(e))
          }
        }
      }
    }
  }
  if (inherits(mleList, "error")) {
    stop("The current error occurs during optimization:\n", mleList$message,
      call. = FALSE)
  }
  # Inverse Hessian + other -------
  mleList$invHessian <- vcovObj(mleObj = mleList$mleObj, hessianType = hessianType,
    method = method, nParm = nParm)
  mleList <- c(mleList, if (method == "ucminf") {
    list(type = "ucminf maximization", nIter = unname(mleList$mleObj$info["neval"]),
      status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$value, gradient = mleList$mleObj$gradient)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
      list(type = substr(mleList$mleObj$type, 1, 27), nIter = mleList$mleObj$iterations,
        status = mleList$mleObj$message, mleLoglik = mleList$mleObj$maximum,
        gradient = mleList$mleObj$gradient)
    } else {
      if (method == "sr1") {
        list(type = "SR1 maximization", nIter = mleList$mleObj$iterations,
          status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
          gradient = -mleList$mleObj$gradient)
      } else {
        if (method == "mla") {
          list(type = "Levenberg-Marquardt maximization", nIter = mleList$mleObj$ni,
          status = switch(mleList$mleObj$istop, `1` = "convergence criteria were satisfied",
            `2` = "maximum number of iterations was reached", `4` = "algorithm encountered a problem in the function computation"),
          mleLoglik = -mleList$mleObj$fn.value, gradient = -mleList$mleObj$grad)
        } else {
          if (method == "sparse") {
          list(type = "Sparse Hessian maximization", nIter = mleList$mleObj$iterations,
            status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
            gradient = -mleList$mleObj$gradient)
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
  # quick renaming (only when start is provided)-------
  if (!is.null(start)) {
    if (udist %in% c("tnormal", "lognormal")) {
      names(mleList$startVal) <- fName_mu_sfapanel1(Xvar = Xvar, udist = udist,
        muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p, modelType = modelType,
        gHvar = if (modelType %in% c("c00", "bc92c", "mols93"))
          gHvar else NULL)
    } else {
      names(mleList$startVal) <- fName_uv_sfapanel1(Xvar = Xvar, udist = udist,
        uHvar = uHvar_p, vHvar = vHvar_p, modelType = modelType, gHvar = if (modelType %in%
          c("c00", "bc92c", "mols93"))
          gHvar else NULL)
    }
    names(mleList$mlParam) <- names(mleList$startVal)
  }
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
  names(mleList$gradient) <- names(mleList$mlParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
  # Return object -------
  mlDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
  dataTable$mlResiduals <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$mlFitted <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  datlogL <- data.frame(levels(pindex[, 1]), logL_OBS = mleList$mleObj$logL_OBS/TT)
  names(datlogL)[1] <- names(pindex)[1]
  dataTable <- merge(dataTable, datlogL, by = names(pindex)[1])
  dataTable <- pdata.frame(dataTable, names(dataTable)[1:2])
  returnObj <- list()
  returnObj$call <- cl
  returnObj$formula <- formula
  returnObj$S <- S
  returnObj$typeSfa <- typeSfa
  returnObj$Nobs <- NT
  returnObj$Nid <- N
  returnObj$Vtime <- TT
  returnObj$Ntime <- mean(TT)
  returnObj$nXvar <- nXvar
  if (udist %in% c("tnormal", "lognormal")) {
    returnObj$nmuZUvar <- nmuZUvar
  }
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  returnObj$logDepVar <- logDepVar
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
  returnObj$modelType <- modelType
  if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00", "mols93")) {
    returnObj$gHvar <- gHvar
    returnObj$ngZGvar <- ncol(gHvar)
  }
  returnObj$invariance <- invariance
  returnObj$isWeights <- !all.equal(wHvar_p, rep(1, N))
  returnObj$optType <- mleList$type
  returnObj$nIter <- mleList$nIter
  returnObj$optStatus <- mleList$status
  returnObj$startLoglik <- mleList$startLoglik
  returnObj$mlLoglik <- mleList$mleLoglik
  returnObj$mlParam <- mleList$mlParam
  returnObj$gradient <- mleList$gradient
  datGradL <- as.data.frame(mleList$mleObj$gradL_OBS)
  datGradL <- cbind(levels(pindex[, 1]), datGradL)
  colnames(datGradL)[1] <- names(pindex)[1]
  returnObj$gradL_OBS <- datGradL
  returnObj$gradientNorm <- sqrt(sum(mleList$gradient^2))
  returnObj$invHessian <- mleList$invHessian
  returnObj$conditionNums <- condiNum(mleObj = mleList$mleObj, method = method,
    nParm = nParm)
  returnObj$hessianType <- if (hessianType == 1) {
    "Analytic/Numeric Hessian"
  } else {
    if (hessianType == 2) {
      "BHHH Hessian"
    }
  }
  returnObj$mlDate <- mlDate
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    returnObj$simDist <- simDist
    returnObj$Nsim <- Nsim
    returnObj$FiMat <- FiMat_N
  }
  rm(mleList)
  class(returnObj) <- "sfapanel1"
  return(returnObj)
}

# print for sfapanel1 ----------
#' @rdname sfapanel1
#' @export
print.sfapanel1 <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat(sfapaneldist(x$udist), "\n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}

# Bread for Sandwich Estimator ----------
#' @rdname sfapanel1
#' @export
bread.sfapanel1 <- function(x, ...) {
  if (x$hessianType == "Analytic/Numeric Hessian") {
    return(x$invHessian * x$Nid)
  } else {
    cat("Computing Analytical/Numerical Hessian \n")
    Yvar <- model.response(model.frame(x$formula, data = x$dataTable))
    Xvar <- model.matrix(x$formula, rhs = 1, data = x$dataTable)
    if (x$modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00", "mols93")) {
      gHvar <- x$gHvar
      ngZGvar <- dim(gHvar)[2]
    }
    pindex <- x$dataTable[, 1:2]
    invariance <- x$invariance
    if (x$udist %in% c("tnormal", "lognormal")) {
      muHvar_c <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
      uHvar_c <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
      vHvar_c <- model.matrix(x$formula, rhs = 4, data = x$dataTable)
      if (invariance == 1) {
        muHvar_p <- apply(muHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[1])
        })
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[1])
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[1])
        })
      } else {
        if (invariance == 2) {
          muHvar_p <- apply(muHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
          })
          uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
          })
        } else {
          if (invariance == 3) {
          muHvar_p <- apply(muHvar_c, 2, function(x) {
            tapply(x, pindex[, 1], mean)
          })
          uHvar_p <- apply(uHvar_c, 2, function(x) {
            tapply(x, pindex[, 1], mean)
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
            tapply(x, pindex[, 1], mean)
          })
          }
        }
      }
    } else {
      uHvar_c <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
      vHvar_c <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
      if (invariance == 1) {
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[1])
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[1])
        })
      } else {
        if (invariance == 2) {
          uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
          })
        } else {
          if (invariance == 3) {
          uHvar_p <- apply(uHvar_c, 2, function(x) {
            tapply(x, pindex[, 1], mean)
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
            tapply(x, pindex[, 1], mean)
          })
          }
        }
      }
    }
    wHvar_c <- x$dataTable$weights
    varw_p <- if (invariance == 1) {
      as.vector(tapply(wHvar_c, pindex[, 1], function(u) u[1]))
    } else {
      if (invariance == 2) {
        as.vector(tapply(wHvar_c, pindex[, 1], function(u) u[length(u)]))
      } else {
        if (invariance == 3) {
          as.vector(tapply(wHvar_c, pindex[, 1], mean))
        }
      }
    }
    wHvar_p <- x$Nid * varw_p/sum(varw_p)
    if (x$modelType %in% c("pl81", "mbc92", "k90")) {
      if (x$udist == "hnormal") {
        hessAnalytical <- eval(parse(text = paste0("phesshalfnormlike_",
          x$modelType, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar, 
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)")))
      } else {
        if (x$udist == "exponential") {
          hessAnalytical <- eval(parse(text = paste0("phessexponormlike_",
          x$modelType, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)")))
        } else {
          if (x$udist == "tnormal") {
          hessAnalytical <- eval(parse(text = paste0("phesstruncnormlike_",
            x$modelType, "(parm = x$mlParam, nXvar = nXvar, 
            nmuZUvar = ncol(muHvar_p), nuZUvar = nuZUvar,nvZVvar = nvZVvar, 
            muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, 
            Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)")))
          } else {
          if (x$udist == "rayleigh") {
            hessAnalytical <- eval(parse(text = paste0("phessraynormlike_",
            x$modelType, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)")))
          } else {
            if (x$udist == "uniform") {
            hessAnalytical <- eval(parse(text = paste0("phessuninormlike_",
              x$modelType, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)")))
            } else {
            if (x$udist == "genexponential") {
              hessAnalytical <- eval(parse(text = paste0("phessgenexponormlike_",
              x$modelType, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)")))
            } else {
              if (x$udist == "tslaplace") {
              hessAnalytical <- eval(parse(text = paste0("phesstslnormlike_",
                x$modelType, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)")))
              } else {
              if (x$udist == "gamma") {
                hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradgammanormlike_",
                x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, N = x$Nid, FiMat = x$FiMat)), unname(x$mlParam))")))
              } else {
                if (x$udist == "weibull") {
                hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradweibullnormlike_",
                  x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, N = x$Nid, FiMat = x$FiMat)), unname(x$mlParam))")))
                } else {
                if (x$udist == "lognormal") {
                  hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                  function(parm) -colSums(pgradlognormlike_",
                  x$modelType, "(parm, 
                  nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar_p), 
                  nuZUvar = ncol(uHvar_p), nvZVvar = ncol(vHvar_p), 
                  muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p, 
                  Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, 
                  S = x$S, wHvar = wHvar_p, N = x$Nid, FiMat = x$FiMat)), 
                  unname(x$mlParam))")))
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
      if (x$modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00", "mols93")) {
        toPaste <- if (x$modelType == "mols93")
          "mols93" else "gzit"
        if (x$udist == "hnormal") {
          hessAnalytical <- eval(parse(text = paste0("phesshalfnormlike_",
          toPaste, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar, 
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, 
          ngZGvar = ngZGvar, gHvar = gHvar)")))
        } else {
          if (x$udist == "exponential") {
          hessAnalytical <- eval(parse(text = paste0("phessexponormlike_",
            toPaste, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar, 
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, 
          ngZGvar = ngZGvar, gHvar = gHvar)")))
          } else {
          if (x$udist == "tnormal") {
            hessAnalytical <- eval(parse(text = paste0("phesstruncnormlike_",
            toPaste, "(parm = x$mlParam, nXvar = nXvar, 
            nmuZUvar = ncol(muHvar_p), nuZUvar = nuZUvar,nvZVvar = nvZVvar, 
            muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, 
            Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, 
            ngZGvar = ngZGvar, gHvar = gHvar)")))
          } else {
            if (x$udist == "rayleigh") {
            hessAnalytical <- eval(parse(text = paste0("phessraynormlike_",
              toPaste, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar, 
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, 
          ngZGvar = ngZGvar, gHvar = gHvar)")))
            } else {
            if (x$udist == "uniform") {
              hessAnalytical <- eval(parse(text = paste0("phessuninormlike_",
              toPaste, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar, 
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, 
          ngZGvar = ngZGvar, gHvar = gHvar)")))
            } else {
              if (x$udist == "genexponential") {
              hessAnalytical <- eval(parse(text = paste0("phessgenexponormlike_",
                toPaste, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar, 
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, 
          ngZGvar = ngZGvar, gHvar = gHvar)")))
              } else {
              if (x$udist == "tslaplace") {
                hessAnalytical <- eval(parse(text = paste0("phesstslnormlike_",
                toPaste, "(parm = x$mlParam, nXvar = nXvar, nuZUvar = nuZUvar, 
          nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, 
          Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, 
          ngZGvar = ngZGvar, gHvar = gHvar)")))
              } else {
                if (x$udist == "gamma") {
                hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(function(parm) -colSums(pgradgammanormlike_",
                  toPaste, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                  nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                  Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, 
                  S = x$S, wHvar = wHvar_p, N = x$Nid, FiMat = x$FiMat, 
                  ngZGvar = ngZGvar, gHvar = gHvar)), unname(x$mlParam))")))
                } else {
                if (x$udist == "weibull") {
                  hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(function(parm) -colSums(pgradweibullnormlike_",
                  toPaste, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                  nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                  Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, 
                  S = x$S, wHvar = wHvar_p, N = x$Nid, FiMat = x$FiMat, 
                  ngZGvar = ngZGvar, gHvar = gHvar)), unname(x$mlParam))")))
                } else {
                  if (x$udist == "lognormal") {
                  hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(function(parm) -colSums(pgradlognormlike_",
                    toPaste, "(parm, nXvar = ncol(Xvar), 
                    nmuZUvar = ncol(muHvar_p), nuZUvar = ncol(uHvar_p), 
                    nvZVvar = ncol(vHvar_p), muHvar = muHvar_p, uHvar = uHvar_p, 
                    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, 
                    TT = x$Vtime, S = x$S, wHvar = wHvar_p, N = x$Nid, 
                    FiMat = x$FiMat, ngZGvar = ngZGvar, gHvar = gHvar)), 
                    unname(x$mlParam))")))
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
    return(invHess * x$Nid)
  }
}

# Gradients Evaluated at each cross section ----------
#' @rdname sfapanel1
#' @export
estfun.sfapanel1 <- function(x, ...) {
  return(x$gradL_OBS)
}
