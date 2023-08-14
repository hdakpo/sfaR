################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Model: Latent Stochastic Frontier Analysis                                   #
# First generation panel models                                                #
# Data: Panel data                                                             #
#------------------------------------------------------------------------------#

#' Latent class stochastic frontier using panel data
#' 
#' \code{\link{sfalcmpanel}} is a symbolic formula based function for the
#' estimation of the latent class stochastic frontier model (LCM) in the case
#' of panel data. Several panel models are implemented: the time invariant 
#' inefficiency model discussed in Pitt and Lee (1981), the time varying 
#' efficiency models suggested in Kumbhakar (1990), Battese and Coelli (1992), 
#' Cuesta (2000), Cuesta and Orea (2002), Kumbhakar and Wang (2005), 
#' Alvarez \emph{et al.} (2006), Feng and Serletis (2009). We also suggested a 
#' modified version of the model by Lee and Schmidt (1993) - see 
#' \sQuote{Details} section.The model is estimated using maximum likelihood (ML). 
#' See Orea and Kumbhakar (2004), Parmeter and Kumbhakar (2014, p282).
#'
#' Only the half-normal distribution is possible for the one-sided error term.
#' Eleven optimization algorithms are available.
#'
#' Depending on the specification, the function accounts for heteroscedasticity 
#' in both one-sided and two-sided error terms as in 
#' Reifschneider and Stevenson (1991), Caudill and Ford (1993), 
#' Caudill \emph{et al.} (1995) and Hadri (1999).
#' Alvarez \emph{et al.} (2006) implements a version of the time varying 
#' inefficiency using the scaling property as in Wang and Schmidt (2002).
#'
#' The model can estimate up to five classes.
#'
#' @aliases sfalcmpanel bread.sfalcmpanel estfun.sfalcmpanel print.sfalcmpanel
#'
#' @param formula A symbolic description of the model to be estimated based on
#' the generic function \code{formula} (see section \sQuote{Details}).
#' @param uhet A one-part formula to account for heteroscedasticity in the
#' one-sided error variance (see section \sQuote{Details}).
#' @param vhet A one-part formula to account for heteroscedasticity in the
#' two-sided error variance (see section \sQuote{Details}).
#' @param thet A one-part formula to account for technological heterogeneity in
#' the construction of the classes.
#' @param logDepVar Logical. Informs whether the dependent variable is logged
#' (\code{TRUE}) or not (\code{FALSE}). Default = \code{TRUE}.
#' @param data The data frame containing the data.
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
#' @param udist Character string. Distribution specification for the one-sided
#' error term. Only the half normal distribution \code{'hnormal'} (Aigner
#' \emph{et al.}) is currently implemented. 
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
#' @param lcmClasses Number of classes to be estimated (default = \code{2}). A
#' maximum of five classes can be estimated.
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
#' @param x an object of class sfalcmpanel (returned by the function 
#' \code{\link{sfalcmpanel}}).
#' @param ... additional arguments of frontier are passed to sfalcmpanel; 
#' additional arguments of the print, bread, estfun, nobs methods are currently
#'  ignored.
#'  
#' @details
#' LCM is an estimation of a finite mixture of production functions. In the case 
#' of the time invariant inefficiency model (Pitt and Lee (1981): \code{'pl81'}) 
#' the production frontier is defined as: 
#'
#'\deqn{y_{it} = \alpha_j + \mathbf{x_{it}^{\prime}} 
#' \bm{\beta_j} + v_{it|j} - Su_{i|j}}
#' 
#' \deqn{\epsilon_{it|j} = v_{it|j} - Su_{i|j}}
#' 
#' The contribution of cross-section \eqn{i} to the likelihood conditional on
#' class \eqn{j} is defined as: 
#' 
#' The prior probability of using a particular technology can depend on some
#' covariates (namely the variables separating the observations into classes)
#' using a logit specification: 
#' 
#' \deqn{\pi(i,j) = \frac{\exp{(\bm{\theta}_j'\mathbf{Z}_{hi})}}{
#' \sum_{m=1}^{J}\exp{(\bm{\theta}_m'\mathbf{Z}_{hi})}}}
#' 
#' with \eqn{\mathbf{Z}_h} the covariates, \eqn{\bm{\theta}} the coefficients 
#' estimated for the covariates, and \eqn{\exp(\bm{\theta}_J'\mathbf{Z}_h)=1}.
#' 
#' The particularity of the prior probability is that it is time-invariant. In 
#' the case of \code{'sfalcmpanel'}, the \eqn{\mathbf{Z}_h} are averaged for 
#' each cross-section.
#'
#' @return \code{\link{sfalcmpanel}} returns a list of class 
#' \code{'sfalcmpanel'} containing the following elements:
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
#' \item{nZHvar}{Number of variables in the logit specification of the finite
#' mixture model (i.e. number of covariates).}
#' 
#'  \item{nuZUvar}{Number of variables explaining heteroscedasticity in the
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
#'  \item{dataTable}{A data frame containing information on data
#' used for optimization along with residuals and fitted values of the OLS and
#' ML estimations, and the individual observation log-likelihood. When
#' \code{weights} is specified an additional variable is also provided in 
#' \code{dataTable}.}
#' 
#'  \item{modelType}{Spefication used for the panel model.}
#' 
#' \item{gHvar}{Matrix of variables used in the case of time varying 
#' inefficiency when \code{'modelType %in% c('bc92a', 'bc92b', 'bc92c', 'kw05',
#' 'c00', 'mols93')'} (for internal use).}
#' 
#' \item{invariance}{Methodology used to obtain time invariant exogenous 
#' variables (in the case of heteroscedasticity).}
#' 
#' \item{initHalf}{When \code{start = NULL} and \code{whichStart == 2L}. 
#' Initial ML estimation with half normal distribution for the one-sided error 
#' term. Model to construct the starting values for 
#' the latent class estimation. Object of class \code{'maxLik'} and 
#' \code{'maxim'} returned.}
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
#' \item{nClasses}{The number of classes estimated.}
#' 
#' \item{mlLoglik}{Log-likelihood value of the M(S)L estimation.}
#'
#' \item{mlParam}{Parameters obtained from M(S)L estimation.}
#' 
#' \item{mlParamMatrix}{Double. Matrix of ML parameters by class.}
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
# @author K Hervé Dakpo
#'
#' @seealso \code{\link[=print.sfalcmpanel]{print}} for printing 
#' \code{sfalcmpanel} object.
#' 
#' \code{\link[=summary.sfalcmpanel]{summary}} for creating and printing
#' summary results.
#'
#' \code{\link[=coef.sfalcmpanel]{coef}} for extracting coefficients of the
#' estimation.
#'
#' \code{\link[=efficiencies.sfalcmpanel]{efficiencies}} for computing
#' (in-)efficiency estimates.
#'
#' \code{\link[=fitted.sfalcmpanel]{fitted}} for extracting the fitted frontier
#' values.
#'
#' \code{\link[=ic.sfalcmpanel]{ic}} for extracting information criteria.
#'
#' \code{\link[=logLik.sfalcmpanel]{logLik}} for extracting log-likelihood
#' value(s) of the estimation.
#'
#' \code{\link[=residuals.sfalcmpanel]{residuals}} for extracting residuals of 
#' the estimation.
#'
#' \code{\link[=vcov.sfalcmpanel]{vcov}} for computing the variance-covariance
#' matrix of the coefficients.
#' 
#' \code{\link[=bread.sfalcmpanel]{bread}} for bread for sandwich estimator.
#' 
#' \code{\link[=estfun.sfalcmpanel]{estfun}} for gradient extraction for each 
#' observation.
#'
#' @references Aigner, D., Lovell, C. A. K., and P. Schmidt. 1977. Formulation
#' and estimation of stochastic frontier production function models.
#' \emph{Journal of Econometrics}, \bold{6}(1), 21--37.
#'
#' Caudill, S. Aigner, D., Lovell, C. A. K., and Schmidt, P. 1977. Formulation
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
#' Hadri, K. 1999. Estimation of a doubly heteroscedastic stochastic frontier
#' cost function. \emph{Journal of Business & Economic Statistics},
#' \bold{17}(3), 359--363.
#' 
#' Kumbhakar, S. C. (1990). Production Frontiers, Panel Data, and Time-Varying 
#' Technical Inefficiency. \emph{Journal of Econometrics}, \bold{46}(1-2), 
#' 201-211. doi: Doi 10.1016/0304-4076(90)90055-X.
#' 
#' Kumbhakar, S. C., & Wang, H.-J. (2005). Estimation of growth convergence 
#' using a stochastic production frontier approach. 
#' \emph{Economics Letters}, \bold{88}(3), 300-305. 
#' doi: https://doi.org/10.1016/j.econlet.2005.01.023.
#' 
#' Orea, L., and S.C. Kumbhakar. 2004. Efficiency measurement using a latent
#' class stochastic frontier model. \emph{Empirical Economics}, \bold{29},
#' 169--183.
#' 
#' Parmeter, C.F., and S.C. Kumbhakar. 2014. Efficiency analysis: A primer on
#' recent advances. \emph{Foundations and Trends in Econometrics}, \bold{7},
#' 191--385.
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
#' Wang, H.J., and Schmidt, P. 2002. One-step and two-step estimation of the
#' effects of exogenous variables on technical efficiency levels. \emph{Journal
#' of Productivity Analysis}, \bold{18}:129--144.
#' 
#'
#' @keywords models optimize panel-data latent-class likelihood
#'
#'
#' @examples
#' 
#' 
#' 
#' @export
sfalcmpanel <- function(formula, uhet, vhet, thet, logDepVar = TRUE, data, idVar = NULL,
  timeVar = NULL, subset, weights, wscale = TRUE, S = 1L, modelType = "bc92a",
  udist = "hnormal", start = NULL, randStart = FALSE, whichStart = 2L, initAlg = "nm",
  initIter = 500, invariance = 2L, lcmClasses = 2, method = "bfgs", hessianType = 1,
  itermax = 2000L, printInfo = FALSE, tol = 1e-12, gradtol = 1e-06, stepmax = 0.1,
  qac = "marquardt") {
  # panel model check -------
  modelType <- tolower(modelType)
  if (!(modelType %in% c("pl81", "bc92a", "bc92b", "bc92c", "mbc92", "k90", "kw05",
    "c00", "mols93"))) {
    stop("Unknown SFA panel model: ", paste(modelType), call. = FALSE)
  }
  # u distribution check -------
  udist <- tolower(udist)
  if (udist != "hnormal") {
    stop("Currently latent class model only handles half-normal distribution ... ",
      call. = FALSE)
  }
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
    if (length(attr(terms(uhet), "term.labels")) == 0) {
      stop("at least one exogeneous variable must be provided for the scaling option for modelType = 'bc92c'",
           call. = FALSE)
    } else {
      ghet <- plhsCheck_u_bc92c(formula = uhet)
      uhet <- ~1
    }
  }
  if (!missing(thet)) {
    thet <- clhsCheck_t(formula = thet)
  } else {
    thet <- ~1
  }
  if (modelType == "bc92c") {
    formula <- formDist_sfalcmpanel_bc92c(formula = formula, uhet = uhet, vhet = vhet,
      ghet = ghet, thet = thet)
  } else {
    formula <- formDist_sfalcmcross(formula = formula, uhet = uhet, vhet = vhet,
      thet = thet)
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
    mtZ <- delete.response(terms(formula, data = data, rhs = 5))
    Zvar <- model.matrix(mtZ, mc)
    Zvar <- Zvar[validObs, , drop = FALSE]
    nZHvar <- ncol(Zvar)
  } else {
    mtZ <- delete.response(terms(formula, data = data, rhs = 4))
    Zvar <- model.matrix(mtZ, mc)
    Zvar <- Zvar[validObs, , drop = FALSE]
    nZHvar <- ncol(Zvar)
  }
  # for logit specification use mean Zit
  Zvar <- apply(Zvar, 2, function(x) tapply(x, index(data)[[1]], mean))
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
    "Latent Class Production/Profit Frontier, e = v - u"
  } else {
    "Latent Class Cost Frontier, e = v + u"
  }
  if (length(logDepVar) != 1 || !is.logical(logDepVar[1])) {
    stop("argument 'logDepVar' must be a single logical value", call. = FALSE)
  }
  if (!(lcmClasses %in% 2:5)) {
    stop("argument 'lcmClasses' must be comprised between 2 and 5", call. = FALSE)
  }
  # Number of parameters -------
  if (modelType == "pl81") {
    nParm <- lcmClasses * (nXvar + nuZUvar + nvZVvar) + (lcmClasses - 1) * nZHvar
  } else {
    if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00")) {
      nParm <- lcmClasses * (nXvar + ngZGvar + nvZVvar + 1) + (lcmClasses -
        1) * nZHvar
    } else {
      if (modelType %in% c("k90", "mbc92")) {
        nParm <- lcmClasses * (nXvar + nuZUvar + nvZVvar + 2) + (lcmClasses -
          1) * nZHvar
      } else {
        if (modelType == "mols93") {
          nParm <- lcmClasses * (nXvar + ngZGvar + nvZVvar) + (lcmClasses -
          1) * nZHvar
        }
      }
    }
  }
  # Checking starting values when provided -------
  if (!is.null(start)) {
    if (length(start) != nParm) {
      stop("Wrong number of initial values: model has ", nParm, " parameters",
        call. = FALSE)
    }
  }
  if (nParm > N) {
    stop("Model has more parameters than observations", call. = FALSE)
  }
  # check whichStart
  if (length(whichStart) != 1 || !(whichStart %in% c(1L, 2L))) {
    stop("argument 'whichStart' must equal either 1 or 2", call. = FALSE)
  }
  # Set std. Error when random start is allowed -------
  sdStart <- if (randStart)
    0.01 else NULL
  # Check algorithms -------
  method <- tolower(method)
  if (!(method %in% c("ucminf", "bfgs", "bhhh", "nr", "nm", "cg", "sann", "sr1",
    "mla", "sparse", "nlminb"))) {
    stop("Unknown or non-available optimization algorithm: ", paste(method),
      call. = FALSE)
  }
  # Check hessian type -------
  if (length(hessianType) != 1 || !(hessianType %in% c(1L, 2L))) {
    stop("argument 'hessianType' must equal either 1 or 2", call. = FALSE)
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
  if (!is.numeric(itermax) || length(itermax) != 1) {
    stop("argument 'itermax' must be a single numeric scalar", call. = FALSE)
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
  dataTable <- data[, names(pindex)]
  dataTable <- cbind(dataTable, data[, all.vars(terms(formula))], weights = wHvar_c)
  dataTable <- cbind(dataTable, olsResiduals = residuals(olsRes), olsFitted = fitted(olsRes))
  # possibility to have duplicated columns if ID or TIME appears in ols
  dataTable <- dataTable[!duplicated(as.list(dataTable))]
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
  if (modelType %in% c("pl81", "k90", "mbc92")) {
    FunArgs <- list(start = start, randStart = randStart, sdStart = sdStart,
      olsParam = olsParam, dataTable = dataTable, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar_c = uHvar_c, vHvar_c = vHvar_c, uHvar_p = uHvar_p,
      vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar, Zvar = Zvar, nZHvar = nZHvar,
      wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S, pindex = pindex, TT = TT,
      method = method, printInfo = printInfo, itermax = itermax, stepmax = stepmax,
      tol = tol, gradtol = gradtol, hessianType = hessianType, initAlg = initAlg,
      initIter = initIter, whichStart = whichStart, qac = qac)
  } else {
    if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00")) {
      FunArgs <- list(start = start, randStart = randStart, sdStart = sdStart,
        olsParam = olsParam, dataTable = dataTable, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar_c = uHvar_c, vHvar_c = vHvar_c, uHvar_p = uHvar_p,
        vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar, Zvar = Zvar, nZHvar = nZHvar,
        wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S, modelType = modelType,
        ngZGvar = ngZGvar, gHvar = gHvar, pindex = pindex, TT = TT, method = method,
        printInfo = printInfo, itermax = itermax, stepmax = stepmax, tol = tol,
        gradtol = gradtol, hessianType = hessianType, initAlg = initAlg,
        initIter = initIter, whichStart = whichStart, qac = qac)

    } else {
      if (modelType == "mols93") {
        FunArgs <- list(start = start, randStart = randStart, sdStart = sdStart,
          olsParam = olsParam, dataTable = dataTable, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar_c = uHvar_c, vHvar_c = vHvar_c, uHvar_p = uHvar_p,
          vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar, Zvar = Zvar, nZHvar = nZHvar,
          wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S, ngZGvar = ngZGvar,
          gHvar = gHvar, pindex = pindex, TT = TT, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, initAlg = initAlg, initIter = initIter,
          whichStart = whichStart, qac = qac)

      }
    }
  }
  ## MLE run -------
  if (modelType == "pl81") {
    mleList <- tryCatch(switch(as.character(lcmClasses), `2` = do.call(LCM2ChnormAlgOpt_pl81,
      FunArgs), `3` = do.call(LCM3ChnormAlgOpt_pl81, FunArgs), `4` = do.call(LCM4ChnormAlgOpt_pl81,
      FunArgs), `5` = do.call(LCM5ChnormAlgOpt_pl81, FunArgs)), error = function(e) print(e))
  } else {
    if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00")) {
      mleList <- tryCatch(switch(as.character(lcmClasses), `2` = do.call(LCM2ChnormAlgOpt_gzit,
        FunArgs), `3` = do.call(LCM3ChnormAlgOpt_gzit, FunArgs), `4` = do.call(LCM4ChnormAlgOpt_gzit,
        FunArgs), `5` = do.call(LCM5ChnormAlgOpt_gzit, FunArgs)), error = function(e) print(e))
    } else {
      if (modelType == "k90") {
        mleList <- tryCatch(switch(as.character(lcmClasses), `2` = do.call(LCM2ChnormAlgOpt_k90,
          FunArgs), `3` = do.call(LCM3ChnormAlgOpt_k90, FunArgs), `4` = do.call(LCM4ChnormAlgOpt_k90,
          FunArgs), `5` = do.call(LCM5ChnormAlgOpt_k90, FunArgs)), error = function(e) print(e))
      } else {
        if (modelType == "mbc92") {
          mleList <- tryCatch(switch(as.character(lcmClasses), `2` = do.call(LCM2ChnormAlgOpt_mbc92,
          FunArgs), `3` = do.call(LCM3ChnormAlgOpt_mbc92, FunArgs), `4` = do.call(LCM4ChnormAlgOpt_mbc92,
          FunArgs), `5` = do.call(LCM5ChnormAlgOpt_mbc92, FunArgs)), error = function(e) print(e))
        } else {
          if (modelType == "mols93") {
          mleList <- tryCatch(switch(as.character(lcmClasses), `2` = do.call(LCM2ChnormAlgOpt_mols93,
            FunArgs), `3` = do.call(LCM3ChnormAlgOpt_mols93, FunArgs),
            `4` = do.call(LCM4ChnormAlgOpt_mols93, FunArgs), `5` = do.call(LCM5ChnormAlgOpt_mols93,
            FunArgs)), error = function(e) print(e))
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
    names(mleList$startVal) <- fName_sfalcmpanel(Xvar = Xvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Zvar = Zvar, nZHvar = nZHvar, gHvar = if (modelType %in%
        c("c00", "bc92c", "mols93"))
        gHvar else NULL, modelType = modelType, lcmClasses = lcmClasses)
    names(mleList$mlParam) <- names(mleList$startVal)
  }
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
  names(mleList$gradient) <- names(mleList$mlParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
  # Return object -------
  mleDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
  dataTable$mlResiduals_c1 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$mlFitted_c1 <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$mlResiduals_c2 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(nXvar +
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]), t(Xvar)))
  dataTable$mlFitted_c2 <- as.numeric(crossprod(matrix(mleList$mlParam[(nXvar +
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]), t(Xvar)))
  if (lcmClasses == 3) {
    dataTable$mlResiduals_c3 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(2 *
      nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 *
      nvZVvar)]), t(Xvar)))
    dataTable$mlFitted_c3 <- as.numeric(crossprod(matrix(mleList$mlParam[(2 *
      nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 *
      nvZVvar)]), t(Xvar)))
  } else {
    if (lcmClasses == 4) {
      dataTable$mlResiduals_c3 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(2 *
        nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
        2 * nvZVvar)]), t(Xvar)))
      dataTable$mlFitted_c3 <- as.numeric(crossprod(matrix(mleList$mlParam[(2 *
        nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
        2 * nvZVvar)]), t(Xvar)))
      dataTable$mlResiduals_c4 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(3 *
        nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
        3 * nvZVvar)]), t(Xvar)))
      dataTable$mlFitted_c4 <- as.numeric(crossprod(matrix(mleList$mlParam[(3 *
        nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
        3 * nvZVvar)]), t(Xvar)))
    } else {
      if (lcmClasses == 5) {
        dataTable$mlResiduals_c3 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(2 *
          nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
          2 * nvZVvar)]), t(Xvar)))
        dataTable$mlFitted_c3 <- as.numeric(crossprod(matrix(mleList$mlParam[(2 *
          nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
          2 * nvZVvar)]), t(Xvar)))
        dataTable$mlResiduals_c4 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(3 *
          nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
          3 * nvZVvar)]), t(Xvar)))
        dataTable$mlFitted_c4 <- as.numeric(crossprod(matrix(mleList$mlParam[(3 *
          nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
          3 * nvZVvar)]), t(Xvar)))
        dataTable$mlResiduals_c5 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(4 *
          nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
          4 * nvZVvar)]), t(Xvar)))
        dataTable$mlFitted_c5 <- as.numeric(crossprod(matrix(mleList$mlParam[(4 *
          nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
          4 * nvZVvar)]), t(Xvar)))
      }
    }
  }
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
  returnObj$nZHvar <- nZHvar
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  returnObj$logDepVar <- logDepVar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- mleList$startVal
  returnObj$dataTable <- dataTable
  returnObj$modelType <- modelType
  if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00", "mols93")) {
    returnObj$gHvar <- gHvar
  }
  returnObj$invariance <- invariance
  if (is.null(start)) {
    if (whichStart == 2L) {
      returnObj$initHalf <- mleList$initHalf
    }
  }
  returnObj$isWeights <- !all.equal(wHvar_p, rep(1, N))
  returnObj$optType <- mleList$type
  returnObj$nIter <- mleList$nIter
  returnObj$optStatus <- mleList$status
  returnObj$startLoglik <- mleList$startLoglik
  returnObj$nClasses <- lcmClasses
  returnObj$mlLoglik <- mleList$mleLoglik
  returnObj$mlParam <- mleList$mlParam
  returnObj$mlParamMatrix <- rbind(matrix(mleList$mlParam[1:(lcmClasses * (nXvar +
    nuZUvar + nvZVvar))], ncol = lcmClasses), cbind(matrix(mleList$mlParam[(lcmClasses *
    (nXvar + nuZUvar + nvZVvar) + 1):(lcmClasses * (nXvar + nuZUvar + nvZVvar) +
    (lcmClasses - 1) * nZHvar)], ncol = lcmClasses - 1), NA))
  colnames(returnObj$mlParamMatrix) <- paste0("Class", 1:lcmClasses)
  d1names <- names(mleList$mlParam)[c(1:(nXvar + nuZUvar + nvZVvar), (lcmClasses *
    (nXvar + nuZUvar + nvZVvar) + 1):(lcmClasses * (nXvar + nuZUvar + nvZVvar) +
    nZHvar))]
  rownames(returnObj$mlParamMatrix) <- gsub("Cl1", "Cl", d1names)
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
  returnObj$mlDate <- mleDate
  rm(mleList)
  class(returnObj) <- "sfalcmpanel"
  return(returnObj)
}

# print for sfalcmpanel ----------
#' @rdname sfalcmpanel
#' @export
print.sfalcmpanel <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call, width.cutoff = 500))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat("Normal-Half Normal Panel Latent Class Stochastic Frontier Model", "\n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}

# Bread for Sandwich Estimator ----------
#' @rdname sfalcmpanel
#' @export
bread.sfalcmpanel <- function(x, ...) {
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
    if (x$modelType == "bc92c") {
      Zvar <- model.matrix(x$formula, rhs = 5, data = x$dataTable)
    } else {
      Zvar <- model.matrix(x$formula, rhs = 4, data = x$dataTable)
    }
    Zvar <- apply(Zvar, 2, function(x) tapply(x, pindex[, 1], mean))
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
    if (x$nClasses == 2) {
      if (x$modelType %in% c("pl81", "mbc92", "k90")) {
        hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradLCMhalfnormlike2C_",
          x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = ncol(Zvar))), unname(x$mlParam))")))
      } else {
        if (x$modelType %in% c("bc92a", "bc92b", "bc92c", "kw05", "c00",
          "mols93")) {
          hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradLCMhalfnormlike2C_",
          x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar, Zvar = Zvar,
        nZHvar = ncol(Zvar))), unname(x$mlParam))")))
        }
      }
    } else {
      if (x$nClasses == 3) {
        if (x$modelType %in% c("pl81", "mbc92", "k90")) {
          hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradLCMhalfnormlike3C_",
          x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = ncol(Zvar))), unname(x$mlParam))")))
        } else {
          if (x$modelType %in% c("bc92a", "bc92b", "bc93C", "kw05", "c00",
          "mols93")) {
          hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradLCMhalfnormlike3C_",
            x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar, Zvar = Zvar,
        nZHvar = ncol(Zvar))), unname(x$mlParam))")))
          }
        }
      } else {
        if (x$nClasses == 4) {
          if (x$modelType %in% c("pl81", "mbc92", "k90")) {
          hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradLCMhalfnormlike4C_",
            x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = ncol(Zvar))), unname(x$mlParam))")))
          } else {
          if (x$modelType %in% c("bc92a", "bc92b", "bc94C", "kw05", "c00",
            "mols93")) {
            hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradLCMhalfnormlike4C_",
            x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar, Zvar = Zvar,
        nZHvar = ncol(Zvar))), unname(x$mlParam))")))
          }
          }
        } else {
          if (x$nClasses == 5) {
          if (x$modelType %in% c("pl81", "mbc92", "k90")) {
            hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradLCMhalfnormlike5C_",
            x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = ncol(Zvar))), unname(x$mlParam))")))
          } else {
            if (x$modelType %in% c("bc92a", "bc92b", "bc95C", "kw05", "c00",
            "mols93")) {
            hessAnalytical <- eval(parse(text = paste0("calculus::jacobian(
                function(parm) -colSums(pgradLCMhalfnormlike5C_",
              x$modelType, "(parm, nXvar = ncol(Xvar), nuZUvar = ncol(uHvar_p), 
                nvZVvar = ncol(vHvar_p), uHvar = uHvar_p, vHvar = vHvar_p, 
                Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = x$Vtime, S = x$S, 
                wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar, Zvar = Zvar,
        nZHvar = ncol(Zvar))), unname(x$mlParam))")))
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

# Gradients Evaluated at each Observation ----------
#' @rdname sfalcmpanel
#' @export
estfun.sfalcmpanel <- function(x, ...) {
  return(x$gradL_OBS)
}
