################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Model: Latent Class Stochastic Frontier Analysis                             #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Latent class stochastic frontier using cross-sectional data
#'
#' @description
#' \code{\link{sfalcmcross}} is a symbolic formula based function for the
#' estimation of the latent class stochastic frontier model (LCM) in the case
#' of cross-sectional or pooled cross-sectional data. The model is estimated
#' using maximum likelihood (ML). See Orea and Kumbhakar (2004), Parmeter and
#' Kumbhakar (2014, p282).
#'
#' Only the half-normal distribution is possible for the one-sided error term.
#' Eleven optimization algorithms are available.
#'
#' The function also accounts for heteroscedasticity in both one-sided and
#' two-sided error terms, as in Reifschneider and Stevenson (1991), Caudill and
#' Ford (1993), Caudill \emph{et al.} (1995) and Hadri (1999).
#'
#' The model can estimate up to five classes.
#'
#' @aliases sfalcmcross print.sfalcmcross
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
#' @param subset An optional vector specifying a subset of observations to be
#' used in the optimization process.
#' @param weights An optional vector of weights to be used for weighted 
#' log-likelihood. Should be \code{NULL} or numeric vector with positive values. 
#' When \code{NULL}, a numeric vector of 1 is used.
#' @param wscale Logical. When \code{weights} is not \code{NULL}, a scaling
#' transformation is used such that the \code{weights} sums to the sample 
#' size. Default \code{TRUE}. When \code{FALSE} no scaling is used.
#' @param S If \code{S = 1} (default), a production (profit) frontier is
#' estimated: \eqn{\epsilon_i = v_i-u_i}. If \code{S = -1}, a cost frontier is
#' estimated: \eqn{\epsilon_i = v_i+u_i}.
#' @param udist Character string. Distribution specification for the one-sided
#' error term. Only the half normal distribution \code{'hnormal'} (Aigner
#' \emph{et al.}, 1977, Meeusen and Vandenbroeck, 1977) is currently
#' implemented.
#' @param start Numeric vector. Optional starting values for the maximum
#' likelihood (ML) estimation.
#' @param whichStart Integer. If \code{'whichStart = 1'}, the starting values 
#' are obtained from the method of moments. When \code{'whichStart = 2'}
#' (Default), the model is initialized by solving the homoscedastic pooled 
#' cross section SFA model.
#' @param initAlg Character string specifying the algorithm used for 
#' initialization and obtain the starting values (when \code{'whichStart = 2'}).
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
#' @param initIter Maximum number of iterations for initialization algorithm.
#' Default \code{100}.
#' @param lcmClasses Number of classes to be estimated (default = \code{2}). A
#' maximum of five classes can be estimated.
##' @param method Optimization algorithm used for the estimation.  Default =
#' \code{'bfgs'}. 11 algorithms are available: \itemize{ \item \code{'bfgs'},
#' for Broyden-Fletcher-Goldfarb-Shanno (see
#' \code{\link[maxLik:maxBFGS]{maxBFGS}}) \item \code{'bhhh'}, for
#' Berndt-Hall-Hall-Hausman (see \code{\link[maxLik:maxBHHH]{maxBHHH}}) \item
#' \code{'nr'}, for Newton-Raphson (see \code{\link[maxLik:maxNR]{maxNR}}) 
#' \item \code{'nm'}, for Nelder-Mead (see \code{\link[maxLik:maxNM]{maxNM}}) 
#' \item \code{'cg'}, for Conjugate Gradient 
#' (see \code{\link[maxLik:maxCG]{maxCG}}) \item \code{'sann'}, for Simulated 
#' Annealing (see \code{\link[maxLik:maxSANN]{maxSANN}})
#' \item \code{'ucminf'}, for a quasi-Newton type optimization with BFGS updating of 
#' the inverse Hessian and soft line search with a trust region type monitoring 
#' of the input to the line search algorithm 
#' (see \code{\link[ucminf:ucminf]{ucminf}})
#' \item \code{'mla'}, for general-purpose optimization based on
#' Marquardt-Levenberg algorithm (see \code{\link[marqLevAlg:mla]{mla}})
#' \item \code{'sr1'}, for Symmetric Rank 1 (see
#' \code{\link[trustOptim:trust.optim]{trust.optim}}) \item \code{'sparse'}, 
#' for trust regions and sparse Hessian 
#' (see \code{\link[trustOptim:trust.optim]{trust.optim}}) \item
#' \code{'nlminb'}, for optimization using PORT routines (see
#' \code{\link[stats:nlminb]{nlminb}})}
#' @param hessianType Integer. If \code{1} (default), analytic Hessian is
#' returned. If \code{2}, bhhh Hessian is estimated (\eqn{g'g}).
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
#' and \code{'nr'} algorithms. If \code{'qac = stephalving'}, the step length
#' is decreased but the direction is kept. If \code{'qac = marquardt'}
#' (default), the step length is decreased while also moving closer to the pure
#' gradient direction. See \code{\link[maxLik:maxBHHH]{maxBHHH}} and
#' \code{\link[maxLik:maxNR]{maxNR}}.
#' @param x an object of class sfalcmcross (returned by the function
#' \code{\link{sfalcmcross}}).
#' @param ... additional arguments of frontier are passed to sfalcmcross; 
#' additional arguments of the print, bread, estfun, nobs methods are currently 
#' ignored.
#'
#' @details
#' LCM is an estimation of a finite mixture of production functions:
#'
#'\deqn{y_i = \alpha_j + \mathbf{x_i^{\prime}} 
#' \bm{\beta_j} + v_{i|j} - Su_{i|j}}
#' 
#' \deqn{\epsilon_{i|j} = v_{i|j} - Su_{i|j}}
#'
#' where \eqn{i} is the observation, \eqn{j} is the class, \eqn{y} is the
#' output (cost, revenue, profit), \eqn{x} is the vector of main explanatory
#' variables (inputs and other control variables), \eqn{u} is the one-sided
#' error term with variance \eqn{\sigma_{u}^2}, and \eqn{v} is the two-sided
#' error term with variance \eqn{\sigma_{v}^2}.
#'
#' \code{S = 1} in the case of production (profit) frontier function and
#' \code{S = -1} in the case of cost frontier function.
#'
#' The contribution of observation \eqn{i} to the likelihood conditional on
#' class \eqn{j} is defined as: 
#' 
#' \deqn{P(i|j) = \frac{2}{\sqrt{\sigma_{u|j}^2 + 
#' \sigma_{v|j}^2}}\phi\left(\frac{S\epsilon_{i|j}}{\sqrt{
#' \sigma_{u|j}^2 +\sigma_{v|j}^2}}\right)\Phi\left(\frac{
#' \mu_{i*|j}}{\sigma_{*|j}}\right)}
#'
#' where 
#' 
#' \deqn{\mu_{i*|j}=\frac{- S\epsilon_{i|j}
#' \sigma_{u|j}^2}{\sigma_{u|j}^2 + \sigma_{v|j}^2}}
#'
#' and 
#' 
#' \deqn{\sigma_*^2 = \frac{\sigma_{u|j}^2 
#' \sigma_{v|j}^2}{\sigma_{u|j}^2 + \sigma_{v|j}^2}}
#' 
#' The prior probability of using a particular technology can depend on some
#' covariates (namely the variables separating the observations into classes)
#' using a logit specification: 
#' 
#' \deqn{\pi(i,j) = \frac{\exp{(\bm{\theta}_j'\mathbf{Z}_{hi})}}{
#' \sum_{m=1}^{J}\exp{(\bm{\theta}_m'\mathbf{Z}_{hi})}}}
#' 
#' with \eqn{\mathbf{Z}_h} the covariates, \eqn{\bm{\theta}} the coefficients estimated for
#' the covariates, and \eqn{\exp(\bm{\theta}_J'\mathbf{Z}_h)=1}.
#'
#' The unconditional likelihood of observation \eqn{i} is simply the average
#' over the \eqn{J} classes:
#'
#' \deqn{P(i) = \sum_{m=1}^{J}\pi(i,m)P(i|m)}
#'
#' The number of classes to retain can be based on information criterion (see
#' for instance \code{\link[=ic.sfalcmcross]{ic}}).
#'
#' Class assignment is based on the largest posterior probability. This
#' probability is obtained using Bayes' rule, as follows for class \eqn{j}:
#' 
#' \deqn{w\left(j|i\right)=\frac{P\left(i|j\right)
#' \pi\left(i,j\right)}{\sum_{m=1}^JP\left(i|m\right)
#' \pi\left(i, m\right)}}
#'
#' To accommodate heteroscedasticity in the variance parameters of the error
#' terms, a single part (right) formula can also be specified. To impose the
#' positivity on these parameters, the variances are modelled respectively as:
#' \eqn{\sigma^2_{u|j} = \exp{(\bm{\delta}_j'\mathbf{Z}_u)}} and \eqn{\sigma^2_{v|j} =
#' \exp{(\bm{\phi}_j'\mathbf{Z}_v)}}, where \eqn{Z_u} and \eqn{Z_v} are the
#' heteroscedasticity variables (inefficiency drivers in the case of \eqn{\mathbf{Z}_u})
#' and \eqn{\bm{\delta}} and \eqn{\bm{\phi}} the coefficients. \code{'sfalcmcross'} only 
#' supports the half-normal distribution for the one-sided error term.
#' 
#' \code{sfalcmcross} allows for the maximization of weighted log-likelihood.
#' When option \code{weights} is specified and \code{wscale = TRUE}, the weights
#' are scaled as: 
#' 
#' \deqn{new_{weights} = sample_{size} \times 
#' \frac{old_{weights}}{\sum(old_{weights})}}
#' 
#' For complex problems, non-gradient methods (e.g. \code{nm} or 
#' \code{sann}) can be used to warm start the optimization and zoom in the 
#' neighborhood of the solution. Then a gradient-based methods is recommended 
#' in the second step. In the case of \code{sann}, we recommend to significantly 
#' increase the iteration limit (e.g. \code{itermax = 20000}). The Conjugate 
#' Gradient (\code{cg}) can also be used in the first stage.
#' 
#' A set of extractor functions for fitted model objects is available for 
#' objects of class \code{'sfalcmcross'} including methods to the generic functions 
#' \code{\link[=print.sfalcmcross]{print}}, 
#' \code{\link[=summary.sfalcmcross]{summary}}, 
#' \code{\link[=coef.sfalcmcross]{coef}}, 
#' \code{\link[=fitted.sfalcmcross]{fitted}}, 
#' \code{\link[=logLik.sfalcmcross]{logLik}}, 
#' \code{\link[=residuals.sfalcmcross]{residuals}}, 
#' \code{\link[=vcov.sfalcmcross]{vcov}}, 
#' \code{\link[=efficiencies.sfalcmcross]{efficiencies}}, 
#' \code{\link[=ic.sfalcmcross]{ic}}, 
#' \code{\link[=marginal.sfalcmcross]{marginal}}, 
#' \code{\link[=estfun.sfalcmcross]{estfun}} and 
#' \code{\link[=bread.sfalcmcross]{bread}} (from the \CRANpkg{sandwich} package), 
#' [lmtest::coeftest()] (from the \CRANpkg{lmtest} package).
#'
#' @return \code{\link{sfalcmcross}} returns a list of class \code{'sfalcmcross'}
#' containing the following elements:
#'
#' \item{call}{The matched call.}
#'
#' \item{formula}{Multi parts formula describing the estimated model.}
#'
#' \item{S}{The argument \code{'S'}. See the section \sQuote{Arguments}.}
#'
#' \item{typeSfa}{Character string. 'Latent Class Production/Profit Frontier, e
#' = v - u' when \code{S = 1} and 'Latent Class Cost Frontier, e = v + u' when
#' \code{S = -1}.}
#'
#' \item{Nobs}{Number of observations used for optimization.}
#'
#' \item{nXvar}{Number of main explanatory variables.}
#'
#' \item{nZHvar}{Number of variables in the logit specification of the finite
#' mixture model (i.e. number of covariates).}
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
#' \item{startVal}{Numeric vector. Starting value for ML estimation.}
#'
#' \item{dataTable}{A data frame (tibble format) containing information on data
#' used for optimization along with residuals and fitted values of the OLS and
#' ML estimations, and the individual observation log-likelihood. When
#' \code{weights} is specified an additional variable is also provided in 
#' \code{dataTable}.}
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
#' \item{optType}{The optimization algorithm used.}
#'
#' \item{nIter}{Number of iterations of the ML estimation.}
#'
#' \item{optStatus}{An optimization algorithm termination message.}
#'
#' \item{startLoglik}{Log-likelihood at the starting values.}
#'
#' \item{nClasses}{The number of classes estimated.}
#'
#' \item{mlLoglik}{Log-likelihood value of the ML estimation.}
#'
#' \item{mlParam}{Numeric vector. Parameters obtained from ML estimation.}
#' 
#' \item{mlParamMatrix}{Double. Matrix of ML parameters by class.}
#'
#' \item{gradient}{Numeric vector. Each variable gradient of the ML
#' estimation.}
#'
#' \item{gradL_OBS}{Matrix. Each variable individual observation gradient of
#' the ML estimation.}
#'
#' \item{gradientNorm}{Numeric. Gradient norm of the ML estimation.}
#'
#' \item{invHessian}{The covariance matrix of the parameters obtained from the
#' ML estimation.}
#'
#' \item{hessianType}{The argument \code{'hessianType'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{mlDate}{Date and time of the estimated model.}
#'
#' @note In the case of panel data, \code{\link{sfalcmcross}} estimates a pooled
#' cross-section where the probability of belonging to a class a priori is not
#' permanent (not fixed over time).
#'
# @author K Hervé Dakpo
#'
#' @seealso \code{\link[=print.sfalcmcross]{print}} for printing \code{sfalcmcross} 
#' object.
#' 
#' \code{\link[=summary.sfalcmcross]{summary}} for creating and printing
#' summary results.
#'
#' \code{\link[=coef.sfalcmcross]{coef}} for extracting coefficients of the
#' estimation.
#'
#' \code{\link[=efficiencies.sfalcmcross]{efficiencies}} for computing
#' (in-)efficiency estimates.
#'
#' \code{\link[=fitted.sfalcmcross]{fitted}} for extracting the fitted frontier
#' values.
#'
#' \code{\link[=ic.sfalcmcross]{ic}} for extracting information criteria.
#'
#' \code{\link[=logLik.sfalcmcross]{logLik}} for extracting log-likelihood
#' value(s) of the estimation.
#'
#' \code{\link[=marginal.sfalcmcross]{marginal}} for computing marginal effects of
#' inefficiency drivers.
#'
#' \code{\link[=residuals.sfalcmcross]{residuals}} for extracting residuals of the
#' estimation.
#'
#' \code{\link[=vcov.sfalcmcross]{vcov}} for computing the variance-covariance
#' matrix of the coefficients.
#' 
#' \code{\link[=bread.sfalcmcross]{bread}} for bread for sandwich estimator.
#' 
#' \code{\link[=estfun.sfalcmcross]{estfun}} for gradient extraction for each 
#' observation.
#'
#' @references Aigner, D., Lovell, C. A. K., and P. Schmidt. 1977. Formulation
#' and estimation of stochastic frontier production function models.
#' \emph{Journal of Econometrics}, \bold{6}(1), 21--37.
#'
#' Caudill, S. B., and J. M. Ford. 1993. Biases in frontier estimation due to
#' heteroscedasticity. \emph{Economics Letters}, \bold{41}(1), 17--20.
#'
#' Caudill, S. B., Ford, J. M., and D. M. Gropper. 1995. Frontier estimation
#' and firm-specific inefficiency measures in the presence of
#' heteroscedasticity. \emph{Journal of Business & Economic Statistics},
#' \bold{13}(1), 105--111.
#'
#' Hadri, K. 1999. Estimation of a doubly heteroscedastic stochastic frontier
#' cost function. \emph{Journal of Business & Economic Statistics},
#' \bold{17}(3), 359--363.
#'
#' Meeusen, W., and J. Vandenbroeck. 1977. Efficiency estimation from
#' Cobb-Douglas production functions with composed error. \emph{International
#' Economic Review}, \bold{18}(2), 435--445.
#'
#' Orea, L., and S.C. Kumbhakar. 2004. Efficiency measurement using a latent
#' class stochastic frontier model. \emph{Empirical Economics}, \bold{29},
#' 169--183.
#'
#' Parmeter, C.F., and S.C. Kumbhakar. 2014. Efficiency analysis: A primer on
#' recent advances. \emph{Foundations and Trends in Econometrics}, \bold{7},
#' 191--385.
#'
#' Reifschneider, D., and R. Stevenson. 1991. Systematic departures from the
#' frontier: A framework for the analysis of firm inefficiency.
#' \emph{International Economic Review}, \bold{32}(3), 715--723.
#'
#' @keywords models optimize cross-section latent-class likelihood
#'
#' @examples
#'
#' ## Using data on eighty-two countries production (GDP)
#' # LCM Cobb Douglas (production function) half normal distribution
#' # Intercept and initStat used as separating variables
#' cb_2c_h1 <- sfalcmcross(formula = ly ~ lk + ll + yr, thet = ~initStat, 
#' data = worldprod)
#' summary(cb_2c_h1)
#'
#' # summary of the initial ML model
#' summary(cb_2c_h1$InitHalf)
#'
#' # Only the intercept is used as the separating variable
#' # and only variable initStat is used as inefficiency driver
#' cb_2c_h3 <- sfalcmcross(formula = ly ~ lk + ll + yr, uhet = ~initStat, 
#' data = worldprod)
#' summary(cb_2c_h3)
#'
#' @export
sfalcmcross <- function(formula, uhet, vhet, thet, logDepVar = TRUE,
  data, subset, weights, wscale = TRUE, S = 1L, udist = "hnormal",
  start = NULL, whichStart = 2L, initAlg = "nm", initIter = 100,
  lcmClasses = 2, method = "bfgs", hessianType = 1, itermax = 2000L,
  printInfo = FALSE, tol = 1e-12, gradtol = 1e-06, stepmax = 0.1,
  qac = "marquardt") {
  # u distribution check -------
  udist <- tolower(udist)
  if (udist != "hnormal") {
    stop("Currently latent class model only handles half-normal distribution ... ",
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
  if (!missing(thet)) {
    thet <- clhsCheck_t(formula = thet)
  } else {
    thet <- ~1
  }
  formula <- formDist_sfalcmcross(formula = formula, uhet = uhet,
    vhet = vhet, thet = thet)
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
  # if subset is non-missing and there NA, force data to
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
      length(Yvar), ") must be the same to the number of observations of the exogenous variables (",
      nrow(Xvar), ")", sep = ""), call. = FALSE)
  }
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
  # Check other supplied options -------
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
    stop("argument 'logDepVar' must be a single logical value",
      call. = FALSE)
  }
  if (!(lcmClasses %in% 2:5)) {
    stop("argument 'lcmClasses' must be comprised between 2 and 5",
      call. = FALSE)
  }
  # Number of parameters -------
  nParm <- lcmClasses * (nXvar + nuZUvar + nvZVvar) + (lcmClasses -
    1) * nZHvar
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
  # check whichStart
  if (length(whichStart) != 1 || !(whichStart %in% c(1L, 2L))) {
    stop("argument 'whichStart' must equal either 1 or 2",
      call. = FALSE)
  }
  # Check algorithms -------
  method <- tolower(method)
  if (!(method %in% c("ucminf", "bfgs", "bhhh", "nr", "nm",
    "cg", "sann", "sr1", "mla", "sparse", "nlminb"))) {
    stop("Unknown or non-available optimization algorithm: ",
      paste(method), call. = FALSE)
  }
  # Check hessian type -------
  if (length(hessianType) != 1 || !(hessianType %in% c(1L,
    2L))) {
    stop("argument 'hessianType' must equal either 1 or 2",
      call. = FALSE)
  }
  # Other optimization options -------
  if (!is.numeric(initIter) || length(initIter) != 1) {
    stop("argument 'initIter' must be a single numeric scalar",
      call. = FALSE)
  }
  if (initIter != round(initIter)) {
    stop("argument 'initIter' must be an integer", call. = FALSE)
  }
  if (initIter <= 0) {
    stop("argument 'initIter' must be positive", call. = FALSE)
  }
  initIter <- as.integer(initIter)
  if (!is.numeric(itermax) || length(itermax) != 1) {
    stop("argument 'itermax' must be a single numeric scalar",
      call. = FALSE)
  }
  initAlg <- tolower(initAlg)
  if (!(initAlg %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann"))) {
    stop("Unknown or non-available optimization algorithm: ",
      paste(initAlg), call. = FALSE)
  }
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
  names(olsRes$coefficients) <- colnames(Xvar)
  olsParam <- c(olsRes$coefficients)
  if (inherits(data, "pdata.frame")) {
    dataTable <- data[, names(index(data))][validObs, ]
  } else {
    dataTable <- data.frame(IdObs = 1:sum(validObs))
  }
  dataTable <- cbind(dataTable, data[, all.vars(terms(formula))][validObs,
    ], weights = wHvar)
  dataTable <- cbind(dataTable, olsResiduals = residuals(olsRes),
    olsFitted = fitted(olsRes))
  # possibility to have duplicated columns if ID or TIME
  # appears in ols in the case of panel data
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
  # arguments -------
  FunArgs <- list(start = start, olsParam = olsParam, dataTable = dataTable,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Zvar = Zvar,
    nZHvar = nZHvar, Xvar = Xvar, S = S, wHvar = wHvar, method = method,
    initAlg = initAlg, initIter = initIter, whichStart = whichStart,
    printInfo = printInfo, itermax = itermax, stepmax = stepmax,
    tol = tol, gradtol = gradtol, hessianType = hessianType,
    qac = qac)
  ## MLE run -------
  mleList <- tryCatch(switch(as.character(lcmClasses), `2` = do.call(LCM2ChnormAlgOpt,
    FunArgs), `3` = do.call(LCM3ChnormAlgOpt, FunArgs), `4` = do.call(LCM4ChnormAlgOpt,
    FunArgs), `5` = do.call(LCM5ChnormAlgOpt, FunArgs)),
    error = function(e) print(e))
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
          list(type = "Levenberg-Marquardt maximization",
          nIter = mleList$mleObj$ni, status = switch(mleList$mleObj$istop,
            `1` = "convergence criteria were satisfied",
            `2` = "maximum number of iterations was reached",
            `4` = "algorithm encountered a problem in the function computation"),
          mleLoglik = -mleList$mleObj$fn.value, gradient = -mleList$mleObj$grad)
        } else {
          if (method == "sparse") {
          list(type = "Sparse Hessian maximization",
            nIter = mleList$mleObj$iterations, status = mleList$mleObj$status,
            mleLoglik = -mleList$mleObj$fval, gradient = -mleList$mleObj$gradient)
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
  if (!is.null(start)) {
    names(mleList$startVal) <- fName_sfalcmcross(Xvar = Xvar,
      uHvar = uHvar, vHvar = vHvar, Zvar = Zvar, nZHvar = nZHvar,
      lcmClasses = lcmClasses)
    names(mleList$mlParam) <- names(mleList$startVal)
  }
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
  names(mleList$gradient) <- names(mleList$mlParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
  # Return object -------
  mleDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M")
  dataTable$mlResiduals_c1 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$mlFitted_c1 <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$mlResiduals_c2 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(nXvar +
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]),
    t(Xvar)))
  dataTable$mlFitted_c2 <- as.numeric(crossprod(matrix(mleList$mlParam[(nXvar +
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]),
    t(Xvar)))
  if (lcmClasses == 3) {
    dataTable$mlResiduals_c3 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(2 *
      nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
      2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
    dataTable$mlFitted_c3 <- as.numeric(crossprod(matrix(mleList$mlParam[(2 *
      nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
      2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
  } else {
    if (lcmClasses == 4) {
      dataTable$mlResiduals_c3 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(2 *
        nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
        2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
      dataTable$mlFitted_c3 <- as.numeric(crossprod(matrix(mleList$mlParam[(2 *
        nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
        2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
      dataTable$mlResiduals_c4 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(3 *
        nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar +
        3 * nuZUvar + 3 * nvZVvar)]), t(Xvar)))
      dataTable$mlFitted_c4 <- as.numeric(crossprod(matrix(mleList$mlParam[(3 *
        nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar +
        3 * nuZUvar + 3 * nvZVvar)]), t(Xvar)))
    } else {
      if (lcmClasses == 5) {
        dataTable$mlResiduals_c3 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(2 *
          nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
          nXvar + 2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
        dataTable$mlFitted_c3 <- as.numeric(crossprod(matrix(mleList$mlParam[(2 *
          nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
          nXvar + 2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
        dataTable$mlResiduals_c4 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(3 *
          nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 *
          nXvar + 3 * nuZUvar + 3 * nvZVvar)]), t(Xvar)))
        dataTable$mlFitted_c4 <- as.numeric(crossprod(matrix(mleList$mlParam[(3 *
          nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 *
          nXvar + 3 * nuZUvar + 3 * nvZVvar)]), t(Xvar)))
        dataTable$mlResiduals_c5 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(4 *
          nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 *
          nXvar + 4 * nuZUvar + 4 * nvZVvar)]), t(Xvar)))
        dataTable$mlFitted_c5 <- as.numeric(crossprod(matrix(mleList$mlParam[(4 *
          nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 *
          nXvar + 4 * nuZUvar + 4 * nvZVvar)]), t(Xvar)))
      }
    }
  }
  dataTable$logL_OBS <- mleList$mleObj$logL_OBS
  returnObj <- list()
  returnObj$call <- cl
  returnObj$formula <- formula
  returnObj$S <- S
  returnObj$typeSfa <- typeSfa
  returnObj$Nobs <- N
  returnObj$nXvar <- nXvar
  returnObj$nZHvar <- nZHvar
  returnObj$logDepVar <- logDepVar
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- mleList$startVal
  returnObj$dataTable <- dataTable
  if (is.null(start)) {
    if (whichStart == 2L) {
      returnObj$initHalf <- mleList$initHalf
    }
  }
  returnObj$isWeights <- !all.equal(wHvar, rep(1, N))
  returnObj$optType <- mleList$type
  returnObj$nIter <- mleList$nIter
  returnObj$optStatus <- mleList$status
  returnObj$startLoglik <- mleList$startLoglik
  returnObj$nClasses <- lcmClasses
  returnObj$mlLoglik <- mleList$mleLoglik
  returnObj$mlParam <- mleList$mlParam
  returnObj$mlParamMatrix <- rbind(matrix(mleList$mlParam[1:(lcmClasses *
    (nXvar + nuZUvar + nvZVvar))], ncol = lcmClasses), cbind(matrix(mleList$mlParam[(lcmClasses *
    (nXvar + nuZUvar + nvZVvar) + 1):(lcmClasses * (nXvar +
    nuZUvar + nvZVvar) + (lcmClasses - 1) * nZHvar)], ncol = lcmClasses -
    1), NA))
  colnames(returnObj$mlParamMatrix) <- paste0("Class", 1:lcmClasses)
  d1names <- names(mleList$mlParam)[c(1:(nXvar + nuZUvar +
    nvZVvar), (lcmClasses * (nXvar + nuZUvar + nvZVvar) +
    1):(lcmClasses * (nXvar + nuZUvar + nvZVvar) + nZHvar))]
  rownames(returnObj$mlParamMatrix) <- gsub("Cl1", "Cl", d1names)
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
  returnObj$mlDate <- mleDate
  rm(mleList)
  class(returnObj) <- "sfalcmcross"
  return(returnObj)
}

# print for sfalcmcross ----------
#' @rdname sfalcmcross
#' @export
print.sfalcmcross <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call, width.cutoff = 500))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat("Normal-Half Normal Latent Class Stochastic Frontier Model",
    "\n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}

# Bread for Sandwich Estimator ----------
#' @rdname sfalcmcross
#' @export
bread.sfalcmcross <- function(x, ...) {
  if (x$hessianType == "Analytic Hessian") {
    return(x$invHessian * x$Nobs)
  } else {
    cat("Computing Analytical Hessian \n")
    Yvar <- model.response(model.frame(x$formula, data = x$dataTable))
    Xvar <- model.matrix(x$formula, rhs = 1, data = x$dataTable)
    uHvar <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
    vHvar <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
    Zvar <- model.matrix(x$formula, rhs = 4, data = x$dataTable)
    if (x$nClasses == 2) {
      hessAnalytical <- chessLCMhalfnormlike2C(parm = x$mlParam,
        nXvar = ncol(Xvar), nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = x$dataTable$weights, S = x$S, Zvar = Zvar,
        nZHvar = ncol(Zvar))
    } else {
      if (x$nClasses == 3) {
        hessAnalytical <- chessLCMhalfnormlike3C(parm = x$mlParam,
          nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
          nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
          S = x$S, Zvar = Zvar, nZHvar = ncol(Zvar))
      } else {
        if (x$nClasses == 4) {
          hessAnalytical <- chessLCMhalfnormlike4C(parm = x$mlParam,
          nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
          nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
          S = x$S, Zvar = Zvar, nZHvar = ncol(Zvar))
        } else {
          if (x$nClasses == 5) {
          hessAnalytical <- chessLCMhalfnormlike5C(parm = x$mlParam,
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
            S = x$S, Zvar = Zvar, nZHvar = ncol(Zvar))
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
#' @rdname sfalcmcross
#' @export
estfun.sfalcmcross <- function(x, ...) {
  return(x$gradL_OBS)
}
