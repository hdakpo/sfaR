################################################################################
#                                                                              #
# Deprecated functions for the sfaR package                                    #
#                                                                              #
################################################################################

#' Deprecated functions of sfaR
#' 
#' @description
#' These functions are provided for compatibility with older versions of 
#' \sQuote{sfaR} only, and could be defunct at a future release.
#' 
#' @name sfaR-deprecated
# @aliases sfaR-deprecated
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
#' cross section SFA model. \code{'whichStart = 1'} can be fast.
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
#' @param extraPar Logical (default = \code{FALSE}). If \code{TRUE}, additional
#' parameters are returned (see \code{\link{coef}} or \code{\link{vcov}}).
#' @param ci Logical. Default = \code{FALSE}. If \code{TRUE}, the 95%
#' confidence interval for the different parameters (OLS or/and ML estimates) is
#' returned.
#' @param digits Numeric. Number of digits displayed in values.
#' @param grad Logical. Default = \code{FALSE}. If \code{TRUE}, the gradient
#' for the maximum likelihood (ML) estimates of the different parameters is
#' returned.
#' @param IC Character string. Information criterion measure. Three criteria
#' are available: \itemize{ \item \code{'AIC'} for Akaike information criterion
#' (default) \item \code{'BIC'} for Bayesian information criterion \item
#' \code{'HQIC'} for Hannan-Quinn information criterion }.
#' @param individual Logical. If \code{FALSE} (default), the sum of all
#' observations' log-likelihood values is returned. If \code{TRUE}, a vector of
#' each observation's log-likelihood value is returned.
#' @param level A number between between 0 and 0.9999 used for the computation
#' of (in-)efficiency confidence intervals (defaut = \code{0.95}). Not used in the
#' case of \code{lcmcross}.
#' @param object an object of class lcmcross (returned by the function
#' \code{\link{lcmcross}}).
#' @param newData Optional data frame that is used to calculate the efficiency 
#' estimates. If NULL (the default), the efficiency estimates are calculated 
#' for the observations that were used in the estimation.
#' @param x an object of class lcmcross (returned by the function
#' \code{\link{lcmcross}}).
#' @param ... additional arguments of frontier are passed to lcmcross; 
#' additional arguments of the print, bread, estfun, nobs methods are currently 
#' ignored.
#' 
#' @details
#'  The following functions are deprecated and could be removed from \pkg{sfaR} 
#'  in a near future. Use the replacement indicated below:
#'  \itemize{
#'  \item{lcmcross: \code{\link{sfalcmcross}}}
#'  \item{bread.lcmcross: \code{\link{bread.sfalcmcross}}}
#'  \item{coef.lcmcross: \code{\link{coef.sfalcmcross}}}
#'  \item{coef.summary.lcmcross: \code{\link{coef.summary.sfalcmcross}}}
#'  \item{efficiencies.lcmcross: \code{\link{efficiencies.sfalcmcross}}}
#'  \item{estfun.lcmcross: \code{\link{estfun.sfalcmcross}}}
#'  \item{fitted.lcmcross: \code{\link{fitted.sfalcmcross}}}
#'  \item{ic.lcmcross: \code{\link{ic.sfalcmcross}}}
#'  \item{logLik.lcmcross: \code{\link{logLik.sfalcmcross}}}
#'  \item{marginal.lcmcross: \code{\link{marginal.sfalcmcross}}}
#'  \item{nobs.lcmcross: \code{\link{nobs.sfalcmcross}}}
#'  \item{print.lcmcross: \code{\link{print.sfalcmcross}}}
#'  \item{print.summary.lcmcross: \code{\link{print.summary.sfalcmcross}}}
#'  \item{residuals.lcmcross: \code{\link{residuals.sfalcmcross}}}
#'  \item{summary.lcmcross: \code{\link{summary.sfalcmcross}}}
#'  \item{vcov.lcmcross: \code{\link{vcov.sfalcmcross}}}
#'  }
#'
NULL

# lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
lcmcross <- function(formula, uhet, vhet, thet, logDepVar = TRUE,
  data, subset, weights, wscale = TRUE, S = 1L, udist = "hnormal",
  start = NULL, whichStart = 2L, initAlg = "nm", initIter = 100,
  lcmClasses = 2, method = "bfgs", hessianType = 1, itermax = 2000L,
  printInfo = FALSE, tol = 1e-12, gradtol = 1e-06, stepmax = 0.1,
  qac = "marquardt") {
  .Deprecated(new = "sfalcmcross", package = "sfaR", msg = " `lcmcross` is deprecated, use 'sfalcmcross' instead for same functionality",
    old = "lcmcross")
  if (missing(subset)) {
    if (missing(weights)) {
      x <- sfalcmcross(formula = formula, uhet = uhet,
        vhet = vhet, thet = thet, logDepVar = logDepVar,
        data = data, wscale = wscale, S = S, udist = udist,
        start = start, whichStart = whichStart, initAlg = initAlg,
        initIter = initIter, lcmClasses = lcmClasses,
        method = method, hessianType = hessianType, itermax = itermax,
        printInfo = printInfo, tol = tol, gradtol = gradtol,
        stepmax = stepmax, qac = qac)
    } else {
      x <- sfalcmcross(formula = formula, uhet = uhet,
        vhet = vhet, thet = thet, logDepVar = logDepVar,
        data = data, weights = weights, wscale = wscale,
        S = S, udist = udist, start = start, whichStart = whichStart,
        initAlg = initAlg, initIter = initIter, lcmClasses = lcmClasses,
        method = method, hessianType = hessianType, itermax = itermax,
        printInfo = printInfo, tol = tol, gradtol = gradtol,
        stepmax = stepmax, qac = qac)
    }
  } else {
    if (missing(weights)) {
      x <- sfalcmcross(formula = formula, uhet = uhet,
        vhet = vhet, thet = thet, logDepVar = logDepVar,
        data = data, subset = subset, wscale = wscale,
        S = S, udist = udist, start = start, whichStart = whichStart,
        initAlg = initAlg, initIter = initIter, lcmClasses = lcmClasses,
        method = method, hessianType = hessianType, itermax = itermax,
        printInfo = printInfo, tol = tol, gradtol = gradtol,
        stepmax = stepmax, qac = qac)
    } else {
      x <- sfalcmcross(formula = formula, uhet = uhet,
        vhet = vhet, thet = thet, logDepVar = logDepVar,
        data = data, subset = subset, weights = weights,
        wscale = wscale, S = S, udist = udist, start = start,
        whichStart = whichStart, initAlg = initAlg, initIter = initIter,
        lcmClasses = lcmClasses, method = method, hessianType = hessianType,
        itermax = itermax, printInfo = printInfo, tol = tol,
        gradtol = gradtol, stepmax = stepmax, qac = qac)
    }
  }
  class(x) <- c("lcmcross", "sfalcmcross")
  return(x)
}

# print.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
print.lcmcross <- function(x, ...) {
  .Deprecated(new = "print.sfalcmcross", package = "sfaR",
    msg = " `print.lcmcross` is deprecated, use 'print.sfalcmcross' instead for same functionality",
    old = "print.lcmcross")
  print.sfalcmcross(x, ...)
}

# bread.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
bread.lcmcross <- function(x, ...) {
  .Deprecated(new = "bread.sfalcmcross", package = "sfaR",
    msg = " `bread.lcmcross` is deprecated, use 'bread.sfalcmcross' instead for same functionality",
    old = "bread.lcmcross")
  bread.sfalcmcross(x = x, ...)
}

# estfun.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
estfun.lcmcross <- function(x, ...) {
  .Deprecated(new = "estfun.sfalcmcross", package = "sfaR",
    msg = " `estfun.lcmcross` is deprecated, use 'estfun.sfalcmcross' instead for same functionality",
    old = "estfun.lcmcross")
  estfun.sfalcmcross(x = x, ...)
}

# coef.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
coef.lcmcross <- function(object, extraPar = FALSE, ...) {
  .Deprecated(new = "coef.sfalcmcross", package = "sfaR", msg = " `coef.lcmcross` is deprecated, use 'coef.sfalcmcross' instead for same functionality",
    old = "coef.lcmcross")
  coef.sfalcmcross(object = object, extraPar = extraPar, ...)
}

# coef.summary.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
coef.summary.lcmcross <- function(object, ...) {
  .Deprecated(new = "coef.summary.sfalcmcross", package = "sfaR",
    msg = " `coef.summary.lcmcross` is deprecated, use 'coef.summary.sfalcmcross' instead for same functionality",
    old = "coef.summary.lcmcross")
  coef.summary.sfalcmcross(object = object, ...)
}

# fitted.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
fitted.lcmcross <- function(object, ...) {
  .Deprecated(new = "fitted.sfalcmcross", package = "sfaR",
    msg = " `fitted.lcmcross` is deprecated, use 'fitted.sfalcmcross' instead for same functionality",
    old = "fitted.lcmcross")
  fitted.sfalcmcross(object = object, ...)
}

# ic.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
ic.lcmcross <- function(object, IC = "AIC", ...) {
  .Deprecated(new = "ic.sfalcmcross", package = "sfaR", msg = " `ic.lcmcross` is deprecated, use 'ic.sfalcmcross' instead for same functionality",
    old = "ic.lcmcross")
  ic.sfalcmcross(object = object, IC = IC, ...)
}

# logLik.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
logLik.lcmcross <- function(object, individual = FALSE, ...) {
  .Deprecated(new = "logLik.sfalcmcross", package = "sfaR",
    msg = " `logLik.lcmcross` is deprecated, use 'logLik.sfalcmcross' instead for same functionality",
    old = "logLik.lcmcross")
  logLik.sfalcmcross(object = object, individual = individual,
    ...)
}

# marginal.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
marginal.lcmcross <- function(object, newData = NULL, ...) {
  .Deprecated(new = "marginal.sfalcmcross", package = "sfaR",
    msg = " `marginal.lcmcross` is deprecated, use 'marginal.sfalcmcross' instead for same functionality",
    old = "marginal.lcmcross")
  marginal.sfalcmcross(object = object, newData = newData)
}

# nobs.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
nobs.lcmcross <- function(object, ...) {
  .Deprecated(new = "nobs.sfalcmcross", package = "sfaR", msg = " `nobs.lcmcross` is deprecated, use 'nobs.sfalcmcross' instead for same functionality",
    old = "nobs.lcmcross")
  nobs.sfalcmcross(object = object, ...)
}

# residuals from sfalcmcross ----------
#' @rdname sfaR-deprecated
#' @export
residuals.lcmcross <- function(object, ...) {
  .Deprecated(new = "residuals.sfalcmcross", package = "sfaR",
    msg = " `residuals.lcmcross` is deprecated, use 'residuals.sfalcmcross' instead for same functionality",
    old = "residuals.lcmcross")
  residuals.sfalcmcross(object = object, ...)
}

# summary.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
summary.lcmcross <- function(object, grad = FALSE, ci = FALSE,
  ...) {
  .Deprecated(new = "summary.sfalcmcross", package = "sfaR",
    msg = " `summary.lcmcross` is deprecated, use 'summary.sfalcmcross' instead for same functionality",
    old = "summary.lcmcross")
  # class here
  su <- summary.sfalcmcross(object = object, grad = grad, ci = ci,
    ...)
  class(su) <- c("summary.lcmcross", "summary.sfalcmcross")
  su
}

# print.summary.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
print.summary.lcmcross <- function(x, digits = max(3, getOption("digits") -
  2), ...) {
  .Deprecated(new = "print.summary.sfalcmcross", package = "sfaR",
    msg = " `print.summary.lcmcross` is deprecated, use 'print.summary.sfalcmcross' instead for same functionality",
    old = "print.summary.lcmcross")
  print.summary.sfalcmcross(x = x, digits = digits, ...)
}

# efficiencies.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
efficiencies.lcmcross <- function(object, level = 0.95, newData = NULL,
  ...) {
  .Deprecated(new = "efficiencies.sfalcmcross", package = "sfaR",
    msg = " `efficiencies.lcmcross` is deprecated, use 'efficiencies.sfalcmcross' instead for same functionality",
    old = "efficiencies.lcmcross")
  efficiencies.sfalcmcross(object = object, level = level,
    newData = newData, ...)
}

# vcov.lcmcross ----------
#' @rdname sfaR-deprecated
#' @export
vcov.lcmcross <- function(object, ...) {
  .Deprecated(new = "vcov.sfalcmcross", package = "sfaR", msg = " `vcov.lcmcross` is deprecated, use 'vcov.sfalcmcross' instead for same functionality",
    old = "vcov.lcmcross")
  vcov.sfalcmcross(object = object, ...)
}
