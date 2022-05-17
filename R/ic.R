################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Information Criteria extraction                                              #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
#         -Sample selection correction                                         #
#         -Zero inefficiency stochastic frontier                               #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract information criteria of stochastic frontier models
#'
#' \code{\link{ic}} returns information criterion from classic or latent class
#' stochastic frontier models estimated with \code{\link{sfacross}},
#' \code{\link{lcmcross}}, \code{\link{sfaselectioncross}} or 
#' \code{\link{zisfcross}}.
#'
#' The different information criteria are computed as follows: \itemize{ \item
#' AIC: \eqn{-2 \log{LL} + 2 * K} \item BIC: \eqn{-2 \log{LL} + \log{N} * K}
#' \item HQIC: \eqn{-2 \log{LL} + 2 \log{\left[\log{N}\right]} * K} } where
#' \eqn{LL} is the maximum likelihood value, \eqn{K} the number of parameters
#' estimated and \eqn{N} the number of observations.
#'
#' @name ic
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{sfacross}}, \code{\link{lcmcross}}, 
#' \code{\link{sfaselectioncross}} or \code{\link{zisfcross}}.
#' @param IC Character string. Information criterion measure. Three criteria
#' are available: \itemize{ \item \code{'AIC'} for Akaike information criterion
#' (default) \item \code{'BIC'} for Bayesian information criterion \item
#' \code{'HQIC'} for Hannan-Quinn information criterion }
#' @param ... Currently ignored.
#'
#' @return \code{\link{ic}} returns the value of the information criterion
#' (AIC, BIC or HQIC) of the maximum likelihood coefficients.
#'
# @author K Hervé Dakpo
#'
#' @seealso \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function.
#'
#' \code{\link{lcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function.
#' 
#' \code{\link{zisfcross}} for zero inefficiency in stochastic frontier model
#' fitting function.
#'
#' @keywords methods AIC BIC HQIC
#'
#' @examples
#'
#' ## Using data on Swiss railway
#' # LCM (cost function) half normal distribution
#' cb_2c_u <- lcmcross(formula = LNCT ~ LNQ2 + LNQ3 + LNNET + LNPK + LNPL,
#' udist = 'hnormal', uhet = ~ 1, data = swissrailways, S = -1, method='ucminf')
#' ic(cb_2c_u)
#' ic(cb_2c_u, IC = 'BIC')
#' ic(cb_2c_u, IC = 'HQIC')
#'
#' @aliases ic.sfacross
#' @export
#' @export ic
# information criteria for sfacross ----------
ic.sfacross <- function(object, IC = "AIC", ...) {
  if (!(IC %in% c("AIC", "BIC", "HQIC"))) {
    stop("Unknown information criteria: ", paste(IC), call. = FALSE)
  }
  if (IC == "AIC") {
    obj <- -2 * object$mlLoglik + 2 * object$nParm
  } else {
    if (IC == "BIC") {
      obj <- -2 * object$mlLoglik + log(object$Nobs) *
        object$nParm
    } else {
      if (IC == "HQIC") {
        obj <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) *
          object$nParm
      }
    }
  }
  message(IC, ": ", prettyNum(obj), sep = "")
}

# information criteria for lcmcross ----------
#' @rdname ic
#' @aliases ic.lcmcross
#' @export
ic.lcmcross <- function(object, IC = "AIC", ...) {
  if (!(IC %in% c("AIC", "BIC", "HQIC"))) {
    stop("Unknown information criteria: ", paste(IC), call. = FALSE)
  }
  if (IC == "AIC") {
    obj <- -2 * object$mlLoglik + 2 * object$nParm
  } else {
    if (IC == "BIC") {
      obj <- -2 * object$mlLoglik + log(object$Nobs) *
        object$nParm
    } else {
      if (IC == "HQIC") {
        obj <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) *
          object$nParm
      }
    }
  }
  message(IC, ": ", prettyNum(obj), sep = "")
}

# information criteria for sfaselectioncross ----------
#' @rdname ic
#' @aliases ic.sfaselectioncross
#' @export
ic.sfaselectioncross <- function(object, IC = "AIC", ...) {
  if (!(IC %in% c("AIC", "BIC", "HQIC"))) {
    stop("Unknown information criteria: ", paste(IC), call. = FALSE)
  }
  if (IC == "AIC") {
    obj <- -2 * object$mlLoglik + 2 * object$nParm
  } else {
    if (IC == "BIC") {
      obj <- -2 * object$mlLoglik + log(object$Nobs) *
        object$nParm
    } else {
      if (IC == "HQIC") {
        obj <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) *
          object$nParm
      }
    }
  }
  message(IC, ": ", prettyNum(obj), sep = "")
}

# information criteria for zisfcross ----------
#' @rdname ic
#' @aliases ic.zisfcross
#' @export
ic.zisfcross <- function(object, IC = "AIC", ...) {
  if (!(IC %in% c("AIC", "BIC", "HQIC"))) {
    stop("Unknown information criteria: ", paste(IC), call. = FALSE)
  }
  if (IC == "AIC") {
    obj <- -2 * object$mlLoglik + 2 * object$nParm
  } else {
    if (IC == "BIC") {
      obj <- -2 * object$mlLoglik + log(object$Nobs) *
        object$nParm
    } else {
      if (IC == "HQIC") {
        obj <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) *
          object$nParm
      }
    }
  }
  message(IC, ": ", prettyNum(obj), sep = "")
}
