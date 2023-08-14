################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Information Criteria extraction                                              #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Generalized Zero Inefficiency Stochastic Frontier Analysis        #
#           -Zero inefficiency Stochastic Frontier                             #
#           -Contaminated noise Stochastic Frontier                            #
#           -Multi-Modal Inefficiency Stochastic Frontier Analysis             #
#           -Stochastic/Deterministic Metafrontier Analysis                    #
#           -Sample selection correction for Stochastic Frontier Model         #
#         + Panel data                                                         #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
# Data: Cross sectional data & Pooled data & Panel data                        #
#------------------------------------------------------------------------------#

#' Extract information criteria of stochastic frontier models
#'
#' \code{\link{ic}} returns information criterion from stochastic 
#' frontier models estimated with \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, 
#' \code{\link{sfaselectioncross}}, \code{\link{sfapanel1}}, or 
#' \code{\link{sfalcmpanel}}.
#'
#' The different information criteria are computed as follows: \itemize{ \item
#' AIC: \eqn{-2 \log{LL} + 2 * K} \item BIC: \eqn{-2 \log{LL} + \log{N} * K}
#' \item HQIC: \eqn{-2 \log{LL} + 2 \log{\left[\log{N}\right]} * K} } where
#' \eqn{LL} is the maximum likelihood value, \eqn{K} the number of parameters
#' estimated and \eqn{N} the number of observations.
#'
#' @name ic
#' @aliases ic.sfacross ic.sfalcmcross ic.sfagzisfcross ic.sfacnsfcross
#' ic.sfamisfcross ic.sfazisfcross ic.sfaselectioncross ic.sfametacross
#' ic.sfapanel1 ic.sfalcmpanel
#'
#' @param object A stochastic frontier model returned by \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, 
#' \code{\link{sfaselectioncross}}, \code{\link{sfapanel1}}, or 
#' \code{\link{sfalcmpanel}}.
#' @param IC Character string. Information criterion measure. Three criteria
#' are available: \itemize{ \item \code{'AIC'} for Akaike information criterion
#' (default) \item \code{'BIC'} for Bayesian information criterion \item
#' \code{'HQIC'} for Hannan-Quinn information criterion }.
#' @param ... Currently ignored.
#'
#' @return \code{\link{ic}} returns the value of the information criterion
#' (AIC, BIC or HQIC) of the maximum likelihood coefficients.
#'
# @author K Hervé Dakpo
#'
#' @seealso \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfalcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfagzisfcross}}, for the generalized zero inefficiency
#'  stochastic frontier analysis model fitting function using cross-sectional or 
#'  pooled data.
#' 
#' \code{\link{sfacnsfcross}}, for the contaminated noise stochastic 
#' frontier analysis model fitting function using cross-sectional data.
#' 
#' \code{\link{sfamisfcross}}, for the multi-modal inefficiency stochastic 
#' frontier analysis model fitting function using cross-sectional data.
#' 
#' \code{\link{sfazisfcross}} for zero inefficiency in stochastic frontier model
#' fitting function using cross-sectional data.
#' 
#' \code{\link{sfametacross}}, for fitting different metafrontier models
#' using cross-sectional or pooled data.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfapanel1}}, for the first generation stochastic frontier 
#' analysis model fitting function using panel data.
#' 
#' \code{\link{sfalcmpanel}}, for the latent class stochastic frontier analysis
#' model fitting function using panel data.
#'
#' @keywords methods AIC BIC HQIC
#'
#' @examples
#'
#' \dontrun{
#' # Using data on Swiss railway
#' ## LCM (cost function) half normal distribution
#' cb_2c_u <- sfalcmcross(formula = LNCT ~ LNQ2 + LNQ3 + LNNET + LNPK + LNPL,
#' udist = 'hnormal', uhet = ~ 1, data = swissrailways, S = -1, method='ucminf')
#' ic(cb_2c_u)
#' ic(cb_2c_u, IC = 'BIC')
#' ic(cb_2c_u, IC = 'HQIC')
#' }
#'
#' @export
#' @export ic
# @exportS3Method ic sfacross
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

# information criteria for sfalcmcross ----------
#' @rdname ic
#' @export
# @exportS3Method ic sfalcmcross
ic.sfalcmcross <- function(...) {
  ic.sfacross(...)
}

# information criteria for sfagzisfcross ----------
#' @rdname ic
#' @export
# @exportS3Method ic sfagzisfcross
ic.sfagzisfcross <- function(...) {
  ic.sfacross(...)
}

# information criteria for sfacnsfcross ----------
#' @rdname ic
#' @export
# @exportS3Method ic sfacnsfcross
ic.sfacnsfcross <- function(...) {
  ic.sfacross(...)
}

# information criteria for sfamisfcross ----------
#' @rdname ic
#' @export
# @exportS3Method ic sfamisfcross
ic.sfamisfcross <- function(...) {
  ic.sfacross(...)
}

# information criteria for sfazisfcross ----------
#' @rdname ic
#' @export
# @exportS3Method ic sfazisfcross
ic.sfazisfcross <- function(...) {
  ic.sfacross(...)
}

# information criteria for sfaselectioncross ----------
#' @rdname ic
#' @export
# @exportS3Method ic sfaselectioncross
ic.sfaselectioncross <- function(...) {
  ic.sfacross(...)
}

# information criteria for sfametacross ----------
#' @rdname ic
#' @export
# @exportS3Method ic sfametacross
ic.sfametacross <- function(object, IC = "AIC", ...) {
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
  obj <- rbind(obj)
  row.names(obj) <- IC
  message(paste0(captureIC(obj), collapse = "\n"))
}

# information criteria for sfapanel1 ----------
#' @rdname ic
#' @aliases ic.sfapanel1
#' @export
ic.sfapanel1 <- function(...) {
  ic.sfacross(...)
}

# information criteria for sfalcmpanel ----------
#' @rdname ic
#' @aliases ic.sfalcmpanel
#' @export
ic.sfalcmpanel <- function(...) {
  ic.sfacross(...)
}
