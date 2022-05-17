################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Variance - Covariance Matrix of estimates                                    #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
#         -Sample selection correction                                         #
#         -Zero inefficiency stochastic frontier                               #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract Total Number of Observations Used in frontier models
#' 
#' This function extracts the total number of 'observations' from a
#' fitted frontier model.
#' 
#' `nobs` gives the number of observations actually
#' used by the estimation procedure. It is not necessarily the number
#' of observations of the model frame (number of rows in the model
#' frame), because sometimes the model frame is further reduced by the
#' estimation procedure especially in the presence of NA. In the case of 
#' `sfaselectioncross`, `nobs` returns the number of observations used in the 
#' frontier equation.
#'
#' @name nobs
#' 
#' @param object a `sfacross`, `lcmcross`, `sfaselectioncross` or `zisfcross` 
#' object for which the number of total observations is to be extracted,
#' @param \dots further arguments.
#' 
#' @return A single number, normally an integer.
#' 
# @author K Hervé Dakpo
#' 
#' @keywords attribute
#' 
#' @examples
#' 
#' ## Using data on fossil fuel fired steam electric power generation plants in U.S.
#' # Translog (cost function) half normal with heteroscedasticity
#' tl_u_h <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'hnormal', uhet = ~ regu, data = utility, S = -1, method = 'bfgs')
#'     
#'  nobs(tl_u_h)
#' 
#' @aliases nobs.sfacross
#' @export
# Extract number of observations for sfacros ----------
nobs.sfacross <- function(object, ...) {
  return(object$Nobs)
}

# Extract number of observations for lcmcros ----------
#' @rdname nobs
#' @aliases nobs.lcmcross
#' @export
nobs.lcmcross <- function(object, ...) {
  return(object$Nobs)
}

# Extract number of observations for sfaselectioncros ----------
#' @rdname nobs
#' @aliases nobs.sfaselectioncross
#' @export
nobs.sfaselectioncross <- function(object, ...) {
  return(object$Nobs)
}

# Extract number of observations for zisfcros ----------
#' @rdname nobs
#' @aliases nobs.zisfcross
#' @export
nobs.zisfcross <- function(object, ...) {
  return(object$Nobs)
}
