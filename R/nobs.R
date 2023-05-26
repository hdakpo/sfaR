################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Variance - Covariance Matrix of estimates                                    #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Sample selection correction for Stochastic Frontier Model         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract total number of observations used in frontier models
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
#' @param object a `sfacross`, `sfalcmcross`, or `sfaselectioncross`
#' object for which the number of total observations is to be extracted.
#' @param \dots Currently ignored.
#' 
#' @return A single number, normally an integer.
#' 
# @author K Hervé Dakpo
#' 
#' @seealso \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfalcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function using cross-sectional or pooled data.
#'
#' @keywords attribute
#' 
#' @examples
#' 
#' \dontrun{
#' ## Using data on fossil fuel fired steam electric power generation plants in the U.S.
#' # Translog (cost function) half normal with heteroscedasticity
#' tl_u_h <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'hnormal', uhet = ~ regu, data = utility, S = -1, method = 'bfgs')
#' nobs(tl_u_h)
#' }
#' 
#' @aliases nobs.sfacross
#' @export
# Extract number of observations for sfacross ----------
nobs.sfacross <- function(object, ...) {
  return(object$Nobs)
}

# Extract number of observations for sfalcmcross ----------
#' @rdname nobs
#' @aliases nobs.sfalcmcross
#' @export
nobs.sfalcmcross <- function(object, ...) {
  return(object$Nobs)
}

# Extract number of observations for sfaselectioncross ----------
#' @rdname nobs
#' @aliases nobs.sfaselectioncross
#' @export
nobs.sfaselectioncross <- function(object, ...) {
  return(object$Nobs)
}
