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

#' Extract total number of observations used in frontier models
#' 
#' This function extracts the total number of 'observations' from a
#' fitted frontier model.
#' 
#' `nobs` gives the number of observations actually
#' used by the estimation procedure. It is not necessarily the number
#' of observations of the model frame (number of rows in the model
#' frame), because sometimes the model frame is further reduced by the
#' estimation procedure especially in the presence of 'NA' or 'Inf'. In the case 
#' of `sfaselectioncross`, `nobs` returns the number of observations used in the 
#' frontier equation.
#'
#' @name nobs
#' 
#' @aliases nobs.sfacross nobs.sfalcmcross nobs.sfagzisfcross nobs.sfacnsfcross 
#' nobs.sfamisfcross nobs.sfazisfcross nobs.sfaselectioncross nobs.sfametacross
#' nobs.sfapanel1 nobs.sfalcmpanel
#' 
#' @param object a `sfacross`, `sfalcmcross`, `sfagzisfcross`, `sfacnsfcross`, 
#' `sfamisfcross`, `sfazisfcross`, `sfaselectioncross`, `sfapanel1`, or 
#' `sfalcmpanel` object for which the number of total observations is to be 
#' extracted. For object `sfametacross` the number of observations per group is 
#' also returned.
#' @param \dots Currently ignored.
#' 
#' @return A single number, normally an integer, except for the metafrontier 
#' case where of vector of integers is returned.
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
#' @keywords methods nobs
#' 
#' @examples
#' 
#' \dontrun{
#' # Using data on fossil fuel fired steam electric power generation plants in 
#' # the U.S.
#' ## Translog (cost function) half normal with heteroscedasticity
#' tl_u_h <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'hnormal', uhet = ~ regu, data = utility, S = -1, method = 'bfgs')
#' nobs(tl_u_h)
#' }
#' 
#' @export
# @exportS3Method nobs sfacross
# Extract number of observations for sfacross ----------
nobs.sfacross <- function(object, ...) {
  return(object$Nobs)
}

# Extract number of observations for sfalcmcross ----------
#' @rdname nobs
#' @export
# @exportS3Method nobs sfalcmcross
nobs.sfalcmcross <- function(...) {
  nobs.sfacross(...)
}

# Extract number of observations for sfagzisfcross ----------
#' @rdname nobs
#' @export
# @exportS3Method nobs sfagzisfcross
nobs.sfagzisfcross <- function(...) {
  nobs.sfacross(...)
}

# Extract number of observations for sfacnsfcross ----------
#' @rdname nobs
#' @export
# @exportS3Method nobs sfacnsfcross
nobs.sfacnsfcross <- function(...) {
  nobs.sfacross(...)
}

# Extract number of observations for sfamisfcross ----------
#' @rdname nobs
#' @export
# @exportS3Method nobs sfamisfcross
nobs.sfamisfcross <- function(...) {
  nobs.sfacross(...)
}

# Extract number of observations for sfazisfcross ----------
#' @rdname nobs
#' @export
# @exportS3Method nobs sfazisfcross
nobs.sfazisfcross <- function(...) {
  nobs.sfacross(...)
}

# Extract number of observations for sfaselectioncross ----------
#' @rdname nobs
#' @export
# @exportS3Method nobs sfaselectioncross
nobs.sfaselectioncross <- function(...) {
  nobs.sfacross(...)
}

# Extract number of observations for sfametacross ----------
#' @rdname nobs
#' @export
# @exportS3Method nobs sfametacross
nobs.sfametacross <- function(...) {
  nobs.sfacross(...)
}

# Extract number of observations for sfapanel1 ----------
#' @rdname nobs
#' @aliases nobs.sfapanel1
#' @export
nobs.sfapanel1 <- function(object, ...) {
  return(list(Observations = object$Nobs, `Cross-sections` = object$Nid))
}

# Extract number of observations for sfalcmpanel ----------
#' @rdname nobs
#' @aliases nobs.sfalcmpanel
#' @export
nobs.sfalcmpanel <- function(object, ...) {
  return(list(Observations = object$Nobs, `Cross-sections` = object$Nid))
}
