################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Log-likelihood extraction                                                    #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Sample selection correction for Stochastic Frontier Model         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract log-likelihood value of stochastic frontier models
#'
#' \code{\link{logLik}} extracts the log-likelihood value(s) from stochastic 
#' frontier models estimated with \code{\link{sfacross}}, \code{\link{sfalcmcross}}, 
#' or \code{\link{sfaselectioncross}}.
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{sfacross}}, \code{\link{sfalcmcross}}, or  
#' \code{\link{sfaselectioncross}}.
#' @param individual Logical. If \code{FALSE} (default), the sum of all
#' observations' log-likelihood values is returned. If \code{TRUE}, a vector of
#' each observation's log-likelihood value is returned.
#' @param ... Currently ignored.
#'
#' @name logLik
#'
#' @return \code{\link{logLik}} returns either an object of class 
#' \code{'logLik'}, which is the log-likelihood value with the total number of 
#' observations (\code{nobs}) and the number of free parameters (\code{df}) as 
#' attributes, when \code{individual = FALSE}, or a list of elements, containing 
#' the log-likelihood of each observation (\code{logLik}), the total number of 
#' observations (\code{Nobs}) and the number of free parameters (\code{df}), 
#' when \code{individual = TRUE}.
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
#' @keywords methods likelihood
#'
#' @examples
#'
#' \dontrun{
#' ## Using data on fossil fuel fired steam electric power generation plants in the U.S.
#' # Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = 'mla')
#' logLik(tl_u_ts)
#'
#' ## Using data on eighty-two countries production (GDP)
#' # LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- sfalcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', 
#' data = worldprod, S = 1)
#' logLik(cb_2c_h, individual = TRUE)
#' }
#'
#' @aliases logLik.sfacross
#' @export
#' @export logLik
# log likelihood extraction for sfacross ----------
logLik.sfacross <- function(object, individual = FALSE, ...) {
  if (length(individual) != 1 || !is.logical(individual[1]))
    stop("argument 'individual' must be a single logical value",
      call. = FALSE)
  if (individual) {
    LL <- list()
    LL[["logLik"]] <- object$dataTable$logL_OBS
    LL[["Nobs"]] <- object$Nobs
    LL[["df"]] <- object$nParm
    return(LL)
  } else {
    LL <- object$mlLoglik
    attributes(LL)$nobs <- object$Nobs
    attributes(LL)$df <- object$nParm
    class(LL) <- "logLik"
    return(LL)
  }
}

# log likelihood (LL) extraction for sfalcmcross ----------
#' @rdname logLik
#' @aliases logLik.sfalcmcross
#' @export
logLik.sfalcmcross <- function(object, individual = FALSE, ...) {
  if (length(individual) != 1 || !is.logical(individual[1]))
    stop("argument 'individual' must be a single logical value",
      call. = FALSE)
  if (individual) {
    LL <- list()
    LL[["logLik"]] <- object$dataTable$logL_OBS
    LL[["Nobs"]] <- object$Nobs
    LL[["df"]] <- object$nParm
    return(LL)
  } else {
    LL <- object$mlLoglik
    attributes(LL)$nobs <- object$Nobs
    attributes(LL)$df <- object$nParm
    class(LL) <- "logLik"
    return(LL)
  }
}

# LL extraction for sfaselectioncross ----------
#' @rdname logLik
#' @aliases logLik.sfaselectioncross
#' @export
logLik.sfaselectioncross <- function(object, individual = FALSE,
  ...) {
  if (length(individual) != 1 || !is.logical(individual[1]))
    stop("argument 'individual' must be a single logical value",
      call. = FALSE)
  if (individual) {
    LL <- list()
    LL[["logLik"]] <- object$dataTable$logL_OBS
    LL[["Nobs"]] <- object$Nobs
    LL[["df"]] <- object$nParm
    return(LL)
  } else {
    LL <- object$mlLoglik
    attributes(LL)$nobs <- object$Nobs
    attributes(LL)$df <- object$nParm
    class(LL) <- "logLik"
    return(LL)
  }
}
