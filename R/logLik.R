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
#' Extract log-likelihood value of stochastic frontier models
#'
#' \code{\link{logLik}} extracts the log-likelihood value(s) from stochastic 
#' frontier models estimated with \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, 
#' \code{\link{sfaselectioncross}}, \code{\link{sfapanel1}}, or 
#' \code{\link{sfalcmpanel}}.
#'
#' @param object A stochastic frontier model returned by \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, 
#' \code{\link{sfaselectioncross}}, \code{\link{sfapanel1}}, or 
#' \code{\link{sfalcmpanel}}.
#' @param individual Logical. If \code{FALSE} (default), the sum of all
#' observations' log-likelihood values is returned. If \code{TRUE}, a vector of
#' each observation's log-likelihood value is returned.
#' @param ... Currently ignored.
#'
#' @name logLik
#' 
#' @aliases logLik.sfacross logLik.sfalcmcross logLik.sfagzisfcross
#' logLik.sfacnsfcross logLik.sfamisfcross logLik.sfazisfcross
#' logLik.sfaselectioncross logLik.sfametacross logLik.sfapanel1 
#' logLik.sfalcmpanel
#'
#' @return \code{\link{logLik}} returns either an object of class 
#' \code{'logLik'}, which is the log-likelihood value with the total number of 
#' observations (\code{nobs}) and the number of free parameters (\code{df}) as 
#' attributes, when \code{individual = FALSE}, or a list of elements, containing 
#' the log-likelihood of each observation (\code{logLik}), the total number of 
#' observations (\code{Nobs}) and the number of free parameters (\code{df}), 
#' when \code{individual = TRUE}.
#' 
#' For object of class \code{'sfametacross'}, the groups log-likelihood values 
#' are returned, and in the case of \code{'modelType'} = 'hhl14', the 
#' metafrontier log-likelihood value is also returned.
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
#' @keywords methods logLik
#'
#' @examples
#'
#' \dontrun{
#' # Using data on fossil fuel fired steam electric power generation plants in 
#' # the U.S.
#' ## Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = 'mla')
#' logLik(tl_u_ts)
#'
#' # Using data on eighty-two countries production (GDP)
#' ## LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- sfalcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', 
#' data = worldprod, S = 1)
#' logLik(cb_2c_h, individual = TRUE)
#' }
#'
#' @export
# @exportS3Method logLik sfacross
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

# log likelihood extraction for sfalcmcross ----------
#' @rdname logLik
#' @export
# @exportS3Method logLik sfalcmcross
logLik.sfalcmcross <- function(...) {
  logLik.sfacross(...)
}

# log likelihood extraction for sfagzisfcross ----------
#' @rdname logLik
#' @export
# @exportS3Method logLik sfagzisfcross
logLik.sfagzisfcross <- function(...) {
  logLik.sfacross(...)
}

# log likelihood extraction for sfacnsfcross ----------
#' @rdname logLik
#' @export
# @exportS3Method logLik sfacnsfcross
logLik.sfacnsfcross <- function(...) {
  logLik.sfacross(...)
}

# log likelihood extraction for sfamisfcross ----------
#' @rdname logLik
#' @export
# @exportS3Method logLik sfamisfcross
logLik.sfamisfcross <- function(...) {
  logLik.sfacross(...)
}

# log likelihood extraction for sfazisfcross ----------
#' @rdname logLik
#' @export
# @exportS3Method logLik sfazisfcross
logLik.sfazisfcross <- function(...) {
  logLik.sfacross(...)
}

# LL extraction for sfaselectioncross ----------
#' @rdname logLik
#' @export
# @exportS3Method logLik sfaselectioncross
logLik.sfaselectioncross <- function(...) {
  logLik.sfacross(...)
}

# log likelihood extraction for sfametacross ----------
#' @rdname logLik
#' @export
# @exportS3Method logLik sfametacross
logLik.sfametacross <- function(...) {
  logLik.sfacross(...)
}

# log likelihood extraction for sfapanel1 ----------
#' @rdname logLik
#' @aliases logLik.sfapanel1
#' @export
logLik.sfapanel1 <- function(...) {
  logLik.sfacross(...)
}

# log likelihood extraction for sfalcmpanel ----------
#' @rdname logLik
#' @aliases logLik.sfalcmpanel
#' @export
logLik.sfalcmpanel <- function(...) {
  logLik.sfacross(...)
}
