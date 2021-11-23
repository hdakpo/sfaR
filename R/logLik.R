################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Log-likelihood extraction                                                    #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract log-likelihood value of classic or latent class stochastic models
#'
#' \code{\link{logLik}} extracts the log-likelihood value(s) from classic or
#' latent class stochastic frontier models estimated with
#' \code{\link{sfacross}} or \code{\link{lcmcross}}.
#'
#' @param object A classic or latent class stochastic frontier model returned
#' by \code{\link{sfacross}} or \code{\link{lcmcross}}.
#' @param individual Logical. If \code{FALSE} (default), the sum of all
#' observations' log-likelihood values is returned. If \code{TRUE}, a vector of
#' each observation's log-likelihood value is returned.
#' @param ... Currently ignored.
#'
#' @name logLik
#'
#' @return \code{\link{logLik}} returns an object of class \code{'logLik'},
#' which is either a numeric matrix with the log-likelihood value
#' (\code{logLik}), the total number of observations (\code{Nobs}) and the
#' number of free parameters (\code{df}), when \code{individual = FALSE},
#' or a list of elements, containing the log-likelihood of each observation
#' (\code{logLik}), the total number of observations (\code{Nobs}) and the
#' number of free parameters (\code{df}), when \code{individual = TRUE}.
#'
#' @author K Hervé Dakpo, Yann Desjeux and Laure Latruffe
#'
#' @seealso \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function.
#'
#' \code{\link{lcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function.
#'
#' @keywords methods likelihood
#'
#' @examples
#'
#' ## Using data on fossil fuel fired steam electric power generation plants in the U.S.
#' # Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#'     log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#'     I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#'     udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#'     scaling = TRUE, method = 'mla')
#'   logLik(tl_u_ts)
#'
#' ## Using data on eighty-two countries production (DGP)
#' # LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', data = worldprod, S = 1)
#'   logLik(cb_2c_h, individual = TRUE)
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
    cat("'log Lik.'", round(object$mlLoglik, 5), paste0("(df=", object$nParm, ")"))
    # LL <- rbind(`logLik: ` = object$mlLoglik, `Nobs: ` = object$Nobs,
    #   `df: ` = object$nParm)
  }
}

# log likelihood extraction for lcmcross ----------
#' @rdname logLik
#' @aliases logLik.lcmcross
#' @export
logLik.lcmcross <- function(object, individual = FALSE, ...) {
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
    cat("'log Lik.'", round(object$mlLoglik, 5), paste0("(df=", object$nParm, ")"))
    # LL <- rbind(`logLik: ` = object$mlLoglik, `Nobs: ` = object$Nobs,
    #   `df: ` = object$nParm)
  }
}
