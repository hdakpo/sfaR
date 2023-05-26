################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Marginal impact of Z variables on inefficiency                               #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Sample selection correction for Stochastic Frontier Model         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Marginal effects of the inefficiency drivers in stochastic frontier models
#'
#' This function returns marginal effects of the inefficiency drivers from stochastic 
#' frontier models estimated with \code{\link{sfacross}}, \code{\link{sfalcmcross}}, 
#' or \code{\link{sfaselectioncross}}.
#'
#' @details \code{\link{marginal}} operates in the presence of exogenous 
#' variables that explain inefficiency, namely the inefficiency drivers
#'  (\eqn{uhet = ~ Z_u} or \eqn{muhet = ~ Z_{mu}}).
#'
#' Two components are computed for each variable: the marginal effects on the
#' expected inefficiency (\eqn{\frac{\partial E[u]}{\partial Z_{mu}}}) and
#' the marginal effects on the variance of inefficiency (\eqn{\frac{\partial
#' V[u]}{\partial Z_{mu}}}).
#'
#' The model also allows the Wang (2002) parametrization of \eqn{\mu} and
#' \eqn{\sigma_u^2} by the same vector of exogenous variables. This double
#' parameterization accounts for non-monotonic relationships between the
#' inefficiency and its drivers.
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{sfacross}}, \code{\link{sfalcmcross}}, or 
#' \code{\link{sfaselectioncross}}.
#' @param newData Optional data frame that is used to calculate the marginal 
#' effect of \eqn{Z} variables on inefficiency. If NULL (the default), the 
#' marginal estimates are calculated for the observations that were used in the 
#' estimation.
#' @param ... Currently ignored.
#'
#' @name marginal
#'
#' @return \code{\link{marginal}} returns a data frame containing the marginal
#' effects of the \eqn{Z_u} variables on the expected inefficiency (each
#' variable has the prefix \code{'Eu_'}) and on the variance of the
#' inefficiency (each variable has the prefix \code{'Vu_'}).
#'
#' In the case of the latent class stochastic frontier (LCM), each variable 
#' ends with \code{'_c#'} where \code{'#'} is the class number.
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
#' @references Wang, H.J. 2002. Heteroscedasticity and non-monotonic efficiency
#' effects of a stochastic frontier model. \emph{Journal of Productivity
#' Analysis}, \bold{18}:241--253.
#'
#' @keywords methods marginal
#'
#' @examples
#'
#' \dontrun{
#' ## Using data on fossil fuel fired steam electric power generation plants in the U.S.
#' # Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu + wl, uhet = ~ regu + wl, data = utility, 
#' S = -1, scaling = TRUE, method = 'mla')
#' marg.tl_u_ts <- marginal(tl_u_ts)
#' summary(marg.tl_u_ts)
#'
#' ## Using data on eighty-two countries production (GDP)
#' # LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- sfalcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal',
#'     data = worldprod, uhet = ~ initStat + h, S = 1, method = 'mla')
#'   marg.cb_2c_h <- marginal(cb_2c_h)
#'   summary(marg.cb_2c_h)
#'   }
#'
#' @aliases marginal.sfacross
#' @export
#' @export marginal
# marginal effects computation sfacross ----------
marginal.sfacross <- function(object, newData = NULL, ...) {
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$udist == "hnormal") {
    if (object$nuZUvar == 1) {
      stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
        call. = FALSE)
    } else {
      EffMarg <- data.frame(cbind(cmarghalfnorm_Eu(object = object),
        cmarghalfnorm_Vu(object = object)))
    }
  } else {
    if (object$udist == "exponential") {
      if (object$nuZUvar == 1) {
        stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
      } else {
        EffMarg <- data.frame(cbind(cmargexponorm_Eu(object = object),
          cmargexponorm_Vu(object = object)))
      }
    } else {
      if (object$udist == "gamma") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- data.frame(cbind(cmarggammanorm_Eu(object = object),
          cmarggammanorm_Vu(object = object)))
        }
      } else {
        if (object$udist == "rayleigh") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- data.frame(cbind(cmargraynorm_Eu(object = object),
            cmargraynorm_Vu(object = object)))
          }
        } else {
          if (object$udist == "uniform") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- data.frame(cbind(cmarguninorm_Eu(object = object),
            cmarguninorm_Vu(object = object)))
          }
          } else {
          if (object$udist == "tnormal") {
            if (object$scaling) {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- data.frame(cbind(cmargtruncnormscal_Eu(object = object),
              cmargtruncnormscal_Vu(object = object)))
            }
            } else {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- data.frame(cbind(cmargtruncnorm_Eu(object = object),
              cmargtruncnorm_Vu(object = object)))
            }
            }
          } else {
            if (object$udist == "lognormal") {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- data.frame(cbind(cmarglognorm_Eu(object = object),
              cmarglognorm_Vu(object = object)))
            }
            } else {
            if (object$udist == "genexponential") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- data.frame(cbind(cmarggenexponorm_Eu(object = object),
                cmarguninorm_Vu(object = object)))
              }
            } else {
              if (object$udist == "tslaplace") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- data.frame(cbind(cmargtslnorm_Eu(object = object),
                cmargtslnorm_Vu(object = object)))
              }
              } else {
              if (object$udist == "weibull") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- data.frame(cbind(cmargweibullnorm_Eu(object = object),
                  cmargweibullnorm_Vu(object = object)))
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  return(EffMarg)
}

# marginal effects computation sfalcmcross ----------
#' @rdname marginal
#' @aliases marginal.sfalcmcross
#' @export
marginal.sfalcmcross <- function(object, newData = NULL, ...) {
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$nuZUvar == 1) {
    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
      call. = FALSE)
  } else {
    if (object$nClasses == 2) {
      EffMarg <- data.frame(cbind(cmargLCM2Chalfnorm_Eu(object = object),
        cmargLCM2Chalfnorm_Vu(object = object)))
    } else {
      if (object$nClasses == 3) {
        EffMarg <- data.frame(cbind(cmargLCM3Chalfnorm_Eu(object = object),
          cmargLCM3Chalfnorm_Vu(object = object)))
      } else {
        if (object$nClasses == 4) {
          EffMarg <- data.frame(cbind(cmargLCM4Chalfnorm_Eu(object = object),
          cmargLCM4Chalfnorm_Vu(object = object)))
        } else {
          if (object$nClasses == 5) {
          EffMarg <- data.frame(cbind(cmargLCM5Chalfnorm_Eu(object = object),
            cmargLCM5Chalfnorm_Vu(object = object)))
          }
        }
      }
    }
  }
  return(EffMarg)
}

# marginal effects computation sfaselectioncross ----------
#' @rdname marginal
#' @aliases marginal.sfaselectioncross
#' @export
marginal.sfaselectioncross <- function(object, newData = NULL,
  ...) {
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$nuZUvar == 1) {
    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
      call. = FALSE)
  } else {
    EffMarg <- data.frame(cbind(cmarghalfnorm_Eu_ss(object = object),
      cmarghalfnorm_Vu_ss(object = object)))
  }
  return(EffMarg)
}
