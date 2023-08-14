################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Coefficients extraction                                                      #
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

#' Extract coefficients of stochastic frontier models
#' 
#' @description
#' From an object of class \code{'summary.sfacross'},
#' \code{'summary.sfalcmcross'}, \code{'summary.sfagzisfcross'}, 
#' \code{'summary.sfacnsfcross'}, \code{'summary.sfamisfcross'}, 
#' \code{'summary.sfazisfcross'}, \code{'summary.sfametacross'},
#' \code{'summary.sfaselectioncross'}, \code{'summary.sfapanel1'}, or
#' \code{'summary.sfalcmpanel'} \code{\link{coef}} extracts the 
#' coefficients, their standard errors, z-values, and (asymptotic) P-values.
#'
#' From on object of class \code{'sfacross'}, \code{'sfalcmcross'}, 
#' \code{'sfagzisfcross'}, \code{'sfacnsfcross'}, \code{'sfamisfcross'}, 
#' \code{'sfazisfcross'}, \code{'sfametacross'}, \code{'sfaselectioncross'}, 
#' \code{'sfapanel1'}, or \code{'sfalcmpanel'} it extracts only the estimated 
#' coefficients.
#'
#' @name coef
#' 
#' @aliases coef.sfacross coef.sfalcmcross coef.sfagzisfcross coef.sfacnsfcross 
#' coef.sfamisfcross coef.sfazisfcross coef.sfaselectioncross coef.sfametacross 
#' coef.sfapanel1 coef.sfalcmpanel coef.summary.sfacross 
#' coef.summary.sfalcmcross coef.summary.sfagzisfcross coef.summary.sfacnsfcross 
#' coef.summary.sfamisfcross coef.summary.sfazisfcross 
#' coef.summary.sfaselectioncross coef.summary.sfametacross
#' coef.summary.sfapanel1 coef.summary.sfalcmpanel
#'
#' @param object A stochastic frontier model returned by \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, 
#' \code{\link{sfaselectioncross}}, \code{\link{sfapanel1}}, or 
#' \code{\link{sfalcmpanel}} or an object of class \code{'summary.sfacross'}, 
#' \code{'summary.sfalcmcross'}, \code{'summary.sfagzisfcross'}, 
#' \code{'summary.sfacnsfcross'}, \code{'summary.sfamisfcross'}, 
#' \code{'summary.sfazisfcross'}, \code{'summary.sfametacross'}, 
#' \code{'summary.sfaselectioncross'}, \code{'summary.sfapanel1'} or 
#' \code{'summary.sfalcmpanel'}.
#' @param extraPar Logical (default = \code{FALSE}). If \code{TRUE}, additional
#' parameters are returned:
#'
#' \code{sigmaSq} = \code{sigmauSq} + \code{sigmavSq}
#'
#' \code{lambdaSq} = \code{sigmauSq}/\code{sigmavSq}
#'
#' \code{sigmauSq} = \eqn{\exp{(Wu)}} = \eqn{\exp{(\bm{\delta}' \mathbf{Z}_u)}}
#'
#' \code{sigmavSq} = \eqn{\exp{(Wv)}} = \eqn{\exp{(\bm{\phi}' \mathbf{Z}_v)}}
#'
#' \code{sigma} = \code{sigmaSq}^0.5
#'
#' \code{lambda} = \code{lambdaSq}^0.5
#'
#' \code{sigmau} = \code{sigmauSq}^0.5
#'
#' \code{sigmav} = \code{sigmavSq}^0.5
#'
#' \code{gamma} = \code{sigmauSq}/(\code{sigmauSq} + \code{sigmavSq})
#' @param ... Currently ignored.
#'
#' @return For objects of class \code{'summary.sfacross'}, 
#' \code{'summary.sfalcmcross'}, \code{'summary.sfagzisfcross'},
#' \code{'summary.sfacnsfcross'}, \code{'summary.sfamisfcross'}, 
#' \code{'summary.sfazisfcross'}, \code{'summary.sfaselectioncross'}, 
#' \code{'summary.sfapanel1'} or \code{'summary.sfalcmpanel'}, 
#' \code{\link{coef}} returns a matrix with four columns, namely, the estimated 
#' coefficients, their standard errors, z-values, and (asymptotic) P-values. For
#' objects of class \code{'summary.metasfacross'}, a list of four-column 
#' matrices is returned.
#'
#' For objects of class \code{'sfacross'}, \code{'sfalcmcross'},
#' \code{'sfagzisfcross'}, \code{'sfacnsfcross'},  \code{'sfamisfcross'}, 
#' \code{'sfazisfcross'}, \code{'sfaselectioncross'}, \code{'sfapanel1'}, 
#' \code{'sfalcmpanel'}, \code{\link{coef}} returns a numeric vector of the 
#' estimated coefficients. If \code{extraPar = TRUE}, additional parameters, 
#' detailed in the section \sQuote{Arguments}, are also returned. In the case 
#' of object of class \code{'sfalcmcross'}, , each additional parameter ends with 
#' \code{'#'} that 
#' represents the class number.
#' 
#' For an object of class \code{'sfametacross'}, \code{\link{coef}} returns a 
#' matrix of the estimated coefficients for each group and the metafrontier (
#' whenever available). If \code{extraPar = TRUE}, additional parameters, 
#' detailed in the section \sQuote{Arguments}, are also returned.
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
#' stochastic frontier analysis model fitting function using cross-sectional or 
#' pooled data.
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
#' @keywords methods coefficients
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
#' coef(tl_u_ts, extraPar = TRUE)
#' coef(summary(tl_u_ts))
#' }
#'
#' @export
# @exportS3Method coef sfacross coefficients from sfacross ----------
coef.sfacross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value", call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    if (object$udist == "tnormal") {
      if (object$scaling) {
        beta <- object$mlParam[1:(object$nXvar)]
        delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + (object$nuZUvar -
          1))]
        tau <- object$mlParam[object$nXvar + (object$nuZUvar - 1) + 1]
        cu <- object$mlParam[object$nXvar + (object$nuZUvar - 1) + 2]
        phi <- object$mlParam[(object$nXvar + (object$nuZUvar - 1) + 2 +
          1):(object$nXvar + (object$nuZUvar - 1) + 2 + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
        Wu <- cu + 2 * as.numeric(crossprod(matrix(delta), t(uHvar[, -1])))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      } else {
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar)]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
        Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      }
    } else {
      if (object$udist == "lognormal") {
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar)]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
        Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      } else {
        delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
        phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
          object$nuZUvar + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
        vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
        Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      }
    }
    if (object$udist == "lognormal" || (object$udist == "tnormal" & object$scaling ==
      FALSE)) {
      if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar > 1)
        cat("Variances averaged over observations     \n\n")
    } else {
      if (object$nuZUvar > 1 || object$nvZVvar > 1)
        cat("Variances averaged over observations     \n\n")
    }
    cRes <- c(cRes, sigmaSq = mean(exp(Wu)) + mean(exp(Wv)), lambdaSq = mean(exp(Wu))/mean(exp(Wv)),
      sigmauSq = mean(exp(Wu)), sigmavSq = mean(exp(Wv)), sigma = sqrt(mean(exp(Wu)) +
        mean(exp(Wv))), lambda = sqrt(mean(exp(Wu))/mean(exp(Wv))), sigmau = sqrt(mean(exp(Wu))),
      sigmav = sqrt(mean(exp(Wv))), gamma = mean(exp(Wu))/(mean(exp(Wu)) +
        mean(exp(Wv))))
  }
  return(cRes)
}

# coefficients from summary.sfacross ----------
#' @rdname coef
#' @export
# @exportS3Method coef summary.sfacross
coef.summary.sfacross <- function(object, ...) {
  object$mlRes
}

# coefficients from sfalcmcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef sfalcmcross
coef.sfalcmcross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value", call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
    vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
    delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
    phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
      object$nuZUvar + object$nvZVvar)]
    delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + object$nvZVvar +
      1):(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
    phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
      1):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
    Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
    Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
    Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
    Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
    if (object$nClasses == 3) {
      delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 2 *
        object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar + 2 *
        object$nvZVvar)]
      phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar +
        1):(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar)]
      Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
      Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
    } else {
      if (object$nClasses == 4) {
        delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
          2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
          2 * object$nvZVvar)]
        phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 *
          object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar + 3 *
          object$nvZVvar)]
        delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar +
          3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
          3 * object$nvZVvar)]
        phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 3 *
          object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 4 *
          object$nvZVvar)]
        Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
        Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
        Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
        Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
      } else {
        if (object$nClasses == 5) {
          delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
          2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
          2 * object$nvZVvar)]
          phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
          2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
          3 * object$nvZVvar)]
          delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar +
          3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
          3 * object$nvZVvar)]
          phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
          3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
          4 * object$nvZVvar)]
          delta5 <- object$mlParam[(5 * object$nXvar + 4 * object$nuZUvar +
          4 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
          4 * object$nvZVvar)]
          phi5 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
          4 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
          5 * object$nvZVvar)]
          Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
          Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
          Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
          Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
          Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar)))
          Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
        }
      }
    }
    if (object$nuZUvar > 1 || object$nvZVvar > 1)
      cat("Variances averaged over observations     \n\n")
    cRes <- c(cRes, sigmaSq1 = mean(exp(Wu1)) + mean(exp(Wv1)), lambdaSq1 = mean(exp(Wu1))/mean(exp(Wv1)),
      sigmauSq1 = mean(exp(Wu1)), sigmavSq1 = mean(exp(Wv1)), sigma1 = sqrt(mean(exp(Wu1)) +
        mean(exp(Wv1))), lambda1 = sqrt(mean(exp(Wu1))/mean(exp(Wv1))), sigmau1 = sqrt(mean(exp(Wu1))),
      sigmav1 = sqrt(mean(exp(Wv1))), gamma1 = mean(exp(Wu1))/(mean(exp(Wu1)) +
        mean(exp(Wv1))), sigmaSq2 = mean(exp(Wu2)) + mean(exp(Wv2)), lambdaSq2 = mean(exp(Wu2))/mean(exp(Wv2)),
      sigmauSq2 = mean(exp(Wu2)), sigmavSq2 = mean(exp(Wv2)), sigma2 = sqrt(mean(exp(Wu2)) +
        mean(exp(Wv2))), lambda2 = sqrt(mean(exp(Wu2))/mean(exp(Wv2))), sigmau2 = sqrt(mean(exp(Wu2))),
      sigmav2 = sqrt(mean(exp(Wv2))), gamma2 = mean(exp(Wu2))/(mean(exp(Wu2)) +
        mean(exp(Wv2))))
    if (object$nClasses == 3) {
      cRes <- c(cRes, sigmaSq3 = mean(exp(Wu3)) + mean(exp(Wv3)), lambdaSq3 = mean(exp(Wu3))/mean(exp(Wv3)),
        sigmauSq3 = mean(exp(Wu3)), sigmavSq3 = mean(exp(Wv3)), sigma3 = sqrt(mean(exp(Wu3)) +
          mean(exp(Wv3))), lambda3 = sqrt(mean(exp(Wu3))/mean(exp(Wv3))),
        sigmau3 = sqrt(mean(exp(Wu3))), sigmav3 = sqrt(mean(exp(Wv3))), gamma3 = mean(exp(Wu3))/(mean(exp(Wu3)) +
          mean(exp(Wv3))))
    } else {
      if (object$nClasses == 4) {
        cRes <- c(cRes, sigmaSq3 = mean(exp(Wu3)) + mean(exp(Wv3)), lambdaSq3 = mean(exp(Wu3))/mean(exp(Wv3)),
          sigmauSq3 = mean(exp(Wu3)), sigmavSq3 = mean(exp(Wv3)), sigma3 = sqrt(mean(exp(Wu3)) +
          mean(exp(Wv3))), lambda3 = sqrt(mean(exp(Wu3))/mean(exp(Wv3))),
          sigmau3 = sqrt(mean(exp(Wu3))), sigmav3 = sqrt(mean(exp(Wv3))),
          gamma3 = mean(exp(Wu3))/(mean(exp(Wu3)) + mean(exp(Wv3))), sigmaSq4 = mean(exp(Wu4)) +
          mean(exp(Wv4)), lambdaSq4 = mean(exp(Wu4))/mean(exp(Wv4)), sigmauSq4 = mean(exp(Wu4)),
          sigmavSq4 = mean(exp(Wv4)), sigma4 = sqrt(mean(exp(Wu4)) + mean(exp(Wv4))),
          lambda4 = sqrt(mean(exp(Wu4))/mean(exp(Wv4))), sigmau4 = sqrt(mean(exp(Wu4))),
          sigmav4 = sqrt(mean(exp(Wv4))), gamma4 = mean(exp(Wu4))/(mean(exp(Wu4)) +
          mean(exp(Wv4))))
      } else {
        if (object$nClasses == 5) {
          cRes <- c(cRes, sigmaSq3 = mean(exp(Wu3)) + mean(exp(Wv3)), lambdaSq3 = mean(exp(Wu3))/mean(exp(Wv3)),
          sigmauSq3 = mean(exp(Wu3)), sigmavSq3 = mean(exp(Wv3)), sigma3 = sqrt(mean(exp(Wu3)) +
            mean(exp(Wv3))), lambda3 = sqrt(mean(exp(Wu3))/mean(exp(Wv3))),
          sigmau3 = sqrt(mean(exp(Wu3))), sigmav3 = sqrt(mean(exp(Wv3))),
          gamma3 = mean(exp(Wu3))/(mean(exp(Wu3)) + mean(exp(Wv3))), sigmaSq4 = mean(exp(Wu4)) +
            mean(exp(Wv4)), lambdaSq4 = mean(exp(Wu4))/mean(exp(Wv4)),
          sigmauSq4 = mean(exp(Wu4)), sigmavSq4 = mean(exp(Wv4)), sigma4 = sqrt(mean(exp(Wu4)) +
            mean(exp(Wv4))), lambda4 = sqrt(mean(exp(Wu4))/mean(exp(Wv4))),
          sigmau4 = sqrt(mean(exp(Wu4))), sigmav4 = sqrt(mean(exp(Wv4))),
          gamma4 = mean(exp(Wu4))/(mean(exp(Wu4)) + mean(exp(Wv4))), sigmaSq5 = mean(exp(Wu5)) +
            mean(exp(Wv5)), lambdaSq5 = mean(exp(Wu5))/mean(exp(Wv5)),
          sigmauSq5 = mean(exp(Wu5)), sigmavSq5 = mean(exp(Wv5)), sigma5 = sqrt(mean(exp(Wu5)) +
            mean(exp(Wv5))), lambda5 = sqrt(mean(exp(Wu5))/mean(exp(Wv5))),
          sigmau5 = sqrt(mean(exp(Wu5))), sigmav5 = sqrt(mean(exp(Wv5))),
          gamma5 = mean(exp(Wu5))/(mean(exp(Wu5)) + mean(exp(Wv5))))
        }
      }
    }
  }
  return(cRes)
}

# coefficients from summary.sfalcmcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef summary.sfalcmcross
coef.summary.sfalcmcross <- function(...) {
  coef.summary.sfacross(...)
}

# coefficients from sfagzisfcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef sfagzisfcross
coef.sfagzisfcross <- function(...) {
  coef.sfalcmcross(...)
}

# coefficients from summary.sfagzisfcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef summary.sfagzisfcross
coef.summary.sfagzisfcross <- function(...) {
  coef.summary.sfacross(...)
}

# coefficients from sfacnsfcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef sfacnsfcross
coef.sfacnsfcross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value", call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    GROUPS <- efficiencies.sfacnsfcross(object = object)$Group_c
    if (object$sigmauType == "common") {
      if (object$udist == "tnormal") {
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar)]
        phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
        phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
          2 * object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
        Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
        Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
      } else {
        if (object$udist == "lognormal") {
          delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar)]
          phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
          phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
          2 * object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 4)
          Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
          Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
          Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
        } else {
          delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
          phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
          object$nuZUvar + object$nvZVvar)]
          phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
          1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 2)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
          Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
          Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
          Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
        }
      }
      if (object$udist == "lognormal" || object$udist == "tnormal") {
        if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar >
          1)
          cat("Variances averaged over observations     \n\n")
      } else {
        if (object$nuZUvar > 1 || object$nvZVvar > 1)
          cat("Variances averaged over observations     \n\n")
      }
      cRes <- c(cRes, sigmaSq1 = mean(exp(Wu[GROUPS == 1])) + mean(exp(Wv1[GROUPS ==
        1])), sigmaSq2 = mean(exp(Wu[GROUPS == 1])) + mean(exp(Wv2[GROUPS ==
        1])), lambdaSq1 = mean(exp(Wu[GROUPS == 1]))/mean(exp(Wv1[GROUPS ==
        1])), lambdaSq2 = mean(exp(Wu[GROUPS == 1]))/mean(exp(Wv2[GROUPS ==
        1])), sigmauSq = mean(exp(Wu[GROUPS == 1])), sigmavSq1 = mean(exp(Wv1[GROUPS ==
        1])), sigmavSq2 = mean(exp(Wv2[GROUPS == 2])), sigma1 = sqrt(mean(exp(Wu[GROUPS ==
        1])) + mean(exp(Wv1[GROUPS == 1]))), sigma2 = sqrt(mean(exp(Wu[GROUPS ==
        1])) + mean(exp(Wv2[GROUPS == 1]))), lambda1 = sqrt(mean(exp(Wu[GROUPS ==
        1]))/mean(exp(Wv1[GROUPS == 1]))), lambda2 = sqrt(mean(exp(Wu[GROUPS ==
        1]))/mean(exp(Wv2[GROUPS == 1]))), sigmau = sqrt(mean(exp(Wu[GROUPS ==
        1]))), sigmav1 = sqrt(mean(exp(Wv1[GROUPS == 1]))), sigmav2 = sqrt(mean(exp(Wv2[GROUPS ==
        2]))), gamma1 = mean(exp(Wu[GROUPS == 1]))/(mean(exp(Wu[GROUPS ==
        1])) + mean(exp(Wv1[GROUPS == 1]))), gamma2 = mean(exp(Wu[GROUPS ==
        1]))/(mean(exp(Wu[GROUPS == 1])) + mean(exp(Wv2[GROUPS == 1]))))
    } else {
      if (object$sigmauType == "different") {
        if (object$udist == "tnormal") {
          delta1 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar)]
          delta2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar)]
          phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
          phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar +
          object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar +
          2 * object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 4)
          Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
          Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
          Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
          Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
        } else {
          if (object$udist == "lognormal") {
          delta1 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
            object$nmuZUvar + object$nuZUvar)]
          delta2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
            1):(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar)]
          phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + 2 *
            object$nuZUvar + 1):(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar +
            object$nvZVvar)]
          phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 2 *
            object$nuZUvar + object$nvZVvar + 1):(object$nXvar + object$nmuZUvar +
            2 * object$nuZUvar + 2 * object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
            rhs = 3)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
            rhs = 4)
          Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
          Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
          Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
          Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
          } else {
          delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
          delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
            2 * object$nuZUvar)]
          phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
            2 * object$nuZUvar + object$nvZVvar)]
          phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
            1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
            rhs = 2)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
            rhs = 3)
          Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
          Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
          Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
          Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
          }
        }
        if (object$udist == "lognormal" || object$udist == "tnormal") {
          if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar >
          1)
          cat("Variances averaged over observations     \n\n")
        } else {
          if (object$nuZUvar > 1 || object$nvZVvar > 1)
          cat("Variances averaged over observations     \n\n")
        }
        cRes <- c(cRes, sigmaSq1 = mean(exp(Wu1[GROUPS == 1])) + mean(exp(Wv1[GROUPS ==
          1])), sigmaSq2 = mean(exp(Wu2[GROUPS == 1])) + mean(exp(Wv2[GROUPS ==
          1])), lambdaSq1 = mean(exp(Wu1[GROUPS == 1]))/mean(exp(Wv1[GROUPS ==
          1])), lambdaSq2 = mean(exp(Wu2[GROUPS == 1]))/mean(exp(Wv2[GROUPS ==
          1])), sigmauSq1 = mean(exp(Wu1[GROUPS == 1])), sigmauSq2 = mean(exp(Wu2[GROUPS ==
          1])), sigmavSq1 = mean(exp(Wv1[GROUPS == 1])), sigmavSq2 = mean(exp(Wv2[GROUPS ==
          2])), sigma1 = sqrt(mean(exp(Wu1[GROUPS == 1])) + mean(exp(Wv1[GROUPS ==
          1]))), sigma2 = sqrt(mean(exp(Wu2[GROUPS == 1])) + mean(exp(Wv2[GROUPS ==
          1]))), lambda1 = sqrt(mean(exp(Wu1[GROUPS == 1]))/mean(exp(Wv1[GROUPS ==
          1]))), lambda2 = sqrt(mean(exp(Wu2[GROUPS == 1]))/mean(exp(Wv2[GROUPS ==
          1]))), sigmau1 = sqrt(mean(exp(Wu1[GROUPS == 1]))), sigmau2 = sqrt(mean(exp(Wu2[GROUPS ==
          1]))), sigmav1 = sqrt(mean(exp(Wv1[GROUPS == 1]))), sigmav2 = sqrt(mean(exp(Wv2[GROUPS ==
          2]))), gamma1 = mean(exp(Wu1[GROUPS == 1]))/(mean(exp(Wu1[GROUPS ==
          1])) + mean(exp(Wv1[GROUPS == 1]))), gamma2 = mean(exp(Wu2[GROUPS ==
          1]))/(mean(exp(Wu2[GROUPS == 1])) + mean(exp(Wv2[GROUPS == 1]))))
      }
    }
  }
  return(cRes)
}

# coefficients from summary.sfacnsfcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef summary.sfacnsfcross
coef.summary.sfacnsfcross <- function(...) {
  coef.summary.sfacross(...)
}

# coefficients from sfamisfcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef sfamisfcross
coef.sfamisfcross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value", call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    if (object$udist == "tnormal") {
      delta1 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
        object$nmuZUvar + object$nuZUvar)]
      delta2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
        1):(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar +
        1):(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
      uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
      vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
      Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
      Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
      Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
    } else {
      if (object$udist == "lognormal") {
        delta1 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar)]
        delta2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar)]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
        Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
        Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      } else {
        delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
        delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
          2 * object$nuZUvar)]
        phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
          2 * object$nuZUvar + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
        vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
        Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
        Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      }
    }
    if (object$udist == "lognormal" || object$udist == "tnormal") {
      if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar > 1)
        cat("Variances averaged over observations     \n\n")
    } else {
      if (object$nuZUvar > 1 || object$nvZVvar > 1)
        cat("Variances averaged over observations     \n\n")
    }
    GROUPS <- efficiencies.sfamisfcross(object = object)$Group_c
    cRes <- c(cRes, sigmaSq1 = mean(exp(Wu1[GROUPS == 1])) + mean(exp(Wv[GROUPS ==
      1])), sigmaSq2 = mean(exp(Wu2[GROUPS == 1])) + mean(exp(Wv[GROUPS ==
      1])), lambdaSq1 = mean(exp(Wu1[GROUPS == 1]))/mean(exp(Wv[GROUPS == 1])),
      lambdaSq2 = mean(exp(Wu2[GROUPS == 1]))/mean(exp(Wv[GROUPS == 1])), sigmauSq1 = mean(exp(Wu1[GROUPS ==
        1])), sigmauSq2 = mean(exp(Wu2[GROUPS == 1])), sigmavSq1 = mean(exp(Wv[GROUPS ==
        1])), sigmavSq2 = mean(exp(Wv[GROUPS == 2])), sigma1 = sqrt(mean(exp(Wu1[GROUPS ==
        1])) + mean(exp(Wv[GROUPS == 1]))), sigma2 = sqrt(mean(exp(Wu2[GROUPS ==
        1])) + mean(exp(Wv[GROUPS == 1]))), lambda1 = sqrt(mean(exp(Wu1[GROUPS ==
        1]))/mean(exp(Wv[GROUPS == 1]))), lambda2 = sqrt(mean(exp(Wu2[GROUPS ==
        1]))/mean(exp(Wv[GROUPS == 1]))), sigmau1 = sqrt(mean(exp(Wu1[GROUPS ==
        1]))), sigmau2 = sqrt(mean(exp(Wu2[GROUPS == 1]))), sigmav1 = sqrt(mean(exp(Wv[GROUPS ==
        1]))), sigmav2 = sqrt(mean(exp(Wv[GROUPS == 2]))), gamma1 = mean(exp(Wu1[GROUPS ==
        1]))/(mean(exp(Wu1[GROUPS == 1])) + mean(exp(Wv[GROUPS == 1]))),
      gamma2 = mean(exp(Wu2[GROUPS == 1]))/(mean(exp(Wu2[GROUPS == 1])) + mean(exp(Wv[GROUPS ==
        1]))))
  }
  return(cRes)
}

# coefficients from summary.sfamisfcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef summary.sfamisfcross
coef.summary.sfamisfcross <- function(...) {
  coef.summary.sfacross(...)
}

# coefficients from sfazisfcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef sfazisfcross
coef.sfazisfcross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value", call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    GROUPS <- efficiencies.sfazisfcross(object = object)$Group_c
    if (object$sigmavType == "common") {
      if (object$udist == "tnormal") {
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar)]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
        Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      } else {
        if (object$udist == "lognormal") {
          delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar)]
          phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 4)
          Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
          Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
        } else {
          delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
          phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
          object$nuZUvar + object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 2)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
          Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
          Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
        }
      }
      if (object$udist == "lognormal" || object$udist == "tnormal") {
        if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar >
          1)
          cat("Variances averaged over observations     \n\n")
      } else {
        if (object$nuZUvar > 1 || object$nvZVvar > 1)
          cat("Variances averaged over observations     \n\n")
      }
      cRes <- c(cRes, sigmaSq1 = mean(exp(Wu[GROUPS == 1])) + mean(exp(Wv[GROUPS ==
        1])), sigmaSq2 = mean(exp(Wv[GROUPS == 2])), lambdaSq1 = mean(exp(Wu[GROUPS ==
        1]))/mean(exp(Wv[GROUPS == 1])), lambdaSq2 = 0, sigmauSq1 = mean(exp(Wu[GROUPS ==
        1])), sigmauSq2 = 0, sigmavSq1 = mean(exp(Wv[GROUPS == 1])), sigmavSq2 = mean(exp(Wv[GROUPS ==
        2])), sigma1 = sqrt(mean(exp(Wu[GROUPS == 1])) + mean(exp(Wv[GROUPS ==
        1]))), sigma2 = sqrt(mean(exp(Wv[GROUPS == 2]))), lambda1 = sqrt(mean(exp(Wu[GROUPS ==
        1]))/mean(exp(Wv[GROUPS == 1]))), lambda2 = 0, sigmau1 = sqrt(mean(exp(Wu[GROUPS ==
        1]))), sigmau2 = 0, sigmav1 = sqrt(mean(exp(Wv[GROUPS == 1]))), sigmav2 = sqrt(mean(exp(Wv[GROUPS ==
        2]))), gamma1 = mean(exp(Wu[GROUPS == 1]))/(mean(exp(Wu[GROUPS ==
        1])) + mean(exp(Wv[GROUPS == 1]))), gamma2 = 0)
    } else {
      if (object$sigmavType == "different") {
        if (object$udist == "tnormal") {
          delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar)]
          phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
          phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
          2 * object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 4)
          Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
          Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
          Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
        } else {
          if (object$udist == "lognormal") {
          delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
            object$nmuZUvar + object$nuZUvar)]
          phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
            1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
          phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
            object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
            2 * object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
            rhs = 3)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
            rhs = 4)
          Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
          Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
          Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
          } else {
          delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
          phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
            1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
          phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
            object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
            2 * object$nvZVvar)]
          uHvar <- model.matrix(object$formula, data = object$dataTable,
            rhs = 2)
          vHvar <- model.matrix(object$formula, data = object$dataTable,
            rhs = 3)
          Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
          Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
          Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
          }
        }
        if (object$udist == "lognormal" || object$udist == "tnormal") {
          if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar >
          1)
          cat("Variances averaged over observations     \n\n")
        } else {
          if (object$nuZUvar > 1 || object$nvZVvar > 1)
          cat("Variances averaged over observations     \n\n")
        }
        cRes <- c(cRes, sigmaSq1 = mean(exp(Wu[GROUPS == 1])) + mean(exp(Wv1[GROUPS ==
          1])), sigmaSq2 = mean(exp(Wv2[GROUPS == 2])), lambdaSq1 = mean(exp(Wu[GROUPS ==
          1]))/mean(exp(Wv1[GROUPS == 1])), lambdaSq2 = 0, sigmauSq1 = mean(exp(Wu[GROUPS ==
          1])), sigmauSq2 = 0, sigmavSq1 = mean(exp(Wv1[GROUPS == 1])), sigmavSq2 = mean(exp(Wv2[GROUPS ==
          2])), sigma1 = sqrt(mean(exp(Wu[GROUPS == 1])) + mean(exp(Wv1[GROUPS ==
          1]))), sigma2 = sqrt(mean(exp(Wv2[GROUPS == 2]))), lambda1 = sqrt(mean(exp(Wu[GROUPS ==
          1]))/mean(exp(Wv1[GROUPS == 1]))), lambda2 = 0, sigmau1 = sqrt(mean(exp(Wu[GROUPS ==
          1]))), sigmau2 = 0, sigmav1 = sqrt(mean(exp(Wv1[GROUPS == 1]))),
          sigmav2 = sqrt(mean(exp(Wv2[GROUPS == 2]))), gamma1 = mean(exp(Wu[GROUPS ==
          1]))/(mean(exp(Wu[GROUPS == 1])) + mean(exp(Wv1[GROUPS == 1]))),
          gamma2 = 0)
      }
    }
  }
  return(cRes)
}

# coefficients from summary.sfazisfcross ----------
#' @rdname coef
#' @export
# @exportS3Method coef summary.sfazisfcross
coef.summary.sfazisfcross <- function(...) {
  coef.summary.sfacross(...)
}

# coefficients from sfaselectioncross ----------
#' @rdname coef
#' @export
# @exportS3Method coef sfaselectioncross
coef.sfaselectioncross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value", call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
    phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
      object$nuZUvar + object$nvZVvar)]
    uHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable$selectDum ==
      1, ], rhs = 2)
    vHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable$selectDum ==
      1, ], rhs = 3)
    Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
    Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
    if (object$nuZUvar > 1 || object$nvZVvar > 1)
      cat("Variances averaged over observations     \n\n")
    cRes <- c(cRes, sigmaSq = mean(exp(Wu)) + mean(exp(Wv)), lambdaSq = mean(exp(Wu))/mean(exp(Wv)),
      sigmauSq = mean(exp(Wu)), sigmavSq = mean(exp(Wv)), sigma = sqrt(mean(exp(Wu)) +
        mean(exp(Wv))), lambda = sqrt(mean(exp(Wu))/mean(exp(Wv))), sigmau = sqrt(mean(exp(Wu))),
      sigmav = sqrt(mean(exp(Wv))), gamma = mean(exp(Wu))/(mean(exp(Wu)) +
        mean(exp(Wv))))
  }
  return(cRes)
}

# coefficients from summary.sfaselectioncross ----------
#' @rdname coef
#' @export
# @exportS3Method coef summary.sfaselectioncross
coef.summary.sfaselectioncross <- function(...) {
  coef.summary.sfacross(...)
}

# coefficients from sfametacross ----------
#' @rdname coef
#' @export
# @exportS3Method coef sfametacross
coef.sfametacross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value", call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    group_var <- object$dataTable[[object$Ngroup + 1]][object$name_meta_var][,
      1]
    group_var_list <- sort(unique(group_var))
    Wu <- list()
    Wv <- list()
    for (g in group_var_list) {
      if (object$udist %in% c("tnormal", "lognormal")) {
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar), which(group_var_list == g)]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar),
          which(group_var_list == g)]
        uHvar <- model.matrix(object$formula, data = object$dataTable[[which(group_var_list ==
          g)]], rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable[[which(group_var_list ==
          g)]], rhs = 4)
        Wu[[which(group_var_list == g)]] <- as.numeric(crossprod(matrix(delta),
          t(uHvar)))
        Wv[[which(group_var_list == g)]] <- as.numeric(crossprod(matrix(phi),
          t(vHvar)))
        if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar >
          1)
          cat("Variances averaged over observations     \n\n")
      } else {
        delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar),
          which(group_var_list == g)]
        phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
          object$nuZUvar + object$nvZVvar), which(group_var_list == g)]
        uHvar <- model.matrix(object$formula, data = object$dataTable[[which(group_var_list ==
          g)]], rhs = 2)
        vHvar <- model.matrix(object$formula, data = object$dataTable[[which(group_var_list ==
          g)]], rhs = 3)
        Wu[[which(group_var_list == g)]] <- as.numeric(crossprod(matrix(delta),
          t(uHvar)))
        Wv[[which(group_var_list == g)]] <- as.numeric(crossprod(matrix(phi),
          t(vHvar)))
        if (object$nuZUvar > 1 || object$nvZVvar > 1)
          cat("Variances averaged over observations     \n\n")
      }
    }
    cInt <- rbind(sapply(Wu, FUN = function(x) mean(exp(x))) + sapply(Wv, FUN = function(x) mean(exp(x))),
      sapply(Wu, FUN = function(x) mean(exp(x)))/sapply(Wv, FUN = function(x) mean(exp(x))),
      sapply(Wu, FUN = function(x) mean(exp(x))), sapply(Wv, FUN = function(x) mean(exp(x))),
      sqrt(sapply(Wu, FUN = function(x) mean(exp(x))) + sapply(Wv, FUN = function(x) mean(exp(x)))),
      sqrt(sapply(Wu, FUN = function(x) mean(exp(x)))/sapply(Wv, FUN = function(x) mean(exp(x)))),
      sqrt(sapply(Wu, FUN = function(x) mean(exp(x)))), sqrt(sapply(Wv, FUN = function(x) mean(exp(x)))),
      sapply(Wu, FUN = function(x) mean(exp(x)))/(sapply(Wu, FUN = function(x) mean(exp(x))) +
        sapply(Wv, FUN = function(x) mean(exp(x)))))
    if (object$modelType %in% c("hhl14")) {
      if (object$udist %in% c("tnormal", "lognormal")) {
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
          object$nmuZUvar + object$nuZUvar), object$Ngroup + 1]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar),
          object$Ngroup + 1]
        uHvar <- model.matrix(object$formula, data = object$dataTable[[object$Ngroup +
          1]], rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable[[object$Ngroup +
          1]], rhs = 4)
        Wu[[object$Ngroup + 1]] <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv[[object$Ngroup + 1]] <- as.numeric(crossprod(matrix(phi), t(vHvar)))
        if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar >
          1)
          cat("Variances averaged over observations     \n\n")
      } else {
        delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar),
          object$Ngroup + 1]
        phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
          object$nuZUvar + object$nvZVvar), object$Ngroup + 1]
        uHvar <- model.matrix(object$formula, data = object$dataTable[[object$Ngroup +
          1]], rhs = 2)
        vHvar <- model.matrix(object$formula, data = object$dataTable[[object$Ngroup +
          1]], rhs = 3)
        Wu[[object$Ngroup + 1]] <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv[[object$Ngroup + 1]] <- as.numeric(crossprod(matrix(phi), t(vHvar)))
        if (object$nuZUvar > 1 || object$nvZVvar > 1)
          cat("Variances averaged over observations     \n\n")
      }
      cInt <- cbind(cInt, mean(exp(Wu[[object$Ngroup + 1]])) + mean(exp(Wv[[object$Ngroup +
        1]])), mean(exp(Wu[[object$Ngroup + 1]]))/mean(exp(Wv[[object$Ngroup +
        1]])), mean(exp(Wu[[object$Ngroup + 1]])), mean(exp(Wv[[object$Ngroup +
        1]])), sqrt(mean(exp(Wu[[object$Ngroup + 1]])) + mean(exp(Wv[[object$Ngroup +
        1]]))), sqrt(mean(exp(Wu[[object$Ngroup + 1]]))/mean(exp(Wv[[object$Ngroup +
        1]]))), sqrt(mean(exp(Wu[[object$Ngroup + 1]]))), sqrt(mean(exp(Wv[[object$Ngroup +
        1]]))), mean(exp(Wu[[object$Ngroup + 1]]))/(mean(exp(Wu[[object$Ngroup +
        1]])) + mean(exp(Wv[[object$Ngroup + 1]]))))
    } else {
      cInt <- cbind(cInt, NA)
    }
    row.names(cInt) <- c("sigmaSq", "lambdaSq", "sigmauSq", "sigmavSq", "sigma",
      "lambda", "sigmau", "sigmav", "gamma")
    cRes <- rbind(cRes, cInt)
  }
  return(cRes)
}

# coefficients from summary.sfametacross ----------
#' @rdname coef
#' @export
# @exportS3Method coef summary.sfametacross
coef.summary.sfametacross <- function(...) {
  coef.summary.sfacross(...)
}

# coefficients from sfapanel1 ----------
#' @rdname coef
#' @aliases coef.sfapanel1
#' @export
coef.sfapanel1 <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value", call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    if (object$udist %in% c("tnormal", "lognormal")) {
      delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
        object$nmuZUvar + object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
        1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
      uHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
      vHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
    } else {
      delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
        object$nuZUvar + object$nvZVvar)]
      uHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
      vHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
    }
    pindex <- object$dataTable[, 1:2]
    invariance <- object$invariance
    if (invariance == 1) {
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
    } else {
      if (invariance == 2) {
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
      } else {
        if (invariance == 3) {
          uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
        }
      }
    }
    Wu <- as.numeric(crossprod(matrix(delta), t(uHvar_p)))
    Wv <- as.numeric(crossprod(matrix(phi), t(vHvar_p)))
    if (object$udist == "lognormal" || (object$udist == "tnormal" & object$scaling ==
      FALSE)) {
      if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar > 1)
        cat("Variances averaged over observations     \n\n")
    } else {
      if (object$nuZUvar > 1 || object$nvZVvar > 1)
        cat("Variances averaged over observations     \n\n")
    }
    cRes <- c(cRes, sigmaSq = mean(exp(Wu)) + mean(exp(Wv)), lambdaSq = mean(exp(Wu))/mean(exp(Wv)),
      sigmauSq = mean(exp(Wu)), sigmavSq = mean(exp(Wv)), sigma = sqrt(mean(exp(Wu)) +
        mean(exp(Wv))), lambda = sqrt(mean(exp(Wu))/mean(exp(Wv))), sigmau = sqrt(mean(exp(Wu))),
      sigmav = sqrt(mean(exp(Wv))), gamma = mean(exp(Wu))/(mean(exp(Wu)) +
        mean(exp(Wv))))
  }
  return(cRes)
}

# coefficients from summary.sfapanel1 ----------
#' @rdname coef
#' @aliases coef.summary.sfapanel1
#' @export
coef.summary.sfapanel1 <- function(...) {
  coef.summary.sfacross(...)
}

# coefficients from sfalcmpanel ----------
#' @rdname coef
#' @aliases coef.sfalcmpanel
#' @export
coef.sfalcmpanel <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value", call. = FALSE)
  cRes <- object$mlParam
  nParamExtra <- if (object$modelType %in% c("mbc92", "k90")) {
    2
  } else {
    if (object$modelType %in% c("bc92a", "bc92b", "bc92c", "mbc92", "k90", "kw05")) {
      ncol(object$gHvar)
    } else {
      0
    }
  }
  if (extraPar) {
    uHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
    vHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
    delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
    phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
      object$nuZUvar + object$nvZVvar)]
    delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + object$nvZVvar +
      nParamExtra + 1):(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
      nParamExtra)]
    phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
      nParamExtra + 1):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
      nParamExtra)]
    pindex <- object$dataTable[, 1:2]
    invariance <- object$invariance
    if (invariance == 1) {
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
    } else {
      if (invariance == 2) {
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
      } else {
        if (invariance == 3) {
          uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
        }
      }
    }
    Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar_p)))
    Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar_p)))
    Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar_p)))
    Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar_p)))
    if (object$nClasses == 3) {
      delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 2 *
        object$nvZVvar + 2 * nParamExtra + 1):(3 * object$nXvar + 3 * object$nuZUvar +
        2 * object$nvZVvar + 2 * nParamExtra)]
      phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar +
        2 * nParamExtra + 1):(3 * object$nXvar + 3 * object$nuZUvar + 3 *
        object$nvZVvar + 2 * nParamExtra)]
      Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar_p)))
      Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar_p)))
    } else {
      if (object$nClasses == 4) {
        delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
          2 * object$nvZVvar + 2 * nParamExtra + 1):(3 * object$nXvar + 3 *
          object$nuZUvar + 2 * object$nvZVvar + 2 * nParamExtra)]
        phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 *
          object$nvZVvar + 2 * nParamExtra + 1):(3 * object$nXvar + 3 * object$nuZUvar +
          3 * object$nvZVvar + 2 * nParamExtra)]
        delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar +
          3 * object$nvZVvar + 3 * nParamExtra + 1):(4 * object$nXvar + 4 *
          object$nuZUvar + 3 * object$nvZVvar + 3 * nParamExtra)]
        phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 3 *
          object$nvZVvar + 3 * nParamExtra + 1):(4 * object$nXvar + 4 * object$nuZUvar +
          4 * object$nvZVvar + 3 * nParamExtra)]
        Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar_p)))
        Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar_p)))
        Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar_p)))
        Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar_p)))
      } else {
        if (object$nClasses == 5) {
          delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
          2 * object$nvZVvar + 2 * nParamExtra + 1):(3 * object$nXvar +
          3 * object$nuZUvar + 2 * object$nvZVvar + 2 * nParamExtra)]
          phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
          2 * object$nvZVvar + 2 * nParamExtra + 1):(3 * object$nXvar +
          3 * object$nuZUvar + 3 * object$nvZVvar + 2 * nParamExtra)]
          delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar +
          3 * object$nvZVvar + 3 * nParamExtra + 1):(4 * object$nXvar +
          4 * object$nuZUvar + 3 * object$nvZVvar + 3 * nParamExtra)]
          phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
          3 * object$nvZVvar + 3 * nParamExtra + 1):(4 * object$nXvar +
          4 * object$nuZUvar + 4 * object$nvZVvar + 3 * nParamExtra)]
          delta5 <- object$mlParam[(5 * object$nXvar + 4 * object$nuZUvar +
          4 * object$nvZVvar + 4 * nParamExtra + 1):(5 * object$nXvar +
          5 * object$nuZUvar + 4 * object$nvZVvar + 4 * nParamExtra)]
          phi5 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
          4 * object$nvZVvar + 4 * nParamExtra + 1):(5 * object$nXvar +
          5 * object$nuZUvar + 5 * object$nvZVvar + 4 * nParamExtra)]
          Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar_p)))
          Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar_p)))
          Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar_p)))
          Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar_p)))
          Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar_p)))
          Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar_p)))
        }
      }
    }
    if (object$nuZUvar > 1 || object$nvZVvar > 1)
      cat("Variances averaged over observations     \n\n")
    cRes <- c(cRes, sigmaSq1 = mean(exp(Wu1)) + mean(exp(Wv1)), lambdaSq1 = mean(exp(Wu1))/mean(exp(Wv1)),
      sigmauSq1 = mean(exp(Wu1)), sigmavSq1 = mean(exp(Wv1)), sigma1 = sqrt(mean(exp(Wu1)) +
        mean(exp(Wv1))), lambda1 = sqrt(mean(exp(Wu1))/mean(exp(Wv1))), sigmau1 = sqrt(mean(exp(Wu1))),
      sigmav1 = sqrt(mean(exp(Wv1))), gamma1 = mean(exp(Wu1))/(mean(exp(Wu1)) +
        mean(exp(Wv1))), sigmaSq2 = mean(exp(Wu2)) + mean(exp(Wv2)), lambdaSq2 = mean(exp(Wu2))/mean(exp(Wv2)),
      sigmauSq2 = mean(exp(Wu2)), sigmavSq2 = mean(exp(Wv2)), sigma2 = sqrt(mean(exp(Wu2)) +
        mean(exp(Wv2))), lambda2 = sqrt(mean(exp(Wu2))/mean(exp(Wv2))), sigmau2 = sqrt(mean(exp(Wu2))),
      sigmav2 = sqrt(mean(exp(Wv2))), gamma2 = mean(exp(Wu2))/(mean(exp(Wu2)) +
        mean(exp(Wv2))))
    if (object$nClasses == 3) {
      cRes <- c(cRes, sigmaSq3 = mean(exp(Wu3)) + mean(exp(Wv3)), lambdaSq3 = mean(exp(Wu3))/mean(exp(Wv3)),
        sigmauSq3 = mean(exp(Wu3)), sigmavSq3 = mean(exp(Wv3)), sigma3 = sqrt(mean(exp(Wu3)) +
          mean(exp(Wv3))), lambda3 = sqrt(mean(exp(Wu3))/mean(exp(Wv3))),
        sigmau3 = sqrt(mean(exp(Wu3))), sigmav3 = sqrt(mean(exp(Wv3))), gamma3 = mean(exp(Wu3))/(mean(exp(Wu3)) +
          mean(exp(Wv3))))
    } else {
      if (object$nClasses == 4) {
        cRes <- c(cRes, sigmaSq3 = mean(exp(Wu3)) + mean(exp(Wv3)), lambdaSq3 = mean(exp(Wu3))/mean(exp(Wv3)),
          sigmauSq3 = mean(exp(Wu3)), sigmavSq3 = mean(exp(Wv3)), sigma3 = sqrt(mean(exp(Wu3)) +
          mean(exp(Wv3))), lambda3 = sqrt(mean(exp(Wu3))/mean(exp(Wv3))),
          sigmau3 = sqrt(mean(exp(Wu3))), sigmav3 = sqrt(mean(exp(Wv3))),
          gamma3 = mean(exp(Wu3))/(mean(exp(Wu3)) + mean(exp(Wv3))), sigmaSq4 = mean(exp(Wu4)) +
          mean(exp(Wv4)), lambdaSq4 = mean(exp(Wu4))/mean(exp(Wv4)), sigmauSq4 = mean(exp(Wu4)),
          sigmavSq4 = mean(exp(Wv4)), sigma4 = sqrt(mean(exp(Wu4)) + mean(exp(Wv4))),
          lambda4 = sqrt(mean(exp(Wu4))/mean(exp(Wv4))), sigmau4 = sqrt(mean(exp(Wu4))),
          sigmav4 = sqrt(mean(exp(Wv4))), gamma4 = mean(exp(Wu4))/(mean(exp(Wu4)) +
          mean(exp(Wv4))))
      } else {
        if (object$nClasses == 5) {
          cRes <- c(cRes, sigmaSq3 = mean(exp(Wu3)) + mean(exp(Wv3)), lambdaSq3 = mean(exp(Wu3))/mean(exp(Wv3)),
          sigmauSq3 = mean(exp(Wu3)), sigmavSq3 = mean(exp(Wv3)), sigma3 = sqrt(mean(exp(Wu3)) +
            mean(exp(Wv3))), lambda3 = sqrt(mean(exp(Wu3))/mean(exp(Wv3))),
          sigmau3 = sqrt(mean(exp(Wu3))), sigmav3 = sqrt(mean(exp(Wv3))),
          gamma3 = mean(exp(Wu3))/(mean(exp(Wu3)) + mean(exp(Wv3))), sigmaSq4 = mean(exp(Wu4)) +
            mean(exp(Wv4)), lambdaSq4 = mean(exp(Wu4))/mean(exp(Wv4)),
          sigmauSq4 = mean(exp(Wu4)), sigmavSq4 = mean(exp(Wv4)), sigma4 = sqrt(mean(exp(Wu4)) +
            mean(exp(Wv4))), lambda4 = sqrt(mean(exp(Wu4))/mean(exp(Wv4))),
          sigmau4 = sqrt(mean(exp(Wu4))), sigmav4 = sqrt(mean(exp(Wv4))),
          gamma4 = mean(exp(Wu4))/(mean(exp(Wu4)) + mean(exp(Wv4))), sigmaSq5 = mean(exp(Wu5)) +
            mean(exp(Wv5)), lambdaSq5 = mean(exp(Wu5))/mean(exp(Wv5)),
          sigmauSq5 = mean(exp(Wu5)), sigmavSq5 = mean(exp(Wv5)), sigma5 = sqrt(mean(exp(Wu5)) +
            mean(exp(Wv5))), lambda5 = sqrt(mean(exp(Wu5))/mean(exp(Wv5))),
          sigmau5 = sqrt(mean(exp(Wu5))), sigmav5 = sqrt(mean(exp(Wv5))),
          gamma5 = mean(exp(Wu5))/(mean(exp(Wu5)) + mean(exp(Wv5))))
        }
      }
    }
  }
  return(cRes)
}

# coefficients from summary.sfalcmpanel ----------
#' @rdname coef
#' @aliases coef.summary.sfalcmpanel
#' @export
coef.summary.sfalcmpanel <- function(...) {
  coef.summary.sfacross(...)
}
