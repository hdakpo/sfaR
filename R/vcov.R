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

#' Compute variance-covariance matrix of stochastic frontier models
#'
#' \code{\link{vcov}} computes the variance-covariance matrix of the maximum
#' likelihood (ML) coefficients from stochastic frontier models estimated with 
#' \code{\link{sfacross}}, \code{\link{sfalcmcross}}, 
#' \code{\link{sfagzisfcross}}, \code{\link{sfacnsfcross}}, 
#' \code{\link{sfamisfcross}}, \code{\link{sfazisfcross}}, 
#' \code{\link{sfametacross}}, \code{\link{sfaselectioncross}}, 
#' \code{\link{sfapanel1}}, or \code{\link{sfalcmpanel}}.
#'
#' @details The variance-covariance matrix is obtained by the inversion of the 
#' negative Hessian matrix. Depending on the distribution and the
#' \code{'hessianType'} option, the analytical/numeric Hessian or the bhhh 
#' Hessian is evaluated. For the metafrontier, when the `modelType` is 'hhl14', 
#' a sandwich-form of the covariance matrix is computed.
#'
#' The argument \code{extraPar}, is currently available only for objects of 
#' class \code{'sfacross'}, \code{'sfametacross'}, and 
#' \code{'sfaselectioncross'}. When \code{'extraPar = TRUE'}, the 
#' variance-covariance of the additional parameters is obtained using the 
#' delta method.
#'
#' @param object A stochastic frontier model returned by \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, 
#' \code{\link{sfaselectioncross}}, \code{\link{sfapanel1}}, or 
#' \code{\link{sfalcmpanel}}.
#' @param extraPar Logical. Only available for non heteroscedastic models
#' returned by \code{\link{sfacross}}, \code{\link{sfametacross}}, 
#' and \code{\link{sfaselectioncross}}. Default = \code{FALSE}. If \code{TRUE}, 
#' variances and covariances of additional parameters are returned:
#'
#' \code{sigmaSq} = \code{sigmauSq} + \code{sigmavSq}
#'
#' \code{lambdaSq} = \code{sigmauSq}/\code{sigmavSq}
#'
#' \code{sigmauSq} = \eqn{\exp{(Wu)}} = \eqn{\exp{(\bm{\delta} \mathbf{Z}_u)}}
#'
#' \code{sigmavSq} = \eqn{\exp{(Wv)}} = \eqn{\exp{(\bm{\phi} \mathbf{Z}_v)}}
#'
#' \code{sigma} = \code{sigmaSq}^0.5
#'
#' \code{lambda} = \code{lambdaSq}^0.5
#'
#' \code{sigmau} = \code{sigmauSq}^0.5
#'
#' \code{sigmav} = \code{sigmavSq}^0.5
#'
#' \code{gamma} = \code{sigmauSq}/(\code{sigmauSq} + \code{sigmavSq}).
#' 
#' For object of class \code{'sfametacross'}, 
#' when \code{'modelType'} is 'bpo04a', 'bpo04b', 'aos17a', or 'aos17b', the
#' additional parameters are only returned for the group frontiers. In the case 
#' of 'aos17a', or 'aos17b', no parameters are estimated for the metafrontier.
#' 
#' @param ... Currently ignored
#'
#' @name vcov
#' 
#' @aliases vcov.sfacross vcov.sfalcmcross vcov.sfagzisfcross vcov.sfacnsfcross 
#' vcov.sfamisfcross vcov.sfazisfcross vcov.sfaselectioncross vcov.sfametacross
#' vcov.sfapanel1 vcov.sfalcmpanel
#'
#' @return The variance-covariance matrix of the maximum likelihood
#' coefficients is returned. For object of class \code{'sfametacross'}, a list
#' containing the variance-covariance matrices for each group, and in some cases
#' the metafrontier are returned.
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
#' model fitting function using cross-sectional data.
#' 
#' \code{\link{sfapanel1}}, for the first generation stochastic frontier 
#' analysis model fitting function using panel data.
#' 
#' \code{\link{sfalcmpanel}}, for the latent class stochastic frontier analysis
#' model fitting function using panel data.
#' 
#' @keywords methods vcov
#'
#' @examples
#'
#' # Using data on Spanish dairy farms
#' ## Cobb Douglas (production function) half normal distribution
#' cb_s_h <- sfacross(formula = YIT ~ X1 + X2 + X3 + X4, udist = 'hnormal',
#' data = dairyspain, S = 1, method = 'bfgs')
#' vcov(cb_s_h)
#' vcov(cb_s_h, extraPar = TRUE)
#'  
#' # Other variance-covariance matrices can be obtained using the 
#' # sandwich package
#'  
#' # Robust variance-covariance matrix
#'  
#' requireNamespace('sandwich', quietly = TRUE)
#'  
#' sandwich::vcovCL(cb_s_h)
#'  
#' # Coefficients and standard errors can be obtained using lmtest package
#'  
#' requireNamespace('lmtest', quietly = TRUE)
#'  
#' lmtest::coeftest(cb_s_h, vcov. = sandwich::vcovCL)
#'  
#' # Clustered standard errors
#'  
#' lmtest::coeftest(cb_s_h, vcov. = sandwich::vcovCL, cluster = ~ FARM)
#'  
#' # Doubly clustered standard errors
#'  
#' lmtest::coeftest(cb_s_h, vcov. = sandwich::vcovCL, cluster = ~ FARM + YEAR)
#'  
#' # BHHH standard errors
#'  
#' lmtest::coeftest(cb_s_h, vcov. = sandwich::vcovOPG)
#'  
#' # Adjusted BHHH standard errors
#'  
#' lmtest::coeftest(cb_s_h, vcov. = sandwich::vcovOPG, adjust = TRUE)
#'
#' # Using data on eighty-two countries production (GDP)
#' ## LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- sfalcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal',
#' data = worldprod, uhet = ~ initStat, S = 1)
#' vcov(cb_2c_h)
#'
#' @export
# @exportS3Method vcov sfacross
# variance covariance matrix for sfacross ----------
vcov.sfacross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1])) {
    stop("argument 'extraPar' must be a single logical value",
      call. = FALSE)
  }
  resCov <- object$invHessian
  if (extraPar) {
    if (object$udist %in% c("tnormal", "lognormal")) {
      if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar >
        1) {
        stop("argument 'extraPar' is not available for heteroscedasctic models",
          call. = FALSE)
      }
      delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
        1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
        object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
        object$nuZUvar + object$nvZVvar)]
      uHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 3)
      vHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 4)
    } else {
      if (object$nuZUvar > 1 || object$nvZVvar > 1) {
        stop("argument 'extraPar' is not available for heteroscedasctic models",
          call. = FALSE)
      }
      delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nuZUvar +
        1):(object$nXvar + object$nuZUvar + object$nvZVvar)]
      uHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 2)
      vHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 3)
    }
    Wu <- mean(as.numeric(crossprod(matrix(delta), t(uHvar))))
    Wv <- mean(as.numeric(crossprod(matrix(phi), t(vHvar))))
    jac <- diag(nrow(resCov))
    jac <- rbind(jac, matrix(0, nrow = 9, ncol = ncol(resCov)))
    rownames(jac) <- c(rownames(resCov), "sigmaSq", "lambdaSq",
      "sigmauSq", "sigmavSq", "sigma", "lambda", "sigmau",
      "sigmav", "gamma")
    colnames(jac) <- colnames(resCov)
    jac["sigmaSq", "Zu_(Intercept)"] <- exp(Wu)
    jac["sigmaSq", "Zv_(Intercept)"] <- exp(Wv)
    jac["lambdaSq", "Zu_(Intercept)"] <- exp(Wu)/exp(Wv)
    jac["lambdaSq", "Zv_(Intercept)"] <- -exp(Wu + Wv)/exp(2 *
      Wv)
    jac["sigmauSq", "Zu_(Intercept)"] <- exp(Wu)
    jac["sigmavSq", "Zv_(Intercept)"] <- exp(Wv)
    jac["sigma", "Zu_(Intercept)"] <- 1/2 * exp(Wu) * (exp(Wu) +
      exp(Wv))^(-1/2)
    jac["sigma", "Zv_(Intercept)"] <- 1/2 * exp(Wv) * (exp(Wu) +
      exp(Wv))^(-1/2)
    jac["lambda", "Zu_(Intercept)"] <- 1/2 * exp(Wu/2)/exp(Wv/2)
    jac["lambda", "Zv_(Intercept)"] <- -1/2 * exp(Wu/2 +
      Wv/2)/exp(Wv)
    jac["sigmau", "Zu_(Intercept)"] <- 1/2 * exp(Wu/2)
    jac["sigmav", "Zv_(Intercept)"] <- 1/2 * exp(Wv/2)
    jac["gamma", "Zu_(Intercept)"] <- (exp(Wu) * (exp(Wu) +
      exp(Wv)) - exp(2 * Wu))/(exp(Wu) + exp(Wv))^2
    jac["gamma", "Zv_(Intercept)"] <- -exp(Wu + Wv)/(exp(Wu) +
      exp(Wv))^2
    resCov <- jac %*% resCov %*% t(jac)
  }
  return(resCov)
}

# variance covariance matrix for sfalcmcross ----------
#' @rdname vcov
#' @export
# @exportS3Method vcov sfalcmcross
vcov.sfalcmcross <- function(object, ...) {
  resCov <- object$invHessian
  return(resCov)
}

# variance covariance matrix for sfagzisfcross ----------
#' @rdname vcov
#' @export
# @exportS3Method vcov sfagzisfcross
vcov.sfagzisfcross <- function(...) {
  vcov.sfalcmcross(...)
}

# variance covariance matrix for sfacnsfcross ----------
#' @rdname vcov
#' @export
# @exportS3Method vcov sfacnsfcross
vcov.sfacnsfcross <- function(...) {
  vcov.sfalcmcross(...)
}

# variance covariance matrix for sfamisfcross ----------
#' @rdname vcov
#' @export
# @exportS3Method vcov sfamisfcross
vcov.sfamisfcross <- function(...) {
  vcov.sfalcmcross(...)
}

# variance covariance matrix for sfazisfcross ----------
#' @rdname vcov
#' @export
# @exportS3Method vcov sfazisfcross
vcov.sfazisfcross <- function(...) {
  vcov.sfalcmcross(...)
}

# vcov matrix for sfaselectioncross ----------
#' @rdname vcov
#' @export
# @exportS3Method vcov sfaselectioncross
vcov.sfaselectioncross <- function(object, extraPar = FALSE,
  ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1])) {
    stop("argument 'extraPar' must be a single logical value",
      call. = FALSE)
  }
  resCov <- object$invHessian
  if (extraPar) {
    if (object$nuZUvar > 1 || object$nvZVvar > 1) {
      stop("argument 'extraPar' is not available for heteroscedasctic models",
        call. = FALSE)
    }
    delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
      object$nuZUvar)]
    phi <- object$mlParam[(object$nXvar + object$nuZUvar +
      1):(object$nXvar + object$nuZUvar + object$nvZVvar)]
    uHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable$selectDum ==
      1, ], rhs = 2)
    vHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable$selectDum ==
      1, ], rhs = 3)
    Wu <- mean(as.numeric(crossprod(matrix(delta), t(uHvar))))
    Wv <- mean(as.numeric(crossprod(matrix(phi), t(vHvar))))
    jac <- diag(nrow(resCov))
    jac <- rbind(jac, matrix(0, nrow = 9, ncol = ncol(resCov)))
    rownames(jac) <- c(rownames(resCov), "sigmaSq", "lambdaSq",
      "sigmauSq", "sigmavSq", "sigma", "lambda", "sigmau",
      "sigmav", "gamma")
    colnames(jac) <- colnames(resCov)
    jac["sigmaSq", "Zu_(Intercept)"] <- exp(Wu)
    jac["sigmaSq", "Zv_(Intercept)"] <- exp(Wv)
    jac["lambdaSq", "Zu_(Intercept)"] <- exp(Wu)/exp(Wv)
    jac["lambdaSq", "Zv_(Intercept)"] <- -exp(Wu + Wv)/exp(2 *
      Wv)
    jac["sigmauSq", "Zu_(Intercept)"] <- exp(Wu)
    jac["sigmavSq", "Zv_(Intercept)"] <- exp(Wv)
    jac["sigma", "Zu_(Intercept)"] <- 1/2 * exp(Wu) * (exp(Wu) +
      exp(Wv))^(-1/2)
    jac["sigma", "Zv_(Intercept)"] <- 1/2 * exp(Wv) * (exp(Wu) +
      exp(Wv))^(-1/2)
    jac["lambda", "Zu_(Intercept)"] <- 1/2 * exp(Wu/2)/exp(Wv/2)
    jac["lambda", "Zv_(Intercept)"] <- -1/2 * exp(Wu/2 +
      Wv/2)/exp(Wv)
    jac["sigmau", "Zu_(Intercept)"] <- 1/2 * exp(Wu/2)
    jac["sigmav", "Zv_(Intercept)"] <- 1/2 * exp(Wv/2)
    jac["gamma", "Zu_(Intercept)"] <- (exp(Wu) * (exp(Wu) +
      exp(Wv)) - exp(2 * Wu))/(exp(Wu) + exp(Wv))^2
    jac["gamma", "Zv_(Intercept)"] <- -exp(Wu + Wv)/(exp(Wu) +
      exp(Wv))^2
    resCov <- jac %*% resCov %*% t(jac)
  }
  return(resCov)
}

# vcov matrix for sfametacross ----------
#' @rdname vcov
#' @export
# @exportS3Method vcov sfametacross
vcov.sfametacross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1])) {
    stop("argument 'extraPar' must be a single logical value",
      call. = FALSE)
  }
  resCov <- object$invHessian
  if (extraPar) {
    group_var <- object$dataTable[[object$Ngroup + 1]][object$name_meta_var][,
      1]
    group_var_list <- sort(unique(group_var))
    for (g in group_var_list) {
      if (object$udist %in% c("tnormal", "lognormal")) {
        if (object$nuZUvar > 1 || object$nvZVvar > 1 ||
          object$nmuZUvar > 1) {
          stop("argument 'extraPar' is not available for heteroscedasctic models",
          call. = FALSE)
        }
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar),
          which(group_var_list == g)]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
          object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
          object$nuZUvar + object$nvZVvar), which(group_var_list ==
          g)]
        uHvar <- model.matrix(object$formula, data = object$dataTable[[which(group_var_list ==
          g)]], rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable[[which(group_var_list ==
          g)]], rhs = 4)
        Wu <- mean(as.numeric(crossprod(matrix(delta),
          t(uHvar))))
        Wv <- mean(as.numeric(crossprod(matrix(phi),
          t(vHvar))))
      } else {
        if (object$nuZUvar > 1 || object$nvZVvar > 1) {
          stop("argument 'extraPar' is not available for heteroscedasctic models",
          call. = FALSE)
        }
        delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
          object$nuZUvar), which(group_var_list == g)]
        phi <- object$mlParam[(object$nXvar + object$nuZUvar +
          1):(object$nXvar + object$nuZUvar + object$nvZVvar),
          which(group_var_list == g)]
        uHvar <- model.matrix(object$formula, data = object$dataTable[[which(group_var_list ==
          g)]], rhs = 2)
        vHvar <- model.matrix(object$formula, data = object$dataTable[[which(group_var_list ==
          g)]], rhs = 3)
        Wu <- mean(as.numeric(crossprod(matrix(delta),
          t(uHvar))))
        Wv <- mean(as.numeric(crossprod(matrix(phi),
          t(vHvar))))
      }
      jac <- diag(nrow(resCov[[which(group_var_list ==
        g)]]))
      jac <- rbind(jac, matrix(0, nrow = 9, ncol = ncol(resCov[[which(group_var_list ==
        g)]])))
      rownames(jac) <- c(rownames(resCov[[which(group_var_list ==
        g)]]), "sigmaSq", "lambdaSq", "sigmauSq", "sigmavSq",
        "sigma", "lambda", "sigmau", "sigmav", "gamma")
      colnames(jac) <- colnames(resCov[[which(group_var_list ==
        g)]])
      jac["sigmaSq", "Zu_(Intercept)"] <- exp(Wu)
      jac["sigmaSq", "Zv_(Intercept)"] <- exp(Wv)
      jac["lambdaSq", "Zu_(Intercept)"] <- (exp(Wu)/exp(Wv))
      jac["lambdaSq", "Zv_(Intercept)"] <- (-exp(Wu + Wv)/exp(2 *
        Wv))
      jac["sigmauSq", "Zu_(Intercept)"] <- exp(Wu)
      jac["sigmavSq", "Zv_(Intercept)"] <- exp(Wv)
      jac["sigma", "Zu_(Intercept)"] <- (1/2 * exp(Wu) *
        (exp(Wu) + exp(Wv))^(-1/2))
      jac["sigma", "Zv_(Intercept)"] <- (1/2 * exp(Wv) *
        (exp(Wu) + exp(Wv))^(-1/2))
      jac["lambda", "Zu_(Intercept)"] <- (1/2 * exp(Wu/2)/exp(Wv/2))
      jac["lambda", "Zv_(Intercept)"] <- (-1/2 * exp(Wu/2 +
        Wv/2)/exp(Wv))
      jac["sigmau", "Zu_(Intercept)"] <- (1/2 * exp(Wu/2))
      jac["sigmav", "Zv_(Intercept)"] <- (1/2 * exp(Wv/2))
      jac["gamma", "Zu_(Intercept)"] <- ((exp(Wu) * (exp(Wu) +
        exp(Wv)) - exp(2 * Wu))/(exp(Wu) + exp(Wv))^2)
      jac["gamma", "Zv_(Intercept)"] <- (-exp(Wu + Wv)/(exp(Wu) +
        exp(Wv))^2)
      resCov[[which(group_var_list == g)]] <- jac %*% resCov[[which(group_var_list ==
        g)]] %*% t(jac)
    }
    if (object$modelType %in% c("hhl14")) {
      if (object$udist %in% c("tnormal", "lognormal")) {
        if (object$nuZUvar > 1 || object$nvZVvar > 1 ||
          object$nmuZUvar > 1) {
          stop("argument 'extraPar' is not available for heteroscedasctic models",
          call. = FALSE)
        }
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar),
          object$Ngroup + 1]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
          object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
          object$nuZUvar + object$nvZVvar), object$Ngroup +
          1]
        uHvar <- model.matrix(object$formula, data = object$dataTable[[object$Ngroup +
          1]], rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable[[object$Ngroup +
          1]], rhs = 4)
        Wu <- mean(as.numeric(crossprod(matrix(delta),
          t(uHvar))))
        Wv <- mean(as.numeric(crossprod(matrix(phi),
          t(vHvar))))
      } else {
        if (object$nuZUvar > 1 || object$nvZVvar > 1) {
          stop("argument 'extraPar' is not available for heteroscedasctic models",
          call. = FALSE)
        }
        delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
          object$nuZUvar), object$Ngroup + 1]
        phi <- object$mlParam[(object$nXvar + object$nuZUvar +
          1):(object$nXvar + object$nuZUvar + object$nvZVvar),
          object$Ngroup + 1]
        uHvar <- model.matrix(object$formula, data = object$dataTable[[object$Ngroup +
          1]], rhs = 2)
        vHvar <- model.matrix(object$formula, data = object$dataTable[[object$Ngroup +
          1]], rhs = 3)
        Wu <- mean(as.numeric(crossprod(matrix(delta),
          t(uHvar))))
        Wv <- mean(as.numeric(crossprod(matrix(phi),
          t(vHvar))))
      }
      jac <- diag(nrow(resCov[[object$Ngroup + 1]]))
      jac <- rbind(jac, matrix(0, nrow = 9, ncol = ncol(resCov[[object$Ngroup +
        1]])))
      rownames(jac) <- c(rownames(resCov[[object$Ngroup +
        1]]), "sigmaSq", "lambdaSq", "sigmauSq", "sigmavSq",
        "sigma", "lambda", "sigmau", "sigmav", "gamma")
      colnames(jac) <- colnames(resCov[[object$Ngroup +
        1]])
      jac["sigmaSq", "Zu_(Intercept)"] <- exp(Wu)
      jac["sigmaSq", "Zv_(Intercept)"] <- exp(Wv)
      jac["lambdaSq", "Zu_(Intercept)"] <- (exp(Wu)/exp(Wv))
      jac["lambdaSq", "Zv_(Intercept)"] <- (-exp(Wu + Wv)/exp(2 *
        Wv))
      jac["sigmauSq", "Zu_(Intercept)"] <- exp(Wu)
      jac["sigmavSq", "Zv_(Intercept)"] <- exp(Wv)
      jac["sigma", "Zu_(Intercept)"] <- (1/2 * exp(Wu) *
        (exp(Wu) + exp(Wv))^(-1/2))
      jac["sigma", "Zv_(Intercept)"] <- (1/2 * exp(Wv) *
        (exp(Wu) + exp(Wv))^(-1/2))
      jac["lambda", "Zu_(Intercept)"] <- (1/2 * exp(Wu/2)/exp(Wv/2))
      jac["lambda", "Zv_(Intercept)"] <- (-1/2 * exp(Wu/2 +
        Wv/2)/exp(Wv))
      jac["sigmau", "Zu_(Intercept)"] <- (1/2 * exp(Wu/2))
      jac["sigmav", "Zv_(Intercept)"] <- (1/2 * exp(Wv/2))
      jac["gamma", "Zu_(Intercept)"] <- ((exp(Wu) * (exp(Wu) +
        exp(Wv)) - exp(2 * Wu))/(exp(Wu) + exp(Wv))^2)
      jac["gamma", "Zv_(Intercept)"] <- (-exp(Wu + Wv)/(exp(Wu) +
        exp(Wv))^2)
      resCov[[object$Ngroup + 1]] <- jac %*% resCov[[object$Ngroup +
        1]] %*% t(jac)
    }
  }
  return(resCov)
}

# variance covariance matrix for sfapanel1 ----------
#' @rdname vcov
#' @aliases vcov.sfapanel1
#' @export
vcov.sfapanel1 <- function(...) {
  vcov.sfalcmcross(...)
}

# variance covariance matrix for sfalcmpanel ----------
#' @rdname vcov
#' @aliases vcov.sfalcmpanel
#' @export
vcov.sfalcmpanel <- function(...) {
  vcov.sfalcmcross(...)
}
