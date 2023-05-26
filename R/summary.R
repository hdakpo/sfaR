################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Summary of optimization objects                                              #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Sample selection correction for Stochastic Frontier Model         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Summary of results for stochastic frontier models
#'
#' Create and print summary results for stochastic frontier models returned by 
#' \code{\link{sfacross}}, \code{\link{sfalcmcross}}, or 
#' \code{\link{sfaselectioncross}}.
#'
#' @param object An object of either class \code{'sfacross'} returned by the
#' function \code{\link{sfacross}}, or \code{'sfalcmcross'} returned by the
#' function \code{\link{sfalcmcross}}, or class \code{'sfaselectioncross'} returned 
#' by the function \code{\link{sfaselectioncross}}.
#' @param grad Logical. Default = \code{FALSE}. If \code{TRUE}, the gradient
#' for the maximum likelihood (ML) estimates of the different parameters is
#' returned.
#' @param ci Logical. Default = \code{FALSE}. If \code{TRUE}, the 95%
#' confidence interval for the different parameters (OLS or/and ML estimates) is
#' returned.
#' @param ... Currently ignored.
#' @param x An object of either class \code{'summary.sfacross'}, \code{'summary.sfalcmcross'}, or\cr
#' \code{'summary.sfaselectioncross'}.
#' @param digits Numeric. Number of digits displayed in values.
#'
#' @name summary
#'
#' @return The \code{\link{summary}} method returns a list of class
#' \code{'summary.sfacross'}, \code{'summary.sfalcmcross'}, or\cr
#' \code{'summary.sfaselectioncross'} 
#' that contains the same elements as an object returned by \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, or \code{\link{sfaselectioncross}} with the 
#' following additional elements:
#'
#' \item{AIC}{Akaike information criterion.}
#'
#' \item{BIC}{Bayesian information criterion.}
#'
#' \item{HQIC}{Hannan-Quinn information criterion.}
#'
#' \item{sigmavSq}{For \code{object} of class \code{'sfacross'} or 
#' \code{'sfaselectioncross'}. Variance of
#' the two-sided error term (\eqn{\sigma_v^2}).}
#'
#' \item{sigmauSq}{For \code{object} of class \code{'sfacross'} or 
#' \code{'sfaselectioncross'}. Parametrization of the variance of the one-sided
#'  error term (\eqn{\sigma_u^2}).}
#'
#' \item{Varu}{For \code{object} of class \code{'sfacross'} or 
#' \code{'sfaselectioncross'}. Variance of the one-sided error term.}
#'
#' \item{theta}{For \code{object} of class \code{'sfacross'} with \code{'udist
#' = uniform'}.  \eqn{\Theta} value in the case the uniform distribution is
#' defined as: \eqn{u_i \in [0, \Theta]}.}
#'
#' \item{Eu}{For \code{object} of class \code{'sfacross'} or 
#' \code{'sfaselectioncross'}. Expected unconditional inefficiency
#' (\eqn{E[u]}).}
#'
#' \item{Expu}{For \code{object} of class \code{'sfacross'} or 
#' \code{'sfaselectioncross'}. Expected unconditional efficiency 
#' (\eqn{E[\exp(u)]}).}
#'
#' \item{olsRes}{For \code{object} of class \code{'sfacross'}. Matrix of OLS
#' estimates, their standard errors, t-values, P-values, and when \code{ci =
#' TRUE} their confidence intervals.}
#' 
#' \item{ols2StepRes}{For \code{object} of class \code{'sfaselectioncross'}. 
#' Matrix of OLS 2 step estimates, their standard errors, t-values, P-values,
#' and when \code{ci = TRUE} their confidence intervals.}
#'
#' \item{mlRes}{Matrix of ML estimates, their standard errors, z-values,
#' asymptotic P-values, and when \code{grad = TRUE} their gradient, \code{ci =
#' TRUE} their confidence intervals.}
#'
#' \item{chisq}{For \code{object} of class \code{'sfacross'}. Chi-square
#' statistics of the difference between the stochastic frontier and the OLS.}
#'
#' \item{df}{Degree of freedom for the inefficiency model.}
#' 
# @author K Hervé Dakpo
#' 
#' @seealso \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function for cross-sectional or pooled data.
#' 
#' \code{\link{sfalcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function for cross-sectional or pooled data.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function for cross-sectional or pooled data.
#' 
#' \code{\link[=print.sfacross]{print}} for printing \code{sfacross} object.
#'
#' \code{\link[=coef.sfacross]{coef}} for extracting coefficients of the
#' estimation.
#'
#' \code{\link[=efficiencies.sfacross]{efficiencies}} for computing
#' (in-)efficiency estimates.
#'
#' \code{\link[=fitted.sfacross]{fitted}} for extracting the fitted frontier
#' values.
#'
#' \code{\link[=ic.sfacross]{ic}} for extracting information criteria.
#'
#' \code{\link[=logLik.sfacross]{logLik}} for extracting log-likelihood
#' value(s) of the estimation.
#'
#' \code{\link[=marginal.sfacross]{marginal}} for computing marginal effects of
#' inefficiency drivers.
#'
#' \code{\link[=residuals.sfacross]{residuals}} for extracting residuals of the
#' estimation.
#'
#' \code{\link[=vcov.sfacross]{vcov}} for computing the variance-covariance
#' matrix of the coefficients.
#' 
#'  \code{\link[=bread.sfacross]{bread}} for bread for sandwich estimator.
#' 
#' \code{\link[=estfun.sfacross]{estfun}} for gradient extraction for each 
#' observation.
#'
#' \code{\link{skewnessTest}} for implementing skewness test.
#'
#' @keywords methods summary
#'
#' @examples
#'
#' ## Using data on fossil fuel fired steam electric power generation plants in the U.S.
#' # Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = 'mla')
#' summary(tl_u_ts, grad = TRUE, ci = TRUE)
#'
#' @aliases summary.sfacross
#' @export
# summary for sfacross ----------
summary.sfacross <- function(object, grad = FALSE, ci = FALSE,
  ...) {
  if (length(grad) != 1 || !is.logical(grad[1])) {
    stop("argument 'grad' must be a single logical value",
      call. = FALSE)
  }
  if (length(ci) != 1 || !is.logical(ci[1])) {
    stop("argument 'ci' must be a single logical value",
      call. = FALSE)
  }
  object$AIC <- -2 * object$mlLoglik + 2 * object$nParm
  object$BIC <- -2 * object$mlLoglik + log(object$Nobs) * object$nParm
  object$HQIC <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) *
    object$nParm
  if (object$udist == "tnormal") {
    if (object$scaling) {
      muHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 2)
      uHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 3)
      vHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 4)
      delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        (object$nuZUvar - 1))]
      tau <- object$mlParam[object$nXvar + (object$nuZUvar -
        1) + 1]
      cu <- object$mlParam[object$nXvar + (object$nuZUvar -
        1) + 2]
      phi <- object$mlParam[(object$nXvar + (object$nuZUvar -
        1) + 2 + 1):(object$nXvar + (object$nuZUvar -
        1) + 2 + object$nvZVvar)]
    } else {
      omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        object$nmuZUvar)]
      delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
        1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
        object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
        object$nuZUvar + object$nvZVvar)]
      muHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 2)
      uHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 3)
      vHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 4)
    }
  } else {
    if (object$udist == "lognormal") {
      omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        object$nmuZUvar)]
      delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
        1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
        object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
        object$nuZUvar + object$nvZVvar)]
      muHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 2)
      uHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 3)
      vHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 4)
    } else {
      delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nuZUvar +
        1):(object$nXvar + object$nuZUvar + object$nvZVvar)]
      uHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 2)
      vHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 3)
    }
  }
  mu <- if (object$udist == "tnormal") {
    if (object$scaling) {
      exp(as.numeric(crossprod(matrix(delta), t(uHvar[,
        -1, drop = FALSE])))) * tau
    } else {
      as.numeric(crossprod(matrix(omega), t(muHvar)))
    }
  } else {
    if (object$udist == "lognormal") {
      as.numeric(crossprod(matrix(omega), t(muHvar)))
    } else {
      NULL
    }
  }
  P <- if (object$udist == "gamma") {
    object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
      1]
  } else {
    NULL
  }
  k <- if (object$udist == "weibull") {
    object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
      1]
  } else {
    NULL
  }
  lambda <- if (object$udist == "tslaplace") {
    object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
      1]
  } else {
    NULL
  }
  Wu <- if (object$udist == "tnormal" & object$scaling == TRUE) {
    cu + 2 * as.numeric(crossprod(matrix(delta), t(uHvar[,
      -1, drop = FALSE])))
  } else {
    as.numeric(crossprod(matrix(delta), t(uHvar)))
  }
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  object$sigmavSq <- mean(exp(Wv))
  object$sigmauSq <- mean(exp(Wu))
  if (object$udist == "uniform") {
    object$theta <- sqrt(12 * object$sigmauSq)
  }
  object$Varu <- varuFun(object = object, mu = mu, P = P, k = k,
    lambda = lambda)
  object$Eu <- euFun(object = object, mu = mu, P = P, k = k,
    lambda = lambda)
  if (object$logDepVar == TRUE) {
    object$Expu <- eExpuFun(object = object, mu = mu, P = P,
      k = k, lambda = lambda)
  }
  # OLS estimates and stder, p-values, CI
  dfOLS <- object$Nobs - object$nXvar
  if (ci) {
    olsRes <- matrix(nrow = object$nXvar + 1, ncol = 6)
    colnames(olsRes) <- c("Coefficient", "Std. Error", "binf",
      "bsup", "t value", "Pr(>|t|)")
    olsRes[, 1] <- c(object$olsParam, object$olsSigmasq)
    olsRes[, 2] <- c(object$olsStder, NA)
    olsRes[, 3] <- olsRes[, 1] - qt(0.975, df = dfOLS) *
      olsRes[, 2]
    olsRes[, 4] <- olsRes[, 1] + qt(0.975, df = dfOLS) *
      olsRes[, 2]
    olsRes[, 5] <- olsRes[, 1]/olsRes[, 2]
    olsRes[, 6] <- 2 * pt(-abs(olsRes[, 5]), df = dfOLS)
  } else {
    olsRes <- matrix(nrow = object$nXvar + 1, ncol = 4)
    colnames(olsRes) <- c("Coefficient", "Std. Error", "t value",
      "Pr(>|t|)")
    olsRes[, 1] <- c(object$olsParam, object$olsSigmasq)
    olsRes[, 2] <- c(object$olsStder, NA)
    olsRes[, 3] <- olsRes[, 1]/olsRes[, 2]
    olsRes[, 4] <- 2 * pt(-abs(olsRes[, 3]), df = dfOLS)
  }
  row.names(olsRes) <- c(names(object$olsParam), "sigmaSq")
  object$olsRes <- olsRes
  # MLE estimates and stder, p-values, CI, Gradient
  if (grad && ci) {
    mlRes <- matrix(nrow = object$nParm, ncol = 7)
    colnames(mlRes) <- c("Coefficient", "Std. Error", "binf",
      "bsup", "gradient", "z-value", "Pr(>|z|)")
    mlRes[, 1] <- object$mlParam
    mlRes[, 2] <- sqrt(diag(object$invHessian))
    mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[, 2]
    mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[, 2]
    mlRes[, 5] <- object$gradient
    mlRes[, 6] <- mlRes[, 1]/mlRes[, 2]
    mlRes[, 7] <- 2 * pnorm(-abs(mlRes[, 6]))
  } else {
    if (grad == TRUE && ci == FALSE) {
      mlRes <- matrix(nrow = object$nParm, ncol = 5)
      colnames(mlRes) <- c("Coefficient", "Std. Error",
        "gradient", "z-value", "Pr(>|z|)")
      mlRes[, 1] <- object$mlParam
      mlRes[, 2] <- sqrt(diag(object$invHessian))
      mlRes[, 3] <- object$gradient
      mlRes[, 4] <- mlRes[, 1]/mlRes[, 2]
      mlRes[, 5] <- 2 * pnorm(-abs(mlRes[, 4]))
    } else {
      if (grad == FALSE && ci == TRUE) {
        mlRes <- matrix(nrow = object$nParm, ncol = 6)
        colnames(mlRes) <- c("Coefficient", "Std. Error",
          "binf", "bsup", "z-value", "Pr(>|z|)")
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[,
          2]
        mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[,
          2]
        mlRes[, 5] <- mlRes[, 1]/mlRes[, 2]
        mlRes[, 6] <- 2 * pnorm(-abs(mlRes[, 5]))
      } else {
        mlRes <- matrix(nrow = object$nParm, ncol = 4)
        colnames(mlRes) <- c("Coefficient", "Std. Error",
          "z value", "Pr(>|z|)")
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1]/mlRes[, 2]
        mlRes[, 4] <- 2 * pnorm(-abs(mlRes[, 3]))
      }
    }
  }
  row.names(mlRes) <- names(object$startVal)
  object$mlRes <- mlRes
  object$chisq <- 2 * (object$mlLoglik - object$olsLoglik)
  object$df <- object$nParm - object$nXvar - object$nvZVvar
  class(object) <- "summary.sfacross"
  return(object)
}

# print summary for sfacross ----------
#' @rdname summary
#' @aliases print.summary.sfacross
#' @export
print.summary.sfacross <- function(x, digits = max(3, getOption("digits") -
  2), ...) {
  mlRes <- x$mlRes
  if (dim(mlRes)[2] == 4) {
    mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1], digits = digits,
      format = "f"))
    mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2], digits = digits,
      format = "f"))
    mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3], digits = digits,
      format = "f"))
    mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4], digits = digits,
      format = "e"))
  } else {
    if (dim(mlRes)[2] == 5) {
      mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1], digits = digits,
        format = "f"))
      mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2], digits = digits,
        format = "f"))
      mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3], digits = digits,
        format = "e"))
      mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4], digits = digits,
        format = "f"))
      mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5], digits = digits,
        format = "e"))
    } else {
      if (dim(mlRes)[2] == 6) {
        mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
          digits = digits, format = "f"))
        mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
          digits = digits, format = "f"))
        mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
          digits = digits, format = "f"))
        mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
          digits = digits, format = "f"))
        mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5],
          digits = digits, format = "f"))
        mlRes[, 6] <- as.numeric(formatC(x$mlRes[, 6],
          digits = digits, format = "e"))
      } else {
        if (dim(mlRes)[2] == 7) {
          mlRes[, 1] <- as.numeric(formatC(x$mlRes[,
          1], digits = digits, format = "f"))
          mlRes[, 2] <- as.numeric(formatC(x$mlRes[,
          2], digits = digits, format = "f"))
          mlRes[, 3] <- as.numeric(formatC(x$mlRes[,
          3], digits = digits, format = "f"))
          mlRes[, 4] <- as.numeric(formatC(x$mlRes[,
          4], digits = digits, format = "f"))
          mlRes[, 5] <- as.numeric(formatC(x$mlRes[,
          5], digits = digits, format = "e"))
          mlRes[, 6] <- as.numeric(formatC(x$mlRes[,
          6], digits = digits, format = "f"))
          mlRes[, 7] <- as.numeric(formatC(x$mlRes[,
          7], digits = digits, format = "e"))
        }
      }
    }
  }
  row.names(mlRes) <- formatC(row.names(mlRes), width = max(nchar(row.names(mlRes))),
    flag = "-")
  mlRes1 <- mlRes[1:x$nXvar, , drop = FALSE]
  if (x$udist == "tnormal") {
    if (x$scaling) {
      uHvar <- model.matrix(x$formula, data = x$dataTable,
        rhs = 3)
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + (x$nuZUvar -
        1)), , drop = FALSE]
      mlRes3 <- mlRes[x$nXvar + (x$nuZUvar - 1) + 1, ,
        drop = FALSE]
      mlRes4 <- mlRes[x$nXvar + (x$nuZUvar - 1) + 2, ,
        drop = FALSE]
      mlRes5 <- mlRes[(x$nXvar + (x$nuZUvar - 1) + 2 +
        1):(x$nXvar + (x$nuZUvar - 1) + 2 + x$nvZVvar),
        , drop = FALSE]
    } else {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nmuZUvar),
        , drop = FALSE]
      mlRes3 <- mlRes[(x$nXvar + x$nmuZUvar + 1):(x$nXvar +
        x$nmuZUvar + x$nuZUvar), , drop = FALSE]
      mlRes4 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
        1):(x$nXvar + x$nmuZUvar + x$nuZUvar + x$nvZVvar),
        , drop = FALSE]
    }
  } else {
    if (x$udist == "lognormal") {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nmuZUvar),
        , drop = FALSE]
      mlRes3 <- mlRes[(x$nXvar + x$nmuZUvar + 1):(x$nXvar +
        x$nmuZUvar + x$nuZUvar), , drop = FALSE]
      mlRes4 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
        1):(x$nXvar + x$nmuZUvar + x$nuZUvar + x$nvZVvar),
        , drop = FALSE]
    } else {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar),
        , drop = FALSE]
      mlRes3 <- mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar +
        x$nuZUvar + x$nvZVvar), , drop = FALSE]
      if (x$udist %in% c("gamma", "weibull", "tslaplace")) {
        mlRes4 <- mlRes[x$nXvar + x$nuZUvar + x$nvZVvar +
          1, , drop = FALSE]
      }
    }
  }
  lengthSum <- 80
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(sfacrossdist(x$udist), "\n")
  cat("Dependent Variable:", paste0(rep(" ", lengthSum - 2 -
    nchar("Dependent Variable:") - nchar(paste0(attr(x$formula,
    "lhs")))), collapse = ""), paste0(attr(x$formula, "lhs")),
    "\n")
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    2 - nchar("Log likelihood solver:") - nchar(x$optType)),
    collapse = ""), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - 2 -
    nchar("Log likelihood iter:") - nchar(x$nIter)), collapse = ""),
    x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    2 - nchar("Log likelihood value:") - nchar(if (x$mlLoglik >
    1e+12) formatC(x$mlLoglik, digits = digits, format = "e") else formatC(x$mlLoglik,
    digits = digits, format = "f"))), collapse = ""), if (x$mlLoglik >
    1e+12)
    formatC(x$mlLoglik, digits = digits, format = "e") else formatC(x$mlLoglik, digits = digits, format = "f"),
    "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    2 - nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
    digits = digits, format = "e"))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"), "\n")
  cat("Estimation based on:", paste0(rep(" ", lengthSum - 2 -
    nchar("Estimation based on:") - nchar(x$Nobs) - nchar(x$nParm) -
    nchar("N = ") - nchar("and K = ") - 3), collapse = ""),
    "N = ", x$Nobs, "and K = ", x$nParm, "\n")
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - 2 - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(if (abs(x$AIC) > 1e+12) formatC(x$AIC,
    digits = 1, format = "e") else formatC(x$AIC, digits = 1,
    format = "f")) - nchar("AIC/N  = ") - nchar(if (abs(x$AIC/x$Nobs) >
    1e+12) formatC(x$AIC/x$Nobs, digits = 3, format = "e") else formatC(x$AIC/x$Nobs,
    digits = 3, format = "f")) - 3), collapse = ""), "AIC  = ",
    if (abs(x$AIC) > 1e+12)
      formatC(x$AIC, digits = 1, format = "e") else formatC(x$AIC, digits = 1, format = "f"), "AIC/N  = ",
    if (abs(x$AIC/x$Nobs) > 1e+12)
      formatC(x$AIC/x$Nobs, digits = 3, format = "e") else formatC(x$AIC/x$Nobs, digits = 3, format = "f"),
    "\n")
  cat(paste0(rep(" ", lengthSum - 2 - nchar("BIC  = ") - nchar(if (abs(x$BIC) >
    1e+12) formatC(x$BIC, digits = 1, format = "e") else formatC(x$BIC,
    digits = 1, format = "f")) - nchar("BIC/N  = ") - nchar(if (abs(x$BIC/x$Nobs) >
    1e+12) formatC(x$BIC/x$Nobs, digits = 3, format = "e") else formatC(x$BIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "BIC  = ",
    if (abs(x$BIC) > 1e+12)
      formatC(x$BIC, digits = 1, format = "e") else formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    if (abs(x$BIC/x$Nobs) > 1e+12)
      formatC(x$BIC/x$Nobs, digits = 3, format = "e") else formatC(x$BIC/x$Nobs, digits = 3, format = "f"),
    "\n")
  cat(paste0(rep(" ", lengthSum - 2 - nchar("HQIC = ") - nchar(if (abs(x$HQIC) >
    1e+12) formatC(x$HQIC, digits = 1, format = "e") else formatC(x$HQIC,
    digits = 1, format = "f")) - nchar("HQIC/N = ") - nchar(if (abs(x$HQIC/x$Nobs) >
    1e+12) formatC(x$HQIC/x$Nobs, digits = 3, format = "e") else formatC(x$HQIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "HQIC = ",
    if (abs(x$HQIC) > 1e+12)
      formatC(x$HQIC, digits = 1, format = "e") else formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    if (abs(x$HQIC/x$Nobs) > 1e+12)
      formatC(x$HQIC/x$Nobs, digits = 3, format = "e") else formatC(x$HQIC/x$Nobs, digits = 3, format = "f"),
    "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat("Variances: Sigma-squared(v)   = ", paste0(rep(" ", lengthSum -
    2 - nchar("Variances: Sigma-squared(v)   = ") - nchar(if (x$sigmavSq >
    1e+12) formatC(x$sigmavSq, digits = digits, format = "e") else formatC(x$sigmavSq,
    digits = digits, format = "f"))), collapse = ""), if (x$sigmavSq >
    1e+12)
    formatC(x$sigmavSq, digits = digits, format = "e") else formatC(x$sigmavSq, digits = digits, format = "f"),
    "\n")
  cat("           Sigma(v)           = ", paste0(rep(" ", lengthSum -
    2 - nchar("Variances: Sigma-squared(v)   = ") - nchar(if (x$sigmavSq >
    1e+12) formatC(x$sigmavSq, digits = digits, format = "e") else formatC(x$sigmavSq,
    digits = digits, format = "f"))), collapse = ""), if (x$sigmavSq >
    1e+12)
    formatC(x$sigmavSq, digits = digits, format = "e") else formatC(x$sigmavSq, digits = digits, format = "f"),
    "\n")
  cat("           Sigma-squared(u)   = ", paste0(rep(" ", lengthSum -
    2 - nchar("Variances: Sigma-squared(u)   = ") - nchar(if (x$sigmauSq >
    1e+12) formatC(x$sigmauSq, digits = digits, format = "e") else formatC(x$sigmauSq,
    digits = digits, format = "f"))), collapse = ""), if (x$sigmauSq >
    1e+12)
    formatC(x$sigmauSq, digits = digits, format = "e") else formatC(x$sigmauSq, digits = digits, format = "f"),
    "\n")
  cat("           Sigma(u)           = ", paste0(rep(" ", lengthSum -
    2 - nchar("Variances: Sigma-squared(u)   = ") - nchar(if (x$sigmauSq >
    1e+12) formatC(x$sigmauSq, digits = digits, format = "e") else formatC(x$sigmauSq,
    digits = digits, format = "f"))), collapse = ""), if (x$sigmauSq >
    1e+12)
    formatC(x$sigmauSq, digits = digits, format = "e") else formatC(x$sigmauSq, digits = digits, format = "f"),
    "\n")
  cat("Sigma = Sqrt[(s^2(u)+s^2(v))] = ", paste0(rep(" ", lengthSum -
    2 - nchar("Sigma = Sqrt[(s^2(u)+s^2(v))] = ") - nchar(if (sqrt(x$sigmavSq +
    x$sigmauSq) > 1e+12) formatC(sqrt(x$sigmavSq + x$sigmauSq),
    digits = digits, format = "e") else formatC(sqrt(x$sigmavSq +
    x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    if (sqrt(x$sigmavSq + x$sigmauSq) > 1e+12)
      formatC(sqrt(x$sigmavSq + x$sigmauSq), digits = digits,
        format = "e") else formatC(sqrt(x$sigmavSq + x$sigmauSq), digits = digits,
      format = "f"), "\n")
  cat("Gamma = sigma(u)^2/sigma^2    = ", paste0(rep(" ", lengthSum -
    2 - nchar("Gamma = sigma(u)^2/sigma^2    = ") - nchar(if (x$sigmauSq/(x$sigmavSq +
    x$sigmauSq) > 1e+12 || is.na(x$sigmauSq/(x$sigmavSq +
    x$sigmauSq))) formatC(x$sigmauSq/(x$sigmavSq + x$sigmauSq),
    digits = digits, format = "e") else formatC(x$sigmauSq/(x$sigmavSq +
    x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    if (x$sigmauSq/(x$sigmavSq + x$sigmauSq) > 1e+12 || is.na(x$sigmauSq/(x$sigmavSq +
      x$sigmauSq)))
      formatC(x$sigmauSq/(x$sigmavSq + x$sigmauSq), digits = digits,
        format = "e") else formatC(x$sigmauSq/(x$sigmavSq + x$sigmauSq), digits = digits,
      format = "f"), "\n")
  cat("Lambda = sigma(u)/sigma(v)    = ", paste0(rep(" ", lengthSum -
    2 - nchar("Lambda = sigma(u)/sigma(v)    = ") - nchar(if (sqrt(x$sigmauSq/x$sigmavSq) >
    1e+12) formatC(sqrt(x$sigmauSq/x$sigmavSq), digits = digits,
    format = "e") else formatC(sqrt(x$sigmauSq/x$sigmavSq),
    digits = digits, format = "f"))), collapse = ""), if (sqrt(x$sigmauSq/x$sigmavSq) >
    1e+12)
    formatC(sqrt(x$sigmauSq/x$sigmavSq), digits = digits,
      format = "e") else formatC(sqrt(x$sigmauSq/x$sigmavSq), digits = digits,
    format = "f"), "\n")
  cat("Var[u]/{Var[u]+Var[v]}        = ", paste0(rep(" ", lengthSum -
    2 - nchar("Var[u]/{Var[u]+Var[v]}        = ") - nchar(if (x$Varu/(x$Varu +
    x$sigmavSq) > 1e+12 || is.na(x$Varu/(x$Varu + x$sigmavSq))) formatC(x$Varu/(x$Varu +
    x$sigmavSq), digits = digits, format = "e") else formatC(x$Varu/(x$Varu +
    x$sigmavSq), digits = digits, format = "f"))), collapse = ""),
    if (x$Varu/(x$Varu + x$sigmavSq) > 1e+12 || is.na(x$Varu/(x$Varu +
      x$sigmavSq)))
      formatC(x$Varu/(x$Varu + x$sigmavSq), digits = digits,
        format = "e") else formatC(x$Varu/(x$Varu + x$sigmavSq), digits = digits,
      format = "f"), "\n")
  if (x$udist == "uniform") {
    if (x$nuZUvar > 1) {
      cat("theta averaged over observations \n")
    }
    cat("theta                         = ", paste0(rep(" ",
      lengthSum - 2 - nchar("theta                         = ") -
        nchar(if (x$theta > 1e+12) formatC(x$theta, digits = digits,
          format = "e") else formatC(x$theta, digits = digits,
          format = "f"))), collapse = ""), if (x$theta >
      1e+12)
      formatC(x$theta, digits = digits, format = "e") else formatC(x$theta, digits = digits, format = "f"),
      "\n")
  }
  if (x$nuZUvar > 1 || x$nvZVvar > 1) {
    cat("Variances averaged over observations \n")
  }
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat("Average inefficiency E[ui]     = ", paste0(rep(" ",
    lengthSum - 2 - nchar("Average inefficiency E[ui]     = ") -
      nchar(if (x$Eu > 1e+12) formatC(x$Eu, digits = digits,
        format = "e") else formatC(x$Eu, digits = digits,
        format = "f"))), collapse = ""), if (x$Eu > 1e+12)
    formatC(x$Eu, digits = digits, format = "e") else formatC(x$Eu, digits = digits, format = "f"), "\n")
  if (x$logDepVar == TRUE) {
    cat("Average efficiency E[exp(-ui)] = ", paste0(rep(" ",
      lengthSum - 2 - nchar("Average efficiency E[exp(-ui)] = ") -
        nchar(if (is.na(x$Expu)) {
          x$Expu
        } else {
          if (x$Expu < 1e-12) formatC(x$Expu, digits = digits,
          format = "e") else formatC(x$Expu, digits = digits,
          format = "f")
        })), collapse = ""), if (is.na(x$Expu)) {
      x$Expu
    } else {
      if (x$Expu < 1e-12)
        formatC(x$Expu, digits = digits, format = "e") else formatC(x$Expu, digits = digits, format = "f")
    }, "\n")
  }
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("-----[ Tests vs. No Inefficiency ]-----\n")
  cat("Likelihood Ratio Test of Inefficiency\n")
  cat("Deg. freedom for inefficiency model", paste0(rep(" ",
    lengthSum - 2 - nchar("Deg. freedom for inefficiency model") -
      nchar(formatC(x$df, digits = digits, format = "d"))),
    collapse = ""), formatC(x$df, digits = digits, format = "d"),
    "\n")
  cat("Log Likelihood for OLS Log(H0) = ", paste0(rep(" ",
    lengthSum - 2 - nchar("Log Likelihood for OLS Log(H0) = ") -
      nchar(formatC(x$olsLoglik, digits = digits, format = "f"))),
    collapse = ""), formatC(x$olsLoglik, digits = digits,
    format = "f"), "\n")
  cat("LR statistic: ", "\n")
  cat("Chisq = 2*[LogL(H0)-LogL(H1)]  = ", paste0(rep(" ",
    lengthSum - 2 - nchar("Chisq = 2*[LogL(H0)-LogL(H1)]  = ") -
      nchar(formatC(x$chisq, digits = digits, format = "f"))),
    collapse = ""), formatC(x$chisq, digits = digits, format = "f"),
    "\n")
  cat("Kodde-Palm C*:       95%:", formatC(qchibarsq(0.95,
    df = x$df), digits = digits, format = "f"), paste0(rep(" ",
    lengthSum - 2 - nchar("Kodde-Palm C*:       95%:") -
      nchar(formatC(qchibarsq(0.95, df = x$df), digits = digits,
        format = "f")) - nchar(formatC(qchibarsq(0.99,
      df = x$df), digits = digits, format = "f")) - nchar("99%") -
      3), collapse = ""), "99%:", formatC(qchibarsq(0.99,
    df = x$df), digits = digits, format = "f"), "\n")
  cat("Coelli (1995) skewness test on OLS residuals\n")
  cat("M3T: z                         = ", paste0(rep(" ",
    lengthSum - 2 - nchar("M3T: z                         = ") -
      nchar(formatC(x$CoelliM3Test[1], digits = digits,
        format = "f"))), collapse = ""), formatC(x$CoelliM3Test[1],
    digits = digits, format = "f"), "\n")
  cat("M3T: p.value                   = ", paste0(rep(" ",
    lengthSum - 2 - nchar("M3T: p.value                   = ") -
      nchar(formatC(x$CoelliM3Test[2], digits = digits,
        format = "f"))), collapse = ""), formatC(x$CoelliM3Test[2],
    digits = digits, format = "f"), "\n")
  cat("Final maximum likelihood estimates \n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Deterministic Component of SFA", width = lengthSum),
    "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes1, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  if (x$udist == "tnormal") {
    if (x$scaling) {
      cat(centerText("Scaling property parameters", width = lengthSum),
        "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes2, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes3, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes4, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes5, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
      if (all(mlRes5[, 4] > 0.1) || is.na(all(mlRes5[,
        4] > 0.1)) || all(is.na(mlRes5[, 4]))) {
        cat("---\n")
        cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
      }
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
    } else {
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes2, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes3, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes4, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
      if (all(mlRes4[, 4] > 0.1) || is.na(all(mlRes4[,
        4] > 0.1)) || all(is.na(mlRes4[, 4]))) {
        cat("---\n")
        cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
      }
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
    }
  } else {
    if (x$udist == "lognormal") {
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes2, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes3, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes4, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
      if (all(mlRes4[, 4] > 0.1) || is.na(all(mlRes4[,
        4] > 0.1)) || all(is.na(mlRes4[, 4]))) {
        cat("---\n")
        cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
      }
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
    } else {
      if (x$udist == "gamma") {
        cat(centerText("Parameter in variance of u (one-sided error)",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error)",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        cat(centerText("Location parameter P in u (one-sided error)",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        if (all(mlRes4[, 4] > 0.1) || is.na(all(mlRes4[,
          4] > 0.1)) || all(is.na(mlRes4[, 4]))) {
          cat("---\n")
          cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
        }
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
      } else {
        if (x$udist == "weibull") {
          cat(centerText("Parameter in variance of u (one-sided error)",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Parameters in variance of v (two-sided error)",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Shape parameter k in u (one-sided error)",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          if (all(mlRes4[, 4] > 0.1) || is.na(all(mlRes4[,
          4] > 0.1)) || all(is.na(mlRes4[, 4]))) {
          cat("---\n")
          cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
          }
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        } else {
          if (x$udist == "tslaplace") {
          cat(centerText("Parameter in variance of u (one-sided error)",
            width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          cat(centerText("Parameters in variance of v (two-sided error)",
            width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          cat(centerText("Skewness parameter 'lambda' in u (one-sided error)",
            width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          if (all(mlRes4[, 4] > 0.1) || is.na(all(mlRes4[,
            4] > 0.1)) || all(is.na(mlRes4[, 4]))) {
            cat("---\n")
            cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
          }
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          } else {
          cat(centerText("Parameter in variance of u (one-sided error)",
            width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          cat(centerText("Parameters in variance of v (two-sided error)",
            width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          if (all(mlRes3[, 4] > 0.1) || is.na(all(mlRes3[,
            4] > 0.1)) || all(is.na(mlRes3[, 4]))) {
            cat("---\n")
            cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
          }
          cat(paste0(rep("-", lengthSum), collapse = ""),
            "\n")
          }
        }
      }
    }
  }
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  invisible(x)
}

# summary for sfalcmcross ----------
#' @rdname summary
#' @aliases summary.sfalcmcross
#' @export
summary.sfalcmcross <- function(object, grad = FALSE, ci = FALSE,
  ...) {
  if (length(grad) != 1 || !is.logical(grad[1])) {
    stop("argument 'grad' must be a single logical value",
      call. = FALSE)
  }
  if (length(ci) != 1 || !is.logical(ci[1])) {
    stop("argument 'ci' must be a single logical value",
      call. = FALSE)
  }
  object$AIC <- -2 * object$mlLoglik + 2 * object$nParm
  object$BIC <- -2 * object$mlLoglik + log(object$Nobs) * object$nParm
  object$HQIC <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) *
    object$nParm
  # MLE estimates and stder, p-values, CI, Gradient
  if (grad && ci) {
    mlRes <- matrix(nrow = object$nParm, ncol = 7)
    colnames(mlRes) <- c("Coefficient", "Std. Error", "binf",
      "bsup", "gradient", "z-value", "Pr(>|z|)")
    mlRes[, 1] <- object$mlParam
    mlRes[, 2] <- sqrt(diag(object$invHessian))
    mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[, 2]
    mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[, 2]
    mlRes[, 5] <- object$gradient
    mlRes[, 6] <- mlRes[, 1]/mlRes[, 2]
    mlRes[, 7] <- 2 * pnorm(-abs(mlRes[, 6]))
  } else {
    if (grad == TRUE && ci == FALSE) {
      mlRes <- matrix(nrow = object$nParm, ncol = 5)
      colnames(mlRes) <- c("Coefficient", "Std. Error",
        "gradient", "z-value", "Pr(>|z|)")
      mlRes[, 1] <- object$mlParam
      mlRes[, 2] <- sqrt(diag(object$invHessian))
      mlRes[, 3] <- object$gradient
      mlRes[, 4] <- mlRes[, 1]/mlRes[, 2]
      mlRes[, 5] <- 2 * pnorm(-abs(mlRes[, 4]))
    } else {
      if (grad == FALSE && ci == TRUE) {
        mlRes <- matrix(nrow = object$nParm, ncol = 6)
        colnames(mlRes) <- c("Coefficient", "Std. Error",
          "binf", "bsup", "z-value", "Pr(>|z|)")
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[,
          2]
        mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[,
          2]
        mlRes[, 5] <- mlRes[, 1]/mlRes[, 2]
        mlRes[, 6] <- 2 * pnorm(-abs(mlRes[, 5]))
      } else {
        mlRes <- matrix(nrow = object$nParm, ncol = 4)
        colnames(mlRes) <- c("Coefficient", "Std. Error",
          "z value", "Pr(>|z|)")
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1]/mlRes[, 2]
        mlRes[, 4] <- 2 * pnorm(-abs(mlRes[, 3]))
      }
    }
  }
  row.names(mlRes) <- names(object$startVal)
  object$mlRes <- mlRes
  # object$chisq <- 2 * (object$mlLoglik -
  # object$olsLoglik)
  object$df <- object$nParm - object$nClasses * object$nXvar -
    object$nClasses * object$nvZVvar - object$nZHvar * (object$nClasses -
    1)
  class(object) <- "summary.sfalcmcross"
  return(object)
}

# print summary for sfalcmcross ----------
#' @rdname summary
#' @aliases print.summary.sfalcmcross
#' @export
print.summary.sfalcmcross <- function(x, digits = max(3, getOption("digits") -
  2), ...) {
  mlRes <- x$mlRes
  if (dim(mlRes)[2] == 4) {
    mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1], digits = digits,
      format = "f"))
    mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2], digits = digits,
      format = "f"))
    mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3], digits = digits,
      format = "f"))
    mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4], digits = digits,
      format = "e"))
  } else {
    if (dim(mlRes)[2] == 5) {
      mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1], digits = digits,
        format = "f"))
      mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2], digits = digits,
        format = "f"))
      mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3], digits = digits,
        format = "e"))
      mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4], digits = digits,
        format = "f"))
      mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5], digits = digits,
        format = "e"))
    } else {
      if (dim(mlRes)[2] == 6) {
        mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
          digits = digits, format = "f"))
        mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
          digits = digits, format = "f"))
        mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
          digits = digits, format = "f"))
        mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
          digits = digits, format = "f"))
        mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5],
          digits = digits, format = "f"))
        mlRes[, 6] <- as.numeric(formatC(x$mlRes[, 6],
          digits = digits, format = "e"))
      } else {
        if (dim(mlRes)[2] == 7) {
          mlRes[, 1] <- as.numeric(formatC(x$mlRes[,
          1], digits = digits, format = "f"))
          mlRes[, 2] <- as.numeric(formatC(x$mlRes[,
          2], digits = digits, format = "f"))
          mlRes[, 3] <- as.numeric(formatC(x$mlRes[,
          3], digits = digits, format = "f"))
          mlRes[, 4] <- as.numeric(formatC(x$mlRes[,
          4], digits = digits, format = "f"))
          mlRes[, 5] <- as.numeric(formatC(x$mlRes[,
          5], digits = digits, format = "e"))
          mlRes[, 6] <- as.numeric(formatC(x$mlRes[,
          6], digits = digits, format = "f"))
          mlRes[, 7] <- as.numeric(formatC(x$mlRes[,
          7], digits = digits, format = "e"))
        }
      }
    }
  }
  row.names(mlRes) <- formatC(row.names(mlRes), width = max(nchar(row.names(mlRes))),
    flag = "-")
  sfaModel <- "Normal-Half Normal Latent Class Stochastic Frontier Model"
  lengthSum <- 80
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(sfaModel, "\n")
  cat("Dependent Variable:", paste0(rep(" ", lengthSum - 2 -
    nchar("Dependent Variable:") - nchar(paste0(attr(x$formula,
    "lhs")))), collapse = ""), paste0(attr(x$formula, "lhs")),
    "\n")
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    2 - nchar("Log likelihood solver:") - nchar(x$optType)),
    collapse = ""), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - 2 -
    nchar("Log likelihood iter:") - nchar(x$nIter)), collapse = ""),
    x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    2 - nchar("Log likelihood value:") - nchar(if (x$mlLoglik >
    1e+12) formatC(x$mlLoglik, digits = digits, format = "e") else formatC(x$mlLoglik,
    digits = digits, format = "f"))), collapse = ""), if (x$mlLoglik >
    1e+12)
    formatC(x$mlLoglik, digits = digits, format = "e") else formatC(x$mlLoglik, digits = digits, format = "f"),
    "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    2 - nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
    digits = digits, format = "e"))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"), "\n")
  cat("Estimation based on:", paste0(rep(" ", lengthSum - 2 -
    nchar("Estimation based on:") - nchar(x$Nobs) - nchar(x$nParm) -
    nchar("N = ") - nchar("and K = ") - 3), collapse = ""),
    "N = ", x$Nobs, "and K = ", x$nParm, "\n")
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - 2 - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(if (abs(x$AIC) > 1e+12) formatC(x$AIC,
    digits = 1, format = "e") else formatC(x$AIC, digits = 1,
    format = "f")) - nchar("AIC/N  = ") - nchar(if (abs(x$AIC/x$Nobs) >
    1e+12) formatC(x$AIC/x$Nobs, digits = 3, format = "e") else formatC(x$AIC/x$Nobs,
    digits = 3, format = "f")) - 3), collapse = ""), "AIC  = ",
    if (abs(x$AIC) > 1e+12)
      formatC(x$AIC, digits = 1, format = "e") else formatC(x$AIC, digits = 1, format = "f"), "AIC/N  = ",
    if (abs(x$AIC/x$Nobs) > 1e+12)
      formatC(x$AIC/x$Nobs, digits = 3, format = "e") else formatC(x$AIC/x$Nobs, digits = 3, format = "f"),
    "\n")
  cat(paste0(rep(" ", lengthSum - 2 - nchar("BIC  = ") - nchar(if (abs(x$BIC) >
    1e+12) formatC(x$BIC, digits = 1, format = "e") else formatC(x$BIC,
    digits = 1, format = "f")) - nchar("BIC/N  = ") - nchar(if (abs(x$BIC/x$Nobs) >
    1e+12) formatC(x$BIC/x$Nobs, digits = 3, format = "e") else formatC(x$BIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "BIC  = ",
    if (abs(x$BIC) > 1e+12)
      formatC(x$BIC, digits = 1, format = "e") else formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    if (abs(x$BIC/x$Nobs) > 1e+12)
      formatC(x$BIC/x$Nobs, digits = 3, format = "e") else formatC(x$BIC/x$Nobs, digits = 3, format = "f"),
    "\n")
  cat(paste0(rep(" ", lengthSum - 2 - nchar("HQIC = ") - nchar(if (abs(x$HQIC) >
    1e+12) formatC(x$HQIC, digits = 1, format = "e") else formatC(x$HQIC,
    digits = 1, format = "f")) - nchar("HQIC/N = ") - nchar(if (abs(x$HQIC/x$Nobs) >
    1e+12) formatC(x$HQIC/x$Nobs, digits = 3, format = "e") else formatC(x$HQIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "HQIC = ",
    if (abs(x$HQIC) > 1e+12)
      formatC(x$HQIC, digits = 1, format = "e") else formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    if (abs(x$HQIC/x$Nobs) > 1e+12)
      formatC(x$HQIC/x$Nobs, digits = 3, format = "e") else formatC(x$HQIC/x$Nobs, digits = 3, format = "f"),
    "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("Latent class model with", x$nClasses, "latent classes \n")
  cat("Final maximum likelihood estimates \n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Deterministic Component of SFA for latent class 1",
    width = lengthSum), "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes[1:x$nXvar, , drop = FALSE], P.values = TRUE,
    digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Parameter in variance of u (one-sided error) for latent class 1",
    width = lengthSum), "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar), ,
    drop = FALSE], P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Parameters in variance of v (two-sided error) for latent class 1",
    width = lengthSum), "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar + x$nuZUvar +
    x$nvZVvar), , drop = FALSE], P.values = TRUE, digits = digits,
    signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Deterministic Component of SFA for latent class 2",
    width = lengthSum), "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes[(x$nXvar + x$nuZUvar + x$nvZVvar + 1):(2 *
    x$nXvar + x$nuZUvar + x$nvZVvar), , drop = FALSE], P.values = TRUE,
    digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Parameter in variance of u (one-sided error) for latent class 2",
    width = lengthSum), "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes[(2 * x$nXvar + x$nuZUvar + x$nvZVvar +
    1):(2 * x$nXvar + 2 * x$nuZUvar + x$nvZVvar), , drop = FALSE],
    P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Parameters in variance of v (two-sided error) for latent class 2",
    width = lengthSum), "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar + x$nvZVvar +
    1):(2 * x$nXvar + 2 * x$nuZUvar + 2 * x$nvZVvar), , drop = FALSE],
    P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  if (x$nClasses == 2) {
    cat(centerText("Estimated prior probabilities for class membership",
      width = lengthSum), "\n")
    cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
    printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar + 2 *
      x$nvZVvar + 1):(2 * x$nXvar + 2 * x$nuZUvar + 2 *
      x$nvZVvar + x$nZHvar), , drop = FALSE], P.values = TRUE,
      digits = digits, signif.legend = TRUE)
    if (all(mlRes[(2 * x$nXvar + 2 * x$nuZUvar + 2 * x$nvZVvar +
      1):(2 * x$nXvar + 2 * x$nuZUvar + 2 * x$nvZVvar +
      x$nZHvar), 4] > 0.1) || is.na(all(mlRes[(2 * x$nXvar +
      2 * x$nuZUvar + 2 * x$nvZVvar + 1):(2 * x$nXvar +
      2 * x$nuZUvar + 2 * x$nvZVvar + x$nZHvar), 4] > 0.1)) ||
      all(is.na(mlRes[(2 * x$nXvar + 2 * x$nuZUvar + 2 *
        x$nvZVvar + 1):(2 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar + x$nZHvar), 4]))) {
      cat("---\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    }
    cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  } else {
    if (x$nClasses == 3) {
      cat(centerText("Deterministic Component of SFA for latent class 3",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar + 1):(3 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
        digits = digits, signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Parameter in variance of u (one-sided error) for latent class 3",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes[(3 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
        digits = digits, signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Parameters in variance of v (two-sided error) for latent class 3",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
        2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
        digits = digits, signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      cat(centerText("Estimated prior probabilities for class membership",
        width = lengthSum), "\n")
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
      printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar + 2 * x$nZHvar), , drop = FALSE],
        P.values = TRUE, digits = digits, signif.legend = TRUE)
      if (all(mlRes[(3 * x$nXvar + 3 * x$nuZUvar + 3 *
        x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar + 2 * x$nZHvar), 4] > 0.1) || is.na(all(mlRes[(3 *
        x$nXvar + 3 * x$nuZUvar + 3 * x$nvZVvar + 1):(3 *
        x$nXvar + 3 * x$nuZUvar + 3 * x$nvZVvar + 2 *
        x$nZHvar), 4] > 0.1)) || all(is.na(mlRes[(3 *
        x$nXvar + 3 * x$nuZUvar + 3 * x$nvZVvar + 1):(3 *
        x$nXvar + 3 * x$nuZUvar + 3 * x$nvZVvar + 2 *
        x$nZHvar), 4]))) {
        cat("---\n")
        cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
      }
      cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
    } else {
      if (x$nClasses == 4) {
        cat(centerText("Deterministic Component of SFA for latent class 3",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        cat(centerText("Parameter in variance of u (one-sided error) for latent class 3",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes[(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error) for latent class 3",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        cat(centerText("Deterministic Component of SFA for latent class 4",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        cat(centerText("Parameter in variance of u (one-sided error) for latent class 4",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes[(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error) for latent class 4",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        cat(centerText("Estimated prior probabilities for class membership",
          width = lengthSum), "\n")
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 3 * x$nZHvar), , drop = FALSE],
          P.values = TRUE, digits = digits, signif.legend = TRUE)
        if (all(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 3 * x$nZHvar), 4] > 0.1) ||
          is.na(all(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 3 * x$nZHvar), 4] > 0.1)) ||
          all(is.na(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 3 * x$nZHvar), 4]))) {
          cat("---\n")
          cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
        }
        cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
      } else {
        if (x$nClasses == 5) {
          cat(centerText("Deterministic Component of SFA for latent class 3",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Parameter in variance of u (one-sided error) for latent class 3",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Parameters in variance of v (two-sided error) for latent class 3",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Deterministic Component of SFA for latent class 4",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Parameter in variance of u (one-sided error) for latent class 4",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Parameters in variance of v (two-sided error) for latent class 4",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Deterministic Component of SFA for latent class 5",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(5 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Parameter in variance of u (one-sided error) for latent class 5",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(5 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
          4 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Parameters in variance of v (two-sided error) for latent class 5",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(5 * x$nXvar + 5 * x$nuZUvar +
          4 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
          5 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          cat(centerText("Estimated prior probabilities for class membership",
          width = lengthSum), "\n")
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
          printCoefmat(mlRes[(5 * x$nXvar + 5 * x$nuZUvar +
          5 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
          5 * x$nvZVvar + 4 * x$nZHvar), , drop = FALSE],
          P.values = TRUE, digits = digits, signif.legend = TRUE)
          if (all(mlRes[(5 * x$nXvar + 5 * x$nuZUvar +
          5 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
          5 * x$nvZVvar + 4 * x$nZHvar), 4] > 0.1) ||
          is.na(all(mlRes[(5 * x$nXvar + 5 * x$nuZUvar +
            5 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
            5 * x$nvZVvar + 4 * x$nZHvar), 4] > 0.1)) ||
          all(is.na(mlRes[(5 * x$nXvar + 5 * x$nuZUvar +
            5 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
            5 * x$nvZVvar + 4 * x$nZHvar), 4]))) {
          cat("---\n")
          cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
          }
          cat(paste0(rep("-", lengthSum), collapse = ""),
          "\n")
        }
      }
    }
  }
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  invisible(x)
}

# summary for sfaselectioncross ----------
#' @rdname summary
#' @aliases summary.sfaselectioncross
#' @export
summary.sfaselectioncross <- function(object, grad = FALSE, ci = FALSE,
  ...) {
  if (length(grad) != 1 || !is.logical(grad[1])) {
    stop("argument 'grad' must be a single logical value",
      call. = FALSE)
  }
  if (length(ci) != 1 || !is.logical(ci[1])) {
    stop("argument 'ci' must be a single logical value",
      call. = FALSE)
  }
  object$AIC <- -2 * object$mlLoglik + 2 * object$nParm
  object$BIC <- -2 * object$mlLoglik + log(object$Nobs) * object$nParm
  object$HQIC <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) *
    object$nParm
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable[all.vars(object$selectionF)[1]] ==
    1, ], rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable[all.vars(object$selectionF)[1]] ==
    1, ], rhs = 3)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  object$sigmavSq <- mean(exp(Wv))
  object$sigmauSq <- mean(exp(Wu))
  object$Varu <- varuFun(object = object, mu = NULL, P = NULL,
    k = NULL, lambda = NULL)
  object$Eu <- euFun(object = object, mu = NULL, P = NULL,
    k = NULL, lambda = NULL)
  if (object$logDepVar == TRUE) {
    object$Expu <- eExpuFun(object = object, mu = NULL, P = NULL,
      k = NULL, lambda = NULL)
  }
  # OLS 2 Step estimates and stder, p-values, CI
  dfOLS <- object$Nobs - object$nXvar
  if (ci) {
    ols2StepRes <- matrix(nrow = object$nXvar + 2, ncol = 6)
    colnames(ols2StepRes) <- c("Coefficient", "Std. Error",
      "binf", "bsup", "t value", "Pr(>|t|)")
    ols2StepRes[, 1] <- c(object$ols2stepParam, object$ols2stepSigmasq)
    ols2StepRes[, 2] <- c(object$ols2stepStder, NA)
    ols2StepRes[, 3] <- ols2StepRes[, 1] - qt(0.975, df = dfOLS) *
      ols2StepRes[, 2]
    ols2StepRes[, 4] <- ols2StepRes[, 1] + qt(0.975, df = dfOLS) *
      ols2StepRes[, 2]
    ols2StepRes[, 5] <- ols2StepRes[, 1]/ols2StepRes[, 2]
    ols2StepRes[, 6] <- 2 * pt(-abs(ols2StepRes[, 5]), df = dfOLS)
  } else {
    ols2StepRes <- matrix(nrow = object$nXvar + 2, ncol = 4)
    colnames(ols2StepRes) <- c("Coefficient", "Std. Error",
      "t value", "Pr(>|t|)")
    ols2StepRes[, 1] <- c(object$ols2stepParam, object$ols2stepSigmasq)
    ols2StepRes[, 2] <- c(object$ols2stepStder, NA)
    ols2StepRes[, 3] <- ols2StepRes[, 1]/ols2StepRes[, 2]
    ols2StepRes[, 4] <- 2 * pt(-abs(ols2StepRes[, 3]), df = dfOLS)
  }
  row.names(ols2StepRes) <- c(names(object$ols2stepParam),
    "sigmaSq")
  object$ols2StepRes <- ols2StepRes
  # MLE estimates and stder, p-values, CI, Gradient
  if (grad && ci) {
    mlRes <- matrix(nrow = object$nParm, ncol = 7)
    colnames(mlRes) <- c("Coefficient", "Std. Error", "binf",
      "bsup", "gradient", "z-value", "Pr(>|z|)")
    mlRes[, 1] <- object$mlParam
    mlRes[, 2] <- sqrt(diag(object$invHessian))
    mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[, 2]
    mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[, 2]
    mlRes[, 5] <- object$gradient
    mlRes[, 6] <- mlRes[, 1]/mlRes[, 2]
    mlRes[, 7] <- 2 * pnorm(-abs(mlRes[, 6]))
  } else {
    if (grad == TRUE && ci == FALSE) {
      mlRes <- matrix(nrow = object$nParm, ncol = 5)
      colnames(mlRes) <- c("Coefficient", "Std. Error",
        "gradient", "z-value", "Pr(>|z|)")
      mlRes[, 1] <- object$mlParam
      mlRes[, 2] <- sqrt(diag(object$invHessian))
      mlRes[, 3] <- object$gradient
      mlRes[, 4] <- mlRes[, 1]/mlRes[, 2]
      mlRes[, 5] <- 2 * pnorm(-abs(mlRes[, 4]))
    } else {
      if (grad == FALSE && ci == TRUE) {
        mlRes <- matrix(nrow = object$nParm, ncol = 6)
        colnames(mlRes) <- c("Coefficient", "Std. Error",
          "binf", "bsup", "z-value", "Pr(>|z|)")
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[,
          2]
        mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[,
          2]
        mlRes[, 5] <- mlRes[, 1]/mlRes[, 2]
        mlRes[, 6] <- 2 * pnorm(-abs(mlRes[, 5]))
      } else {
        mlRes <- matrix(nrow = object$nParm, ncol = 4)
        colnames(mlRes) <- c("Coefficient", "Std. Error",
          "z value", "Pr(>|z|)")
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1]/mlRes[, 2]
        mlRes[, 4] <- 2 * pnorm(-abs(mlRes[, 3]))
      }
    }
  }
  row.names(mlRes) <- names(object$startVal)
  object$mlRes <- mlRes
  object$df <- object$nParm - object$nXvar - object$nvZVvar -
    1
  class(object) <- "summary.sfaselectioncross"
  return(object)
}

# print summary for sfaselectioncross ----------
#' @rdname summary
#' @aliases print.summary.sfaselectioncross
#' @export
print.summary.sfaselectioncross <- function(x, digits = max(3,
  getOption("digits") - 2), ...) {
  mlRes <- x$mlRes
  if (dim(mlRes)[2] == 4) {
    mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1], digits = digits,
      format = "f"))
    mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2], digits = digits,
      format = "f"))
    mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3], digits = digits,
      format = "f"))
    mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4], digits = digits,
      format = "e"))
  } else {
    if (dim(mlRes)[2] == 5) {
      mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1], digits = digits,
        format = "f"))
      mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2], digits = digits,
        format = "f"))
      mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3], digits = digits,
        format = "e"))
      mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4], digits = digits,
        format = "f"))
      mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5], digits = digits,
        format = "e"))
    } else {
      if (dim(mlRes)[2] == 6) {
        mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
          digits = digits, format = "f"))
        mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
          digits = digits, format = "f"))
        mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
          digits = digits, format = "f"))
        mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
          digits = digits, format = "f"))
        mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5],
          digits = digits, format = "f"))
        mlRes[, 6] <- as.numeric(formatC(x$mlRes[, 6],
          digits = digits, format = "e"))
      } else {
        if (dim(mlRes)[2] == 7) {
          mlRes[, 1] <- as.numeric(formatC(x$mlRes[,
          1], digits = digits, format = "f"))
          mlRes[, 2] <- as.numeric(formatC(x$mlRes[,
          2], digits = digits, format = "f"))
          mlRes[, 3] <- as.numeric(formatC(x$mlRes[,
          3], digits = digits, format = "f"))
          mlRes[, 4] <- as.numeric(formatC(x$mlRes[,
          4], digits = digits, format = "f"))
          mlRes[, 5] <- as.numeric(formatC(x$mlRes[,
          5], digits = digits, format = "e"))
          mlRes[, 6] <- as.numeric(formatC(x$mlRes[,
          6], digits = digits, format = "f"))
          mlRes[, 7] <- as.numeric(formatC(x$mlRes[,
          7], digits = digits, format = "e"))
        }
      }
    }
  }
  row.names(mlRes) <- formatC(row.names(mlRes), width = max(nchar(row.names(mlRes))),
    flag = "-")
  mlRes1 <- mlRes[1:x$nXvar, ]
  mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar), , drop = FALSE]
  mlRes3 <- mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar + x$nuZUvar +
    x$nvZVvar), , drop = FALSE]
  mlRes4 <- mlRes[x$nXvar + x$nuZUvar + x$nvZVvar + 1, , drop = FALSE]
  sfaModel <- "Sample Selection Correction Stochastic Frontier Model"
  lengthSum <- 80
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(sfaModel, "\n")
  cat("Dependent Variable:", paste0(rep(" ", lengthSum - 2 -
    nchar("Dependent Variable:") - nchar(paste0(attr(x$formula,
    "lhs")))), collapse = ""), paste0(attr(x$formula, "lhs")),
    "\n")
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    2 - nchar("Log likelihood solver:") - nchar(x$optType)),
    collapse = ""), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - 2 -
    nchar("Log likelihood iter:") - nchar(x$nIter)), collapse = ""),
    x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    2 - nchar("Log likelihood value:") - nchar(if (x$mlLoglik >
    1e+12) formatC(x$mlLoglik, digits = digits, format = "e") else formatC(x$mlLoglik,
    digits = digits, format = "f"))), collapse = ""), if (x$mlLoglik >
    1e+12)
    formatC(x$mlLoglik, digits = digits, format = "e") else formatC(x$mlLoglik, digits = digits, format = "f"),
    "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    2 - nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
    digits = digits, format = "e"))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"), "\n")
  cat("Estimation based on:", paste0(rep(" ", lengthSum - 2 -
    nchar("Estimation based on:") - nchar(x$Nobs) - nchar(x$nParm) -
    nchar("N = ") - nchar("of") - nchar(x$Ninit) - nchar("obs.") -
    nchar("and K = ") - 6), collapse = ""), "N = ", x$Nobs,
    "of", x$Ninit, "obs.", "and K = ", x$nParm, "\n")
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - 2 - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(if (abs(x$AIC) > 1e+12) formatC(x$AIC,
    digits = 1, format = "e") else formatC(x$AIC, digits = 1,
    format = "f")) - nchar("AIC/N  = ") - nchar(if (abs(x$AIC/x$Nobs) >
    1e+12) formatC(x$AIC/x$Nobs, digits = 3, format = "e") else formatC(x$AIC/x$Nobs,
    digits = 3, format = "f")) - 3), collapse = ""), "AIC  = ",
    if (abs(x$AIC) > 1e+12)
      formatC(x$AIC, digits = 1, format = "e") else formatC(x$AIC, digits = 1, format = "f"), "AIC/N  = ",
    if (abs(x$AIC/x$Nobs) > 1e+12)
      formatC(x$AIC/x$Nobs, digits = 3, format = "e") else formatC(x$AIC/x$Nobs, digits = 3, format = "f"),
    "\n")
  cat(paste0(rep(" ", lengthSum - 2 - nchar("BIC  = ") - nchar(if (abs(x$BIC) >
    1e+12) formatC(x$BIC, digits = 1, format = "e") else formatC(x$BIC,
    digits = 1, format = "f")) - nchar("BIC/N  = ") - nchar(if (abs(x$BIC/x$Nobs) >
    1e+12) formatC(x$BIC/x$Nobs, digits = 3, format = "e") else formatC(x$BIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "BIC  = ",
    if (abs(x$BIC) > 1e+12)
      formatC(x$BIC, digits = 1, format = "e") else formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    if (abs(x$BIC/x$Nobs) > 1e+12)
      formatC(x$BIC/x$Nobs, digits = 3, format = "e") else formatC(x$BIC/x$Nobs, digits = 3, format = "f"),
    "\n")
  cat(paste0(rep(" ", lengthSum - 2 - nchar("HQIC = ") - nchar(if (abs(x$HQIC) >
    1e+12) formatC(x$HQIC, digits = 1, format = "e") else formatC(x$HQIC,
    digits = 1, format = "f")) - nchar("HQIC/N = ") - nchar(if (abs(x$HQIC/x$Nobs) >
    1e+12) formatC(x$HQIC/x$Nobs, digits = 3, format = "e") else formatC(x$HQIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "HQIC = ",
    if (abs(x$HQIC) > 1e+12)
      formatC(x$HQIC, digits = 1, format = "e") else formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    if (abs(x$HQIC/x$Nobs) > 1e+12)
      formatC(x$HQIC/x$Nobs, digits = 3, format = "e") else formatC(x$HQIC/x$Nobs, digits = 3, format = "f"),
    "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat("Variances: Sigma-squared(v)   = ", paste0(rep(" ", lengthSum -
    2 - nchar("Variances: Sigma-squared(v)   = ") - nchar(if (x$sigmavSq >
    1e+12) formatC(x$sigmavSq, digits = digits, format = "e") else formatC(x$sigmavSq,
    digits = digits, format = "f"))), collapse = ""), if (x$sigmavSq >
    1e+12)
    formatC(x$sigmavSq, digits = digits, format = "e") else formatC(x$sigmavSq, digits = digits, format = "f"),
    "\n")
  cat("           Sigma(v)           = ", paste0(rep(" ", lengthSum -
    2 - nchar("Variances: Sigma-squared(v)   = ") - nchar(if (x$sigmavSq >
    1e+12) formatC(x$sigmavSq, digits = digits, format = "e") else formatC(x$sigmavSq,
    digits = digits, format = "f"))), collapse = ""), if (x$sigmavSq >
    1e+12)
    formatC(x$sigmavSq, digits = digits, format = "e") else formatC(x$sigmavSq, digits = digits, format = "f"),
    "\n")
  cat("           Sigma-squared(u)   = ", paste0(rep(" ", lengthSum -
    2 - nchar("Variances: Sigma-squared(u)   = ") - nchar(if (x$sigmauSq >
    1e+12) formatC(x$sigmauSq, digits = digits, format = "e") else formatC(x$sigmauSq,
    digits = digits, format = "f"))), collapse = ""), if (x$sigmauSq >
    1e+12)
    formatC(x$sigmauSq, digits = digits, format = "e") else formatC(x$sigmauSq, digits = digits, format = "f"),
    "\n")
  cat("           Sigma(u)           = ", paste0(rep(" ", lengthSum -
    2 - nchar("Variances: Sigma-squared(u)   = ") - nchar(if (x$sigmauSq >
    1e+12) formatC(x$sigmauSq, digits = digits, format = "e") else formatC(x$sigmauSq,
    digits = digits, format = "f"))), collapse = ""), if (x$sigmauSq >
    1e+12)
    formatC(x$sigmauSq, digits = digits, format = "e") else formatC(x$sigmauSq, digits = digits, format = "f"),
    "\n")
  cat("Sigma = Sqrt[(s^2(u)+s^2(v))] = ", paste0(rep(" ", lengthSum -
    2 - nchar("Sigma = Sqrt[(s^2(u)+s^2(v))] = ") - nchar(if (sqrt(x$sigmavSq +
    x$sigmauSq) > 1e+12) formatC(sqrt(x$sigmavSq + x$sigmauSq),
    digits = digits, format = "e") else formatC(sqrt(x$sigmavSq +
    x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    if (sqrt(x$sigmavSq + x$sigmauSq) > 1e+12)
      formatC(sqrt(x$sigmavSq + x$sigmauSq), digits = digits,
        format = "e") else formatC(sqrt(x$sigmavSq + x$sigmauSq), digits = digits,
      format = "f"), "\n")
  cat("Gamma = sigma(u)^2/sigma^2    = ", paste0(rep(" ", lengthSum -
    2 - nchar("Gamma = sigma(u)^2/sigma^2    = ") - nchar(if (x$sigmauSq/(x$sigmavSq +
    x$sigmauSq) > 1e+12) formatC(x$sigmauSq/(x$sigmavSq +
    x$sigmauSq), digits = digits, format = "e") else formatC(x$sigmauSq/(x$sigmavSq +
    x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    if (x$sigmauSq/(x$sigmavSq + x$sigmauSq) > 1e+12)
      formatC(x$sigmauSq/(x$sigmavSq + x$sigmauSq), digits = digits,
        format = "e") else formatC(x$sigmauSq/(x$sigmavSq + x$sigmauSq), digits = digits,
      format = "f"), "\n")
  cat("Lambda = sigma(u)/sigma(v)    = ", paste0(rep(" ", lengthSum -
    2 - nchar("Lambda = sigma(u)/sigma(v)    = ") - nchar(if (sqrt(x$sigmauSq/x$sigmavSq) >
    1e+12) formatC(sqrt(x$sigmauSq/x$sigmavSq), digits = digits,
    format = "e") else formatC(sqrt(x$sigmauSq/x$sigmavSq),
    digits = digits, format = "f"))), collapse = ""), if (sqrt(x$sigmauSq/x$sigmavSq) >
    1e+12)
    formatC(sqrt(x$sigmauSq/x$sigmavSq), digits = digits,
      format = "e") else formatC(sqrt(x$sigmauSq/x$sigmavSq), digits = digits,
    format = "f"), "\n")
  cat("Var[u]/{Var[u]+Var[v]}        = ", paste0(rep(" ", lengthSum -
    2 - nchar("Var[u]/{Var[u]+Var[v]}        = ") - nchar(if (x$Varu/(x$Varu +
    x$sigmavSq) > 1e+12) formatC(x$Varu/(x$Varu + x$sigmavSq),
    digits = digits, format = "e") else formatC(x$Varu/(x$Varu +
    x$sigmavSq), digits = digits, format = "f"))), collapse = ""),
    if (x$Varu/(x$Varu + x$sigmavSq) > 1e+12)
      formatC(x$Varu/(x$Varu + x$sigmavSq), digits = digits,
        format = "e") else formatC(x$Varu/(x$Varu + x$sigmavSq), digits = digits,
      format = "f"), "\n")
  if (x$nuZUvar > 1 || x$nvZVvar > 1) {
    cat("Variances averaged over observations \n")
  }
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat("Average inefficiency E[ui]     = ", paste0(rep(" ",
    lengthSum - 2 - nchar("Average inefficiency E[ui]     = ") -
      nchar(if (x$Eu > 1e+12) formatC(x$Eu, digits = digits,
        format = "e") else formatC(x$Eu, digits = digits,
        format = "f"))), collapse = ""), if (x$Eu > 1e+12)
    formatC(x$Eu, digits = digits, format = "e") else formatC(x$Eu, digits = digits, format = "f"), "\n")
  if (x$logDepVar == TRUE) {
    cat("Average efficiency E[exp(-ui)] = ", paste0(rep(" ",
      lengthSum - 2 - nchar("Average efficiency E[exp(-ui)] = ") -
        nchar(if (is.na(x$Expu)) {
          x$Expu
        } else {
          if (x$Expu < 1e-12) formatC(x$Expu, digits = digits,
          format = "e") else formatC(x$Expu, digits = digits,
          format = "f")
        })), collapse = ""), if (is.na(x$Expu)) {
      x$Expu
    } else {
      if (x$Expu < 1e-12)
        formatC(x$Expu, digits = digits, format = "e") else formatC(x$Expu, digits = digits, format = "f")
    }, "\n")
  }
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("Estimator is 2 step Maximum Likelihood \n")
  cat("Final maximum likelihood estimates \n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Deterministic Component of SFA", width = lengthSum),
    "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes1, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Parameter in variance of u (one-sided error)",
    width = lengthSum), "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes2, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Parameters in variance of v (two-sided error)",
    width = lengthSum), "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes3, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(centerText("Selection bias parameter", width = lengthSum),
    "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  printCoefmat(mlRes4, P.values = TRUE, digits = digits, signif.legend = TRUE)
  if (all(mlRes4[, 4] > 0.1) || is.na(all(mlRes4[, 4] > 0.1)) ||
    all(is.na(mlRes4[, 4]))) {
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(paste0(rep("-", lengthSum), collapse = ""), "\n")
  invisible(x)
}
