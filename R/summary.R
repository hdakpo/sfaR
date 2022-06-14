################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Summary of optimization objects                                              #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
#         -Sample selection correction                                         #
#         -Zero inefficiency stochastic frontier                               #
#         -Contaminated noise stochastic frontier                              #
#         -Multi-Modal Inefficiency Stochastic Frontier Analysis               #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Summary of results for stochastic frontier models
#'
#' Create and print summary results for stochastic frontier models returned by 
#' \code{\link{cnsfcross}}, \code{\link{lcmcross}}, \code{\link{misfcross}}, 
#' \code{\link{sfacross}}, \code{\link{sfaselectioncross}} or 
#' \code{\link{zisfcross}}.
#'
#' @param object An object of either class \code{'cnsfcross'} returned by the
#' function \code{\link{cnsfcross}}, or class \code{'lcmcross'} returned by the
#' function \code{\link{lcmcross}}, or \code{'misfcross'} returned by 
#' \code{\link{cnsfcross}}, or \code{'sfacross'} returned by the
#' function \code{\link{sfacross}}, or class \code{'sfaselectioncross'} returned 
#' by the function \code{\link{sfaselectioncross}}, or class \code{'zisfcross'} 
#' returned by the function \code{\link{zisfcross}}.
#' @param grad Logical. Default = \code{FALSE}. If \code{TRUE}, the gradient
#' for the maximum likelihood (ML) estimates of the different parameters is
#' returned.
#' @param ci Logical. Default = \code{FALSE}. If \code{TRUE}, the 95\%
#' confidence interval for the different parameters (OLS and ML estimates) is
#' returned.
#' @param ... Currently ignored.
#' @param x An object of either class \code{'summary.cnsfcross'},
#' \code{'summary.lcmcross'}, \code{'summary.misfcross'}, 
#' \code{'summary.sfacross'}, \code{'summary.sfaselectioncross'} or 
#' \code{'summary.zisfcross'}.
#' @param digits Numeric. Number of digits displayed in values.
#'
#' @name summary
#'
#' @return The \code{\link{summary}} method returns a list of class
#' \code{'summary.cnsfcross'}, \code{'summary.lcmcross'}, 
#' \code{'summary.misfcross'}, \code{'summary.sfacross'}, 
#' \code{'summary.sfaselectioncross'}, \code{'summary.zisfcross'} that contains
#' the same elements as an object returned by \code{\link{cnsfcross}}, 
#' \code{\link{lcmcross}},  \code{\link{misfcross}}, \code{\link{sfacross}}, 
#' \code{\link{sfaselectioncross}} or \code{\link{zisfcross}} with the 
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
#' \item{THETA}{For \code{object} of class \code{'sfacross'} with \code{'udist
#' = uniform'}.  \eqn{\Theta} value in the case the uniform distribution is
#' defined as: \eqn{u_i \in [0, \Theta]}.}
#'
#' \item{Eu}{For \code{object} of class \code{'sfacross'} or 
#' \code{'sfaselectioncross'}. Expected unconditional inefficiency.}
#'
#' \item{Expu}{For \code{object} of class \code{'sfacross'} or 
#' \code{'sfaselectioncross'}. Expected unconditional efficiency.}
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
#' @seealso \code{\link{cnsfcross}}, for the contaminated noise stochastic 
#' frontier analysis model fitting function.
#'
#' \code{\link{lcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function.
#' 
#' \code{\link{misfcross}}, for the multi-modal inefficiency stochastic frontier 
#' analysis model fitting function.
#' 
#' \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function.
#' 
#' \code{\link{zisfcross}} for zero inefficiency in stochastic frontier model
#' fitting function.
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
      delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        (object$nuZUvar - 1))]
      tau <- object$mlParam[object$nXvar + (object$nuZUvar -
        1) + 1]
      cu <- object$mlParam[object$nXvar + (object$nuZUvar -
        1) + 2]
      phi <- object$mlParam[(object$nXvar + (object$nuZUvar -
        1) + 2 + 1):(object$nXvar + (object$nuZUvar -
        1) + 2 + object$nvZVvar)]
      muHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 2)
      uHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 3)
      vHvar <- model.matrix(object$formula, data = object$dataTable,
        rhs = 4)
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
        -1])))) * tau
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
      -1])))
  } else {
    as.numeric(crossprod(matrix(delta), t(uHvar)))
  }
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  object$sigmavSq <- mean(exp(Wv))
  object$sigmauSq <- mean(exp(Wu))
  object$Varu <- varuFun(object = object, mu = mu, P = P, k = k,
    lambda = lambda)
  if (object$udist == "uniform") {
    object$THETA <- sqrt(12 * object$sigmauSq)
  }
  object$Eu <- euFun(object = object, mu = mu, P = P, k = k,
    lambda = lambda)
  object$Expu <- eExpuFun(object = object, mu = mu, P = P,
    k = k, lambda = lambda)
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
  lengthSum <- nchar(sfadist(x$udist)) + 27
  dimCoefTable <- as.character(dim(x$mlRes)[2])
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(sfadist(x$udist), "\n")
  cat("Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
    nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n")
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
    collapse = ""), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mlLoglik,
    digits = digits, format = "f"))), collapse = ""), formatC(x$mlLoglik,
    digits = digits, format = "f"), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
    digits = digits, format = "e"))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"), "\n")
  cat("Estimation based on:", if (lengthSum - nchar("Estimation based on:") -
    nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("and K = ") -
    3 > 0)
    paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") -
      nchar("and K = ") - 3), collapse = ""), "N = ", x$Nobs,
    "and K = ", x$nParm, "\n")
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC/x$Nobs, digits = 3,
    format = "f")) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"), "AIC/N  = ", formatC(x$AIC/x$Nobs,
    digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
    digits = 1, format = "f")) - nchar("BIC/N  = ") - nchar(formatC(x$BIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
    digits = 1, format = "f")) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat("Variances: Sigma-squared(v)   = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(v)   = ") - nchar(formatC(x$sigmavSq,
    digits = digits, format = "f"))), collapse = ""), formatC(x$sigmavSq,
    digits = digits, format = "f"), "\n")
  cat("           Sigma(v)           = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(v)   = ") - nchar(formatC(sqrt(x$sigmavSq),
    digits = digits, format = "f"))), collapse = ""), formatC(sqrt(x$sigmavSq),
    digits = digits, format = "f"), "\n")
  cat("           Sigma-squared(u)   = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(u)   = ") - nchar(formatC(x$sigmauSq,
    digits = digits, format = "f"))), collapse = ""), formatC(x$sigmauSq,
    digits = digits, format = "f"), "\n")
  cat("           Sigma(u)           = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(u)   = ") - nchar(formatC(sqrt(x$sigmauSq),
    digits = digits, format = "f"))), collapse = ""), formatC(sqrt(x$sigmauSq),
    digits = digits, format = "f"), "\n")
  cat("Sigma = Sqrt[(s^2(u)+s^2(v))] = ", paste0(rep(" ", lengthSum -
    nchar("Sigma = Sqrt[(s^2(u)+s^2(v))] = ") - nchar(formatC(sqrt(x$sigmavSq +
    x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    formatC(sqrt(x$sigmavSq + x$sigmauSq), digits = digits,
      format = "f"), "\n")
  cat("Gamma = sigma(u)^2/sigma^2    = ", paste0(rep(" ", lengthSum -
    nchar("Gamma = sigma(u)^2/sigma^2    = ") - nchar(formatC(x$sigmauSq/(x$sigmavSq +
    x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    formatC(x$sigmauSq/(x$sigmavSq + x$sigmauSq), digits = digits,
      format = "f"), "\n")
  cat("Lambda = sigma(u)/sigma(v)    = ", paste0(rep(" ", lengthSum -
    nchar("Lambda = sigma(u)/sigma(v)    = ") - nchar(formatC(sqrt(x$sigmauSq/x$sigmavSq),
    digits = digits, format = "f"))), collapse = ""), formatC(sqrt(x$sigmauSq/x$sigmavSq),
    digits = digits, format = "f"), "\n")
  cat("Var[u]/{Var[u]+Var[v]}        = ", paste0(rep(" ", lengthSum -
    nchar("Var[u]/{Var[u]+Var[v]}        = ") - nchar(formatC(x$Varu/(x$Varu +
    x$sigmavSq), digits = digits, format = "f"))), collapse = ""),
    formatC(x$Varu/(x$Varu + x$sigmavSq), digits = digits,
      format = "f"), "\n")
  cat("Var[e]                        = ", paste0(rep(" ", lengthSum -
    nchar("Var[e]                        = ") - nchar(formatC(x$Varu +
    x$sigmavSq, digits = digits, format = "f"))), collapse = ""),
    formatC(x$Varu + x$sigmavSq, digits = digits, format = "f"),
    "\n")
  if (x$udist == "uniform") {
    cat("THETA                         = ", paste0(rep(" ",
      lengthSum - nchar("THETA                         = ") -
        nchar(formatC(x$theta, digits = digits, format = "f"))),
      collapse = ""), formatC(x$theta, digits = digits,
      format = "f"), "\n")
  }
  if (x$nuZUvar > 1 || x$nvZVvar > 1) {
    cat("Variances averaged over observations \n")
  }
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat("Average inefficiency E[u]     = ", paste0(rep(" ", lengthSum -
    nchar("Average inefficiency E[u]     = ") - nchar(formatC(x$Eu,
    digits = digits, format = "f"))), collapse = ""), formatC(x$Eu,
    digits = digits, format = "f"), "\n")
  cat("Average efficiency E[exp(-u)] = ", paste0(rep(" ", lengthSum -
    nchar("Average efficiency E[exp(-u)] = ") - nchar(formatC(x$Expu,
    digits = digits, format = "f"))), collapse = ""), formatC(x$Expu,
    digits = digits, format = "f"), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("-----[ Tests vs. No Inefficiency ]-----\n")
  cat("Likelihood Ratio Test of Inefficiency\n")
  cat("Deg. freedom for inefficiency model", paste0(rep(" ",
    lengthSum - nchar("Deg. freedom for inefficiency model") -
      nchar(formatC(x$df, digits = digits, format = "d"))),
    collapse = ""), formatC(x$df, digits = digits, format = "d"),
    "\n")
  cat("Log Likelihood for OLS Log(H0) = ", paste0(rep(" ",
    lengthSum - nchar("Log Likelihood for OLS Log(H0) = ") -
      nchar(formatC(x$olsLoglik, digits = digits, format = "f"))),
    collapse = ""), formatC(x$olsLoglik, digits = digits,
    format = "f"), "\n")
  cat("LR statistic: ", "\n")
  cat("Chisq = 2*[LogL(H0)-LogL(H1)]  = ", paste0(rep(" ",
    lengthSum - nchar("Chisq = 2*[LogL(H0)-LogL(H1)]  = ") -
      nchar(formatC(x$chisq, digits = digits, format = "f"))),
    collapse = ""), formatC(x$chisq, digits = digits, format = "f"),
    "\n")
  cat("Kodde-Palm C*:       95%:", formatC(qchibarsq(0.95,
    df = x$df), digits = digits, format = "f"), paste0(rep(" ",
    lengthSum - nchar("Kodde-Palm C*:       95%:") - nchar(formatC(qchibarsq(0.95,
      df = x$df), digits = digits, format = "f")) - nchar(formatC(qchibarsq(0.99,
      df = x$df), digits = digits, format = "f")) - nchar("99%") -
      3), collapse = ""), "99%:", formatC(qchibarsq(0.99,
    df = x$df), digits = digits, format = "f"), "\n")
  cat("Coelli (1995) skewness test on OLS residuals\n")
  cat("M3T                            = ", paste0(rep(" ",
    lengthSum - nchar("M3T                            = ") -
      nchar(formatC(x$CoelliM3Test[1], digits = digits,
        format = "f"))), collapse = ""), formatC(x$CoelliM3Test[1],
    digits = digits, format = "f"), "\n")
  cat("final maximum likelihood estimates \n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Deterministic Component of SFA", width = lengthSum +
    2 + switch(dimCoefTable, `4` = 18, `5` = 31, `6` = 43,
    `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes1, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  if (x$udist == "tnormal") {
    if (x$scaling) {
      cat(centerText("Scaling property parameters", width = lengthSum +
        2 + switch(dimCoefTable, `4` = 18, `5` = 31,
        `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes2, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes3, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes4, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes5, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
    } else {
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes2, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes3, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes4, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
    }
  } else {
    if (x$udist == "lognormal") {
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes2, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes3, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes4, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
    } else {
      if (x$udist == "gamma") {
        cat(centerText("Parameter in variance of u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Location parameter P in u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
      } else {
        if (x$udist == "weibull") {
          cat(centerText("Parameter in variance of u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Shape parameter k in u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
        } else {
          if (x$udist == "tslaplace") {
          cat(centerText("Parameter in variance of u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Skewness parameter 'lambda' in u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          } else {
          cat(centerText("Parameter in variance of u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          }
        }
      }
    }
  }
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  invisible(x)
}

# summary for lcmcross ----------
#' @rdname summary
#' @aliases summary.lcmcross
#' @export
summary.lcmcross <- function(object, grad = FALSE, ci = FALSE,
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
  class(object) <- "summary.lcmcross"
  return(object)
}

# print summary for lcmcross ----------
#' @rdname summary
#' @aliases print.summary.lcmcross
#' @export
print.summary.lcmcross <- function(x, digits = max(3, getOption("digits") -
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
  lengthSum <- nchar(sfaModel)  # + 10
  dimCoefTable <- as.character(dim(x$mlRes)[2])
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(sfaModel, "\n")
  cat("Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
    nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n")
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
    collapse = ""), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mlLoglik,
    digits = digits, format = "f"))), collapse = ""), formatC(x$mlLoglik,
    digits = digits, format = "f"), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
    digits = digits, format = "e"))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"), "\n")
  cat("Estimation based on:", if (lengthSum - nchar("Estimation based on:") -
    nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("and K = ") -
    3 > 0)
    paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") -
      nchar("and K = ") - 3), collapse = ""), "N = ", x$Nobs,
    "and K = ", x$nParm, "\n")
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC/x$Nobs, digits = 3,
    format = "f")) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"), "AIC/N  = ", formatC(x$AIC/x$Nobs,
    digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
    digits = 1, format = "f")) - nchar("BIC/N  = ") - nchar(formatC(x$BIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
    digits = 1, format = "f")) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("Latent class model with", x$nClasses, "latent classes \n")
  cat("final maximum likelihood estimates \n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Deterministic Component of SFA for latent class 1",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes[1:x$nXvar, , drop = FALSE], P.values = TRUE,
    digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Parameter in variance of u (one-sided error) for latent class 1",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar), ,
    drop = FALSE], P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Parameters in variance of v (two-sided error) for latent class 1",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar + x$nuZUvar +
    x$nvZVvar), , drop = FALSE], P.values = TRUE, digits = digits,
    signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Deterministic Component of SFA for latent class 2",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes[(x$nXvar + x$nuZUvar + x$nvZVvar + 1):(2 *
    x$nXvar + x$nuZUvar + x$nvZVvar), , drop = FALSE], P.values = TRUE,
    digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Parameter in variance of u (one-sided error) for latent class 2",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes[(2 * x$nXvar + x$nuZUvar + x$nvZVvar +
    1):(2 * x$nXvar + 2 * x$nuZUvar + x$nvZVvar), , drop = FALSE],
    P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Parameters in variance of v (two-sided error) for latent class 2",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar + x$nvZVvar +
    1):(2 * x$nXvar + 2 * x$nuZUvar + 2 * x$nvZVvar), , drop = FALSE],
    P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  if (x$nClasses == 2) {
    cat(centerText("Estimated prior probabilities for class membership",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57)), "\n")
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar + 2 *
      x$nvZVvar + 1):(2 * x$nXvar + 2 * x$nuZUvar + 2 *
      x$nvZVvar + x$nZHvar), , drop = FALSE], P.values = TRUE,
      digits = digits, signif.legend = TRUE)
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
  } else {
    if (x$nClasses == 3) {
      cat(centerText("Deterministic Component of SFA for latent class 3",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar + 1):(3 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
        digits = digits, signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameter in variance of u (one-sided error) for latent class 3",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes[(3 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
        digits = digits, signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameters in variance of v (two-sided error) for latent class 3",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
        2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
        digits = digits, signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Estimated prior probabilities for class membership",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar + 2 * x$nZHvar), , drop = FALSE],
        P.values = TRUE, digits = digits, signif.legend = TRUE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
    } else {
      if (x$nClasses == 4) {
        cat(centerText("Deterministic Component of SFA for latent class 3",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameter in variance of u (one-sided error) for latent class 3",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes[(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error) for latent class 3",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Deterministic Component of SFA for latent class 4",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameter in variance of u (one-sided error) for latent class 4",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes[(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error) for latent class 4",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Estimated prior probabilities for class membership",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 3 * x$nZHvar), , drop = FALSE],
          P.values = TRUE, digits = digits, signif.legend = TRUE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
      } else {
        if (x$nClasses == 5) {
          cat(centerText("Deterministic Component of SFA for latent class 3",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameter in variance of u (one-sided error) for latent class 3",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error) for latent class 3",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Deterministic Component of SFA for latent class 4",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameter in variance of u (one-sided error) for latent class 4",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error) for latent class 4",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Deterministic Component of SFA for latent class 5",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(5 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameter in variance of u (one-sided error) for latent class 5",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(5 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
          4 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error) for latent class 5",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(5 * x$nXvar + 5 * x$nuZUvar +
          4 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
          5 * x$nvZVvar), , drop = FALSE], P.values = TRUE,
          digits = digits, signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for class membership",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes[(5 * x$nXvar + 5 * x$nuZUvar +
          5 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
          5 * x$nvZVvar + 4 * x$nZHvar), , drop = FALSE],
          P.values = TRUE, digits = digits, signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
        }
      }
    }
  }
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
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
  uHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable$selectDum ==
    1, ], rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable$selectDum ==
    1, ], rhs = 3)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  object$sigmavSq <- mean(exp(Wv))
  object$sigmauSq <- mean(exp(Wu))
  object$Varu <- varuFun(object = object, mu = NULL, P = NULL,
    k = NULL, lambda = NULL)
  object$Eu <- euFun(object = object, mu = NULL, P = NULL,
    k = NULL, lambda = NULL)
  object$Expu <- eExpuFun(object = object, mu = NULL, P = NULL,
    k = NULL, lambda = NULL)
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
  object$df <- object$nParm - object$nClasses * object$nXvar -
    object$nClasses * object$nvZVvar - object$nZHvar * (object$nClasses -
    1)
  class(object) <- "summary.selectioncross"
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
  lengthSum <- nchar(sfaModel)
  dimCoefTable <- as.character(dim(x$mlRes)[2])
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(sfaModel, "\n")
  cat("Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
    nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n")
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
    collapse = ""), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mlLoglik,
    digits = digits, format = "f"))), collapse = ""), formatC(x$mlLoglik,
    digits = digits, format = "f"), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
    digits = digits, format = "e"))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"), "\n")
  cat("Estimation based on:", if (lengthSum - nchar("Estimation based on:") -
    nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("of") -
    nchar(x$Ninit) - nchar("obs.") - nchar("and K = ") -
    6 > 0)
    paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") -
      nchar("of") - nchar(x$Ninit) - nchar("obs.") - nchar("and K = ") -
      6), collapse = ""), "N = ", x$Nobs, "of", x$Ninit,
    "obs.", "and K = ", x$nParm, "\n")
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC/x$Nobs, digits = 3,
    format = "f")) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"), "AIC/N  = ", formatC(x$AIC/x$Nobs,
    digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
    digits = 1, format = "f")) - nchar("BIC/N  = ") - nchar(formatC(x$BIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
    digits = 1, format = "f")) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat("Variances: Sigma-squared(v)   = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(v)   = ") - nchar(formatC(x$sigmavSq,
    digits = digits, format = "f"))), collapse = ""), formatC(x$sigmavSq,
    digits = digits, format = "f"), "\n")
  cat("           Sigma(v)           = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(v)   = ") - nchar(formatC(sqrt(x$sigmavSq),
    digits = digits, format = "f"))), collapse = ""), formatC(sqrt(x$sigmavSq),
    digits = digits, format = "f"), "\n")
  cat("           Sigma-squared(u)   = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(u)   = ") - nchar(formatC(x$sigmauSq,
    digits = digits, format = "f"))), collapse = ""), formatC(x$sigmauSq,
    digits = digits, format = "f"), "\n")
  cat("           Sigma(u)           = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(u)   = ") - nchar(formatC(sqrt(x$sigmauSq),
    digits = digits, format = "f"))), collapse = ""), formatC(sqrt(x$sigmauSq),
    digits = digits, format = "f"), "\n")
  cat("Sigma = Sqrt[(s^2(u)+s^2(v))] = ", paste0(rep(" ", lengthSum -
    nchar("Sigma = Sqrt[(s^2(u)+s^2(v))] = ") - nchar(formatC(sqrt(x$sigmavSq +
    x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    formatC(sqrt(x$sigmavSq + x$sigmauSq), digits = digits,
      format = "f"), "\n")
  cat("Gamma = sigma(u)^2/sigma^2    = ", paste0(rep(" ", lengthSum -
    nchar("Gamma = sigma(u)^2/sigma^2    = ") - nchar(formatC(x$sigmauSq/(x$sigmavSq +
    x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    formatC(x$sigmauSq/(x$sigmavSq + x$sigmauSq), digits = digits,
      format = "f"), "\n")
  cat("Lambda = sigma(u)/sigma(v)    = ", paste0(rep(" ", lengthSum -
    nchar("Lambda = sigma(u)/sigma(v)    = ") - nchar(formatC(sqrt(x$sigmauSq/x$sigmavSq),
    digits = digits, format = "f"))), collapse = ""), formatC(sqrt(x$sigmauSq/x$sigmavSq),
    digits = digits, format = "f"), "\n")
  cat("Var[u]/{Var[u]+Var[v]}        = ", paste0(rep(" ", lengthSum -
    nchar("Var[u]/{Var[u]+Var[v]}        = ") - nchar(formatC(x$Varu/(x$Varu +
    x$sigmavSq), digits = digits, format = "f"))), collapse = ""),
    formatC(x$Varu/(x$Varu + x$sigmavSq), digits = digits,
      format = "f"), "\n")
  cat("Var[e]                        = ", paste0(rep(" ", lengthSum -
    nchar("Var[e]                        = ") - nchar(formatC(x$Varu +
    x$sigmavSq, digits = digits, format = "f"))), collapse = ""),
    formatC(x$Varu + x$sigmavSq, digits = digits, format = "f"),
    "\n")
  if (x$nuZUvar > 1 || x$nvZVvar > 1) {
    cat("Variances averaged over observations \n")
  }
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat("Average inefficiency E[u]     = ", paste0(rep(" ", lengthSum -
    nchar("Average inefficiency E[u]     = ") - nchar(formatC(x$Eu,
    digits = digits, format = "f"))), collapse = ""), formatC(x$Eu,
    digits = digits, format = "f"), "\n")
  cat("Average efficiency E[exp(-u)] = ", paste0(rep(" ", lengthSum -
    nchar("Average efficiency E[exp(-u)] = ") - nchar(formatC(x$Expu,
    digits = digits, format = "f"))), collapse = ""), formatC(x$Expu,
    digits = digits, format = "f"), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("Estimator is 2 step Maximum Likelihood \n")
  cat("final maximum likelihood estimates \n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Deterministic Component of SFA", width = lengthSum +
    2 + switch(dimCoefTable, `4` = 18, `5` = 31, `6` = 43,
    `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes1, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Parameter in variance of u (one-sided error)",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes2, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Parameters in variance of v (two-sided error)",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes3, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Selection bias parameter", width = lengthSum +
    2 + switch(dimCoefTable, `4` = 18, `5` = 31, `6` = 43,
    `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes4, P.values = TRUE, digits = digits, signif.legend = TRUE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  invisible(x)
}

# summary for zisfcross ----------
#' @rdname summary
#' @aliases summary.zisfcross
#' @export
summary.zisfcross <- function(object, grad = FALSE, ci = FALSE,
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
  object$chisq <- 2 * (object$mlLoglik - object$olsLoglik)
  object$df <- object$nParm - object$nXvar - object$nvZVvar
  class(object) <- "summary.zisfcross"
  return(object)
}

# print summary for zisfcross ----------
#' @rdname summary
#' @aliases print.summary.zisfcross
#' @export
print.summary.zisfcross <- function(x, digits = max(3, getOption("digits") -
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
  if (x$sigmavType == "common") {
    if (x$udist %in% c("lognormal", "tnormal")) {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nmuZUvar),
        , drop = FALSE]
      mlRes3 <- mlRes[(x$nXvar + x$nmuZUvar + 1):(x$nXvar +
        x$nmuZUvar + x$nuZUvar), , drop = FALSE]
      mlRes4 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
        1):(x$nXvar + x$nmuZUvar + x$nuZUvar + x$nvZVvar),
        , drop = FALSE]
      mlRes5 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
        x$nvZVvar + 1):(x$nXvar + x$nmuZUvar + x$nuZUvar +
        x$nvZVvar + x$nZHvar), , drop = FALSE]
    } else {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar),
        , drop = FALSE]
      mlRes3 <- mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar +
        x$nuZUvar + x$nvZVvar), , drop = FALSE]
      if (x$udist %in% c("gamma", "weibull", "tslaplace")) {
        mlRes4 <- mlRes[x$nXvar + x$nuZUvar + x$nvZVvar +
          1, , drop = FALSE]
        mlRes5 <- mlRes[(x$nXvar + x$nuZUvar + x$nvZVvar +
          2):(x$nXvar + x$nuZUvar + x$nvZVvar + x$nZHvar +
          1), , drop = FALSE]
      } else {
        mlRes4 <- mlRes[(x$nXvar + x$nuZUvar + x$nvZVvar +
          1):(x$nXvar + x$nuZUvar + x$nvZVvar + x$nZHvar),
          , drop = FALSE]
      }
    }
  } else {
    if (x$sigmavType == "different") {
      if (x$udist %in% c("lognormal", "tnormal")) {
        mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nmuZUvar),
          , drop = FALSE]
        mlRes3 <- mlRes[(x$nXvar + x$nmuZUvar + 1):(x$nXvar +
          x$nmuZUvar + x$nuZUvar), , drop = FALSE]
        mlRes4 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
          1):(x$nXvar + x$nmuZUvar + x$nuZUvar + x$nvZVvar),
          , drop = FALSE]
        mlRes5 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
          x$nvZVvar + 1):(x$nXvar + x$nmuZUvar + x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE]
        mlRes6 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
          2 * x$nvZVvar + 1):(x$nXvar + x$nmuZUvar +
          x$nuZUvar + 2 * x$nvZVvar + x$nZHvar), , drop = FALSE]
      } else {
        mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar),
          , drop = FALSE]
        mlRes3 <- mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar +
          x$nuZUvar + x$nvZVvar), , drop = FALSE]
        mlRes4 <- mlRes[(x$nXvar + x$nuZUvar + x$nvZVvar +
          1):(x$nXvar + x$nuZUvar + 2 * x$nvZVvar), ,
          drop = FALSE]
        if (x$udist %in% c("gamma", "weibull", "tslaplace")) {
          mlRes5 <- mlRes[x$nXvar + x$nuZUvar + 2 * x$nvZVvar +
          1, , drop = FALSE]
          mlRes6 <- mlRes[(x$nXvar + x$nuZUvar + 2 *
          x$nvZVvar + 2):(x$nXvar + x$nuZUvar + 2 *
          x$nvZVvar + x$nZHvar + 1), , drop = FALSE]
        } else {
          mlRes5 <- mlRes[(x$nXvar + x$nuZUvar + 2 *
          x$nvZVvar + 1):(x$nXvar + x$nuZUvar + 2 *
          x$nvZVvar + x$nZHvar), , drop = FALSE]
        }
      }
    }
  }
  lengthSum <- nchar(zisfdist(x$udist)) + 27
  dimCoefTable <- as.character(dim(x$mlRes)[2])
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(zisfdist(x$udist), "\n")
  cat("Class Membership link function:", paste0(rep(" ", lengthSum -
    nchar("Class Membership link function:") - nchar(x$linkF)),
    collapse = ""), x$linkF, "\n")
  cat("Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
    nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n")
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
    collapse = ""), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mlLoglik,
    digits = digits, format = "f"))), collapse = ""), formatC(x$mlLoglik,
    digits = digits, format = "f"), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
    digits = digits, format = "e"))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"), "\n")
  cat("Estimation based on:", if (lengthSum - nchar("Estimation based on:") -
    nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("and K = ") -
    3 > 0)
    paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") -
      nchar("and K = ") - 3), collapse = ""), "N = ", x$Nobs,
    "and K = ", x$nParm, "\n")
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC/x$Nobs, digits = 3,
    format = "f")) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"), "AIC/N  = ", formatC(x$AIC/x$Nobs,
    digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
    digits = 1, format = "f")) - nchar("BIC/N  = ") - nchar(formatC(x$BIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
    digits = 1, format = "f")) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("final maximum likelihood estimates \n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Deterministic Component of SFA", width = lengthSum +
    2 + switch(dimCoefTable, `4` = 18, `5` = 31, `6` = 43,
    `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes1, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  if (x$udist %in% c("lognormal", "tnormal")) {
    cat(centerText("Location parameter [offset mu] in u (one-sided error)",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57)), "\n")
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    printCoefmat(mlRes2, P.values = TRUE, digits = digits,
      signif.legend = FALSE)
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    cat(centerText("Parameter in variance of u (one-sided error)",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57)), "\n")
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    printCoefmat(mlRes3, P.values = TRUE, digits = digits,
      signif.legend = FALSE)
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    if (x$sigmavType == "common") {
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes4, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Estimated prior probabilities for class membership",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes5, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
    } else {
      if (x$sigmavType == "different") {
        cat(centerText("Parameters in variance of v (two-sided error) for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error) for efficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Estimated prior probabilities for class membership",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes6, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
      }
    }
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
  } else {
    if (x$udist == "gamma") {
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes2, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      if (x$sigmavType == "common") {
        cat(centerText("Parameters in variance of v (two-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Location parameter P in u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Estimated prior probabilities for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
      } else {
        if (x$sigmavType == "different") {
          cat(centerText("Parameters in variance of v (two-sided error) for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error) for efficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Location parameter P in u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes6, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        }
      }
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
    } else {
      if (x$udist == "weibull") {
        cat(centerText("Parameter in variance of u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        if (x$sigmavType == "common") {
          cat(centerText("Parameters in variance of v (two-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Shape parameter k in u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        } else {
          if (x$sigmavType == "different") {
          cat(centerText("Parameters in variance of v (two-sided error) for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error) for efficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Shape parameter k in u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes6, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          }
        }
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
      } else {
        if (x$udist == "tslaplace") {
          cat(centerText("Parameter in variance of u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          if (x$sigmavType == "common") {
          cat(centerText("Parameters in variance of v (two-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Skewness parameter 'lambda' in u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          } else {
          if (x$sigmavType == "different") {
            cat(centerText("Parameters in variance of v (two-sided error) for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Parameters in variance of v (two-sided error) for efficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Skewness parameter 'lambda' in u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            printCoefmat(mlRes5, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes6, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          }
          }
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
        } else {
          cat(centerText("Parameter in variance of u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          if (x$sigmavType == "common") {
          cat(centerText("Parameters in variance of v (two-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          } else {
          if (x$sigmavType == "different") {
            cat(centerText("Parameters in variance of v (two-sided error) for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Parameters in variance of v (two-sided error) for efficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes5, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          }
          }
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
        }
      }
    }
  }
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  invisible(x)
}

# summary for cnsfcross ----------
#' @rdname summary
#' @aliases summary.cnsfcross
#' @export
summary.cnsfcross <- function(object, grad = FALSE, ci = FALSE,
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
  object$chisq <- 2 * (object$mlLoglik - object$olsLoglik)
  object$df <- object$nParm - object$nXvar - object$nvZVvar
  class(object) <- "summary.cnsfcross"
  return(object)
}

# print summary for cnsfcross ----------
#' @rdname summary
#' @aliases print.summary.cnsfcross
#' @export
print.summary.cnsfcross <- function(x, digits = max(3, getOption("digits") -
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
  if (x$sigmauType == "common") {
    if (x$udist %in% c("lognormal", "tnormal")) {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nmuZUvar),
        , drop = FALSE]
      mlRes3 <- mlRes[(x$nXvar + x$nmuZUvar + 1):(x$nXvar +
        x$nmuZUvar + x$nuZUvar), , drop = FALSE]
      mlRes4 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
        1):(x$nXvar + x$nmuZUvar + x$nuZUvar + x$nvZVvar),
        , drop = FALSE]
      mlRes5 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
        x$nvZVvar + 1):(x$nXvar + x$nmuZUvar + x$nuZUvar +
        2 * x$nvZVvar), , drop = FALSE]
      mlRes6 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
        2 * x$nvZVvar + 1):(x$nXvar + x$nmuZUvar + x$nuZUvar +
        2 * x$nvZVvar + x$nZHvar), , drop = FALSE]
    } else {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar),
        , drop = FALSE]
      mlRes3 <- mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar +
        x$nuZUvar + x$nvZVvar), , drop = FALSE]
      mlRes4 <- mlRes[(x$nXvar + x$nuZUvar + x$nvZVvar +
        1):(x$nXvar + x$nuZUvar + 2 * x$nvZVvar), , drop = FALSE]
      if (x$udist %in% c("gamma", "weibull", "tslaplace")) {
        mlRes5 <- mlRes[x$nXvar + x$nuZUvar + 2 * x$nvZVvar +
          1, , drop = FALSE]
        mlRes6 <- mlRes[(x$nXvar + x$nuZUvar + 2 * x$nvZVvar +
          2):(x$nXvar + x$nuZUvar + 2 * x$nvZVvar + x$nZHvar +
          1), , drop = FALSE]
      } else {
        mlRes5 <- mlRes[(x$nXvar + x$nuZUvar + 2 * x$nvZVvar +
          1):(x$nXvar + x$nuZUvar + 2 * x$nvZVvar + x$nZHvar),
          , drop = FALSE]
      }
    }
  } else {
    if (x$sigmauType == "different") {
      if (x$udist %in% c("lognormal", "tnormal")) {
        mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nmuZUvar),
          , drop = FALSE]
        mlRes3 <- mlRes[(x$nXvar + x$nmuZUvar + 1):(x$nXvar +
          2 * x$nmuZUvar), , drop = FALSE]
        mlRes4 <- mlRes[(x$nXvar + 2 * x$nmuZUvar + 1):(x$nXvar +
          2 * x$nmuZUvar + x$nuZUvar), , drop = FALSE]
        mlRes5 <- mlRes[(x$nXvar + 2 * x$nmuZUvar + x$nuZUvar +
          1):(x$nXvar + 2 * x$nmuZUvar + 2 * x$nuZUvar),
          , drop = FALSE]
        mlRes6 <- mlRes[(x$nXvar + 2 * x$nmuZUvar + 2 *
          x$nuZUvar + 1):(x$nXvar + 2 * x$nmuZUvar +
          2 * x$nuZUvar + x$nvZVvar), , drop = FALSE]
        mlRes7 <- mlRes[(x$nXvar + 2 * x$nmuZUvar + 2 *
          x$nuZUvar + x$nvZVvar + 1):(x$nXvar + 2 * x$nmuZUvar +
          2 * x$nuZUvar + 2 * x$nvZVvar), , drop = FALSE]
        mlRes8 <- mlRes[(x$nXvar + 2 * x$nmuZUvar + 2 *
          x$nuZUvar + 2 * x$nvZVvar + 1):(x$nXvar + 2 *
          x$nmuZUvar + 2 * x$nuZUvar + 2 * x$nvZVvar +
          x$nZHvar), , drop = FALSE]
      } else {
        mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar),
          , drop = FALSE]
        mlRes3 <- mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar +
          2 * x$nuZUvar), , drop = FALSE]
        mlRes4 <- mlRes[(x$nXvar + 2 * x$nuZUvar + 1):(x$nXvar +
          2 * x$nuZUvar + x$nvZVvar), , drop = FALSE]
        mlRes5 <- mlRes[(x$nXvar + 2 * x$nuZUvar + x$nvZVvar +
          1):(x$nXvar + 2 * x$nuZUvar + 2 * x$nvZVvar),
          , drop = FALSE]
        if (x$udist %in% c("gamma", "weibull", "tslaplace")) {
          mlRes6 <- mlRes[x$nXvar + 2 * x$nuZUvar + 2 *
          x$nvZVvar + 1, , drop = FALSE]
          mlRes7 <- mlRes[x$nXvar + 2 * x$nuZUvar + 2 *
          x$nvZVvar + 2, , drop = FALSE]
          mlRes8 <- mlRes[(x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 3):(x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + x$nZHvar + 2), , drop = FALSE]
        } else {
          mlRes6 <- mlRes[(x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + x$nZHvar), , drop = FALSE]
        }
      }
    }
  }
  lengthSum <- nchar(cnsfdist(x$udist)) + 27
  dimCoefTable <- as.character(dim(x$mlRes)[2])
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(cnsfdist(x$udist), "\n")
  cat("Class Membership link function:", paste0(rep(" ", lengthSum -
    nchar("Class Membership link function:") - nchar(x$linkF)),
    collapse = ""), x$linkF, "\n")
  cat("Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
    nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n")
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
    collapse = ""), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mlLoglik,
    digits = digits, format = "f"))), collapse = ""), formatC(x$mlLoglik,
    digits = digits, format = "f"), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
    digits = digits, format = "e"))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"), "\n")
  cat("Estimation based on:", if (lengthSum - nchar("Estimation based on:") -
    nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("and K = ") -
    3 > 0)
    paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") -
      nchar("and K = ") - 3), collapse = ""), "N = ", x$Nobs,
    "and K = ", x$nParm, "\n")
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC/x$Nobs, digits = 3,
    format = "f")) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"), "AIC/N  = ", formatC(x$AIC/x$Nobs,
    digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
    digits = 1, format = "f")) - nchar("BIC/N  = ") - nchar(formatC(x$BIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
    digits = 1, format = "f")) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("final maximum likelihood estimates \n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Deterministic Component of SFA", width = lengthSum +
    2 + switch(dimCoefTable, `4` = 18, `5` = 31, `6` = 43,
    `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes1, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  if (x$udist %in% c("lognormal", "tnormal")) {
    if (x$sigmauType == "common") {
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes2, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes3, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameters in variance of v (two-sided error): class 1",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes4, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameters in variance of v (two-sided error): class 2",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes5, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Estimated prior probabilities for class membership",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes6, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
    } else {
      if (x$sigmauType == "different") {
        cat(centerText("Location parameter [offset mu] in u (one-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Location parameter [offset mu] in u (one-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameter in variance of u (one-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameter in variance of u (one-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes6, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes7, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Estimated prior probabilities for class membership",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes8, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
      }
    }
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
  } else {
    if (x$udist == "gamma") {
      if (x$sigmauType == "common") {
        cat(centerText("Parameter in variance of u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Location parameter P in u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Estimated prior probabilities for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes6, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
      } else {
        if (x$sigmauType == "different") {
          cat(centerText("Parameter in variance of u (one-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameter in variance of u (one-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Location parameter P in u (one-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes6, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(centerText("Location parameter P in u (one-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes7, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes8, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        }
      }
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
    } else {
      if (x$udist == "weibull") {
        if (x$sigmauType == "common") {
          cat(centerText("Parameter in variance of u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Shape parameter k in u (one-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes6, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        } else {
          if (x$sigmauType == "different") {
          cat(centerText("Parameter in variance of u (one-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameter in variance of u (one-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Shape parameter k in u (one-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes6, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Shape parameter k in u (one-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes7, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes8, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          }
        }
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
      } else {
        if (x$udist == "tslaplace") {
          if (x$sigmauType == "common") {
          cat(centerText("Parameter in variance of u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Skewness parameter 'lambda' in u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes6, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          } else {
          if (x$sigmauType == "different") {
            cat(centerText("Parameter in variance of u (one-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes2, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Parameter in variance of u (one-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Parameters in variance of v (two-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Parameters in variance of v (two-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes5, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Skewness parameter 'lambda' in u (one-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            printCoefmat(mlRes6, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Skewness parameter 'lambda' in u (one-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            printCoefmat(mlRes7, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes8, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          }
          }
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
        } else {
          if (x$sigmauType == "common") {
          cat(centerText("Parameter in variance of u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          } else {
          if (x$sigmauType == "different") {
            cat(centerText("Parameter in variance of u (one-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes2, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Parameter in variance of u (one-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes3, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Parameters in variance of v (two-sided error): class 1",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes4, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Parameters in variance of v (two-sided error): class 2",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes5, P.values = TRUE, digits = digits,
            signif.legend = FALSE)
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            cat(centerText("Estimated prior probabilities for inefficient class",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            "\n")
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
            collapse = ""), "\n")
            printCoefmat(mlRes6, P.values = TRUE, digits = digits,
            signif.legend = TRUE)
          }
          }
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
        }
      }
    }
  }
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  invisible(x)
}

# summary for misfcross ----------
#' @rdname summary
#' @aliases summary.misfcross
#' @export
summary.misfcross <- function(object, grad = FALSE, ci = FALSE,
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
  object$chisq <- 2 * (object$mlLoglik - object$olsLoglik)
  object$df <- object$nParm - object$nXvar - object$nvZVvar
  class(object) <- "summary.misfcross"
  return(object)
}

# print summary for misfcross ----------
#' @rdname summary
#' @aliases print.summary.misfcross
#' @export
print.summary.misfcross <- function(x, digits = max(3, getOption("digits") -
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
  if (x$udist %in% c("lognormal", "tnormal")) {
    mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nmuZUvar),
      , drop = FALSE]
    mlRes3 <- mlRes[(x$nXvar + x$nmuZUvar + 1):(x$nXvar +
      2 * x$nmuZUvar), , drop = FALSE]
    mlRes4 <- mlRes[(x$nXvar + 2 * x$nmuZUvar + 1):(x$nXvar +
      2 * x$nmuZUvar + x$nuZUvar), , drop = FALSE]
    mlRes5 <- mlRes[(x$nXvar + 2 * x$nmuZUvar + x$nuZUvar +
      1):(x$nXvar + 2 * x$nmuZUvar + 2 * x$nuZUvar), ,
      drop = FALSE]
    mlRes6 <- mlRes[(x$nXvar + 2 * x$nmuZUvar + 2 * x$nuZUvar +
      1):(x$nXvar + 2 * x$nmuZUvar + 2 * x$nuZUvar + x$nvZVvar),
      , drop = FALSE]
    mlRes7 <- mlRes[(x$nXvar + 2 * x$nmuZUvar + 2 * x$nuZUvar +
      x$nvZVvar + 1):(x$nXvar + 2 * x$nmuZUvar + 2 * x$nuZUvar +
      x$nvZVvar + x$nZHvar), , drop = FALSE]
  } else {
    mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar),
      , drop = FALSE]
    mlRes3 <- mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar +
      2 * x$nuZUvar), , drop = FALSE]
    mlRes4 <- mlRes[(x$nXvar + 2 * x$nuZUvar + 1):(x$nXvar +
      2 * x$nuZUvar + x$nvZVvar), , drop = FALSE]
    if (x$udist %in% c("gamma", "weibull", "tslaplace")) {
      mlRes5 <- mlRes[x$nXvar + 2 * x$nuZUvar + x$nvZVvar +
        1, , drop = FALSE]
      mlRes6 <- mlRes[x$nXvar + 2 * x$nuZUvar + x$nvZVvar +
        2, , drop = FALSE]
      mlRes7 <- mlRes[(x$nXvar + 2 * x$nuZUvar + x$nvZVvar +
        3):(x$nXvar + 2 * x$nuZUvar + x$nvZVvar + x$nZHvar +
        2), , drop = FALSE]
    } else {
      mlRes5 <- mlRes[(x$nXvar + 2 * x$nuZUvar + x$nvZVvar +
        1):(x$nXvar + 2 * x$nuZUvar + x$nvZVvar +
        x$nZHvar), , drop = FALSE]
    }
  }
  lengthSum <- nchar(misfdist(x$udist)) + 27
  dimCoefTable <- as.character(dim(x$mlRes)[2])
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(misfdist(x$udist), "\n")
  cat("Class Membership link function:", paste0(rep(" ", lengthSum -
    nchar("Class Membership link function:") - nchar(x$linkF)),
    collapse = ""), x$linkF, "\n")
  cat("Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
    nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n")
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
    collapse = ""), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mlLoglik,
    digits = digits, format = "f"))), collapse = ""), formatC(x$mlLoglik,
    digits = digits, format = "f"), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
    digits = digits, format = "e"))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"), "\n")
  cat("Estimation based on:", if (lengthSum - nchar("Estimation based on:") -
    nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("and K = ") -
    3 > 0)
    paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") -
      nchar("and K = ") - 3), collapse = ""), "N = ", x$Nobs,
    "and K = ", x$nParm, "\n")
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC/x$Nobs, digits = 3,
    format = "f")) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"), "AIC/N  = ", formatC(x$AIC/x$Nobs,
    digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
    digits = 1, format = "f")) - nchar("BIC/N  = ") - nchar(formatC(x$BIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
    digits = 1, format = "f")) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC/x$Nobs,
    digits = 3, format = "f")) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC/x$Nobs, digits = 3, format = "f"), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("final maximum likelihood estimates \n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  cat(centerText("Deterministic Component of SFA", width = lengthSum +
    2 + switch(dimCoefTable, `4` = 18, `5` = 31, `6` = 43,
    `7` = 57)), "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  printCoefmat(mlRes1, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  if (x$udist %in% c("lognormal", "tnormal")) {
    cat(centerText("Location parameter [offset mu] in u (one-sided error): class 1",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57)), "\n")
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    printCoefmat(mlRes2, P.values = TRUE, digits = digits,
      signif.legend = FALSE)
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    cat(centerText("Location parameter [offset mu] in u (one-sided error): class 2",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57)), "\n")
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    printCoefmat(mlRes3, P.values = TRUE, digits = digits,
      signif.legend = FALSE)
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    cat(centerText("Parameter in variance of u (one-sided error): class 1",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57)), "\n")
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    printCoefmat(mlRes4, P.values = TRUE, digits = digits,
      signif.legend = FALSE)
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    cat(centerText("Parameter in variance of u (one-sided error): class 2",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57)), "\n")
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    printCoefmat(mlRes5, P.values = TRUE, digits = digits,
      signif.legend = FALSE)
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    cat(centerText("Parameters in variance of v (two-sided error)",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57)), "\n")
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    printCoefmat(mlRes6, P.values = TRUE, digits = digits,
      signif.legend = FALSE)
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    cat(centerText("Estimated prior probabilities for class membership",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57)), "\n")
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
    printCoefmat(mlRes7, P.values = TRUE, digits = digits,
      signif.legend = TRUE)
    cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
      "\n")
  } else {
    if (x$udist == "gamma") {
      cat(centerText("Parameter in variance of u (one-sided error): class 1",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes2, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameter in variance of u (one-sided error): class 2",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes3, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes4, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Location parameter P in u (one-sided error): class 1",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes5, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(centerText("Location parameter P in u (one-sided error): class 2",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes6, P.values = TRUE, digits = digits,
        signif.legend = FALSE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      cat(centerText("Estimated prior probabilities for inefficient class",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), "\n")
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
      printCoefmat(mlRes7, P.values = TRUE, digits = digits,
        signif.legend = TRUE)
      cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
        "\n")
    } else {
      if (x$udist == "weibull") {
        cat(centerText("Parameter in variance of u (one-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameter in variance of u (one-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Parameters in variance of v (two-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Shape parameter k in u (one-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Shape parameter k in u (one-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes6, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        cat(centerText("Estimated prior probabilities for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
        printCoefmat(mlRes7, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
        cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
          "\n")
      } else {
        if (x$udist == "tslaplace") {
          cat(centerText("Parameter in variance of u (one-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameter in variance of u (one-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Skewness parameter 'lambda' in u (one-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Skewness parameter 'lambda' in u (one-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          printCoefmat(mlRes6, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes7, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
        } else {
          cat(centerText("Parameter in variance of u (one-sided error): class 1",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes2, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameter in variance of u (one-sided error): class 2",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes3, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Parameters in variance of v (two-sided error)",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes4, P.values = TRUE, digits = digits,
          signif.legend = FALSE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          cat(centerText("Estimated prior probabilities for inefficient class",
          width = lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          "\n")
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
          printCoefmat(mlRes5, P.values = TRUE, digits = digits,
          signif.legend = TRUE)
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57)),
          collapse = ""), "\n")
        }
      }
    }
  }
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
    `4` = 18, `5` = 31, `6` = 43, `7` = 57)), collapse = ""),
    "\n")
  invisible(x)
}
