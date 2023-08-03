################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Zero Inefficiency Stochastic Frontier Model                           #
# Two types:      - Common noise component (sigma_v)                           #
#                 - Different noise component (multimodal noise - mnsf)        #
# Link functions: - logit exp(theta * Z)/(1 + exp(theta * Z))                  #
#                 - cauchit 1/pi * atan(theta * Z) + 1/2                       #
#                 - probit pnorm(theta * Z)                                    #
#                 - cloglog 1 - exp(-exp(theta * Z))                           #
# Convolution: truncated skewed laplace - normal                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf truncated_skewed_laplace-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
# Same sigma_v

## logit specification class membership
czisftslnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  if (lambda < 0)
    return(NA)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cauchit specification class membership
czisftslnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  if (lambda < 0)
    return(NA)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## probit specification class membership
czisftslnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  if (lambda < 0)
    return(NA)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cloglog specification class membership
czisftslnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  if (lambda < 0)
    return(NA)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

# Different sigma_v

## logit specification class membership
cmnsftslnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  if (lambda < 0)
    return(NA)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cauchit specification class membership
cmnsftslnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  if (lambda < 0)
    return(NA)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## probit specification class membership
cmnsftslnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  if (lambda < 0)
    return(NA)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cloglog specification class membership
cmnsftslnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  if (lambda < 0)
    return(NA)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf truncated_skewed_laplace-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param nXvar number of main variables (inputs + env. var)
#' @param wHvar vector of weights (weighted likelihood)
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
# Same sigma_v
cstzisftslnorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- csttslnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initTSL <- NULL
  } else {
    cat("Initialization: SFA + truncated skewed laplace - normal distributions...\n")
    initTSL <- maxLik::maxLik(logLik = ctslnormlike, start = csttslnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradtslnormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initTSL$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), Esti[nXvar + 3],
    rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "lambda", paste0("ZISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initTSL = initTSL))
}

# Different sigma_v
cstmnsftslnorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- csttslnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initTSL <- NULL
  } else {
    cat("Initialization: SFA + truncated skewed laplace - normal distributions...\n")
    initTSL <- maxLik::maxLik(logLik = ctslnormlike, start = csttslnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradtslnormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initTSL$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar +
    2], if (nvZVvar > 1) rep(0, nvZVvar - 1), Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), "lambda", paste0("ZISF_",
    colnames(Zvar)))
  return(list(StartVal = StartVal, initTSL = initTSL))
}

# Gradient of the likelihood function ----------
#' gradient for zisf truncated_skewed_laplace-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
# Same sigma_v

## logit specification class membership
cgradzisftslnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  wvwu <- ewvsr/ewusr
  lwv <- ((1 + lambda) * wvwu + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (wvwu + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wzl <- (lwu * wzdeno)
  lwzu <- (1/wzl - lwu * ewz/wzl^2)
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (wvwu) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * wvwu) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * wvwu) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - (1 + lambda) * pd/ewusr) *
    eB)
  sigx2 <- ((1 + lambda) * sigx1 * ewz/wzl + S * prC * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (prC * dwsr/ewvsr + (1 + lambda) * eaeb * ewz/wzl)
  sigx4 <- ((0.5 * (da * wvwu) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - (0.5 * (dd * wvwu) - sigx5 * pd) * (1 + lambda) * eB)
  sigx7 <- (sigx6/wzl - lambda * wzdeno * eaeb * ewusr/wzl^2)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (0.5 * sigx8 - 0.5 * dwsr) * prC/ewvsr
  sigx10 <- (eA * pawv)
  sigx11 <- ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd)
  sigx12 <- (sigx9 + (1 + lambda) * (2 * sigx10 - sigx11 * eB) * ewz/wzl)
  sigx13 <- (2 * (eA * pa) - (llwusq + pd) * eB)
  sigx14 <- (sigx13/wzl - 2 * (wzdeno * (1 + lambda) * eaeb * ewusr/wzl^2))
  sigx15 <- ((1 + lambda) * lwzu * eaeb - prC * dwsr/(wzdeno * ewvsr))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx7 * (1 + lambda) * ewz/sigx3, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx12/sigx3, FUN = "*"), sigx14 * ewz/sigx3,
    sweep(Zvar, MARGIN = 1, STATS = sigx15 * ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradzisftslnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  lwv <- ((1 + lambda) * ewvsr/ewusr + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * ewvsr/ewusr) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - (1 + lambda) * pd/ewusr) *
    eB)
  sigx2 <- (ewz1 * (1 + lambda) * sigx1/lwu + S * ewz2 * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (ewz2 * dwsr/ewvsr + ewz1 * (1 + lambda) * eaeb/lwu)
  sigx4 <- ((0.5 * (da * ewvsr/ewusr) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - ((0.5 * (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) *
    eB + lambda * eaeb * ewusr/lwu))
  sigx7 <- (2 * (eA * pawv) - ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd) *
    eB)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (ewz2 * (0.5 * sigx8 - 0.5 * dwsr)/ewvsr + ewz1 * (1 + lambda) * sigx7/lwu)
  sigx10 <- (2 * (eA * pa) - ((llwusq + pd) * eB + 2 * ((1 + lambda) * eaeb * ewusr/lwu)))
  sigx11 <- ((1 + lambda) * eaeb/lwu - dwsr/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = ewz1 * (1 + lambda) * sigx6/(sigx3 * lwu),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx9/sigx3, FUN = "*"),
    ewz1 * sigx10/(sigx3 * lwu), sweep(Zvar, MARGIN = 1, STATS = sigx11/(pi *
      sigx3 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradzisftslnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  lwv <- ((1 + lambda) * ewvsr/ewusr + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * ewvsr/ewusr) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - (1 + lambda) * pd/ewusr) *
    eB)
  sigx2 <- ((1 + lambda) * sigx1 * pwZ/lwu + S * (1 - pwZ) * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - pwZ) * dwsr/ewvsr + (1 + lambda) * eaeb * pwZ/lwu)
  sigx4 <- ((0.5 * (da * ewvsr/ewusr) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - ((0.5 * (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) *
    eB + lambda * eaeb * ewusr/lwu))
  sigx7 <- (2 * (eA * pawv) - ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd) *
    eB)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- ((0.5 * sigx8 - 0.5 * dwsr) * (1 - pwZ)/ewvsr + (1 + lambda) * sigx7 *
    pwZ/lwu)
  sigx10 <- (2 * (eA * pa) - ((llwusq + pd) * eB + 2 * ((1 + lambda) * eaeb * ewusr/lwu)))
  sigx11 <- ((1 + lambda) * eaeb/lwu - dwsr/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (1 + lambda) * sigx6 * pwZ/(sigx3 * lwu),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx9/sigx3, FUN = "*"),
    sigx10 * pwZ/(sigx3 * lwu), sweep(Zvar, MARGIN = 1, STATS = sigx11 * dwZ/sigx3,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradzisftslnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  lwv <- ((1 + lambda) * ewvsr/ewusr + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * ewvsr/ewusr) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - (1 + lambda) * pd/ewusr) *
    eB)
  sigx2 <- ((1 - prZ) * (1 + lambda) * sigx1/lwu + S * dwsr * prZ * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - prZ) * (1 + lambda) * eaeb/lwu + dwsr * prZ/ewvsr)
  sigx4 <- ((0.5 * (da * ewvsr/ewusr) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - ((0.5 * (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) *
    eB + lambda * eaeb * ewusr/lwu))
  sigx7 <- (2 * (eA * pawv) - ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd) *
    eB)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- ((0.5 * sigx8 - 0.5 * dwsr) * prZ/ewvsr + (1 - prZ) * (1 + lambda) *
    sigx7/lwu)
  sigx10 <- (2 * (eA * pa) - ((llwusq + pd) * eB + 2 * ((1 + lambda) * eaeb * ewusr/lwu)))
  sigx11 <- ((1 + lambda) * eaeb/lwu - dwsr/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (1 - prZ) * (1 + lambda) * sigx6/(sigx3 *
      lwu), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx9/sigx3, FUN = "*"),
    (1 - prZ) * sigx10/(sigx3 * lwu), sweep(Zvar, MARGIN = 1, STATS = sigx11 *
      prZ * ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v

## logit specification class membership
cgradmnsftslnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  epsiuv <- (ewv1_h/ewu_h + S * epsilon/ewv1_h)
  depsi <- dnorm(-epsiuv, 0, 1)
  pepsi <- pnorm(-epsiuv)
  euv <- exp(ewv1/(2 * ewu) + S * epsilon/ewu_h)
  luvepsi <- ((1 + lambda) * ewv1_h/ewu_h + S * epsilon/ewv1_h)
  dluvepsi <- dnorm(-luvepsi, 0, 1)
  pluvepsi <- pnorm(-luvepsi)
  dv2 <- dnorm(S * epsilon/ewv2_h, 0, 1)
  sigx1 <- (depsi/ewv1_h - pepsi/ewu_h)
  sigx2 <- (dluvepsi/ewv1_h - (1 + lambda) * pluvepsi/ewu_h)
  sigx3 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * epsilon/ewu_h) * (1 + lambda))
  sigx4 <- ((1 + 2 * (lambda * ewu_h)) * wzdeno)
  sigx5 <- (2 * (sigx1 * euv) - sigx2 * sigx3)
  sigx6 <- ((1 + lambda) * sigx5 * ewz/sigx4 + S * prC * dv2 * epsilon/ewv2_h^3)
  sigx7 <- (2 * (euv * pepsi) - sigx3 * pluvepsi)
  sigx8 <- (prC * dv2/ewv2_h + (1 + lambda) * sigx7 * ewz/sigx4)
  sigx9 <- (0.5 * (S * epsilon/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx10 <- (0.5 * (depsi * ewv1_h/ewu_h) - sigx9 * pepsi)
  sigx11 <- (0.5 * (S * epsilon/ewu_h) + 2 * ((1 + lambda) * ewu * ewv1/(2 * ewu)^2))
  sigx12 <- (0.5 * (dluvepsi * ewv1_h/ewu_h) - sigx11 * pluvepsi)
  sigx13 <- (2 * (sigx10 * euv) - sigx12 * (1 + lambda) * sigx3)
  sigx14 <- (sigx13/sigx4 - lambda * wzdeno * sigx7 * ewu_h/sigx4^2)
  sigx15 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h))
  sigx16 <- (euv * (ewv1 * pepsi/(2 * ewu) - sigx15 * depsi))
  sigx17 <- (0.5 * ((1 + lambda) * ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h))
  sigx18 <- ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx17 * dluvepsi)
  sigx19 <- (sigx8 * (1 + 2 * (lambda * ewu_h)) * wzdeno)
  sigx20 <- (S^2 * dv2 * epsilon^2/ewv2_h^2)
  sigx21 <- ((1 + lambda) * ewv1/ewu + S * epsilon/ewu_h)
  sigx22 <- (sigx21 * pluvepsi - dluvepsi * ewv1_h/ewu_h)
  sigx23 <- (2 * (euv * pepsi) - (sigx22 * (1 + lambda) + pluvepsi) * sigx3)
  sigx24 <- (wzdeno * (1 + lambda) * sigx7 * ewu_h/sigx4^2)
  sigx25 <- ((1 + lambda) * (1/sigx4 - (1 + 2 * (lambda * ewu_h)) * ewz/sigx4^2) *
    sigx7 - prC * dv2/(wzdeno * ewv2_h))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx6/sigx8, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx14 * (1 + lambda) * ewz/sigx8, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (1 + lambda) * (2 * sigx16 - sigx18 * sigx3) *
      ewz/sigx19, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (0.5 * sigx20 -
      0.5 * dv2) * prC/(sigx8 * ewv2_h), FUN = "*"), (sigx23/sigx4 - 2 * sigx24) *
      ewz/sigx8, sweep(Zvar, MARGIN = 1, STATS = sigx25 * ewz/sigx8, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmnsftslnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  epsiuv <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  depsi <- dnorm(-epsiuv, 0, 1)
  pepsi <- pnorm(-epsiuv)
  euv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  luvepsi <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  dluvepsi <- dnorm(-luvepsi, 0, 1)
  pluvepsi <- pnorm(-luvepsi)
  dv2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  lwu <- (1 + 2 * (lambda * ewusr))
  sigx1 <- (depsi/ewvsr1 - pepsi/ewusr)
  sigx2 <- (dluvepsi/ewvsr1 - (1 + lambda) * pluvepsi/ewusr)
  sigx3 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  sigx4 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx5 <- (2 * (sigx1 * euv) - sigx2 * sigx3)
  sigx6 <- (ewz1 * (1 + lambda) * sigx5/lwu + S * ewz2 * dv2 * (epsilon)/ewvsr2^3)
  sigx7 <- (2 * (euv * pepsi) - sigx3 * pluvepsi)
  sigx8 <- (ewz2 * dv2/ewvsr2 + ewz1 * (1 + lambda) * sigx7/lwu)
  sigx9 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx10 <- (0.5 * (depsi * ewvsr1/ewusr) - sigx9 * pepsi)
  sigx11 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv1/(2 *
    ewu)^2))
  sigx12 <- (0.5 * (dluvepsi * ewvsr1/ewusr) - sigx11 * pluvepsi)
  sigx13 <- (2 * (sigx10 * euv) - (sigx12 * (1 + lambda) * sigx3 + lambda * sigx7 *
    ewusr/lwu))
  sigx14 <- (euv * (ewv1 * pepsi/(2 * ewu) - sigx4 * depsi))
  sigx15 <- ((1 + lambda) * ewvsr1/ewusr)
  sigx16 <- (0.5 * sigx15 - 0.5 * (S * (epsilon)/ewvsr1))
  sigx17 <- (2 * sigx14 - ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 *
    dluvepsi) * sigx3)
  sigx18 <- (0.5 * (S^2 * dv2 * (epsilon)^2/ewvsr2^2) - 0.5 * dv2)
  sigx19 <- ((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr)
  sigx20 <- (sigx19 * pluvepsi - dluvepsi * ewvsr1/ewusr)
  sigx21 <- (sigx20 * (1 + lambda) + pluvepsi)
  sigx22 <- (sigx21 * sigx3 + 2 * ((1 + lambda) * sigx7 * ewusr/lwu))
  sigx23 <- (2 * (euv * pepsi) - sigx22)
  sigx24 <- ((1 + lambda) * sigx7/lwu - dv2/ewvsr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx6/sigx8, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = ewz1 * (1 + lambda) * sigx13/(sigx8 * lwu),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ewz1 * (1 + lambda) * sigx17/(sigx8 *
      lwu), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ewz2 * sigx18/(sigx8 *
      ewvsr2), FUN = "*"), ewz1 * sigx23/(sigx8 * lwu), sweep(Zvar, MARGIN = 1,
      STATS = sigx24/(pi * sigx8 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmnsftslnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  epsiuv <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  depsi <- dnorm(-epsiuv, 0, 1)
  pepsi <- pnorm(-epsiuv)
  euv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  luvepsi <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  dluvepsi <- dnorm(-luvepsi, 0, 1)
  pluvepsi <- pnorm(-luvepsi)
  dv2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  lwu <- (1 + 2 * (lambda * ewusr))
  sigx1 <- (depsi/ewvsr1 - pepsi/ewusr)
  sigx2 <- (dluvepsi/ewvsr1 - (1 + lambda) * pluvepsi/ewusr)
  sigx3 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  sigx4 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx5 <- (2 * (sigx1 * euv) - sigx2 * sigx3)
  sigx6 <- ((1 + lambda) * sigx5 * pwZ/lwu + S * (1 - pwZ) * dv2 * (epsilon)/ewvsr2^3)
  sigx7 <- (2 * (euv * pepsi) - sigx3 * pluvepsi)
  sigx8 <- ((1 - pwZ) * dv2/ewvsr2 + (1 + lambda) * sigx7 * pwZ/lwu)
  sigx9 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx10 <- (0.5 * (depsi * ewvsr1/ewusr) - sigx9 * pepsi)
  sigx11 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv1/(2 *
    ewu)^2))
  sigx12 <- (0.5 * (dluvepsi * ewvsr1/ewusr) - sigx11 * pluvepsi)
  sigx13 <- (2 * (sigx10 * euv) - (sigx12 * (1 + lambda) * sigx3 + lambda * sigx7 *
    ewusr/lwu))
  sigx14 <- (euv * (ewv1 * pepsi/(2 * ewu) - sigx4 * depsi))
  sigx15 <- ((1 + lambda) * ewvsr1/ewusr)
  sigx16 <- (0.5 * sigx15 - 0.5 * (S * (epsilon)/ewvsr1))
  sigx17 <- (2 * sigx14 - ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 *
    dluvepsi) * sigx3)
  sigx18 <- (0.5 * (S^2 * dv2 * (epsilon)^2/ewvsr2^2) - 0.5 * dv2)
  sigx19 <- ((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr)
  sigx20 <- (sigx19 * pluvepsi - dluvepsi * ewvsr1/ewusr)
  sigx21 <- (sigx20 * (1 + lambda) + pluvepsi)
  sigx22 <- (sigx21 * sigx3 + 2 * ((1 + lambda) * sigx7 * ewusr/lwu))
  sigx23 <- (2 * (euv * pepsi) - sigx22)
  sigx24 <- ((1 + lambda) * sigx7/lwu - dv2/ewvsr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx6/sigx8, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (1 + lambda) * sigx13 * pwZ/(sigx8 * lwu),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (1 + lambda) * sigx17 *
      pwZ/(sigx8 * lwu), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx18 *
      (1 - pwZ)/(sigx8 * ewvsr2), FUN = "*"), sigx23 * pwZ/(sigx8 * lwu), sweep(Zvar,
      MARGIN = 1, STATS = sigx24 * dwZ/sigx8, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmnsftslnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  epsiuv <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  depsi <- dnorm(-epsiuv, 0, 1)
  pepsi <- pnorm(-epsiuv)
  euv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  luvepsi <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  dluvepsi <- dnorm(-luvepsi, 0, 1)
  pluvepsi <- pnorm(-luvepsi)
  dv2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  lwu <- (1 + 2 * (lambda * ewusr))
  sigx1 <- (depsi/ewvsr1 - pepsi/ewusr)
  sigx2 <- (dluvepsi/ewvsr1 - (1 + lambda) * pluvepsi/ewusr)
  sigx3 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  sigx4 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx5 <- (2 * (sigx1 * euv) - sigx2 * sigx3)
  sigx6 <- ((1 - prZ) * (1 + lambda) * sigx5/lwu + S * dv2 * prZ * (epsilon)/ewvsr2^3)
  sigx7 <- (2 * (euv * pepsi) - sigx3 * pluvepsi)
  sigx8 <- ((1 - prZ) * (1 + lambda) * sigx7/lwu + dv2 * prZ/ewvsr2)
  sigx9 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx10 <- (0.5 * (depsi * ewvsr1/ewusr) - sigx9 * pepsi)
  sigx11 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv1/(2 *
    ewu)^2))
  sigx12 <- (0.5 * (dluvepsi * ewvsr1/ewusr) - sigx11 * pluvepsi)
  sigx13 <- (2 * (sigx10 * euv) - (sigx12 * (1 + lambda) * sigx3 + lambda * sigx7 *
    ewusr/lwu))
  sigx14 <- (euv * (ewv1 * pepsi/(2 * ewu) - sigx4 * depsi))
  sigx15 <- ((1 + lambda) * ewvsr1/ewusr)
  sigx16 <- (0.5 * sigx15 - 0.5 * (S * (epsilon)/ewvsr1))
  sigx17 <- (2 * sigx14 - ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 *
    dluvepsi) * sigx3)
  sigx18 <- (0.5 * (S^2 * dv2 * (epsilon)^2/ewvsr2^2) - 0.5 * dv2)
  sigx19 <- ((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr)
  sigx20 <- (sigx19 * pluvepsi - dluvepsi * ewvsr1/ewusr)
  sigx21 <- (sigx20 * (1 + lambda) + pluvepsi)
  sigx22 <- (sigx21 * sigx3 + 2 * ((1 + lambda) * sigx7 * ewusr/lwu))
  sigx23 <- (2 * (euv * pepsi) - sigx22)
  sigx24 <- ((1 + lambda) * sigx7/lwu - dv2/ewvsr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx6/sigx8, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (1 - prZ) * (1 + lambda) * sigx13/(sigx8 *
      lwu), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (1 - prZ) * (1 + lambda) *
      sigx17/(sigx8 * lwu), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx18 *
      prZ/(sigx8 * ewvsr2), FUN = "*"), (1 - prZ) * sigx23/(sigx8 * lwu), sweep(Zvar,
      MARGIN = 1, STATS = sigx24 * prZ * ewz/sigx8, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf truncated_skewed_laplace-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
# Same sigma_v

## logit specification class membership
chesszisftslnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  wvwu <- ewvsr/ewusr
  lwv <- ((1 + lambda) * wvwu + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (wvwu + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wzl <- (lwu * wzdeno)
  lwzu <- (1/wzl - lwu * ewz/wzl^2)
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (wvwu) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * wvwu) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * wvwu) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - (1 + lambda) * pd/ewusr) *
    eB)
  sigx2 <- ((1 + lambda) * sigx1 * ewz/wzl + S * prC * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (prC * dwsr/ewvsr + (1 + lambda) * eaeb * ewz/wzl)
  sigx4 <- ((0.5 * (da * wvwu) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - (0.5 * (dd * wvwu) - sigx5 * pd) * (1 + lambda) * eB)
  sigx7 <- (sigx6/wzl - lambda * wzdeno * eaeb * ewusr/wzl^2)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (0.5 * sigx8 - 0.5 * dwsr) * prC/ewvsr
  sigx10 <- (eA * pawv)
  sigx11 <- ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd)
  sigx12 <- (sigx9 + (1 + lambda) * (2 * sigx10 - sigx11 * eB) * ewz/wzl)
  sigx13 <- (2 * (eA * pa) - (llwusq + pd) * eB)
  sigx14 <- (sigx13/wzl - 2 * (wzdeno * (1 + lambda) * eaeb * ewusr/wzl^2))
  sigx15 <- ((1 + lambda) * lwzu * eaeb - prC * dwsr/(wzdeno * ewvsr))
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (prC * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 + (1 + lambda) *
    (2 * (((epsiwv/ewvsr - 1/ewusr) * da/ewvsr - (da/ewvsr - pa/ewusr)/ewusr) *
      eA) - ((lwv/ewvsr - (1 + lambda)/ewusr) * dd/ewvsr - (1 + lambda) * (dd/ewvsr -
      (1 + lambda) * pd/ewusr)/ewusr) * eB) * ewz/wzl - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv) *
      pa - 0.5 * (da * wvwu))/ewusr + (0.5 * (epsiwv/ewusr) - epsiwu/ewvsr) *
      da) * eA) - ((0.5 * (lwv/ewusr) - sigx5/ewvsr) * dd + (0.5 * pd - (0.5 *
      (dd * wvwu) - sigx5 * pd) * (1 + lambda))/ewusr) * (1 + lambda) * eB)/wzl -
      (sigx2 * sigx7/sigx3 + lambda * wzdeno * sigx1 * ewusr/wzl^2)) * (1 +
      lambda) * ewz/sigx3, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((1 + lambda) * (2 * ((da * (ewv/(2 * ewu) -
      (epsiwr * epsiwv + 0.5))/ewvsr - pawv/ewusr) * eA) - (((1 + lambda)^2 *
      ewv/(2 * ewu) - (lwv * llwvsq + 0.5)) * dd/ewvsr - sigx11 * (1 + lambda)/ewusr) *
      eB) * ewz/wzl + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * prC *
      dwsr * (epsilon)/ewvsr^3 - sigx12 * sigx2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- matrix(colSums(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((da/ewvsr - pa/ewusr) * eA) - (((lwusq/ewvsr -
      lwv/ewusr) * dd - (llwusq + 2 * pd)/ewusr) * (1 + lambda) + dd/ewvsr) *
      eB)/wzl - (sigx2 * sigx14/sigx3 + 2 * (wzdeno * (1 + lambda) * sigx1 *
      ewusr/wzl^2))) * ewz/sigx3, FUN = "*")), ncol = 1)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((1 +
    lambda) * lwzu * sigx1 - (sigx15 * sigx2/sigx3 + S * prC * dwsr * (epsilon)/(wzdeno *
    ewvsr^3))) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr * epsiwv/ewusr) -
      0.5) - 0.5 * epsiwu) * da * wvwu - ((0.5 * (da * wvwu) - epsiwu * pa) *
      epsiwu + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv/(2 * ewu)^2) -
      0.25 * (S * (epsilon)/ewusr)) * pa)) * eA) - ((0.5 * (0.5 * (lwv * (1 +
      lambda) * wvwu) - 0.5) - 0.5 * (sigx5 * (1 + lambda))) * dd * wvwu -
      ((0.5 * (dd * wvwu) - sigx5 * pd) * sigx5 * (1 + lambda) + (2 * ((1 -
        8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) * ewu * ewv/(2 * ewu)^2) -
        0.25 * (S * (epsilon)/ewusr)) * pd)) * (1 + lambda) * eB)/wzl - (sigx7^2 *
      (1 + lambda) * ewz/sigx3 + lambda * ((0.5 - 2 * (lambda * lwu * wzdeno^2 *
      ewusr/wzl^2)) * eaeb + 4 * sigx4 - 2 * ((0.5 * (dd * wvwu) - sigx5 *
      pd) * (1 + lambda) * eB)) * wzdeno * ewusr/wzl^2)) * (1 + lambda) * ewz/sigx3,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((2 * (((da *
    ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pa/(2 * ewu)^2)) * ewv - ((0.5 * (epsiwr *
    epsiwv) - 0.25) * da * wvwu + epsiwu * pawv)) * eA) - (((1 + lambda) * dd *
    ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pd/(2 * ewu)^2)) * (1 + lambda) *
    ewv - (sigx11 * sigx5 + (0.5 * (lwv * llwvsq) - 0.25) * dd * wvwu)) * (1 +
    lambda) * eB)/wzl - (sigx12 * sigx7/sigx3 + lambda * wzdeno * (2 * sigx10 -
    sigx11 * eB) * ewusr/wzl^2)) * (1 + lambda) * ewz/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar + 1] <- matrix(colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * sigx4 - (((0.5 * lwusq - 0.5 * (lwv * wvwu)) *
      (1 + lambda) + 1) * dd * wvwu - ((llwusq + pd) * sigx5 + ((1 + lambda) *
      ewv/ewu + 0.5 * (S * (epsilon)/ewusr)) * pd)) * (1 + lambda) * eB)/wzl -
      (sigx7 * sigx14 * (1 + lambda) * ewz/sigx3 + wzdeno * (2 * (((0.5 - 2 *
        (lambda * lwu * wzdeno^2 * ewusr/wzl^2)) * eaeb + 2 * sigx4 - (0.5 *
        (dd * wvwu) - sigx5 * pd) * (1 + lambda) * eB) * (1 + lambda)) +
        lambda * sigx13) * ewusr/wzl^2)) * ewz/sigx3, FUN = "*")), ncol = 1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (lwzu * sigx6 - (sigx15 * sigx7 * ewz/sigx3 + lambda * ((2 - 2 * (lwu^2 *
      wzdeno^2/wzl^2)) * ewz + 1) * eaeb * ewusr/wzl^2)) * (1 + lambda) * ewz/sigx3,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr *
      (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * sigx8 - 0.5 * dwsr))/ewvsr + (1 +
      lambda) * (2 * (((pawv/2 + (pa - epsiwr * da)/2) * ewv/ewu - (0.25 *
      (wvwu) + 0.25 * (S * (epsilon)/ewvsr) - epsiwr^2 * epsiwv) * da) * eA) -
      ((sigx11/2 + (pd - llwvsq * dd)/2) * (1 + lambda)^2 * ewv/ewu - (0.25 *
        ((1 + lambda) * wvwu) + 0.25 * (S * (epsilon)/ewvsr) - lwv * llwvsq^2) *
        dd) * eB) * ewz/wzl - sigx12^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), nXvar + nuZUvar + nvZVvar +
    1] <- matrix(colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((2 * sigx10 -
    ((((llwusq + pd)/2 + pd) * (1 + lambda) * ewv/ewu - (lwusq * llwvsq + (0.5 -
      lwv * llwvsq) * wvwu) * dd) * (1 + lambda) - llwvsq * dd) * eB)/wzl -
    (sigx12 * sigx14/sigx3 + 2 * (wzdeno * (1 + lambda) * (2 * sigx10 - sigx11 *
      eB) * ewusr/wzl^2))) * ewz/sigx3, FUN = "*")), ncol = 1)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((1 + lambda) * lwzu * (2 * sigx10 - sigx11 *
      eB) - (sigx12 * sigx15/sigx3 + (0.5 * (S^2 * dwsr * (epsilon)^2/(wzdeno *
      ewvsr^3)) - 0.5 * (wzdeno * dwsr * ewvsr/(wzdeno * ewvsr)^2)) * prC)) *
      ewz/sigx3, FUN = "*"), Zvar)
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar + 1] <- sum(-((((((lwv *
    ewvsr - S * (epsilon))/ewusr - (1 + lambda) * ewv/ewu) * dd * wvwu + ewv *
    pd/ewu) * (1 + lambda) + (llwusq + 2 * pd) * lwusq - 2 * (dd * wvwu)) * eB/wzl +
    sigx14^2 * ewz/sigx3 + wzdeno * (2 * (2 * (eA * pa) - ((llwusq + pd) * eB +
    4 * (lwu * wzdeno^2 * (1 + lambda) * eaeb * ewusr/wzl^2))) + 2 * sigx13) *
    ewusr/wzl^2) * ewz/sigx3) * wHvar)
  hessll[nXvar + nuZUvar + nvZVvar + 1, (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- matrix(colSums(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * ((1/wzl - (((2 - 4 * (lwu^2 * wzdeno^2/wzl^2)) * ewz + 2 *
      wzdeno) * (1 + lambda) * ewusr + lwu * ewz)/wzl^2) * eaeb - (llwusq *
      lwzu * eB + sigx15 * sigx14 * ewz/sigx3)) * ewz/sigx3, FUN = "*")), ncol = 1)
  hessll[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 *
    ewvsr) + ewvsr/(wzdeno * ewvsr)^2) * dwsr - (sigx15^2/sigx3 + lwu * (1 +
    lambda) * (2 - 2 * (lwu^2 * wzdeno * ewz/wzl^2)) * eaeb/wzl^2)) * ewz + (1 +
    lambda) * lwzu * eaeb - prC * dwsr/(wzdeno * ewvsr)) * ewz/sigx3, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesszisftslnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  lwv <- ((1 + lambda) * ewvsr/ewusr + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * ewvsr/ewusr) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - (1 + lambda) * pd/ewusr) *
    eB)
  sigx2 <- (ewz1 * (1 + lambda) * sigx1/lwu + S * ewz2 * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (ewz2 * dwsr/ewvsr + ewz1 * (1 + lambda) * eaeb/lwu)
  sigx4 <- ((0.5 * (da * ewvsr/ewusr) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - ((0.5 * (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) *
    eB + lambda * eaeb * ewusr/lwu))
  sigx7 <- (2 * (eA * pawv) - ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd) *
    eB)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (ewz2 * (0.5 * sigx8 - 0.5 * dwsr)/ewvsr + ewz1 * (1 + lambda) * sigx7/lwu)
  sigx10 <- (2 * (eA * pa) - ((llwusq + pd) * eB + 2 * ((1 + lambda) * eaeb * ewusr/lwu)))
  sigx11 <- ((1 + lambda) * eaeb/lwu - dwsr/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (ewz2 * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 + ewz1 * (1 + lambda) *
      (2 * (((epsiwv/ewvsr - 1/ewusr) * da/ewvsr - (da/ewvsr - pa/ewusr)/ewusr) *
        eA) - ((lwv/ewvsr - (1 + lambda)/ewusr) * dd/ewvsr - (1 + lambda) *
        (dd/ewvsr - (1 + lambda) * pd/ewusr)/ewusr) * eB)/lwu - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv) *
      pa - 0.5 * (da * ewvsr/ewusr))/ewusr + (0.5 * (epsiwv/ewusr) - epsiwu/ewvsr) *
      da) * eA) - (((0.5 * (lwv/ewusr) - sigx5/ewvsr) * dd + (0.5 * pd - (0.5 *
      (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda))/ewusr) * (1 + lambda) *
      eB + lambda * sigx1 * ewusr/lwu))/(sigx3 * lwu) - sigx2 * lwu * sigx6/(sigx3 *
      lwu)^2) * ewz1 * (1 + lambda), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (ewz1 * (1 + lambda) * (2 * ((da * (ewv/(2 *
      ewu) - (epsiwr * epsiwv + 0.5))/ewvsr - pawv/ewusr) * eA) - (((1 + lambda)^2 *
      ewv/(2 * ewu) - (lwv * llwvsq + 0.5)) * dd/ewvsr - ((1 + lambda)^2 *
      ewv * pd/(2 * ewu) - llwvsq * dd) * (1 + lambda)/ewusr) * eB)/lwu + S *
      ewz2 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * dwsr * (epsilon)/ewvsr^3 -
      sigx9 * sigx2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((da/ewvsr - pa/ewusr) * eA) - ((((lwusq/ewvsr -
      lwv/ewusr) * dd - (llwusq + 2 * pd)/ewusr) * (1 + lambda) + dd/ewvsr) *
      eB + 2 * ((1 + lambda) * sigx1 * ewusr/lwu)))/(sigx3 * lwu) - sigx2 *
      lwu * sigx10/(sigx3 * lwu)^2) * ewz1, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (((1 +
    lambda) * sigx1/lwu - S * dwsr * (epsilon)/ewvsr^3)/(pi * sigx3 * ((Wz)^2 +
    1)) - pi * sigx2 * sigx11 * ((Wz)^2 + 1)/(pi * sigx3 * ((Wz)^2 + 1))^2),
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr * epsiwv/ewusr) -
      0.5) - 0.5 * epsiwu) * da * ewvsr/ewusr - ((0.5 * (da * ewvsr/ewusr) -
      epsiwu * pa) * epsiwu + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv/(2 *
      ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pa)) * eA) - (((0.5 * (0.5 *
      (lwv * (1 + lambda) * ewvsr/ewusr) - 0.5) - 0.5 * (sigx5 * (1 + lambda))) *
      dd * ewvsr/ewusr - ((0.5 * (dd * ewvsr/ewusr) - sigx5 * pd) * sigx5 *
      (1 + lambda) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) * ewu *
      ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pd)) * (1 + lambda) *
      eB + lambda * ((0.5 - lambda * ewusr/lwu) * eaeb + 2 * sigx4 - (0.5 *
      (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) * eB) * ewusr/lwu))/(sigx3 *
      lwu) - (ewz1 * (1 + lambda) * sigx6 + lambda * sigx3 * ewusr) * sigx6/(sigx3 *
      lwu)^2) * ewz1 * (1 + lambda), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ewz1 * (1 +
    lambda) * (2 * (((da * ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pa/(2 * ewu)^2)) *
    ewv - ((0.5 * (epsiwr * epsiwv) - 0.25) * da * ewvsr/ewusr + epsiwu * pawv)) *
    eA) - ((((1 + lambda) * dd * ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pd/(2 *
    ewu)^2)) * (1 + lambda) * ewv - (((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq *
    dd) * sigx5 + (0.5 * (lwv * llwvsq) - 0.25) * dd * ewvsr/ewusr)) * (1 + lambda) *
    eB + sigx9 * sigx6/sigx3 + lambda * sigx7 * ewusr/lwu))/(sigx3 * lwu), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * sigx4 - ((((0.5 * lwusq - 0.5 * (lwv *
      ewvsr/ewusr)) * (1 + lambda) + 1) * dd * ewvsr/ewusr - ((llwusq + pd) *
      sigx5 + ((1 + lambda) * ewv/ewu + 0.5 * (S * (epsilon)/ewusr)) * pd)) *
      eB + 2 * (((0.5 - lambda * ewusr/lwu) * eaeb + 2 * sigx4 - (0.5 * (dd *
      ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) * eB) * ewusr/lwu)) * (1 +
      lambda))/(sigx3 * lwu) - (ewz1 * (1 + lambda) * sigx6 + lambda * sigx3 *
      ewusr) * sigx10/(sigx3 * lwu)^2) * ewz1, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 + lambda) * (1/(pi * sigx3 * ((Wz)^2 + 1)) - pi * sigx11 * ((Wz)^2 + 1) *
    ewz1/(pi * sigx3 * ((Wz)^2 + 1))^2) * sigx6/lwu, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (ewz2 * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr *
      (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * sigx8 - 0.5 * dwsr))/ewvsr + ewz1 *
      (1 + lambda) * (2 * (((pawv/2 + (pa - epsiwr * da)/2) * ewv/ewu - (0.25 *
      (ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) - epsiwr^2 * epsiwv) * da) *
      eA) - ((((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd)/2 + (pd -
      llwvsq * dd)/2) * (1 + lambda)^2 * ewv/ewu - (0.25 * ((1 + lambda) *
      ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) - lwv * llwvsq^2) * dd) *
      eB)/lwu - sigx9^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((2 * (eA *
    pawv) - (((((llwusq + pd)/2 + pd) * (1 + lambda) * ewv/ewu - (lwusq * llwvsq +
    (0.5 - lwv * llwvsq) * ewvsr/ewusr) * dd) * (1 + lambda) - llwvsq * dd) *
    eB + 2 * ((1 + lambda) * sigx7 * ewusr/lwu)))/(sigx3 * lwu) - sigx9 * lwu *
    sigx10/(sigx3 * lwu)^2) * ewz1, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((1 + lambda) * sigx7/lwu - (0.5 * sigx8 - 0.5 *
      dwsr)/ewvsr)/(pi * sigx3 * ((Wz)^2 + 1)) - pi * sigx9 * sigx11 * ((Wz)^2 +
      1)/(pi * sigx3 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 1)] <- sum(wHvar *
    (-(((((((lwv * ewvsr - S * (epsilon))/ewusr - (1 + lambda) * ewv/ewu) * dd *
      ewvsr/ewusr + ewv * pd/ewu) * (1 + lambda) + (llwusq + 2 * pd) * lwusq -
      2 * (dd * ewvsr/ewusr)) * eB + 2 * (sigx10 * ewusr/lwu))/(sigx3 * lwu) +
      (ewz1 * sigx10 + 2 * (sigx3 * ewusr)) * sigx10/(sigx3 * lwu)^2) * ewz1)))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    (1/(pi * sigx3 * ((Wz)^2 + 1)) - pi * sigx11 * ((Wz)^2 + 1) * ewz1/(pi *
      sigx3 * ((Wz)^2 + 1))^2) * sigx10/lwu, FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar * (sigx11 * ((1 +
    lambda) * eaeb/lwu + 2 * (pi * Wz * sigx3) - dwsr/ewvsr)/(pi * sigx3 * ((Wz)^2 +
    1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesszisftslnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  lwv <- ((1 + lambda) * ewvsr/ewusr + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * ewvsr/ewusr) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - (1 + lambda) * pd/ewusr) *
    eB)
  sigx2 <- ((1 + lambda) * sigx1 * pwZ/lwu + S * (1 - pwZ) * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - pwZ) * dwsr/ewvsr + (1 + lambda) * eaeb * pwZ/lwu)
  sigx4 <- ((0.5 * (da * ewvsr/ewusr) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - ((0.5 * (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) *
    eB + lambda * eaeb * ewusr/lwu))
  sigx7 <- (2 * (eA * pawv) - ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd) *
    eB)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- ((0.5 * sigx8 - 0.5 * dwsr) * (1 - pwZ)/ewvsr + (1 + lambda) * sigx7 *
    pwZ/lwu)
  sigx10 <- (2 * (eA * pa) - ((llwusq + pd) * eB + 2 * ((1 + lambda) * eaeb * ewusr/lwu)))
  sigx11 <- ((1 + lambda) * eaeb/lwu - dwsr/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((1 - pwZ) * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 + (1 + lambda) *
      (2 * (((epsiwv/ewvsr - 1/ewusr) * da/ewvsr - (da/ewvsr - pa/ewusr)/ewusr) *
        eA) - ((lwv/ewvsr - (1 + lambda)/ewusr) * dd/ewvsr - (1 + lambda) *
        (dd/ewvsr - (1 + lambda) * pd/ewusr)/ewusr) * eB) * pwZ/lwu - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv) *
      pa - 0.5 * (da * ewvsr/ewusr))/ewusr + (0.5 * (epsiwv/ewusr) - epsiwu/ewvsr) *
      da) * eA) - (((0.5 * (lwv/ewusr) - sigx5/ewvsr) * dd + (0.5 * pd - (0.5 *
      (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda))/ewusr) * (1 + lambda) *
      eB + lambda * sigx1 * ewusr/lwu))/(sigx3 * lwu) - sigx2 * lwu * sigx6/(sigx3 *
      lwu)^2) * (1 + lambda) * pwZ, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((1 + lambda) * (2 * ((da * (ewv/(2 * ewu) -
      (epsiwr * epsiwv + 0.5))/ewvsr - pawv/ewusr) * eA) - (((1 + lambda)^2 *
      ewv/(2 * ewu) - (lwv * llwvsq + 0.5)) * dd/ewvsr - ((1 + lambda)^2 *
      ewv * pd/(2 * ewu) - llwvsq * dd) * (1 + lambda)/ewusr) * eB) * pwZ/lwu +
      S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * (1 - pwZ) * dwsr *
        (epsilon)/ewvsr^3 - sigx9 * sigx2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((da/ewvsr - pa/ewusr) * eA) - ((((lwusq/ewvsr -
      lwv/ewusr) * dd - (llwusq + 2 * pd)/ewusr) * (1 + lambda) + dd/ewvsr) *
      eB + 2 * ((1 + lambda) * sigx1 * ewusr/lwu)))/(sigx3 * lwu) - sigx2 *
      lwu * sigx10/(sigx3 * lwu)^2) * pwZ, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((1 +
    lambda) * sigx1/lwu - (sigx2 * sigx11/sigx3 + S * dwsr * (epsilon)/ewvsr^3)) *
    dwZ/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr * epsiwv/ewusr) -
      0.5) - 0.5 * epsiwu) * da * ewvsr/ewusr - ((0.5 * (da * ewvsr/ewusr) -
      epsiwu * pa) * epsiwu + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv/(2 *
      ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pa)) * eA) - (((0.5 * (0.5 *
      (lwv * (1 + lambda) * ewvsr/ewusr) - 0.5) - 0.5 * (sigx5 * (1 + lambda))) *
      dd * ewvsr/ewusr - ((0.5 * (dd * ewvsr/ewusr) - sigx5 * pd) * sigx5 *
      (1 + lambda) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) * ewu *
      ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pd)) * (1 + lambda) *
      eB + lambda * ((0.5 - lambda * ewusr/lwu) * eaeb + 2 * sigx4 - (0.5 *
      (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) * eB) * ewusr/lwu))/(sigx3 *
      lwu) - ((1 + lambda) * sigx6 * pwZ + lambda * sigx3 * ewusr) * sigx6/(sigx3 *
      lwu)^2) * (1 + lambda) * pwZ, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (1 + lambda) *
    (2 * (((da * ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pa/(2 * ewu)^2)) * ewv -
      ((0.5 * (epsiwr * epsiwv) - 0.25) * da * ewvsr/ewusr + epsiwu * pawv)) *
      eA) - ((((1 + lambda) * dd * ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pd/(2 *
      ewu)^2)) * (1 + lambda) * ewv - (((1 + lambda)^2 * ewv * pd/(2 * ewu) -
      llwvsq * dd) * sigx5 + (0.5 * (lwv * llwvsq) - 0.25) * dd * ewvsr/ewusr)) *
      (1 + lambda) * eB + sigx9 * sigx6/sigx3 + lambda * sigx7 * ewusr/lwu)) *
    pwZ/(sigx3 * lwu), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * sigx4 - ((((0.5 * lwusq - 0.5 * (lwv *
      ewvsr/ewusr)) * (1 + lambda) + 1) * dd * ewvsr/ewusr - ((llwusq + pd) *
      sigx5 + ((1 + lambda) * ewv/ewu + 0.5 * (S * (epsilon)/ewusr)) * pd)) *
      eB + 2 * (((0.5 - lambda * ewusr/lwu) * eaeb + 2 * sigx4 - (0.5 * (dd *
      ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) * eB) * ewusr/lwu)) * (1 +
      lambda))/(sigx3 * lwu) - ((1 + lambda) * sigx6 * pwZ + lambda * sigx3 *
      ewusr) * sigx10/(sigx3 * lwu)^2) * pwZ, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 - sigx11 * pwZ/sigx3) * (1 + lambda) * sigx6 * dwZ/(sigx3 * lwu), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((1 - pwZ) * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) *
      dwsr * (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * sigx8 - 0.5 * dwsr))/ewvsr +
      (1 + lambda) * (2 * (((pawv/2 + (pa - epsiwr * da)/2) * ewv/ewu - (0.25 *
        (ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) - epsiwr^2 * epsiwv) *
        da) * eA) - ((((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd)/2 +
        (pd - llwvsq * dd)/2) * (1 + lambda)^2 * ewv/ewu - (0.25 * ((1 +
        lambda) * ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) - lwv * llwvsq^2) *
        dd) * eB) * pwZ/lwu - sigx9^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((2 * (eA *
    pawv) - (((((llwusq + pd)/2 + pd) * (1 + lambda) * ewv/ewu - (lwusq * llwvsq +
    (0.5 - lwv * llwvsq) * ewvsr/ewusr) * dd) * (1 + lambda) - llwvsq * dd) *
    eB + 2 * ((1 + lambda) * sigx7 * ewusr/lwu)))/(sigx3 * lwu) - sigx9 * lwu *
    sigx10/(sigx3 * lwu)^2) * pwZ, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((1 + lambda) * sigx7/lwu - (sigx9 * sigx11/sigx3 +
      (0.5 * sigx8 - 0.5 * dwsr)/ewvsr)) * dwZ/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 1)] <- sum(wHvar *
    (-(((((((lwv * ewvsr - S * (epsilon))/ewusr - (1 + lambda) * ewv/ewu) * dd *
      ewvsr/ewusr + ewv * pd/ewu) * (1 + lambda) + (llwusq + 2 * pd) * lwusq -
      2 * (dd * ewvsr/ewusr)) * eB + 2 * (sigx10 * ewusr/lwu))/(sigx3 * lwu) +
      (sigx10 * pwZ + 2 * (sigx3 * ewusr)) * sigx10/(sigx3 * lwu)^2) * pwZ)))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    (1 - sigx11 * pwZ/sigx3) * sigx10 * dwZ/(sigx3 * lwu), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar * ((sigx11 * dwZ/sigx3 +
    Wz) * sigx11 * dwZ/sigx3), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesszisftslnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  lwv <- ((1 + lambda) * ewvsr/ewusr + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * ewvsr/ewusr) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - (1 + lambda) * pd/ewusr) *
    eB)
  sigx2 <- ((1 - prZ) * (1 + lambda) * sigx1/lwu + S * dwsr * prZ * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - prZ) * (1 + lambda) * eaeb/lwu + dwsr * prZ/ewvsr)
  sigx4 <- ((0.5 * (da * ewvsr/ewusr) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - ((0.5 * (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) *
    eB + lambda * eaeb * ewusr/lwu))
  sigx7 <- (2 * (eA * pawv) - ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd) *
    eB)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- ((0.5 * sigx8 - 0.5 * dwsr) * prZ/ewvsr + (1 - prZ) * (1 + lambda) *
    sigx7/lwu)
  sigx10 <- (2 * (eA * pa) - ((llwusq + pd) * eB + 2 * ((1 + lambda) * eaeb * ewusr/lwu)))
  sigx11 <- ((1 + lambda) * eaeb/lwu - dwsr/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((1 - prZ) * (1 + lambda) * (2 * (((epsiwv/ewvsr - 1/ewusr) * da/ewvsr -
      (da/ewvsr - pa/ewusr)/ewusr) * eA) - ((lwv/ewvsr - (1 + lambda)/ewusr) *
      dd/ewvsr - (1 + lambda) * (dd/ewvsr - (1 + lambda) * pd/ewusr)/ewusr) *
      eB)/lwu + dwsr * prZ * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv) *
      pa - 0.5 * (da * ewvsr/ewusr))/ewusr + (0.5 * (epsiwv/ewusr) - epsiwu/ewvsr) *
      da) * eA) - (((0.5 * (lwv/ewusr) - sigx5/ewvsr) * dd + (0.5 * pd - (0.5 *
      (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda))/ewusr) * (1 + lambda) *
      eB + lambda * sigx1 * ewusr/lwu))/(sigx3 * lwu) - sigx2 * lwu * sigx6/(sigx3 *
      lwu)^2) * (1 - prZ) * (1 + lambda), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((1 - prZ) * (1 + lambda) * (2 * ((da * (ewv/(2 *
      ewu) - (epsiwr * epsiwv + 0.5))/ewvsr - pawv/ewusr) * eA) - (((1 + lambda)^2 *
      ewv/(2 * ewu) - (lwv * llwvsq + 0.5)) * dd/ewvsr - ((1 + lambda)^2 *
      ewv * pd/(2 * ewu) - llwvsq * dd) * (1 + lambda)/ewusr) * eB)/lwu + S *
      (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * dwsr * prZ * (epsilon)/ewvsr^3 -
      sigx9 * sigx2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((da/ewvsr - pa/ewusr) * eA) - ((((lwusq/ewvsr -
      lwv/ewusr) * dd - (llwusq + 2 * pd)/ewusr) * (1 + lambda) + dd/ewvsr) *
      eB + 2 * ((1 + lambda) * sigx1 * ewusr/lwu)))/(sigx3 * lwu) - sigx2 *
      lwu * sigx10/(sigx3 * lwu)^2) * (1 - prZ), FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((1 +
    lambda) * sigx1/lwu - (sigx2 * sigx11/sigx3 + S * dwsr * (epsilon)/ewvsr^3)) *
    prZ * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr * epsiwv/ewusr) -
      0.5) - 0.5 * epsiwu) * da * ewvsr/ewusr - ((0.5 * (da * ewvsr/ewusr) -
      epsiwu * pa) * epsiwu + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv/(2 *
      ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pa)) * eA) - (((0.5 * (0.5 *
      (lwv * (1 + lambda) * ewvsr/ewusr) - 0.5) - 0.5 * (sigx5 * (1 + lambda))) *
      dd * ewvsr/ewusr - ((0.5 * (dd * ewvsr/ewusr) - sigx5 * pd) * sigx5 *
      (1 + lambda) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) * ewu *
      ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pd)) * (1 + lambda) *
      eB + lambda * ((0.5 - lambda * ewusr/lwu) * eaeb + 2 * sigx4 - (0.5 *
      (dd * ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) * eB) * ewusr/lwu))/(sigx3 *
      lwu) - ((1 - prZ) * (1 + lambda) * sigx6 + lambda * sigx3 * ewusr) *
      sigx6/(sigx3 * lwu)^2) * (1 - prZ) * (1 + lambda), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (1 - prZ) *
    (1 + lambda) * (2 * (((da * ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pa/(2 *
    ewu)^2)) * ewv - ((0.5 * (epsiwr * epsiwv) - 0.25) * da * ewvsr/ewusr + epsiwu *
    pawv)) * eA) - ((((1 + lambda) * dd * ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu *
    pd/(2 * ewu)^2)) * (1 + lambda) * ewv - (((1 + lambda)^2 * ewv * pd/(2 *
    ewu) - llwvsq * dd) * sigx5 + (0.5 * (lwv * llwvsq) - 0.25) * dd * ewvsr/ewusr)) *
    (1 + lambda) * eB + sigx9 * sigx6/sigx3 + lambda * sigx7 * ewusr/lwu))/(sigx3 *
    lwu), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * sigx4 - ((((0.5 * lwusq - 0.5 * (lwv *
      ewvsr/ewusr)) * (1 + lambda) + 1) * dd * ewvsr/ewusr - ((llwusq + pd) *
      sigx5 + ((1 + lambda) * ewv/ewu + 0.5 * (S * (epsilon)/ewusr)) * pd)) *
      eB + 2 * (((0.5 - lambda * ewusr/lwu) * eaeb + 2 * sigx4 - (0.5 * (dd *
      ewvsr/ewusr) - sigx5 * pd) * (1 + lambda) * eB) * ewusr/lwu)) * (1 +
      lambda))/(sigx3 * lwu) - ((1 - prZ) * (1 + lambda) * sigx6 + lambda *
      sigx3 * ewusr) * sigx10/(sigx3 * lwu)^2) * (1 - prZ), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 - sigx11 * (1 - prZ)/sigx3) * (1 + lambda) * sigx6 * prZ * ewz/(sigx3 *
    lwu), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((1 - prZ) * (1 + lambda) * (2 * (((pawv/2 + (pa - epsiwr * da)/2) * ewv/ewu -
      (0.25 * (ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) - epsiwr^2 * epsiwv) *
        da) * eA) - ((((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * dd)/2 +
      (pd - llwvsq * dd)/2) * (1 + lambda)^2 * ewv/ewu - (0.25 * ((1 + lambda) *
      ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) - lwv * llwvsq^2) * dd) *
      eB)/lwu + prZ * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) -
      0.25) * dwsr * (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * sigx8 - 0.5 * dwsr))/ewvsr -
      sigx9^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((2 * (eA *
    pawv) - (((((llwusq + pd)/2 + pd) * (1 + lambda) * ewv/ewu - (lwusq * llwvsq +
    (0.5 - lwv * llwvsq) * ewvsr/ewusr) * dd) * (1 + lambda) - llwvsq * dd) *
    eB + 2 * ((1 + lambda) * sigx7 * ewusr/lwu)))/(sigx3 * lwu) - sigx9 * lwu *
    sigx10/(sigx3 * lwu)^2) * (1 - prZ), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((1 + lambda) * sigx7/lwu - (sigx9 * sigx11/sigx3 +
      (0.5 * sigx8 - 0.5 * dwsr)/ewvsr)) * prZ * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 1)] <- sum(wHvar *
    (-(((((((lwv * ewvsr - S * (epsilon))/ewusr - (1 + lambda) * ewv/ewu) * dd *
      ewvsr/ewusr + ewv * pd/ewu) * (1 + lambda) + (llwusq + 2 * pd) * lwusq -
      2 * (dd * ewvsr/ewusr)) * eB + 2 * (sigx10 * ewusr/lwu))/(sigx3 * lwu) +
      ((1 - prZ) * sigx10 + 2 * (sigx3 * ewusr)) * sigx10/(sigx3 * lwu)^2) *
      (1 - prZ))))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    (1 - sigx11 * (1 - prZ)/sigx3) * sigx10 * prZ * ewz/(sigx3 * lwu), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar +
    1)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar * sigx11 * (1 - (sigx11 *
    prZ/sigx3 + 1) * ewz) * prZ * ewz/sigx3, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Different sigma_v

## logit specification class membership
chessmnsftslnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  epsiuv <- (ewv1_h/ewu_h + S * epsilon/ewv1_h)
  depsi <- dnorm(-epsiuv, 0, 1)
  pepsi <- pnorm(-epsiuv)
  euv <- exp(ewv1/(2 * ewu) + S * epsilon/ewu_h)
  luvepsi <- ((1 + lambda) * ewv1_h/ewu_h + S * epsilon/ewv1_h)
  dluvepsi <- dnorm(-luvepsi, 0, 1)
  pluvepsi <- pnorm(-luvepsi)
  dv2 <- dnorm(S * epsilon/ewv2_h, 0, 1)
  sigx1 <- (depsi/ewv1_h - pepsi/ewu_h)
  sigx2 <- (dluvepsi/ewv1_h - (1 + lambda) * pluvepsi/ewu_h)
  sigx3 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * epsilon/ewu_h) * (1 + lambda))
  sigx4 <- ((1 + 2 * (lambda * ewu_h)) * wzdeno)
  sigx5 <- (2 * (sigx1 * euv) - sigx2 * sigx3)
  sigx6 <- ((1 + lambda) * sigx5 * ewz/sigx4 + S * prC * dv2 * epsilon/ewv2_h^3)
  sigx7 <- (2 * (euv * pepsi) - sigx3 * pluvepsi)
  sigx8 <- (prC * dv2/ewv2_h + (1 + lambda) * sigx7 * ewz/sigx4)
  sigx9 <- (0.5 * (S * epsilon/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx10 <- (0.5 * (depsi * ewv1_h/ewu_h) - sigx9 * pepsi)
  sigx11 <- (0.5 * (S * epsilon/ewu_h) + 2 * ((1 + lambda) * ewu * ewv1/(2 * ewu)^2))
  sigx12 <- (0.5 * (dluvepsi * ewv1_h/ewu_h) - sigx11 * pluvepsi)
  sigx13 <- (2 * (sigx10 * euv) - sigx12 * (1 + lambda) * sigx3)
  sigx14 <- (sigx13/sigx4 - lambda * wzdeno * sigx7 * ewu_h/sigx4^2)
  sigx15 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h))
  sigx16 <- (euv * (ewv1 * pepsi/(2 * ewu) - sigx15 * depsi))
  sigx17 <- (0.5 * ((1 + lambda) * ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h))
  sigx18 <- ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx17 * dluvepsi)
  sigx19 <- (sigx8 * (1 + 2 * (lambda * ewu_h)) * wzdeno)
  sigx20 <- (S^2 * dv2 * epsilon^2/ewv2_h^2)
  sigx21 <- ((1 + lambda) * ewv1/ewu + S * epsilon/ewu_h)
  sigx22 <- (sigx21 * pluvepsi - dluvepsi * ewv1_h/ewu_h)
  sigx23 <- (2 * (euv * pepsi) - (sigx22 * (1 + lambda) + pluvepsi) * sigx3)
  sigx24 <- (wzdeno * (1 + lambda) * sigx7 * ewu_h/sigx4^2)
  sigx25 <- ((1 + lambda) * (1/sigx4 - (1 + 2 * (lambda * ewu_h)) * ewz/sigx4^2) *
    sigx7 - prC * dv2/(wzdeno * ewv2_h))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (prC * dv2 * (S^2 * epsilon^2/ewv2_h^2 - 1)/ewv2_h^3 + (1 + lambda) *
    (2 * (((epsiuv/ewv1_h - 1/ewu_h) * depsi/ewv1_h - sigx1/ewu_h) * euv) - ((luvepsi/ewv1_h -
      (1 + lambda)/ewu_h) * dluvepsi/ewv1_h - (1 + lambda) * sigx2/ewu_h) *
      sigx3) * ewz/sigx4 - sigx6^2/sigx8)/sigx8, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((((0.5 + 0.5 * (S * epsilon/ewu_h) + 2 * (ewu *
      ewv1/(2 * ewu)^2)) * pepsi - 0.5 * (depsi * ewv1_h/ewu_h))/ewu_h + (0.5 *
      (epsiuv/ewu_h) - sigx9/ewv1_h) * depsi) * euv) - ((0.5 * (luvepsi/ewu_h) -
      sigx11/ewv1_h) * dluvepsi + (0.5 * pluvepsi - sigx12 * (1 + lambda))/ewu_h) *
      (1 + lambda) * sigx3)/sigx4 - (sigx6 * sigx14/sigx8 + lambda * wzdeno *
      sigx5 * ewu_h/sigx4^2)) * (1 + lambda) * ewz/sigx8, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((depsi * (ewv1/(2 * ewu) - (sigx15 *
      epsiuv + 0.5))/ewv1_h - (ewv1 * pepsi/(2 * ewu) - sigx15 * depsi)/ewu_h) *
      euv) - (((1 + lambda)^2 * ewv1/(2 * ewu) - (luvepsi * sigx17 + 0.5)) *
      dluvepsi/ewv1_h - sigx18 * (1 + lambda)/ewu_h) * sigx3)/sigx19 - sigx6 *
      (1 + 2 * (lambda * ewu_h)) * wzdeno * (2 * sigx16 - sigx18 * sigx3)/sigx19^2) *
      (1 + lambda) * ewz, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * prC * (S * (0.5 * (S^2 * epsilon^2/ewv2_h^2 -
      2) - 0.5) * dv2 * epsilon/(sigx8 * ewv2_h^3) - sigx6 * (0.5 * sigx20 -
      0.5 * dv2) * ewv2_h/(sigx8 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * (sigx1 * euv) - (((sigx21/ewv1_h - luvepsi/ewu_h) *
      dluvepsi - (sigx22 * (1 + lambda) + 2 * pluvepsi)/ewu_h) * (1 + lambda) +
      dluvepsi/ewv1_h) * sigx3)/sigx4 - (sigx6 * (sigx23/sigx4 - 2 * sigx24)/sigx8 +
      2 * (wzdeno * (1 + lambda) * sigx5 * ewu_h/sigx4^2))) * ewz/sigx8, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((1 +
    lambda) * (1/sigx4 - (1 + 2 * (lambda * ewu_h)) * ewz/sigx4^2) * sigx5 -
    (sigx25 * sigx6/sigx8 + S * prC * dv2 * epsilon/(wzdeno * ewv2_h^3))) * ewz/sigx8,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewv1_h * epsiuv/ewu_h) -
      0.5) - 0.5 * sigx9) * depsi * ewv1_h/ewu_h - (sigx10 * sigx9 + (2 * ((1 -
      8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * epsilon/ewu_h)) *
      pepsi)) * euv) - ((0.5 * (0.5 * (luvepsi * (1 + lambda) * ewv1_h/ewu_h) -
      0.5) - 0.5 * (sigx11 * (1 + lambda))) * dluvepsi * ewv1_h/ewu_h - (sigx12 *
      sigx11 * (1 + lambda) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) *
      ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * epsilon/ewu_h)) * pluvepsi)) *
      (1 + lambda) * sigx3)/sigx4 - (sigx14^2 * (1 + lambda) * ewz/sigx8 +
      lambda * ((0.5 - 2 * (lambda * (1 + 2 * (lambda * ewu_h)) * wzdeno^2 *
        ewu_h/sigx4^2)) * sigx7 + 4 * (sigx10 * euv) - 2 * (sigx12 * (1 +
        lambda) * sigx3)) * wzdeno * ewu_h/sigx4^2)) * (1 + lambda) * ewz/sigx8,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((2 * (((depsi *
    ewv1_h/(4 * (ewu * ewu_h)) - 2 * (ewu * pepsi/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx15 * epsiuv) - 0.25) * depsi * ewv1_h/ewu_h + sigx9 * (ewv1 * pepsi/(2 *
    ewu) - sigx15 * depsi))) * euv) - (((1 + lambda) * dluvepsi * ewv1_h/(4 *
    (ewu * ewu_h)) - 2 * (ewu * pluvepsi/(2 * ewu)^2)) * (1 + lambda) * ewv1 -
    (sigx18 * sigx11 + (0.5 * (luvepsi * sigx17) - 0.25) * dluvepsi * ewv1_h/ewu_h)) *
    (1 + lambda) * sigx3)/sigx19 - (sigx14 * (1 + 2 * (lambda * ewu_h)) * (1 +
    lambda) * ewz + lambda * sigx8 * ewu_h) * wzdeno * (2 * sigx16 - sigx18 *
    sigx3)/sigx19^2) * (1 + lambda) * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx14 * (0.5 * sigx20 - 0.5 * dv2) * prC * (1 + lambda) * ewv2_h * ewz/(sigx8 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx10 * euv) - (((0.5 * sigx21 - 0.5 *
      (luvepsi * ewv1_h/ewu_h)) * (1 + lambda) + 1) * dluvepsi * ewv1_h/ewu_h -
      ((sigx22 * (1 + lambda) + pluvepsi) * sigx11 + ((1 + lambda) * ewv1/ewu +
        0.5 * (S * epsilon/ewu_h)) * pluvepsi)) * (1 + lambda) * sigx3)/sigx4 -
      (sigx14 * (sigx23/sigx4 - 2 * sigx24) * (1 + lambda) * ewz/sigx8 + wzdeno *
        (2 * (((0.5 - 2 * (lambda * (1 + 2 * (lambda * ewu_h)) * wzdeno^2 *
          ewu_h/sigx4^2)) * sigx7 + 2 * (sigx10 * euv) - sigx12 * (1 + lambda) *
          sigx3) * (1 + lambda)) + lambda * sigx23) * ewu_h/sigx4^2)) * ewz/sigx8,
    FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ((1/sigx4 - (1 + 2 * (lambda * ewu_h)) * ewz/sigx4^2) * sigx13 -
      (sigx25 * sigx14 * ewz/sigx8 + lambda * ((2 - 2 * ((1 + 2 * (lambda *
        ewu_h))^2 * wzdeno^2/sigx4^2)) * ewz + 1) * sigx7 * ewu_h/sigx4^2)) *
      (1 + lambda) * ewz/sigx8, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((2 * ((((ewv1 * pepsi/(2 * ewu) - sigx15 * depsi)/2 + (pepsi - sigx15 *
      depsi)/2) * ewv1/ewu - (0.25 * (ewv1_h/ewu_h) + 0.25 * (S * epsilon/ewv1_h) -
      sigx15^2 * epsiuv) * depsi) * euv) - ((sigx18/2 + (pluvepsi - sigx17 *
      dluvepsi)/2) * (1 + lambda)^2 * ewv1/ewu - (0.25 * ((1 + lambda) * ewv1_h/ewu_h) +
      0.25 * (S * epsilon/ewv1_h) - luvepsi * sigx17^2) * dluvepsi) * sigx3)/sigx19 -
      (1 + lambda) * (2 * sigx16 - sigx18 * sigx3)^2 * ewz/sigx19^2) * (1 +
    lambda) * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * ((0.5 * sigx20 - 0.5 * dv2) * prC * (1 + lambda) * (2 *
      sigx16 - sigx18 * sigx3) * ewv2_h * ewz/((sigx8 * ewv2_h)^2 * (1 + 2 *
      (lambda * ewu_h)) * wzdeno)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((2 *
    sigx16 - ((((sigx22 * (1 + lambda) + pluvepsi)/2 + pluvepsi) * (1 + lambda) *
    ewv1/ewu - (sigx21 * sigx17 + (0.5 - luvepsi * sigx17) * ewv1_h/ewu_h) *
    dluvepsi) * (1 + lambda) - sigx17 * dluvepsi) * sigx3)/sigx4 - ((sigx23/sigx4 -
    2 * sigx24) * ewz/sigx19 + 2 * (wzdeno * ewu_h/sigx4^2)) * (1 + lambda) *
    (2 * sigx16 - sigx18 * sigx3)) * ewz/sigx8, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 + lambda) * (1/sigx4 - (sigx25/sigx19 + (1 +
      2 * (lambda * ewu_h))/sigx4^2) * ewz) * (2 * sigx16 - sigx18 * sigx3) *
      ewz/sigx8, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (S^2 * (0.5 * (0.5 * (S^2 * epsilon^2/ewv2_h^2) -
      1) - 0.25) * dv2 * epsilon^2/(sigx8 * ewv2_h^3) - ((0.5 * sigx20 - 0.5 *
      dv2) * prC + 0.5 * (sigx8 * ewv2_h)) * (0.5 * sigx20 - 0.5 * dv2)/(sigx8 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx23/sigx4 - 2 * sigx24) * (0.5 * sigx20 - 0.5 * dv2) * prC * ewz/(sigx8^2 *
      ewv2_h)), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx25 * (0.5 * sigx20 - 0.5 * dv2)/sigx8 +
      0.5 * (S^2 * dv2 * epsilon^2/(wzdeno * ewv2_h^2)))/ewv2_h - 0.5 * (wzdeno *
      dv2 * ewv2_h/(wzdeno * ewv2_h)^2)) * prC * ewz/sigx8), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    1)] <- sum(wHvar * (-((((((luvepsi * ewv1_h - S * epsilon)/ewu_h - (1 + lambda) *
    ewv1/ewu) * dluvepsi * ewv1_h/ewu_h + ewv1 * pluvepsi/ewu) * (1 + lambda) +
    (sigx22 * (1 + lambda) + 2 * pluvepsi) * sigx21 - 2 * (dluvepsi * ewv1_h/ewu_h)) *
    sigx3/sigx4 + (sigx23/sigx4 - 2 * sigx24)^2 * ewz/sigx8 + wzdeno * (2 * (2 *
    (euv * pepsi) - ((sigx22 * (1 + lambda) + pluvepsi) * sigx3 + 4 * ((1 + 2 *
    (lambda * ewu_h)) * wzdeno^2 * (1 + lambda) * sigx7 * ewu_h/sigx4^2))) +
    2 * sigx23) * ewu_h/sigx4^2) * ewz/sigx8)))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * ((1/sigx4 - (((2 - 4 * ((1 + 2 * (lambda * ewu_h))^2 * wzdeno^2/sigx4^2)) *
      ewz + 2 * wzdeno) * (1 + lambda) * ewu_h + (1 + 2 * (lambda * ewu_h)) *
      ewz)/sigx4^2) * sigx7 - (sigx22 * (1 + lambda) * (1/sigx4 - (1 + 2 *
      (lambda * ewu_h)) * ewz/sigx4^2) * sigx3 + sigx25 * (sigx23/sigx4 - 2 *
      sigx24) * ewz/sigx8)) * ewz/sigx8, FUN = "*"))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar * ((prC *
    (1/(wzdeno^2 * ewv2_h) + ewv2_h/(wzdeno * ewv2_h)^2) * dv2 - (sigx25^2/sigx8 +
    (1 + 2 * (lambda * ewu_h)) * (1 + lambda) * (2 - 2 * ((1 + 2 * (lambda *
      ewu_h))^2 * wzdeno * ewz/sigx4^2)) * sigx7/sigx4^2)) * ewz + (1 + lambda) *
    (1/sigx4 - (1 + 2 * (lambda * ewu_h)) * ewz/sigx4^2) * sigx7 - prC * dv2/(wzdeno *
    ewv2_h)) * ewz/sigx8, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmnsftslnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  epsiuv <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  depsi <- dnorm(-epsiuv, 0, 1)
  pepsi <- pnorm(-epsiuv)
  euv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  luvepsi <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  dluvepsi <- dnorm(-luvepsi, 0, 1)
  pluvepsi <- pnorm(-luvepsi)
  dv2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  lwu <- (1 + 2 * (lambda * ewusr))
  sigx1 <- (depsi/ewvsr1 - pepsi/ewusr)
  sigx2 <- (dluvepsi/ewvsr1 - (1 + lambda) * pluvepsi/ewusr)
  sigx3 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  sigx4 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx5 <- (2 * (sigx1 * euv) - sigx2 * sigx3)
  sigx6 <- (ewz1 * (1 + lambda) * sigx5/lwu + S * ewz2 * dv2 * (epsilon)/ewvsr2^3)
  sigx7 <- (2 * (euv * pepsi) - sigx3 * pluvepsi)
  sigx8 <- (ewz2 * dv2/ewvsr2 + ewz1 * (1 + lambda) * sigx7/lwu)
  sigx9 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx10 <- (0.5 * (depsi * ewvsr1/ewusr) - sigx9 * pepsi)
  sigx11 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv1/(2 *
    ewu)^2))
  sigx12 <- (0.5 * (dluvepsi * ewvsr1/ewusr) - sigx11 * pluvepsi)
  sigx13 <- (2 * (sigx10 * euv) - (sigx12 * (1 + lambda) * sigx3 + lambda * sigx7 *
    ewusr/lwu))
  sigx14 <- (euv * (ewv1 * pepsi/(2 * ewu) - sigx4 * depsi))
  sigx15 <- ((1 + lambda) * ewvsr1/ewusr)
  sigx16 <- (0.5 * sigx15 - 0.5 * (S * (epsilon)/ewvsr1))
  sigx17 <- (2 * sigx14 - ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 *
    dluvepsi) * sigx3)
  sigx18 <- (0.5 * (S^2 * dv2 * (epsilon)^2/ewvsr2^2) - 0.5 * dv2)
  sigx19 <- ((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr)
  sigx20 <- (sigx19 * pluvepsi - dluvepsi * ewvsr1/ewusr)
  sigx21 <- (sigx20 * (1 + lambda) + pluvepsi)
  sigx22 <- (sigx21 * sigx3 + 2 * ((1 + lambda) * sigx7 * ewusr/lwu))
  sigx23 <- (2 * (euv * pepsi) - sigx22)
  sigx24 <- ((1 + lambda) * sigx7/lwu - dv2/ewvsr2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (ewz2 * dv2 * (S^2 * (epsilon)^2/ewvsr2^2 - 1)/ewvsr2^3 + ewz1 * (1 + lambda) *
      (2 * (((epsiuv/ewvsr1 - 1/ewusr) * depsi/ewvsr1 - sigx1/ewusr) * euv) -
        ((luvepsi/ewvsr1 - (1 + lambda)/ewusr) * dluvepsi/ewvsr1 - (1 + lambda) *
          sigx2/ewusr) * sigx3)/lwu - sigx6^2/sigx8)/sigx8, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv1/(2 * ewu)^2)) * pepsi - 0.5 * (depsi * ewvsr1/ewusr))/ewusr + (0.5 *
      (epsiuv/ewusr) - sigx9/ewvsr1) * depsi) * euv) - (((0.5 * (luvepsi/ewusr) -
      sigx11/ewvsr1) * dluvepsi + (0.5 * pluvepsi - sigx12 * (1 + lambda))/ewusr) *
      (1 + lambda) * sigx3 + lambda * sigx5 * ewusr/lwu))/(sigx8 * lwu) - sigx6 *
      lwu * sigx13/(sigx8 * lwu)^2) * ewz1 * (1 + lambda), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((depsi * (ewv1/(2 * ewu) - (sigx4 *
      epsiuv + 0.5))/ewvsr1 - (ewv1 * pepsi/(2 * ewu) - sigx4 * depsi)/ewusr) *
      euv) - (((1 + lambda)^2 * ewv1/(2 * ewu) - (luvepsi * sigx16 + 0.5)) *
      dluvepsi/ewvsr1 - ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 *
      dluvepsi) * (1 + lambda)/ewusr) * sigx3)/(sigx8 * lwu) - sigx6 * lwu *
      sigx17/(sigx8 * lwu)^2) * ewz1 * (1 + lambda), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ewz2 * (S * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2 -
      2) - 0.5) * dv2 * (epsilon)/(sigx8 * ewvsr2^3) - sigx6 * sigx18 * ewvsr2/(sigx8 *
      ewvsr2)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * (sigx1 * euv) - ((((sigx19/ewvsr1 - luvepsi/ewusr) *
      dluvepsi - (sigx20 * (1 + lambda) + 2 * pluvepsi)/ewusr) * (1 + lambda) +
      dluvepsi/ewvsr1) * sigx3 + 2 * ((1 + lambda) * sigx5 * ewusr/lwu)))/(sigx8 *
      lwu) - sigx6 * lwu * sigx23/(sigx8 * lwu)^2) * ewz1, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (((1 +
    lambda) * sigx5/lwu - S * dv2 * (epsilon)/ewvsr2^3)/(pi * sigx8 * ((Wz)^2 +
    1)) - pi * sigx6 * sigx24 * ((Wz)^2 + 1)/(pi * sigx8 * ((Wz)^2 + 1))^2),
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr1 * epsiuv/ewusr) -
      0.5) - 0.5 * sigx9) * depsi * ewvsr1/ewusr - (sigx10 * sigx9 + (2 * ((1 -
      8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) *
      pepsi)) * euv) - (((0.5 * (0.5 * (luvepsi * (1 + lambda) * ewvsr1/ewusr) -
      0.5) - 0.5 * (sigx11 * (1 + lambda))) * dluvepsi * ewvsr1/ewusr - (sigx12 *
      sigx11 * (1 + lambda) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) *
      ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pluvepsi)) *
      (1 + lambda) * sigx3 + lambda * ((0.5 - lambda * ewusr/lwu) * sigx7 +
      2 * (sigx10 * euv) - sigx12 * (1 + lambda) * sigx3) * ewusr/lwu))/(sigx8 *
      lwu) - (ewz1 * (1 + lambda) * sigx13 + lambda * sigx8 * ewusr) * sigx13/(sigx8 *
      lwu)^2) * ewz1 * (1 + lambda), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((2 * (((depsi *
    ewvsr1/(4 * (ewu * ewusr)) - 2 * (ewu * pepsi/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx4 * epsiuv) - 0.25) * depsi * ewvsr1/ewusr + sigx9 * (ewv1 * pepsi/(2 *
    ewu) - sigx4 * depsi))) * euv) - (((1 + lambda) * dluvepsi * ewvsr1/(4 *
    (ewu * ewusr)) - 2 * (ewu * pluvepsi/(2 * ewu)^2)) * (1 + lambda) * ewv1 -
    (((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 * dluvepsi) * sigx11 +
      (0.5 * (luvepsi * sigx16) - 0.25) * dluvepsi * ewvsr1/ewusr)) * (1 +
    lambda) * sigx3)/(sigx8 * lwu) - (ewz1 * (1 + lambda) * sigx13 + lambda *
    sigx8 * ewusr) * sigx17/(sigx8 * lwu)^2) * ewz1 * (1 + lambda), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (ewz2 * sigx18 * ewz1 * (1 + lambda) * sigx13 * ewvsr2/((sigx8 * ewvsr2)^2 *
      lwu)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx10 * euv) - ((((0.5 * sigx19 - 0.5 *
      (luvepsi * ewvsr1/ewusr)) * (1 + lambda) + 1) * dluvepsi * ewvsr1/ewusr -
      (sigx21 * sigx11 + ((1 + lambda) * ewv1/ewu + 0.5 * (S * (epsilon)/ewusr)) *
        pluvepsi)) * sigx3 + 2 * (((0.5 - lambda * ewusr/lwu) * sigx7 + 2 *
      (sigx10 * euv) - sigx12 * (1 + lambda) * sigx3) * ewusr/lwu)) * (1 +
      lambda))/(sigx8 * lwu) - (ewz1 * (1 + lambda) * sigx13 + lambda * sigx8 *
      ewusr) * sigx23/(sigx8 * lwu)^2) * ewz1, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * (1 + lambda) * (1/(pi * sigx8 * ((Wz)^2 + 1)) - pi * sigx24 *
      ((Wz)^2 + 1) * ewz1/(pi * sigx8 * ((Wz)^2 + 1))^2) * sigx13/lwu, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((2 * ((((ewv1 * pepsi/(2 * ewu) - sigx4 * depsi)/2 + (pepsi - sigx4 * depsi)/2) *
      ewv1/ewu - (0.25 * (ewvsr1/ewusr) + 0.25 * (S * (epsilon)/ewvsr1) - sigx4^2 *
      epsiuv) * depsi) * euv) - ((((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) -
      sigx16 * dluvepsi)/2 + (pluvepsi - sigx16 * dluvepsi)/2) * (1 + lambda)^2 *
      ewv1/ewu - (0.25 * sigx15 + 0.25 * (S * (epsilon)/ewvsr1) - luvepsi *
      sigx16^2) * dluvepsi) * sigx3)/(sigx8 * lwu) - ewz1 * (1 + lambda) *
      sigx17^2/(sigx8 * lwu)^2) * ewz1 * (1 + lambda), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (ewz2 * sigx18 * ewz1 * (1 + lambda) * sigx17 * ewvsr2/((sigx8 *
      ewvsr2)^2 * lwu)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((2 *
    sigx14 - ((((sigx21/2 + pluvepsi) * (1 + lambda) * ewv1/ewu - (sigx19 * sigx16 +
    (0.5 - luvepsi * sigx16) * ewvsr1/ewusr) * dluvepsi) * (1 + lambda) - sigx16 *
    dluvepsi) * sigx3 + 2 * ((1 + lambda) * sigx17 * ewusr/lwu)))/(sigx8 * lwu) -
    ewz1 * (1 + lambda) * sigx17 * sigx23/(sigx8 * lwu)^2) * ewz1, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 + lambda) * (1/(pi * sigx8 * ((Wz)^2 + 1)) -
      pi * sigx24 * ((Wz)^2 + 1) * ewz1/(pi * sigx8 * ((Wz)^2 + 1))^2) * sigx17/lwu,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ewz2 * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2) -
      1) - 0.25) * dv2 * (epsilon)^2/(sigx8 * ewvsr2^3) - (ewz2 * sigx18 +
      0.5 * (sigx8 * ewvsr2)) * sigx18/(sigx8 * ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (ewz2 * sigx18 * ewz1 * lwu * sigx23/((sigx8 * lwu)^2 * ewvsr2)), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx18 * (1/(pi * sigx8 * ((Wz)^2 + 1)) + pi *
      sigx24 * ((Wz)^2 + 1) * ewz2/(pi * sigx8 * ((Wz)^2 + 1))^2)/ewvsr2),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    1)] <- sum(wHvar * (-(((((((luvepsi * ewvsr1 - S * (epsilon))/ewusr - (1 +
    lambda) * ewv1/ewu) * dluvepsi * ewvsr1/ewusr + ewv1 * pluvepsi/ewu) * (1 +
    lambda) + (sigx20 * (1 + lambda) + 2 * pluvepsi) * sigx19 - 2 * (dluvepsi *
    ewvsr1/ewusr)) * sigx3 + 2 * (sigx23 * ewusr/lwu))/(sigx8 * lwu) + (ewz1 *
    sigx23 + 2 * (sigx8 * ewusr)) * sigx23/(sigx8 * lwu)^2) * ewz1)))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * (1/(pi * sigx8 * ((Wz)^2 + 1)) - pi * sigx24 * ((Wz)^2 +
      1) * ewz1/(pi * sigx8 * ((Wz)^2 + 1))^2) * sigx23/lwu, FUN = "*"))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar * (sigx24 *
    ((1 + lambda) * sigx7/lwu + 2 * (pi * Wz * sigx8) - dv2/ewvsr2)/(pi * sigx8 *
    ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmnsftslnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  epsiuv <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  depsi <- dnorm(-epsiuv, 0, 1)
  pepsi <- pnorm(-epsiuv)
  euv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  luvepsi <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  dluvepsi <- dnorm(-luvepsi, 0, 1)
  pluvepsi <- pnorm(-luvepsi)
  dv2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  lwu <- (1 + 2 * (lambda * ewusr))
  sigx1 <- (depsi/ewvsr1 - pepsi/ewusr)
  sigx2 <- (dluvepsi/ewvsr1 - (1 + lambda) * pluvepsi/ewusr)
  sigx3 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  sigx4 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx5 <- (2 * (sigx1 * euv) - sigx2 * sigx3)
  sigx6 <- ((1 + lambda) * sigx5 * pwZ/lwu + S * (1 - pwZ) * dv2 * (epsilon)/ewvsr2^3)
  sigx7 <- (2 * (euv * pepsi) - sigx3 * pluvepsi)
  sigx8 <- ((1 - pwZ) * dv2/ewvsr2 + (1 + lambda) * sigx7 * pwZ/lwu)
  sigx9 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx10 <- (0.5 * (depsi * ewvsr1/ewusr) - sigx9 * pepsi)
  sigx11 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv1/(2 *
    ewu)^2))
  sigx12 <- (0.5 * (dluvepsi * ewvsr1/ewusr) - sigx11 * pluvepsi)
  sigx13 <- (2 * (sigx10 * euv) - (sigx12 * (1 + lambda) * sigx3 + lambda * sigx7 *
    ewusr/lwu))
  sigx14 <- (euv * (ewv1 * pepsi/(2 * ewu) - sigx4 * depsi))
  sigx15 <- ((1 + lambda) * ewvsr1/ewusr)
  sigx16 <- (0.5 * sigx15 - 0.5 * (S * (epsilon)/ewvsr1))
  sigx17 <- (2 * sigx14 - ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 *
    dluvepsi) * sigx3)
  sigx18 <- (0.5 * (S^2 * dv2 * (epsilon)^2/ewvsr2^2) - 0.5 * dv2)
  sigx19 <- ((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr)
  sigx20 <- (sigx19 * pluvepsi - dluvepsi * ewvsr1/ewusr)
  sigx21 <- (sigx20 * (1 + lambda) + pluvepsi)
  sigx22 <- (sigx21 * sigx3 + 2 * ((1 + lambda) * sigx7 * ewusr/lwu))
  sigx23 <- (2 * (euv * pepsi) - sigx22)
  sigx24 <- ((1 + lambda) * sigx7/lwu - dv2/ewvsr2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((1 - pwZ) * dv2 * (S^2 * (epsilon)^2/ewvsr2^2 - 1)/ewvsr2^3 + (1 + lambda) *
      (2 * (((epsiuv/ewvsr1 - 1/ewusr) * depsi/ewvsr1 - sigx1/ewusr) * euv) -
        ((luvepsi/ewvsr1 - (1 + lambda)/ewusr) * dluvepsi/ewvsr1 - (1 + lambda) *
          sigx2/ewusr) * sigx3) * pwZ/lwu - sigx6^2/sigx8)/sigx8, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv1/(2 * ewu)^2)) * pepsi - 0.5 * (depsi * ewvsr1/ewusr))/ewusr + (0.5 *
      (epsiuv/ewusr) - sigx9/ewvsr1) * depsi) * euv) - (((0.5 * (luvepsi/ewusr) -
      sigx11/ewvsr1) * dluvepsi + (0.5 * pluvepsi - sigx12 * (1 + lambda))/ewusr) *
      (1 + lambda) * sigx3 + lambda * sigx5 * ewusr/lwu))/(sigx8 * lwu) - sigx6 *
      lwu * sigx13/(sigx8 * lwu)^2) * (1 + lambda) * pwZ, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((depsi * (ewv1/(2 * ewu) - (sigx4 *
      epsiuv + 0.5))/ewvsr1 - (ewv1 * pepsi/(2 * ewu) - sigx4 * depsi)/ewusr) *
      euv) - (((1 + lambda)^2 * ewv1/(2 * ewu) - (luvepsi * sigx16 + 0.5)) *
      dluvepsi/ewvsr1 - ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 *
      dluvepsi) * (1 + lambda)/ewusr) * sigx3)/(sigx8 * lwu) - sigx6 * lwu *
      sigx17/(sigx8 * lwu)^2) * (1 + lambda) * pwZ, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (1 - pwZ) * (S * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2 -
      2) - 0.5) * dv2 * (epsilon)/(sigx8 * ewvsr2^3) - sigx6 * sigx18 * ewvsr2/(sigx8 *
      ewvsr2)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * (sigx1 * euv) - ((((sigx19/ewvsr1 - luvepsi/ewusr) *
      dluvepsi - (sigx20 * (1 + lambda) + 2 * pluvepsi)/ewusr) * (1 + lambda) +
      dluvepsi/ewvsr1) * sigx3 + 2 * ((1 + lambda) * sigx5 * ewusr/lwu)))/(sigx8 *
      lwu) - sigx6 * lwu * sigx23/(sigx8 * lwu)^2) * pwZ, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((1 +
    lambda) * sigx5/lwu - (sigx6 * sigx24/sigx8 + S * dv2 * (epsilon)/ewvsr2^3)) *
    dwZ/sigx8, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr1 * epsiuv/ewusr) -
      0.5) - 0.5 * sigx9) * depsi * ewvsr1/ewusr - (sigx10 * sigx9 + (2 * ((1 -
      8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) *
      pepsi)) * euv) - (((0.5 * (0.5 * (luvepsi * (1 + lambda) * ewvsr1/ewusr) -
      0.5) - 0.5 * (sigx11 * (1 + lambda))) * dluvepsi * ewvsr1/ewusr - (sigx12 *
      sigx11 * (1 + lambda) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) *
      ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pluvepsi)) *
      (1 + lambda) * sigx3 + lambda * ((0.5 - lambda * ewusr/lwu) * sigx7 +
      2 * (sigx10 * euv) - sigx12 * (1 + lambda) * sigx3) * ewusr/lwu))/(sigx8 *
      lwu) - ((1 + lambda) * sigx13 * pwZ + lambda * sigx8 * ewusr) * sigx13/(sigx8 *
      lwu)^2) * (1 + lambda) * pwZ, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((2 * (((depsi *
    ewvsr1/(4 * (ewu * ewusr)) - 2 * (ewu * pepsi/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx4 * epsiuv) - 0.25) * depsi * ewvsr1/ewusr + sigx9 * (ewv1 * pepsi/(2 *
    ewu) - sigx4 * depsi))) * euv) - (((1 + lambda) * dluvepsi * ewvsr1/(4 *
    (ewu * ewusr)) - 2 * (ewu * pluvepsi/(2 * ewu)^2)) * (1 + lambda) * ewv1 -
    (((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 * dluvepsi) * sigx11 +
      (0.5 * (luvepsi * sigx16) - 0.25) * dluvepsi * ewvsr1/ewusr)) * (1 +
    lambda) * sigx3)/(sigx8 * lwu) - ((1 + lambda) * sigx13 * pwZ + lambda *
    sigx8 * ewusr) * sigx17/(sigx8 * lwu)^2) * (1 + lambda) * pwZ, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx18 * (1 - pwZ) * (1 + lambda) * sigx13 * ewvsr2 * pwZ/((sigx8 * ewvsr2)^2 *
      lwu)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx10 * euv) - ((((0.5 * sigx19 - 0.5 *
      (luvepsi * ewvsr1/ewusr)) * (1 + lambda) + 1) * dluvepsi * ewvsr1/ewusr -
      (sigx21 * sigx11 + ((1 + lambda) * ewv1/ewu + 0.5 * (S * (epsilon)/ewusr)) *
        pluvepsi)) * sigx3 + 2 * (((0.5 - lambda * ewusr/lwu) * sigx7 + 2 *
      (sigx10 * euv) - sigx12 * (1 + lambda) * sigx3) * ewusr/lwu)) * (1 +
      lambda))/(sigx8 * lwu) - ((1 + lambda) * sigx13 * pwZ + lambda * sigx8 *
      ewusr) * sigx23/(sigx8 * lwu)^2) * pwZ, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * (1 - sigx24 * pwZ/sigx8) * (1 + lambda) * sigx13 * dwZ/(sigx8 *
      lwu), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((2 * ((((ewv1 * pepsi/(2 * ewu) - sigx4 * depsi)/2 + (pepsi - sigx4 * depsi)/2) *
      ewv1/ewu - (0.25 * (ewvsr1/ewusr) + 0.25 * (S * (epsilon)/ewvsr1) - sigx4^2 *
      epsiuv) * depsi) * euv) - ((((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) -
      sigx16 * dluvepsi)/2 + (pluvepsi - sigx16 * dluvepsi)/2) * (1 + lambda)^2 *
      ewv1/ewu - (0.25 * sigx15 + 0.25 * (S * (epsilon)/ewvsr1) - luvepsi *
      sigx16^2) * dluvepsi) * sigx3)/(sigx8 * lwu) - (1 + lambda) * sigx17^2 *
      pwZ/(sigx8 * lwu)^2) * (1 + lambda) * pwZ, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx18 * (1 - pwZ) * (1 + lambda) * sigx17 * ewvsr2 * pwZ/((sigx8 *
      ewvsr2)^2 * lwu)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((2 *
    sigx14 - ((((sigx21/2 + pluvepsi) * (1 + lambda) * ewv1/ewu - (sigx19 * sigx16 +
    (0.5 - luvepsi * sigx16) * ewvsr1/ewusr) * dluvepsi) * (1 + lambda) - sigx16 *
    dluvepsi) * sigx3 + 2 * ((1 + lambda) * sigx17 * ewusr/lwu)))/(sigx8 * lwu) -
    (1 + lambda) * sigx17 * sigx23 * pwZ/(sigx8 * lwu)^2) * pwZ, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx24 * pwZ/sigx8) * (1 + lambda) * sigx17 *
      dwZ/(sigx8 * lwu), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - pwZ) * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2) -
      1) - 0.25) * dv2 * (epsilon)^2/(sigx8 * ewvsr2^3) - (sigx18 * (1 - pwZ) +
      0.5 * (sigx8 * ewvsr2)) * sigx18/(sigx8 * ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx18 * (1 - pwZ) * lwu * sigx23 * pwZ/((sigx8 * lwu)^2 * ewvsr2)), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx24 * (1 - pwZ)/sigx8 + 1) * sigx18 * dwZ/(sigx8 *
      ewvsr2)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    1)] <- sum(wHvar * (-(((((((luvepsi * ewvsr1 - S * (epsilon))/ewusr - (1 +
    lambda) * ewv1/ewu) * dluvepsi * ewvsr1/ewusr + ewv1 * pluvepsi/ewu) * (1 +
    lambda) + (sigx20 * (1 + lambda) + 2 * pluvepsi) * sigx19 - 2 * (dluvepsi *
    ewvsr1/ewusr)) * sigx3 + 2 * (sigx23 * ewusr/lwu))/(sigx8 * lwu) + (sigx23 *
    pwZ + 2 * (sigx8 * ewusr)) * sigx23/(sigx8 * lwu)^2) * pwZ)))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * (1 - sigx24 * pwZ/sigx8) * sigx23 * dwZ/(sigx8 * lwu), FUN = "*"))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar * ((sigx24 *
    dwZ/sigx8 + Wz) * sigx24 * dwZ/sigx8), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmnsftslnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  epsiuv <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  depsi <- dnorm(-epsiuv, 0, 1)
  pepsi <- pnorm(-epsiuv)
  euv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  luvepsi <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  dluvepsi <- dnorm(-luvepsi, 0, 1)
  pluvepsi <- pnorm(-luvepsi)
  dv2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  lwu <- (1 + 2 * (lambda * ewusr))
  sigx1 <- (depsi/ewvsr1 - pepsi/ewusr)
  sigx2 <- (dluvepsi/ewvsr1 - (1 + lambda) * pluvepsi/ewusr)
  sigx3 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) * (1 + lambda))
  sigx4 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx5 <- (2 * (sigx1 * euv) - sigx2 * sigx3)
  sigx6 <- ((1 - prZ) * (1 + lambda) * sigx5/lwu + S * dv2 * prZ * (epsilon)/ewvsr2^3)
  sigx7 <- (2 * (euv * pepsi) - sigx3 * pluvepsi)
  sigx8 <- ((1 - prZ) * (1 + lambda) * sigx7/lwu + dv2 * prZ/ewvsr2)
  sigx9 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx10 <- (0.5 * (depsi * ewvsr1/ewusr) - sigx9 * pepsi)
  sigx11 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * ewu * ewv1/(2 *
    ewu)^2))
  sigx12 <- (0.5 * (dluvepsi * ewvsr1/ewusr) - sigx11 * pluvepsi)
  sigx13 <- (2 * (sigx10 * euv) - (sigx12 * (1 + lambda) * sigx3 + lambda * sigx7 *
    ewusr/lwu))
  sigx14 <- (euv * (ewv1 * pepsi/(2 * ewu) - sigx4 * depsi))
  sigx15 <- ((1 + lambda) * ewvsr1/ewusr)
  sigx16 <- (0.5 * sigx15 - 0.5 * (S * (epsilon)/ewvsr1))
  sigx17 <- (2 * sigx14 - ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 *
    dluvepsi) * sigx3)
  sigx18 <- (0.5 * (S^2 * dv2 * (epsilon)^2/ewvsr2^2) - 0.5 * dv2)
  sigx19 <- ((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr)
  sigx20 <- (sigx19 * pluvepsi - dluvepsi * ewvsr1/ewusr)
  sigx21 <- (sigx20 * (1 + lambda) + pluvepsi)
  sigx22 <- (sigx21 * sigx3 + 2 * ((1 + lambda) * sigx7 * ewusr/lwu))
  sigx23 <- (2 * (euv * pepsi) - sigx22)
  sigx24 <- ((1 + lambda) * sigx7/lwu - dv2/ewvsr2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((1 - prZ) * (1 + lambda) * (2 * (((epsiuv/ewvsr1 - 1/ewusr) * depsi/ewvsr1 -
      sigx1/ewusr) * euv) - ((luvepsi/ewvsr1 - (1 + lambda)/ewusr) * dluvepsi/ewvsr1 -
      (1 + lambda) * sigx2/ewusr) * sigx3)/lwu + dv2 * prZ * (S^2 * (epsilon)^2/ewvsr2^2 -
      1)/ewvsr2^3 - sigx6^2/sigx8)/sigx8, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * ((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv1/(2 * ewu)^2)) * pepsi - 0.5 * (depsi * ewvsr1/ewusr))/ewusr + (0.5 *
      (epsiuv/ewusr) - sigx9/ewvsr1) * depsi) * euv) - (((0.5 * (luvepsi/ewusr) -
      sigx11/ewvsr1) * dluvepsi + (0.5 * pluvepsi - sigx12 * (1 + lambda))/ewusr) *
      (1 + lambda) * sigx3 + lambda * sigx5 * ewusr/lwu))/(sigx8 * lwu) - sigx6 *
      lwu * sigx13/(sigx8 * lwu)^2) * (1 - prZ) * (1 + lambda), FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((depsi * (ewv1/(2 * ewu) - (sigx4 *
      epsiuv + 0.5))/ewvsr1 - (ewv1 * pepsi/(2 * ewu) - sigx4 * depsi)/ewusr) *
      euv) - (((1 + lambda)^2 * ewv1/(2 * ewu) - (luvepsi * sigx16 + 0.5)) *
      dluvepsi/ewvsr1 - ((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 *
      dluvepsi) * (1 + lambda)/ewusr) * sigx3)/(sigx8 * lwu) - sigx6 * lwu *
      sigx17/(sigx8 * lwu)^2) * (1 - prZ) * (1 + lambda), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * prZ * (S * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2 -
      2) - 0.5) * dv2 * (epsilon)/(sigx8 * ewvsr2^3) - sigx6 * sigx18 * ewvsr2/(sigx8 *
      ewvsr2)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((2 * (sigx1 * euv) - ((((sigx19/ewvsr1 - luvepsi/ewusr) *
      dluvepsi - (sigx20 * (1 + lambda) + 2 * pluvepsi)/ewusr) * (1 + lambda) +
      dluvepsi/ewvsr1) * sigx3 + 2 * ((1 + lambda) * sigx5 * ewusr/lwu)))/(sigx8 *
      lwu) - sigx6 * lwu * sigx23/(sigx8 * lwu)^2) * (1 - prZ), FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((1 +
    lambda) * sigx5/lwu - (sigx6 * sigx24/sigx8 + S * dv2 * (epsilon)/ewvsr2^3)) *
    prZ * ewz/sigx8, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr1 * epsiuv/ewusr) -
      0.5) - 0.5 * sigx9) * depsi * ewvsr1/ewusr - (sigx10 * sigx9 + (2 * ((1 -
      8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) *
      pepsi)) * euv) - (((0.5 * (0.5 * (luvepsi * (1 + lambda) * ewvsr1/ewusr) -
      0.5) - 0.5 * (sigx11 * (1 + lambda))) * dluvepsi * ewvsr1/ewusr - (sigx12 *
      sigx11 * (1 + lambda) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) *
      ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pluvepsi)) *
      (1 + lambda) * sigx3 + lambda * ((0.5 - lambda * ewusr/lwu) * sigx7 +
      2 * (sigx10 * euv) - sigx12 * (1 + lambda) * sigx3) * ewusr/lwu))/(sigx8 *
      lwu) - ((1 - prZ) * (1 + lambda) * sigx13 + lambda * sigx8 * ewusr) *
      sigx13/(sigx8 * lwu)^2) * (1 - prZ) * (1 + lambda), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((2 * (((depsi *
    ewvsr1/(4 * (ewu * ewusr)) - 2 * (ewu * pepsi/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx4 * epsiuv) - 0.25) * depsi * ewvsr1/ewusr + sigx9 * (ewv1 * pepsi/(2 *
    ewu) - sigx4 * depsi))) * euv) - (((1 + lambda) * dluvepsi * ewvsr1/(4 *
    (ewu * ewusr)) - 2 * (ewu * pluvepsi/(2 * ewu)^2)) * (1 + lambda) * ewv1 -
    (((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) - sigx16 * dluvepsi) * sigx11 +
      (0.5 * (luvepsi * sigx16) - 0.25) * dluvepsi * ewvsr1/ewusr)) * (1 +
    lambda) * sigx3)/(sigx8 * lwu) - ((1 - prZ) * (1 + lambda) * sigx13 + lambda *
    sigx8 * ewusr) * sigx17/(sigx8 * lwu)^2) * (1 - prZ) * (1 + lambda), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx18 * (1 - prZ) * (1 + lambda) * sigx13 * prZ * ewvsr2/((sigx8 * ewvsr2)^2 *
      lwu)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx10 * euv) - ((((0.5 * sigx19 - 0.5 *
      (luvepsi * ewvsr1/ewusr)) * (1 + lambda) + 1) * dluvepsi * ewvsr1/ewusr -
      (sigx21 * sigx11 + ((1 + lambda) * ewv1/ewu + 0.5 * (S * (epsilon)/ewusr)) *
        pluvepsi)) * sigx3 + 2 * (((0.5 - lambda * ewusr/lwu) * sigx7 + 2 *
      (sigx10 * euv) - sigx12 * (1 + lambda) * sigx3) * ewusr/lwu)) * (1 +
      lambda))/(sigx8 * lwu) - ((1 - prZ) * (1 + lambda) * sigx13 + lambda *
      sigx8 * ewusr) * sigx23/(sigx8 * lwu)^2) * (1 - prZ), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * (1 - sigx24 * (1 - prZ)/sigx8) * (1 + lambda) * sigx13 *
      prZ * ewz/(sigx8 * lwu), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((2 * ((((ewv1 * pepsi/(2 * ewu) - sigx4 * depsi)/2 + (pepsi - sigx4 * depsi)/2) *
      ewv1/ewu - (0.25 * (ewvsr1/ewusr) + 0.25 * (S * (epsilon)/ewvsr1) - sigx4^2 *
      epsiuv) * depsi) * euv) - ((((1 + lambda)^2 * ewv1 * pluvepsi/(2 * ewu) -
      sigx16 * dluvepsi)/2 + (pluvepsi - sigx16 * dluvepsi)/2) * (1 + lambda)^2 *
      ewv1/ewu - (0.25 * sigx15 + 0.25 * (S * (epsilon)/ewvsr1) - luvepsi *
      sigx16^2) * dluvepsi) * sigx3)/(sigx8 * lwu) - (1 - prZ) * (1 + lambda) *
      sigx17^2/(sigx8 * lwu)^2) * (1 - prZ) * (1 + lambda), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx18 * (1 - prZ) * (1 + lambda) * sigx17 * prZ * ewvsr2/((sigx8 *
      ewvsr2)^2 * lwu)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((2 *
    sigx14 - ((((sigx21/2 + pluvepsi) * (1 + lambda) * ewv1/ewu - (sigx19 * sigx16 +
    (0.5 - luvepsi * sigx16) * ewvsr1/ewusr) * dluvepsi) * (1 + lambda) - sigx16 *
    dluvepsi) * sigx3 + 2 * ((1 + lambda) * sigx17 * ewusr/lwu)))/(sigx8 * lwu) -
    (1 - prZ) * (1 + lambda) * sigx17 * sigx23/(sigx8 * lwu)^2) * (1 - prZ),
    FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx24 * (1 - prZ)/sigx8) * (1 + lambda) *
      sigx17 * prZ * ewz/(sigx8 * lwu), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prZ * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2) -
      1) - 0.25) * dv2 * (epsilon)^2/(sigx8 * ewvsr2^3) - (sigx18 * prZ + 0.5 *
      (sigx8 * ewvsr2)) * sigx18/(sigx8 * ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx18 * (1 - prZ) * lwu * sigx23 * prZ/((sigx8 * lwu)^2 * ewvsr2)), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx24 * prZ/sigx8 + 1) * sigx18 * prZ * ewz/(sigx8 *
      ewvsr2)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    1)] <- sum(wHvar * (-(((((((luvepsi * ewvsr1 - S * (epsilon))/ewusr - (1 +
    lambda) * ewv1/ewu) * dluvepsi * ewvsr1/ewusr + ewv1 * pluvepsi/ewu) * (1 +
    lambda) + (sigx20 * (1 + lambda) + 2 * pluvepsi) * sigx19 - 2 * (dluvepsi *
    ewvsr1/ewusr)) * sigx3 + 2 * (sigx23 * ewusr/lwu))/(sigx8 * lwu) + ((1 -
    prZ) * sigx23 + 2 * (sigx8 * ewusr)) * sigx23/(sigx8 * lwu)^2) * (1 - prZ))))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * (1 - sigx24 * (1 - prZ)/sigx8) * sigx23 * prZ * ewz/(sigx8 *
      lwu), FUN = "*"))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar + 1)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar * sigx24 *
    (1 - (sigx24 * prZ/sigx8 + 1) * ewz) * prZ * ewz/sigx8, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf truncated_skewed_laplace-normal distribution
#' @param start starting value for optimization
#' @param randStart if random starting values should be used
#' @param sdStart std. Error for random draws for starting values
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
# Same sigma_v

## logit specification class membership
zisftslnormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisftslnormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisftslnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisftslnormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisftslnormlike_logit,
    grad = cgradzisftslnormlike_logit, hess = chesszisftslnormlike_logit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisftslnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo, eval.max = itermax,
        rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisftslnormlike_logit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chesszisftslnormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisftslnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisftslnormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisftslnormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTSL = initTSL))
}

## cauchit specification class membership
zisftslnormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisftslnormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisftslnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisftslnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisftslnormlike_cauchit,
    grad = cgradzisftslnormlike_cauchit, hess = chesszisftslnormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisftslnormlike_cauchit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisftslnormlike_cauchit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chesszisftslnormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisftslnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisftslnormlike_cauchit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisftslnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTSL = initTSL))
}

## probit specification class membership
zisftslnormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisftslnormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisftslnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisftslnormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisftslnormlike_probit,
    grad = cgradzisftslnormlike_probit, hess = chesszisftslnormlike_probit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisftslnormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisftslnormlike_probit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chesszisftslnormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisftslnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisftslnormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisftslnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTSL = initTSL))
}

## cloglog specification class membership
zisftslnormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisftslnormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisftslnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisftslnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisftslnormlike_cloglog,
    grad = cgradzisftslnormlike_cloglog, hess = chesszisftslnormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisftslnormlike_cloglog(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisftslnormlike_cloglog(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chesszisftslnormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisftslnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisftslnormlike_cloglog(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisftslnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTSL = initTSL))
}

# Different sigma_v

## logit specification class membership
mnsftslnormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsftslnormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsftslnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsftslnormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsftslnormlike_logit,
    grad = cgradmnsftslnormlike_logit, hess = chessmnsftslnormlike_logit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsftslnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsftslnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo, eval.max = itermax,
        rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsftslnormlike_logit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessmnsftslnormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsftslnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsftslnormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsftslnormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTSL = initTSL))
}

## cauchit specification class membership
mnsftslnormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsftslnormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsftslnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsftslnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsftslnormlike_cauchit,
    grad = cgradmnsftslnormlike_cauchit, hess = chessmnsftslnormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsftslnormlike_cauchit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsftslnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsftslnormlike_cauchit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessmnsftslnormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsftslnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsftslnormlike_cauchit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsftslnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTSL = initTSL))
}

## probit specification class membership
mnsftslnormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsftslnormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsftslnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsftslnormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsftslnormlike_probit,
    grad = cgradmnsftslnormlike_probit, hess = chessmnsftslnormlike_probit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsftslnormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsftslnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsftslnormlike_probit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessmnsftslnormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsftslnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsftslnormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsftslnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTSL = initTSL))
}

## cloglog specification class membership
mnsftslnormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsftslnormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsftslnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsftslnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsftslnormlike_cloglog,
    grad = cgradmnsftslnormlike_cloglog, hess = chessmnsftslnormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsftslnormlike_cloglog(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsftslnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsftslnormlike_cloglog(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessmnsftslnormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsftslnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsftslnormlike_cloglog(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsftslnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTSL = initTSL))
}

# Conditional efficiencies estimation ----------
#' efficiencies for truncated skewed laplace-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v

## logit specification class membership
czisftslnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 2):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar + 1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a - exp(Wv/2)) -
      exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) * pnorm(a +
      exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cauchit specification class membership
czisftslnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 2):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar + 1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a - exp(Wv/2)) -
      exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) * pnorm(a +
      exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## probit specification class membership
czisftslnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 2):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar + 1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a - exp(Wv/2)) -
      exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) * pnorm(a +
      exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cloglog specification class membership
czisftslnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 2):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar + 1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a - exp(Wv/2)) -
      exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) * pnorm(a +
      exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

# Different sigma_v

## logit specification class membership
cmnsftslnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A) * exp(-a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a - exp(Wv1/2)) -
      exp(B) * exp(-b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b - exp(Wv1/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A) * exp(a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a +
      exp(Wv1/2)) - exp(B) * exp(b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b + exp(Wv1/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cauchit specification class membership
cmnsftslnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A) * exp(-a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a - exp(Wv1/2)) -
      exp(B) * exp(-b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b - exp(Wv1/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A) * exp(a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a +
      exp(Wv1/2)) - exp(B) * exp(b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b + exp(Wv1/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## probit specification class membership
cmnsftslnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A) * exp(-a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a - exp(Wv1/2)) -
      exp(B) * exp(-b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b - exp(Wv1/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A) * exp(a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a +
      exp(Wv1/2)) - exp(B) * exp(b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b + exp(Wv1/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cloglog specification class membership
cmnsftslnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A) * exp(-a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a - exp(Wv1/2)) -
      exp(B) * exp(-b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b - exp(Wv1/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A) * exp(a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a +
      exp(Wv1/2)) - exp(B) * exp(b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b + exp(Wv1/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, odRatio = odRatio,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for zisf tsl-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v

## logit specification class membership
czisfmargtslnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 4 * lambda + 2 *
    lambda^2)/((1 + lambda) * (1 + 2 * lambda)), nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargtslnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 8 * lambda + 16 *
    lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 + lambda)^2 * (1 + 2 * lambda)^2),
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
czisfmargtslnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 4 * lambda + 2 *
    lambda^2)/((1 + lambda) * (1 + 2 * lambda)), nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargtslnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 8 * lambda + 16 *
    lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 + lambda)^2 * (1 + 2 * lambda)^2),
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
czisfmargtslnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 4 * lambda + 2 *
    lambda^2)/((1 + lambda) * (1 + 2 * lambda)), nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargtslnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 8 * lambda + 16 *
    lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 + lambda)^2 * (1 + 2 * lambda)^2),
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
czisfmargtslnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 4 * lambda + 2 *
    lambda^2)/((1 + lambda) * (1 + 2 * lambda)), nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargtslnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 8 * lambda + 16 *
    lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 + lambda)^2 * (1 + 2 * lambda)^2),
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# Different sigma_v

## logit specification class membership
cmnsfmargtslnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 4 * lambda + 2 *
    lambda^2)/((1 + lambda) * (1 + 2 * lambda)), nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargtslnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 8 * lambda + 16 *
    lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 + lambda)^2 * (1 + 2 * lambda)^2),
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmnsfmargtslnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 4 * lambda + 2 *
    lambda^2)/((1 + lambda) * (1 + 2 * lambda)), nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargtslnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 8 * lambda + 16 *
    lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 + lambda)^2 * (1 + 2 * lambda)^2),
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmnsfmargtslnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 4 * lambda + 2 *
    lambda^2)/((1 + lambda) * (1 + 2 * lambda)), nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargtslnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 8 * lambda + 16 *
    lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 + lambda)^2 * (1 + 2 * lambda)^2),
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmnsfmargtslnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 4 * lambda + 2 *
    lambda^2)/((1 + lambda) * (1 + 2 * lambda)), nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargtslnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 + 8 * lambda + 16 *
    lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 + lambda)^2 * (1 + 2 * lambda)^2),
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
