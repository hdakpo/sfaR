################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Contaminated Noise Stochastic Frontier Model                          #
# Two types: - Common inefficiency component (sigma_u)                         #
#            - Mixture composed error (mcesf)                                  # 
# Link functions: - logit exp(theta * Z)/(1 + exp(theta * Z))                  #
#                 - cauchit 1/pi * atan(theta * Z) + 1/2                       #
#                 - probit pnorm(theta * Z)                                    #
#                 - cloglog 1 - exp(-exp(theta * Z))                           #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf halfnormal-normal distribution
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
# same sigma_u

## logit specification class membership
ccnsfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
ccnsfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
ccnsfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
ccnsfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# different sigma_u

## logit specification class membership
cmcesfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu1) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu2) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
cmcesfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu1) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu2) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
cmcesfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu1) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu2) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
cmcesfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu1) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu2) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for cnsf halfnormal-normal distribution
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
# same sigma_u
cstcnsfhalfnorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- csthalfnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initHalf <- NULL
  } else {
    cat("Initialization: SFA halfnormal - normal distribution...\n")
    initHalf <- maxLik::maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradhalfnormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initHalf$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), 0.95 * Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 1.05 *
    Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("CNSF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# different sigma_u
cstmcesfhalfnorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- csthalfnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initHalf <- NULL
  } else {
    cat("Initialization: SFA halfnormal - normal distribution...\n")
    initHalf <- maxLik::maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradhalfnormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initHalf$estimate
  }
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0,
    nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1),
    0.95 * Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar +
      2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_", colnames(vHvar)),
    paste0("CNSF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for cnsf halfnormal-normal distribution
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
# same sigma_u

## logit specification class membership
cgradcnsfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewz <- exp(Wz)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  sigma_sq1 <- ewu + ewv1
  sigma_sq2 <- ewu + ewv2
  sigmastar1 <- sqrt(ewu * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu * (epsilon)/ssq1)
  musig2 <- (S * ewu * (epsilon)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  depsi1 <- dnorm(S * (epsilon)/sqrt(sigma_sq1))
  depsi2 <- dnorm(S * (epsilon)/sqrt(sigma_sq2))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewvsq1 <- (1 - ewv1/(sigma_sq1))
  ewvsq2 <- (1 - ewv2/(sigma_sq2))
  sigx1_1 <- (0.5 * (ewvsq1 * ewu/sigmastar1) + sigmastar1)
  sigx1_2 <- (0.5 * (ewvsq2 * ewu/sigmastar2) + sigmastar2)
  sigx2_1 <- (S * depsi1 * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx2_2 <- (S * depsi2 * pmusig2 * (epsilon)/(sigma_sq2)^2)
  sigx3_1 <- (sigx1_1 * dmusig1 * depsi1 * ewu/ssq1^2 + 0.5 * sigx2_1)
  sigx3_2 <- (sigx1_2 * dmusig2 * depsi2 * ewu/ssq2^2 + 0.5 * sigx2_2)
  sigx4_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + S * depsi1 * pmusig1 * (epsilon))
  sigx4_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + S * depsi2 * pmusig2 * (epsilon))
  sisi1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sisi2 <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx5_1 <- (sigx4_1 * ewz/sisi1)
  sigx5_2 <- (prC * sigx4_2/sisi2)
  wzsq1 <- (wzdeno * sqrt(sigma_sq1))
  wzsq2 <- (wzdeno * sqrt(sigma_sq2))
  sigx6_1 <- (depsi1 * ewz * pmusig1/wzsq1)
  sigx6_2 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2))
  prU1 <- (1 - ewu/(sigma_sq1))
  prU2 <- (1 - ewu/(sigma_sq2))
  sigx7_1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (wzdeno * depsi1 * pmusig1/wzsq1^2)
  sigx8_2 <- (prC * depsi2 * pmusig2/wzsq2)
  sigx9_1 <- (1/wzsq1 - ewz * sqrt(sigma_sq1)/wzsq1^2)
  sigx9_2 <- ((2 * sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq2))
  sigx10_1 <- (sigx9_1 * depsi1 * pmusig1)
  sigx10_2 <- (depsi2 * pmusig2/(sigma_sq2))
  ssqx1 <- (1/ssq1 - sigx7_1 * ewu/ssq1^2)
  ssqx2 <- (1/ssq2 - sigx7_2 * ewu/ssq2^2)
  sigx11_1 <- (0.5 * sigx2_1 - ssqx1 * dmusig1 * depsi1)
  sigx11_2 <- (0.5 * sigx2_2 - ssqx2 * dmusig2 * depsi2)
  sigx12_1 <- (S * sigx11_1 * (epsilon)/wzdeno - 0.5 * sigx8_1)
  sigx12_2 <- (S * sigx11_2 * (epsilon) - 0.5 * sigx10_2)
  sigx13_1 <- (S * sigx3_1 * (epsilon)/wzdeno - 0.5 * sigx8_1)
  sigx13_2 <- (S * sigx3_2 * (epsilon) - 0.5 * sigx10_2)
  sigx14_1 <- (ewz * sigx12_1/sqrt(sigma_sq1))
  sigx14_2 <- (prC * sigx12_2/sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (2 * sigx5_2 + 2 * sigx5_1)/(2 *
    sigx6_2 + 2 * sigx6_1), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (2 *
    sigx14_2 + 2 * sigx14_1) * ewu/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = 2 * (ewv1 * ewz * sigx13_1/((2 * sigx6_2 + 2 * sigx6_1) *
      sqrt(sigma_sq1))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * (prC *
    ewv2 * sigx13_2/sigx9_2), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = (2 *
    sigx10_1 - 2 * sigx8_2) * ewz/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradcnsfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu <- 0
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/sqrt(sigma_sq2) + ewz1 * depsi1 * pmusig1/sqrt(sigma_sq1))
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- (ewz2 * sigx4_2/ssqq2 + ewz1 * sigx4_1/ssqq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (ewz2 * sigx11_2/sqrt(sigma_sq2) + sigx11_1 * ewz1/sqrt(sigma_sq1))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (sigx14 * ewu/sigx2 + 0.5 * (mu * dmu/(ewusr *
      pmu))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_1 * ewz1 *
      ewv1/sigx17_1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_2 *
      ewz2 * ewv2/sigx17_2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx9/(pi *
      sigx2 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradcnsfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu <- 0
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/sqrt(sigma_sq2) + depsi1 * pmusig1 * pwZ/sqrt(sigma_sq1))
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- ((1 - pwZ) * sigx4_2/ssqq2 + sigx4_1 * pwZ/ssqq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (sigx11_1 * pwZ/sqrt(sigma_sq1) + sigx11_2 * (1 - pwZ)/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (sigx14 * ewu/sigx2 + 0.5 * (mu * dmu/(ewusr *
      pmu))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_1 * pwZ *
      ewv1/sigx17_1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_2 *
      (1 - pwZ) * ewv2/sigx17_2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx9 *
      dwZ/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradcnsfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu <- 0
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/sqrt(sigma_sq1) + depsi2 * prZ * pmusig2/sqrt(sigma_sq2))
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- ((1 - prZ) * sigx4_1/ssqq1 + sigx4_2 * prZ/ssqq2)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (sigx11_1 * (1 - prZ)/sqrt(sigma_sq1) + sigx11_2 * prZ/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (sigx14 * ewu/sigx2 + 0.5 * (mu * dmu/(ewusr *
      pmu))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_1 * (1 -
      prZ) * ewv1/sigx17_1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_2 *
      prZ * ewv2/sigx17_2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx9 *
      prZ * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u

## logit specification class membership
cgradmcesfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu1 * (epsilon)/ssq1)
  musig2 <- (S * ewu2 * (epsilon)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsi1 <- S * (epsilon)/sqrt(sigma_sq1)
  epsi2 <- S * (epsilon)/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  sxq <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 * pmusig1 * (epsilon))
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + S * depsi2 * pmusig2 * (epsilon))
  wzsq1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  wzxsq1 <- (wzdeno * sqrt(sigma_sq1))
  sigx2 <- (2 * (prC * sigx1_2/sxq) + 2 * (sigx1_1 * ewz/wzsq1))
  sigx3_1 <- (depsi1 * ewz * pmusig1/wzxsq1)
  sigx3_2 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx4_1 <- (S * depsi1 * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx4_2 <- (S * depsi2 * pmusig2 * (epsilon)/(sigma_sq2)^2)
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  usq1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  usq2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
  sigx5_1 <- (1/ssq1 - usq1 * ewu1/ssq1^2)
  sigx5_2 <- (1/ssq2 - usq2 * ewu2/ssq2^2)
  sigx6_1 <- (0.5 * sigx4_1 - sigx5_1 * dmusig1 * depsi1)
  sigx6_2 <- (0.5 * sigx4_2 - sigx5_2 * dmusig2 * depsi2)
  s7sq2 <- (depsi2 * pmusig2/(sigma_sq2))
  sigx11_1 <- (wzdeno * depsi1 * pmusig1/wzxsq1^2)
  sigx11_2 <- (prC * depsi2 * pmusig2/(wzdeno * sqrt(sigma_sq2)))
  sigx7_1 <- (S * sigx6_1 * (epsilon)/wzdeno - 0.5 * sigx11_1)
  sigx7_2 <- (S * sigx6_2 * (epsilon) - 0.5 * s7sq2)
  sigx8_1 <- ((2 * sigx3_2 + 2 * sigx3_1) * sqrt(sigma_sq1))
  sigx8_2 <- ((2 * sigx3_2 + 2 * sigx3_1) * sqrt(sigma_sq2))
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  sigx9_1 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  s9sq2 <- (0.5 * (prV2 * ewu2/sigmastar2) + sigmastar2)
  sigx9_2 <- (s9sq2 * dmusig2 * depsi2 * ewu2/ssq2^2 + 0.5 * sigx4_2)
  sigx10_1 <- (sigx9_1 * dmusig1 * depsi1 * ewu1/ssq1^2 + 0.5 * sigx4_1)
  sigx10_2 <- (S * sigx9_2 * (epsilon) - 0.5 * s7sq2)
  s12sq1 <- (1/wzxsq1 - ewz * sqrt(sigma_sq1)/wzxsq1^2)
  sigx12 <- (2 * (s12sq1 * depsi1 * pmusig1) - 2 * sigx11_2)
  sigx18_1 <- (S * sigx10_1 * (epsilon)/wzdeno - 0.5 * sigx11_1)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/(2 * sigx3_2 + 2 *
    sigx3_1), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (ewu1 * ewz *
    sigx7_1/sigx8_1), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (prC *
    ewu2 * sigx7_2/sigx8_2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 *
    (ewv1 * ewz * sigx18_1/sigx8_1), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 *
    (prC * ewv2 * sigx10_2/sigx8_2), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx12 *
    ewz/(2 * sigx3_2 + 2 * sigx3_1), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmcesfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/psq2 + ewz1 * depsi1 * pmusig1/psq1)
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx18 <- (pi * sigx2 * ((Wz)^2 + 1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx11_1 * ewz1/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = sigx11_2 * ewz2/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx16_1 * ewz1 * ewv1/sigx17_1, FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx16_2 * ewz2 * ewv2/sigx17_2, FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx9/sigx18, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmcesfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/psq2 + depsi1 * pmusig1 * pwZ/psq1)
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx11_1 * pwZ/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = sigx11_2 * (1 - pwZ)/sigx2, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx16_1 * pwZ * ewv1/sigx17_1, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx16_2 * (1 - pwZ) * ewv2/sigx17_2, FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx9 * dwZ/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmcesfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/psq1 + depsi2 * prZ * pmusig2/psq2)
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx11_1 * (1 - prZ)/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx11_2 * prZ/sigx2, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx16_1 * (1 - prZ) * ewv1/sigx17_1, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx16_2 * prZ * ewv2/sigx17_2, FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx9 * prZ * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for cnsf halfnormal-normal distribution
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
# same sigma_u

## logit specification class membership
chesscnsfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewz <- exp(Wz)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  sigma_sq1 <- ewu + ewv1
  sigma_sq2 <- ewu + ewv2
  sigmastar1 <- sqrt(ewu * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu * (epsilon)/ssq1)
  musig2 <- (S * ewu * (epsilon)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  depsi1 <- dnorm(S * (epsilon)/sqrt(sigma_sq1))
  depsi2 <- dnorm(S * (epsilon)/sqrt(sigma_sq2))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewvsq1 <- (1 - ewv1/(sigma_sq1))
  ewvsq2 <- (1 - ewv2/(sigma_sq2))
  sigx1_1 <- (0.5 * (ewvsq1 * ewu/sigmastar1) + sigmastar1)
  sigx1_2 <- (0.5 * (ewvsq2 * ewu/sigmastar2) + sigmastar2)
  sigx2_1 <- (S * depsi1 * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx2_2 <- (S * depsi2 * pmusig2 * (epsilon)/(sigma_sq2)^2)
  sigx3_1 <- (sigx1_1 * dmusig1 * depsi1 * ewu/ssq1^2 + 0.5 * sigx2_1)
  sigx3_2 <- (sigx1_2 * dmusig2 * depsi2 * ewu/ssq2^2 + 0.5 * sigx2_2)
  sigx4_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + S * depsi1 * pmusig1 * (epsilon))
  sigx4_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + S * depsi2 * pmusig2 * (epsilon))
  sisi1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sisi2 <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx5_1 <- (sigx4_1 * ewz/sisi1)
  sigx5_2 <- (prC * sigx4_2/sisi2)
  wzsq1 <- (wzdeno * sqrt(sigma_sq1))
  wzsq2 <- (wzdeno * sqrt(sigma_sq2))
  sigx6_1 <- (depsi1 * ewz * pmusig1/wzsq1)
  sigx6_2 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2))
  prU1 <- (1 - ewu/(sigma_sq1))
  prU2 <- (1 - ewu/(sigma_sq2))
  sigx7_1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (wzdeno * depsi1 * pmusig1/wzsq1^2)
  sigx8_2 <- (prC * depsi2 * pmusig2/wzsq2)
  sigx9_1 <- (1/wzsq1 - ewz * sqrt(sigma_sq1)/wzsq1^2)
  sigx9_2 <- ((2 * sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq2))
  sigx10_1 <- (sigx9_1 * depsi1 * pmusig1)
  sigx10_2 <- (depsi2 * pmusig2/(sigma_sq2))
  ssqx1 <- (1/ssq1 - sigx7_1 * ewu/ssq1^2)
  ssqx2 <- (1/ssq2 - sigx7_2 * ewu/ssq2^2)
  sigx11_1 <- (0.5 * sigx2_1 - ssqx1 * dmusig1 * depsi1)
  sigx11_2 <- (0.5 * sigx2_2 - ssqx2 * dmusig2 * depsi2)
  sigx12_1 <- (S * sigx11_1 * (epsilon)/wzdeno - 0.5 * sigx8_1)
  sigx12_2 <- (S * sigx11_2 * (epsilon) - 0.5 * sigx10_2)
  sigx13_1 <- (S * sigx3_1 * (epsilon)/wzdeno - 0.5 * sigx8_1)
  sigx13_2 <- (S * sigx3_2 * (epsilon) - 0.5 * sigx10_2)
  sigx14_1 <- (ewz * sigx12_1/sqrt(sigma_sq1))
  sigx14_2 <- (prC * sigx12_2/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * ewu/sigmastar1 + S * pmusig1 * (epsilon))
  sigx15_2 <- (dmusig2 * ewu/sigmastar2 + S * pmusig2 * (epsilon))
  ssqxx1 <- (S * sigx15_1 * (epsilon)/(sigma_sq1) - pmusig1)
  ssqxx2 <- (S * sigx15_2 * (epsilon)/(sigma_sq2) - pmusig2)
  dd1 <- (depsi1 * ewu/ewv1 + depsi1)
  dd2 <- (depsi2 * ewu/ewv2 + depsi2)
  sigx17_1 <- (0.5 * ssqxx1 - 0.5 * pmusig1)
  sigx17_2 <- (0.5 * ssqxx2 - 0.5 * pmusig2)
  sigx18_1 <- (wzdeno * sigx4_1/(wzsq1^2 * (sigma_sq1)))
  sigx18_2 <- (sigx4_2/(sigma_sq2))
  pesig1 <- (S * pmusig1 * (epsilon)/(sigma_sq1)^2)
  pesig2 <- (S * pmusig2 * (epsilon)/(sigma_sq2)^2)
  sigx19_1 <- (0.5 * pesig1 - ssqx1 * dmusig1)
  sigx19_2 <- (0.5 * pesig2 - ssqx2 * dmusig2)
  sigx20_1 <- (S * sigx19_1 * (epsilon) - 2 * (pmusig1/(sigma_sq1)))
  sigx20_2 <- (S * sigx19_2 * (epsilon) - 2 * (pmusig2/(sigma_sq2)))
  sigx21_1 <- (S * depsi1 * sigx20_1 * (epsilon)/(sigma_sq1)^2)
  sigx21_2 <- (S * depsi2 * sigx20_2 * (epsilon)/(sigma_sq2)^2)
  sigx22_1 <- (S * sigx11_1 * (epsilon) - wzdeno^2 * depsi1 * pmusig1/wzsq1^2)
  sigx22_2 <- (S * sigx11_2 * (epsilon) - depsi2 * pmusig2/(sigma_sq2))
  sigx23_1 <- (wzdeno * sigx22_1/wzsq1^2)
  sigx23_2 <- (wzdeno * depsi2 * pmusig2/wzsq2^2)
  sigx24_1 <- (0.5/sqrt(sigma_sq1) - wzdeno^2 * sqrt(sigma_sq1)/wzsq1^2)
  sigx25 <- (2 * sigx10_1 - 2 * sigx8_2)/(2 * sigx6_2 + 2 * sigx6_1)
  sigx26_1 <- (0.5 * (S^2 * sigx9_1 * depsi1 * (epsilon)^2/(sigma_sq1)^2) - (sigx24_1 *
    ewz + 0.5 * (wzdeno/sqrt(sigma_sq1))) * depsi1/wzsq1^2)
  sigx27 <- ((2 * sigx6_2 + 2 * sigx6_1)/sqrt(sigma_sq1))
  sigx28 <- ((2 * sigx6_2 + 2 * sigx6_1)/sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (2 * (prC * (depsi2 * ssqxx2 + S * dmusig2 * dd2 * ewu * (epsilon)/ssq2)/((sigma_sq2) *
    sqrt(sigma_sq2))) + 2 * ((depsi1 * ssqxx1 + S * dmusig1 * dd1 * ewu * (epsilon)/ssq1) *
    ewz/(wzdeno * (sigma_sq1) * sqrt(sigma_sq1))) - (2 * sigx5_2 + 2 * sigx5_1)^2/(2 *
    sigx6_2 + 2 * sigx6_1))/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (2 * (((ssqx1 * dmusig1 * depsi1 + S * (sigx17_1 * depsi1/(sigma_sq1) -
      S * ssqx1 * dmusig1 * dd1 * (epsilon)) * (epsilon)/(sigma_sq1))/wzdeno -
      0.5 * sigx18_1) * ewz/sqrt(sigma_sq1)) + 2 * ((ssqx2 * dmusig2 * depsi2 +
      (S * (sigx17_2 * depsi2/(sigma_sq2) - S * ssqx2 * dmusig2 * dd2 * (epsilon)) *
        (epsilon) - 0.5 * sigx18_2)/(sigma_sq2)) * prC/sqrt(sigma_sq2)) -
      (2 * sigx5_2 + 2 * sigx5_1) * (2 * sigx14_2 + 2 * sigx14_1)/(2 * sigx6_2 +
        2 * sigx6_1)) * ewu/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((S * (sigx17_1 * depsi1/(sigma_sq1) +
      S * sigx1_1 * dmusig1 * dd1 * ewu * (epsilon)/ssq1^2) * (epsilon)/(sigma_sq1) -
      sigx1_1 * dmusig1 * depsi1 * ewu/ssq1^2)/wzdeno - 0.5 * sigx18_1)/((2 *
      sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq1)) - (2 * sigx5_2 + 2 * sigx5_1) *
      sigx13_1 * sqrt(sigma_sq1)/((2 * sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq1))^2) *
      ewv1 * ewz), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((S * (sigx17_2 * depsi2/(sigma_sq2) +
      S * sigx1_2 * dmusig2 * dd2 * ewu * (epsilon)/ssq2^2) * (epsilon) - 0.5 *
      sigx18_2)/(sigma_sq2) - sigx1_2 * dmusig2 * depsi2 * ewu/ssq2^2)/sigx9_2 -
      (2 * sigx5_2 + 2 * sigx5_1) * sigx13_2 * sqrt(sigma_sq2)/sigx9_2^2) *
      prC * ewv2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * (sigx9_1 *
    sigx4_1/(sigma_sq1)) - ((2 * sigx5_2 + 2 * sigx5_1) * sigx25 + 2 * (prC *
    sigx4_2/(wzdeno * (sigma_sq2) * sqrt(sigma_sq2))))) * ewz/(2 * sigx6_2 +
    2 * sigx6_1), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (prC * (S * (0.5 * sigx21_2 - (0.5 * (S^2 *
      ssqx2 * depsi2 * (epsilon)^2/(sigma_sq2)^2) - (((0.5 * (ewu/(sigma_sq2)) +
      1 - 0.5 * (0.5 * prU2 + ewu/(sigma_sq2))) * prU2 * ewv2/sigmastar2 +
      (2 - 2 * (sigx7_2^2 * ewu * (sigma_sq2)/ssq2^2)) * sigmastar2)/ssq2^2 +
      S^2 * ssqx2^2 * ewu * (epsilon)^2/ssq2) * depsi2) * dmusig2) * (epsilon) -
      (0.5 * sigx12_2 + 0.5 * sigx22_2)/(sigma_sq2))/sqrt(sigma_sq2)) + 2 *
      (ewz * (S * (0.5 * sigx21_1 - (0.5 * (S^2 * ssqx1 * depsi1 * (epsilon)^2/(sigma_sq1)^2) -
        (((0.5 * (ewu/(sigma_sq1)) + 1 - 0.5 * (0.5 * prU1 + ewu/(sigma_sq1))) *
          prU1 * ewv1/sigmastar1 + (2 - 2 * (sigx7_1^2 * ewu * (sigma_sq1)/ssq1^2)) *
          sigmastar1)/ssq1^2 + S^2 * ssqx1^2 * ewu * (epsilon)^2/ssq1) *
          depsi1) * dmusig1) * (epsilon)/wzdeno - (0.5 * sigx23_1 + 0.5 *
        (sigx12_1/(sigma_sq1))))/sqrt(sigma_sq1)) - (2 * sigx14_2 + 2 * sigx14_1)^2/(2 *
      sigx6_2 + 2 * sigx6_1)) * ewu + 2 * sigx14_2 + 2 * sigx14_1) * ewu/(2 *
      sigx6_2 + 2 * sigx6_1), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((S *
    (((((0.5 * (prU1 * ewv1) - S^2 * sigx1_1 * ssqx1 * ewu * (epsilon)^2)/(sigma_sq1) +
      0.5 * ((ewu/(sigma_sq1) - 1) * ewv1/(sigma_sq1) + 1 - 0.5 * (prU1 * ewvsq1))) *
      depsi1/sigmastar1 + 0.5 * (S^2 * sigx1_1 * depsi1 * (epsilon)^2/(sigma_sq1)^2)) *
      ewu + sigx1_1 * (1 - 2 * (sigx7_1 * ewu * (sigma_sq1) * sigmastar1/ssq1^2)) *
      depsi1) * dmusig1/ssq1^2 + 0.5 * sigx21_1) * (epsilon)/wzdeno - 0.5 *
    sigx23_1)/((2 * sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq1)) - ((2 * sigx14_2 +
    2 * sigx14_1) * sqrt(sigma_sq1) + 0.5 * sigx27) * sigx13_1/((2 * sigx6_2 +
    2 * sigx6_1) * sqrt(sigma_sq1))^2) * ewu * ewv1 * ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((S * (((((0.5 * (prU2 * ewv2) - S^2 * sigx1_2 * ssqx2 * ewu * (epsilon)^2)/(sigma_sq2) +
    0.5 * ((ewu/(sigma_sq2) - 1) * ewv2/(sigma_sq2) + 1 - 0.5 * (prU2 * ewvsq2))) *
    depsi2/sigmastar2 + 0.5 * (S^2 * sigx1_2 * depsi2 * (epsilon)^2/(sigma_sq2)^2)) *
    ewu + sigx1_2 * (1 - 2 * (sigx7_2 * ewu * (sigma_sq2) * sigmastar2/ssq2^2)) *
    depsi2) * dmusig2/ssq2^2 + 0.5 * sigx21_2) * (epsilon) - 0.5 * (sigx22_2/(sigma_sq2)))/sigx9_2 -
    ((2 * sigx14_2 + 2 * sigx14_1) * sqrt(sigma_sq2) + 0.5 * sigx28) * sigx13_2/sigx9_2^2) *
    prC * ewu * ewv2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (2 * (sigx26_1 * pmusig1 - S * sigx9_1 * ssqx1 * dmusig1 * depsi1 * (epsilon)) -
      ((2 * sigx14_2 + 2 * sigx14_1) * sigx25 + 2 * (prC * (S * sigx11_2 *
        (epsilon)/wzdeno - 0.5 * sigx23_2)/sqrt(sigma_sq2)))) * ewu * ewz/(2 *
    sigx6_2 + 2 * sigx6_1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    2 * (((S * ((((0.5 * (ewv1/(sigma_sq1)) - 0.5 * (0.5 * ewvsq1 + ewv1/(sigma_sq1))) *
    ewvsq1 + S^2 * sigx1_1^2 * ewu * ewv1 * (epsilon)^2/(ssq1^2 * (sigma_sq1))) *
    depsi1 * ewu/sigmastar1 + ((0.5 * (S^2 * depsi1 * (epsilon)^2/(sigma_sq1)^2) -
    2 * (sigx1_1 * depsi1 * (sigma_sq1) * sigmastar1/ssq1^2)) * ewv1 + depsi1) *
    sigx1_1) * dmusig1 * ewu/ssq1^2 + S * (0.5 * (ewv1 * (S * (sigx1_1 * dmusig1 *
    ewu/ssq1^2 + 0.5 * pesig1) * (epsilon) - 2 * (pmusig1/(sigma_sq1)))) + 0.5 *
    pmusig1) * depsi1 * (epsilon)/(sigma_sq1)^2) * (epsilon)/wzdeno - (0.5 *
    (depsi1 * pmusig1) + 0.5 * (ewv1 * (S * sigx3_1 * (epsilon) - wzdeno^2 *
    depsi1 * pmusig1/wzsq1^2))) * wzdeno/wzsq1^2)/((2 * sigx6_2 + 2 * sigx6_1) *
    sqrt(sigma_sq1)) - (0.5 * sigx27 + 2 * (ewz * sigx13_1)) * ewv1 * sigx13_1/((2 *
    sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq1))^2) * ewv1 * ewz), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (prC * ewv1 * ewv2 * ewz * sigx13_1 * sigx13_2 * sqrt(sigma_sq2)/(sigx9_2^2 *
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx26_1 * pmusig1 + S * sigx1_1 * sigx9_1 *
      dmusig1 * depsi1 * ewu * (epsilon)/ssq1^2) - 2 * ((2 * sigx10_1 - 2 *
      sigx8_2) * ewz * sigx13_1/((2 * sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq1)))) *
      ewv1 * ewz/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * ((((0.5 * (ewv2/(sigma_sq2)) - 0.5 *
      (0.5 * ewvsq2 + ewv2/(sigma_sq2))) * ewvsq2 + S^2 * sigx1_2^2 * ewu *
      ewv2 * (epsilon)^2/(ssq2^2 * (sigma_sq2))) * depsi2 * ewu/sigmastar2 +
      ((0.5 * (S^2 * depsi2 * (epsilon)^2/(sigma_sq2)^2) - 2 * (sigx1_2 * depsi2 *
        (sigma_sq2) * sigmastar2/ssq2^2)) * ewv2 + depsi2) * sigx1_2) * dmusig2 *
      ewu/ssq2^2 + S * (0.5 * (ewv2 * (S * (sigx1_2 * dmusig2 * ewu/ssq2^2 +
      0.5 * pesig2) * (epsilon) - 2 * (pmusig2/(sigma_sq2)))) + 0.5 * pmusig2) *
      depsi2 * (epsilon)/(sigma_sq2)^2) * (epsilon) - (0.5 * (depsi2 * pmusig2) +
      0.5 * (ewv2 * (S * sigx3_2 * (epsilon) - depsi2 * pmusig2/(sigma_sq2))))/(sigma_sq2))/sigx9_2 -
      (0.5 * sigx28 + 2 * (prC * sigx13_2)) * ewv2 * sigx13_2/sigx9_2^2) *
      prC * ewv2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (prC * (2 * ((2 * sigx10_1 - 2 * sigx8_2) *
      sigx13_2/(2 * sigx6_2 + 2 * sigx6_1)) + 2 * (S * sigx3_2 * (epsilon)/wzdeno -
      0.5 * sigx23_2)) * ewv2 * ewz/sigx9_2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((2 * (prC * (1/(wzdeno^2 * sqrt(sigma_sq2)) +
      sqrt(sigma_sq2)/wzsq2^2) * depsi2 * pmusig2) - ((2 * sigx10_1 - 2 * sigx8_2)^2/(2 *
      sigx6_2 + 2 * sigx6_1) + 2 * ((2 - 2 * (wzdeno * (sigma_sq1) * ewz/wzsq1^2)) *
      depsi1 * pmusig1 * sqrt(sigma_sq1)/wzsq1^2))) * ewz + 2 * sigx10_1 -
      2 * sigx8_2) * ewz/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesscnsfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu <- 0
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/sqrt(sigma_sq2) + ewz1 * depsi1 * pmusig1/sqrt(sigma_sq1))
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- (ewz2 * sigx4_2/ssqq2 + ewz1 * sigx4_1/ssqq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (ewz2 * sigx11_2/sqrt(sigma_sq2) + sigx11_1 * ewz1/sqrt(sigma_sq1))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) * depsi1 +
      dmusig1 * (depsi1 * mupsi - depsi1 * musi1/ewv1) * ewu/starsq1) * ewz1/ssqq1 +
      (((mupsi^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) *
        depsi2 + dmusig2 * (depsi2 * mupsi - depsi2 * musi2/ewv2) * ewu/starsq2) *
        ewz2/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq1 - 2) * pmusig1 + dmusig1 *
      ewu * mupsi/starsq1) - 0.5 * pmusig1)/sigma_sq1 - sigx8_1 * dmusig1) *
      depsi1 * mupsi/sigma_sq1 - ((sigx7_1/starsq1^2 + 0.5/(sigma_sq1^2 * sigmastar1)) *
      ewu - (sigx8_1 * musi1/ewv1 + 1/sigmastar1)/sigma_sq1) * dmusig1 * depsi1) *
      ewz1/sqrt(sigma_sq1) + (((0.5 * ((mupsi^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu * mupsi/starsq2) - 0.5 * pmusig2)/sigma_sq2 - sigx8_2 *
      dmusig2) * depsi2 * mupsi/sigma_sq2 - ((sigx7_2/starsq2^2 + 0.5/(sigma_sq2^2 *
      sigmastar2)) * ewu - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/sigma_sq2) *
      dmusig2 * depsi2) * ewz2/sqrt(sigma_sq2) - sigx14 * sigx3/sigx2) * ewu/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq1 - 2) * pmusig1 +
      dmusig1 * ewu * mupsi/starsq1) - 0.5 * pmusig1)/sigma_sq1 + dmusig1 *
      sigx13_1) * depsi1 * mupsi/sigma_sq1 - ((sigx12_1/starsq1^2 + 0.5/(sigma_sq1^2 *
      sigmastar1)) * ewu + musi1 * sigx13_1/(sigma_sq1 * ewv1)) * dmusig1 *
      depsi1)/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) * ewz1 *
      ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu * mupsi/starsq2) - 0.5 * pmusig2)/sigma_sq2 + dmusig2 *
      sigx13_2) * depsi2 * mupsi/sigma_sq2 - ((sigx12_2/starsq2^2 + 0.5/(sigma_sq2^2 *
      sigmastar2)) * ewu + musi2 * sigx13_2/(sigma_sq2 * ewv2)) * dmusig2 *
      depsi2)/sigx17_2 - sigx3 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) * ewz2 *
      ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((sigx1_1/ssqq1 -
    sigx1_2/ssqq2)/(pi * sigx2 * ((Wz)^2 + 1)) - pi * sigx3 * ((Wz)^2 + 1) *
    sigx9/(pi * sigx2 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) - 0.5 * sigx10_1) * depsi1 * mupsi^2/sigma_sq1 -
      0.5 * sigx11_1) * ewu + 0.5 * (depsi1 * mupsi^2 * pmusig1/sigma_sq1))/sigma_sq1 -
      (((sigx8_1^2 * ewu/starsq1 + sigx7_1/starsq1^2) * musi1 + ((0.5 * (ewu/sigma_sq1) -
        0.5 * (0.5 * usq1 + ewu/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
        sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
          2 * (S * (epsilon))) * ewu)/starsq1^2 + S * (epsilon)/starsq1) *
        dmusig1 + (0.5 * pmusig1 - 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1) *
        ewu))/sigma_sq1) * depsi1) * ewz1/sqrt(sigma_sq1) + ((((0.5 * ((0.5 *
      (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) - 0.5 *
      sigx10_2) * depsi2 * mupsi^2/sigma_sq2 - 0.5 * sigx11_2) * ewu + 0.5 *
      (depsi2 * mupsi^2 * pmusig2/sigma_sq2))/sigma_sq2 - (((sigx8_2^2 * ewu/starsq2 +
      sigx7_2/starsq2^2) * musi2 + ((0.5 * (ewu/sigma_sq2) - 0.5 * (0.5 * usq2 +
      ewu/sigma_sq2)) * usq2 * ewv2 * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 *
      sigma_sq2 * musi2 * sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu)/starsq2^2 +
      S * (epsilon)/starsq2) * dmusig2 + (0.5 * pmusig2 - 0.5 * ((sigx8_2 *
      dmusig2 + pmusig2/sigma_sq2) * ewu))/sigma_sq2) * depsi2) * ewz2/sqrt(sigma_sq2) -
      sigx14^2 * ewu/sigx2) * ewu/sigx2 + 0.5 * (mu * (0.5 * (mu^2/(ewusr^3 *
      pmu)) - (0.5 * (ewusr * pmu) - 0.5 * (mu * dmu))/(ewusr * pmu)^2) * dmu)),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((sigx8_1 *
    musi1 * sigx13_1/starsq1 - ((0.5 * (usq1 * ewv1/sigma_sq1) + 0.5 * ((ewu/sigma_sq1 -
    1) * ewv1/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu *
    sigx7_1 - sigx12_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
    S * (epsilon)))/starsq1^2) * dmusig1 + 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1)/sigma_sq1)) *
    depsi1 + (0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 - sigx8_1 *
    dmusig1) + 0.5 * sigx15_1) * depsi1 * mupsi^2/sigma_sq1^2)/sigx17_1 - (sigx14 *
    sqrt(sigma_sq1) + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1/sigx17_1^2) *
    ewz1 * ewu * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((sigx8_2 * musi2 * sigx13_2/starsq2 - ((0.5 * (usq2 * ewv2/sigma_sq2) +
      0.5 * ((ewu/sigma_sq2 - 1) * ewv2/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) *
      musi2/sigmastar2 + mu * sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 *
      musi2 * sigmastar2/starsq2^2) + S * (epsilon)))/starsq2^2) * dmusig2 +
      0.5 * ((sigx8_2 * dmusig2 + pmusig2/sigma_sq2)/sigma_sq2)) * depsi2 +
      (0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 *
        dmusig2) + 0.5 * sigx15_2) * depsi2 * mupsi^2/sigma_sq2^2)/sigx17_2 -
      (sigx14 * sqrt(sigma_sq2) + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2/sigx17_2^2) *
    ewz2 * ewu * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((sigx11_1/sqrt(sigma_sq1) - sigx11_2/sqrt(sigma_sq2))/(pi * sigx2 * ((Wz)^2 +
      1)) - pi * sigx14 * ((Wz)^2 + 1) * sigx9/(pi * sigx2 * ((Wz)^2 + 1))^2) *
    ewu, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 + dmusig1 *
      sigx13_1) + 0.5 * sigx15_1) * ewv1 + 0.5 * pmusig1) * depsi1 * mupsi^2/sigma_sq1^2 +
      (dmusig1 * (mu/starsq1 - ((sigx12_1 * (2 * (mu) - 2 * (sigx12_1 * sigma_sq1 *
        musi1 * sigmastar1/starsq1^2)) * ewv1 + (0.5 * (ewv1/sigma_sq1) -
        0.5 * (0.5 * vsq1 + ewv1/sigma_sq1)) * vsq1 * ewu * musi1/sigmastar1)/starsq1^2 +
        (sigx12_1/starsq1^2 + ewv1 * sigx13_1^2/starsq1) * musi1)) - (0.5 *
        ((dmusig1 * sigx13_1 - pmusig1/sigma_sq1) * ewv1) + 0.5 * pmusig1)/sigma_sq1) *
        depsi1)/sigx17_1 - (sigx16_1 * ewz1 + 0.5 * (sigx2/sqrt(sigma_sq1))) *
      sigx16_1 * ewv1/sigx17_1^2) * ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx16_1 * sigx16_2 * ewz2 * ewz1 * ewv1 * ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx16_1 * (1/(pi * sigx2 * ((Wz)^2 + 1)) - pi *
      ((Wz)^2 + 1) * ewz1 * sigx9/(pi * sigx2 * ((Wz)^2 + 1))^2) * ewv1/sqrt(sigma_sq1),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) *
      pmusig2/sigma_sq2 + dmusig2 * sigx13_2) + 0.5 * sigx15_2) * ewv2 + 0.5 *
      pmusig2) * depsi2 * mupsi^2/sigma_sq2^2 + (dmusig2 * (mu/starsq2 - ((sigx12_2 *
      (2 * (mu) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
      ewv2 + (0.5 * (ewv2/sigma_sq2) - 0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) *
      vsq2 * ewu * musi2/sigmastar2)/starsq2^2 + (sigx12_2/starsq2^2 + ewv2 *
      sigx13_2^2/starsq2) * musi2)) - (0.5 * ((dmusig2 * sigx13_2 - pmusig2/sigma_sq2) *
      ewv2) + 0.5 * pmusig2)/sigma_sq2) * depsi2)/sigx17_2 - (sigx16_2 * ewz2 +
      0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2 * ewv2/sigx17_2^2) * ewz2 *
      ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_2 * (1/(pi * sigx2 * ((Wz)^2 + 1)) +
      pi * ((Wz)^2 + 1) * ewz2 * sigx9/(pi * sigx2 * ((Wz)^2 + 1))^2) * ewv2/sqrt(sigma_sq2)),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (pi * Wz * sigx2) + depsi1 * pmusig1/sqrt(sigma_sq1) -
      depsi2 * pmusig2/sqrt(sigma_sq2)) * sigx9/(pi * sigx2 * ((Wz)^2 + 1))^2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesscnsfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu <- 0
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/sqrt(sigma_sq2) + depsi1 * pmusig1 * pwZ/sqrt(sigma_sq1))
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- ((1 - pwZ) * sigx4_2/ssqq2 + sigx4_1 * pwZ/ssqq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (sigx11_1 * pwZ/sqrt(sigma_sq1) + sigx11_2 * (1 - pwZ)/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) * depsi1 +
      dmusig1 * (depsi1 * mupsi - depsi1 * musi1/ewv1) * ewu/starsq1) * pwZ/ssqq1 +
      (((mupsi^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) *
        depsi2 + dmusig2 * (depsi2 * mupsi - depsi2 * musi2/ewv2) * ewu/starsq2) *
        (1 - pwZ)/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq1 - 2) * pmusig1 + dmusig1 *
      ewu * mupsi/starsq1) - 0.5 * pmusig1)/sigma_sq1 - sigx8_1 * dmusig1) *
      depsi1 * mupsi/sigma_sq1 - ((sigx7_1/starsq1^2 + 0.5/(sigma_sq1^2 * sigmastar1)) *
      ewu - (sigx8_1 * musi1/ewv1 + 1/sigmastar1)/sigma_sq1) * dmusig1 * depsi1) *
      pwZ/sqrt(sigma_sq1) + (((0.5 * ((mupsi^2/sigma_sq2 - 2) * pmusig2 + dmusig2 *
      ewu * mupsi/starsq2) - 0.5 * pmusig2)/sigma_sq2 - sigx8_2 * dmusig2) *
      depsi2 * mupsi/sigma_sq2 - ((sigx7_2/starsq2^2 + 0.5/(sigma_sq2^2 * sigmastar2)) *
      ewu - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/sigma_sq2) * dmusig2 * depsi2) *
      (1 - pwZ)/sqrt(sigma_sq2) - sigx14 * sigx3/sigx2) * ewu/sigx2, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq1 - 2) * pmusig1 +
      dmusig1 * ewu * mupsi/starsq1) - 0.5 * pmusig1)/sigma_sq1 + dmusig1 *
      sigx13_1) * depsi1 * mupsi/sigma_sq1 - ((sigx12_1/starsq1^2 + 0.5/(sigma_sq1^2 *
      sigmastar1)) * ewu + musi1 * sigx13_1/(sigma_sq1 * ewv1)) * dmusig1 *
      depsi1)/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) * pwZ *
      ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu * mupsi/starsq2) - 0.5 * pmusig2)/sigma_sq2 + dmusig2 *
      sigx13_2) * depsi2 * mupsi/sigma_sq2 - ((sigx12_2/starsq2^2 + 0.5/(sigma_sq2^2 *
      sigmastar2)) * ewu + musi2 * sigx13_2/(sigma_sq2 * ewv2)) * dmusig2 *
      depsi2)/sigx17_2 - sigx3 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) * (1 -
      pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx1_1/ssqq1 -
    (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) * dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) - 0.5 * sigx10_1) * depsi1 * mupsi^2/sigma_sq1 -
      0.5 * sigx11_1) * ewu + 0.5 * (depsi1 * mupsi^2 * pmusig1/sigma_sq1))/sigma_sq1 -
      (((sigx8_1^2 * ewu/starsq1 + sigx7_1/starsq1^2) * musi1 + ((0.5 * (ewu/sigma_sq1) -
        0.5 * (0.5 * usq1 + ewu/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
        sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
          2 * (S * (epsilon))) * ewu)/starsq1^2 + S * (epsilon)/starsq1) *
        dmusig1 + (0.5 * pmusig1 - 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1) *
        ewu))/sigma_sq1) * depsi1) * pwZ/sqrt(sigma_sq1) + ((((0.5 * ((0.5 *
      (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) - 0.5 *
      sigx10_2) * depsi2 * mupsi^2/sigma_sq2 - 0.5 * sigx11_2) * ewu + 0.5 *
      (depsi2 * mupsi^2 * pmusig2/sigma_sq2))/sigma_sq2 - (((sigx8_2^2 * ewu/starsq2 +
      sigx7_2/starsq2^2) * musi2 + ((0.5 * (ewu/sigma_sq2) - 0.5 * (0.5 * usq2 +
      ewu/sigma_sq2)) * usq2 * ewv2 * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 *
      sigma_sq2 * musi2 * sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu)/starsq2^2 +
      S * (epsilon)/starsq2) * dmusig2 + (0.5 * pmusig2 - 0.5 * ((sigx8_2 *
      dmusig2 + pmusig2/sigma_sq2) * ewu))/sigma_sq2) * depsi2) * (1 - pwZ)/sqrt(sigma_sq2) -
      sigx14^2 * ewu/sigx2) * ewu/sigx2 + 0.5 * (mu * (0.5 * (mu^2/(ewusr^3 *
      pmu)) - (0.5 * (ewusr * pmu) - 0.5 * (mu * dmu))/(ewusr * pmu)^2) * dmu)),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((sigx8_1 *
    musi1 * sigx13_1/starsq1 - ((0.5 * (usq1 * ewv1/sigma_sq1) + 0.5 * ((ewu/sigma_sq1 -
    1) * ewv1/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu *
    sigx7_1 - sigx12_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
    S * (epsilon)))/starsq1^2) * dmusig1 + 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1)/sigma_sq1)) *
    depsi1 + (0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 - sigx8_1 *
    dmusig1) + 0.5 * sigx15_1) * depsi1 * mupsi^2/sigma_sq1^2)/sigx17_1 - (sigx14 *
    sqrt(sigma_sq1) + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1/sigx17_1^2) *
    pwZ * ewu * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((sigx8_2 * musi2 * sigx13_2/starsq2 - ((0.5 * (usq2 * ewv2/sigma_sq2) +
      0.5 * ((ewu/sigma_sq2 - 1) * ewv2/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) *
      musi2/sigmastar2 + mu * sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 *
      musi2 * sigmastar2/starsq2^2) + S * (epsilon)))/starsq2^2) * dmusig2 +
      0.5 * ((sigx8_2 * dmusig2 + pmusig2/sigma_sq2)/sigma_sq2)) * depsi2 +
      (0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 *
        dmusig2) + 0.5 * sigx15_2) * depsi2 * mupsi^2/sigma_sq2^2)/sigx17_2 -
      (sigx14 * sqrt(sigma_sq2) + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2/sigx17_2^2) *
    (1 - pwZ) * ewu * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (sigx11_1/sqrt(sigma_sq1) - (sigx14 * sigx9/sigx2 + sigx11_2/sqrt(sigma_sq2))) *
    dwZ * ewu/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 + dmusig1 *
      sigx13_1) + 0.5 * sigx15_1) * ewv1 + 0.5 * pmusig1) * depsi1 * mupsi^2/sigma_sq1^2 +
      (dmusig1 * (mu/starsq1 - ((sigx12_1 * (2 * (mu) - 2 * (sigx12_1 * sigma_sq1 *
        musi1 * sigmastar1/starsq1^2)) * ewv1 + (0.5 * (ewv1/sigma_sq1) -
        0.5 * (0.5 * vsq1 + ewv1/sigma_sq1)) * vsq1 * ewu * musi1/sigmastar1)/starsq1^2 +
        (sigx12_1/starsq1^2 + ewv1 * sigx13_1^2/starsq1) * musi1)) - (0.5 *
        ((dmusig1 * sigx13_1 - pmusig1/sigma_sq1) * ewv1) + 0.5 * pmusig1)/sigma_sq1) *
        depsi1)/sigx17_1 - (sigx16_1 * pwZ + 0.5 * (sigx2/sqrt(sigma_sq1))) *
      sigx16_1 * ewv1/sigx17_1^2) * pwZ * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx16_1 * sigx16_2 * (1 - pwZ) * pwZ * ewv1 * ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx16_1 * (1 - sigx9 * pwZ/sigx2) * dwZ * ewv1/(sigx2 *
      sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) *
      pmusig2/sigma_sq2 + dmusig2 * sigx13_2) + 0.5 * sigx15_2) * ewv2 + 0.5 *
      pmusig2) * depsi2 * mupsi^2/sigma_sq2^2 + (dmusig2 * (mu/starsq2 - ((sigx12_2 *
      (2 * (mu) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
      ewv2 + (0.5 * (ewv2/sigma_sq2) - 0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) *
      vsq2 * ewu * musi2/sigmastar2)/starsq2^2 + (sigx12_2/starsq2^2 + ewv2 *
      sigx13_2^2/starsq2) * musi2)) - (0.5 * ((dmusig2 * sigx13_2 - pmusig2/sigma_sq2) *
      ewv2) + 0.5 * pmusig2)/sigma_sq2) * depsi2)/sigx17_2 - (sigx16_2 * (1 -
      pwZ) + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2 * ewv2/sigx17_2^2) *
      (1 - pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((1 - pwZ) * sigx9/sigx2 + 1) * sigx16_2 *
      dwZ * ewv2/(sigx2 * sqrt(sigma_sq2))), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx9 * dwZ/sigx2 + Wz) * sigx9 * dwZ/sigx2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesscnsfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu <- 0
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/sqrt(sigma_sq1) + depsi2 * prZ * pmusig2/sqrt(sigma_sq2))
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- ((1 - prZ) * sigx4_1/ssqq1 + sigx4_2 * prZ/ssqq2)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (sigx11_1 * (1 - prZ)/sqrt(sigma_sq1) + sigx11_2 * prZ/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) * depsi1 +
      dmusig1 * (depsi1 * mupsi - depsi1 * musi1/ewv1) * ewu/starsq1) * (1 -
      prZ)/ssqq1 + (((mupsi^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) *
      depsi2 + dmusig2 * (depsi2 * mupsi - depsi2 * musi2/ewv2) * ewu/starsq2) *
      prZ/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq1 - 2) * pmusig1 + dmusig1 *
      ewu * mupsi/starsq1) - 0.5 * pmusig1)/sigma_sq1 - sigx8_1 * dmusig1) *
      depsi1 * mupsi/sigma_sq1 - ((sigx7_1/starsq1^2 + 0.5/(sigma_sq1^2 * sigmastar1)) *
      ewu - (sigx8_1 * musi1/ewv1 + 1/sigmastar1)/sigma_sq1) * dmusig1 * depsi1) *
      (1 - prZ)/sqrt(sigma_sq1) + (((0.5 * ((mupsi^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu * mupsi/starsq2) - 0.5 * pmusig2)/sigma_sq2 - sigx8_2 *
      dmusig2) * depsi2 * mupsi/sigma_sq2 - ((sigx7_2/starsq2^2 + 0.5/(sigma_sq2^2 *
      sigmastar2)) * ewu - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/sigma_sq2) *
      dmusig2 * depsi2) * prZ/sqrt(sigma_sq2) - sigx14 * sigx3/sigx2) * ewu/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq1 - 2) * pmusig1 +
      dmusig1 * ewu * mupsi/starsq1) - 0.5 * pmusig1)/sigma_sq1 + dmusig1 *
      sigx13_1) * depsi1 * mupsi/sigma_sq1 - ((sigx12_1/starsq1^2 + 0.5/(sigma_sq1^2 *
      sigmastar1)) * ewu + musi1 * sigx13_1/(sigma_sq1 * ewv1)) * dmusig1 *
      depsi1)/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) * (1 -
      prZ) * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu * mupsi/starsq2) - 0.5 * pmusig2)/sigma_sq2 + dmusig2 *
      sigx13_2) * depsi2 * mupsi/sigma_sq2 - ((sigx12_2/starsq2^2 + 0.5/(sigma_sq2^2 *
      sigmastar2)) * ewu + musi2 * sigx13_2/(sigma_sq2 * ewv2)) * dmusig2 *
      depsi2)/sigx17_2 - sigx3 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) * prZ *
      ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx1_1/ssqq1 -
    (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) * prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) - 0.5 * sigx10_1) * depsi1 * mupsi^2/sigma_sq1 -
      0.5 * sigx11_1) * ewu + 0.5 * (depsi1 * mupsi^2 * pmusig1/sigma_sq1))/sigma_sq1 -
      (((sigx8_1^2 * ewu/starsq1 + sigx7_1/starsq1^2) * musi1 + ((0.5 * (ewu/sigma_sq1) -
        0.5 * (0.5 * usq1 + ewu/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
        sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
          2 * (S * (epsilon))) * ewu)/starsq1^2 + S * (epsilon)/starsq1) *
        dmusig1 + (0.5 * pmusig1 - 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1) *
        ewu))/sigma_sq1) * depsi1) * (1 - prZ)/sqrt(sigma_sq1) + ((((0.5 *
      ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) -
      0.5 * sigx10_2) * depsi2 * mupsi^2/sigma_sq2 - 0.5 * sigx11_2) * ewu +
      0.5 * (depsi2 * mupsi^2 * pmusig2/sigma_sq2))/sigma_sq2 - (((sigx8_2^2 *
      ewu/starsq2 + sigx7_2/starsq2^2) * musi2 + ((0.5 * (ewu/sigma_sq2) -
      0.5 * (0.5 * usq2 + ewu/sigma_sq2)) * usq2 * ewv2 * musi2/sigmastar2 -
      sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
        2 * (S * (epsilon))) * ewu)/starsq2^2 + S * (epsilon)/starsq2) *
      dmusig2 + (0.5 * pmusig2 - 0.5 * ((sigx8_2 * dmusig2 + pmusig2/sigma_sq2) *
      ewu))/sigma_sq2) * depsi2) * prZ/sqrt(sigma_sq2) - sigx14^2 * ewu/sigx2) *
      ewu/sigx2 + 0.5 * (mu * (0.5 * (mu^2/(ewusr^3 * pmu)) - (0.5 * (ewusr *
      pmu) - 0.5 * (mu * dmu))/(ewusr * pmu)^2) * dmu)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((sigx8_1 *
    musi1 * sigx13_1/starsq1 - ((0.5 * (usq1 * ewv1/sigma_sq1) + 0.5 * ((ewu/sigma_sq1 -
    1) * ewv1/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu *
    sigx7_1 - sigx12_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
    S * (epsilon)))/starsq1^2) * dmusig1 + 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1)/sigma_sq1)) *
    depsi1 + (0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 - sigx8_1 *
    dmusig1) + 0.5 * sigx15_1) * depsi1 * mupsi^2/sigma_sq1^2)/sigx17_1 - (sigx14 *
    sqrt(sigma_sq1) + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1/sigx17_1^2) *
    (1 - prZ) * ewu * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((sigx8_2 * musi2 * sigx13_2/starsq2 - ((0.5 * (usq2 * ewv2/sigma_sq2) +
      0.5 * ((ewu/sigma_sq2 - 1) * ewv2/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) *
      musi2/sigmastar2 + mu * sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 *
      musi2 * sigmastar2/starsq2^2) + S * (epsilon)))/starsq2^2) * dmusig2 +
      0.5 * ((sigx8_2 * dmusig2 + pmusig2/sigma_sq2)/sigma_sq2)) * depsi2 +
      (0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 *
        dmusig2) + 0.5 * sigx15_2) * depsi2 * mupsi^2/sigma_sq2^2)/sigx17_2 -
      (sigx14 * sqrt(sigma_sq2) + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2/sigx17_2^2) *
    prZ * ewu * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (sigx11_1/sqrt(sigma_sq1) - (sigx14 * sigx9/sigx2 + sigx11_2/sqrt(sigma_sq2))) *
    prZ * ewu * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 + dmusig1 *
      sigx13_1) + 0.5 * sigx15_1) * ewv1 + 0.5 * pmusig1) * depsi1 * mupsi^2/sigma_sq1^2 +
      (dmusig1 * (mu/starsq1 - ((sigx12_1 * (2 * (mu) - 2 * (sigx12_1 * sigma_sq1 *
        musi1 * sigmastar1/starsq1^2)) * ewv1 + (0.5 * (ewv1/sigma_sq1) -
        0.5 * (0.5 * vsq1 + ewv1/sigma_sq1)) * vsq1 * ewu * musi1/sigmastar1)/starsq1^2 +
        (sigx12_1/starsq1^2 + ewv1 * sigx13_1^2/starsq1) * musi1)) - (0.5 *
        ((dmusig1 * sigx13_1 - pmusig1/sigma_sq1) * ewv1) + 0.5 * pmusig1)/sigma_sq1) *
        depsi1)/sigx17_1 - (sigx16_1 * (1 - prZ) + 0.5 * (sigx2/sqrt(sigma_sq1))) *
      sigx16_1 * ewv1/sigx17_1^2) * (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx16_1 * sigx16_2 * prZ * (1 - prZ) * ewv1 * ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx16_1 * (1 - (1 - prZ) * sigx9/sigx2) * prZ *
      ewv1 * ewz/(sigx2 * sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) *
      pmusig2/sigma_sq2 + dmusig2 * sigx13_2) + 0.5 * sigx15_2) * ewv2 + 0.5 *
      pmusig2) * depsi2 * mupsi^2/sigma_sq2^2 + (dmusig2 * (mu/starsq2 - ((sigx12_2 *
      (2 * (mu) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
      ewv2 + (0.5 * (ewv2/sigma_sq2) - 0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) *
      vsq2 * ewu * musi2/sigmastar2)/starsq2^2 + (sigx12_2/starsq2^2 + ewv2 *
      sigx13_2^2/starsq2) * musi2)) - (0.5 * ((dmusig2 * sigx13_2 - pmusig2/sigma_sq2) *
      ewv2) + 0.5 * pmusig2)/sigma_sq2) * depsi2)/sigx17_2 - (sigx16_2 * prZ +
      0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2 * ewv2/sigx17_2^2) * prZ *
      ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_2 * (sigx9 * prZ/sigx2 + 1) * prZ *
      ewv2 * ewz/(sigx2 * sqrt(sigma_sq2))), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (sigx9 * prZ/sigx2 + 1) * ewz) * sigx9 *
      prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# different sigma_u

## logit specification class membership
chessmcesfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu1 * (epsilon)/ssq1)
  musig2 <- (S * ewu2 * (epsilon)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsi1 <- S * (epsilon)/sqrt(sigma_sq1)
  epsi2 <- S * (epsilon)/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  sxq <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 * pmusig1 * (epsilon))
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + S * depsi2 * pmusig2 * (epsilon))
  wzsq1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  wzxsq1 <- (wzdeno * sqrt(sigma_sq1))
  sigx2 <- (2 * (prC * sigx1_2/sxq) + 2 * (sigx1_1 * ewz/wzsq1))
  sigx3_1 <- (depsi1 * ewz * pmusig1/wzxsq1)
  sigx3_2 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx4_1 <- (S * depsi1 * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx4_2 <- (S * depsi2 * pmusig2 * (epsilon)/(sigma_sq2)^2)
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  usq1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  usq2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
  sigx5_1 <- (1/ssq1 - usq1 * ewu1/ssq1^2)
  sigx5_2 <- (1/ssq2 - usq2 * ewu2/ssq2^2)
  sigx6_1 <- (0.5 * sigx4_1 - sigx5_1 * dmusig1 * depsi1)
  sigx6_2 <- (0.5 * sigx4_2 - sigx5_2 * dmusig2 * depsi2)
  s7sq2 <- (depsi2 * pmusig2/(sigma_sq2))
  sigx11_1 <- (wzdeno * depsi1 * pmusig1/wzxsq1^2)
  sigx11_2 <- (prC * depsi2 * pmusig2/(wzdeno * sqrt(sigma_sq2)))
  sigx7_1 <- (S * sigx6_1 * (epsilon)/wzdeno - 0.5 * sigx11_1)
  sigx7_2 <- (S * sigx6_2 * (epsilon) - 0.5 * s7sq2)
  sigx8_1 <- ((2 * sigx3_2 + 2 * sigx3_1) * sqrt(sigma_sq1))
  sigx8_2 <- ((2 * sigx3_2 + 2 * sigx3_1) * sqrt(sigma_sq2))
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  sigx9_1 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  s9sq2 <- (0.5 * (prV2 * ewu2/sigmastar2) + sigmastar2)
  sigx9_2 <- (s9sq2 * dmusig2 * depsi2 * ewu2/ssq2^2 + 0.5 * sigx4_2)
  sigx10_1 <- (sigx9_1 * dmusig1 * depsi1 * ewu1/ssq1^2 + 0.5 * sigx4_1)
  sigx10_2 <- (S * sigx9_2 * (epsilon) - 0.5 * s7sq2)
  s12sq1 <- (1/wzxsq1 - ewz * sqrt(sigma_sq1)/wzxsq1^2)
  sigx12 <- (2 * (s12sq1 * depsi1 * pmusig1) - 2 * sigx11_2)
  sigx18_1 <- (S * sigx10_1 * (epsilon)/wzdeno - 0.5 * sigx11_1)
  sigx13_1 <- (dmusig1 * ewu1/sigmastar1 + S * pmusig1 * (epsilon))
  sigx13_2 <- (dmusig2 * ewu2/sigmastar2 + S * pmusig2 * (epsilon))
  sigx14_1 <- (S * sigx13_1 * (epsilon)/(sigma_sq1) - pmusig1)
  sigx14_2 <- (S * sigx13_2 * (epsilon)/(sigma_sq2) - pmusig2)
  sigx15_1 <- (0.5 * sigx14_1 - 0.5 * pmusig1) * depsi1/(sigma_sq1)
  sigx15_2 <- (0.5 * sigx14_2 - 0.5 * pmusig2) * depsi2/(sigma_sq2)
  sigx16_1 <- (depsi1 * ewu1/ewv1 + depsi1)
  sigx16_2 <- (depsi2 * ewu2/ewv2 + depsi2)
  sigx17_1 <- (wzdeno * sigx1_1/(wzxsq1^2 * (sigma_sq1)))
  sigx17_2 <- (sigx1_2/(sigma_sq2))
  sigx19_1 <- (S^2 * s12sq1 * depsi1 * (epsilon)^2/(sigma_sq1)^2)
  sigx20_1 <- ((0.5/sqrt(sigma_sq1) - wzdeno^2 * sqrt(sigma_sq1)/wzxsq1^2) * ewz +
    0.5 * (wzdeno/sqrt(sigma_sq1)))
  sigx21_1 <- (0.5 * sigx19_1 - sigx20_1 * depsi1/wzxsq1^2) * pmusig1
  sigx22_1 <- (S * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx22_2 <- (S * pmusig2 * (epsilon)/(sigma_sq2)^2)
  sigx23_1 <- (S * (0.5 * sigx22_1 - sigx5_1 * dmusig1) * (epsilon) - 2 * (pmusig1/(sigma_sq1)))
  sigx23_2 <- (S * (0.5 * sigx22_2 - sigx5_2 * dmusig2) * (epsilon) - 2 * (pmusig2/(sigma_sq2)))
  sigx24_1 <- (S * depsi1 * sigx23_1 * (epsilon)/(sigma_sq1)^2)
  sigx24_2 <- (S * depsi2 * sigx23_2 * (epsilon)/(sigma_sq2)^2)
  sigx25_1 <- (wzdeno * (S * sigx6_1 * (epsilon) - wzdeno^2 * depsi1 * pmusig1/wzxsq1^2)/wzxsq1^2)
  sigx25_2 <- ((S * sigx6_2 * (epsilon) - depsi2 * pmusig2/(sigma_sq2))/(sigma_sq2))
  sigx26_1 <- ((2 * sigx3_2 + 2 * sigx3_1)/sqrt(sigma_sq1))
  sigx26_2 <- ((2 * sigx3_2 + 2 * sigx3_1)/sqrt(sigma_sq2))
  sigx27_1 <- (sigx8_2^2 * sqrt(sigma_sq1))
  sigx27_2 <- (wzdeno * depsi2 * pmusig2/(wzdeno * sqrt(sigma_sq2))^2)
  sigx28_1 <- (0.5 * sigx26_1 + 2 * (ewz * sigx7_1))
  sigx28_2 <- (0.5 * sigx26_2 + 2 * (prC * sigx7_2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (2 * (prC * (depsi2 * sigx14_2 + S * dmusig2 * sigx16_2 * ewu2 * (epsilon)/ssq2)/sxq) +
    2 * ((depsi1 * sigx14_1 + S * dmusig1 * sigx16_1 * ewu1 * (epsilon)/ssq1) *
      ewz/wzsq1) - sigx2^2/(2 * sigx3_2 + 2 * sigx3_1))/(2 * sigx3_2 + 2 *
    sigx3_1), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((sigx5_1 * dmusig1 * depsi1 + S * (sigx15_1 -
      S * sigx5_1 * dmusig1 * sigx16_1 * (epsilon)) * (epsilon)/(sigma_sq1))/wzdeno -
      0.5 * sigx17_1)/sigx8_1 - sigx2 * sigx7_1 * sqrt(sigma_sq1)/sigx8_1^2) *
      ewu1 * ewz), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * ((sigx5_2 * dmusig2 * depsi2 + (S *
      (sigx15_2 - S * sigx5_2 * dmusig2 * sigx16_2 * (epsilon)) * (epsilon) -
      0.5 * sigx17_2)/(sigma_sq2))/sigx8_2 - sigx2 * sigx7_2 * sqrt(sigma_sq2)/sigx8_2^2) *
      prC * ewu2), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((S * (sigx15_1 + S * sigx9_1 * dmusig1 *
      sigx16_1 * ewu1 * (epsilon)/ssq1^2) * (epsilon)/(sigma_sq1) - sigx9_1 *
      dmusig1 * depsi1 * ewu1/ssq1^2)/wzdeno - 0.5 * sigx17_1)/sigx8_1 - sigx2 *
      sigx18_1 * sqrt(sigma_sq1)/sigx8_1^2) * ewv1 * ewz), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * 2 * (S * (((S *
    (sigx15_2 + S * s9sq2 * dmusig2 * sigx16_2 * ewu2 * (epsilon)/ssq2^2) * (epsilon) -
    0.5 * sigx17_2)/(sigma_sq2) - s9sq2 * dmusig2 * depsi2 * ewu2/ssq2^2)/sigx8_2 -
    sigx2 * sigx10_2 * sqrt(sigma_sq2)/sigx8_2^2) * prC * ewv2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (2 * (s12sq1 * sigx1_1/(sigma_sq1)) - (sigx2 * sigx12/(2 * sigx3_2 +
    2 * sigx3_1) + 2 * (prC * sigx1_2/(wzdeno * (sigma_sq2) * sqrt(sigma_sq2))))) *
    ewz/(2 * sigx3_2 + 2 * sigx3_1), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((ewu1 * (S * (0.5 * sigx24_1 - (0.5 * (S^2 *
      sigx5_1 * depsi1 * (epsilon)^2/(sigma_sq1)^2) - (((0.5 * (ewu1/(sigma_sq1)) +
      1 - 0.5 * (0.5 * prU1 + ewu1/(sigma_sq1))) * prU1 * ewv1/sigmastar1 +
      (2 - 2 * (usq1^2 * ewu1 * (sigma_sq1)/ssq1^2)) * sigmastar1)/ssq1^2 +
      S^2 * sigx5_1^2 * ewu1 * (epsilon)^2/ssq1) * depsi1) * dmusig1) * (epsilon)/wzdeno -
      0.5 * sigx25_1) + S * sigx6_1 * (epsilon)/wzdeno - 0.5 * sigx11_1)/sigx8_1 -
      sigx28_1 * ewu1 * sigx7_1/sigx8_1^2) * ewu1 * ewz), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (prC * ewu1 * ewu2 * ewz * sigx7_1 * sigx7_2 *
      sqrt(sigma_sq2)/sigx27_1)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((S * (((((0.5 * (prU1 * ewv1) - S^2 * sigx9_1 * sigx5_1 * ewu1 * (epsilon)^2)/(sigma_sq1) +
    0.5 * ((ewu1/(sigma_sq1) - 1) * ewv1/(sigma_sq1) + 1 - 0.5 * (prU1 * prV1))) *
    depsi1/sigmastar1 + 0.5 * (S^2 * sigx9_1 * depsi1 * (epsilon)^2/(sigma_sq1)^2)) *
    ewu1 + sigx9_1 * (1 - 2 * (usq1 * ewu1 * (sigma_sq1) * sigmastar1/ssq1^2)) *
    depsi1) * dmusig1/ssq1^2 + 0.5 * sigx24_1) * (epsilon)/wzdeno - 0.5 * sigx25_1)/sigx8_1 -
    sigx28_1 * sigx18_1/sigx8_1^2) * ewu1 * ewv1 * ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (4 * (prC * ewu1 * ewv2 * ewz * sigx10_2 * sigx7_1 * sqrt(sigma_sq2)/sigx27_1)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * (2 * (sigx21_1 - S * s12sq1 * sigx5_1 * dmusig1 * depsi1 *
      (epsilon)) - 2 * (sigx12 * ewz * sigx7_1/sigx8_1)) * ewu1 * ewz/(2 *
      sigx3_2 + 2 * sigx3_1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((ewu2 *
    (S * (0.5 * sigx24_2 - (0.5 * (S^2 * sigx5_2 * depsi2 * (epsilon)^2/(sigma_sq2)^2) -
      (((0.5 * (ewu2/(sigma_sq2)) + 1 - 0.5 * (0.5 * prU2 + ewu2/(sigma_sq2))) *
        prU2 * ewv2/sigmastar2 + (2 - 2 * (usq2^2 * ewu2 * (sigma_sq2)/ssq2^2)) *
        sigmastar2)/ssq2^2 + S^2 * sigx5_2^2 * ewu2 * (epsilon)^2/ssq2) *
        depsi2) * dmusig2) * (epsilon) - 0.5 * sigx25_2) + S * sigx6_2 *
    (epsilon) - 0.5 * s7sq2)/sigx8_2 - sigx28_2 * ewu2 * sigx7_2/sigx8_2^2) *
    prC * ewu2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (4 * (prC * ewu2 * ewv1 * ewz * sigx18_1 * sigx7_2 * sqrt(sigma_sq1)/(sigx8_1^2 *
      sqrt(sigma_sq2)))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * 2 * (((S * (((((0.5 * (prU2 * ewv2) - S^2 * s9sq2 * sigx5_2 *
      ewu2 * (epsilon)^2)/(sigma_sq2) + 0.5 * ((ewu2/(sigma_sq2) - 1) * ewv2/(sigma_sq2) +
      1 - 0.5 * (prU2 * prV2))) * depsi2/sigmastar2 + 0.5 * (S^2 * s9sq2 *
      depsi2 * (epsilon)^2/(sigma_sq2)^2)) * ewu2 + s9sq2 * (1 - 2 * (usq2 *
      ewu2 * (sigma_sq2) * sigmastar2/ssq2^2)) * depsi2) * dmusig2/ssq2^2 +
      0.5 * sigx24_2) * (epsilon) - 0.5 * sigx25_2)/sigx8_2 - sigx28_2 * sigx10_2/sigx8_2^2) *
      prC * ewu2 * ewv2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 2 *
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (prC * (2 * (sigx12 * sigx7_2/(2 * sigx3_2 +
      2 * sigx3_1)) + 2 * (S * sigx6_2 * (epsilon)/wzdeno - 0.5 * sigx27_2)) *
      ewu2 * ewz/sigx8_2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * 2 * (((S * ((((0.5 * (ewv1/(sigma_sq1)) - 0.5 * (0.5 * prV1 +
      ewv1/(sigma_sq1))) * prV1 + S^2 * sigx9_1^2 * ewu1 * ewv1 * (epsilon)^2/(ssq1^2 *
      (sigma_sq1))) * depsi1 * ewu1/sigmastar1 + ((0.5 * (S^2 * depsi1 * (epsilon)^2/(sigma_sq1)^2) -
      2 * (sigx9_1 * depsi1 * (sigma_sq1) * sigmastar1/ssq1^2)) * ewv1 + depsi1) *
      sigx9_1) * dmusig1 * ewu1/ssq1^2 + S * (0.5 * (ewv1 * (S * (sigx9_1 *
      dmusig1 * ewu1/ssq1^2 + 0.5 * sigx22_1) * (epsilon) - 2 * (pmusig1/(sigma_sq1)))) +
      0.5 * pmusig1) * depsi1 * (epsilon)/(sigma_sq1)^2) * (epsilon)/wzdeno -
      (0.5 * (depsi1 * pmusig1) + 0.5 * (ewv1 * (S * sigx10_1 * (epsilon) -
        wzdeno^2 * depsi1 * pmusig1/wzxsq1^2))) * wzdeno/wzxsq1^2)/sigx8_1 -
      (0.5 * sigx26_1 + 2 * (ewz * sigx18_1)) * ewv1 * sigx18_1/sigx8_1^2) *
      ewv1 * ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (prC * ewv1 * ewv2 * ewz * sigx18_1 * sigx10_2 *
      sqrt(sigma_sq2)/sigx27_1)), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx21_1 + S * sigx9_1 * s12sq1 * dmusig1 *
      depsi1 * ewu1 * (epsilon)/ssq1^2) - 2 * (sigx12 * ewz * sigx18_1/sigx8_1)) *
      ewv1 * ewz/(2 * sigx3_2 + 2 * sigx3_1), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * ((((0.5 * (ewv2/(sigma_sq2)) - 0.5 *
      (0.5 * prV2 + ewv2/(sigma_sq2))) * prV2 + S^2 * s9sq2^2 * ewu2 * ewv2 *
      (epsilon)^2/(ssq2^2 * (sigma_sq2))) * depsi2 * ewu2/sigmastar2 + ((0.5 *
      (S^2 * depsi2 * (epsilon)^2/(sigma_sq2)^2) - 2 * (s9sq2 * depsi2 * (sigma_sq2) *
      sigmastar2/ssq2^2)) * ewv2 + depsi2) * s9sq2) * dmusig2 * ewu2/ssq2^2 +
      S * (0.5 * (ewv2 * (S * (s9sq2 * dmusig2 * ewu2/ssq2^2 + 0.5 * sigx22_2) *
        (epsilon) - 2 * (pmusig2/(sigma_sq2)))) + 0.5 * pmusig2) * depsi2 *
        (epsilon)/(sigma_sq2)^2) * (epsilon) - (0.5 * (depsi2 * pmusig2) +
      0.5 * (ewv2 * (S * sigx9_2 * (epsilon) - depsi2 * pmusig2/(sigma_sq2))))/(sigma_sq2))/sigx8_2 -
      (0.5 * sigx26_2 + 2 * (prC * sigx10_2)) * ewv2 * sigx10_2/sigx8_2^2) *
      prC * ewv2), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (prC *
    (2 * (sigx12 * sigx10_2/(2 * sigx3_2 + 2 * sigx3_1)) + 2 * (S * sigx9_2 *
      (epsilon)/wzdeno - 0.5 * sigx27_2)) * ewv2 * ewz/sigx8_2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar * ((2 *
    (prC * (1/(wzdeno^2 * sqrt(sigma_sq2)) + sqrt(sigma_sq2)/(wzdeno * sqrt(sigma_sq2))^2) *
      depsi2 * pmusig2) - (sigx12^2/(2 * sigx3_2 + 2 * sigx3_1) + 2 * ((2 -
    2 * (wzdeno * (sigma_sq1) * ewz/wzxsq1^2)) * depsi1 * pmusig1 * sqrt(sigma_sq1)/wzxsq1^2))) *
    ewz + 2 * (s12sq1 * depsi1 * pmusig1) - 2 * sigx11_2) * ewz/(2 * sigx3_2 +
    2 * sigx3_1), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmcesfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/psq2 + ewz1 * depsi1 * pmusig1/psq1)
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx18 <- (pi * sigx2 * ((Wz)^2 + 1))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
      depsi1 + dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv1) * ewu1/starsq1) *
      ewz1/ssqq1 + (((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu2 *
      mupsi2/starsq2) * depsi2 + dmusig2 * (depsi2 * mupsi2 - depsi2 * musi2/ewv2) *
      ewu2/starsq2) * ewz2/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 *
      ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - (sigx8_1 * depsi1 *
      mupsi1/sigma_sq1 + (sigx7_1 * ewu1/starsq1^2 - (sigx8_1 * musi1/ewv1 +
      1/sigmastar1)/sigma_sq1) * depsi1) * dmusig1) * ewu1/psq1 - (sigx3 *
      sigx11_1/sigx2 + (0.5 * wup1 - 0.5 * mud1) * sigx1_1/(sigma_sq1 * psq1^2))) *
      ewz1/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - (sigx8_2 *
      depsi2 * mupsi2/sigma_sq2 + (sigx7_2 * ewu2/starsq2^2 - (sigx8_2 * musi2/ewv2 +
      1/sigmastar2)/sigma_sq2) * depsi2) * dmusig2) * ewu2/psq2 - (sigx3 *
      sigx11_2/sigx2 + (0.5 * wup2 - 0.5 * mud2) * sigx1_2/(sigma_sq2 * psq2^2))) *
      ewz2/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((((depsi1 * mupsi1 - depsi1 * musi1/ewv1) *
      sigx13_1/sigma_sq1 - sigx12_1 * depsi1 * ewu1/starsq1^2) * dmusig1 +
      0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
        depsi1 * mupsi1/sigma_sq1^2))/pmu1 - 0.5 * (sigx1_1 * pmu1/(sigma_sq1 *
      psq1^2)))/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) *
      ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (((((depsi2 *
    mupsi2 - depsi2 * musi2/ewv2) * sigx13_2/sigma_sq2 - sigx12_2 * depsi2 *
    ewu2/starsq2^2) * dmusig2 + 0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 +
    dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2))/pmu2 - 0.5 *
    (sigx1_2 * pmu2/(sigma_sq2 * psq2^2)))/sigx17_2 - sigx3 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) *
    ewz2 * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((sigx1_1/ssqq1 - sigx1_2/ssqq2)/sigx18 - pi * sigx3 * ((Wz)^2 + 1) *
    sigx9/sigx18^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      mupsi1^2/sigma_sq1^2 - ((sigx8_1^2 * ewu1 * musi1/starsq1 + ((0.5 * (ewu1/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu1/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu1)/starsq1^2) * depsi1 + sigx8_1 * (0.5 *
      (depsi1 * ewu1 * mupsi1^2/sigma_sq1^2) + depsi1)) * dmusig1)/psq1 - sigx10_1 *
      (0.5 * wup1 - 0.5 * mud1)/psq1^2) * ewu1 - ((((0.5 * (((1 - 0.5 * (ewu1/sigma_sq1)) *
      pmu1 - 0.5 * (mu1 * dmu1/ewusr1)) * ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 *
      ((0.5 * (mu1^2/ewusr1^2) - 0.5) * sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) *
      dmu1/ewusr1)) * depsi1 + 0.5 * ((0.5 * wup1 - 0.5 * mud1) * depsi1 *
      ewu1 * mupsi1^2/sigma_sq1^2)) * pmusig1 - (sigx8_1 * dmusig1 * ewu1 +
      2 * ((0.5 * wup1 - 0.5 * mud1) * pmusig1 * pmu1 * sqrt(sigma_sq1)/psq1^2)) *
      (0.5 * wup1 - 0.5 * mud1) * depsi1)/psq1^2 + sigx11_1^2 * ewz1/sigx2)) *
      ewz1/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * ewz2 * ewz1/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 * (depsi1 * mupsi1^2/sigma_sq1)) *
      sigx13_1/sigma_sq1 - ((0.5 * (usq1 * ewv1/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 -
      1) * ewv1/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 +
      mu1 * sigx7_1 - sigx12_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
      S * (epsilon))) * depsi1/starsq1^2) * dmusig1 + 0.5 * (((0.5 * (mupsi1^2/sigma_sq1) -
      2) * pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)) *
      ewu1 + 0.5 * (mu1 * sigx15_1 * dmu1/(ewusr1 * pmu1)))/pmu1 - 0.5 * ((sigx10_1 *
      ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewusr1) + 2 * ((0.5 * wup1 - 0.5 * mud1) *
      pmu1^2 * sqrt(sigma_sq1)/psq1^2)) * depsi1 * pmusig1)/psq1^2))/sigx17_1 -
      (sigx11_1 * ewz1 * sqrt(sigma_sq1) + 0.5 * (sigx2 * ewu1/sqrt(sigma_sq1))) *
        sigx16_1/sigx17_1^2) * ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx11_1 * sigx16_2 * ewz2 * ewz1 * ewv2 * sqrt(sigma_sq2)/sigx17_2^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * sigx11_1 * (1/sigx18 - pi * ((Wz)^2 + 1) * ewz1 * sigx9/sigx18^2),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((0.5 *
    (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) *
      ewu2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2 - ((sigx8_2^2 *
    ewu2 * musi2/starsq2 + ((0.5 * (ewu2/sigma_sq2) - 0.5 * (0.5 * usq2 + ewu2/sigma_sq2)) *
    usq2 * ewv2 * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 *
    sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu2)/starsq2^2) * depsi2 +
    sigx8_2 * (0.5 * (depsi2 * ewu2 * mupsi2^2/sigma_sq2^2) + depsi2)) * dmusig2)/psq2 -
    sigx10_2 * (0.5 * wup2 - 0.5 * mud2)/psq2^2) * ewu2 - ((((0.5 * (((1 - 0.5 *
    (ewu2/sigma_sq2)) * pmu2 - 0.5 * (mu2 * dmu2/ewusr2)) * ewu2/sqrt(sigma_sq2)) -
    0.5 * (mu2 * ((0.5 * (mu2^2/ewusr2^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewusr2)) * depsi2 + 0.5 * ((0.5 * wup2 - 0.5 * mud2) * depsi2 *
    ewu2 * mupsi2^2/sigma_sq2^2)) * pmusig2 - (sigx8_2 * dmusig2 * ewu2 + 2 *
    ((0.5 * wup2 - 0.5 * mud2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/psq2^2)) *
    (0.5 * wup2 - 0.5 * mud2) * depsi2)/psq2^2 + sigx11_2^2 * ewz2/sigx2)) *
    ewz2/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx16_1 * sigx11_2 * ewz2 * ewz1 * ewv1 * sqrt(sigma_sq1)/sigx17_1^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ((((((sigx8_2 * depsi2 * musi2/sigmastar2 + 0.5 * (depsi2 *
      mupsi2^2/sigma_sq2)) * sigx13_2/sigma_sq2 - ((0.5 * (usq2 * ewv2/sigma_sq2) +
      0.5 * ((ewu2/sigma_sq2 - 1) * ewv2/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) *
      musi2/sigmastar2 + mu2 * sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 *
      musi2 * sigmastar2/starsq2^2) + S * (epsilon))) * depsi2/starsq2^2) *
      dmusig2 + 0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 -
      sigx8_2 * dmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)) * ewu2 + 0.5 * (mu2 *
      sigx15_2 * dmu2/(ewusr2 * pmu2)))/pmu2 - 0.5 * ((sigx10_2 * ewu2 * pmu2 -
      (0.5 * (mu2 * dmu2/ewusr2) + 2 * ((0.5 * wup2 - 0.5 * mud2) * pmu2^2 *
        sqrt(sigma_sq2)/psq2^2)) * depsi2 * pmusig2)/psq2^2))/sigx17_2 -
      (sigx11_2 * ewz2 * sqrt(sigma_sq2) + 0.5 * (sigx2 * ewu2/sqrt(sigma_sq2))) *
        sigx16_2/sigx17_2^2) * ewz2 * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 2 *
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_2 * (1/sigx18 + pi * ((Wz)^2 + 1) *
      ewz2 * sigx9/sigx18^2)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((((0.5 * (depsi1 * mupsi1^2/sigma_sq1) - depsi1 * musi1 *
      sigx13_1/sigmastar1) * ewv1 * sigx13_1/sigma_sq1 + depsi1 * (mu1/starsq1 -
      (((3 * (mu1) - 2 * (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) *
        ewv1 - S * ewu1 * (epsilon)) * sigx12_1 + (0.5 * (ewv1/sigma_sq1) -
        0.5 * (0.5 * vsq1 + ewv1/sigma_sq1)) * vsq1 * ewu1 * musi1/sigmastar1)/starsq1^2)) *
      dmusig1 + (0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 +
      dmusig1 * sigx13_1) * ewv1) + 0.5 * pmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)/pmu1 -
      (0.5 * (((dmusig1 * sigx13_1 - pmusig1 * pmu1^2/psq1^2) * depsi1 + 0.5 *
        sigx6_1) * ewv1) + 0.5 * (depsi1 * pmusig1)) * pmu1/psq1^2)/sigx17_1 -
      (sigx16_1 * ewz1 + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1 * ewv1/sigx17_1^2) *
      ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_1 * sigx16_2 * ewz2 * ewz1 * ewv1 *
      ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 * sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx16_1 * (1/sigx18 - pi * ((Wz)^2 + 1) * ewz1 *
      sigx9/sigx18^2) * ewv1/sqrt(sigma_sq1), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (depsi2 * mupsi2^2/sigma_sq2) - depsi2 *
      musi2 * sigx13_2/sigmastar2) * ewv2 * sigx13_2/sigma_sq2 + depsi2 * (mu2/starsq2 -
      (((3 * (mu2) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
        ewv2 - S * ewu2 * (epsilon)) * sigx12_2 + (0.5 * (ewv2/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) * vsq2 * ewu2 * musi2/sigmastar2)/starsq2^2)) *
      dmusig2 + (0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 +
      dmusig2 * sigx13_2) * ewv2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)/pmu2 -
      (0.5 * (((dmusig2 * sigx13_2 - pmusig2 * pmu2^2/psq2^2) * depsi2 + 0.5 *
        sigx6_2) * ewv2) + 0.5 * (depsi2 * pmusig2)) * pmu2/psq2^2)/sigx17_2 -
      (sigx16_2 * ewz2 + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2 * ewv2/sigx17_2^2) *
      ewz2 * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (sigx16_2 *
    (1/sigx18 + pi * ((Wz)^2 + 1) * ewz2 * sigx9/sigx18^2) * ewv2/sqrt(sigma_sq2)),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar *
    ((2 * (pi * Wz * sigx2) + depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2) *
      sigx9/sigx18^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmcesfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/psq2 + depsi1 * pmusig1 * pwZ/psq1)
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
      depsi1 + dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv1) * ewu1/starsq1) *
      pwZ/ssqq1 + (((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) *
      depsi2 + dmusig2 * (depsi2 * mupsi2 - depsi2 * musi2/ewv2) * ewu2/starsq2) *
      (1 - pwZ)/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 *
      ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - (sigx8_1 * depsi1 *
      mupsi1/sigma_sq1 + (sigx7_1 * ewu1/starsq1^2 - (sigx8_1 * musi1/ewv1 +
      1/sigmastar1)/sigma_sq1) * depsi1) * dmusig1) * ewu1/psq1 - (sigx3 *
      sigx11_1/sigx2 + (0.5 * wup1 - 0.5 * mud1) * sigx1_1/(sigma_sq1 * psq1^2))) *
      pwZ/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - (sigx8_2 *
      depsi2 * mupsi2/sigma_sq2 + (sigx7_2 * ewu2/starsq2^2 - (sigx8_2 * musi2/ewv2 +
      1/sigmastar2)/sigma_sq2) * depsi2) * dmusig2) * ewu2/psq2 - (sigx3 *
      sigx11_2/sigx2 + (0.5 * wup2 - 0.5 * mud2) * sigx1_2/(sigma_sq2 * psq2^2))) *
      (1 - pwZ)/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((((depsi1 * mupsi1 - depsi1 * musi1/ewv1) *
      sigx13_1/sigma_sq1 - sigx12_1 * depsi1 * ewu1/starsq1^2) * dmusig1 +
      0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
        depsi1 * mupsi1/sigma_sq1^2))/pmu1 - 0.5 * (sigx1_1 * pmu1/(sigma_sq1 *
      psq1^2)))/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) *
      pwZ * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (((((depsi2 *
    mupsi2 - depsi2 * musi2/ewv2) * sigx13_2/sigma_sq2 - sigx12_2 * depsi2 *
    ewu2/starsq2^2) * dmusig2 + 0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 +
    dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2))/pmu2 - 0.5 *
    (sigx1_2 * pmu2/(sigma_sq2 * psq2^2)))/sigx17_2 - sigx3 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) *
    (1 - pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (sigx1_1/ssqq1 - (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) * dwZ/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      mupsi1^2/sigma_sq1^2 - ((sigx8_1^2 * ewu1 * musi1/starsq1 + ((0.5 * (ewu1/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu1/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu1)/starsq1^2) * depsi1 + sigx8_1 * (0.5 *
      (depsi1 * ewu1 * mupsi1^2/sigma_sq1^2) + depsi1)) * dmusig1)/psq1 - sigx10_1 *
      (0.5 * wup1 - 0.5 * mud1)/psq1^2) * ewu1 - ((((0.5 * (((1 - 0.5 * (ewu1/sigma_sq1)) *
      pmu1 - 0.5 * (mu1 * dmu1/ewusr1)) * ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 *
      ((0.5 * (mu1^2/ewusr1^2) - 0.5) * sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) *
      dmu1/ewusr1)) * depsi1 + 0.5 * ((0.5 * wup1 - 0.5 * mud1) * depsi1 *
      ewu1 * mupsi1^2/sigma_sq1^2)) * pmusig1 - (sigx8_1 * dmusig1 * ewu1 +
      2 * ((0.5 * wup1 - 0.5 * mud1) * pmusig1 * pmu1 * sqrt(sigma_sq1)/psq1^2)) *
      (0.5 * wup1 - 0.5 * mud1) * depsi1)/psq1^2 + sigx11_1^2 * pwZ/sigx2)) *
      pwZ/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * (1 - pwZ) * pwZ/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 * (depsi1 * mupsi1^2/sigma_sq1)) *
      sigx13_1/sigma_sq1 - ((0.5 * (usq1 * ewv1/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 -
      1) * ewv1/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 +
      mu1 * sigx7_1 - sigx12_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
      S * (epsilon))) * depsi1/starsq1^2) * dmusig1 + 0.5 * (((0.5 * (mupsi1^2/sigma_sq1) -
      2) * pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)) *
      ewu1 + 0.5 * (mu1 * sigx15_1 * dmu1/(ewusr1 * pmu1)))/pmu1 - 0.5 * ((sigx10_1 *
      ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewusr1) + 2 * ((0.5 * wup1 - 0.5 * mud1) *
      pmu1^2 * sqrt(sigma_sq1)/psq1^2)) * depsi1 * pmusig1)/psq1^2))/sigx17_1 -
      (sigx11_1 * pwZ * sqrt(sigma_sq1) + 0.5 * (sigx2 * ewu1/sqrt(sigma_sq1))) *
        sigx16_1/sigx17_1^2) * pwZ * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx11_1 * sigx16_2 * (1 - pwZ) * pwZ * ewv2 * sqrt(sigma_sq2)/sigx17_2^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * sigx11_1 * (1 - sigx9 * pwZ/sigx2) * dwZ/sigx2, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((0.5 *
    (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) *
      ewu2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2 - ((sigx8_2^2 *
    ewu2 * musi2/starsq2 + ((0.5 * (ewu2/sigma_sq2) - 0.5 * (0.5 * usq2 + ewu2/sigma_sq2)) *
    usq2 * ewv2 * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 *
    sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu2)/starsq2^2) * depsi2 +
    sigx8_2 * (0.5 * (depsi2 * ewu2 * mupsi2^2/sigma_sq2^2) + depsi2)) * dmusig2)/psq2 -
    sigx10_2 * (0.5 * wup2 - 0.5 * mud2)/psq2^2) * ewu2 - ((((0.5 * (((1 - 0.5 *
    (ewu2/sigma_sq2)) * pmu2 - 0.5 * (mu2 * dmu2/ewusr2)) * ewu2/sqrt(sigma_sq2)) -
    0.5 * (mu2 * ((0.5 * (mu2^2/ewusr2^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewusr2)) * depsi2 + 0.5 * ((0.5 * wup2 - 0.5 * mud2) * depsi2 *
    ewu2 * mupsi2^2/sigma_sq2^2)) * pmusig2 - (sigx8_2 * dmusig2 * ewu2 + 2 *
    ((0.5 * wup2 - 0.5 * mud2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/psq2^2)) *
    (0.5 * wup2 - 0.5 * mud2) * depsi2)/psq2^2 + sigx11_2^2 * (1 - pwZ)/sigx2)) *
    (1 - pwZ)/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx16_1 * sigx11_2 * (1 - pwZ) * pwZ * ewv1 * sqrt(sigma_sq1)/sigx17_1^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ((((((sigx8_2 * depsi2 * musi2/sigmastar2 + 0.5 * (depsi2 *
      mupsi2^2/sigma_sq2)) * sigx13_2/sigma_sq2 - ((0.5 * (usq2 * ewv2/sigma_sq2) +
      0.5 * ((ewu2/sigma_sq2 - 1) * ewv2/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) *
      musi2/sigmastar2 + mu2 * sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 *
      musi2 * sigmastar2/starsq2^2) + S * (epsilon))) * depsi2/starsq2^2) *
      dmusig2 + 0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 -
      sigx8_2 * dmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)) * ewu2 + 0.5 * (mu2 *
      sigx15_2 * dmu2/(ewusr2 * pmu2)))/pmu2 - 0.5 * ((sigx10_2 * ewu2 * pmu2 -
      (0.5 * (mu2 * dmu2/ewusr2) + 2 * ((0.5 * wup2 - 0.5 * mud2) * pmu2^2 *
        sqrt(sigma_sq2)/psq2^2)) * depsi2 * pmusig2)/psq2^2))/sigx17_2 -
      (sigx11_2 * (1 - pwZ) * sqrt(sigma_sq2) + 0.5 * (sigx2 * ewu2/sqrt(sigma_sq2))) *
        sigx16_2/sigx17_2^2) * (1 - pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 2 *
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_2 * ((1 - pwZ) * sigx9/sigx2 + 1) *
      dwZ/sigx2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((((0.5 * (depsi1 * mupsi1^2/sigma_sq1) - depsi1 * musi1 *
      sigx13_1/sigmastar1) * ewv1 * sigx13_1/sigma_sq1 + depsi1 * (mu1/starsq1 -
      (((3 * (mu1) - 2 * (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) *
        ewv1 - S * ewu1 * (epsilon)) * sigx12_1 + (0.5 * (ewv1/sigma_sq1) -
        0.5 * (0.5 * vsq1 + ewv1/sigma_sq1)) * vsq1 * ewu1 * musi1/sigmastar1)/starsq1^2)) *
      dmusig1 + (0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 +
      dmusig1 * sigx13_1) * ewv1) + 0.5 * pmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)/pmu1 -
      (0.5 * (((dmusig1 * sigx13_1 - pmusig1 * pmu1^2/psq1^2) * depsi1 + 0.5 *
        sigx6_1) * ewv1) + 0.5 * (depsi1 * pmusig1)) * pmu1/psq1^2)/sigx17_1 -
      (sigx16_1 * pwZ + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1 * ewv1/sigx17_1^2) *
      pwZ * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_1 * sigx16_2 * (1 - pwZ) * pwZ * ewv1 *
      ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 * sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx16_1 * (1 - sigx9 * pwZ/sigx2) * dwZ * ewv1/(sigx2 *
      sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (depsi2 * mupsi2^2/sigma_sq2) - depsi2 *
      musi2 * sigx13_2/sigmastar2) * ewv2 * sigx13_2/sigma_sq2 + depsi2 * (mu2/starsq2 -
      (((3 * (mu2) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
        ewv2 - S * ewu2 * (epsilon)) * sigx12_2 + (0.5 * (ewv2/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) * vsq2 * ewu2 * musi2/sigmastar2)/starsq2^2)) *
      dmusig2 + (0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 +
      dmusig2 * sigx13_2) * ewv2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)/pmu2 -
      (0.5 * (((dmusig2 * sigx13_2 - pmusig2 * pmu2^2/psq2^2) * depsi2 + 0.5 *
        sigx6_2) * ewv2) + 0.5 * (depsi2 * pmusig2)) * pmu2/psq2^2)/sigx17_2 -
      (sigx16_2 * (1 - pwZ) + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2 * ewv2/sigx17_2^2) *
      (1 - pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (sigx16_2 *
    ((1 - pwZ) * sigx9/sigx2 + 1) * dwZ * ewv2/(sigx2 * sqrt(sigma_sq2))), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar *
    ((sigx9 * dwZ/sigx2 + Wz) * sigx9 * dwZ/sigx2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmcesfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/psq1 + depsi2 * prZ * pmusig2/psq2)
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
      depsi1 + dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv1) * ewu1/starsq1) *
      (1 - prZ)/ssqq1 + (((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu2 *
      mupsi2/starsq2) * depsi2 + dmusig2 * (depsi2 * mupsi2 - depsi2 * musi2/ewv2) *
      ewu2/starsq2) * prZ/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 *
      ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - (sigx8_1 * depsi1 *
      mupsi1/sigma_sq1 + (sigx7_1 * ewu1/starsq1^2 - (sigx8_1 * musi1/ewv1 +
      1/sigmastar1)/sigma_sq1) * depsi1) * dmusig1) * ewu1/psq1 - (sigx3 *
      sigx11_1/sigx2 + (0.5 * wup1 - 0.5 * mud1) * sigx1_1/(sigma_sq1 * psq1^2))) *
      (1 - prZ)/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - (sigx8_2 *
      depsi2 * mupsi2/sigma_sq2 + (sigx7_2 * ewu2/starsq2^2 - (sigx8_2 * musi2/ewv2 +
      1/sigmastar2)/sigma_sq2) * depsi2) * dmusig2) * ewu2/psq2 - (sigx3 *
      sigx11_2/sigx2 + (0.5 * wup2 - 0.5 * mud2) * sigx1_2/(sigma_sq2 * psq2^2))) *
      prZ/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((((depsi1 * mupsi1 - depsi1 * musi1/ewv1) *
      sigx13_1/sigma_sq1 - sigx12_1 * depsi1 * ewu1/starsq1^2) * dmusig1 +
      0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
        depsi1 * mupsi1/sigma_sq1^2))/pmu1 - 0.5 * (sigx1_1 * pmu1/(sigma_sq1 *
      psq1^2)))/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) *
      (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (((((depsi2 *
    mupsi2 - depsi2 * musi2/ewv2) * sigx13_2/sigma_sq2 - sigx12_2 * depsi2 *
    ewu2/starsq2^2) * dmusig2 + 0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 +
    dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2))/pmu2 - 0.5 *
    (sigx1_2 * pmu2/(sigma_sq2 * psq2^2)))/sigx17_2 - sigx3 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) *
    prZ * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (sigx1_1/ssqq1 - (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) * prZ * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      mupsi1^2/sigma_sq1^2 - ((sigx8_1^2 * ewu1 * musi1/starsq1 + ((0.5 * (ewu1/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu1/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu1)/starsq1^2) * depsi1 + sigx8_1 * (0.5 *
      (depsi1 * ewu1 * mupsi1^2/sigma_sq1^2) + depsi1)) * dmusig1)/psq1 - sigx10_1 *
      (0.5 * wup1 - 0.5 * mud1)/psq1^2) * ewu1 - ((((0.5 * (((1 - 0.5 * (ewu1/sigma_sq1)) *
      pmu1 - 0.5 * (mu1 * dmu1/ewusr1)) * ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 *
      ((0.5 * (mu1^2/ewusr1^2) - 0.5) * sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) *
      dmu1/ewusr1)) * depsi1 + 0.5 * ((0.5 * wup1 - 0.5 * mud1) * depsi1 *
      ewu1 * mupsi1^2/sigma_sq1^2)) * pmusig1 - (sigx8_1 * dmusig1 * ewu1 +
      2 * ((0.5 * wup1 - 0.5 * mud1) * pmusig1 * pmu1 * sqrt(sigma_sq1)/psq1^2)) *
      (0.5 * wup1 - 0.5 * mud1) * depsi1)/psq1^2 + sigx11_1^2 * (1 - prZ)/sigx2)) *
      (1 - prZ)/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * prZ * (1 - prZ)/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 * (depsi1 * mupsi1^2/sigma_sq1)) *
      sigx13_1/sigma_sq1 - ((0.5 * (usq1 * ewv1/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 -
      1) * ewv1/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 +
      mu1 * sigx7_1 - sigx12_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
      S * (epsilon))) * depsi1/starsq1^2) * dmusig1 + 0.5 * (((0.5 * (mupsi1^2/sigma_sq1) -
      2) * pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)) *
      ewu1 + 0.5 * (mu1 * sigx15_1 * dmu1/(ewusr1 * pmu1)))/pmu1 - 0.5 * ((sigx10_1 *
      ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewusr1) + 2 * ((0.5 * wup1 - 0.5 * mud1) *
      pmu1^2 * sqrt(sigma_sq1)/psq1^2)) * depsi1 * pmusig1)/psq1^2))/sigx17_1 -
      (sigx11_1 * (1 - prZ) * sqrt(sigma_sq1) + 0.5 * (sigx2 * ewu1/sqrt(sigma_sq1))) *
        sigx16_1/sigx17_1^2) * (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx11_1 * sigx16_2 * prZ * (1 - prZ) * ewv2 * sqrt(sigma_sq2)/sigx17_2^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * sigx11_1 * (1 - (1 - prZ) * sigx9/sigx2) * prZ * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((0.5 *
    (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) *
      ewu2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2 - ((sigx8_2^2 *
    ewu2 * musi2/starsq2 + ((0.5 * (ewu2/sigma_sq2) - 0.5 * (0.5 * usq2 + ewu2/sigma_sq2)) *
    usq2 * ewv2 * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 *
    sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu2)/starsq2^2) * depsi2 +
    sigx8_2 * (0.5 * (depsi2 * ewu2 * mupsi2^2/sigma_sq2^2) + depsi2)) * dmusig2)/psq2 -
    sigx10_2 * (0.5 * wup2 - 0.5 * mud2)/psq2^2) * ewu2 - ((((0.5 * (((1 - 0.5 *
    (ewu2/sigma_sq2)) * pmu2 - 0.5 * (mu2 * dmu2/ewusr2)) * ewu2/sqrt(sigma_sq2)) -
    0.5 * (mu2 * ((0.5 * (mu2^2/ewusr2^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewusr2)) * depsi2 + 0.5 * ((0.5 * wup2 - 0.5 * mud2) * depsi2 *
    ewu2 * mupsi2^2/sigma_sq2^2)) * pmusig2 - (sigx8_2 * dmusig2 * ewu2 + 2 *
    ((0.5 * wup2 - 0.5 * mud2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/psq2^2)) *
    (0.5 * wup2 - 0.5 * mud2) * depsi2)/psq2^2 + sigx11_2^2 * prZ/sigx2)) * prZ/sigx2,
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx16_1 * sigx11_2 * prZ * (1 - prZ) * ewv1 * sqrt(sigma_sq1)/sigx17_1^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ((((((sigx8_2 * depsi2 * musi2/sigmastar2 + 0.5 * (depsi2 *
      mupsi2^2/sigma_sq2)) * sigx13_2/sigma_sq2 - ((0.5 * (usq2 * ewv2/sigma_sq2) +
      0.5 * ((ewu2/sigma_sq2 - 1) * ewv2/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) *
      musi2/sigmastar2 + mu2 * sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 *
      musi2 * sigmastar2/starsq2^2) + S * (epsilon))) * depsi2/starsq2^2) *
      dmusig2 + 0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 -
      sigx8_2 * dmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)) * ewu2 + 0.5 * (mu2 *
      sigx15_2 * dmu2/(ewusr2 * pmu2)))/pmu2 - 0.5 * ((sigx10_2 * ewu2 * pmu2 -
      (0.5 * (mu2 * dmu2/ewusr2) + 2 * ((0.5 * wup2 - 0.5 * mud2) * pmu2^2 *
        sqrt(sigma_sq2)/psq2^2)) * depsi2 * pmusig2)/psq2^2))/sigx17_2 -
      (sigx11_2 * prZ * sqrt(sigma_sq2) + 0.5 * (sigx2 * ewu2/sqrt(sigma_sq2))) *
        sigx16_2/sigx17_2^2) * prZ * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 2 *
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_2 * (sigx9 * prZ/sigx2 + 1) * prZ *
      ewz/sigx2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((((0.5 * (depsi1 * mupsi1^2/sigma_sq1) - depsi1 * musi1 *
      sigx13_1/sigmastar1) * ewv1 * sigx13_1/sigma_sq1 + depsi1 * (mu1/starsq1 -
      (((3 * (mu1) - 2 * (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) *
        ewv1 - S * ewu1 * (epsilon)) * sigx12_1 + (0.5 * (ewv1/sigma_sq1) -
        0.5 * (0.5 * vsq1 + ewv1/sigma_sq1)) * vsq1 * ewu1 * musi1/sigmastar1)/starsq1^2)) *
      dmusig1 + (0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 +
      dmusig1 * sigx13_1) * ewv1) + 0.5 * pmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)/pmu1 -
      (0.5 * (((dmusig1 * sigx13_1 - pmusig1 * pmu1^2/psq1^2) * depsi1 + 0.5 *
        sigx6_1) * ewv1) + 0.5 * (depsi1 * pmusig1)) * pmu1/psq1^2)/sigx17_1 -
      (sigx16_1 * (1 - prZ) + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1 * ewv1/sigx17_1^2) *
      (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_1 * sigx16_2 * prZ * (1 - prZ) * ewv1 *
      ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 * sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx16_1 * (1 - (1 - prZ) * sigx9/sigx2) * prZ *
      ewv1 * ewz/(sigx2 * sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (depsi2 * mupsi2^2/sigma_sq2) - depsi2 *
      musi2 * sigx13_2/sigmastar2) * ewv2 * sigx13_2/sigma_sq2 + depsi2 * (mu2/starsq2 -
      (((3 * (mu2) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
        ewv2 - S * ewu2 * (epsilon)) * sigx12_2 + (0.5 * (ewv2/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) * vsq2 * ewu2 * musi2/sigmastar2)/starsq2^2)) *
      dmusig2 + (0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 +
      dmusig2 * sigx13_2) * ewv2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)/pmu2 -
      (0.5 * (((dmusig2 * sigx13_2 - pmusig2 * pmu2^2/psq2^2) * depsi2 + 0.5 *
        sigx6_2) * ewv2) + 0.5 * (depsi2 * pmusig2)) * pmu2/psq2^2)/sigx17_2 -
      (sigx16_2 * prZ + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2 * ewv2/sigx17_2^2) *
      prZ * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (sigx16_2 *
    (sigx9 * prZ/sigx2 + 1) * prZ * ewv2 * ewz/(sigx2 * sqrt(sigma_sq2))), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar * (1 -
    (sigx9 * prZ/sigx2 + 1) * ewz) * sigx9 * prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf halfnormal-normal distribution
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
# same sigma_u

## logit specification class membership
cnsfhalfnormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ccnsfhalfnormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ccnsfhalfnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsfhalfnormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ccnsfhalfnormlike_logit,
    grad = cgradcnsfhalfnormlike_logit, hess = chesscnsfhalfnormlike_logit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsfhalfnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(ccnsfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(ccnsfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesscnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfhalfnormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesscnsfhalfnormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfhalfnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfhalfnormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfhalfnormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

## cauchit specification class membership
cnsfhalfnormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ccnsfhalfnormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ccnsfhalfnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsfhalfnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ccnsfhalfnormlike_cauchit,
    grad = cgradcnsfhalfnormlike_cauchit, hess = chesscnsfhalfnormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsfhalfnormlike_cauchit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(ccnsfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(ccnsfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesscnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfhalfnormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesscnsfhalfnormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfhalfnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfhalfnormlike_cauchit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfhalfnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

## probit specification class membership
cnsfhalfnormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ccnsfhalfnormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ccnsfhalfnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsfhalfnormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ccnsfhalfnormlike_probit,
    grad = cgradcnsfhalfnormlike_probit, hess = chesscnsfhalfnormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsfhalfnormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(ccnsfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(ccnsfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesscnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfhalfnormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesscnsfhalfnormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfhalfnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfhalfnormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfhalfnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

## cloglog specification class membership
cnsfhalfnormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ccnsfhalfnormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ccnsfhalfnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsfhalfnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ccnsfhalfnormlike_cloglog,
    grad = cgradcnsfhalfnormlike_cloglog, hess = chesscnsfhalfnormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsfhalfnormlike_cloglog(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(ccnsfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(ccnsfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesscnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfhalfnormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesscnsfhalfnormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfhalfnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfhalfnormlike_cloglog(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfhalfnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

# Different sigma_u

## logit specification class membership
mcesfhalfnormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmcesfhalfnormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmcesfhalfnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesfhalfnormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmcesfhalfnormlike_logit,
    grad = cgradmcesfhalfnormlike_logit, hess = chessmcesfhalfnormlike_logit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesfhalfnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmcesfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmcesfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmcesfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfhalfnormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmcesfhalfnormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfhalfnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfhalfnormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfhalfnormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

## cauchit specification class membership
mcesfhalfnormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmcesfhalfnormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmcesfhalfnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesfhalfnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmcesfhalfnormlike_cauchit,
    grad = cgradmcesfhalfnormlike_cauchit, hess = chessmcesfhalfnormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesfhalfnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmcesfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmcesfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmcesfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfhalfnormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmcesfhalfnormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfhalfnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfhalfnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfhalfnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

## probit specification class membership
mcesfhalfnormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmcesfhalfnormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmcesfhalfnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesfhalfnormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmcesfhalfnormlike_probit,
    grad = cgradmcesfhalfnormlike_probit, hess = chessmcesfhalfnormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesfhalfnormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmcesfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmcesfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmcesfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfhalfnormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmcesfhalfnormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfhalfnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfhalfnormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfhalfnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

## cloglog specification class membership
mcesfhalfnormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmcesfhalfnormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmcesfhalfnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesfhalfnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmcesfhalfnormlike_cloglog,
    grad = cgradmcesfhalfnormlike_cloglog, hess = chessmcesfhalfnormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesfhalfnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmcesfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmcesfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmcesfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfhalfnormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmcesfhalfnormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfhalfnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfhalfnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfhalfnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

# Conditional efficiencies estimation ----------
#' efficiencies for cnsf halfnormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u

## logit specification class membership
ccnsfhalfnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

## cauchit specification class membership
ccnsfhalfnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

## probit specification class membership
ccnsfhalfnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

## cloglog specification class membership
ccnsfhalfnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

# different sigma_u

## logit specification class membership
cmcesfhalfnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

## cauchit specification class membership
cmcesfhalfnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

## probit specification class membership
cmcesfhalfnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

## cloglog specification class membership
cmcesfhalfnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for cnsf halfnormal-normal distribution
#' @param object object of class sfacross
#' @noRd
# same sigma_u

## logit specification class membership
ccnsfmarghalfnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmarghalfnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
ccnsfmarghalfnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmarghalfnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
ccnsfmarghalfnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmarghalfnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
ccnsfmarghalfnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmarghalfnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# different sigma_u

## logit specification class membership
cmcesfmarghalfnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmarghalfnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmcesfmarghalfnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmarghalfnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmcesfmarghalfnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmarghalfnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmcesfmarghalfnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmarghalfnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

