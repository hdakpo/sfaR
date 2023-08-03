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
# Convolution: rayleigh - normal                                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf rayleigh-normal distribution
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
czisfraynormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cauchit specification class membership
czisfraynormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## probit specification class membership
czisfraynormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cloglog specification class membership
czisfraynormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

# Different sigma_v

## logit specification class membership
cmnsfraynormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, S, wHvar, Zvar, nZHvar) {
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
cmnsfraynormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
cmnsfraynormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
cmnsfraynormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf rayleigh-normal distribution
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
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
# Same sigma_v
cstzisfraynorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- cstraynorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initRay <- NULL

  } else {
    cat("Initialization: SFA + rayleigh - normal distributions...\n")
    initRay <- maxLik::maxLik(logLik = craynormlike, start = cstraynorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradraynormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initRay$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("ZISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initRay = initRay))
}

# Different sigma_v
cstmnsfraynorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- cstraynorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initRay <- NULL

  } else {
    cat("Initialization: SFA + rayleigh - normal distributions...\n")
    initRay <- maxLik::maxLik(logLik = craynormlike, start = cstraynorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradraynormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initRay$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar +
    2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("ZISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initRay = initRay))
}

# Gradient of the likelihood function ----------
#' gradient for zisf rayleigh-normal distribution
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
cgradzisfraynormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  dmusig <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  svsu <- (sigmastar/ewv - ewu/(ssq))
  prU <- (1 - ewu/(sigma_sq))
  prV <- (1 - ewv/(sigma_sq))
  dvs <- (dmusig * ewv/sigmastar)
  wzsq <- (wzdeno * (sigma_sq))
  ssg <- (0.5 * (prU * ewv/sigmastar) + sigmastar)
  sgsq <- (1/(ssq) - ssg * ewu/(ssq)^2)
  prvsig <- (prV * dmusig/sigmastar)
  prusig <- (prV * ewu/sigmastar)
  wzwz <- (1/(wzdeno * ewu) - ewu * ewz/(wzdeno * ewu)^2)
  euepsi <- exp(0.5 * (-(S * ewu * (epsilon)/(ssq)))^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1 <- (pmusig + S * dmusig * svsu * (epsilon))
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * (epsilon)/(sigma_sq))
  sigx3 <- (ewu * sigx1/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv)
  sigx4 <- (euepsi * sigx3 * ewz * sigmastar/(wzdeno * ewu) + S * prC * dwsr *
    (epsilon)/ewvsr^2)
  sigx5 <- (prC * dwsr + sigx2 * euepsi * ewz * sigmastar/(wzdeno * ewu))
  sigx6 <- (0.5 * dvs - S * pmusig * (epsilon))
  sigx7 <- (sigx6 * prU + S^2 * sgsq * sigx2 * ewu * (epsilon)^2/sigmastar)
  sigx8 <- (sigx7 * sigmastar + 0.5 * (prU * sigx2 * ewv/sigmastar))
  sigx9 <- (sigx8/wzsq - wzdeno * sigx2 * ewu * sigmastar/(wzdeno * ewu)^2)
  sigx10 <- (0.5 * prvsig + S * pmusig * (epsilon)/(sigma_sq))
  sigx11 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx12 <- (2 * sigx11 - S^2 * (0.5 * prusig + sigmastar) * ewu^2 * (epsilon)^2/((ssq)^2 *
    ssq))
  sigx13 <- (sigx10 * ewu/(sigma_sq) + sigx12 * sigx2)
  sigx14 <- (sigx13 * sigmastar + 0.5 * (prV * sigx2 * ewu/(ssq)))
  sigx15 <- (sigx14 * euepsi * ewv * ewz/(wzdeno * ewu) + 0.5 * (S^2 * prC * dwsr *
    (epsilon)^2/ewvsr^2))
  sigx16 <- (wzwz * sigx2 * euepsi * sigmastar - prC * dwsr/wzdeno)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx5, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx9 * euepsi * ewz/sigx5, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (sigx15/sigx5 - 0.5), FUN = "*"), sweep(Zvar,
      MARGIN = 1, STATS = sigx16 * ewz/sigx5, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradzisfraynormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  dmusig <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  svsu <- (sigmastar/ewv - ewu/(ssq))
  prU <- (1 - ewu/(sigma_sq))
  prV <- (1 - ewv/(sigma_sq))
  dvs <- (dmusig * ewv/sigmastar)
  ssg <- (0.5 * (prU * ewv/sigmastar) + sigmastar)
  sgsq <- (1/(ssq) - ssg * ewu/(ssq)^2)
  prvsig <- (prV * dmusig/sigmastar)
  prusig <- (prV * ewu/sigmastar)
  euepsi <- exp(0.5 * (-(S * ewu * (epsilon)/(ssq)))^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1 <- (pmusig + S * dmusig * svsu * (epsilon))
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * (epsilon)/(sigma_sq))
  sigx3 <- (ewu * sigx1/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv)
  sigx4 <- (ewz1 * euepsi * sigx3 * sigmastar/ewu + S * ewz2 * dwsr * (epsilon)/ewvsr^2)
  sigx5 <- (ewz2 * dwsr + ewz1 * sigx2 * euepsi * sigmastar/ewu)
  sigx6 <- (0.5 * dvs - S * pmusig * (epsilon))
  sigx7 <- (sigx6 * prU + S^2 * sgsq * sigx2 * ewu * (epsilon)^2/sigmastar)
  sigx8 <- (sigx7 * sigmastar + 0.5 * (prU * sigx2 * ewv/sigmastar))
  sigx9 <- (sigx8/(sigma_sq) - sigx2 * sigmastar/ewu)
  sigx10 <- (0.5 * prvsig + S * pmusig * (epsilon)/(sigma_sq))
  sigx11 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx12 <- (2 * sigx11 - S^2 * (0.5 * prusig + sigmastar) * ewu^2 * (epsilon)^2/((ssq)^2 *
    ssq))
  sigx13 <- (sigx10 * ewu/(sigma_sq) + sigx12 * sigx2)
  sigx14 <- (sigx13 * sigmastar + 0.5 * (prV * sigx2 * ewu/(ssq)))
  sigx15 <- (sigx14 * ewz1 * euepsi * ewv/ewu + 0.5 * (S^2 * ewz2 * dwsr * (epsilon)^2/ewvsr^2))
  sigx16 <- (sigx2 * euepsi * sigmastar/ewu - dwsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx5, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx9 * ewz1 * euepsi/sigx5, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (sigx15/sigx5 - 0.5), FUN = "*"), sweep(Zvar,
      MARGIN = 1, STATS = sigx16/(pi * sigx5 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradzisfraynormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  dmusig <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  svsu <- (sigmastar/ewv - ewu/(ssq))
  prU <- (1 - ewu/(sigma_sq))
  prV <- (1 - ewv/(sigma_sq))
  dvs <- (dmusig * ewv/sigmastar)
  ssg <- (0.5 * (prU * ewv/sigmastar) + sigmastar)
  sgsq <- (1/(ssq) - ssg * ewu/(ssq)^2)
  prvsig <- (prV * dmusig/sigmastar)
  prusig <- (prV * ewu/sigmastar)
  euepsi <- exp(0.5 * (-(S * ewu * (epsilon)/(ssq)))^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1 <- (pmusig + S * dmusig * svsu * (epsilon))
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * (epsilon)/(sigma_sq))
  sigx3 <- (ewu * sigx1/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv)
  sigx4 <- (euepsi * sigx3 * pwZ * sigmastar/ewu + S * (1 - pwZ) * dwsr * (epsilon)/ewvsr^2)
  sigx5 <- ((1 - pwZ) * dwsr + sigx2 * euepsi * pwZ * sigmastar/ewu)
  sigx6 <- (0.5 * dvs - S * pmusig * (epsilon))
  sigx7 <- (sigx6 * prU + S^2 * sgsq * sigx2 * ewu * (epsilon)^2/sigmastar)
  sigx8 <- (sigx7 * sigmastar + 0.5 * (prU * sigx2 * ewv/sigmastar))
  sigx9 <- (sigx8/(sigma_sq) - sigx2 * sigmastar/ewu)
  sigx10 <- (0.5 * prvsig + S * pmusig * (epsilon)/(sigma_sq))
  sigx11 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx12 <- (2 * sigx11 - S^2 * (0.5 * prusig + sigmastar) * ewu^2 * (epsilon)^2/((ssq)^2 *
    ssq))
  sigx13 <- (sigx10 * ewu/(sigma_sq) + sigx12 * sigx2)
  sigx14 <- (sigx13 * sigmastar + 0.5 * (prV * sigx2 * ewu/(ssq)))
  sigx15 <- (sigx14 * euepsi * ewv * pwZ/ewu + 0.5 * (S^2 * (1 - pwZ) * dwsr *
    (epsilon)^2/ewvsr^2))
  sigx16 <- (sigx2 * euepsi * sigmastar/ewu - dwsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx5, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx9 * euepsi * pwZ/sigx5, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (sigx15/sigx5 - 0.5), FUN = "*"), sweep(Zvar,
      MARGIN = 1, STATS = sigx16 * dwZ/sigx5, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradzisfraynormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  dmusig <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  svsu <- (sigmastar/ewv - ewu/(ssq))
  prU <- (1 - ewu/(sigma_sq))
  prV <- (1 - ewv/(sigma_sq))
  dvs <- (dmusig * ewv/sigmastar)
  ssg <- (0.5 * (prU * ewv/sigmastar) + sigmastar)
  sgsq <- (1/(ssq) - ssg * ewu/(ssq)^2)
  prvsig <- (prV * dmusig/sigmastar)
  prusig <- (prV * ewu/sigmastar)
  euepsi <- exp(0.5 * (-(S * ewu * (epsilon)/(ssq)))^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1 <- (pmusig + S * dmusig * svsu * (epsilon))
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * (epsilon)/(sigma_sq))
  sigx3 <- (ewu * sigx1/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv)
  sigx4 <- ((1 - prZ) * euepsi * sigx3 * sigmastar/ewu + S * dwsr * prZ * (epsilon)/ewvsr^2)
  sigx5 <- ((1 - prZ) * sigx2 * euepsi * sigmastar/ewu + dwsr * prZ)
  sigx6 <- (0.5 * dvs - S * pmusig * (epsilon))
  sigx7 <- (sigx6 * prU + S^2 * sgsq * sigx2 * ewu * (epsilon)^2/sigmastar)
  sigx8 <- (sigx7 * sigmastar + 0.5 * (prU * sigx2 * ewv/sigmastar))
  sigx9 <- (sigx8/(sigma_sq) - sigx2 * sigmastar/ewu)
  sigx10 <- (0.5 * prvsig + S * pmusig * (epsilon)/(sigma_sq))
  sigx11 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx12 <- (2 * sigx11 - S^2 * (0.5 * prusig + sigmastar) * ewu^2 * (epsilon)^2/((ssq)^2 *
    ssq))
  sigx13 <- (sigx10 * ewu/(sigma_sq) + sigx12 * sigx2)
  sigx14 <- (sigx13 * sigmastar + 0.5 * (prV * sigx2 * ewu/(ssq)))
  sigx15 <- (sigx14 * (1 - prZ) * euepsi * ewv/ewu + 0.5 * (S^2 * dwsr * prZ *
    (epsilon)^2/ewvsr^2))
  sigx16 <- (sigx2 * euepsi * sigmastar/ewu - dwsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx5, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx9 * (1 - prZ) * euepsi/sigx5, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (sigx15/sigx5 - 0.5), FUN = "*"), sweep(Zvar,
      MARGIN = 1, STATS = sigx16 * prZ * ewz/sigx5, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v

## logit specification class membership
cgradmnsfraynormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  sqsq <- ((sigma_sq) * sigmastar)
  mustar <- (S * ewu * (epsilon)/sqsq)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wzuv <- (wzdeno * ewu * ewv1_h)
  sigx1 <- exp(0.5 * (-mustar)^2 - (S * (epsilon))^2/(2 * ewv1))
  sigx2 <- (sigmastar/ewv1 - ewu/sqsq)
  sigx3 <- (pmustar + S * dmustar * sigx2 * (epsilon))
  sigx4 <- (dmustar * sigmastar - S * ewu * pmustar * (epsilon)/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx5 <- (ewu * sigx3/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1)
  sigx6 <- (sigx1 * sigx5 * ewz * sigmastar/wzuv + S * prC * depsi * (epsilon)/ewv2_h^3)
  sigx7 <- (prC * depsi/ewv2_h + sigx4 * sigx1 * ewz * sigmastar/wzuv)
  sigx8 <- (0.5 * (dmustar * ewv1/sigmastar) - S * pmustar * (epsilon))
  sigx9 <- (prU * ewv1/sigmastar)
  sigx10 <- (1/sqsq - (0.5 * sigx9 + sigmastar) * ewu/sqsq^2)
  sigx11 <- (sigx8 * prU + S^2 * sigx10 * sigx4 * ewu * (epsilon)^2/sigmastar)
  sigx12 <- (prU * sigx4 * ewv1/sigmastar)
  sigx13 <- (sigx11 * sigmastar + 0.5 * sigx12)
  wsqv <- (wzdeno * (sigma_sq) * ewv1_h)
  sigx14 <- (sigx13/wsqv - wzdeno * sigx4 * ewu * ewv1_h * sigmastar/wzuv^2)
  prV <- (1 - ewv1/(sigma_sq))
  sigx15 <- (0.5 * (prV * dmustar/sigmastar) + S * pmustar * (epsilon)/(sigma_sq))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx17 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx18 <- (sqsq^2 * (sigma_sq) * sigmastar)
  sigx19 <- (2 * sigx16 - S^2 * sigx17 * ewu^2 * (epsilon)^2/sigx18)
  sigx20 <- (sigx15 * ewu/(sigma_sq) + sigx19 * sigx4)
  sigx21 <- (sigx20 * sigmastar + 0.5 * (prV * sigx4 * ewu/sqsq))
  sigx22 <- (wzdeno * sigx4 * ewu * ewv1_h * sigmastar/wzuv^2)
  sigx23 <- (sigx21 * ewv1/wzuv - 0.5 * sigx22)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  sigx25 <- ((1/wzuv - ewu * ewv1_h * ewz/wzuv^2) * sigx4 * sigx1 * sigmastar -
    prC * depsi/(wzdeno * ewv2_h))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx6/sigx7, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx14 * sigx1 * ewz/sigx7, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx23 * sigx1 * ewz/sigx7, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx24 * prC/(sigx7 * ewv2_h), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx25 * ewz/sigx7, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmnsfraynormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  sqsq <- ((sigma_sq) * sigmastar)
  mustar <- (S * ewu * (epsilon)/sqsq)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigx1 <- exp(0.5 * (-mustar)^2 - (S * (epsilon))^2/(2 * ewv1))
  sigx2 <- (sigmastar/ewv1 - ewu/sqsq)
  sigx3 <- (pmustar + S * dmustar * sigx2 * (epsilon))
  sigx4 <- (dmustar * sigmastar - S * ewu * pmustar * (epsilon)/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx5 <- (ewu * sigx3/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1)
  sigx6 <- (ewz1 * sigx1 * sigx5 * sigmastar/(ewu * ewv1_h) + S * ewz2 * depsi *
    (epsilon)/ewv2_h^3)
  sigx7 <- (ewz2 * depsi/ewv2_h + ewz1 * sigx4 * sigx1 * sigmastar/(ewu * ewv1_h))
  sigx8 <- (0.5 * (dmustar * ewv1/sigmastar) - S * pmustar * (epsilon))
  sigx9 <- (prU * ewv1/sigmastar)
  sigx10 <- (1/sqsq - (0.5 * sigx9 + sigmastar) * ewu/sqsq^2)
  sigx11 <- (sigx8 * prU + S^2 * sigx10 * sigx4 * ewu * (epsilon)^2/sigmastar)
  sigx12 <- (prU * sigx4 * ewv1/sigmastar)
  sigx13 <- (sigx11 * sigmastar + 0.5 * sigx12)
  sigx14 <- (sigx13/((sigma_sq) * ewv1_h) - sigx4 * ewu * ewv1_h * sigmastar/(ewu *
    ewv1_h)^2)
  prV <- (1 - ewv1/(sigma_sq))
  sigx15 <- (0.5 * (prV * dmustar/sigmastar) + S * pmustar * (epsilon)/(sigma_sq))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx17 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx18 <- (sqsq^2 * (sigma_sq) * sigmastar)
  sigx19 <- (2 * sigx16 - S^2 * sigx17 * ewu^2 * (epsilon)^2/sigx18)
  sigx20 <- (sigx15 * ewu/(sigma_sq) + sigx19 * sigx4)
  sigx21 <- (sigx20 * sigmastar + 0.5 * (prV * sigx4 * ewu/sqsq))
  sigx22 <- (sigx4 * ewu * ewv1_h * sigmastar/(ewu * ewv1_h)^2)
  sigx23 <- (sigx21 * ewv1/(ewu * ewv1_h) - 0.5 * sigx22)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  sigx25 <- (sigx4 * sigx1 * sigmastar/(ewu * ewv1_h) - depsi/ewv2_h)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx6/sigx7, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx14 * ewz1 * sigx1/sigx7, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx23 * ewz1 * sigx1/sigx7, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = ewz2 * sigx24/(sigx7 * ewv2_h), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx25/(pi * sigx7 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmnsfraynormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  sqsq <- ((sigma_sq) * sigmastar)
  mustar <- (S * ewu * (epsilon)/sqsq)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigx1 <- exp(0.5 * (-mustar)^2 - (S * (epsilon))^2/(2 * ewv1))
  sigx2 <- (sigmastar/ewv1 - ewu/sqsq)
  sigx3 <- (pmustar + S * dmustar * sigx2 * (epsilon))
  sigx4 <- (dmustar * sigmastar - S * ewu * pmustar * (epsilon)/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx5 <- (ewu * sigx3/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1)
  sigx6 <- (sigx1 * sigx5 * pwZ * sigmastar/(ewu * ewv1_h) + S * (1 - pwZ) * depsi *
    (epsilon)/ewv2_h^3)
  sigx7 <- ((1 - pwZ) * depsi/ewv2_h + sigx4 * sigx1 * pwZ * sigmastar/(ewu * ewv1_h))
  sigx8 <- (0.5 * (dmustar * ewv1/sigmastar) - S * pmustar * (epsilon))
  sigx9 <- (prU * ewv1/sigmastar)
  sigx10 <- (1/sqsq - (0.5 * sigx9 + sigmastar) * ewu/sqsq^2)
  sigx11 <- (sigx8 * prU + S^2 * sigx10 * sigx4 * ewu * (epsilon)^2/sigmastar)
  sigx12 <- (prU * sigx4 * ewv1/sigmastar)
  sigx13 <- (sigx11 * sigmastar + 0.5 * sigx12)
  sigx14 <- (sigx13/((sigma_sq) * ewv1_h) - sigx4 * ewu * ewv1_h * sigmastar/(ewu *
    ewv1_h)^2)
  prV <- (1 - ewv1/(sigma_sq))
  sigx15 <- (0.5 * (prV * dmustar/sigmastar) + S * pmustar * (epsilon)/(sigma_sq))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx17 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx18 <- (sqsq^2 * (sigma_sq) * sigmastar)
  sigx19 <- (2 * sigx16 - S^2 * sigx17 * ewu^2 * (epsilon)^2/sigx18)
  sigx20 <- (sigx15 * ewu/(sigma_sq) + sigx19 * sigx4)
  sigx21 <- (sigx20 * sigmastar + 0.5 * (prV * sigx4 * ewu/sqsq))
  sigx22 <- (sigx4 * ewu * ewv1_h * sigmastar/(ewu * ewv1_h)^2)
  sigx23 <- (sigx21 * ewv1/(ewu * ewv1_h) - 0.5 * sigx22)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  sigx25 <- (sigx4 * sigx1 * sigmastar/(ewu * ewv1_h) - depsi/ewv2_h)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx6/sigx7, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx14 * sigx1 * pwZ/sigx7, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx23 * sigx1 * pwZ/sigx7, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx24 * (1 - pwZ)/(sigx7 * ewv2_h), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx25 * dwZ/sigx7, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmnsfraynormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  sqsq <- ((sigma_sq) * sigmastar)
  mustar <- (S * ewu * (epsilon)/sqsq)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigx1 <- exp(0.5 * (-mustar)^2 - (S * (epsilon))^2/(2 * ewv1))
  sigx2 <- (sigmastar/ewv1 - ewu/sqsq)
  sigx3 <- (pmustar + S * dmustar * sigx2 * (epsilon))
  sigx4 <- (dmustar * sigmastar - S * ewu * pmustar * (epsilon)/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx5 <- (ewu * sigx3/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1)
  sigx6 <- ((1 - prZ) * sigx1 * sigx5 * sigmastar/(ewu * ewv1_h) + S * depsi *
    prZ * (epsilon)/ewv2_h^3)
  sigx7 <- ((1 - prZ) * sigx4 * sigx1 * sigmastar/(ewu * ewv1_h) + depsi * prZ/ewv2_h)
  sigx8 <- (0.5 * (dmustar * ewv1/sigmastar) - S * pmustar * (epsilon))
  sigx9 <- (prU * ewv1/sigmastar)
  sigx10 <- (1/sqsq - (0.5 * sigx9 + sigmastar) * ewu/sqsq^2)
  sigx11 <- (sigx8 * prU + S^2 * sigx10 * sigx4 * ewu * (epsilon)^2/sigmastar)
  sigx12 <- (prU * sigx4 * ewv1/sigmastar)
  sigx13 <- (sigx11 * sigmastar + 0.5 * sigx12)
  sigx14 <- (sigx13/((sigma_sq) * ewv1_h) - sigx4 * ewu * ewv1_h * sigmastar/(ewu *
    ewv1_h)^2)
  prV <- (1 - ewv1/(sigma_sq))
  sigx15 <- (0.5 * (prV * dmustar/sigmastar) + S * pmustar * (epsilon)/(sigma_sq))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx17 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx18 <- (sqsq^2 * (sigma_sq) * sigmastar)
  sigx19 <- (2 * sigx16 - S^2 * sigx17 * ewu^2 * (epsilon)^2/sigx18)
  sigx20 <- (sigx15 * ewu/(sigma_sq) + sigx19 * sigx4)
  sigx21 <- (sigx20 * sigmastar + 0.5 * (prV * sigx4 * ewu/sqsq))
  sigx22 <- (sigx4 * ewu * ewv1_h * sigmastar/(ewu * ewv1_h)^2)
  sigx23 <- (sigx21 * ewv1/(ewu * ewv1_h) - 0.5 * sigx22)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  sigx25 <- (sigx4 * sigx1 * sigmastar/(ewu * ewv1_h) - depsi/ewv2_h)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx6/sigx7, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx14 * (1 - prZ) * sigx1/sigx7, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx23 * (1 - prZ) * sigx1/sigx7, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx24 * prZ/(sigx7 * ewv2_h), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx25 * prZ * ewz/sigx7, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf rayleigh-normal distribution
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
chesszisfraynormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  dmusig <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  svsu <- (sigmastar/ewv - ewu/(ssq))
  prU <- (1 - ewu/(sigma_sq))
  prV <- (1 - ewv/(sigma_sq))
  dvs <- (dmusig * ewv/sigmastar)
  wzsq <- (wzdeno * (sigma_sq))
  ssg <- (0.5 * (prU * ewv/sigmastar) + sigmastar)
  sgsq <- (1/(ssq) - ssg * ewu/(ssq)^2)
  prvsig <- (prV * dmusig/sigmastar)
  prusig <- (prV * ewu/sigmastar)
  wzwz <- (1/(wzdeno * ewu) - ewu * ewz/(wzdeno * ewu)^2)
  euepsi <- exp(0.5 * (-(S * ewu * (epsilon)/(ssq)))^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1 <- (pmusig + S * dmusig * svsu * (epsilon))
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * (epsilon)/(sigma_sq))
  sigx3 <- (ewu * sigx1/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv)
  sigx4 <- (euepsi * sigx3 * ewz * sigmastar/(wzdeno * ewu) + S * prC * dwsr *
    (epsilon)/ewvsr^2)
  sigx5 <- (prC * dwsr + sigx2 * euepsi * ewz * sigmastar/(wzdeno * ewu))
  sigx6 <- (0.5 * dvs - S * pmusig * (epsilon))
  sigx7 <- (sigx6 * prU + S^2 * sgsq * sigx2 * ewu * (epsilon)^2/sigmastar)
  sigx8 <- (sigx7 * sigmastar + 0.5 * (prU * sigx2 * ewv/sigmastar))
  sigx9 <- (sigx8/wzsq - wzdeno * sigx2 * ewu * sigmastar/(wzdeno * ewu)^2)
  sigx10 <- (0.5 * prvsig + S * pmusig * (epsilon)/(sigma_sq))
  sigx11 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sp1 <- (0.5 * prusig + sigmastar)
  sigx12 <- (2 * sigx11 - S^2 * sp1 * ewu^2 * (epsilon)^2/((ssq)^2 * ssq))
  sigx13 <- (sigx10 * ewu/(sigma_sq) + sigx12 * sigx2)
  sigx14 <- (sigx13 * sigmastar + 0.5 * (prV * sigx2 * ewu/(ssq)))
  sigx15 <- (sigx14 * euepsi * ewv * ewz/(wzdeno * ewu) + 0.5 * (S^2 * prC * dwsr *
    (epsilon)^2/ewvsr^2))
  sigx16 <- (wzwz * sigx2 * euepsi * sigmastar - prC * dwsr/wzdeno)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * ((((S^2 * ewu * (epsilon)^2/((sigma_sq) * ewv) - 1) * svsu + ewu/(ssq)) *
    dmusig * ewu/(sigma_sq) + prU * (S * ((2 * (S * dmusig * svsu * (epsilon)) +
    3 * pmusig) * ewu/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv) * (epsilon) -
    dmusig * sigmastar)/ewv) * euepsi * ewz * sigmastar/(wzdeno * ewu) + prC *
    dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^2 - sigx4^2/sigx5)/sigx5, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (sigx9 * (S * prU * (epsilon)/ewv - sigx4/sigx5) + (((prU *
      (pmusig - 0.5 * (S * dmusig * ewu * (epsilon)/(ssq))) + S * sgsq * ewu *
      (S * ewu * sigx1 * (epsilon)/(sigma_sq) - 2 * sigx2) * (epsilon)/sigmastar) *
      sigmastar + 0.5 * (prU * ewu * ewv * sigx1/(ssq)))/wzdeno - wzdeno *
      ewu^2 * sigx1 * sigmastar/(wzdeno * ewu)^2)/(sigma_sq)) * euepsi * ewz/sigx5,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((sigx12 * sigx1 + (S * (0.5 * (prV/ewv) +
      1/(sigma_sq)) * dmusig * ewu * (epsilon)/sigmastar - pmusig)/(sigma_sq)) *
      ewu/(sigma_sq) + S * (2 * (sp1 * ewu^2/((ssq)^2 * ssq)) - 4/(2 * ewv)^2) *
      sigx2 * (epsilon)) * sigmastar + 0.5 * (prV * ewu^2 * sigx1/((sigma_sq)^2 *
      sigmastar)) + S * sigx14 * prU * (epsilon)/ewv) * euepsi * ewv * ewz/(wzdeno *
      ewu) + 0.5 * (S * prC * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 2) * (epsilon)/ewvsr^2) -
      sigx15 * sigx4/sigx5)/sigx5, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (wzwz *
    euepsi * sigx3 * sigmastar - (sigx16 * sigx4/sigx5 + S * prC * dwsr * (epsilon)/(wzdeno *
    ewvsr^2))) * ewz/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((prU * (ewu * (S^2 * sgsq * dmusig * (epsilon)^2 -
      sigx6/(sigma_sq)) - 0.5 * (0.5 * (prU * dmusig * ewv/sigmastar) + S^2 *
      sgsq * dmusig * ewu * (epsilon)^2)) + S^2 * ((sigx6 * prU * sgsq/(sigma_sq) -
      ((0.5 * (ewu/(sigma_sq)) + 1 - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv/sigmastar + (2 - 2 * (ssg^2 * ewu * (sigma_sq)/(ssq)^2)) *
        sigmastar) * sigx2/(ssq)^2) * ewu + (1 - 0.5 * prU) * sgsq * sigx2) *
      ewu * (epsilon)^2/sigmastar) * sigmastar + (0.5 * ((sigx6 * prU + S *
      ewu * pmusig * (epsilon)/(sigma_sq) - dmusig * sigmastar) * ewu/(sigma_sq) -
      0.5 * (prU * sigx2)) + 0.5 * (sigx7 * ewu/(sigma_sq))) * prU * ewv/sigmastar)/wzsq +
      ewu * (S^2 * sigx9 * sgsq * ewu * (epsilon)^2/(ssq) - ((((sigx6 * prU -
        S * pmusig * (epsilon)) * ewu/(sigma_sq) + dmusig * sigmastar) *
        sigmastar + (0.5 * (prU * ewv/(ssq)) - 2 * (wzdeno^2 * ewu * sigmastar/(wzdeno *
        ewu)^2)) * sigx2 * ewu)/(wzdeno * ewu)^2 + sigx8/wzsq^2) * wzdeno) -
      sigx9^2 * euepsi * ewz/sigx5) * euepsi * ewz/sigx5, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((((0.5 *
    (prV * dmusig) + 0.5 * ((dmusig * ewv/(sigma_sq) - S^2 * prV * sgsq * dmusig *
    ewu * (epsilon)^2/sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * prV * dmusig)))/sigmastar +
    sigx6 * prU * sigx12 + (S * (pmusig - ewu * (pmusig/(sigma_sq) + S * sgsq *
    dmusig * (epsilon))) * (epsilon) - sigx10 * ewu)/(sigma_sq))/(sigma_sq) -
    S^2 * (((0.5 * (prU * ewv/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) - 1) * ewv/(sigma_sq) +
      1 - 0.5 * (prU * prV))) * ewu/sigmastar + 2 * sp1)/((ssq)^2 * ssq) -
      (((ssq)^2 + 2 * (ssg * (sigma_sq)^2 * sigmastar)) * sigmastar + 0.5 *
        ((ssq)^2 * prU * ewv/sigmastar)) * sp1 * ewu/((ssq)^2 * ssq)^2) *
      sigx2 * ewu * (epsilon)^2) * sigmastar + (0.5 * (sigx13 * prU * ewv) +
    S^2 * sigx14 * sgsq * ewu * (epsilon)^2)/(ssq) + 0.5 * (((sigx6 * prU * prV +
    sigx2 * ewv/(sigma_sq)) * ewu/(sigma_sq) + prV * sigx2)/(ssq) - ssg * prV *
    sigx2 * ewu/(ssq)^2))/wzdeno - sigx14 * wzdeno * ewu/(wzdeno * ewu)^2) *
    ewv - sigx15 * sigx9/sigx5) * euepsi * ewz/sigx5, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((sigx7 * wzwz/(sigma_sq) - ((2 - 2 * (wzdeno^2 * ewu^2/(wzdeno * ewu)^2)) *
      ewz + 1) * sigx2/(wzdeno * ewu)^2) * sigmastar + 0.5 * (prU * wzwz *
      sigx2 * ewv/(ssq))) * ewu - sigx9 * sigx16 * ewz/sigx5) * euepsi * ewz/sigx5,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((sigx10 * sigx12 + (S * (S * sp1 * dmusig * ewu * (epsilon)/(ssq)^2 -
      2 * (pmusig/(sigma_sq))) * (epsilon) - 0.5 * prvsig)/(sigma_sq)) * ewv +
      (0.5 * (ewv * (S^2 * sp1 * dmusig * ewu^2 * (epsilon)^2/((ssq)^2 * sigmastar) -
        dmusig)/(sigma_sq) - 0.5 * (prV * dmusig)) + 0.5 * dmusig) * prV/sigmastar +
      S * pmusig * (epsilon)/(sigma_sq)) * ewu/(sigma_sq) + ((2 - 16 * (ewv^2/(2 *
      ewv)^2)) * (S * (epsilon))^2/(2 * ewv)^2 - S^2 * (sp1 * (1/((ssq)^2 *
      ssq) - (((ssq)^2 + 2 * (sp1 * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * ((ssq)^2 * prV * ewu/sigmastar)) * ewv/((ssq)^2 * ssq)^2) + (0.5 *
      (ewv/(sigma_sq)) - 0.5 * (0.5 * prV + ewv/(sigma_sq))) * prV/((ssq)^2 *
      ewv)) * ewu^2 * (epsilon)^2) * sigx2) * sigmastar + (sigx14 * sigx12 +
      (0.5 * (((0.5 * prvsig + 2 * (S * pmusig * (epsilon)/(sigma_sq))) * ewu -
        dmusig * sigmastar)/((sigma_sq)^2 * sigmastar) - sp1 * sigx2/(ssq)^2) +
        0.5 * (sigx13/(ssq))) * prV * ewu) * ewv + 0.5 * (prV * sigx2 * ewu/(ssq))) *
      euepsi * ewv * ewz/(wzdeno * ewu) + 0.5 * (S^2 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
      1) * prC * dwsr * (epsilon)^2/ewvsr^2) - sigx15^2/sigx5)/sigx5, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx14 * wzwz * euepsi * ewv - (sigx15 * sigx16/sigx5 +
      0.5 * (S^2 * prC * dwsr * (epsilon)^2/(wzdeno * ewvsr^2)))) * ewz/sigx5,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (wzwz * sigx2 * euepsi * sigmastar + (2 * (prC *
      dwsr/wzdeno^2) - (sigx16^2/sigx5 + (2 - 2 * (wzdeno * ewu^2 * ewz/(wzdeno *
      ewu)^2)) * sigx2 * euepsi * ewu * sigmastar/(wzdeno * ewu)^2)) * ewz -
      prC * dwsr/wzdeno) * ewz/sigx5, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesszisfraynormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  dmusig <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  svsu <- (sigmastar/ewv - ewu/(ssq))
  prU <- (1 - ewu/(sigma_sq))
  prV <- (1 - ewv/(sigma_sq))
  dvs <- (dmusig * ewv/sigmastar)
  ssg <- (0.5 * (prU * ewv/sigmastar) + sigmastar)
  sgsq <- (1/(ssq) - ssg * ewu/(ssq)^2)
  prvsig <- (prV * dmusig/sigmastar)
  prusig <- (prV * ewu/sigmastar)
  euepsi <- exp(0.5 * (-(S * ewu * (epsilon)/(ssq)))^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1 <- (pmusig + S * dmusig * svsu * (epsilon))
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * (epsilon)/(sigma_sq))
  sigx3 <- (ewu * sigx1/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv)
  sigx4 <- (ewz1 * euepsi * sigx3 * sigmastar/ewu + S * ewz2 * dwsr * (epsilon)/ewvsr^2)
  sigx5 <- (ewz2 * dwsr + ewz1 * sigx2 * euepsi * sigmastar/ewu)
  sigx6 <- (0.5 * dvs - S * pmusig * (epsilon))
  sigx7 <- (sigx6 * prU + S^2 * sgsq * sigx2 * ewu * (epsilon)^2/sigmastar)
  sigx8 <- (sigx7 * sigmastar + 0.5 * (prU * sigx2 * ewv/sigmastar))
  sigx9 <- (sigx8/(sigma_sq) - sigx2 * sigmastar/ewu)
  sigx10 <- (0.5 * prvsig + S * pmusig * (epsilon)/(sigma_sq))
  sigx11 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx12 <- (2 * sigx11 - S^2 * (0.5 * prusig + sigmastar) * ewu^2 * (epsilon)^2/((ssq)^2 *
    ssq))
  sigx13 <- (sigx10 * ewu/(sigma_sq) + sigx12 * sigx2)
  sigx14 <- (sigx13 * sigmastar + 0.5 * (prV * sigx2 * ewu/(ssq)))
  sigx15 <- (sigx14 * ewz1 * euepsi * ewv/ewu + 0.5 * (S^2 * ewz2 * dwsr * (epsilon)^2/ewvsr^2))
  sigx16 <- (sigx2 * euepsi * sigmastar/ewu - dwsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((S^2 * ewu * (epsilon)^2/((sigma_sq) * ewv) - 1) * svsu + ewu/(ssq)) *
      dmusig * ewu/(sigma_sq) + prU * (S * ((2 * (S * dmusig * svsu * (epsilon)) +
      3 * pmusig) * ewu/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv) * (epsilon) -
      dmusig * sigmastar)/ewv) * ewz1 * euepsi * sigmastar/ewu + ewz2 * dwsr *
      (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^2 - sigx4^2/sigx5)/sigx5, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (sigx9 * (S * prU * (epsilon)/ewv - sigx4/sigx5) + ((prU *
      (pmusig - 0.5 * (S * dmusig * ewu * (epsilon)/(ssq))) + S * sgsq * ewu *
      (S * ewu * sigx1 * (epsilon)/(sigma_sq) - 2 * sigx2) * (epsilon)/sigmastar) *
      sigmastar + (0.5 * (prU * ewu * ewv/(ssq)) - sigmastar) * sigx1)/(sigma_sq)) *
      ewz1 * euepsi/sigx5, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((sigx12 * sigx1 + (S * (0.5 * (prV/ewv) +
      1/(sigma_sq)) * dmusig * ewu * (epsilon)/sigmastar - pmusig)/(sigma_sq)) *
      ewu/(sigma_sq) + S * (2 * ((0.5 * prusig + sigmastar) * ewu^2/((ssq)^2 *
      ssq)) - 4/(2 * ewv)^2) * sigx2 * (epsilon)) * sigmastar + 0.5 * (prV *
      ewu^2 * sigx1/((sigma_sq)^2 * sigmastar)) + S * sigx14 * prU * (epsilon)/ewv) *
      ewz1 * euepsi * ewv/ewu + 0.5 * (S * ewz2 * dwsr * (S^2 * (epsilon)^2/ewvsr^2 -
      2) * (epsilon)/ewvsr^2) - sigx15 * sigx4/sigx5)/sigx5, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((euepsi *
    sigx3 * sigmastar/ewu - S * dwsr * (epsilon)/ewvsr^2)/(pi * sigx5 * ((Wz)^2 +
    1)) - pi * sigx4 * sigx16 * ((Wz)^2 + 1)/(pi * sigx5 * ((Wz)^2 + 1))^2),
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((sigx6 * prU + S * ewu * pmusig *
      (epsilon)/(sigma_sq) - dmusig * sigmastar) * ewu/(sigma_sq) - 0.5 * (prU *
      sigx2)) + 0.5 * (sigx7 * ewu/(sigma_sq)) - 0.5 * sigx2) * ewv/sigmastar -
      sigx6 * sigmastar) * prU + (prU * (ewu * (S^2 * sgsq * dmusig * (epsilon)^2 -
      sigx6/(sigma_sq)) - 0.5 * (0.5 * (prU * dmusig * ewv/sigmastar) + S^2 *
      sgsq * dmusig * ewu * (epsilon)^2)) + S^2 * ((sigx6 * prU * sgsq/(sigma_sq) -
      ((0.5 * (ewu/(sigma_sq)) + 1 - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv/sigmastar + (2 - 2 * (ssg^2 * ewu * (sigma_sq)/(ssq)^2)) *
        sigmastar) * sigx2/(ssq)^2) * ewu + (1 - 0.5 * prU) * sgsq * sigx2) *
      ewu * (epsilon)^2/sigmastar) * sigmastar + ewu * (S^2 * sigx9 * sgsq *
      ewu * (epsilon)^2/sigmastar - sigx8/(sigma_sq)))/(sigma_sq) + sigx2 *
      sigmastar/ewu - sigx9^2 * ewz1 * euepsi/sigx5) * ewz1 * euepsi/sigx5,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((((0.5 *
    (prV * dmusig) + 0.5 * ((dmusig * ewv/(sigma_sq) - S^2 * prV * sgsq * dmusig *
    ewu * (epsilon)^2/sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * prV * dmusig)))/sigmastar +
    sigx6 * prU * sigx12 + (S * (pmusig - ewu * (pmusig/(sigma_sq) + S * sgsq *
    dmusig * (epsilon))) * (epsilon) - sigx10 * ewu)/(sigma_sq))/(sigma_sq) -
    S^2 * (((0.5 * (prU * ewv/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) - 1) * ewv/(sigma_sq) +
      1 - 0.5 * (prU * prV))) * ewu/sigmastar + 2 * (0.5 * prusig + sigmastar))/((ssq)^2 *
      ssq) - (((ssq)^2 + 2 * (ssg * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * ((ssq)^2 * prU * ewv/sigmastar)) * (0.5 * prusig + sigmastar) *
      ewu/((ssq)^2 * ssq)^2) * sigx2 * ewu * (epsilon)^2) * sigmastar + (0.5 *
    (sigx13 * prU * ewv) + S^2 * sigx14 * sgsq * ewu * (epsilon)^2)/(ssq) + 0.5 *
    (((sigx6 * prU * prV + sigx2 * ewv/(sigma_sq)) * ewu/(sigma_sq) + prV * sigx2)/(ssq) -
      ssg * prV * sigx2 * ewu/(ssq)^2) - sigx14/ewu) * ewv - sigx15 * sigx9/sigx5) *
    ewz1 * euepsi/sigx5, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx9 * (1/(pi * sigx5 * ((Wz)^2 + 1)) - pi * sigx16 * ((Wz)^2 + 1) * ewz1/(pi *
    sigx5 * ((Wz)^2 + 1))^2) * euepsi, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((sigx10 * sigx12 + (S * (S * (0.5 * prusig + sigmastar) * dmusig * ewu *
      (epsilon)/(ssq)^2 - 2 * (pmusig/(sigma_sq))) * (epsilon) - 0.5 * prvsig)/(sigma_sq)) *
      ewv + (0.5 * (ewv * (S^2 * (0.5 * prusig + sigmastar) * dmusig * ewu^2 *
      (epsilon)^2/((ssq)^2 * sigmastar) - dmusig)/(sigma_sq) - 0.5 * (prV *
      dmusig)) + 0.5 * dmusig) * prV/sigmastar + S * pmusig * (epsilon)/(sigma_sq)) *
      ewu/(sigma_sq) + ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
      ewv)^2 - S^2 * ((0.5 * prusig + sigmastar) * (1/((ssq)^2 * ssq) - (((ssq)^2 +
      2 * ((0.5 * prusig + sigmastar) * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * ((ssq)^2 * prV * ewu/sigmastar)) * ewv/((ssq)^2 * ssq)^2) + (0.5 *
      (ewv/(sigma_sq)) - 0.5 * (0.5 * prV + ewv/(sigma_sq))) * prV/((ssq)^2 *
      ewv)) * ewu^2 * (epsilon)^2) * sigx2) * sigmastar + (sigx14 * sigx12 +
      (0.5 * (((0.5 * prvsig + 2 * (S * pmusig * (epsilon)/(sigma_sq))) * ewu -
        dmusig * sigmastar)/((sigma_sq)^2 * sigmastar) - (0.5 * prusig +
        sigmastar) * sigx2/(ssq)^2) + 0.5 * (sigx13/(ssq))) * prV * ewu) *
      ewv + 0.5 * (prV * sigx2 * ewu/(ssq))) * ewz1 * euepsi * ewv/ewu + 0.5 *
      (S^2 * ewz2 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) * dwsr * (epsilon)^2/ewvsr^2) -
      sigx15^2/sigx5)/sigx5, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx14 * euepsi * ewv/ewu - 0.5 * (S^2 * dwsr *
      (epsilon)^2/ewvsr^2))/(pi * sigx5 * ((Wz)^2 + 1)) - pi * sigx15 * sigx16 *
      ((Wz)^2 + 1)/(pi * sigx5 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (sigx16 * (sigx2 * euepsi * sigmastar/ewu +
      2 * (pi * Wz * sigx5) - dwsr)/(pi * sigx5 * ((Wz)^2 + 1))^2), FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesszisfraynormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  dmusig <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  svsu <- (sigmastar/ewv - ewu/(ssq))
  prU <- (1 - ewu/(sigma_sq))
  prV <- (1 - ewv/(sigma_sq))
  dvs <- (dmusig * ewv/sigmastar)
  ssg <- (0.5 * (prU * ewv/sigmastar) + sigmastar)
  sgsq <- (1/(ssq) - ssg * ewu/(ssq)^2)
  prvsig <- (prV * dmusig/sigmastar)
  prusig <- (prV * ewu/sigmastar)
  euepsi <- exp(0.5 * (-(S * ewu * (epsilon)/(ssq)))^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1 <- (pmusig + S * dmusig * svsu * (epsilon))
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * (epsilon)/(sigma_sq))
  sigx3 <- (ewu * sigx1/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv)
  sigx4 <- (euepsi * sigx3 * pwZ * sigmastar/ewu + S * (1 - pwZ) * dwsr * (epsilon)/ewvsr^2)
  sigx5 <- ((1 - pwZ) * dwsr + sigx2 * euepsi * pwZ * sigmastar/ewu)
  sigx6 <- (0.5 * dvs - S * pmusig * (epsilon))
  sigx7 <- (sigx6 * prU + S^2 * sgsq * sigx2 * ewu * (epsilon)^2/sigmastar)
  sigx8 <- (sigx7 * sigmastar + 0.5 * (prU * sigx2 * ewv/sigmastar))
  sigx9 <- (sigx8/(sigma_sq) - sigx2 * sigmastar/ewu)
  sigx10 <- (0.5 * prvsig + S * pmusig * (epsilon)/(sigma_sq))
  sigx11 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx12 <- (2 * sigx11 - S^2 * (0.5 * prusig + sigmastar) * ewu^2 * (epsilon)^2/((ssq)^2 *
    ssq))
  sigx13 <- (sigx10 * ewu/(sigma_sq) + sigx12 * sigx2)
  sigx14 <- (sigx13 * sigmastar + 0.5 * (prV * sigx2 * ewu/(ssq)))
  sigx15 <- (sigx14 * euepsi * ewv * pwZ/ewu + 0.5 * (S^2 * (1 - pwZ) * dwsr *
    (epsilon)^2/ewvsr^2))
  sigx16 <- (sigx2 * euepsi * sigmastar/ewu - dwsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((S^2 * ewu * (epsilon)^2/((sigma_sq) * ewv) - 1) * svsu + ewu/(ssq)) *
      dmusig * ewu/(sigma_sq) + prU * (S * ((2 * (S * dmusig * svsu * (epsilon)) +
      3 * pmusig) * ewu/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv) * (epsilon) -
      dmusig * sigmastar)/ewv) * euepsi * pwZ * sigmastar/ewu + (1 - pwZ) *
      dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^2 - sigx4^2/sigx5)/sigx5,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (sigx9 * (S * prU * (epsilon)/ewv - sigx4/sigx5) + ((prU *
      (pmusig - 0.5 * (S * dmusig * ewu * (epsilon)/(ssq))) + S * sgsq * ewu *
      (S * ewu * sigx1 * (epsilon)/(sigma_sq) - 2 * sigx2) * (epsilon)/sigmastar) *
      sigmastar + (0.5 * (prU * ewu * ewv/(ssq)) - sigmastar) * sigx1)/(sigma_sq)) *
      euepsi * pwZ/sigx5, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((sigx12 * sigx1 + (S * (0.5 * (prV/ewv) +
      1/(sigma_sq)) * dmusig * ewu * (epsilon)/sigmastar - pmusig)/(sigma_sq)) *
      ewu/(sigma_sq) + S * (2 * ((0.5 * prusig + sigmastar) * ewu^2/((ssq)^2 *
      ssq)) - 4/(2 * ewv)^2) * sigx2 * (epsilon)) * sigmastar + 0.5 * (prV *
      ewu^2 * sigx1/((sigma_sq)^2 * sigmastar)) + S * sigx14 * prU * (epsilon)/ewv) *
      euepsi * ewv * pwZ/ewu + 0.5 * (S * (1 - pwZ) * dwsr * (S^2 * (epsilon)^2/ewvsr^2 -
      2) * (epsilon)/ewvsr^2) - sigx15 * sigx4/sigx5)/sigx5, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * dwZ * (euepsi *
    sigx3 * sigmastar/ewu - (sigx16 * sigx4/sigx5 + S * dwsr * (epsilon)/ewvsr^2))/sigx5,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((sigx6 * prU + S * ewu * pmusig *
      (epsilon)/(sigma_sq) - dmusig * sigmastar) * ewu/(sigma_sq) - 0.5 * (prU *
      sigx2)) + 0.5 * (sigx7 * ewu/(sigma_sq)) - 0.5 * sigx2) * ewv/sigmastar -
      sigx6 * sigmastar) * prU + (prU * (ewu * (S^2 * sgsq * dmusig * (epsilon)^2 -
      sigx6/(sigma_sq)) - 0.5 * (0.5 * (prU * dmusig * ewv/sigmastar) + S^2 *
      sgsq * dmusig * ewu * (epsilon)^2)) + S^2 * ((sigx6 * prU * sgsq/(sigma_sq) -
      ((0.5 * (ewu/(sigma_sq)) + 1 - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv/sigmastar + (2 - 2 * (ssg^2 * ewu * (sigma_sq)/(ssq)^2)) *
        sigmastar) * sigx2/(ssq)^2) * ewu + (1 - 0.5 * prU) * sgsq * sigx2) *
      ewu * (epsilon)^2/sigmastar) * sigmastar + ewu * (S^2 * sigx9 * sgsq *
      ewu * (epsilon)^2/sigmastar - sigx8/(sigma_sq)))/(sigma_sq) + sigx2 *
      sigmastar/ewu - sigx9^2 * euepsi * pwZ/sigx5) * euepsi * pwZ/sigx5, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((((0.5 *
    (prV * dmusig) + 0.5 * ((dmusig * ewv/(sigma_sq) - S^2 * prV * sgsq * dmusig *
    ewu * (epsilon)^2/sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * prV * dmusig)))/sigmastar +
    sigx6 * prU * sigx12 + (S * (pmusig - ewu * (pmusig/(sigma_sq) + S * sgsq *
    dmusig * (epsilon))) * (epsilon) - sigx10 * ewu)/(sigma_sq))/(sigma_sq) -
    S^2 * (((0.5 * (prU * ewv/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) - 1) * ewv/(sigma_sq) +
      1 - 0.5 * (prU * prV))) * ewu/sigmastar + 2 * (0.5 * prusig + sigmastar))/((ssq)^2 *
      ssq) - (((ssq)^2 + 2 * (ssg * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * ((ssq)^2 * prU * ewv/sigmastar)) * (0.5 * prusig + sigmastar) *
      ewu/((ssq)^2 * ssq)^2) * sigx2 * ewu * (epsilon)^2) * sigmastar + (0.5 *
    (sigx13 * prU * ewv) + S^2 * sigx14 * sgsq * ewu * (epsilon)^2)/(ssq) + 0.5 *
    (((sigx6 * prU * prV + sigx2 * ewv/(sigma_sq)) * ewu/(sigma_sq) + prV * sigx2)/(ssq) -
      ssg * prV * sigx2 * ewu/(ssq)^2) - sigx14/ewu) * ewv - sigx15 * sigx9/sigx5) *
    euepsi * pwZ/sigx5, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx9 * (1 - sigx16 * pwZ/sigx5) * dwZ * euepsi/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((sigx10 * sigx12 + (S * (S * (0.5 * prusig + sigmastar) * dmusig * ewu *
      (epsilon)/(ssq)^2 - 2 * (pmusig/(sigma_sq))) * (epsilon) - 0.5 * prvsig)/(sigma_sq)) *
      ewv + (0.5 * (ewv * (S^2 * (0.5 * prusig + sigmastar) * dmusig * ewu^2 *
      (epsilon)^2/((ssq)^2 * sigmastar) - dmusig)/(sigma_sq) - 0.5 * (prV *
      dmusig)) + 0.5 * dmusig) * prV/sigmastar + S * pmusig * (epsilon)/(sigma_sq)) *
      ewu/(sigma_sq) + ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
      ewv)^2 - S^2 * ((0.5 * prusig + sigmastar) * (1/((ssq)^2 * ssq) - (((ssq)^2 +
      2 * ((0.5 * prusig + sigmastar) * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * ((ssq)^2 * prV * ewu/sigmastar)) * ewv/((ssq)^2 * ssq)^2) + (0.5 *
      (ewv/(sigma_sq)) - 0.5 * (0.5 * prV + ewv/(sigma_sq))) * prV/((ssq)^2 *
      ewv)) * ewu^2 * (epsilon)^2) * sigx2) * sigmastar + (sigx14 * sigx12 +
      (0.5 * (((0.5 * prvsig + 2 * (S * pmusig * (epsilon)/(sigma_sq))) * ewu -
        dmusig * sigmastar)/((sigma_sq)^2 * sigmastar) - (0.5 * prusig +
        sigmastar) * sigx2/(ssq)^2) + 0.5 * (sigx13/(ssq))) * prV * ewu) *
      ewv + 0.5 * (prV * sigx2 * ewu/(ssq))) * euepsi * ewv * pwZ/ewu + 0.5 *
      (S^2 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) * (1 - pwZ) * dwsr * (epsilon)^2/ewvsr^2) -
      sigx15^2/sigx5)/sigx5, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx14 * euepsi * ewv/ewu - (sigx15 * sigx16/sigx5 +
      0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2))) * dwZ/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx16 * dwZ/sigx5 + Wz) * sigx16 * dwZ/sigx5),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesszisfraynormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  dmusig <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  svsu <- (sigmastar/ewv - ewu/(ssq))
  prU <- (1 - ewu/(sigma_sq))
  prV <- (1 - ewv/(sigma_sq))
  dvs <- (dmusig * ewv/sigmastar)
  ssg <- (0.5 * (prU * ewv/sigmastar) + sigmastar)
  sgsq <- (1/(ssq) - ssg * ewu/(ssq)^2)
  prvsig <- (prV * dmusig/sigmastar)
  prusig <- (prV * ewu/sigmastar)
  euepsi <- exp(0.5 * (-(S * ewu * (epsilon)/(ssq)))^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1 <- (pmusig + S * dmusig * svsu * (epsilon))
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * (epsilon)/(sigma_sq))
  sigx3 <- (ewu * sigx1/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv)
  sigx4 <- ((1 - prZ) * euepsi * sigx3 * sigmastar/ewu + S * dwsr * prZ * (epsilon)/ewvsr^2)
  sigx5 <- ((1 - prZ) * sigx2 * euepsi * sigmastar/ewu + dwsr * prZ)
  sigx6 <- (0.5 * dvs - S * pmusig * (epsilon))
  sigx7 <- (sigx6 * prU + S^2 * sgsq * sigx2 * ewu * (epsilon)^2/sigmastar)
  sigx8 <- (sigx7 * sigmastar + 0.5 * (prU * sigx2 * ewv/sigmastar))
  sigx9 <- (sigx8/(sigma_sq) - sigx2 * sigmastar/ewu)
  sigx10 <- (0.5 * prvsig + S * pmusig * (epsilon)/(sigma_sq))
  sigx11 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx12 <- (2 * sigx11 - S^2 * (0.5 * prusig + sigmastar) * ewu^2 * (epsilon)^2/((ssq)^2 *
    ssq))
  sigx13 <- (sigx10 * ewu/(sigma_sq) + sigx12 * sigx2)
  sigx14 <- (sigx13 * sigmastar + 0.5 * (prV * sigx2 * ewu/(ssq)))
  sigx15 <- (sigx14 * (1 - prZ) * euepsi * ewv/ewu + 0.5 * (S^2 * dwsr * prZ *
    (epsilon)^2/ewvsr^2))
  sigx16 <- (sigx2 * euepsi * sigmastar/ewu - dwsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((S^2 * ewu * (epsilon)^2/((sigma_sq) * ewv) - 1) * svsu + ewu/(ssq)) *
      dmusig * ewu/(sigma_sq) + prU * (S * ((2 * (S * dmusig * svsu * (epsilon)) +
      3 * pmusig) * ewu/(sigma_sq) + S * prU * sigx2 * (epsilon)/ewv) * (epsilon) -
      dmusig * sigmastar)/ewv) * (1 - prZ) * euepsi * sigmastar/ewu + dwsr *
      prZ * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^2 - sigx4^2/sigx5)/sigx5,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (sigx9 * (S * prU * (epsilon)/ewv - sigx4/sigx5) + ((prU *
      (pmusig - 0.5 * (S * dmusig * ewu * (epsilon)/(ssq))) + S * sgsq * ewu *
      (S * ewu * sigx1 * (epsilon)/(sigma_sq) - 2 * sigx2) * (epsilon)/sigmastar) *
      sigmastar + (0.5 * (prU * ewu * ewv/(ssq)) - sigmastar) * sigx1)/(sigma_sq)) *
      (1 - prZ) * euepsi/sigx5, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((sigx12 * sigx1 + (S * (0.5 * (prV/ewv) +
      1/(sigma_sq)) * dmusig * ewu * (epsilon)/sigmastar - pmusig)/(sigma_sq)) *
      ewu/(sigma_sq) + S * (2 * ((0.5 * prusig + sigmastar) * ewu^2/((ssq)^2 *
      ssq)) - 4/(2 * ewv)^2) * sigx2 * (epsilon)) * sigmastar + 0.5 * (prV *
      ewu^2 * sigx1/((sigma_sq)^2 * sigmastar)) + S * sigx14 * prU * (epsilon)/ewv) *
      (1 - prZ) * euepsi * ewv/ewu + 0.5 * (S * dwsr * prZ * (S^2 * (epsilon)^2/ewvsr^2 -
      2) * (epsilon)/ewvsr^2) - sigx15 * sigx4/sigx5)/sigx5, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * prZ * (euepsi *
    sigx3 * sigmastar/ewu - (sigx4 * sigx16/sigx5 + S * dwsr * (epsilon)/ewvsr^2)) *
    ewz/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((sigx6 * prU + S * ewu * pmusig *
      (epsilon)/(sigma_sq) - dmusig * sigmastar) * ewu/(sigma_sq) - 0.5 * (prU *
      sigx2)) + 0.5 * (sigx7 * ewu/(sigma_sq)) - 0.5 * sigx2) * ewv/sigmastar -
      sigx6 * sigmastar) * prU + (prU * (ewu * (S^2 * sgsq * dmusig * (epsilon)^2 -
      sigx6/(sigma_sq)) - 0.5 * (0.5 * (prU * dmusig * ewv/sigmastar) + S^2 *
      sgsq * dmusig * ewu * (epsilon)^2)) + S^2 * ((sigx6 * prU * sgsq/(sigma_sq) -
      ((0.5 * (ewu/(sigma_sq)) + 1 - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv/sigmastar + (2 - 2 * (ssg^2 * ewu * (sigma_sq)/(ssq)^2)) *
        sigmastar) * sigx2/(ssq)^2) * ewu + (1 - 0.5 * prU) * sgsq * sigx2) *
      ewu * (epsilon)^2/sigmastar) * sigmastar + ewu * (S^2 * sigx9 * sgsq *
      ewu * (epsilon)^2/sigmastar - sigx8/(sigma_sq)))/(sigma_sq) + sigx2 *
      sigmastar/ewu - sigx9^2 * (1 - prZ) * euepsi/sigx5) * (1 - prZ) * euepsi/sigx5,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((((0.5 *
    (prV * dmusig) + 0.5 * ((dmusig * ewv/(sigma_sq) - S^2 * prV * sgsq * dmusig *
    ewu * (epsilon)^2/sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * prV * dmusig)))/sigmastar +
    sigx6 * prU * sigx12 + (S * (pmusig - ewu * (pmusig/(sigma_sq) + S * sgsq *
    dmusig * (epsilon))) * (epsilon) - sigx10 * ewu)/(sigma_sq))/(sigma_sq) -
    S^2 * (((0.5 * (prU * ewv/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) - 1) * ewv/(sigma_sq) +
      1 - 0.5 * (prU * prV))) * ewu/sigmastar + 2 * (0.5 * prusig + sigmastar))/((ssq)^2 *
      ssq) - (((ssq)^2 + 2 * (ssg * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * ((ssq)^2 * prU * ewv/sigmastar)) * (0.5 * prusig + sigmastar) *
      ewu/((ssq)^2 * ssq)^2) * sigx2 * ewu * (epsilon)^2) * sigmastar + (0.5 *
    (sigx13 * prU * ewv) + S^2 * sigx14 * sgsq * ewu * (epsilon)^2)/(ssq) + 0.5 *
    (((sigx6 * prU * prV + sigx2 * ewv/(sigma_sq)) * ewu/(sigma_sq) + prV * sigx2)/(ssq) -
      ssg * prV * sigx2 * ewu/(ssq)^2) - sigx14/ewu) * ewv - sigx15 * sigx9/sigx5) *
    (1 - prZ) * euepsi/sigx5, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx9 * (1 - sigx16 * (1 - prZ)/sigx5) * prZ * euepsi * ewz/sigx5, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((sigx10 * sigx12 + (S * (S * (0.5 * prusig + sigmastar) * dmusig * ewu *
      (epsilon)/(ssq)^2 - 2 * (pmusig/(sigma_sq))) * (epsilon) - 0.5 * prvsig)/(sigma_sq)) *
      ewv + (0.5 * (ewv * (S^2 * (0.5 * prusig + sigmastar) * dmusig * ewu^2 *
      (epsilon)^2/((ssq)^2 * sigmastar) - dmusig)/(sigma_sq) - 0.5 * (prV *
      dmusig)) + 0.5 * dmusig) * prV/sigmastar + S * pmusig * (epsilon)/(sigma_sq)) *
      ewu/(sigma_sq) + ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
      ewv)^2 - S^2 * ((0.5 * prusig + sigmastar) * (1/((ssq)^2 * ssq) - (((ssq)^2 +
      2 * ((0.5 * prusig + sigmastar) * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * ((ssq)^2 * prV * ewu/sigmastar)) * ewv/((ssq)^2 * ssq)^2) + (0.5 *
      (ewv/(sigma_sq)) - 0.5 * (0.5 * prV + ewv/(sigma_sq))) * prV/((ssq)^2 *
      ewv)) * ewu^2 * (epsilon)^2) * sigx2) * sigmastar + (sigx14 * sigx12 +
      (0.5 * (((0.5 * prvsig + 2 * (S * pmusig * (epsilon)/(sigma_sq))) * ewu -
        dmusig * sigmastar)/((sigma_sq)^2 * sigmastar) - (0.5 * prusig +
        sigmastar) * sigx2/(ssq)^2) + 0.5 * (sigx13/(ssq))) * prV * ewu) *
      ewv + 0.5 * (prV * sigx2 * ewu/(ssq))) * (1 - prZ) * euepsi * ewv/ewu +
      0.5 * (S^2 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) * dwsr * prZ * (epsilon)^2/ewvsr^2) -
      sigx15^2/sigx5)/sigx5, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx14 * euepsi * ewv/ewu - (sigx15 * sigx16/sigx5 +
      0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2))) * prZ * ewz/sigx5, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * sigx16 * (1 - (sigx16 * prZ/sigx5 + 1) * ewz) *
      prZ * ewz/sigx5, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Different sigma_v

## logit specification class membership
chessmnsfraynormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  sqsq <- ((sigma_sq) * sigmastar)
  mustar <- (S * ewu * (epsilon)/sqsq)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wzuv <- (wzdeno * ewu * ewv1_h)
  sigx1 <- exp(0.5 * (-mustar)^2 - (S * (epsilon))^2/(2 * ewv1))
  sigx2 <- (sigmastar/ewv1 - ewu/sqsq)
  sigx3 <- (pmustar + S * dmustar * sigx2 * (epsilon))
  sigx4 <- (dmustar * sigmastar - S * ewu * pmustar * (epsilon)/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx5 <- (ewu * sigx3/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1)
  sigx6 <- (sigx1 * sigx5 * ewz * sigmastar/wzuv + S * prC * depsi * (epsilon)/ewv2_h^3)
  sigx7 <- (prC * depsi/ewv2_h + sigx4 * sigx1 * ewz * sigmastar/wzuv)
  sigx8 <- (0.5 * (dmustar * ewv1/sigmastar) - S * pmustar * (epsilon))
  sigx9 <- (prU * ewv1/sigmastar)
  sigx10 <- (1/sqsq - (0.5 * sigx9 + sigmastar) * ewu/sqsq^2)
  sigx11 <- (sigx8 * prU + S^2 * sigx10 * sigx4 * ewu * (epsilon)^2/sigmastar)
  sigx12 <- (prU * sigx4 * ewv1/sigmastar)
  sigx13 <- (sigx11 * sigmastar + 0.5 * sigx12)
  wsqv <- (wzdeno * (sigma_sq) * ewv1_h)
  sigx14 <- (sigx13/wsqv - wzdeno * sigx4 * ewu * ewv1_h * sigmastar/wzuv^2)
  prV <- (1 - ewv1/(sigma_sq))
  sigx15 <- (0.5 * (prV * dmustar/sigmastar) + S * pmustar * (epsilon)/(sigma_sq))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx17 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx18 <- (sqsq^2 * (sigma_sq) * sigmastar)
  sigx19 <- (2 * sigx16 - S^2 * sigx17 * ewu^2 * (epsilon)^2/sigx18)
  sigx20 <- (sigx15 * ewu/(sigma_sq) + sigx19 * sigx4)
  sigx21 <- (sigx20 * sigmastar + 0.5 * (prV * sigx4 * ewu/sqsq))
  sigx22 <- (wzdeno * sigx4 * ewu * ewv1_h * sigmastar/wzuv^2)
  sigx23 <- (sigx21 * ewv1/wzuv - 0.5 * sigx22)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  sigx25 <- ((1/wzuv - ewu * ewv1_h * ewz/wzuv^2) * sigx4 * sigx1 * sigmastar -
    prC * depsi/(wzdeno * ewv2_h))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * ((((S^2 * ewu * (epsilon)^2/((sigma_sq) * ewv1) - 1) * sigx2 + ewu/sqsq) *
    dmustar * ewu/(sigma_sq) + prU * (S * ((2 * (S * dmustar * sigx2 * (epsilon)) +
    3 * pmustar) * ewu/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1) * (epsilon) -
    dmustar * sigmastar)/ewv1) * sigx1 * ewz * sigmastar/wzuv + prC * depsi *
    (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 - sigx6^2/sigx7)/sigx7, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (sigx14 * (S * prU * (epsilon)/ewv1 - sigx6/sigx7) +
      (((prU * (pmustar - 0.5 * (S * dmustar * ewu * (epsilon)/sqsq)) + S *
        sigx10 * ewu * (S * ewu * sigx3 * (epsilon)/(sigma_sq) - 2 * sigx4) *
        (epsilon)/sigmastar) * sigmastar + 0.5 * (prU * ewu * ewv1 * sigx3/sqsq))/(wzdeno *
        ewv1_h) - wzdeno * ewu^2 * ewv1_h * sigx3 * sigmastar/wzuv^2)/(sigma_sq)) *
      sigx1 * ewz/sigx7, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx23 * (S * prU * (epsilon)/ewv1 - sigx6/sigx7) +
      (((sigx19 * sigx3 + (S * (0.5 * (prV/ewv1) + 1/(sigma_sq)) * dmustar *
        ewu * (epsilon)/sigmastar - pmustar)/(sigma_sq)) * ewu/(sigma_sq) +
        S * (2 * (sigx17 * ewu^2/sigx18) - 4/(2 * ewv1)^2) * sigx4 * (epsilon)) *
        sigmastar + 0.5 * (prV * ewu^2 * sigx3/((sigma_sq)^2 * sigmastar))) *
        ewv1/wzuv - 0.5 * (wzdeno * ewu^2 * ewv1_h * sigx3 * sigmastar/(wzuv^2 *
      (sigma_sq)))) * sigx1 * ewz/sigx7, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * prC * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx7 * ewv2_h^3) - sigx24 * sigx6 * ewv2_h/(sigx7 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((1/wzuv -
    ewu * ewv1_h * ewz/wzuv^2) * sigx1 * sigx5 * sigmastar - (sigx25 * sigx6/sigx7 +
    S * prC * depsi * (epsilon)/(wzdeno * ewv2_h^3))) * ewz/sigx7, FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((prU * (ewu * (S^2 * sigx10 * dmustar * (epsilon)^2 -
      sigx8/(sigma_sq)) - 0.5 * (0.5 * (prU * dmustar * ewv1/sigmastar) + S^2 *
      sigx10 * dmustar * ewu * (epsilon)^2)) + S^2 * ((sigx8 * prU * sigx10/(sigma_sq) -
      ((0.5 * (ewu/(sigma_sq)) + 1 - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv1/sigmastar + (2 - 2 * ((0.5 * sigx9 + sigmastar)^2 * ewu *
        (sigma_sq)/sqsq^2)) * sigmastar) * sigx4/sqsq^2) * ewu + (1 - 0.5 *
      prU) * sigx10 * sigx4) * ewu * (epsilon)^2/sigmastar) * sigmastar + (0.5 *
      ((sigx8 * prU + S * ewu * pmustar * (epsilon)/(sigma_sq) - dmustar *
        sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * sigx4)) + 0.5 * (sigx11 *
      ewu/(sigma_sq))) * prU * ewv1/sigmastar)/wsqv + ewu * (S^2 * sigx14 *
      sigx10 * ewu * (epsilon)^2/sqsq - ((((sigx8 * prU - S * pmustar * (epsilon)) *
      ewu/(sigma_sq) + dmustar * sigmastar) * sigmastar + (0.5 * (prU * ewv1/sqsq) -
      2 * (wzdeno^2 * ewu * ewv1_h^2 * sigmastar/wzuv^2)) * sigx4 * ewu)/wzuv^2 +
      sigx13/wsqv^2) * wzdeno * ewv1_h) - sigx14^2 * sigx1 * ewz/sigx7) * sigx1 *
      ewz/sigx7, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((((0.5 *
    (prV * dmustar) + 0.5 * ((dmustar * ewv1/(sigma_sq) - S^2 * prV * sigx10 *
    dmustar * ewu * (epsilon)^2/sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * prV *
    dmustar)))/sigmastar + sigx8 * prU * sigx19 + (S * (pmustar - ewu * (pmustar/(sigma_sq) +
    S * sigx10 * dmustar * (epsilon))) * (epsilon) - sigx15 * ewu)/(sigma_sq))/(sigma_sq) -
    S^2 * (((0.5 * (prU * ewv1/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) - 1) * ewv1/(sigma_sq) +
      1 - 0.5 * (prU * prV))) * ewu/sigmastar + 2 * sigx17)/sigx18 - ((sqsq^2 +
      2 * ((0.5 * sigx9 + sigmastar) * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * (sqsq^2 * prU * ewv1/sigmastar)) * sigx17 * ewu/sigx18^2) * sigx4 *
      ewu * (epsilon)^2) * sigmastar + 0.5 * (((sigx8 * prU * prV + sigx4 *
    ewv1/(sigma_sq)) * ewu/(sigma_sq) + prV * sigx4)/sqsq - (0.5 * sigx9 + sigmastar) *
    prV * sigx4 * ewu/sqsq^2) + 0.5 * (sigx20 * prU * ewv1/sqsq))/(wzdeno * ewv1_h) -
    sigx21 * wzdeno * ewu * ewv1_h/wzuv^2) * ewv1 + ewu * (S^2 * sigx23 * sigx10 *
    ewu * (epsilon)^2/sqsq - 0.5 * ((((sigx8 * prU - S * pmustar * (epsilon)) *
    ewu/(sigma_sq) + dmustar * sigmastar) * sigmastar + (0.5 * (prU * ewv1/sqsq) -
    2 * (wzdeno^2 * ewu * ewv1_h^2 * sigmastar/wzuv^2)) * sigx4 * ewu) * wzdeno *
    ewv1_h/wzuv^2)) - sigx23 * sigx14 * sigx1 * ewz/sigx7) * sigx1 * ewz/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx14 * sigx24 * prC * sigx1 * ewv2_h * ewz/(sigx7 * ewv2_h)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((sigx11 * (1/wzuv - ewu * ewv1_h * ewz/wzuv^2)/(sigma_sq) - ((2 - 2 * (wzdeno^2 *
      ewu^2 * ewv1_h^2/wzuv^2)) * ewz + 1) * sigx4 * ewv1_h/wzuv^2) * sigmastar +
      0.5 * (prU * (1/wzuv - ewu * ewv1_h * ewz/wzuv^2) * sigx4 * ewv1/sqsq)) *
      ewu - sigx14 * sigx25 * ewz/sigx7) * sigx1 * ewz/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx15 * sigx19 + (S * (S * sigx17 * dmustar * ewu * (epsilon)/sqsq^2 -
      2 * (pmustar/(sigma_sq))) * (epsilon) - 0.5 * (prV * dmustar/sigmastar))/(sigma_sq)) *
      ewv1 + (0.5 * (ewv1 * (S^2 * sigx17 * dmustar * ewu^2 * (epsilon)^2/(sqsq^2 *
      sigmastar) - dmustar)/(sigma_sq) - 0.5 * (prV * dmustar)) + 0.5 * dmustar) *
      prV/sigmastar + S * pmustar * (epsilon)/(sigma_sq)) * ewu/(sigma_sq) +
      ((2 - 16 * (ewv1^2/(2 * ewv1)^2)) * (S * (epsilon))^2/(2 * ewv1)^2 -
        S^2 * (sigx17 * (1/sigx18 - ((sqsq^2 + 2 * (sigx17 * (sigma_sq)^2 *
          sigmastar)) * sigmastar + 0.5 * (sqsq^2 * prV * ewu/sigmastar)) *
          ewv1/sigx18^2) + (0.5 * (ewv1/(sigma_sq)) - 0.5 * (0.5 * prV +
          ewv1/(sigma_sq))) * prV/(sqsq^2 * ewv1)) * ewu^2 * (epsilon)^2) *
        sigx4) * sigmastar + ((0.5 * (((0.5 * (prV * dmustar/sigmastar) +
      2 * (S * pmustar * (epsilon)/(sigma_sq))) * ewu - dmustar * sigmastar)/((sigma_sq)^2 *
      sigmastar) - sigx17 * sigx4/sqsq^2) + 0.5 * (sigx20/sqsq)) * ewv1 + 0.5 *
      (sigx4/sqsq)) * prV * ewu)/wzuv + sigx23 * sigx19 - 0.5 * (sigx21 * wzdeno *
      ewu * ewv1_h/wzuv^2)) * ewv1 - (sigx23^2 * sigx1 * ewz/sigx7 + 0.5 *
      (((sigx15 * ewu * ewv1/(sigma_sq) + 0.5 * sigx4) * sigmastar + (0.5 *
        (prV * ewv1/sqsq) - wzdeno^2 * ewu * ewv1_h^2 * sigmastar/wzuv^2) *
        sigx4 * ewu) * wzdeno * ewu * ewv1_h/wzuv^2))) * sigx1 * ewz/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx23 * sigx24 * prC * sigx1 * ewv2_h * ewz/(sigx7 * ewv2_h)^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx15 * (1/wzuv - ewu * ewv1_h * ewz/wzuv^2) *
      ewv1/(sigma_sq) - ((0.5 - wzdeno^2 * ewu^2 * ewv1_h^2/wzuv^2) * ewz +
      0.5 * wzdeno) * sigx4 * ewv1_h/wzuv^2) * ewu + (1/wzuv - ewu * ewv1_h *
      ewz/wzuv^2) * sigx19 * sigx4 * ewv1) * sigmastar + 0.5 * (prV * (1/wzuv -
      ewu * ewv1_h * ewz/wzuv^2) * sigx4 * ewu * ewv1/sqsq) - sigx23 * sigx25 *
      ewz/sigx7) * sigx1 * ewz/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx7 * ewv2_h^3) - (sigx24 * prC +
      0.5 * (sigx7 * ewv2_h)) * sigx24/(sigx7 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx25 * sigx24/sigx7 + 0.5 * (S^2 * depsi *
      (epsilon)^2/(wzdeno * ewv2_h^2)))/ewv2_h - 0.5 * (wzdeno * depsi * ewv2_h/(wzdeno *
      ewv2_h)^2)) * prC * ewz/sigx7), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv2_h) + ewv2_h/(wzdeno *
      ewv2_h)^2) * depsi - (sigx25^2/sigx7 + (2 - 2 * (wzdeno * ewu^2 * ewv1_h^2 *
      ewz/wzuv^2)) * sigx4 * sigx1 * ewu * ewv1_h * sigmastar/wzuv^2)) * ewz +
      (1/wzuv - ewu * ewv1_h * ewz/wzuv^2) * sigx4 * sigx1 * sigmastar - prC *
      depsi/(wzdeno * ewv2_h)) * ewz/sigx7, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmnsfraynormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  sqsq <- ((sigma_sq) * sigmastar)
  mustar <- (S * ewu * (epsilon)/sqsq)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigx1 <- exp(0.5 * (-mustar)^2 - (S * (epsilon))^2/(2 * ewv1))
  sigx2 <- (sigmastar/ewv1 - ewu/sqsq)
  sigx3 <- (pmustar + S * dmustar * sigx2 * (epsilon))
  sigx4 <- (dmustar * sigmastar - S * ewu * pmustar * (epsilon)/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx5 <- (ewu * sigx3/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1)
  sigx6 <- (ewz1 * sigx1 * sigx5 * sigmastar/(ewu * ewv1_h) + S * ewz2 * depsi *
    (epsilon)/ewv2_h^3)
  sigx7 <- (ewz2 * depsi/ewv2_h + ewz1 * sigx4 * sigx1 * sigmastar/(ewu * ewv1_h))
  sigx8 <- (0.5 * (dmustar * ewv1/sigmastar) - S * pmustar * (epsilon))
  sigx9 <- (prU * ewv1/sigmastar)
  sigx10 <- (1/sqsq - (0.5 * sigx9 + sigmastar) * ewu/sqsq^2)
  sigx11 <- (sigx8 * prU + S^2 * sigx10 * sigx4 * ewu * (epsilon)^2/sigmastar)
  sigx12 <- (prU * sigx4 * ewv1/sigmastar)
  sigx13 <- (sigx11 * sigmastar + 0.5 * sigx12)
  sigx14 <- (sigx13/((sigma_sq) * ewv1_h) - sigx4 * ewu * ewv1_h * sigmastar/(ewu *
    ewv1_h)^2)
  prV <- (1 - ewv1/(sigma_sq))
  sigx15 <- (0.5 * (prV * dmustar/sigmastar) + S * pmustar * (epsilon)/(sigma_sq))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx17 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx18 <- (sqsq^2 * (sigma_sq) * sigmastar)
  sigx19 <- (2 * sigx16 - S^2 * sigx17 * ewu^2 * (epsilon)^2/sigx18)
  sigx20 <- (sigx15 * ewu/(sigma_sq) + sigx19 * sigx4)
  sigx21 <- (sigx20 * sigmastar + 0.5 * (prV * sigx4 * ewu/sqsq))
  sigx22 <- (sigx4 * ewu * ewv1_h * sigmastar/(ewu * ewv1_h)^2)
  sigx23 <- (sigx21 * ewv1/(ewu * ewv1_h) - 0.5 * sigx22)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  sigx25 <- (sigx4 * sigx1 * sigmastar/(ewu * ewv1_h) - depsi/ewv2_h)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((S^2 * ewu * (epsilon)^2/((sigma_sq) * ewv1) - 1) * sigx2 + ewu/sqsq) *
      dmustar * ewu/(sigma_sq) + prU * (S * ((2 * (S * dmustar * sigx2 * (epsilon)) +
      3 * pmustar) * ewu/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1) * (epsilon) -
      dmustar * sigmastar)/ewv1) * ewz1 * sigx1 * sigmastar/(ewu * ewv1_h) +
      ewz2 * depsi * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 - sigx6^2/sigx7)/sigx7,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (sigx14 * (S * prU * (epsilon)/ewv1 - sigx6/sigx7) +
      (((prU * (pmustar - 0.5 * (S * dmustar * ewu * (epsilon)/sqsq)) + S *
        sigx10 * ewu * (S * ewu * sigx3 * (epsilon)/(sigma_sq) - 2 * sigx4) *
        (epsilon)/sigmastar) * sigmastar + 0.5 * (prU * ewu * ewv1 * sigx3/sqsq))/ewv1_h -
        ewu^2 * ewv1_h * sigx3 * sigmastar/(ewu * ewv1_h)^2)/(sigma_sq)) *
      ewz1 * sigx1/sigx7, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx23 * (S * prU * (epsilon)/ewv1 - sigx6/sigx7) +
      (((sigx19 * sigx3 + (S * (0.5 * (prV/ewv1) + 1/(sigma_sq)) * dmustar *
        ewu * (epsilon)/sigmastar - pmustar)/(sigma_sq)) * ewu/(sigma_sq) +
        S * (2 * (sigx17 * ewu^2/sigx18) - 4/(2 * ewv1)^2) * sigx4 * (epsilon)) *
        sigmastar + 0.5 * (prV * ewu^2 * sigx3/((sigma_sq)^2 * sigmastar))) *
        ewv1/(ewu * ewv1_h) - 0.5 * (ewu^2 * ewv1_h * sigx3 * sigmastar/((ewu *
      ewv1_h)^2 * (sigma_sq)))) * ewz1 * sigx1/sigx7, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ewz2 * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx7 * ewv2_h^3) - sigx6 * sigx24 * ewv2_h/(sigx7 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((sigx1 *
    sigx5 * sigmastar/(ewu * ewv1_h) - S * depsi * (epsilon)/ewv2_h^3)/(pi *
    sigx7 * ((Wz)^2 + 1)) - pi * sigx6 * sigx25 * ((Wz)^2 + 1)/(pi * sigx7 *
    ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((prU * (ewu * (S^2 * sigx10 * dmustar * (epsilon)^2 -
      sigx8/(sigma_sq)) - 0.5 * (0.5 * (prU * dmustar * ewv1/sigmastar) + S^2 *
      sigx10 * dmustar * ewu * (epsilon)^2)) + S^2 * ((sigx8 * prU * sigx10/(sigma_sq) -
      ((0.5 * (ewu/(sigma_sq)) + 1 - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv1/sigmastar + (2 - 2 * ((0.5 * sigx9 + sigmastar)^2 * ewu *
        (sigma_sq)/sqsq^2)) * sigmastar) * sigx4/sqsq^2) * ewu + (1 - 0.5 *
      prU) * sigx10 * sigx4) * ewu * (epsilon)^2/sigmastar) * sigmastar + (0.5 *
      ((sigx8 * prU + S * ewu * pmustar * (epsilon)/(sigma_sq) - dmustar *
        sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * sigx4)) + 0.5 * (sigx11 *
      ewu/(sigma_sq))) * prU * ewv1/sigmastar)/((sigma_sq) * ewv1_h) + ewu *
      (S^2 * sigx14 * sigx10 * ewu * (epsilon)^2/sqsq - ((((sigx8 * prU - S *
        pmustar * (epsilon)) * ewu/(sigma_sq) + dmustar * sigmastar) * sigmastar +
        (0.5 * (prU * ewv1/sqsq) - 2 * (ewu * ewv1_h^2 * sigmastar/(ewu *
          ewv1_h)^2)) * sigx4 * ewu)/(ewu * ewv1_h)^2 + sigx13/((sigma_sq) *
        ewv1_h)^2) * ewv1_h) - sigx14^2 * ewz1 * sigx1/sigx7) * ewz1 * sigx1/sigx7,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((((0.5 *
    (prV * dmustar) + 0.5 * ((dmustar * ewv1/(sigma_sq) - S^2 * prV * sigx10 *
    dmustar * ewu * (epsilon)^2/sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * prV *
    dmustar)))/sigmastar + sigx8 * prU * sigx19 + (S * (pmustar - ewu * (pmustar/(sigma_sq) +
    S * sigx10 * dmustar * (epsilon))) * (epsilon) - sigx15 * ewu)/(sigma_sq))/(sigma_sq) -
    S^2 * (((0.5 * (prU * ewv1/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) - 1) * ewv1/(sigma_sq) +
      1 - 0.5 * (prU * prV))) * ewu/sigmastar + 2 * sigx17)/sigx18 - ((sqsq^2 +
      2 * ((0.5 * sigx9 + sigmastar) * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * (sqsq^2 * prU * ewv1/sigmastar)) * sigx17 * ewu/sigx18^2) * sigx4 *
      ewu * (epsilon)^2) * sigmastar + 0.5 * (((sigx8 * prU * prV + sigx4 *
    ewv1/(sigma_sq)) * ewu/(sigma_sq) + prV * sigx4)/sqsq - (0.5 * sigx9 + sigmastar) *
    prV * sigx4 * ewu/sqsq^2) + 0.5 * (sigx20 * prU * ewv1/sqsq))/ewv1_h - sigx21 *
    ewu * ewv1_h/(ewu * ewv1_h)^2) * ewv1 + ewu * (S^2 * sigx23 * sigx10 * ewu *
    (epsilon)^2/sqsq - 0.5 * ((((sigx8 * prU - S * pmustar * (epsilon)) * ewu/(sigma_sq) +
    dmustar * sigmastar) * sigmastar + (0.5 * (prU * ewv1/sqsq) - 2 * (ewu *
    ewv1_h^2 * sigmastar/(ewu * ewv1_h)^2)) * sigx4 * ewu) * ewv1_h/(ewu * ewv1_h)^2)) -
    sigx23 * sigx14 * ewz1 * sigx1/sigx7) * ewz1 * sigx1/sigx7, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx14 * ewz2 * sigx24 * ewz1 * sigx1 * ewv2_h/(sigx7 * ewv2_h)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx14 * (1/(pi * sigx7 * ((Wz)^2 + 1)) - pi * sigx25 * ((Wz)^2 + 1) * ewz1/(pi *
    sigx7 * ((Wz)^2 + 1))^2) * sigx1, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx15 * sigx19 + (S * (S * sigx17 * dmustar * ewu * (epsilon)/sqsq^2 -
      2 * (pmustar/(sigma_sq))) * (epsilon) - 0.5 * (prV * dmustar/sigmastar))/(sigma_sq)) *
      ewv1 + (0.5 * (ewv1 * (S^2 * sigx17 * dmustar * ewu^2 * (epsilon)^2/(sqsq^2 *
      sigmastar) - dmustar)/(sigma_sq) - 0.5 * (prV * dmustar)) + 0.5 * dmustar) *
      prV/sigmastar + S * pmustar * (epsilon)/(sigma_sq)) * ewu/(sigma_sq) +
      ((2 - 16 * (ewv1^2/(2 * ewv1)^2)) * (S * (epsilon))^2/(2 * ewv1)^2 -
        S^2 * (sigx17 * (1/sigx18 - ((sqsq^2 + 2 * (sigx17 * (sigma_sq)^2 *
          sigmastar)) * sigmastar + 0.5 * (sqsq^2 * prV * ewu/sigmastar)) *
          ewv1/sigx18^2) + (0.5 * (ewv1/(sigma_sq)) - 0.5 * (0.5 * prV +
          ewv1/(sigma_sq))) * prV/(sqsq^2 * ewv1)) * ewu^2 * (epsilon)^2) *
        sigx4) * sigmastar + ((0.5 * (((0.5 * (prV * dmustar/sigmastar) +
      2 * (S * pmustar * (epsilon)/(sigma_sq))) * ewu - dmustar * sigmastar)/((sigma_sq)^2 *
      sigmastar) - sigx17 * sigx4/sqsq^2) + 0.5 * (sigx20/sqsq)) * ewv1 + 0.5 *
      (sigx4/sqsq)) * prV * ewu)/(ewu * ewv1_h) + sigx23 * sigx19 - 0.5 * (sigx21 *
      ewu * ewv1_h/(ewu * ewv1_h)^2)) * ewv1 - (sigx23^2 * ewz1 * sigx1/sigx7 +
      0.5 * (((sigx15 * ewu * ewv1/(sigma_sq) + 0.5 * sigx4) * sigmastar +
        (0.5 * (prV * ewv1/sqsq) - ewu * ewv1_h^2 * sigmastar/(ewu * ewv1_h)^2) *
          sigx4 * ewu) * ewu * ewv1_h/(ewu * ewv1_h)^2))) * ewz1 * sigx1/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx23 * ewz2 * sigx24 * ewz1 * sigx1 * ewv2_h/(sigx7 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx23 * (1/(pi * sigx7 * ((Wz)^2 + 1)) - pi *
      sigx25 * ((Wz)^2 + 1) * ewz1/(pi * sigx7 * ((Wz)^2 + 1))^2) * sigx1,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ewz2 * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx7 * ewv2_h^3) - (ewz2 * sigx24 +
      0.5 * (sigx7 * ewv2_h)) * sigx24/(sigx7 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx24 * (1/(pi * sigx7 * ((Wz)^2 + 1)) + pi *
      sigx25 * ((Wz)^2 + 1) * ewz2/(pi * sigx7 * ((Wz)^2 + 1))^2)/ewv2_h),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (sigx25 * (sigx4 * sigx1 * sigmastar/(ewu *
      ewv1_h) + 2 * (pi * Wz * sigx7) - depsi/ewv2_h)/(pi * sigx7 * ((Wz)^2 +
      1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmnsfraynormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  sqsq <- ((sigma_sq) * sigmastar)
  mustar <- (S * ewu * (epsilon)/sqsq)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigx1 <- exp(0.5 * (-mustar)^2 - (S * (epsilon))^2/(2 * ewv1))
  sigx2 <- (sigmastar/ewv1 - ewu/sqsq)
  sigx3 <- (pmustar + S * dmustar * sigx2 * (epsilon))
  sigx4 <- (dmustar * sigmastar - S * ewu * pmustar * (epsilon)/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx5 <- (ewu * sigx3/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1)
  sigx6 <- (sigx1 * sigx5 * pwZ * sigmastar/(ewu * ewv1_h) + S * (1 - pwZ) * depsi *
    (epsilon)/ewv2_h^3)
  sigx7 <- ((1 - pwZ) * depsi/ewv2_h + sigx4 * sigx1 * pwZ * sigmastar/(ewu * ewv1_h))
  sigx8 <- (0.5 * (dmustar * ewv1/sigmastar) - S * pmustar * (epsilon))
  sigx9 <- (prU * ewv1/sigmastar)
  sigx10 <- (1/sqsq - (0.5 * sigx9 + sigmastar) * ewu/sqsq^2)
  sigx11 <- (sigx8 * prU + S^2 * sigx10 * sigx4 * ewu * (epsilon)^2/sigmastar)
  sigx12 <- (prU * sigx4 * ewv1/sigmastar)
  sigx13 <- (sigx11 * sigmastar + 0.5 * sigx12)
  sigx14 <- (sigx13/((sigma_sq) * ewv1_h) - sigx4 * ewu * ewv1_h * sigmastar/(ewu *
    ewv1_h)^2)
  prV <- (1 - ewv1/(sigma_sq))
  sigx15 <- (0.5 * (prV * dmustar/sigmastar) + S * pmustar * (epsilon)/(sigma_sq))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx17 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx18 <- (sqsq^2 * (sigma_sq) * sigmastar)
  sigx19 <- (2 * sigx16 - S^2 * sigx17 * ewu^2 * (epsilon)^2/sigx18)
  sigx20 <- (sigx15 * ewu/(sigma_sq) + sigx19 * sigx4)
  sigx21 <- (sigx20 * sigmastar + 0.5 * (prV * sigx4 * ewu/sqsq))
  sigx22 <- (sigx4 * ewu * ewv1_h * sigmastar/(ewu * ewv1_h)^2)
  sigx23 <- (sigx21 * ewv1/(ewu * ewv1_h) - 0.5 * sigx22)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  sigx25 <- (sigx4 * sigx1 * sigmastar/(ewu * ewv1_h) - depsi/ewv2_h)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((S^2 * ewu * (epsilon)^2/((sigma_sq) * ewv1) - 1) * sigx2 + ewu/sqsq) *
      dmustar * ewu/(sigma_sq) + prU * (S * ((2 * (S * dmustar * sigx2 * (epsilon)) +
      3 * pmustar) * ewu/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1) * (epsilon) -
      dmustar * sigmastar)/ewv1) * sigx1 * pwZ * sigmastar/(ewu * ewv1_h) +
      (1 - pwZ) * depsi * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 - sigx6^2/sigx7)/sigx7,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (sigx14 * (S * prU * (epsilon)/ewv1 - sigx6/sigx7) +
      (((prU * (pmustar - 0.5 * (S * dmustar * ewu * (epsilon)/sqsq)) + S *
        sigx10 * ewu * (S * ewu * sigx3 * (epsilon)/(sigma_sq) - 2 * sigx4) *
        (epsilon)/sigmastar) * sigmastar + 0.5 * (prU * ewu * ewv1 * sigx3/sqsq))/ewv1_h -
        ewu^2 * ewv1_h * sigx3 * sigmastar/(ewu * ewv1_h)^2)/(sigma_sq)) *
      sigx1 * pwZ/sigx7, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx23 * (S * prU * (epsilon)/ewv1 - sigx6/sigx7) +
      (((sigx19 * sigx3 + (S * (0.5 * (prV/ewv1) + 1/(sigma_sq)) * dmustar *
        ewu * (epsilon)/sigmastar - pmustar)/(sigma_sq)) * ewu/(sigma_sq) +
        S * (2 * (sigx17 * ewu^2/sigx18) - 4/(2 * ewv1)^2) * sigx4 * (epsilon)) *
        sigmastar + 0.5 * (prV * ewu^2 * sigx3/((sigma_sq)^2 * sigmastar))) *
        ewv1/(ewu * ewv1_h) - 0.5 * (ewu^2 * ewv1_h * sigx3 * sigmastar/((ewu *
      ewv1_h)^2 * (sigma_sq)))) * sigx1 * pwZ/sigx7, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (1 - pwZ) * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx7 * ewv2_h^3) - sigx24 * sigx6 * ewv2_h/(sigx7 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * dwZ * (sigx1 *
    sigx5 * sigmastar/(ewu * ewv1_h) - (sigx25 * sigx6/sigx7 + S * depsi * (epsilon)/ewv2_h^3))/sigx7,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((prU * (ewu * (S^2 * sigx10 * dmustar * (epsilon)^2 -
      sigx8/(sigma_sq)) - 0.5 * (0.5 * (prU * dmustar * ewv1/sigmastar) + S^2 *
      sigx10 * dmustar * ewu * (epsilon)^2)) + S^2 * ((sigx8 * prU * sigx10/(sigma_sq) -
      ((0.5 * (ewu/(sigma_sq)) + 1 - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv1/sigmastar + (2 - 2 * ((0.5 * sigx9 + sigmastar)^2 * ewu *
        (sigma_sq)/sqsq^2)) * sigmastar) * sigx4/sqsq^2) * ewu + (1 - 0.5 *
      prU) * sigx10 * sigx4) * ewu * (epsilon)^2/sigmastar) * sigmastar + (0.5 *
      ((sigx8 * prU + S * ewu * pmustar * (epsilon)/(sigma_sq) - dmustar *
        sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * sigx4)) + 0.5 * (sigx11 *
      ewu/(sigma_sq))) * prU * ewv1/sigmastar)/((sigma_sq) * ewv1_h) + ewu *
      (S^2 * sigx14 * sigx10 * ewu * (epsilon)^2/sqsq - ((((sigx8 * prU - S *
        pmustar * (epsilon)) * ewu/(sigma_sq) + dmustar * sigmastar) * sigmastar +
        (0.5 * (prU * ewv1/sqsq) - 2 * (ewu * ewv1_h^2 * sigmastar/(ewu *
          ewv1_h)^2)) * sigx4 * ewu)/(ewu * ewv1_h)^2 + sigx13/((sigma_sq) *
        ewv1_h)^2) * ewv1_h) - sigx14^2 * sigx1 * pwZ/sigx7) * sigx1 * pwZ/sigx7,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((((0.5 *
    (prV * dmustar) + 0.5 * ((dmustar * ewv1/(sigma_sq) - S^2 * prV * sigx10 *
    dmustar * ewu * (epsilon)^2/sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * prV *
    dmustar)))/sigmastar + sigx8 * prU * sigx19 + (S * (pmustar - ewu * (pmustar/(sigma_sq) +
    S * sigx10 * dmustar * (epsilon))) * (epsilon) - sigx15 * ewu)/(sigma_sq))/(sigma_sq) -
    S^2 * (((0.5 * (prU * ewv1/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) - 1) * ewv1/(sigma_sq) +
      1 - 0.5 * (prU * prV))) * ewu/sigmastar + 2 * sigx17)/sigx18 - ((sqsq^2 +
      2 * ((0.5 * sigx9 + sigmastar) * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * (sqsq^2 * prU * ewv1/sigmastar)) * sigx17 * ewu/sigx18^2) * sigx4 *
      ewu * (epsilon)^2) * sigmastar + 0.5 * (((sigx8 * prU * prV + sigx4 *
    ewv1/(sigma_sq)) * ewu/(sigma_sq) + prV * sigx4)/sqsq - (0.5 * sigx9 + sigmastar) *
    prV * sigx4 * ewu/sqsq^2) + 0.5 * (sigx20 * prU * ewv1/sqsq))/ewv1_h - sigx21 *
    ewu * ewv1_h/(ewu * ewv1_h)^2) * ewv1 + ewu * (S^2 * sigx23 * sigx10 * ewu *
    (epsilon)^2/sqsq - 0.5 * ((((sigx8 * prU - S * pmustar * (epsilon)) * ewu/(sigma_sq) +
    dmustar * sigmastar) * sigmastar + (0.5 * (prU * ewv1/sqsq) - 2 * (ewu *
    ewv1_h^2 * sigmastar/(ewu * ewv1_h)^2)) * sigx4 * ewu) * ewv1_h/(ewu * ewv1_h)^2)) -
    sigx23 * sigx14 * sigx1 * pwZ/sigx7) * sigx1 * pwZ/sigx7, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx14 * sigx24 * (1 - pwZ) * sigx1 * ewv2_h * pwZ/(sigx7 * ewv2_h)^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx14 * (1 - sigx25 * pwZ/sigx7) * dwZ * sigx1/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx15 * sigx19 + (S * (S * sigx17 * dmustar * ewu * (epsilon)/sqsq^2 -
      2 * (pmustar/(sigma_sq))) * (epsilon) - 0.5 * (prV * dmustar/sigmastar))/(sigma_sq)) *
      ewv1 + (0.5 * (ewv1 * (S^2 * sigx17 * dmustar * ewu^2 * (epsilon)^2/(sqsq^2 *
      sigmastar) - dmustar)/(sigma_sq) - 0.5 * (prV * dmustar)) + 0.5 * dmustar) *
      prV/sigmastar + S * pmustar * (epsilon)/(sigma_sq)) * ewu/(sigma_sq) +
      ((2 - 16 * (ewv1^2/(2 * ewv1)^2)) * (S * (epsilon))^2/(2 * ewv1)^2 -
        S^2 * (sigx17 * (1/sigx18 - ((sqsq^2 + 2 * (sigx17 * (sigma_sq)^2 *
          sigmastar)) * sigmastar + 0.5 * (sqsq^2 * prV * ewu/sigmastar)) *
          ewv1/sigx18^2) + (0.5 * (ewv1/(sigma_sq)) - 0.5 * (0.5 * prV +
          ewv1/(sigma_sq))) * prV/(sqsq^2 * ewv1)) * ewu^2 * (epsilon)^2) *
        sigx4) * sigmastar + ((0.5 * (((0.5 * (prV * dmustar/sigmastar) +
      2 * (S * pmustar * (epsilon)/(sigma_sq))) * ewu - dmustar * sigmastar)/((sigma_sq)^2 *
      sigmastar) - sigx17 * sigx4/sqsq^2) + 0.5 * (sigx20/sqsq)) * ewv1 + 0.5 *
      (sigx4/sqsq)) * prV * ewu)/(ewu * ewv1_h) + sigx23 * sigx19 - 0.5 * (sigx21 *
      ewu * ewv1_h/(ewu * ewv1_h)^2)) * ewv1 - (sigx23^2 * sigx1 * pwZ/sigx7 +
      0.5 * (((sigx15 * ewu * ewv1/(sigma_sq) + 0.5 * sigx4) * sigmastar +
        (0.5 * (prV * ewv1/sqsq) - ewu * ewv1_h^2 * sigmastar/(ewu * ewv1_h)^2) *
          sigx4 * ewu) * ewu * ewv1_h/(ewu * ewv1_h)^2))) * sigx1 * pwZ/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx23 * sigx24 * (1 - pwZ) * sigx1 * ewv2_h * pwZ/(sigx7 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx23 * (1 - sigx25 * pwZ/sigx7) * dwZ * sigx1/sigx7,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - pwZ) * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx7 * ewv2_h^3) - (sigx24 * (1 -
      pwZ) + 0.5 * (sigx7 * ewv2_h)) * sigx24/(sigx7 * ewv2_h)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx25 * (1 - pwZ)/sigx7 + 1) * sigx24 * dwZ/(sigx7 *
      ewv2_h)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx25 * dwZ/sigx7 + Wz) * sigx25 * dwZ/sigx7),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmnsfraynormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  sqsq <- ((sigma_sq) * sigmastar)
  mustar <- (S * ewu * (epsilon)/sqsq)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigx1 <- exp(0.5 * (-mustar)^2 - (S * (epsilon))^2/(2 * ewv1))
  sigx2 <- (sigmastar/ewv1 - ewu/sqsq)
  sigx3 <- (pmustar + S * dmustar * sigx2 * (epsilon))
  sigx4 <- (dmustar * sigmastar - S * ewu * pmustar * (epsilon)/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx5 <- (ewu * sigx3/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1)
  sigx6 <- ((1 - prZ) * sigx1 * sigx5 * sigmastar/(ewu * ewv1_h) + S * depsi *
    prZ * (epsilon)/ewv2_h^3)
  sigx7 <- ((1 - prZ) * sigx4 * sigx1 * sigmastar/(ewu * ewv1_h) + depsi * prZ/ewv2_h)
  sigx8 <- (0.5 * (dmustar * ewv1/sigmastar) - S * pmustar * (epsilon))
  sigx9 <- (prU * ewv1/sigmastar)
  sigx10 <- (1/sqsq - (0.5 * sigx9 + sigmastar) * ewu/sqsq^2)
  sigx11 <- (sigx8 * prU + S^2 * sigx10 * sigx4 * ewu * (epsilon)^2/sigmastar)
  sigx12 <- (prU * sigx4 * ewv1/sigmastar)
  sigx13 <- (sigx11 * sigmastar + 0.5 * sigx12)
  sigx14 <- (sigx13/((sigma_sq) * ewv1_h) - sigx4 * ewu * ewv1_h * sigmastar/(ewu *
    ewv1_h)^2)
  prV <- (1 - ewv1/(sigma_sq))
  sigx15 <- (0.5 * (prV * dmustar/sigmastar) + S * pmustar * (epsilon)/(sigma_sq))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx17 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx18 <- (sqsq^2 * (sigma_sq) * sigmastar)
  sigx19 <- (2 * sigx16 - S^2 * sigx17 * ewu^2 * (epsilon)^2/sigx18)
  sigx20 <- (sigx15 * ewu/(sigma_sq) + sigx19 * sigx4)
  sigx21 <- (sigx20 * sigmastar + 0.5 * (prV * sigx4 * ewu/sqsq))
  sigx22 <- (sigx4 * ewu * ewv1_h * sigmastar/(ewu * ewv1_h)^2)
  sigx23 <- (sigx21 * ewv1/(ewu * ewv1_h) - 0.5 * sigx22)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  sigx25 <- (sigx4 * sigx1 * sigmastar/(ewu * ewv1_h) - depsi/ewv2_h)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((S^2 * ewu * (epsilon)^2/((sigma_sq) * ewv1) - 1) * sigx2 + ewu/sqsq) *
      dmustar * ewu/(sigma_sq) + prU * (S * ((2 * (S * dmustar * sigx2 * (epsilon)) +
      3 * pmustar) * ewu/(sigma_sq) + S * prU * sigx4 * (epsilon)/ewv1) * (epsilon) -
      dmustar * sigmastar)/ewv1) * (1 - prZ) * sigx1 * sigmastar/(ewu * ewv1_h) +
      depsi * prZ * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 - sigx6^2/sigx7)/sigx7,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (sigx14 * (S * prU * (epsilon)/ewv1 - sigx6/sigx7) +
      (((prU * (pmustar - 0.5 * (S * dmustar * ewu * (epsilon)/sqsq)) + S *
        sigx10 * ewu * (S * ewu * sigx3 * (epsilon)/(sigma_sq) - 2 * sigx4) *
        (epsilon)/sigmastar) * sigmastar + 0.5 * (prU * ewu * ewv1 * sigx3/sqsq))/ewv1_h -
        ewu^2 * ewv1_h * sigx3 * sigmastar/(ewu * ewv1_h)^2)/(sigma_sq)) *
      (1 - prZ) * sigx1/sigx7, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx23 * (S * prU * (epsilon)/ewv1 - sigx6/sigx7) +
      (((sigx19 * sigx3 + (S * (0.5 * (prV/ewv1) + 1/(sigma_sq)) * dmustar *
        ewu * (epsilon)/sigmastar - pmustar)/(sigma_sq)) * ewu/(sigma_sq) +
        S * (2 * (sigx17 * ewu^2/sigx18) - 4/(2 * ewv1)^2) * sigx4 * (epsilon)) *
        sigmastar + 0.5 * (prV * ewu^2 * sigx3/((sigma_sq)^2 * sigmastar))) *
        ewv1/(ewu * ewv1_h) - 0.5 * (ewu^2 * ewv1_h * sigx3 * sigmastar/((ewu *
      ewv1_h)^2 * (sigma_sq)))) * (1 - prZ) * sigx1/sigx7, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * prZ * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx7 * ewv2_h^3) - sigx6 * sigx24 * ewv2_h/(sigx7 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * prZ * (sigx1 *
    sigx5 * sigmastar/(ewu * ewv1_h) - (sigx6 * sigx25/sigx7 + S * depsi * (epsilon)/ewv2_h^3)) *
    ewz/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((prU * (ewu * (S^2 * sigx10 * dmustar * (epsilon)^2 -
      sigx8/(sigma_sq)) - 0.5 * (0.5 * (prU * dmustar * ewv1/sigmastar) + S^2 *
      sigx10 * dmustar * ewu * (epsilon)^2)) + S^2 * ((sigx8 * prU * sigx10/(sigma_sq) -
      ((0.5 * (ewu/(sigma_sq)) + 1 - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv1/sigmastar + (2 - 2 * ((0.5 * sigx9 + sigmastar)^2 * ewu *
        (sigma_sq)/sqsq^2)) * sigmastar) * sigx4/sqsq^2) * ewu + (1 - 0.5 *
      prU) * sigx10 * sigx4) * ewu * (epsilon)^2/sigmastar) * sigmastar + (0.5 *
      ((sigx8 * prU + S * ewu * pmustar * (epsilon)/(sigma_sq) - dmustar *
        sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * sigx4)) + 0.5 * (sigx11 *
      ewu/(sigma_sq))) * prU * ewv1/sigmastar)/((sigma_sq) * ewv1_h) + ewu *
      (S^2 * sigx14 * sigx10 * ewu * (epsilon)^2/sqsq - ((((sigx8 * prU - S *
        pmustar * (epsilon)) * ewu/(sigma_sq) + dmustar * sigmastar) * sigmastar +
        (0.5 * (prU * ewv1/sqsq) - 2 * (ewu * ewv1_h^2 * sigmastar/(ewu *
          ewv1_h)^2)) * sigx4 * ewu)/(ewu * ewv1_h)^2 + sigx13/((sigma_sq) *
        ewv1_h)^2) * ewv1_h) - sigx14^2 * (1 - prZ) * sigx1/sigx7) * (1 -
      prZ) * sigx1/sigx7, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((((0.5 *
    (prV * dmustar) + 0.5 * ((dmustar * ewv1/(sigma_sq) - S^2 * prV * sigx10 *
    dmustar * ewu * (epsilon)^2/sigmastar) * ewu/(sigma_sq) - 0.5 * (prU * prV *
    dmustar)))/sigmastar + sigx8 * prU * sigx19 + (S * (pmustar - ewu * (pmustar/(sigma_sq) +
    S * sigx10 * dmustar * (epsilon))) * (epsilon) - sigx15 * ewu)/(sigma_sq))/(sigma_sq) -
    S^2 * (((0.5 * (prU * ewv1/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) - 1) * ewv1/(sigma_sq) +
      1 - 0.5 * (prU * prV))) * ewu/sigmastar + 2 * sigx17)/sigx18 - ((sqsq^2 +
      2 * ((0.5 * sigx9 + sigmastar) * (sigma_sq)^2 * sigmastar)) * sigmastar +
      0.5 * (sqsq^2 * prU * ewv1/sigmastar)) * sigx17 * ewu/sigx18^2) * sigx4 *
      ewu * (epsilon)^2) * sigmastar + 0.5 * (((sigx8 * prU * prV + sigx4 *
    ewv1/(sigma_sq)) * ewu/(sigma_sq) + prV * sigx4)/sqsq - (0.5 * sigx9 + sigmastar) *
    prV * sigx4 * ewu/sqsq^2) + 0.5 * (sigx20 * prU * ewv1/sqsq))/ewv1_h - sigx21 *
    ewu * ewv1_h/(ewu * ewv1_h)^2) * ewv1 + ewu * (S^2 * sigx23 * sigx10 * ewu *
    (epsilon)^2/sqsq - 0.5 * ((((sigx8 * prU - S * pmustar * (epsilon)) * ewu/(sigma_sq) +
    dmustar * sigmastar) * sigmastar + (0.5 * (prU * ewv1/sqsq) - 2 * (ewu *
    ewv1_h^2 * sigmastar/(ewu * ewv1_h)^2)) * sigx4 * ewu) * ewv1_h/(ewu * ewv1_h)^2)) -
    sigx23 * sigx14 * (1 - prZ) * sigx1/sigx7) * (1 - prZ) * sigx1/sigx7, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx14 * sigx24 * (1 - prZ) * prZ * sigx1 * ewv2_h/(sigx7 * ewv2_h)^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx14 * (1 - sigx25 * (1 - prZ)/sigx7) * prZ * sigx1 * ewz/sigx7, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx15 * sigx19 + (S * (S * sigx17 * dmustar * ewu * (epsilon)/sqsq^2 -
      2 * (pmustar/(sigma_sq))) * (epsilon) - 0.5 * (prV * dmustar/sigmastar))/(sigma_sq)) *
      ewv1 + (0.5 * (ewv1 * (S^2 * sigx17 * dmustar * ewu^2 * (epsilon)^2/(sqsq^2 *
      sigmastar) - dmustar)/(sigma_sq) - 0.5 * (prV * dmustar)) + 0.5 * dmustar) *
      prV/sigmastar + S * pmustar * (epsilon)/(sigma_sq)) * ewu/(sigma_sq) +
      ((2 - 16 * (ewv1^2/(2 * ewv1)^2)) * (S * (epsilon))^2/(2 * ewv1)^2 -
        S^2 * (sigx17 * (1/sigx18 - ((sqsq^2 + 2 * (sigx17 * (sigma_sq)^2 *
          sigmastar)) * sigmastar + 0.5 * (sqsq^2 * prV * ewu/sigmastar)) *
          ewv1/sigx18^2) + (0.5 * (ewv1/(sigma_sq)) - 0.5 * (0.5 * prV +
          ewv1/(sigma_sq))) * prV/(sqsq^2 * ewv1)) * ewu^2 * (epsilon)^2) *
        sigx4) * sigmastar + ((0.5 * (((0.5 * (prV * dmustar/sigmastar) +
      2 * (S * pmustar * (epsilon)/(sigma_sq))) * ewu - dmustar * sigmastar)/((sigma_sq)^2 *
      sigmastar) - sigx17 * sigx4/sqsq^2) + 0.5 * (sigx20/sqsq)) * ewv1 + 0.5 *
      (sigx4/sqsq)) * prV * ewu)/(ewu * ewv1_h) + sigx23 * sigx19 - 0.5 * (sigx21 *
      ewu * ewv1_h/(ewu * ewv1_h)^2)) * ewv1 - (sigx23^2 * (1 - prZ) * sigx1/sigx7 +
      0.5 * (((sigx15 * ewu * ewv1/(sigma_sq) + 0.5 * sigx4) * sigmastar +
        (0.5 * (prV * ewv1/sqsq) - ewu * ewv1_h^2 * sigmastar/(ewu * ewv1_h)^2) *
          sigx4 * ewu) * ewu * ewv1_h/(ewu * ewv1_h)^2))) * (1 - prZ) * sigx1/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx23 * sigx24 * (1 - prZ) * prZ * sigx1 * ewv2_h/(sigx7 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx23 * (1 - sigx25 * (1 - prZ)/sigx7) * prZ *
      sigx1 * ewz/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prZ * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx7 * ewv2_h^3) - (sigx24 * prZ +
      0.5 * (sigx7 * ewv2_h)) * sigx24/(sigx7 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx25 * prZ/sigx7 + 1) * sigx24 * prZ * ewz/(sigx7 *
      ewv2_h)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * sigx25 * (1 - (sigx25 * prZ/sigx7 + 1) * ewz) *
      prZ * ewz/sigx7, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf rayleigh-normal distribution
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
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
# Same sigma_v

## logit specification class membership
zisfraynormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfraynormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfraynormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfraynormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfraynormlike_logit,
    grad = cgradzisfraynormlike_logit, hess = chesszisfraynormlike_logit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfraynormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfraynormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfraynormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfraynormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo, eval.max = itermax,
        rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfraynormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfraynormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfraynormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfraynormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfraynormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initRay = initRay))
}

## cauchit specification class membership
zisfraynormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfraynormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfraynormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfraynormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfraynormlike_cauchit,
    grad = cgradzisfraynormlike_cauchit, hess = chesszisfraynormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfraynormlike_cauchit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfraynormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfraynormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfraynormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfraynormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfraynormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfraynormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfraynormlike_cauchit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfraynormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initRay = initRay))
}

## probit specification class membership
zisfraynormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfraynormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfraynormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfraynormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfraynormlike_probit,
    grad = cgradzisfraynormlike_probit, hess = chesszisfraynormlike_probit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfraynormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfraynormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfraynormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfraynormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfraynormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfraynormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfraynormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfraynormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfraynormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initRay = initRay))
}

## cloglog specification class membership
zisfraynormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfraynormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfraynormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfraynormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfraynormlike_cloglog,
    grad = cgradzisfraynormlike_cloglog, hess = chesszisfraynormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfraynormlike_cloglog(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfraynormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfraynormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfraynormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfraynormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfraynormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfraynormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfraynormlike_cloglog(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfraynormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initRay = initRay))
}

# Different sigma_v

## logit specification class membership
mnsfraynormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfraynormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfraynormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfraynormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfraynormlike_logit,
    grad = cgradmnsfraynormlike_logit, hess = chessmnsfraynormlike_logit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfraynormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfraynormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfraynormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfraynormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfraynormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo, eval.max = itermax,
        rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfraynormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfraynormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfraynormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfraynormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfraynormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initRay = initRay))
}

## cauchit specification class membership
mnsfraynormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfraynormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfraynormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfraynormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfraynormlike_cauchit,
    grad = cgradmnsfraynormlike_cauchit, hess = chessmnsfraynormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfraynormlike_cauchit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfraynormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfraynormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfraynormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfraynormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfraynormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfraynormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfraynormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfraynormlike_cauchit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfraynormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initRay = initRay))
}

## probit specification class membership
mnsfraynormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfraynormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfraynormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfraynormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfraynormlike_probit,
    grad = cgradmnsfraynormlike_probit, hess = chessmnsfraynormlike_probit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfraynormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfraynormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfraynormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfraynormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfraynormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfraynormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfraynormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfraynormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfraynormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfraynormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initRay = initRay))
}

## cloglog specification class membership
mnsfraynormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfraynormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfraynormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfraynormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfraynormlike_cloglog,
    grad = cgradmnsfraynormlike_cloglog, hess = chessmnsfraynormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfraynormlike_cloglog(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfraynormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfraynormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfraynormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfraynormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfraynormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfraynormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfraynormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfraynormlike_cloglog(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfraynormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initRay = initRay))
}

# Conditional efficiencies estimation ----------
#' efficiencies for zisf rayleigh-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v

## logit specification class membership
czisfraynormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar - sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar +
      sigmastar) + (mustar + sigmastar^2) * pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
czisfraynormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar - sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar +
      sigmastar) + (mustar + sigmastar^2) * pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
czisfraynormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar - sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar +
      sigmastar) + (mustar + sigmastar^2) * pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
czisfraynormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar - sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar +
      sigmastar) + (mustar + sigmastar^2) * pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
cmnsfraynormeff_logit <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar - sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar +
      sigmastar) + (mustar + sigmastar^2) * pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
cmnsfraynormeff_cauchit <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar - sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar +
      sigmastar) + (mustar + sigmastar^2) * pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
cmnsfraynormeff_probit <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar - sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar +
      sigmastar) + (mustar + sigmastar^2) * pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
cmnsfraynormeff_cloglog <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar - sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar +
      sigmastar) + (mustar + sigmastar^2) * pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
#' marginal impact on efficiencies for zisf rayleigh-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v

## logit specification class membership
czisfmargraynorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargraynorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (4 - pi)/2, ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
czisfmargraynorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargraynorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (4 - pi)/2, ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
czisfmargraynorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargraynorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (4 - pi)/2, ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
czisfmargraynorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargraynorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv)))/(exp(Wv/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (4 - pi)/2, ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# Different sigma_v

## logit specification class membership
cmnsfmargraynorm_Eu_logit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargraynorm_Vu_logit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (4 - pi)/2, ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmnsfmargraynorm_Eu_cauchit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargraynorm_Vu_cauchit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (4 - pi)/2, ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmnsfmargraynorm_Eu_probit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargraynorm_Vu_probit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (4 - pi)/2, ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmnsfmargraynorm_Eu_cloglog <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) *
    1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargraynorm_Vu_cloglog <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- exp(1/2 * (mustar/sigmastar)^2 - (object$S * epsilon)^2/(2 * exp(Wv1)))/(exp(Wv1/2) *
    exp(Wu)) * sigmastar * (sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (4 - pi)/2, ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
