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
# Convolution: generalized genexponential - normal                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf generalized_genexponential-normal distribution
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
czisfgenexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cauchit specification class membership
czisfgenexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## probit specification class membership
czisfgenexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cloglog specification class membership
czisfgenexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}
# Different sigma_v

## logit specification class membership
cmnsfgenexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  A <- S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cauchit specification class membership
cmnsfgenexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  A <- S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## probit specification class membership
cmnsfgenexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  A <- S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cloglog specification class membership
cmnsfgenexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  A <- S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf generalized_genexponential-normal distribution
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
cstzisfgenexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- cstgenexponorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGenExpo <- NULL
  } else {
    cat("Initialization: SFA + generalized genexponential - normal distributions...\n")
    initGenExpo <- maxLik::maxLik(logLik = cgenexponormlike, start = cstgenexponorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradgenexponormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initGenExpo$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("ZISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
}

# Different sigma_v
cstmnsfgenexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- cstgenexponorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGenExpo <- NULL
  } else {
    cat("Initialization: SFA + generalized genexponential - normal distributions...\n")
    initGenExpo <- maxLik::maxLik(logLik = cgenexponormlike, start = cstgenexponorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradgenexponormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initGenExpo$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar +
    2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("ZISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
}

# Gradient of the likelihood function ----------
#' gradient for zisf generalized_genexponential-normal distribution
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
cgradzisfgenexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) * eB)
  sigx2 <- (2 * (sigx1 * ewz/(wzdeno * ewusr)) + S * prC * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (prC * dwsr/ewvsr + 2 * ((eA * pa - eB * pb) * ewz/(wzdeno * ewusr)))
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx6 <- (sigx5 * eA - (db * ewvsr/ewusr - eC * pb) * eB)
  sigx7 <- (wzdeno * ewusr * (eA * pa - eB * pb)/(wzdeno * ewusr)^2)
  sigx8 <- (sigx6/(wzdeno * ewusr) - 0.5 * sigx7)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)))
  sigx12 <- ((0.5 * depsisr - 0.5 * dwsr) * prC/ewvsr + 2 * ((eA * sigx10 - sigx11 *
    eB) * ewz/(wzdeno * ewusr)))
  sigx13 <- ((1/(wzdeno * ewusr) - ewusr * ewz/(wzdeno * ewusr)^2) * (eA * pa -
    eB * pb))
  sigx14 <- (2 * sigx13 - prC * dwsr/(wzdeno * ewvsr))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx8 * ewz/sigx3), FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx12/sigx3, FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx14 * ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradzisfgenexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) * eB)
  sigx2 <- (2 * (sigx1 * ewz1/ewusr) + S * ewz2 * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (ewz2 * dwsr/ewvsr + 2 * (ewz1 * (eA * pa - eB * pb)/ewusr))
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx7 <- (sigx5 * eA - ((db * ewvsr/ewusr - eC * pb) * eB + 0.5 * (eA * pa -
    eB * pb)))
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)))
  sigx8 <- (ewz1 * (eA * sigx10 - sigx11 * eB)/ewusr)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx12 <- (ewz2 * (0.5 * depsisr - 0.5 * dwsr)/ewvsr + 2 * sigx8)
  sigx13 <- (2 * ((eA * pa - eB * pb)/ewusr) - dwsr/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx7 * ewz1/(sigx3 * ewusr)), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx12/sigx3, FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx13/(pi * sigx3 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradzisfgenexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) * eB)
  sigx2 <- (2 * (sigx1 * pwZ/ewusr) + S * (1 - pwZ) * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - pwZ) * dwsr/ewvsr + 2 * ((eA * pa - eB * pb) * pwZ/ewusr))
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx7 <- (sigx5 * eA - ((db * ewvsr/ewusr - eC * pb) * eB + 0.5 * (eA * pa -
    eB * pb)))
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)))
  sigx8 <- ((eA * sigx10 - sigx11 * eB) * pwZ/ewusr)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx12 <- ((0.5 * depsisr - 0.5 * dwsr) * (1 - pwZ)/ewvsr + 2 * sigx8)
  sigx13 <- (2 * ((eA * pa - eB * pb)/ewusr) - dwsr/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx7 * pwZ/(sigx3 * ewusr)), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx12/sigx3, FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx13 * dwZ/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradzisfgenexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) * eB)
  sigx2 <- (2 * (sigx1 * (1 - prZ)/ewusr) + S * dwsr * prZ * (epsilon)/ewvsr^3)
  sigx3 <- (2 * ((1 - prZ) * (eA * pa - eB * pb)/ewusr) + dwsr * prZ/ewvsr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx7 <- (sigx5 * eA - ((db * ewvsr/ewusr - eC * pb) * eB + 0.5 * (eA * pa -
    eB * pb)))
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)))
  sigx8 <- ((1 - prZ) * (eA * sigx10 - sigx11 * eB)/ewusr)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx12 <- ((0.5 * depsisr - 0.5 * dwsr) * prZ/ewvsr + 2 * sigx8)
  sigx13 <- (2 * ((eA * pa - eB * pb)/ewusr) - dwsr/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx7 * (1 - prZ)/(sigx3 * ewusr)),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx12/sigx3, FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx13 * prZ * ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v

## logit specification class membership
cgradmnsfgenexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewu_h <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  mustar1 <- (ewv1_h/ewu_h + S * epsilon/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(S * epsilon/ewv2_h, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  eusq1 <- exp(ewv1/(2 * ewu) + S * epsilon/ewu_h)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * epsilon/ewu_h))
  evsq1 <- (2 * (ewv1_h/ewu_h) + S * epsilon/ewv1_h)
  depsi <- dnorm(-evsq1, 0, 1)
  pepsi <- pnorm(-evsq1)
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewu_h)
  dpepsi <- (depsi/ewv1_h - 2 * (pepsi/ewu_h))
  epepe <- (eusq1 * pmustar1 - eusq2 * pepsi)
  exusq <- (sigx1 * eusq1 - dpepsi * eusq2)
  ewvu <- (ewu * ewv1/(2 * ewu)^2)
  sigx2 <- (exusq * ewz/(wzdeno * ewu_h))
  sigx3 <- (2 * sigx2 + S * prC * dmustar2 * epsilon/ewv2_h^3)
  sigx4 <- (epepe * ewz/(wzdeno * ewu_h))
  sigx5 <- (prC * dmustar2/ewv2_h + 2 * sigx4)
  sigx6 <- (0.5 * (S * epsilon/ewu_h) + 2 * ewvu)
  sigx7 <- (0.5 * (dmustar1 * ewv1_h/ewu_h) - sigx6 * pmustar1)
  sigx8 <- (depsi * ewv1_h/ewu_h - (2 * (ewv1/ewu) + S * epsilon/ewu_h) * pepsi)
  sigx9 <- (sigx7 * eusq1 - sigx8 * eusq2)
  sigx10 <- (wzdeno * ewu_h * epepe/(wzdeno * ewu_h)^2)
  sigx11 <- (sigx9/(wzdeno * ewu_h) - 0.5 * sigx10)
  sigx12 <- (ewv1 * pmustar1/(2 * ewu) - (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h)) *
    dmustar1)
  sigx13 <- (2 * (ewv1 * pepsi/ewu) - depsi * (ewv1_h/ewu_h - 0.5 * (S * epsilon/ewv1_h)))
  sigx14 <- (eusq1 * sigx12 - sigx13 * eusq2)
  sigx15 <- (0.5 * (S^2 * dmustar2 * epsilon^2/ewv2_h^2) - 0.5 * dmustar2)
  sigx16 <- (1/(wzdeno * ewu_h) - ewu_h * ewz/(wzdeno * ewu_h)^2)
  sigx17 <- (sigx16 * epepe)
  sigx18 <- (2 * sigx17 - prC * dmustar2/(wzdeno * ewv2_h))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx5, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx11 * ewz/sigx5), FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = 2 * (sigx14 * ewz/(sigx5 * wzdeno * ewu_h)), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx15 * prC/(sigx5 * ewv2_h), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx18 * ewz/sigx5, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmnsfgenexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  mustar1 <- (ewv1_h/ewusr + S * (epsilon)/ewv1_h)
  mustar2 <- (2 * (ewv1_h/ewusr) + S * (epsilon)/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(-mustar2, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  pmustar2 <- pnorm(-mustar2)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  eusq1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * (epsilon)/ewusr))
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewusr)
  sigx2 <- (dmustar2/ewv1_h - 2 * (pmustar2/ewusr))
  sigx3 <- ((sigx1 * eusq1 - sigx2 * eusq2) * ewz1/ewusr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar1 * ewv1_h/ewusr) - sigx4 * pmustar1)
  sigx6 <- (dmustar2 * ewv1_h/ewusr - (2 * (ewv1/ewu) + S * (epsilon)/ewusr) *
    pmustar2)
  sigx7 <- (sigx6 * eusq2 + 0.5 * (eusq1 * pmustar1 - eusq2 * pmustar2))
  sigx8 <- (ewz2 * depsi/ewv2_h + 2 * (ewz1 * (eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr))
  sigx9 <- (2 * sigx3 + S * ewz2 * depsi * (epsilon)/ewv2_h^3)
  sigx10 <- (0.5 * (ewv1_h/ewusr) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx11 <- (ewv1 * pmustar1/(2 * ewu) - sigx10 * dmustar1)
  sigx12 <- (ewv1_h/ewusr - 0.5 * (S * (epsilon)/ewv1_h))
  sigx13 <- (eusq1 * sigx11 - (2 * (ewv1 * pmustar2/ewu) - dmustar2 * sigx12) *
    eusq2)
  sigx14 <- (2 * ((eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr) - depsi/ewv2_h)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx9/sigx8, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * ((sigx5 * eusq1 - sigx7) * ewz1/(sigx8 *
      ewusr)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * (ewz1 * sigx13/(sigx8 *
      ewusr)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ewz2 * sigx24/(sigx8 *
      ewv2_h), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx14/(pi * sigx8 *
      ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmnsfgenexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  mustar1 <- (ewv1_h/ewusr + S * (epsilon)/ewv1_h)
  mustar2 <- (2 * (ewv1_h/ewusr) + S * (epsilon)/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(-mustar2, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  pmustar2 <- pnorm(-mustar2)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  eusq1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * (epsilon)/ewusr))
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewusr)
  sigx2 <- (dmustar2/ewv1_h - 2 * (pmustar2/ewusr))
  sigx3 <- ((sigx1 * eusq1 - sigx2 * eusq2) * pwZ/ewusr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar1 * ewv1_h/ewusr) - sigx4 * pmustar1)
  sigx6 <- (dmustar2 * ewv1_h/ewusr - (2 * (ewv1/ewu) + S * (epsilon)/ewusr) *
    pmustar2)
  sigx7 <- (sigx6 * eusq2 + 0.5 * (eusq1 * pmustar1 - eusq2 * pmustar2))
  sigx9 <- (2 * sigx3 + S * (1 - pwZ) * depsi * (epsilon)/ewv2_h^3)
  sigx10 <- (0.5 * (ewv1_h/ewusr) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx11 <- (ewv1 * pmustar1/(2 * ewu) - sigx10 * dmustar1)
  sigx12 <- (ewv1_h/ewusr - 0.5 * (S * (epsilon)/ewv1_h))
  sigx13 <- (eusq1 * sigx11 - (2 * (ewv1 * pmustar2/ewu) - dmustar2 * sigx12) *
    eusq2)
  sigx14 <- (2 * ((eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr) - depsi/ewv2_h)
  sigx15 <- ((1 - pwZ) * depsi/ewv2_h + 2 * ((eusq1 * pmustar1 - eusq2 * pmustar2) *
    pwZ/ewusr))
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx9/sigx15, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * ((sigx5 * eusq1 - sigx7) * pwZ/(sigx15 *
      ewusr)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * (sigx13 * pwZ/(sigx15 *
      ewusr)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx24 * (1 - pwZ)/(sigx15 *
      ewv2_h), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx14 * dwZ/sigx15,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmnsfgenexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  mustar1 <- (ewv1_h/ewusr + S * (epsilon)/ewv1_h)
  mustar2 <- (2 * (ewv1_h/ewusr) + S * (epsilon)/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(-mustar2, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  pmustar2 <- pnorm(-mustar2)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  eusq1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * (epsilon)/ewusr))
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewusr)
  sigx2 <- (dmustar2/ewv1_h - 2 * (pmustar2/ewusr))
  sigx3 <- ((sigx1 * eusq1 - sigx2 * eusq2) * (1 - prZ)/ewusr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar1 * ewv1_h/ewusr) - sigx4 * pmustar1)
  sigx6 <- (dmustar2 * ewv1_h/ewusr - (2 * (ewv1/ewu) + S * (epsilon)/ewusr) *
    pmustar2)
  sigx7 <- (sigx6 * eusq2 + 0.5 * (eusq1 * pmustar1 - eusq2 * pmustar2))
  sigx8 <- (2 * ((1 - prZ) * (eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr) + depsi *
    prZ/ewv2_h)
  sigx9 <- (2 * sigx3 + S * depsi * prZ * (epsilon)/ewv2_h^3)
  sigx10 <- (0.5 * (ewv1_h/ewusr) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx11 <- (ewv1 * pmustar1/(2 * ewu) - sigx10 * dmustar1)
  sigx12 <- (ewv1_h/ewusr - 0.5 * (S * (epsilon)/ewv1_h))
  sigx13 <- (eusq1 * sigx11 - (2 * (ewv1 * pmustar2/ewu) - dmustar2 * sigx12) *
    eusq2)
  sigx14 <- (2 * ((eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr) - depsi/ewv2_h)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx9/sigx8, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * ((sigx5 * eusq1 - sigx7) * (1 - prZ)/(sigx8 *
      ewusr)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * ((1 - prZ) *
      sigx13/(sigx8 * ewusr)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx24 *
      prZ/(sigx8 * ewv2_h), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx14 *
      prZ * ewz/sigx8, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf generalized_genexponential-normal distribution
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
chesszisfgenexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) * eB)
  sigx2 <- (2 * (sigx1 * ewz/(wzdeno * ewusr)) + S * prC * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (prC * dwsr/ewvsr + 2 * ((eA * pa - eB * pb) * ewz/(wzdeno * ewusr)))
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx6 <- (sigx5 * eA - (db * ewvsr/ewusr - eC * pb) * eB)
  sigx7 <- (wzdeno * ewusr * (eA * pa - eB * pb)/(wzdeno * ewusr)^2)
  sigx8 <- (sigx6/(wzdeno * ewusr) - 0.5 * sigx7)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)))
  sigx12 <- ((0.5 * depsisr - 0.5 * dwsr) * prC/ewvsr + 2 * ((eA * sigx10 - sigx11 *
    eB) * ewz/(wzdeno * ewusr)))
  sigx13 <- ((1/(wzdeno * ewusr) - ewusr * ewz/(wzdeno * ewusr)^2) * (eA * pa -
    eB * pb))
  sigx14 <- (2 * sigx13 - prC * dwsr/(wzdeno * ewvsr))
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (prC * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 + 2 * (((((ewvsr/ewusr +
    S * (epsilon)/ewvsr)/ewvsr - 1/ewusr) * da/ewvsr - (da/ewvsr - pa/ewusr)/ewusr) *
    eA - (((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr)/ewvsr - 2/ewusr) * db/ewvsr -
    2 * ((db/ewvsr - 2 * (pb/ewusr))/ewusr)) * eB) * ewz/(wzdeno * ewusr)) -
    sigx2^2/sigx3)/sigx3, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv/(2 * ewu)^2)) * pa - 0.5 * (da * ewvsr/ewusr))/ewusr + (0.5 * ((ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewusr) - sigx4/ewvsr) * da) * eA - (((2 * (ewvsr/ewusr) +
      S * (epsilon)/ewvsr)/ewusr - eC/ewvsr) * db + (pb - 2 * (db * ewvsr/ewusr -
      eC * pb))/ewusr) * eB)/(wzdeno * ewusr) - (sigx8 * sigx2/sigx3 + 0.5 *
      (sigx1 * wzdeno * ewusr/(wzdeno * ewusr)^2))) * ewz/sigx3), FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((da * (ewv/(2 * ewu) - (sigx9 * (ewvsr/ewusr +
      S * (epsilon)/ewvsr) + 0.5))/ewvsr - sigx10/ewusr) * eA - ((2 * (ewv/ewu) -
      ((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr - 0.5 * (S *
        (epsilon)/ewvsr)) + 0.5)) * db/ewvsr - 2 * (sigx11/ewusr)) * eB) *
      ewz/(wzdeno * ewusr)) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) -
      0.5) * prC * dwsr * (epsilon)/ewvsr^3 - sigx12 * sigx2/sigx3)/sigx3,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * (sigx1 *
    (1/(wzdeno * ewusr) - ewusr * ewz/(wzdeno * ewusr)^2)) - (sigx2 * sigx14/sigx3 +
    S * prC * dwsr * (epsilon)/(wzdeno * ewvsr^3))) * ewz/sigx3, FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewvsr * (ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewusr) - 0.5) - 0.5 * sigx4) * da * ewvsr/ewusr -
      (sigx5 * sigx4 + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv/(2 *
        ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pa)) * eA - ((((2 * (ewvsr/ewusr) +
      S * (epsilon)/ewvsr) * ewvsr - S * (epsilon))/ewusr - (0.5 + 2 * (ewv/ewu))) *
      db * ewvsr/ewusr + (0.5 * (S * (epsilon)/ewusr) + 2 * (ewv/ewu)) * pb -
      eC * (db * ewvsr/ewusr - eC * pb)) * eB)/(wzdeno * ewusr) - ((0.5 * ((0.5 -
      wzdeno^2 * ewusr^2/(wzdeno * ewusr)^2) * (eA * pa - eB * pb) + sigx5 *
      eA - (db * ewvsr/ewusr - eC * pb) * eB) + 0.5 * sigx6) * wzdeno * ewusr/(wzdeno *
      ewusr)^2 + 2 * (sigx8^2 * ewz/sigx3))) * ewz/sigx3), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (2 * ((((da *
    ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pa/(2 * ewu)^2)) * ewv - ((0.5 * (sigx9 *
    (ewvsr/ewusr + S * (epsilon)/ewvsr)) - 0.25) * da * ewvsr/ewusr + sigx4 *
    sigx10)) * eA - (2 * ((db * ewvsr/ewusr - pb) * ewv/ewu) - (((2 * (ewvsr/ewusr) +
    S * (epsilon)/ewvsr) * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)) - 0.5) *
    db * ewvsr/ewusr + sigx11 * eC)) * eB)/(wzdeno * ewusr) - 0.5 * (wzdeno *
    ewusr * (eA * sigx10 - sigx11 * eB)/(wzdeno * ewusr)^2)) - 2 * (sigx8 * sigx12/sigx3)) *
    ewz/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (2 * (sigx6 * (1/(wzdeno * ewusr) - ewusr * ewz/(wzdeno * ewusr)^2) - ((0.5 -
      wzdeno^2 * ewusr^2/(wzdeno * ewusr)^2) * ewz + 0.5 * wzdeno) * ewusr *
      (eA * pa - eB * pb)/(wzdeno * ewusr)^2) - 2 * (sigx8 * sigx14 * ewz/sigx3)) *
    ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr *
      (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * depsisr - 0.5 * dwsr))/ewvsr + 2 *
      ((((sigx10/2 + (pa - sigx9 * da)/2) * ewv/ewu - (0.25 * (ewvsr/ewusr) +
        0.25 * (S * (epsilon)/ewvsr) - sigx9^2 * (ewvsr/ewusr + S * (epsilon)/ewvsr)) *
        da) * eA - ((2 * sigx11 + 2 * (pb - db * (ewvsr/ewusr - 0.5 * (S *
        (epsilon)/ewvsr)))) * ewv/ewu - (0.25 * (S * (epsilon)/ewvsr) + 0.5 *
        (ewvsr/ewusr) - (2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr -
        0.5 * (S * (epsilon)/ewvsr))^2) * db) * eB) * ewz/(wzdeno * ewusr)) -
      sigx12^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((1/(wzdeno * ewusr) - ewusr * ewz/(wzdeno *
      ewusr)^2) * (eA * sigx10 - sigx11 * eB)) - (sigx12 * sigx14/sigx3 + (0.5 *
      (S^2 * dwsr * (epsilon)^2/(wzdeno * ewvsr^3)) - 0.5 * (wzdeno * dwsr *
      ewvsr/(wzdeno * ewvsr)^2)) * prC)) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewvsr) + ewvsr/(wzdeno *
      ewvsr)^2) * dwsr - (sigx14^2/sigx3 + 2 * ((2 - 2 * (wzdeno * ewusr^2 *
      ewz/(wzdeno * ewusr)^2)) * ewusr * (eA * pa - eB * pb)/(wzdeno * ewusr)^2))) *
      ewz + 2 * sigx13 - prC * dwsr/(wzdeno * ewvsr)) * ewz/sigx3, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesszisfgenexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) * eB)
  sigx2 <- (2 * (sigx1 * ewz1/ewusr) + S * ewz2 * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (ewz2 * dwsr/ewvsr + 2 * (ewz1 * (eA * pa - eB * pb)/ewusr))
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx6 <- (sigx5 * eA - (db * ewvsr/ewusr - eC * pb) * eB)
  sigx7 <- (sigx5 * eA - ((db * ewvsr/ewusr - eC * pb) * eB + 0.5 * (eA * pa -
    eB * pb)))
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)))
  sigx8 <- (ewz1 * (eA * sigx10 - sigx11 * eB)/ewusr)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx12 <- (ewz2 * (0.5 * depsisr - 0.5 * dwsr)/ewvsr + 2 * sigx8)
  sigx13 <- (2 * ((eA * pa - eB * pb)/ewusr) - dwsr/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (ewz2 * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 + 2 * (((((ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewvsr - 1/ewusr) * da/ewvsr - (da/ewvsr - pa/ewusr)/ewusr) *
      eA - (((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr)/ewvsr - 2/ewusr) * db/ewvsr -
      2 * ((db/ewvsr - 2 * (pb/ewusr))/ewusr)) * eB) * ewz1/ewusr) - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv/(2 * ewu)^2)) * pa - 0.5 * (da * ewvsr/ewusr))/ewusr + (0.5 * ((ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewusr) - sigx4/ewvsr) * da) * eA - ((((2 * (ewvsr/ewusr) +
      S * (epsilon)/ewvsr)/ewusr - eC/ewvsr) * db + (pb - 2 * (db * ewvsr/ewusr -
      eC * pb))/ewusr) * eB + 0.5 * sigx1))/(sigx3 * ewusr) - sigx7 * sigx2 *
      ewusr/(sigx3 * ewusr)^2) * ewz1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((da * (ewv/(2 * ewu) - (sigx9 * (ewvsr/ewusr +
      S * (epsilon)/ewvsr) + 0.5))/ewvsr - sigx10/ewusr) * eA - ((2 * (ewv/ewu) -
      ((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr - 0.5 * (S *
        (epsilon)/ewvsr)) + 0.5)) * db/ewvsr - 2 * (sigx11/ewusr)) * eB) *
      ewz1/ewusr) + S * ewz2 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) *
      dwsr * (epsilon)/ewvsr^3 - sigx12 * sigx2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((2 * (sigx1/ewusr) -
    S * dwsr * (epsilon)/ewvsr^3)/(pi * sigx3 * ((Wz)^2 + 1)) - pi * ((Wz)^2 +
    1) * sigx2 * sigx13/(pi * sigx3 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewvsr * (ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewusr) - 0.5) - 0.5 * sigx4) * da * ewvsr/ewusr -
      (sigx5 * sigx4 + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv/(2 *
        ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pa)) * eA - (((((2 * (ewvsr/ewusr) +
      S * (epsilon)/ewvsr) * ewvsr - S * (epsilon))/ewusr - (0.5 + 2 * (ewv/ewu))) *
      db * ewvsr/ewusr + (0.5 * (S * (epsilon)/ewusr) + 2 * (ewv/ewu)) * pb -
      eC * (db * ewvsr/ewusr - eC * pb)) * eB + 0.5 * sigx6))/(sigx3 * ewusr) -
      sigx7 * (0.5 * (sigx3 * ewusr) + 2 * (sigx7 * ewz1))/(sigx3 * ewusr)^2) *
      ewz1), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ewz1 * (2 *
    (((da * ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pa/(2 * ewu)^2)) * ewv - ((0.5 *
      (sigx9 * (ewvsr/ewusr + S * (epsilon)/ewvsr)) - 0.25) * da * ewvsr/ewusr +
      sigx4 * sigx10)) * eA - ((2 * ((db * ewvsr/ewusr - pb) * ewv/ewu) - (((2 *
      (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)) -
      0.5) * db * ewvsr/ewusr + sigx11 * eC)) * eB + 0.5 * (eA * sigx10 - sigx11 *
      eB))) - 2 * (sigx12 * sigx7/sigx3))/(sigx3 * ewusr), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx7 * (2/(pi * sigx3 * ((Wz)^2 + 1)) - 2 * (pi * ((Wz)^2 + 1) * ewz1 *
    sigx13/(pi * sigx3 * ((Wz)^2 + 1))^2))/ewusr, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (ewz2 * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr *
      (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * depsisr - 0.5 * dwsr))/ewvsr + 2 *
      ((((sigx10/2 + (pa - sigx9 * da)/2) * ewv/ewu - (0.25 * (ewvsr/ewusr) +
        0.25 * (S * (epsilon)/ewvsr) - sigx9^2 * (ewvsr/ewusr + S * (epsilon)/ewvsr)) *
        da) * eA - ((2 * sigx11 + 2 * (pb - db * (ewvsr/ewusr - 0.5 * (S *
        (epsilon)/ewvsr)))) * ewv/ewu - (0.25 * (S * (epsilon)/ewvsr) + 0.5 *
        (ewvsr/ewusr) - (2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr -
        0.5 * (S * (epsilon)/ewvsr))^2) * db) * eB) * ewz1/ewusr) - sigx12^2/sigx3)/sigx3,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * ((eA * sigx10 - sigx11 * eB)/ewusr) - (0.5 *
      depsisr - 0.5 * dwsr)/ewvsr)/(pi * sigx3 * ((Wz)^2 + 1)) - pi * sigx12 *
      ((Wz)^2 + 1) * sigx13/(pi * sigx3 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (sigx13 * (2 * ((eA * pa - eB * pb)/ewusr) +
      2 * (pi * Wz * sigx3) - dwsr/ewvsr)/(pi * sigx3 * ((Wz)^2 + 1))^2), FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesszisfgenexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) * eB)
  sigx2 <- (2 * (sigx1 * pwZ/ewusr) + S * (1 - pwZ) * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - pwZ) * dwsr/ewvsr + 2 * ((eA * pa - eB * pb) * pwZ/ewusr))
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx6 <- (sigx5 * eA - (db * ewvsr/ewusr - eC * pb) * eB)
  sigx7 <- (sigx5 * eA - ((db * ewvsr/ewusr - eC * pb) * eB + 0.5 * (eA * pa -
    eB * pb)))
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)))
  sigx8 <- ((eA * sigx10 - sigx11 * eB) * pwZ/ewusr)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx12 <- ((0.5 * depsisr - 0.5 * dwsr) * (1 - pwZ)/ewvsr + 2 * sigx8)
  sigx13 <- (2 * ((eA * pa - eB * pb)/ewusr) - dwsr/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((1 - pwZ) * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 + 2 * (((((ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewvsr - 1/ewusr) * da/ewvsr - (da/ewvsr - pa/ewusr)/ewusr) *
      eA - (((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr)/ewvsr - 2/ewusr) * db/ewvsr -
      2 * ((db/ewvsr - 2 * (pb/ewusr))/ewusr)) * eB) * pwZ/ewusr) - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv/(2 * ewu)^2)) * pa - 0.5 * (da * ewvsr/ewusr))/ewusr + (0.5 * ((ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewusr) - sigx4/ewvsr) * da) * eA - ((((2 * (ewvsr/ewusr) +
      S * (epsilon)/ewvsr)/ewusr - eC/ewvsr) * db + (pb - 2 * (db * ewvsr/ewusr -
      eC * pb))/ewusr) * eB + 0.5 * sigx1))/(sigx3 * ewusr) - sigx7 * sigx2 *
      ewusr/(sigx3 * ewusr)^2) * pwZ), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((da * (ewv/(2 * ewu) - (sigx9 * (ewvsr/ewusr +
      S * (epsilon)/ewvsr) + 0.5))/ewvsr - sigx10/ewusr) * eA - ((2 * (ewv/ewu) -
      ((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr - 0.5 * (S *
        (epsilon)/ewvsr)) + 0.5)) * db/ewvsr - 2 * (sigx11/ewusr)) * eB) *
      pwZ/ewusr) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * (1 -
      pwZ) * dwsr * (epsilon)/ewvsr^3 - sigx12 * sigx2/sigx3)/sigx3, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * (sigx1/ewusr) -
    (sigx2 * sigx13/sigx3 + S * dwsr * (epsilon)/ewvsr^3)) * dwZ/sigx3, FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewvsr * (ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewusr) - 0.5) - 0.5 * sigx4) * da * ewvsr/ewusr -
      (sigx5 * sigx4 + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv/(2 *
        ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pa)) * eA - (((((2 * (ewvsr/ewusr) +
      S * (epsilon)/ewvsr) * ewvsr - S * (epsilon))/ewusr - (0.5 + 2 * (ewv/ewu))) *
      db * ewvsr/ewusr + (0.5 * (S * (epsilon)/ewusr) + 2 * (ewv/ewu)) * pb -
      eC * (db * ewvsr/ewusr - eC * pb)) * eB + 0.5 * sigx6))/(sigx3 * ewusr) -
      sigx7 * (0.5 * (sigx3 * ewusr) + 2 * (sigx7 * pwZ))/(sigx3 * ewusr)^2) *
      pwZ), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (2 * (((da *
    ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pa/(2 * ewu)^2)) * ewv - ((0.5 * (sigx9 *
    (ewvsr/ewusr + S * (epsilon)/ewvsr)) - 0.25) * da * ewvsr/ewusr + sigx4 *
    sigx10)) * eA - ((2 * ((db * ewvsr/ewusr - pb) * ewv/ewu) - (((2 * (ewvsr/ewusr) +
    S * (epsilon)/ewvsr) * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)) - 0.5) *
    db * ewvsr/ewusr + sigx11 * eC)) * eB + 0.5 * (eA * sigx10 - sigx11 * eB))) -
    2 * (sigx7 * sigx12/sigx3)) * pwZ/(sigx3 * ewusr), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx7 * (2 - 2 * (sigx13 * pwZ/sigx3)) * dwZ/(sigx3 * ewusr), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((1 - pwZ) * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) *
      dwsr * (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * depsisr - 0.5 * dwsr))/ewvsr +
      2 * ((((sigx10/2 + (pa - sigx9 * da)/2) * ewv/ewu - (0.25 * (ewvsr/ewusr) +
        0.25 * (S * (epsilon)/ewvsr) - sigx9^2 * (ewvsr/ewusr + S * (epsilon)/ewvsr)) *
        da) * eA - ((2 * sigx11 + 2 * (pb - db * (ewvsr/ewusr - 0.5 * (S *
        (epsilon)/ewvsr)))) * ewv/ewu - (0.25 * (S * (epsilon)/ewvsr) + 0.5 *
        (ewvsr/ewusr) - (2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr -
        0.5 * (S * (epsilon)/ewvsr))^2) * db) * eB) * pwZ/ewusr) - sigx12^2/sigx3)/sigx3,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((eA * sigx10 - sigx11 * eB)/ewusr) - (sigx12 *
      sigx13/sigx3 + (0.5 * depsisr - 0.5 * dwsr)/ewvsr)) * dwZ/sigx3, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx13 * dwZ/sigx3 + Wz) * sigx13 * dwZ/sigx3),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesszisfgenexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) * eB)
  sigx2 <- (2 * (sigx1 * (1 - prZ)/ewusr) + S * dwsr * prZ * (epsilon)/ewvsr^3)
  sigx3 <- (2 * ((1 - prZ) * (eA * pa - eB * pb)/ewusr) + dwsr * prZ/ewvsr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx6 <- (sigx5 * eA - (db * ewvsr/ewusr - eC * pb) * eB)
  sigx7 <- (sigx5 * eA - ((db * ewvsr/ewusr - eC * pb) * eB + 0.5 * (eA * pa -
    eB * pb)))
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)))
  sigx8 <- ((1 - prZ) * (eA * sigx10 - sigx11 * eB)/ewusr)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx12 <- ((0.5 * depsisr - 0.5 * dwsr) * prZ/ewvsr + 2 * sigx8)
  sigx13 <- (2 * ((eA * pa - eB * pb)/ewusr) - dwsr/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (2 * (((((ewvsr/ewusr + S * (epsilon)/ewvsr)/ewvsr - 1/ewusr) * da/ewvsr -
      (da/ewvsr - pa/ewusr)/ewusr) * eA - (((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr)/ewvsr -
      2/ewusr) * db/ewvsr - 2 * ((db/ewvsr - 2 * (pb/ewusr))/ewusr)) * eB) *
      (1 - prZ)/ewusr) + dwsr * prZ * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 -
      sigx2^2/sigx3)/sigx3, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv/(2 * ewu)^2)) * pa - 0.5 * (da * ewvsr/ewusr))/ewusr + (0.5 * ((ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewusr) - sigx4/ewvsr) * da) * eA - ((((2 * (ewvsr/ewusr) +
      S * (epsilon)/ewvsr)/ewusr - eC/ewvsr) * db + (pb - 2 * (db * ewvsr/ewusr -
      eC * pb))/ewusr) * eB + 0.5 * sigx1))/(sigx3 * ewusr) - sigx7 * sigx2 *
      ewusr/(sigx3 * ewusr)^2) * (1 - prZ)), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((da * (ewv/(2 * ewu) - (sigx9 * (ewvsr/ewusr +
      S * (epsilon)/ewvsr) + 0.5))/ewvsr - sigx10/ewusr) * eA - ((2 * (ewv/ewu) -
      ((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr - 0.5 * (S *
        (epsilon)/ewvsr)) + 0.5)) * db/ewvsr - 2 * (sigx11/ewusr)) * eB) *
      (1 - prZ)/ewusr) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) *
      dwsr * prZ * (epsilon)/ewvsr^3 - sigx12 * sigx2/sigx3)/sigx3, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * (sigx1/ewusr) -
    (sigx2 * sigx13/sigx3 + S * dwsr * (epsilon)/ewvsr^3)) * prZ * ewz/sigx3,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewvsr * (ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewusr) - 0.5) - 0.5 * sigx4) * da * ewvsr/ewusr -
      (sigx5 * sigx4 + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv/(2 *
        ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pa)) * eA - (((((2 * (ewvsr/ewusr) +
      S * (epsilon)/ewvsr) * ewvsr - S * (epsilon))/ewusr - (0.5 + 2 * (ewv/ewu))) *
      db * ewvsr/ewusr + (0.5 * (S * (epsilon)/ewusr) + 2 * (ewv/ewu)) * pb -
      eC * (db * ewvsr/ewusr - eC * pb)) * eB + 0.5 * sigx6))/(sigx3 * ewusr) -
      sigx7 * (0.5 * (sigx3 * ewusr) + 2 * (sigx7 * (1 - prZ)))/(sigx3 * ewusr)^2) *
      (1 - prZ)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (1 - prZ) *
    (2 * (((da * ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pa/(2 * ewu)^2)) * ewv -
      ((0.5 * (sigx9 * (ewvsr/ewusr + S * (epsilon)/ewvsr)) - 0.25) * da *
        ewvsr/ewusr + sigx4 * sigx10)) * eA - ((2 * ((db * ewvsr/ewusr -
      pb) * ewv/ewu) - (((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr -
      0.5 * (S * (epsilon)/ewvsr)) - 0.5) * db * ewvsr/ewusr + sigx11 * eC)) *
      eB + 0.5 * (eA * sigx10 - sigx11 * eB))) - 2 * (sigx7 * sigx12/sigx3))/(sigx3 *
    ewusr), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx7 * (2 - 2 * ((1 - prZ) * sigx13/sigx3)) * prZ * ewz/(sigx3 * ewusr),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (2 * ((((sigx10/2 + (pa - sigx9 * da)/2) * ewv/ewu - (0.25 * (ewvsr/ewusr) +
      0.25 * (S * (epsilon)/ewvsr) - sigx9^2 * (ewvsr/ewusr + S * (epsilon)/ewvsr)) *
      da) * eA - ((2 * sigx11 + 2 * (pb - db * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)))) *
      ewv/ewu - (0.25 * (S * (epsilon)/ewvsr) + 0.5 * (ewvsr/ewusr) - (2 *
      (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr))^2) *
      db) * eB) * (1 - prZ)/ewusr) + prZ * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
      1) - 0.25) * dwsr * (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * depsisr - 0.5 *
      dwsr))/ewvsr - sigx12^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((eA * sigx10 - sigx11 * eB)/ewusr) - (sigx12 *
      sigx13/sigx3 + (0.5 * depsisr - 0.5 * dwsr)/ewvsr)) * prZ * ewz/sigx3,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (sigx13 * prZ/sigx3 + 1) * ewz) * sigx13 *
      prZ * ewz/sigx3, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Different sigma_v

## logit specification class membership
chessmnsfgenexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewu_h <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  mustar1 <- (ewv1_h/ewu_h + S * epsilon/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(S * epsilon/ewv2_h, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  eusq1 <- exp(ewv1/(2 * ewu) + S * epsilon/ewu_h)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * epsilon/ewu_h))
  evsq1 <- (2 * (ewv1_h/ewu_h) + S * epsilon/ewv1_h)
  depsi <- dnorm(-evsq1, 0, 1)
  pepsi <- pnorm(-evsq1)
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewu_h)
  dpepsi <- (depsi/ewv1_h - 2 * (pepsi/ewu_h))
  epepe <- (eusq1 * pmustar1 - eusq2 * pepsi)
  exusq <- (sigx1 * eusq1 - dpepsi * eusq2)
  ewvu <- (ewu * ewv1/(2 * ewu)^2)
  sigx2 <- (exusq * ewz/(wzdeno * ewu_h))
  sigx3 <- (2 * sigx2 + S * prC * dmustar2 * epsilon/ewv2_h^3)
  sigx4 <- (epepe * ewz/(wzdeno * ewu_h))
  sigx5 <- (prC * dmustar2/ewv2_h + 2 * sigx4)
  sigx6 <- (0.5 * (S * epsilon/ewu_h) + 2 * ewvu)
  sigx7 <- (0.5 * (dmustar1 * ewv1_h/ewu_h) - sigx6 * pmustar1)
  sigx8 <- (depsi * ewv1_h/ewu_h - (2 * (ewv1/ewu) + S * epsilon/ewu_h) * pepsi)
  sigx9 <- (sigx7 * eusq1 - sigx8 * eusq2)
  sigx10 <- (wzdeno * ewu_h * epepe/(wzdeno * ewu_h)^2)
  sigx11 <- (sigx9/(wzdeno * ewu_h) - 0.5 * sigx10)
  sigx12 <- (ewv1 * pmustar1/(2 * ewu) - (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h)) *
    dmustar1)
  sigx13 <- (2 * (ewv1 * pepsi/ewu) - depsi * (ewv1_h/ewu_h - 0.5 * (S * epsilon/ewv1_h)))
  sigx14 <- (eusq1 * sigx12 - sigx13 * eusq2)
  sigx15 <- (0.5 * (S^2 * dmustar2 * epsilon^2/ewv2_h^2) - 0.5 * dmustar2)
  sigx16 <- (1/(wzdeno * ewu_h) - ewu_h * ewz/(wzdeno * ewu_h)^2)
  sigx17 <- (sigx16 * epepe)
  sigx18 <- (2 * sigx17 - prC * dmustar2/(wzdeno * ewv2_h))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (prC * dmustar2 * (S^2 * epsilon^2/ewv2_h^2 - 1)/ewv2_h^3 + 2 * ((((mustar1/ewv1_h -
    1/ewu_h) * dmustar1/ewv1_h - sigx1/ewu_h) * eusq1 - ((evsq1/ewv1_h - 2/ewu_h) *
    depsi/ewv1_h - 2 * (dpepsi/ewu_h)) * eusq2) * ewz/(wzdeno * ewu_h)) - sigx3^2/sigx5)/sigx5,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * epsilon/ewu_h) + 2 * ewvu) *
      pmustar1 - 0.5 * (dmustar1 * ewv1_h/ewu_h))/ewu_h + (0.5 * (mustar1/ewu_h) -
      sigx6/ewv1_h) * dmustar1) * eusq1 - ((evsq1/ewu_h - (2 * (ewv1/ewu) +
      S * epsilon/ewu_h)/ewv1_h) * depsi + (pepsi - 2 * sigx8)/ewu_h) * eusq2)/(wzdeno *
      ewu_h) - (sigx11 * sigx3/sigx5 + 0.5 * (exusq * wzdeno * ewu_h/(wzdeno *
      ewu_h)^2))) * ewz/sigx5), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((dmustar1 * (ewv1/(2 * ewu) - ((0.5 *
      (ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h)) * mustar1 + 0.5))/ewv1_h -
      sigx12/ewu_h) * eusq1 - ((2 * (ewv1/ewu) - (evsq1 * (ewv1_h/ewu_h - 0.5 *
      (S * epsilon/ewv1_h)) + 0.5)) * depsi/ewv1_h - 2 * (sigx13/ewu_h)) *
      eusq2)/(sigx5 * wzdeno * ewu_h) - wzdeno * sigx3 * ewu_h * sigx14/(sigx5 *
      wzdeno * ewu_h)^2) * ewz), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * prC * (S * (0.5 * (S^2 * epsilon^2/ewv2_h^2 -
      2) - 0.5) * dmustar2 * epsilon/(sigx5 * ewv2_h^3) - sigx15 * sigx3 *
      ewv2_h/(sigx5 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * (exusq *
    sigx16) - (sigx3 * sigx18/sigx5 + S * prC * dmustar2 * epsilon/(wzdeno *
    ewv2_h^3))) * ewz/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewv1_h * mustar1/ewu_h) -
      0.5) - 0.5 * sigx6) * dmustar1 * ewv1_h/ewu_h - (sigx7 * sigx6 + (2 *
      ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S *
      epsilon/ewu_h)) * pmustar1)) * eusq1 - (((evsq1 * ewv1_h - S * epsilon)/ewu_h -
      (0.5 + 2 * (ewv1/ewu))) * depsi * ewv1_h/ewu_h + (0.5 * (S * epsilon/ewu_h) +
      2 * (ewv1/ewu)) * pepsi - (2 * (ewv1/ewu) + S * epsilon/ewu_h) * sigx8) *
      eusq2)/(wzdeno * ewu_h) - ((0.5 * ((0.5 - wzdeno^2 * ewu_h^2/(wzdeno *
      ewu_h)^2) * epepe + sigx7 * eusq1 - sigx8 * eusq2) + 0.5 * sigx9) * wzdeno *
      ewu_h/(wzdeno * ewu_h)^2 + 2 * (sigx11^2 * ewz/sigx5))) * ewz/sigx5),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((((dmustar1 *
    ewv1_h/(4 * (ewu * ewu_h)) - 2 * (ewu * pmustar1/(2 * ewu)^2)) * ewv1 - ((0.5 *
    ((0.5 * (ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h)) * mustar1) - 0.25) *
    dmustar1 * ewv1_h/ewu_h + sigx6 * sigx12)) * eusq1 - (2 * ((depsi * ewv1_h/ewu_h -
    pepsi) * ewv1/ewu) - ((evsq1 * (ewv1_h/ewu_h - 0.5 * (S * epsilon/ewv1_h)) -
    0.5) * depsi * ewv1_h/ewu_h + sigx13 * (2 * (ewv1/ewu) + S * epsilon/ewu_h))) *
    eusq2)/(sigx5 * wzdeno * ewu_h) - (0.5 * sigx5 + 2 * (sigx11 * ewz)) * wzdeno *
    ewu_h * sigx14/(sigx5 * wzdeno * ewu_h)^2) * ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (2 * (sigx11 * sigx15 * prC * ewv2_h * ewz/(sigx5 * ewv2_h)^2)), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (2 * (sigx9 * sigx16 - ((0.5 - wzdeno^2 * ewu_h^2/(wzdeno * ewu_h)^2) * ewz +
      0.5 * wzdeno) * ewu_h * epepe/(wzdeno * ewu_h)^2) - 2 * (sigx11 * sigx18 *
      ewz/sigx5)) * ewz/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    2 * (((((sigx12/2 + (pmustar1 - (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h)) *
    dmustar1)/2) * ewv1/ewu - (0.25 * (ewv1_h/ewu_h) + 0.25 * (S * epsilon/ewv1_h) -
    (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h))^2 * mustar1) * dmustar1) *
    eusq1 - ((2 * sigx13 + 2 * (pepsi - depsi * (ewv1_h/ewu_h - 0.5 * (S * epsilon/ewv1_h)))) *
    ewv1/ewu - (0.25 * (S * epsilon/ewv1_h) + 0.5 * (ewv1_h/ewu_h) - evsq1 *
    (ewv1_h/ewu_h - 0.5 * (S * epsilon/ewv1_h))^2) * depsi) * eusq2)/(sigx5 *
    wzdeno * ewu_h) - 2 * (sigx14^2 * ewz/(sigx5 * wzdeno * ewu_h)^2)) * ewz),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (2 * (sigx15 * prC * sigx14 * ewv2_h * ewz/((sigx5 * ewv2_h)^2 *
      wzdeno * ewu_h))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * sigx16 - 2 * (sigx18 * ewz/(sigx5 * wzdeno *
      ewu_h))) * sigx14 * ewz/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (S^2 * (0.5 * (0.5 * (S^2 * epsilon^2/ewv2_h^2) -
      1) - 0.25) * dmustar2 * epsilon^2/(sigx5 * ewv2_h^3) - (sigx15 * prC +
      0.5 * (sigx5 * ewv2_h)) * sigx15/(sigx5 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx15 * sigx18/sigx5 + 0.5 * (S^2 * dmustar2 *
      epsilon^2/(wzdeno * ewv2_h^2)))/ewv2_h - 0.5 * (wzdeno * dmustar2 * ewv2_h/(wzdeno *
      ewv2_h)^2)) * prC * ewz/sigx5), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv2_h) + ewv2_h/(wzdeno *
      ewv2_h)^2) * dmustar2 - (sigx18^2/sigx5 + 2 * ((2 - 2 * (wzdeno * ewu_h^2 *
      ewz/(wzdeno * ewu_h)^2)) * ewu_h * epepe/(wzdeno * ewu_h)^2))) * ewz +
      2 * sigx17 - prC * dmustar2/(wzdeno * ewv2_h)) * ewz/sigx5, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmnsfgenexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  mustar1 <- (ewv1_h/ewusr + S * (epsilon)/ewv1_h)
  mustar2 <- (2 * (ewv1_h/ewusr) + S * (epsilon)/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(-mustar2, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  pmustar2 <- pnorm(-mustar2)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  eusq1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * (epsilon)/ewusr))
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewusr)
  sigx2 <- (dmustar2/ewv1_h - 2 * (pmustar2/ewusr))
  sigx3 <- ((sigx1 * eusq1 - sigx2 * eusq2) * ewz1/ewusr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar1 * ewv1_h/ewusr) - sigx4 * pmustar1)
  sigx6 <- (dmustar2 * ewv1_h/ewusr - (2 * (ewv1/ewu) + S * (epsilon)/ewusr) *
    pmustar2)
  sigx7 <- (sigx6 * eusq2 + 0.5 * (eusq1 * pmustar1 - eusq2 * pmustar2))
  sigx8 <- (ewz2 * depsi/ewv2_h + 2 * (ewz1 * (eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr))
  sigx9 <- (2 * sigx3 + S * ewz2 * depsi * (epsilon)/ewv2_h^3)
  sigx10 <- (0.5 * (ewv1_h/ewusr) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx11 <- (ewv1 * pmustar1/(2 * ewu) - sigx10 * dmustar1)
  sigx12 <- (ewv1_h/ewusr - 0.5 * (S * (epsilon)/ewv1_h))
  sigx13 <- (eusq1 * sigx11 - (2 * (ewv1 * pmustar2/ewu) - dmustar2 * sigx12) *
    eusq2)
  sigx14 <- (2 * ((eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr) - depsi/ewv2_h)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (ewz2 * depsi * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 + 2 * ((((mustar1/ewv1_h -
      1/ewusr) * dmustar1/ewv1_h - sigx1/ewusr) * eusq1 - ((mustar2/ewv1_h -
      2/ewusr) * dmustar2/ewv1_h - 2 * (sigx2/ewusr)) * eusq2) * ewz1/ewusr) -
      sigx9^2/sigx8)/sigx8, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv1/(2 * ewu)^2)) * pmustar1 - 0.5 * (dmustar1 * ewv1_h/ewusr))/ewusr +
      (0.5 * (mustar1/ewusr) - sigx4/ewv1_h) * dmustar1) * eusq1 - (((mustar2/ewusr -
      (2 * (ewv1/ewu) + S * (epsilon)/ewusr)/ewv1_h) * dmustar2 + (pmustar2 -
      2 * sigx6)/ewusr) * eusq2 + 0.5 * (sigx1 * eusq1 - sigx2 * eusq2)))/(sigx8 *
      ewusr) - (sigx5 * eusq1 - sigx7) * sigx9 * ewusr/(sigx8 * ewusr)^2) *
      ewz1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((dmustar1 * (ewv1/(2 * ewu) - (sigx10 *
      mustar1 + 0.5))/ewv1_h - sigx11/ewusr) * eusq1 - ((2 * (ewv1/ewu) - (mustar2 *
      sigx12 + 0.5)) * dmustar2/ewv1_h - 2 * ((2 * (ewv1 * pmustar2/ewu) -
      dmustar2 * sigx12)/ewusr)) * eusq2)/(sigx8 * ewusr) - sigx9 * ewusr *
      sigx13/(sigx8 * ewusr)^2) * ewz1), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ewz2 * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx8 * ewv2_h^3) - sigx24 * sigx9 * ewv2_h/(sigx8 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((2 * ((sigx1 *
    eusq1 - sigx2 * eusq2)/ewusr) - S * depsi * (epsilon)/ewv2_h^3)/(pi * sigx8 *
    ((Wz)^2 + 1)) - pi * ((Wz)^2 + 1) * sigx9 * sigx14/(pi * sigx8 * ((Wz)^2 +
    1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewv1_h * mustar1/ewusr) -
      0.5) - 0.5 * sigx4) * dmustar1 * ewv1_h/ewusr - (sigx5 * sigx4 + (2 *
      ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S *
      (epsilon)/ewusr)) * pmustar1)) * eusq1 - ((((mustar2 * ewv1_h - S * (epsilon))/ewusr -
      (0.5 + 2 * (ewv1/ewu))) * dmustar2 * ewv1_h/ewusr + (0.5 * (S * (epsilon)/ewusr) +
      2 * (ewv1/ewu)) * pmustar2 - (2 * (ewv1/ewu) + S * (epsilon)/ewusr) *
      sigx6) * eusq2 + 0.5 * (sigx5 * eusq1 - sigx6 * eusq2)))/(sigx8 * ewusr) -
      (sigx5 * eusq1 - sigx7) * (0.5 * (sigx8 * ewusr) + 2 * ((sigx5 * eusq1 -
        sigx7) * ewz1))/(sigx8 * ewusr)^2) * ewz1), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((((dmustar1 *
    ewv1_h/(4 * (ewu * ewusr)) - 2 * (ewu * pmustar1/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx10 * mustar1) - 0.25) * dmustar1 * ewv1_h/ewusr + sigx4 * sigx11)) *
    eusq1 - (2 * ((dmustar2 * ewv1_h/ewusr - pmustar2) * ewv1/ewu) - ((mustar2 *
    sigx12 - 0.5) * dmustar2 * ewv1_h/ewusr + (2 * (ewv1 * pmustar2/ewu) - dmustar2 *
    sigx12) * (2 * (ewv1/ewu) + S * (epsilon)/ewusr))) * eusq2)/(sigx8 * ewusr) -
    (0.5 * (sigx8 * ewusr) + 2 * ((sigx5 * eusq1 - sigx7) * ewz1)) * sigx13/(sigx8 *
      ewusr)^2) * ewz1), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (2 * ((sigx5 * eusq1 - sigx7) * ewz2 * sigx24 * ewz1 * ewv2_h/((sigx8 * ewv2_h)^2 *
      ewusr))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (sigx5 * eusq1 - sigx7) * (2/(pi * sigx8 * ((Wz)^2 + 1)) - 2 * (pi * ((Wz)^2 +
    1) * ewz1 * sigx14/(pi * sigx8 * ((Wz)^2 + 1))^2))/ewusr, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    2 * (((((sigx11/2 + (pmustar1 - sigx10 * dmustar1)/2) * ewv1/ewu - (0.25 *
    (ewv1_h/ewusr) + 0.25 * (S * (epsilon)/ewv1_h) - sigx10^2 * mustar1) * dmustar1) *
    eusq1 - ((2 * (2 * (ewv1 * pmustar2/ewu) - dmustar2 * sigx12) + 2 * (pmustar2 -
    dmustar2 * sigx12)) * ewv1/ewu - (0.25 * (S * (epsilon)/ewv1_h) + 0.5 * (ewv1_h/ewusr) -
    mustar2 * sigx12^2) * dmustar2) * eusq2)/(sigx8 * ewusr) - 2 * (ewz1 * sigx13^2/(sigx8 *
    ewusr)^2)) * ewz1), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (2 * (ewz2 * sigx24 * ewz1 * sigx13 * ewv2_h/((sigx8 * ewv2_h)^2 *
      ewusr))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2/(pi * sigx8 * ((Wz)^2 + 1)) - 2 * (pi * ((Wz)^2 +
      1) * ewz1 * sigx14/(pi * sigx8 * ((Wz)^2 + 1))^2)) * sigx13/ewusr, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ewz2 * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx8 * ewv2_h^3) - (ewz2 * sigx24 +
      0.5 * (sigx8 * ewv2_h)) * sigx24/(sigx8 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx24 * (1/(pi * sigx8 * ((Wz)^2 + 1)) + pi *
      ((Wz)^2 + 1) * ewz2 * sigx14/(pi * sigx8 * ((Wz)^2 + 1))^2)/ewv2_h),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (sigx14 * (2 * ((eusq1 * pmustar1 - eusq2 *
      pmustar2)/ewusr) + 2 * (pi * Wz * sigx8) - depsi/ewv2_h)/(pi * sigx8 *
      ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmnsfgenexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  mustar1 <- (ewv1_h/ewusr + S * (epsilon)/ewv1_h)
  mustar2 <- (2 * (ewv1_h/ewusr) + S * (epsilon)/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(-mustar2, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  pmustar2 <- pnorm(-mustar2)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  eusq1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * (epsilon)/ewusr))
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewusr)
  sigx2 <- (dmustar2/ewv1_h - 2 * (pmustar2/ewusr))
  sigx3 <- ((sigx1 * eusq1 - sigx2 * eusq2) * pwZ/ewusr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar1 * ewv1_h/ewusr) - sigx4 * pmustar1)
  sigx6 <- (dmustar2 * ewv1_h/ewusr - (2 * (ewv1/ewu) + S * (epsilon)/ewusr) *
    pmustar2)
  sigx7 <- (sigx6 * eusq2 + 0.5 * (eusq1 * pmustar1 - eusq2 * pmustar2))
  sigx9 <- (2 * sigx3 + S * (1 - pwZ) * depsi * (epsilon)/ewv2_h^3)
  sigx10 <- (0.5 * (ewv1_h/ewusr) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx11 <- (ewv1 * pmustar1/(2 * ewu) - sigx10 * dmustar1)
  sigx12 <- (ewv1_h/ewusr - 0.5 * (S * (epsilon)/ewv1_h))
  sigx13 <- (eusq1 * sigx11 - (2 * (ewv1 * pmustar2/ewu) - dmustar2 * sigx12) *
    eusq2)
  sigx14 <- (2 * ((eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr) - depsi/ewv2_h)
  sigx15 <- ((1 - pwZ) * depsi/ewv2_h + 2 * ((eusq1 * pmustar1 - eusq2 * pmustar2) *
    pwZ/ewusr))
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((1 - pwZ) * depsi * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 + 2 * ((((mustar1/ewv1_h -
      1/ewusr) * dmustar1/ewv1_h - sigx1/ewusr) * eusq1 - ((mustar2/ewv1_h -
      2/ewusr) * dmustar2/ewv1_h - 2 * (sigx2/ewusr)) * eusq2) * pwZ/ewusr) -
      sigx9^2/sigx15)/sigx15, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv1/(2 * ewu)^2)) * pmustar1 - 0.5 * (dmustar1 * ewv1_h/ewusr))/ewusr +
      (0.5 * (mustar1/ewusr) - sigx4/ewv1_h) * dmustar1) * eusq1 - (((mustar2/ewusr -
      (2 * (ewv1/ewu) + S * (epsilon)/ewusr)/ewv1_h) * dmustar2 + (pmustar2 -
      2 * sigx6)/ewusr) * eusq2 + 0.5 * (sigx1 * eusq1 - sigx2 * eusq2)))/(sigx15 *
      ewusr) - (sigx5 * eusq1 - sigx7) * sigx9 * ewusr/(sigx15 * ewusr)^2) *
      pwZ), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((dmustar1 * (ewv1/(2 * ewu) - (sigx10 *
      mustar1 + 0.5))/ewv1_h - sigx11/ewusr) * eusq1 - ((2 * (ewv1/ewu) - (mustar2 *
      sigx12 + 0.5)) * dmustar2/ewv1_h - 2 * ((2 * (ewv1 * pmustar2/ewu) -
      dmustar2 * sigx12)/ewusr)) * eusq2)/(sigx15 * ewusr) - sigx9 * ewusr *
      sigx13/(sigx15 * ewusr)^2) * pwZ), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (1 - pwZ) * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx15 * ewv2_h^3) - sigx24 * sigx9 *
      ewv2_h/(sigx15 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * ((sigx1 *
    eusq1 - sigx2 * eusq2)/ewusr) - (sigx9 * sigx14/sigx15 + S * depsi * (epsilon)/ewv2_h^3)) *
    dwZ/sigx15, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewv1_h * mustar1/ewusr) -
      0.5) - 0.5 * sigx4) * dmustar1 * ewv1_h/ewusr - (sigx5 * sigx4 + (2 *
      ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S *
      (epsilon)/ewusr)) * pmustar1)) * eusq1 - ((((mustar2 * ewv1_h - S * (epsilon))/ewusr -
      (0.5 + 2 * (ewv1/ewu))) * dmustar2 * ewv1_h/ewusr + (0.5 * (S * (epsilon)/ewusr) +
      2 * (ewv1/ewu)) * pmustar2 - (2 * (ewv1/ewu) + S * (epsilon)/ewusr) *
      sigx6) * eusq2 + 0.5 * (sigx5 * eusq1 - sigx6 * eusq2)))/(sigx15 * ewusr) -
      (sigx5 * eusq1 - sigx7) * (0.5 * (sigx15 * ewusr) + 2 * ((sigx5 * eusq1 -
        sigx7) * pwZ))/(sigx15 * ewusr)^2) * pwZ), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((((dmustar1 *
    ewv1_h/(4 * (ewu * ewusr)) - 2 * (ewu * pmustar1/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx10 * mustar1) - 0.25) * dmustar1 * ewv1_h/ewusr + sigx4 * sigx11)) *
    eusq1 - (2 * ((dmustar2 * ewv1_h/ewusr - pmustar2) * ewv1/ewu) - ((mustar2 *
    sigx12 - 0.5) * dmustar2 * ewv1_h/ewusr + (2 * (ewv1 * pmustar2/ewu) - dmustar2 *
    sigx12) * (2 * (ewv1/ewu) + S * (epsilon)/ewusr))) * eusq2)/(sigx15 * ewusr) -
    (0.5 * (sigx15 * ewusr) + 2 * ((sigx5 * eusq1 - sigx7) * pwZ)) * sigx13/(sigx15 *
      ewusr)^2) * pwZ), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (2 * ((sigx5 * eusq1 - sigx7) * sigx24 * (1 - pwZ) * ewv2_h * pwZ/((sigx15 *
      ewv2_h)^2 * ewusr))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (sigx5 * eusq1 - sigx7) * (2 - 2 * (sigx14 * pwZ/sigx15)) * dwZ/(sigx15 *
    ewusr), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    2 * (((((sigx11/2 + (pmustar1 - sigx10 * dmustar1)/2) * ewv1/ewu - (0.25 *
    (ewv1_h/ewusr) + 0.25 * (S * (epsilon)/ewv1_h) - sigx10^2 * mustar1) * dmustar1) *
    eusq1 - ((2 * (2 * (ewv1 * pmustar2/ewu) - dmustar2 * sigx12) + 2 * (pmustar2 -
    dmustar2 * sigx12)) * ewv1/ewu - (0.25 * (S * (epsilon)/ewv1_h) + 0.5 * (ewv1_h/ewusr) -
    mustar2 * sigx12^2) * dmustar2) * eusq2)/(sigx15 * ewusr) - 2 * (sigx13^2 *
    pwZ/(sigx15 * ewusr)^2)) * pwZ), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (2 * (sigx24 * (1 - pwZ) * sigx13 * ewv2_h * pwZ/((sigx15 *
      ewv2_h)^2 * ewusr))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 - 2 * (sigx14 * pwZ/sigx15)) * dwZ * sigx13/(sigx15 *
      ewusr), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - pwZ) * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx15 * ewv2_h^3) - (sigx24 * (1 -
      pwZ) + 0.5 * (sigx15 * ewv2_h)) * sigx24/(sigx15 * ewv2_h)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((1 - pwZ) * sigx14/sigx15 + 1) * sigx24 *
      dwZ/(sigx15 * ewv2_h)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx14 * dwZ/sigx15 + Wz) * sigx14 * dwZ/sigx15),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmnsfgenexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  mustar1 <- (ewv1_h/ewusr + S * (epsilon)/ewv1_h)
  mustar2 <- (2 * (ewv1_h/ewusr) + S * (epsilon)/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(-mustar2, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  pmustar2 <- pnorm(-mustar2)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  eusq1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * (epsilon)/ewusr))
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewusr)
  sigx2 <- (dmustar2/ewv1_h - 2 * (pmustar2/ewusr))
  sigx3 <- ((sigx1 * eusq1 - sigx2 * eusq2) * (1 - prZ)/ewusr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar1 * ewv1_h/ewusr) - sigx4 * pmustar1)
  sigx6 <- (dmustar2 * ewv1_h/ewusr - (2 * (ewv1/ewu) + S * (epsilon)/ewusr) *
    pmustar2)
  sigx7 <- (sigx6 * eusq2 + 0.5 * (eusq1 * pmustar1 - eusq2 * pmustar2))
  sigx8 <- (2 * ((1 - prZ) * (eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr) + depsi *
    prZ/ewv2_h)
  sigx9 <- (2 * sigx3 + S * depsi * prZ * (epsilon)/ewv2_h^3)
  sigx10 <- (0.5 * (ewv1_h/ewusr) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx11 <- (ewv1 * pmustar1/(2 * ewu) - sigx10 * dmustar1)
  sigx12 <- (ewv1_h/ewusr - 0.5 * (S * (epsilon)/ewv1_h))
  sigx13 <- (eusq1 * sigx11 - (2 * (ewv1 * pmustar2/ewu) - dmustar2 * sigx12) *
    eusq2)
  sigx14 <- (2 * ((eusq1 * pmustar1 - eusq2 * pmustar2)/ewusr) - depsi/ewv2_h)
  sigx24 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 * depsi)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (2 * ((((mustar1/ewv1_h - 1/ewusr) * dmustar1/ewv1_h - sigx1/ewusr) * eusq1 -
      ((mustar2/ewv1_h - 2/ewusr) * dmustar2/ewv1_h - 2 * (sigx2/ewusr)) *
        eusq2) * (1 - prZ)/ewusr) + depsi * prZ * (S^2 * (epsilon)^2/ewv2_h^2 -
      1)/ewv2_h^3 - sigx9^2/sigx8)/sigx8, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu *
      ewv1/(2 * ewu)^2)) * pmustar1 - 0.5 * (dmustar1 * ewv1_h/ewusr))/ewusr +
      (0.5 * (mustar1/ewusr) - sigx4/ewv1_h) * dmustar1) * eusq1 - (((mustar2/ewusr -
      (2 * (ewv1/ewu) + S * (epsilon)/ewusr)/ewv1_h) * dmustar2 + (pmustar2 -
      2 * sigx6)/ewusr) * eusq2 + 0.5 * (sigx1 * eusq1 - sigx2 * eusq2)))/(sigx8 *
      ewusr) - (sigx5 * eusq1 - sigx7) * sigx9 * ewusr/(sigx8 * ewusr)^2) *
      (1 - prZ)), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((dmustar1 * (ewv1/(2 * ewu) - (sigx10 *
      mustar1 + 0.5))/ewv1_h - sigx11/ewusr) * eusq1 - ((2 * (ewv1/ewu) - (mustar2 *
      sigx12 + 0.5)) * dmustar2/ewv1_h - 2 * ((2 * (ewv1 * pmustar2/ewu) -
      dmustar2 * sigx12)/ewusr)) * eusq2)/(sigx8 * ewusr) - sigx9 * ewusr *
      sigx13/(sigx8 * ewusr)^2) * (1 - prZ)), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * prZ * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx8 * ewv2_h^3) - sigx24 * sigx9 * ewv2_h/(sigx8 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * ((sigx1 *
    eusq1 - sigx2 * eusq2)/ewusr) - (sigx9 * sigx14/sigx8 + S * depsi * (epsilon)/ewv2_h^3)) *
    prZ * ewz/sigx8, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewv1_h * mustar1/ewusr) -
      0.5) - 0.5 * sigx4) * dmustar1 * ewv1_h/ewusr - (sigx5 * sigx4 + (2 *
      ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S *
      (epsilon)/ewusr)) * pmustar1)) * eusq1 - ((((mustar2 * ewv1_h - S * (epsilon))/ewusr -
      (0.5 + 2 * (ewv1/ewu))) * dmustar2 * ewv1_h/ewusr + (0.5 * (S * (epsilon)/ewusr) +
      2 * (ewv1/ewu)) * pmustar2 - (2 * (ewv1/ewu) + S * (epsilon)/ewusr) *
      sigx6) * eusq2 + 0.5 * (sigx5 * eusq1 - sigx6 * eusq2)))/(sigx8 * ewusr) -
      (sigx5 * eusq1 - sigx7) * (0.5 * (sigx8 * ewusr) + 2 * ((sigx5 * eusq1 -
        sigx7) * (1 - prZ)))/(sigx8 * ewusr)^2) * (1 - prZ)), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((((dmustar1 *
    ewv1_h/(4 * (ewu * ewusr)) - 2 * (ewu * pmustar1/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx10 * mustar1) - 0.25) * dmustar1 * ewv1_h/ewusr + sigx4 * sigx11)) *
    eusq1 - (2 * ((dmustar2 * ewv1_h/ewusr - pmustar2) * ewv1/ewu) - ((mustar2 *
    sigx12 - 0.5) * dmustar2 * ewv1_h/ewusr + (2 * (ewv1 * pmustar2/ewu) - dmustar2 *
    sigx12) * (2 * (ewv1/ewu) + S * (epsilon)/ewusr))) * eusq2)/(sigx8 * ewusr) -
    (0.5 * (sigx8 * ewusr) + 2 * ((sigx5 * eusq1 - sigx7) * (1 - prZ))) * sigx13/(sigx8 *
      ewusr)^2) * (1 - prZ)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (2 * ((sigx5 * eusq1 - sigx7) * sigx24 * (1 - prZ) * prZ * ewv2_h/((sigx8 *
      ewv2_h)^2 * ewusr))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (sigx5 * eusq1 - sigx7) * (2 - 2 * ((1 - prZ) * sigx14/sigx8)) * prZ * ewz/(sigx8 *
    ewusr), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    2 * (((((sigx11/2 + (pmustar1 - sigx10 * dmustar1)/2) * ewv1/ewu - (0.25 *
    (ewv1_h/ewusr) + 0.25 * (S * (epsilon)/ewv1_h) - sigx10^2 * mustar1) * dmustar1) *
    eusq1 - ((2 * (2 * (ewv1 * pmustar2/ewu) - dmustar2 * sigx12) + 2 * (pmustar2 -
    dmustar2 * sigx12)) * ewv1/ewu - (0.25 * (S * (epsilon)/ewv1_h) + 0.5 * (ewv1_h/ewusr) -
    mustar2 * sigx12^2) * dmustar2) * eusq2)/(sigx8 * ewusr) - 2 * ((1 - prZ) *
    sigx13^2/(sigx8 * ewusr)^2)) * (1 - prZ)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (2 * (sigx24 * (1 - prZ) * prZ * sigx13 * ewv2_h/((sigx8 *
      ewv2_h)^2 * ewusr))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 - 2 * ((1 - prZ) * sigx14/sigx8)) * prZ *
      sigx13 * ewz/(sigx8 * ewusr), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prZ * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx8 * ewv2_h^3) - (sigx24 * prZ +
      0.5 * (sigx8 * ewv2_h)) * sigx24/(sigx8 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx14 * prZ/sigx8 + 1) * sigx24 * prZ * ewz/(sigx8 *
      ewv2_h)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (sigx14 * prZ/sigx8 + 1) * ewz) * sigx14 *
      prZ * ewz/sigx8, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf generalized_genexponential-normal distribution
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
zisfgenexponormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfgenexponormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfgenexponormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfgenexponormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfgenexponormlike_logit,
    grad = cgradzisfgenexponormlike_logit, hess = chesszisfgenexponormlike_logit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfgenexponormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgenexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfgenexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfgenexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfgenexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfgenexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgenexponormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfgenexponormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfgenexponormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfgenexponormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgenexponormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

## cauchit specification class membership
zisfgenexponormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfgenexponormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfgenexponormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfgenexponormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfgenexponormlike_cauchit,
    grad = cgradzisfgenexponormlike_cauchit, hess = chesszisfgenexponormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfgenexponormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgenexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfgenexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfgenexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfgenexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfgenexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgenexponormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfgenexponormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfgenexponormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfgenexponormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgenexponormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

## probit specification class membership
zisfgenexponormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfgenexponormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfgenexponormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfgenexponormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfgenexponormlike_probit,
    grad = cgradzisfgenexponormlike_probit, hess = chesszisfgenexponormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfgenexponormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgenexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfgenexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfgenexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfgenexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfgenexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgenexponormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfgenexponormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfgenexponormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfgenexponormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgenexponormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

## cloglog specification class membership
zisfgenexponormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfgenexponormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfgenexponormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfgenexponormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfgenexponormlike_cloglog,
    grad = cgradzisfgenexponormlike_cloglog, hess = chesszisfgenexponormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfgenexponormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgenexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfgenexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfgenexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfgenexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfgenexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgenexponormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfgenexponormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfgenexponormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfgenexponormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgenexponormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

# Different sigma_v

## logit specification class membership
mnsfgenexponormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfgenexponormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfgenexponormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfgenexponormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfgenexponormlike_logit,
    grad = cgradmnsfgenexponormlike_logit, hess = chessmnsfgenexponormlike_logit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfgenexponormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgenexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfgenexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfgenexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfgenexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfgenexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfgenexponormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfgenexponormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfgenexponormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfgenexponormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfgenexponormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

## cauchit specification class membership
mnsfgenexponormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfgenexponormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfgenexponormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfgenexponormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfgenexponormlike_cauchit,
    grad = cgradmnsfgenexponormlike_cauchit, hess = chessmnsfgenexponormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfgenexponormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgenexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfgenexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfgenexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfgenexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfgenexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfgenexponormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfgenexponormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfgenexponormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfgenexponormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfgenexponormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

## probit specification class membership
mnsfgenexponormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfgenexponormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfgenexponormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfgenexponormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfgenexponormlike_probit,
    grad = cgradmnsfgenexponormlike_probit, hess = chessmnsfgenexponormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfgenexponormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgenexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfgenexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfgenexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfgenexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfgenexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfgenexponormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfgenexponormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfgenexponormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfgenexponormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfgenexponormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

## cloglog specification class membership
mnsfgenexponormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfgenexponormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfgenexponormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfgenexponormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfgenexponormlike_cloglog,
    grad = cgradmnsfgenexponormlike_cloglog, hess = chessmnsfgenexponormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfgenexponormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgenexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfgenexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfgenexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfgenexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfgenexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfgenexponormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfgenexponormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfgenexponormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfgenexponormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfgenexponormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

# Conditional efficiencies estimation ----------
#' efficiencies for zisf genexpo-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v

## logit specification class membership
czisfgenexponormeff_logit <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a - exp(Wv/2)) -
      exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) * pnorm(a +
      exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
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
czisfgenexponormeff_cauchit <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a - exp(Wv/2)) -
      exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) * pnorm(a +
      exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
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
czisfgenexponormeff_probit <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a - exp(Wv/2)) -
      exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) * pnorm(a +
      exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
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
czisfgenexponormeff_cloglog <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a - exp(Wv/2)) -
      exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) * pnorm(a +
      exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
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
cmnsfgenexponormeff_logit <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv1/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a - exp(Wv1/2)) -
      exp(B) * exp(-b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b - exp(Wv1/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a +
      exp(Wv1/2)) - exp(B) * exp(b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b + exp(Wv1/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
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
cmnsfgenexponormeff_cauchit <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv1/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a - exp(Wv1/2)) -
      exp(B) * exp(-b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b - exp(Wv1/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a +
      exp(Wv1/2)) - exp(B) * exp(b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b + exp(Wv1/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
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
cmnsfgenexponormeff_probit <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv1/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a - exp(Wv1/2)) -
      exp(B) * exp(-b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b - exp(Wv1/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a +
      exp(Wv1/2)) - exp(B) * exp(b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b + exp(Wv1/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
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
cmnsfgenexponormeff_cloglog <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv1/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) * (dnorm(b) +
    b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a - exp(Wv1/2)) -
      exp(B) * exp(-b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b - exp(Wv1/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv1/2) + exp(Wv1)/2) * pnorm(a +
      exp(Wv1/2)) - exp(B) * exp(b * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b + exp(Wv1/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
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
#' marginal impact on efficiencies for zisf genexpo-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v

## logit specification class membership
czisfmarggenexponorm_Eu_logit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmarggenexponorm_Vu_logit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
czisfmarggenexponorm_Eu_cauchit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmarggenexponorm_Vu_cauchit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
czisfmarggenexponorm_Eu_probit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmarggenexponorm_Vu_probit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
czisfmarggenexponorm_Eu_cloglog <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmarggenexponorm_Vu_cloglog <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# Different sigma_v

## logit specification class membership
cmnsfmarggenexponorm_Eu_logit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmarggenexponorm_Vu_logit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmnsfmarggenexponorm_Eu_cauchit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmarggenexponorm_Vu_cauchit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmnsfmarggenexponorm_Eu_probit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmarggenexponorm_Vu_probit <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmnsfmarggenexponorm_Eu_cloglog <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmarggenexponorm_Vu_cloglog <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
