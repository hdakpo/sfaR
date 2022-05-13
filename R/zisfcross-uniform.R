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
# Convolution: uniform - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf uniform-normal distribution
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
czisfuninormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv/2)) - pnorm(S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cauchit specification class membership
czisfuninormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv/2)) - pnorm(S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## probit specification class membership
czisfuninormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv/2)) - pnorm(S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cloglog specification class membership
czisfuninormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv/2)) - pnorm(S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# Different sigma_v

## logit specification class membership
cmnsfuninormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv1/2)) - pnorm(S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
cmnsfuninormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv1/2)) - pnorm(S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
cmnsfuninormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv1/2)) - pnorm(S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
cmnsfuninormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv1/2)) - pnorm(S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf uniform-normal distribution
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
cstzisfuninorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, itermax,
  printInfo, tol) {
  cat("Initialization: SFA + uniform - normal distributions...\n")
  initUni <- maxLik(logLik = cuninormlike, start = cstuninorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgraduninormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initUni$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("ZI_",
    colnames(Zvar)))
  names(initUni$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initUni = initUni))
}

# Different sigma_v
cstmnsfuninorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, itermax,
  printInfo, tol) {
  cat("Initialization: SFA + uniform - normal distributions...\n")
  initUni <- maxLik(logLik = cuninormlike, start = cstuninorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgraduninormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initUni$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("ZI_", colnames(Zvar)))
  names(initUni$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initUni = initUni))
}

# Gradient of the likelihood function ----------
#' gradient for zisf uniform-normal distribution
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
cgradzisfuninormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  ezusq <- (sqrt(12) * (wzdeno * ewusr))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  ddue <- (depsi - dmusig)
  wzsq <- ezusq^2
  sigx1 <- (ddue * ewz/ezusq + S * prC * depsi * (epsilon)/ewvsr^2)
  sigx2 <- (prC * depsi/ewvsr + ewz * (pmusig - pepsi)/ezusq)
  sigx3 <- (dmusig/(2 * (wzdeno * ewvsr)) - sqrt(12)/2 * (wzdeno *
    ewusr * (pmusig - pepsi)/wzsq))
  sigx4 <- prC * depsi/(wzdeno * ewvsr)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (sigx5 * ewz/ezusq + sigx6 * prC)
  sigx8 <- (ewusr * ewz/wzsq)
  sigx9 <- (1/ezusq - sqrt(12) * sigx8)
  sigx10 <- (sigx9 * (pmusig - pepsi) - sigx4)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/(sigx2 *
    ewvsr), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx3 *
    ewz/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx7/(sigx2 *
    ewvsr), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx10 *
    ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradzisfuninormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  sigx1 <- (ewz1 * (depsi - dmusig)/(sqrt(12) * ewusr) + S *
    ewz2 * depsi * (epsilon)/ewvsr^2)
  sigx2 <- (ewz2 * depsi/ewvsr + ewz1 * (pmusig - pepsi)/(sqrt(12) *
    ewusr))
  sigx3 <- (ewusr * (pmusig - pepsi)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr) - sqrt(12)/2 * sigx3)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (ewz2 * sigx6 + sigx5 * ewz1/(sqrt(12) * ewusr))
  sigx8 <- ((pmusig - pepsi)/(sqrt(12) * ewusr) - depsi/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/(sigx2 *
    ewvsr), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = ewz1 *
    sigx4/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx7/(sigx2 *
    ewvsr), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx8/(pi *
    sigx2 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradzisfuninormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  sigx1 <- ((depsi - dmusig) * pwZ/(sqrt(12) * ewusr) + S *
    (1 - pwZ) * depsi * (epsilon)/ewvsr^2)
  prZ <- ((1 - pwZ) * depsi/ewvsr + (pmusig - pepsi) * pwZ/(sqrt(12) *
    ewusr))
  sigx2 <- (prZ * ewvsr)
  sigx3 <- (ewusr * (pmusig - pepsi)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr) - sqrt(12)/2 * sigx3)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (sigx5 * pwZ/(sqrt(12) * ewusr) + sigx6 * (1 - pwZ))
  sigx8 <- ((pmusig - pepsi)/(sqrt(12) * ewusr) - depsi/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx4 *
    pwZ/prZ, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx7/sigx2,
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx8 * dwZ/prZ,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradzisfuninormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  sigx1 <- ((1 - prZ) * (depsi - dmusig)/(sqrt(12) * ewusr) +
    S * depsi * prZ * (epsilon)/ewvsr^2)
  sigx2 <- ((1 - prZ) * (pmusig - pepsi)/(sqrt(12) * ewusr) +
    depsi * prZ/ewvsr)
  sigx3 <- (ewusr * (pmusig - pepsi)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr) - sqrt(12)/2 * sigx3)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (sigx5 * (1 - prZ)/(sqrt(12) * ewusr) + sigx6 *
    prZ)
  sigx8 <- ((pmusig - pepsi)/(sqrt(12) * ewusr) - depsi/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/(sigx2 *
    ewvsr), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (1 -
    prZ) * sigx4/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = sigx7/(sigx2 * ewvsr), FUN = "*"), sweep(Zvar,
    MARGIN = 1, STATS = sigx8 * prZ * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v

## logit specification class membership
cgradmnsfuninormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu_h <- exp(Wu/2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wzuv <- (sqrt(12) * (wzdeno * ewu_h * ewv1_h))
  wuepsi <- (sqrt(12) * ewu_h + S * epsilon)
  dwu <- dnorm(wuepsi/ewv1_h, 0, 1)
  pwu <- pnorm(wuepsi/ewv1_h)
  depsi1 <- dnorm(S * epsilon/ewv1_h, 0, 1)
  depsi2 <- dnorm(S * epsilon/ewv2_h, 0, 1)
  pepsi1 <- pnorm(S * epsilon/ewv1_h)
  wz12u <- (sqrt(12) * (wzdeno * ewu_h))
  sigx1 <- ((depsi1 - dwu) * ewz/wzuv + S * prC * depsi2 *
    epsilon/ewv2_h^3)
  sigx2 <- (prC * depsi2/ewv2_h + ewz * (pwu - pepsi1)/wz12u)
  sigx3 <- (wzdeno * ewu_h * (pwu - pepsi1)/wz12u^2)
  sigx4 <- (dwu/(2 * (wzdeno * ewv1_h)) - sqrt(12)/2 * sigx3)
  sigx5 <- (0.5 * (S * depsi1 * epsilon) - 0.5 * (wuepsi *
    dwu))
  sigx6 <- (sqrt(12) * (sigx2 * wzdeno * ewu_h * ewv1_h))
  sigx7 <- (0.5 * (S^2 * depsi2 * epsilon^2/ewv2_h^2) - 0.5 *
    depsi2)
  sigx8 <- (1/wz12u - sqrt(12) * (ewu_h * ewz/wz12u^2))
  sigx9 <- (sigx8 * (pwu - pepsi1) - prC * depsi2/(wzdeno *
    ewv2_h))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx4 *
    ewz/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx5 *
    ewz/sigx6, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx7 *
    prC/(sigx2 * ewv2_h), FUN = "*"), sweep(Zvar, MARGIN = 1,
    STATS = sigx9 * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmnsfuninormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr1, 0, 1)
  depsi1 <- dnorm(S * (epsilon)/ewvsr1, 0, 1)
  depsi2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  pmusig <- pnorm(S * (epsilon)/ewvsr1)
  pepsi <- pnorm(eeusq/ewvsr1)
  uv12 <- (sqrt(12) * (ewusr * ewvsr1))
  sigx1 <- (ewz1 * (depsi1 - dmusig)/uv12 + S * ewz2 * depsi2 *
    (epsilon)/ewvsr2^3)
  sigx2 <- (ewz2 * depsi2/ewvsr2 + ewz1 * (pepsi - pmusig)/(sqrt(12) *
    ewusr))
  sigx3 <- (ewusr * (pepsi - pmusig)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr1) - sqrt(12)/2 * sigx3)
  sigx5 <- (sqrt(12) * (sigx2 * ewusr * ewvsr1))
  sigx6 <- (0.5 * (S * depsi1 * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx7 <- (0.5 * (S^2 * depsi2 * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi2)
  sigx8 <- ((pepsi - pmusig)/(sqrt(12) * ewusr) - depsi2/ewvsr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = ewz1 * sigx4/sigx2,
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx6 *
    ewz1/sigx5, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ewz2 *
    sigx7/(sigx2 * ewvsr2), FUN = "*"), sweep(Zvar, MARGIN = 1,
    STATS = sigx8/(pi * sigx2 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmnsfuninormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr1, 0, 1)
  depsi1 <- dnorm(S * (epsilon)/ewvsr1, 0, 1)
  depsi2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  pmusig <- pnorm(S * (epsilon)/ewvsr1)
  pepsi <- pnorm(eeusq/ewvsr1)
  uv12 <- (sqrt(12) * (ewusr * ewvsr1))
  sigx1 <- ((1 - pwZ) * depsi2/ewvsr2 + (pepsi - pmusig) *
    pwZ/(sqrt(12) * ewusr))
  sigx2 <- ((depsi1 - dmusig) * pwZ/uv12 + S * (1 - pwZ) *
    depsi2 * (epsilon)/ewvsr2^3)
  sigx3 <- (ewusr * (pepsi - pmusig)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr1) - sqrt(12)/2 * sigx3)
  sigx5 <- (sqrt(12) * (sigx1 * ewusr * ewvsr1))
  sigx6 <- (0.5 * (S * depsi1 * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx7 <- (0.5 * (S^2 * depsi2 * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi2)
  sigx8 <- ((pepsi - pmusig)/(sqrt(12) * ewusr) - depsi2/ewvsr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx1,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx4 *
    pwZ/sigx1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx6 *
    pwZ/sigx5, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx7 *
    (1 - pwZ)/(sigx1 * ewvsr2), FUN = "*"), sweep(Zvar, MARGIN = 1,
    STATS = sigx8 * dwZ/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmnsfuninormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr1, 0, 1)
  depsi1 <- dnorm(S * (epsilon)/ewvsr1, 0, 1)
  depsi2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  pmusig <- pnorm(S * (epsilon)/ewvsr1)
  pepsi <- pnorm(eeusq/ewvsr1)
  uv12 <- (sqrt(12) * (ewusr * ewvsr1))
  sigx1 <- ((1 - prZ) * (pepsi - pmusig)/(sqrt(12) * ewusr) +
    depsi2 * prZ/ewvsr2)
  sigx2 <- ((1 - prZ) * (depsi1 - dmusig)/uv12 + S * depsi2 *
    prZ * (epsilon)/ewvsr2^3)
  sigx3 <- (ewusr * (pepsi - pmusig)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr1) - sqrt(12)/2 * sigx3)
  sigx5 <- (sqrt(12) * (sigx1 * ewusr * ewvsr1))
  sigx6 <- (0.5 * (S * depsi1 * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx7 <- (0.5 * (S^2 * depsi2 * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi2)
  sigx8 <- ((pepsi - pmusig)/(sqrt(12) * ewusr) - depsi2/ewvsr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx1,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (1 - prZ) *
    sigx4/sigx1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx6 *
    (1 - prZ)/sigx5, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = sigx7 * prZ/(sigx1 * ewvsr2), FUN = "*"), sweep(Zvar,
    MARGIN = 1, STATS = sigx8 * prZ * ewz/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf uniform-normal distribution
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
chesszisfuninormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  ezusq <- (sqrt(12) * (wzdeno * ewusr))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  ddue <- (depsi - dmusig)
  wzsq <- ezusq^2
  sigx1 <- (ddue * ewz/ezusq + S * prC * depsi * (epsilon)/ewvsr^2)
  sigx2 <- (prC * depsi/ewvsr + ewz * (pmusig - pepsi)/ezusq)
  sigx3 <- (dmusig/(2 * (wzdeno * ewvsr)) - sqrt(12)/2 * (wzdeno *
    ewusr * (pmusig - pepsi)/wzsq))
  sigx4 <- prC * depsi/(wzdeno * ewvsr)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (sigx5 * ewz/ezusq + sigx6 * prC)
  sigx8 <- (ewusr * ewz/wzsq)
  sigx9 <- (1/ezusq - sqrt(12) * sigx8)
  sigx10 <- (sigx9 * (pmusig - pepsi) - sigx4)
  sigx11 <- (S^2 * (epsilon)^2/ewvsr^2 - 1)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((prC * depsi * sigx11 + ewz * (S * depsi *
      (epsilon) - eeusq * dmusig)/ezusq)/(sigx2 * ewvsr^3) -
      sigx1^2/(sigx2 * ewvsr)^2), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (eeusq * dmusig/(2 *
      (wzdeno * ewvsr^2)) - (sigx1 * sigx3/sigx2 + sqrt(12)/2 *
      (wzdeno * ddue * ewusr/wzsq))) * ewz/(sigx2 * ewvsr),
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((0.5 * (depsi * sigx11) - 0.5 * ((eeusq^2/ewvsr^2 -
    1) * dmusig)) * ewz/ezusq + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 -
    2) - 0.5) * prC * depsi * (epsilon)/ewvsr^2)/(sigx2 *
    ewvsr) - sigx7 * sigx1/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx9 * ddue - (sigx10 *
      sigx1/sigx2 + S * prC * depsi * (epsilon)/(wzdeno *
      ewvsr^2))) * ewz/(sigx2 * ewvsr), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((sqrt(12)/2 * (((sqrt(12)/2 * (dmusig/ewvsr) - 12 *
      (wzdeno^2 * ewusr * (pmusig - pepsi)/wzsq)) * ewusr +
      0.5 * (pmusig - pepsi)) * wzdeno/wzsq) + sqrt(12)/2 *
      (eeusq * dmusig)/(2 * (wzdeno * ewvsr^3))) * ewusr +
      sigx3^2 * ewz/sigx2) * ewz/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx7 * sigx3 * ewvsr/(sigx2 *
      ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 * (eeusq^2/ewvsr^2)) *
      dmusig)/(sqrt(12) * wzdeno) + sqrt(12)/2 * (sigx5 *
      wzdeno * ewusr/wzsq))/(sigx2 * ewvsr)) * ewz), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sqrt(12)/2 * (sigx9 * dmusig/ewvsr) -
      (sqrt(12)/2 * wzdeno + sqrt(12) * ((0.5 - 12 * (wzdeno^2 *
        ewusr^2/wzsq)) * ewz)) * (pmusig - pepsi)/wzsq) *
      ewusr - sigx10 * sigx3 * ewz/sigx2) * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((0.25 * (S^3 * depsi *
      (epsilon)^3) - 0.25 * (eeusq^3 * dmusig)) * ewz/ezusq +
      S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
        1) - 0.25) * prC * depsi * (epsilon)^2)/(sigx2 *
      ewvsr^3) - sigx7 * (sigx5 * ewz/ezusq + sigx6 * prC +
      0.5 * (sigx2 * ewvsr))/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((sigx5 * sigx9 - sigx7 * sigx10/sigx2)/ewvsr - (0.5 *
      (S^2 * depsi * (epsilon)^2/(wzdeno * ewvsr^3)) -
      0.5 * (wzdeno * depsi * ewvsr/(wzdeno * ewvsr)^2)) *
      prC) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewvsr) +
      ewvsr/(wzdeno * ewvsr)^2) * depsi - (sigx10^2/sigx2 +
      (sqrt(12) * (1 - 24 * (wzdeno * ewusr^2 * ewz/wzsq)) +
        sqrt(12)) * ewusr * (pmusig - pepsi)/wzsq)) *
      ewz + sigx9 * (pmusig - pepsi) - sigx4) * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesszisfuninormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  sigx1 <- (ewz1 * (depsi - dmusig)/(sqrt(12) * ewusr) + S *
    ewz2 * depsi * (epsilon)/ewvsr^2)
  sigx2 <- (ewz2 * depsi/ewvsr + ewz1 * (pmusig - pepsi)/(sqrt(12) *
    ewusr))
  sigx3 <- (ewusr * (pmusig - pepsi)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr) - sqrt(12)/2 * sigx3)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (ewz2 * sigx6 + sigx5 * ewz1/(sqrt(12) * ewusr))
  sigx8 <- ((pmusig - pepsi)/(sqrt(12) * ewusr) - depsi/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((ewz2 * depsi * (S^2 * (epsilon)^2/ewvsr^2 -
      1) + ewz1 * (S * depsi * (epsilon) - eeusq * dmusig)/(sqrt(12) *
      ewusr))/(sigx2 * ewvsr^3) - sigx1^2/(sigx2 * ewvsr)^2),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (eeusq * dmusig/(2 *
      ewvsr^2) - (sigx1 * sigx4/sigx2 + sqrt(12)/2 * ((depsi -
      dmusig) * ewusr/(sqrt(12) * ewusr)^2))) * ewz1/(sigx2 *
      ewvsr), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((0.5 * (depsi * (S^2 * (epsilon)^2/ewvsr^2 - 1)) -
    0.5 * ((eeusq^2/ewvsr^2 - 1) * dmusig)) * ewz1/(sqrt(12) *
    ewusr) + S * ewz2 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 -
    2) - 0.5) * depsi * (epsilon)/ewvsr^2)/(sigx2 * ewvsr) -
    sigx7 * sigx1/(sigx2 * ewvsr)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((depsi - dmusig)/(sqrt(12) *
      ewusr) - S * depsi * (epsilon)/ewvsr^2)/(pi * sigx2 *
      ((Wz)^2 + 1)) - pi * sigx1 * sigx8 * ((Wz)^2 + 1)/(pi *
      sigx2 * ((Wz)^2 + 1))^2)/ewvsr, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((ewz1 * sigx4^2/sigx2 + (sqrt(12)/2 * (((sqrt(12)/2 *
      (dmusig/ewvsr) - 12 * sigx3) * ewusr + 0.5 * (pmusig -
      pepsi))/(sqrt(12) * ewusr)^2) + sqrt(12)/2 * (eeusq *
      dmusig)/(2 * ewvsr^3)) * ewusr) * ewz1/sigx2), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx7 * sigx4 * ewvsr/(sigx2 *
      ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 * (eeusq^2/ewvsr^2)) *
      dmusig)/sqrt(12) + sqrt(12)/2 * (sigx5 * ewusr/(sqrt(12) *
      ewusr)^2))/(sigx2 * ewvsr)) * ewz1), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (1/(pi * sigx2 * ((Wz)^2 +
      1)) - pi * sigx8 * ((Wz)^2 + 1) * ewz1/(pi * sigx2 *
      ((Wz)^2 + 1))^2) * sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((0.25 * (S^3 * depsi *
      (epsilon)^3) - 0.25 * (eeusq^3 * dmusig)) * ewz1/(sqrt(12) *
      ewusr) + S^2 * ewz2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
      1) - 0.25) * depsi * (epsilon)^2)/(sigx2 * ewvsr^3) -
      sigx7 * (ewz2 * sigx6 + sigx5 * ewz1/(sqrt(12) *
        ewusr) + 0.5 * (sigx2 * ewvsr))/(sigx2 * ewvsr)^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((sigx5/(sqrt(12) * ewusr) + 0.5 * depsi - 0.5 * (S^2 *
      depsi * (epsilon)^2/ewvsr^2))/(pi * sigx2 * ((Wz)^2 +
      1)) - pi * sigx7 * sigx8 * ((Wz)^2 + 1)/(pi * sigx2 *
      ((Wz)^2 + 1))^2)/ewvsr, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (sigx8 * ((pmusig - pepsi)/(sqrt(12) *
      ewusr) + 2 * (pi * Wz * sigx2) - depsi/ewvsr)/(pi *
      sigx2 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesszisfuninormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  sigx1 <- ((depsi - dmusig) * pwZ/(sqrt(12) * ewusr) + S *
    (1 - pwZ) * depsi * (epsilon)/ewvsr^2)
  prZ <- ((1 - pwZ) * depsi/ewvsr + (pmusig - pepsi) * pwZ/(sqrt(12) *
    ewusr))
  sigx2 <- (prZ * ewvsr)
  sigx3 <- (ewusr * (pmusig - pepsi)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr) - sqrt(12)/2 * sigx3)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (sigx5 * pwZ/(sqrt(12) * ewusr) + sigx6 * (1 - pwZ))
  sigx8 <- ((pmusig - pepsi)/(sqrt(12) * ewusr) - depsi/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((1 - pwZ) * depsi * (S^2 * (epsilon)^2/ewvsr^2 -
      1) + pwZ * (S * depsi * (epsilon) - eeusq * dmusig)/(sqrt(12) *
      ewusr))/(prZ * ewvsr^3) - sigx1^2/sigx2^2), FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (eeusq * dmusig/(2 *
      ewvsr^2) - (sigx1 * sigx4/prZ + sqrt(12)/2 * ((depsi -
      dmusig) * ewusr/(sqrt(12) * ewusr)^2))) * pwZ/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((0.5 * (depsi * (S^2 * (epsilon)^2/ewvsr^2 - 1)) -
    0.5 * ((eeusq^2/ewvsr^2 - 1) * dmusig)) * pwZ/(sqrt(12) *
    ewusr) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) -
    0.5) * (1 - pwZ) * depsi * (epsilon)/ewvsr^2)/sigx2 -
    sigx7 * sigx1/sigx2^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((depsi - dmusig)/(sqrt(12) *
      ewusr) - (sigx1 * sigx8/prZ + S * depsi * (epsilon)/ewvsr^2)) *
      dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((sqrt(12)/2 * (((sqrt(12)/2 * (dmusig/ewvsr) - 12 *
      sigx3) * ewusr + 0.5 * (pmusig - pepsi))/(sqrt(12) *
      ewusr)^2) + sqrt(12)/2 * (eeusq * dmusig)/(2 * ewvsr^3)) *
      ewusr + sigx4^2 * pwZ/prZ) * pwZ/prZ), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx7 * sigx4 * ewvsr/sigx2^2 +
      (0.5 * ((sqrt(12)/2 - sqrt(12)/2 * (eeusq^2/ewvsr^2)) *
        dmusig)/sqrt(12) + sqrt(12)/2 * (sigx5 * ewusr/(sqrt(12) *
        ewusr)^2))/sigx2) * pwZ), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx8 * pwZ/prZ) * sigx4 *
      dwZ/prZ, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((0.25 * (S^3 * depsi *
      (epsilon)^3) - 0.25 * (eeusq^3 * dmusig)) * pwZ/(sqrt(12) *
      ewusr) + S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
      1) - 0.25) * (1 - pwZ) * depsi * (epsilon)^2)/(prZ *
      ewvsr^3) - sigx7 * (sigx5 * pwZ/(sqrt(12) * ewusr) +
      sigx6 * (1 - pwZ) + 0.5 * sigx2)/sigx2^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (sigx5/(sqrt(12) * ewusr) + 0.5 * depsi - (sigx7 * sigx8/prZ +
      0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2))) * dwZ/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx8 * dwZ/prZ + Wz) *
      sigx8 * dwZ/prZ), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesszisfuninormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  sigx1 <- ((1 - prZ) * (depsi - dmusig)/(sqrt(12) * ewusr) +
    S * depsi * prZ * (epsilon)/ewvsr^2)
  sigx2 <- ((1 - prZ) * (pmusig - pepsi)/(sqrt(12) * ewusr) +
    depsi * prZ/ewvsr)
  sigx3 <- (ewusr * (pmusig - pepsi)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr) - sqrt(12)/2 * sigx3)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (sigx5 * (1 - prZ)/(sqrt(12) * ewusr) + sigx6 *
    prZ)
  sigx8 <- ((pmusig - pepsi)/(sqrt(12) * ewusr) - depsi/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((1 - prZ) * (S * depsi * (epsilon) -
      eeusq * dmusig)/(sqrt(12) * ewusr) + depsi * prZ *
      (S^2 * (epsilon)^2/ewvsr^2 - 1))/(sigx2 * ewvsr^3) -
      sigx1^2/(sigx2 * ewvsr)^2), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (eeusq * dmusig/(2 *
      ewvsr^2) - (sigx1 * sigx4/sigx2 + sqrt(12)/2 * ((depsi -
      dmusig) * ewusr/(sqrt(12) * ewusr)^2))) * (1 - prZ)/(sigx2 *
      ewvsr), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((0.5 * (depsi * (S^2 * (epsilon)^2/ewvsr^2 - 1)) -
    0.5 * ((eeusq^2/ewvsr^2 - 1) * dmusig)) * (1 - prZ)/(sqrt(12) *
    ewusr) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) -
    0.5) * depsi * prZ * (epsilon)/ewvsr^2)/(sigx2 * ewvsr) -
    sigx7 * sigx1/(sigx2 * ewvsr)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((depsi - dmusig)/(sqrt(12) *
      ewusr) - (sigx1 * sigx8/sigx2 + S * depsi * (epsilon)/ewvsr^2)) *
      prZ * ewz/(sigx2 * ewvsr), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((1 - prZ) * sigx4^2/sigx2 + (sqrt(12)/2 * (((sqrt(12)/2 *
      (dmusig/ewvsr) - 12 * sigx3) * ewusr + 0.5 * (pmusig -
      pepsi))/(sqrt(12) * ewusr)^2) + sqrt(12)/2 * (eeusq *
      dmusig)/(2 * ewvsr^3)) * ewusr) * (1 - prZ)/sigx2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx7 * sigx4 * ewvsr/(sigx2 *
      ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 * (eeusq^2/ewvsr^2)) *
      dmusig)/sqrt(12) + sqrt(12)/2 * (sigx5 * ewusr/(sqrt(12) *
      ewusr)^2))/(sigx2 * ewvsr)) * (1 - prZ)), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx8 * (1 - prZ)/sigx2) *
      sigx4 * prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((0.25 * (S^3 * depsi *
      (epsilon)^3) - 0.25 * (eeusq^3 * dmusig)) * (1 -
      prZ)/(sqrt(12) * ewusr) + S^2 * (0.5 * (0.5 * (S^2 *
      (epsilon)^2/ewvsr^2) - 1) - 0.25) * depsi * prZ *
      (epsilon)^2)/(sigx2 * ewvsr^3) - sigx7 * (sigx5 *
      (1 - prZ)/(sqrt(12) * ewusr) + sigx6 * prZ + 0.5 *
      (sigx2 * ewvsr))/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (sigx5/(sqrt(12) * ewusr) + 0.5 * depsi - (sigx7 * sigx8/sigx2 +
      0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2))) * prZ *
    ewz/(sigx2 * ewvsr), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * sigx8 * (1 - (sigx8 * prZ/sigx2 +
      1) * ewz) * prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Different sigma_v

## logit specification class membership
chessmnsfuninormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu_h <- exp(Wu/2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wzuv <- (sqrt(12) * (wzdeno * ewu_h * ewv1_h))
  wuepsi <- (sqrt(12) * ewu_h + S * epsilon)
  dwu <- dnorm(wuepsi/ewv1_h, 0, 1)
  pwu <- pnorm(wuepsi/ewv1_h)
  depsi1 <- dnorm(S * epsilon/ewv1_h, 0, 1)
  depsi2 <- dnorm(S * epsilon/ewv2_h, 0, 1)
  pepsi1 <- pnorm(S * epsilon/ewv1_h)
  wz12u <- (sqrt(12) * (wzdeno * ewu_h))
  sigx1 <- ((depsi1 - dwu) * ewz/wzuv + S * prC * depsi2 *
    epsilon/ewv2_h^3)
  sigx2 <- (prC * depsi2/ewv2_h + ewz * (pwu - pepsi1)/wz12u)
  sigx3 <- (wzdeno * ewu_h * (pwu - pepsi1)/wz12u^2)
  sigx4 <- (dwu/(2 * (wzdeno * ewv1_h)) - sqrt(12)/2 * sigx3)
  sigx5 <- (0.5 * (S * depsi1 * epsilon) - 0.5 * (wuepsi *
    dwu))
  sigx6 <- (sqrt(12) * (sigx2 * wzdeno * ewu_h * ewv1_h))
  sigx7 <- (0.5 * (S^2 * depsi2 * epsilon^2/ewv2_h^2) - 0.5 *
    depsi2)
  sigx8 <- (1/wz12u - sqrt(12) * (ewu_h * ewz/wz12u^2))
  sigx9 <- (sigx8 * (pwu - pepsi1) - prC * depsi2/(wzdeno *
    ewv2_h))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (prC * depsi2 * (S^2 * epsilon^2/ewv2_h^2 -
      1)/ewv2_h^3 + ewz * (S * depsi1 * epsilon - wuepsi *
      dwu)/(sqrt(12) * (wzdeno * ewu_h * ewv1_h^3)) - sigx1^2/sigx2)/sigx2,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((wuepsi * dwu/(2 * (wzdeno *
      ewv1_h^2)) - sqrt(12)/2 * (wzdeno * (depsi1 - dwu) *
      ewu_h/wz12u^2))/ewv1_h - sigx1 * sigx4/sigx2) * ewz/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((0.5 * (depsi1 * (S^2 * epsilon^2/ewv1_h^2 - 1)) -
    0.5 * ((wuepsi^2/ewv1_h^2 - 1) * dwu))/sigx6 - sqrt(12) *
    (sigx1 * sigx5 * wzdeno * ewu_h * ewv1_h/sigx6^2)) *
    ewz, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * prC * (S * (0.5 * (S^2 * epsilon^2/ewv2_h^2 -
      2) - 0.5) * depsi2 * epsilon/(sigx2 * ewv2_h^3) -
      sigx1 * sigx7 * ewv2_h/(sigx2 * ewv2_h)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx8 * (depsi1 - dwu)/ewv1_h -
      (sigx9 * sigx1/sigx2 + S * prC * depsi2 * epsilon/(wzdeno *
        ewv2_h^3))) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((sqrt(12)/2 * (((sqrt(12)/2 * (dwu/ewv1_h) - 12 * (wzdeno^2 *
      ewu_h * (pwu - pepsi1)/wz12u^2)) * ewu_h + 0.5 *
      (pwu - pepsi1)) * wzdeno/wz12u^2) + sqrt(12)/2 *
      (wuepsi * dwu)/(2 * (wzdeno * ewv1_h^3))) * ewu_h +
      sigx4^2 * ewz/sigx2) * ewz/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      (wuepsi^2/ewv1_h^2)) * dwu)/(sqrt(12) * (sigx2 *
      wzdeno * ewv1_h)) + sqrt(12) * ((sigx4 * ewz + 0.5 *
      sigx2) * sigx5 * wzdeno * ewu_h * ewv1_h/sigx6^2)) *
      ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx7 * prC * sigx4 * ewv2_h *
      ewz/(sigx2 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sqrt(12)/2 * (sigx8 * dwu/ewv1_h) -
      (sqrt(12)/2 * wzdeno + sqrt(12) * ((0.5 - 12 * (wzdeno^2 *
        ewu_h^2/wz12u^2)) * ewz)) * (pwu - pepsi1)/wz12u^2) *
      ewu_h - sigx9 * sigx4 * ewz/sigx2) * ewz/sigx2, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.25 * (S^3 * depsi1 *
      epsilon^3) - 0.25 * (wuepsi^3 * dwu))/(sqrt(12) *
      (sigx2 * wzdeno * ewu_h * ewv1_h^3)) - sqrt(12) *
      ((sigx5 * ewz/sqrt(12) + 0.5 * (sigx2 * wzdeno *
        ewu_h * ewv1_h)) * sigx5/sigx6^2)) * ewz, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx5 * sigx7 * prC * ewv2_h * ewz/(sqrt(12) * ((sigx2 *
      ewv2_h)^2 * wzdeno * ewu_h * ewv1_h))), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx5 * (1/wz12u - (sigx9/(sqrt(12) *
      (sigx2 * wzdeno * ewu_h)) + sqrt(12) * (ewu_h/wz12u^2)) *
      ewz) * ewz/(sigx2 * ewv1_h), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prC * (S^2 * (0.5 * (0.5 * (S^2 * epsilon^2/ewv2_h^2) -
      1) - 0.25) * depsi2 * epsilon^2/(sigx2 * ewv2_h^3) -
      (sigx7 * prC + 0.5 * (sigx2 * ewv2_h)) * sigx7/(sigx2 *
        ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx9 * sigx7/sigx2 +
      0.5 * (S^2 * depsi2 * epsilon^2/(wzdeno * ewv2_h^2)))/ewv2_h -
      0.5 * (wzdeno * depsi2 * ewv2_h/(wzdeno * ewv2_h)^2)) *
      prC * ewz/sigx2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv2_h) +
      ewv2_h/(wzdeno * ewv2_h)^2) * depsi2 - (sigx9^2/sigx2 +
      (sqrt(12) * (1 - 24 * (wzdeno * ewu_h^2 * ewz/wz12u^2)) +
        sqrt(12)) * ewu_h * (pwu - pepsi1)/wz12u^2)) *
      ewz + sigx8 * (pwu - pepsi1) - prC * depsi2/(wzdeno *
      ewv2_h)) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmnsfuninormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr1, 0, 1)
  depsi1 <- dnorm(S * (epsilon)/ewvsr1, 0, 1)
  depsi2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  pmusig <- pnorm(S * (epsilon)/ewvsr1)
  pepsi <- pnorm(eeusq/ewvsr1)
  uv12 <- (sqrt(12) * (ewusr * ewvsr1))
  sigx1 <- (ewz1 * (depsi1 - dmusig)/uv12 + S * ewz2 * depsi2 *
    (epsilon)/ewvsr2^3)
  sigx2 <- (ewz2 * depsi2/ewvsr2 + ewz1 * (pepsi - pmusig)/(sqrt(12) *
    ewusr))
  sigx3 <- (ewusr * (pepsi - pmusig)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr1) - sqrt(12)/2 * sigx3)
  sigx5 <- (sqrt(12) * (sigx2 * ewusr * ewvsr1))
  sigx6 <- (0.5 * (S * depsi1 * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx7 <- (0.5 * (S^2 * depsi2 * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi2)
  sigx8 <- ((pepsi - pmusig)/(sqrt(12) * ewusr) - depsi2/ewvsr2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (ewz2 * depsi2 * (S^2 * (epsilon)^2/ewvsr2^2 -
      1)/ewvsr2^3 + ewz1 * (S * depsi1 * (epsilon) - eeusq *
      dmusig)/(sqrt(12) * (ewusr * ewvsr1^3)) - sigx1^2/sigx2)/sigx2,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((eeusq * dmusig/(2 *
      ewvsr1^2) - sqrt(12)/2 * ((depsi1 - dmusig) * ewusr/(sqrt(12) *
      ewusr)^2))/ewvsr1 - sigx1 * sigx4/sigx2) * ewz1/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((0.5 * (depsi1 * (S^2 * (epsilon)^2/ewvsr1^2 - 1)) -
    0.5 * ((eeusq^2/ewvsr1^2 - 1) * dmusig))/sigx5 - sqrt(12) *
    (sigx1 * sigx6 * ewusr * ewvsr1/sigx5^2)) * ewz1, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ewz2 * (S * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2 -
      2) - 0.5) * depsi2 * (epsilon)/(sigx2 * ewvsr2^3) -
      sigx1 * sigx7 * ewvsr2/(sigx2 * ewvsr2)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((depsi1 - dmusig)/uv12 -
      S * depsi2 * (epsilon)/ewvsr2^3)/(pi * sigx2 * ((Wz)^2 +
      1)) - pi * sigx1 * sigx8 * ((Wz)^2 + 1)/(pi * sigx2 *
      ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((ewz1 * sigx4^2/sigx2 + (sqrt(12)/2 * (((sqrt(12)/2 *
      (dmusig/ewvsr1) - 12 * sigx3) * ewusr + 0.5 * (pepsi -
      pmusig))/(sqrt(12) * ewusr)^2) + sqrt(12)/2 * (eeusq *
      dmusig)/(2 * ewvsr1^3)) * ewusr) * ewz1/sigx2), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      (eeusq^2/ewvsr1^2)) * dmusig)/(sqrt(12) * (sigx2 *
      ewvsr1)) + sqrt(12) * ((ewz1 * sigx4 + 0.5 * sigx2) *
      sigx6 * ewusr * ewvsr1/sigx5^2)) * ewz1), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (ewz2 * sigx7 * ewz1 * sigx4 *
      ewvsr2/(sigx2 * ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (1/(pi * sigx2 * ((Wz)^2 +
      1)) - pi * sigx8 * ((Wz)^2 + 1) * ewz1/(pi * sigx2 *
      ((Wz)^2 + 1))^2) * sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.25 * (S^3 * depsi1 *
      (epsilon)^3) - 0.25 * (eeusq^3 * dmusig))/(sqrt(12) *
      (sigx2 * ewusr * ewvsr1^3)) - sqrt(12) * ((sigx6 *
      ewz1/sqrt(12) + 0.5 * (sigx2 * ewusr * ewvsr1)) *
      sigx6/sigx5^2)) * ewz1, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (ewz2 * sigx6 * sigx7 * ewz1 * ewvsr2/(sqrt(12) * ((sigx2 *
      ewvsr2)^2 * ewusr * ewvsr1))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx6 * (1/(sqrt(12) * (pi *
      sigx2 * ((Wz)^2 + 1))) - pi * sigx8 * ((Wz)^2 + 1) *
      ewz1/(sqrt(12) * (pi * sigx2 * ((Wz)^2 + 1))^2))/(ewusr *
      ewvsr1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ewz2 * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2) -
      1) - 0.25) * depsi2 * (epsilon)^2/(sigx2 * ewvsr2^3) -
      (ewz2 * sigx7 + 0.5 * (sigx2 * ewvsr2)) * sigx7/(sigx2 *
        ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx7 * (1/(pi * sigx2 *
      ((Wz)^2 + 1)) + pi * sigx8 * ((Wz)^2 + 1) * ewz2/(pi *
      sigx2 * ((Wz)^2 + 1))^2)/ewvsr2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (sigx8 * ((pepsi - pmusig)/(sqrt(12) *
      ewusr) + 2 * (pi * Wz * sigx2) - depsi2/ewvsr2)/(pi *
      sigx2 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmnsfuninormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr1, 0, 1)
  depsi1 <- dnorm(S * (epsilon)/ewvsr1, 0, 1)
  depsi2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  pmusig <- pnorm(S * (epsilon)/ewvsr1)
  pepsi <- pnorm(eeusq/ewvsr1)
  uv12 <- (sqrt(12) * (ewusr * ewvsr1))
  sigx1 <- ((1 - pwZ) * depsi2/ewvsr2 + (pepsi - pmusig) *
    pwZ/(sqrt(12) * ewusr))
  sigx2 <- ((depsi1 - dmusig) * pwZ/uv12 + S * (1 - pwZ) *
    depsi2 * (epsilon)/ewvsr2^3)
  sigx3 <- (ewusr * (pepsi - pmusig)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr1) - sqrt(12)/2 * sigx3)
  sigx5 <- (sqrt(12) * (sigx1 * ewusr * ewvsr1))
  sigx6 <- (0.5 * (S * depsi1 * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx7 <- (0.5 * (S^2 * depsi2 * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi2)
  sigx8 <- ((pepsi - pmusig)/(sqrt(12) * ewusr) - depsi2/ewvsr2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((1 - pwZ) * depsi2 * (S^2 * (epsilon)^2/ewvsr2^2 -
      1)/ewvsr2^3 + pwZ * (S * depsi1 * (epsilon) - eeusq *
      dmusig)/(sqrt(12) * (ewusr * ewvsr1^3)) - sigx2^2/sigx1)/sigx1,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((eeusq * dmusig/(2 *
      ewvsr1^2) - sqrt(12)/2 * ((depsi1 - dmusig) * ewusr/(sqrt(12) *
      ewusr)^2))/ewvsr1 - sigx2 * sigx4/sigx1) * pwZ/sigx1,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((0.5 * (depsi1 * (S^2 * (epsilon)^2/ewvsr1^2 - 1)) -
    0.5 * ((eeusq^2/ewvsr1^2 - 1) * dmusig))/(sqrt(12) *
    (sigx1 * ewusr * ewvsr1)) - sqrt(12) * (sigx2 * sigx6 *
    ewusr * ewvsr1/(sqrt(12) * (sigx1 * ewusr * ewvsr1))^2)) *
    pwZ, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (1 - pwZ) * (S * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2 -
      2) - 0.5) * depsi2 * (epsilon)/(sigx1 * ewvsr2^3) -
      sigx2 * sigx7 * ewvsr2/(sigx1 * ewvsr2)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((depsi1 - dmusig)/uv12 -
      (sigx2 * sigx8/sigx1 + S * depsi2 * (epsilon)/ewvsr2^3)) *
      dwZ/sigx1, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((sqrt(12)/2 * (((sqrt(12)/2 * (dmusig/ewvsr1) - 12 *
      sigx3) * ewusr + 0.5 * (pepsi - pmusig))/(sqrt(12) *
      ewusr)^2) + sqrt(12)/2 * (eeusq * dmusig)/(2 * ewvsr1^3)) *
      ewusr + sigx4^2 * pwZ/sigx1) * pwZ/sigx1), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      (eeusq^2/ewvsr1^2)) * dmusig)/(sqrt(12) * (sigx1 *
      ewvsr1)) + sqrt(12) * ((sigx4 * pwZ + 0.5 * sigx1) *
      sigx6 * ewusr * ewvsr1/(sqrt(12) * (sigx1 * ewusr *
      ewvsr1))^2)) * pwZ), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx7 * (1 - pwZ) * sigx4 *
      ewvsr2 * pwZ/(sigx1 * ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx8 * pwZ/sigx1) *
      sigx4 * dwZ/sigx1, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.25 * (S^3 * depsi1 *
      (epsilon)^3) - 0.25 * (eeusq^3 * dmusig))/(sqrt(12) *
      (sigx1 * ewusr * ewvsr1^3)) - sqrt(12) * ((sigx6 *
      pwZ/sqrt(12) + 0.5 * (sigx1 * ewusr * ewvsr1)) *
      sigx6/(sqrt(12) * (sigx1 * ewusr * ewvsr1))^2)) *
      pwZ, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx6 * sigx7 * (1 - pwZ) * ewvsr2 * pwZ/(sqrt(12) *
      ((sigx1 * ewvsr2)^2 * ewusr * ewvsr1))), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sqrt(12)/12 - sigx8 * pwZ/(sqrt(12) *
      sigx1)) * sigx6 * dwZ/(sigx1 * ewusr * ewvsr1), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (1 - pwZ) * (S^2 * (0.5 * (0.5 * (S^2 *
      (epsilon)^2/ewvsr2^2) - 1) - 0.25) * depsi2 * (epsilon)^2/(sigx1 *
      ewvsr2^3) - (sigx7 * (1 - pwZ) + 0.5 * (sigx1 * ewvsr2)) *
      sigx7/(sigx1 * ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx8 * (1 - pwZ)/sigx1 +
      1) * sigx7 * dwZ/(sigx1 * ewvsr2)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx8 * dwZ/sigx1 + Wz) *
      sigx8 * dwZ/sigx1), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmnsfuninormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  dmusig <- dnorm(eeusq/ewvsr1, 0, 1)
  depsi1 <- dnorm(S * (epsilon)/ewvsr1, 0, 1)
  depsi2 <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  pmusig <- pnorm(S * (epsilon)/ewvsr1)
  pepsi <- pnorm(eeusq/ewvsr1)
  uv12 <- (sqrt(12) * (ewusr * ewvsr1))
  sigx1 <- ((1 - prZ) * (pepsi - pmusig)/(sqrt(12) * ewusr) +
    depsi2 * prZ/ewvsr2)
  sigx2 <- ((1 - prZ) * (depsi1 - dmusig)/uv12 + S * depsi2 *
    prZ * (epsilon)/ewvsr2^3)
  sigx3 <- (ewusr * (pepsi - pmusig)/(sqrt(12) * ewusr)^2)
  sigx4 <- (dmusig/(2 * ewvsr1) - sqrt(12)/2 * sigx3)
  sigx5 <- (sqrt(12) * (sigx1 * ewusr * ewvsr1))
  sigx6 <- (0.5 * (S * depsi1 * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx7 <- (0.5 * (S^2 * depsi2 * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi2)
  sigx8 <- ((pepsi - pmusig)/(sqrt(12) * ewusr) - depsi2/ewvsr2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((1 - prZ) * (S * depsi1 * (epsilon) -
      eeusq * dmusig)/(sqrt(12) * (ewusr * ewvsr1^3)) +
      depsi2 * prZ * (S^2 * (epsilon)^2/ewvsr2^2 - 1)/ewvsr2^3 -
      sigx2^2/sigx1)/sigx1, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((eeusq * dmusig/(2 *
      ewvsr1^2) - sqrt(12)/2 * ((depsi1 - dmusig) * ewusr/(sqrt(12) *
      ewusr)^2))/ewvsr1 - sigx2 * sigx4/sigx1) * (1 - prZ)/sigx1,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((0.5 * (depsi1 * (S^2 * (epsilon)^2/ewvsr1^2 - 1)) -
    0.5 * ((eeusq^2/ewvsr1^2 - 1) * dmusig))/sigx5 - sqrt(12) *
    (sigx2 * sigx6 * ewusr * ewvsr1/sigx5^2)) * (1 - prZ),
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * prZ * (S * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2 -
      2) - 0.5) * depsi2 * (epsilon)/(sigx1 * ewvsr2^3) -
      sigx2 * sigx7 * ewvsr2/(sigx1 * ewvsr2)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((depsi1 - dmusig)/uv12 -
      (sigx2 * sigx8/sigx1 + S * depsi2 * (epsilon)/ewvsr2^3)) *
      prZ * ewz/sigx1, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((1 - prZ) * sigx4^2/sigx1 + (sqrt(12)/2 * (((sqrt(12)/2 *
      (dmusig/ewvsr1) - 12 * sigx3) * ewusr + 0.5 * (pepsi -
      pmusig))/(sqrt(12) * ewusr)^2) + sqrt(12)/2 * (eeusq *
      dmusig)/(2 * ewvsr1^3)) * ewusr) * (1 - prZ)/sigx1),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      (eeusq^2/ewvsr1^2)) * dmusig)/(sqrt(12) * (sigx1 *
      ewvsr1)) + sqrt(12) * (((1 - prZ) * sigx4 + 0.5 *
      sigx1) * sigx6 * ewusr * ewvsr1/sigx5^2)) * (1 -
      prZ)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx7 * (1 - prZ) * sigx4 *
      prZ * ewvsr2/(sigx1 * ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx8 * (1 - prZ)/sigx1) *
      sigx4 * prZ * ewz/sigx1, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.25 * (S^3 * depsi1 *
      (epsilon)^3) - 0.25 * (eeusq^3 * dmusig))/(sqrt(12) *
      (sigx1 * ewusr * ewvsr1^3)) - sqrt(12) * ((sigx6 *
      (1 - prZ)/sqrt(12) + 0.5 * (sigx1 * ewusr * ewvsr1)) *
      sigx6/sigx5^2)) * (1 - prZ), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx6 * sigx7 * (1 - prZ) * prZ * ewvsr2/(sqrt(12) *
      ((sigx1 * ewvsr2)^2 * ewusr * ewvsr1))), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sqrt(12)/12 - sigx8 * (1 -
      prZ)/(sqrt(12) * sigx1)) * sigx6 * prZ * ewz/(sigx1 *
      ewusr * ewvsr1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prZ * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2) -
      1) - 0.25) * depsi2 * (epsilon)^2/(sigx1 * ewvsr2^3) -
      (sigx7 * prZ + 0.5 * (sigx1 * ewvsr2)) * sigx7/(sigx1 *
        ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx8 * prZ/sigx1 + 1) *
      sigx7 * prZ * ewz/(sigx1 * ewvsr2)), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * sigx8 * (1 - (sigx8 * prZ/sigx1 +
      1) * ewz) * prZ * ewz/sigx1, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf uniform-normal distribution
#' @param start starting value for optimization
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
zisfuninormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfuninormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfuninormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfuninormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfuninormlike_logit,
      grad = cgradzisfuninormlike_logit, hess = chesszisfuninormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfuninormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfuninormlike_logit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- chesszisfuninormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfuninormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfuninormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfuninormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

## cauchit specification class membership
zisfuninormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfuninormlike_cauchit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfuninormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfuninormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfuninormlike_cauchit,
      grad = cgradzisfuninormlike_cauchit, hess = chesszisfuninormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfuninormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfuninormlike_cauchit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- chesszisfuninormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfuninormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfuninormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfuninormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

## probit specification class membership
zisfuninormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfuninormlike_probit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfuninormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfuninormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfuninormlike_probit,
      grad = cgradzisfuninormlike_probit, hess = chesszisfuninormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfuninormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfuninormlike_probit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- chesszisfuninormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfuninormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfuninormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfuninormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

## cloglog specification class membership
zisfuninormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfuninormlike_cloglog(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfuninormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfuninormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfuninormlike_cloglog,
      grad = cgradzisfuninormlike_cloglog, hess = chesszisfuninormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfuninormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfuninormlike_cloglog(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- chesszisfuninormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfuninormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfuninormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfuninormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

# Different sigma_v

## logit specification class membership
mnsfuninormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfuninormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsfuninormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfuninormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfuninormlike_logit,
      grad = cgradmnsfuninormlike_logit, hess = chessmnsfuninormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfuninormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfuninormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfuninormlike_logit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- chessmnsfuninormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfuninormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfuninormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfuninormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

## cauchit specification class membership
mnsfuninormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfuninormlike_cauchit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsfuninormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfuninormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfuninormlike_cauchit,
      grad = cgradmnsfuninormlike_cauchit, hess = chessmnsfuninormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfuninormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfuninormlike_cauchit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- chessmnsfuninormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfuninormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfuninormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfuninormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

## probit specification class membership
mnsfuninormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfuninormlike_probit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsfuninormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfuninormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfuninormlike_probit,
      grad = cgradmnsfuninormlike_probit, hess = chessmnsfuninormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfuninormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfuninormlike_probit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- chessmnsfuninormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfuninormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfuninormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfuninormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

## cloglog specification class membership
mnsfuninormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfuninormlike_cloglog(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsfuninormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfuninormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfuninormlike_cloglog,
      grad = cgradmnsfuninormlike_cloglog, hess = chessmnsfuninormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfuninormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfuninormlike_cloglog(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- chessmnsfuninormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfuninormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfuninormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfuninormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

# Conditional efficiencies estimation ----------
#' efficiencies for zisf uniform-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v

## logit specification class membership
czisfuninormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u1_c1 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv/2)) - dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv/2)) - pnorm(object$S *
    epsilon/exp(Wv/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c2 <- rep(0, object$Nobs)
  u2_c2 <- rep(0, object$Nobs)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + exp(Wu/2) * sqrt(12))/exp(Wv/2) + exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2) + exp(Wv/2)))/(pnorm((exp(Wu/2) *
      sqrt(12) + object$S * epsilon)/exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv/2) + exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c2 <- rep(1, object$Nobs)
    teBC2_c2 <- rep(1, object$Nobs)
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv/2) -
        exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) -
        exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv/2) - exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c2 <- rep(1, object$Nobs)
    teBC2_reciprocal_c2 <- rep(1, object$Nobs)
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

## cauchit specification class membership
czisfuninormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u1_c1 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv/2)) - dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv/2)) - pnorm(object$S *
    epsilon/exp(Wv/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c2 <- rep(0, object$Nobs)
  u2_c2 <- rep(0, object$Nobs)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + exp(Wu/2) * sqrt(12))/exp(Wv/2) + exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2) + exp(Wv/2)))/(pnorm((exp(Wu/2) *
      sqrt(12) + object$S * epsilon)/exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv/2) + exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c2 <- rep(1, object$Nobs)
    teBC2_c2 <- rep(1, object$Nobs)
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv/2) -
        exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) -
        exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv/2) - exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c2 <- rep(1, object$Nobs)
    teBC2_reciprocal_c2 <- rep(1, object$Nobs)
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

## probit specification class membership
czisfuninormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u1_c1 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv/2)) - dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv/2)) - pnorm(object$S *
    epsilon/exp(Wv/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c2 <- rep(0, object$Nobs)
  u2_c2 <- rep(0, object$Nobs)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + exp(Wu/2) * sqrt(12))/exp(Wv/2) + exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2) + exp(Wv/2)))/(pnorm((exp(Wu/2) *
      sqrt(12) + object$S * epsilon)/exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv/2) + exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c2 <- rep(1, object$Nobs)
    teBC2_c2 <- rep(1, object$Nobs)
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv/2) -
        exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) -
        exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv/2) - exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c2 <- rep(1, object$Nobs)
    teBC2_reciprocal_c2 <- rep(1, object$Nobs)
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

## cloglog specification class membership
czisfuninormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u1_c1 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv/2)) - dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv/2)) - pnorm(object$S *
    epsilon/exp(Wv/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c2 <- rep(0, object$Nobs)
  u2_c2 <- rep(0, object$Nobs)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + exp(Wu/2) * sqrt(12))/exp(Wv/2) + exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2) + exp(Wv/2)))/(pnorm((exp(Wu/2) *
      sqrt(12) + object$S * epsilon)/exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv/2) + exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c2 <- rep(1, object$Nobs)
    teBC2_c2 <- rep(1, object$Nobs)
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv/2) -
        exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) -
        exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv/2) - exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c2 <- rep(1, object$Nobs)
    teBC2_reciprocal_c2 <- rep(1, object$Nobs)
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

# Different sigma_v

## logit specification class membership
cmnsfuninormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u1_c1 <- -exp(Wv1/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv1/2)) - dnorm(object$S * epsilon/exp(Wv1/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
    epsilon/exp(Wv1/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv1/2) * (dnorm(object$S * epsilon/exp(Wv1/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv1/2))) - object$S * epsilon/exp(Wv1/2))
  u1_c2 <- rep(0, object$Nobs)
  u2_c2 <- rep(0, object$Nobs)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (pnorm((object$S *
      epsilon + exp(Wu/2) * sqrt(12))/exp(Wv1/2) + exp(Wv1/2)) -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(pnorm((exp(Wu/2) *
      sqrt(12) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
      epsilon/exp(Wv1/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_c2 <- rep(1, object$Nobs)
    teBC2_c2 <- rep(1, object$Nobs)
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv1/2) -
        exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2) -
        exp(Wv1/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_reciprocal_c2 <- rep(1, object$Nobs)
    teBC2_reciprocal_c2 <- rep(1, object$Nobs)
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

## cauchit specification class membership
cmnsfuninormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u1_c1 <- -exp(Wv1/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv1/2)) - dnorm(object$S * epsilon/exp(Wv1/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
    epsilon/exp(Wv1/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv1/2) * (dnorm(object$S * epsilon/exp(Wv1/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv1/2))) - object$S * epsilon/exp(Wv1/2))
  u1_c2 <- rep(0, object$Nobs)
  u2_c2 <- rep(0, object$Nobs)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (pnorm((object$S *
      epsilon + exp(Wu/2) * sqrt(12))/exp(Wv1/2) + exp(Wv1/2)) -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(pnorm((exp(Wu/2) *
      sqrt(12) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
      epsilon/exp(Wv1/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_c2 <- rep(1, object$Nobs)
    teBC2_c2 <- rep(1, object$Nobs)
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv1/2) -
        exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2) -
        exp(Wv1/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_reciprocal_c2 <- rep(1, object$Nobs)
    teBC2_reciprocal_c2 <- rep(1, object$Nobs)
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

## probit specification class membership
cmnsfuninormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u1_c1 <- -exp(Wv1/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv1/2)) - dnorm(object$S * epsilon/exp(Wv1/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
    epsilon/exp(Wv1/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv1/2) * (dnorm(object$S * epsilon/exp(Wv1/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv1/2))) - object$S * epsilon/exp(Wv1/2))
  u1_c2 <- rep(0, object$Nobs)
  u2_c2 <- rep(0, object$Nobs)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (pnorm((object$S *
      epsilon + exp(Wu/2) * sqrt(12))/exp(Wv1/2) + exp(Wv1/2)) -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(pnorm((exp(Wu/2) *
      sqrt(12) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
      epsilon/exp(Wv1/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_c2 <- rep(1, object$Nobs)
    teBC2_c2 <- rep(1, object$Nobs)
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv1/2) -
        exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2) -
        exp(Wv1/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_reciprocal_c2 <- rep(1, object$Nobs)
    teBC2_reciprocal_c2 <- rep(1, object$Nobs)
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

## cloglog specification class membership
cmnsfuninormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u1_c1 <- -exp(Wv1/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv1/2)) - dnorm(object$S * epsilon/exp(Wv1/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
    epsilon/exp(Wv1/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv1/2) * (dnorm(object$S * epsilon/exp(Wv1/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv1/2))) - object$S * epsilon/exp(Wv1/2))
  u1_c2 <- rep(0, object$Nobs)
  u2_c2 <- rep(0, object$Nobs)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (pnorm((object$S *
      epsilon + exp(Wu/2) * sqrt(12))/exp(Wv1/2) + exp(Wv1/2)) -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(pnorm((exp(Wu/2) *
      sqrt(12) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
      epsilon/exp(Wv1/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_c2 <- rep(1, object$Nobs)
    teBC2_c2 <- rep(1, object$Nobs)
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv1/2) -
        exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2) -
        exp(Wv1/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_reciprocal_c2 <- rep(1, object$Nobs)
    teBC2_reciprocal_c2 <- rep(1, object$Nobs)
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for zisf uniform-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v

## logit specification class membership
czisfmarguninorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarguninorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
czisfmarguninorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarguninorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## probit specification class membership
czisfmarguninorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarguninorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
czisfmarguninorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarguninorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# Different sigma_v

## logit specification class membership
cmnsfmarguninorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarguninorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmnsfmarguninorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarguninorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmnsfmarguninorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarguninorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmnsfmarguninorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarguninorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}
