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
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf halfnormal-normal distribution
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
czisfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cauchit specification class membership
czisfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## probit specification class membership
czisfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cloglog specification class membership
czisfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# Different sigma_v

## logit specification class membership
cmnsfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
cmnsfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
cmnsfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
cmnsfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf halfnormal-normal distribution
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
cstzisfhalfnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA halfnormal - normal distribution...\n")
  initHalf <- maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradhalfnormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initHalf$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("ZI_",
    colnames(Zvar)))
  names(initHalf$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Different sigma_v
cstmnsfhalfnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA halfnormal - normal distribution...\n")
  initHalf <- maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradhalfnormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initHalf$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("ZI_", colnames(Zvar)))
  names(initHalf$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for zisf halfnormal-normal distribution
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
cgradzisfhalfnormlike_logit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  dmusig <- dnorm(S * (epsilon)/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  depsi <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- (dmusig * ewz * pmusig/(wzdeno * sqrt(sigma_sq)))
  sigx2 <- (prC * dwsr/ewvsr + 2 * sigx1)
  sigx3 <- (S * dmusig * pmusig * (epsilon)/(sigma_sq)^2)
  sigx4 <- (wzdeno * dmusig * pmusig/(wzdeno * sqrt(sigma_sq))^2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (2 *
    ((depsi * dmusig * ewu/sigmastar + S * dmusig * pmusig *
      (epsilon)) * ewz/(wzdeno * (sigma_sq)^(3/2))) + S *
    prC * dwsr * (epsilon)/ewvsr^3)/sigx2, FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = 2 * (ewu * ewz * (S * (0.5 * sigx3 -
      (1/(ssq) - (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
        sigmastar) * ewu/(ssq)^2) * depsi * dmusig) *
      (epsilon)/wzdeno - 0.5 * sigx4)/(sigx2 * sqrt(sigma_sq))),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ((0.5 *
    (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr) * prC/ewvsr +
    2 * (ewv * ewz * (S * ((0.5 * ((1 - ewv/(sigma_sq)) *
      ewu/sigmastar) + sigmastar) * depsi * dmusig * ewu/(ssq)^2 +
      0.5 * sigx3) * (epsilon)/wzdeno - 0.5 * sigx4)/sqrt(sigma_sq)))/sigx2,
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = (2 * ((1/(wzdeno *
    sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/(wzdeno * sqrt(sigma_sq))^2) *
    dmusig * pmusig) - prC * dwsr/(wzdeno * ewvsr)) * ewz/sigx2,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradzisfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewusr <- exp(Wu/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq <- ewu + ewv
  mu <- 0
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm((mu * ewv - S * ewu * (epsilon))/(ssq))
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  dmu <- dnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  prV <- (1 - ewv/(sigma_sq))
  s2sig <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  dpmu <- (dmusig * muepsi^2 * pmusig/(sigma_sq)^2)
  pssq <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  ewusq <- (ewu * pwu/sqrt(sigma_sq))
  sigx1 <- (mu/(ssq) - s2sig * (mu * ewv - S * ewu * (epsilon))/(ssq)^2)
  sigx2 <- (dmu * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmusig)
  sigx12 <- (pwu * sqrt(sigma_sq))
  sigx3 <- (ewz2 * depsi/ewvsr + ewz1 * dmusig * pmusig/sigx12)
  sigx4 <- (ewz1 * sigx2/pssq + S * ewz2 * depsi * (epsilon)/ewvsr^3)
  sigx5 <- (dmu * dmusig * ewv/sigmastar - dmusig * muepsi *
    pmusig)
  sigx6 <- (ewusr * sigx12^2)
  sigx7 <- (sigx5/pssq - dmusig * dwu * pmusig * sqrt(sigma_sq)/sigx6)
  sigx8 <- (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
    sigmastar)
  sigx9 <- (sigx8 * (mu * ewv - S * ewu * (epsilon))/(ssq)^2 +
    S * (epsilon)/(ssq))
  sigx10 <- (0.5 * ewusq - 0.5 * (mu * dwu * sqrt(sigma_sq)/ewusr))
  sigx11 <- (0.5 * dpmu - sigx9 * dmu * dmusig)
  sigx13 <- (sigx11 * ewu/sigx12 - sigx10 * dmusig * pmusig/sigx12^2)
  sigx14 <- ((0.5 * dpmu + dmu * dmusig * sigx1)/pwu - 0.5 *
    (dmusig * pmusig * pwu/sigx12^2))
  sigx15 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx16 <- (sigx14 * ewz1 * ewv/sqrt(sigma_sq) + ewz2 * sigx15/ewvsr)
  sigx17 <- (dmusig * pmusig/sigx12 - depsi/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx3,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx13 *
    ewz1/sigx3, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16/sigx3,
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx17/(pi *
    sigx3 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradzisfhalfnormlike_probit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  mu <- 0
  ewvsr <- exp(Wv/2)
  ewusr <- exp(Wu/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm((mu * ewv - S * ewu * (epsilon))/(ssq))
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  dmu <- dnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  prV <- (1 - ewv/(sigma_sq))
  s2sig <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  dpmu <- (dmusig * muepsi^2 * pmusig/(sigma_sq)^2)
  pssq <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  ewusq <- (ewu * pwu/sqrt(sigma_sq))
  sigx1 <- (mu/(ssq) - s2sig * (mu * ewv - S * ewu * (epsilon))/(ssq)^2)
  sigx2 <- (dmu * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmusig)
  sigx12 <- (pwu * sqrt(sigma_sq))
  sigx3 <- ((1 - pwZ) * depsi/ewvsr + dmusig * pmusig * pwZ/sigx12)
  sigx4 <- (sigx2 * pwZ/pssq + S * (1 - pwZ) * depsi * (epsilon)/ewvsr^3)
  sigx5 <- (dmu * dmusig * ewv/sigmastar - dmusig * muepsi *
    pmusig)
  sigx6 <- (ewusr * sigx12^2)
  sigx7 <- (sigx5/pssq - dmusig * dwu * pmusig * sqrt(sigma_sq)/sigx6)
  sigx8 <- (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
    sigmastar)
  sigx9 <- (sigx8 * (mu * ewv - S * ewu * (epsilon))/(ssq)^2 +
    S * (epsilon)/(ssq))
  sigx10 <- (0.5 * ewusq - 0.5 * (mu * dwu * sqrt(sigma_sq)/ewusr))
  sigx11 <- (0.5 * dpmu - sigx9 * dmu * dmusig)
  sigx13 <- (sigx11 * ewu/sigx12 - sigx10 * dmusig * pmusig/sigx12^2)
  sigx14 <- ((0.5 * dpmu + dmu * dmusig * sigx1)/pwu - 0.5 *
    (dmusig * pmusig * pwu/sigx12^2))
  sigx15 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx16 <- (sigx14 * ewv * pwZ/sqrt(sigma_sq) + sigx15 * (1 -
    pwZ)/ewvsr)
  sigx17 <- (dmusig * pmusig/sigx12 - depsi/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx3,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx13 *
    pwZ/sigx3, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16/sigx3,
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx17 *
    dwZ/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradzisfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  mu <- 0
  ewvsr <- exp(Wv/2)
  ewusr <- exp(Wu/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm((mu * ewv - S * ewu * (epsilon))/(ssq))
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  dmu <- dnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  prV <- (1 - ewv/(sigma_sq))
  s2sig <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  dpmu <- (dmusig * muepsi^2 * pmusig/(sigma_sq)^2)
  pssq <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  ewusq <- (ewu * pwu/sqrt(sigma_sq))
  sigx1 <- (mu/(ssq) - s2sig * (mu * ewv - S * ewu * (epsilon))/(ssq)^2)
  sigx2 <- (dmu * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmusig)
  sigx12 <- (pwu * sqrt(sigma_sq))
  sigx3 <- ((1 - prZ) * dmusig * pmusig/sigx12 + depsi * prZ/ewvsr)
  sigx4 <- ((1 - prZ) * sigx2/pssq + S * depsi * prZ * (epsilon)/ewvsr^3)
  sigx5 <- (dmu * dmusig * ewv/sigmastar - dmusig * muepsi *
    pmusig)
  sigx6 <- (ewusr * sigx12^2)
  sigx7 <- (sigx5/pssq - dmusig * dwu * pmusig * sqrt(sigma_sq)/sigx6)
  sigx8 <- (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
    sigmastar)
  sigx9 <- (sigx8 * (mu * ewv - S * ewu * (epsilon))/(ssq)^2 +
    S * (epsilon)/(ssq))
  sigx10 <- (0.5 * ewusq - 0.5 * (mu * dwu * sqrt(sigma_sq)/ewusr))
  sigx11 <- (0.5 * dpmu - sigx9 * dmu * dmusig)
  sigx13 <- (sigx11 * ewu/sigx12 - sigx10 * dmusig * pmusig/sigx12^2)
  sigx14 <- ((0.5 * dpmu + dmu * dmusig * sigx1)/pwu - 0.5 *
    (dmusig * pmusig * pwu/sigx12^2))
  sigx15 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx16 <- (sigx14 * (1 - prZ) * ewv/sqrt(sigma_sq) + sigx15 *
    prZ/ewvsr)
  sigx17 <- (dmusig * pmusig/sigx12 - depsi/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx3,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx13 *
    (1 - prZ)/sigx3, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = sigx16/sigx3, FUN = "*"), sweep(Zvar, MARGIN = 1,
    STATS = sigx17 * prZ * ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v

## logit specification class membership
cgradmnsfhalfnormlike_logit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv1 <- exp(Wv1)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  ssq <- ((sigma_sq) * sigmastar)
  dmusig <- dnorm(-(S * ewu * (epsilon)/ssq))
  pmusig <- pnorm(-(S * ewu * (epsilon)/ssq))
  depsi <- dnorm(S * (epsilon)/sqrt(sigma_sq))
  dwsr <- dnorm(S * (epsilon)/ewv2_h)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- (wzdeno * depsi * pmusig/(wzdeno * sqrt(sigma_sq))^2)
  sigx2 <- (S * depsi * pmusig * (epsilon)/(sigma_sq)^2)
  wsqsq <- (wzdeno * (sigma_sq) * sqrt(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx3 <- (0.5 * (prU * ewv1/sigmastar) + sigmastar)
  sigx4 <- (0.5 * sigx2 - (1/ssq - sigx3 * ewu/ssq^2) * dmusig *
    depsi)
  sigx5 <- (S * sigx4 * (epsilon)/wzdeno - 0.5 * sigx1)
  sigx6 <- (depsi * ewz * pmusig/(wzdeno * sqrt(sigma_sq)))
  sigx7 <- (prC * dwsr/ewv2_h + 2 * sigx6)
  sigx8 <- (sigx7 * sqrt(sigma_sq))
  prV <- (1 - ewv1/(sigma_sq))
  sigx9 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx10 <- (sigx9 * dmusig * depsi * ewu/ssq^2 + 0.5 * sigx2)
  sigx11 <- (S * sigx10 * (epsilon)/wzdeno - 0.5 * sigx1)
  sigx12 <- (S^2 * dwsr * (epsilon)^2/ewv2_h^2)
  sigx13 <- (0.5 * sigx12 - 0.5 * dwsr)
  sigx14 <- (1/(wzdeno * sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/(wzdeno *
    sqrt(sigma_sq))^2)
  sigx15 <- (2 * (sigx14 * depsi * pmusig) - prC * dwsr/(wzdeno *
    ewv2_h))
  sigx16 <- (dmusig * depsi * ewu/sigmastar + S * depsi * pmusig *
    (epsilon))
  sigx17 <- (sigx16 * ewz/wsqsq)
  sigx18 <- (2 * sigx17 + S * prC * dwsr * (epsilon)/ewv2_h^3)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx18/sigx7,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (ewu *
    ewz * sigx5/sigx8), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = 2 * (ewv1 * ewz * sigx11/sigx8), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx13 * prC/(sigx7 *
      ewv2_h), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx15 *
      ewz/sigx7, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmnsfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  mu <- 0
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  ssq <- ((sigma_sq) * sigmastar)
  muvu <- (mu * ewv1 - S * ewu * (epsilon))
  mustar <- muvu/ssq
  dmustar <- dnorm(mustar, 0, 1)
  pmustar <- pnorm(mustar)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  depsi <- dnorm(S * (epsilon)/ewvsr2)
  musq <- (mu * dwu * sqrt(sigma_sq)/ewusr)
  prV <- (1 - ewv1/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx1 <- (dmustar * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmustar)
  sigx2 <- (ewz2 * depsi/ewvsr2 + ewz1 * dmusig * pmustar/(pwu *
    sqrt(sigma_sq)))
  sigx3 <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  sigx4 <- (ewz1 * sigx1/sigx3 + S * ewz2 * depsi * (epsilon)/ewvsr2^3)
  sigx5 <- (dmustar * dmusig * ewv1/sigmastar - dmusig * muepsi *
    pmustar)
  sigx6 <- (ewusr * (pwu * sqrt(sigma_sq))^2)
  sigx7 <- (sigx5/sigx3 - dmusig * dwu * pmustar * sqrt(sigma_sq)/sigx6)
  sigx8 <- (dmusig * muepsi^2 * pmustar/(sigma_sq)^2)
  sigx9 <- (0.5 * (prU * ewv1/sigmastar) + sigmastar)
  sigx10 <- (sigx9 * muvu/ssq^2 + S * (epsilon)/ssq)
  sigx11 <- (0.5 * sigx8 - sigx10 * dmustar * dmusig)
  sigx12 <- (0.5 * (ewu * pwu/sqrt(sigma_sq)) - 0.5 * musq)
  sigx13 <- (sigx11 * ewu/(pwu * sqrt(sigma_sq)) - sigx12 *
    dmusig * pmustar/(pwu * sqrt(sigma_sq))^2)
  sigx14 <- (mu/ssq - (0.5 * (prV * ewu/sigmastar) + sigmastar) *
    muvu/ssq^2)
  sigx15 <- (0.5 * sigx8 + dmustar * dmusig * sigx14)
  sigx16 <- (dmusig * pmustar * pwu/(pwu * sqrt(sigma_sq))^2)
  sigx17 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi)
  sigx18 <- (dmusig * pmustar/(pwu * sqrt(sigma_sq)) - depsi/ewvsr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx2,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx13 *
    ewz1/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx15/pwu -
    0.5 * sigx16) * ewz1 * ewv1/(sigx2 * sqrt(sigma_sq)),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ewz2 * sigx17/(sigx2 *
    ewvsr2), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx18/(pi *
    sigx2 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmnsfhalfnormlike_probit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  mu <- 0
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  ssq <- ((sigma_sq) * sigmastar)
  muvu <- (mu * ewv1 - S * ewu * (epsilon))
  mustar <- muvu/ssq
  dmustar <- dnorm(mustar, 0, 1)
  pmustar <- pnorm(mustar)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  depsi <- dnorm(S * (epsilon)/ewvsr2)
  musq <- (mu * dwu * sqrt(sigma_sq)/ewusr)
  prV <- (1 - ewv1/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx1 <- (dmustar * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmustar)
  sigx2 <- ((1 - pwZ) * depsi/ewvsr2 + dmusig * pmustar * pwZ/(pwu *
    sqrt(sigma_sq)))
  sigx3 <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  sigx4 <- (sigx1 * pwZ/sigx3 + S * (1 - pwZ) * depsi * (epsilon)/ewvsr2^3)
  sigx5 <- (dmustar * dmusig * ewv1/sigmastar - dmusig * muepsi *
    pmustar)
  sigx6 <- (ewusr * (pwu * sqrt(sigma_sq))^2)
  sigx7 <- (sigx5/sigx3 - dmusig * dwu * pmustar * sqrt(sigma_sq)/sigx6)
  sigx8 <- (dmusig * muepsi^2 * pmustar/(sigma_sq)^2)
  sigx9 <- (0.5 * (prU * ewv1/sigmastar) + sigmastar)
  sigx10 <- (sigx9 * muvu/ssq^2 + S * (epsilon)/ssq)
  sigx11 <- (0.5 * sigx8 - sigx10 * dmustar * dmusig)
  sigx12 <- (0.5 * (ewu * pwu/sqrt(sigma_sq)) - 0.5 * musq)
  sigx13 <- (sigx11 * ewu/(pwu * sqrt(sigma_sq)) - sigx12 *
    dmusig * pmustar/(pwu * sqrt(sigma_sq))^2)
  sigx14 <- (mu/ssq - (0.5 * (prV * ewu/sigmastar) + sigmastar) *
    muvu/ssq^2)
  sigx15 <- (0.5 * sigx8 + dmustar * dmusig * sigx14)
  sigx16 <- (dmusig * pmustar * pwu/(pwu * sqrt(sigma_sq))^2)
  sigx17 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi)
  sigx18 <- (dmusig * pmustar/(pwu * sqrt(sigma_sq)) - depsi/ewvsr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx2,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx13 *
    pwZ/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx15/pwu -
    0.5 * sigx16) * ewv1 * pwZ/(sigx2 * sqrt(sigma_sq)),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx17 *
    (1 - pwZ)/(sigx2 * ewvsr2), FUN = "*"), sweep(Zvar, MARGIN = 1,
    STATS = sigx18 * dwZ/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmnsfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq <- ewu + ewv1
  mu <- 0
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  ssq <- ((sigma_sq) * sigmastar)
  muvu <- (mu * ewv1 - S * ewu * (epsilon))
  mustar <- muvu/ssq
  dmustar <- dnorm(mustar, 0, 1)
  pmustar <- pnorm(mustar)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  depsi <- dnorm(S * (epsilon)/ewvsr2)
  musq <- (mu * dwu * sqrt(sigma_sq)/ewusr)
  prV <- (1 - ewv1/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx1 <- (dmustar * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmustar)
  sigx2 <- ((1 - prZ) * dmusig * pmustar/(pwu * sqrt(sigma_sq)) +
    depsi * prZ/ewvsr2)
  sigx3 <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  sigx4 <- ((1 - prZ) * sigx1/sigx3 + S * depsi * prZ * (epsilon)/ewvsr2^3)
  sigx5 <- (dmustar * dmusig * ewv1/sigmastar - dmusig * muepsi *
    pmustar)
  sigx6 <- (ewusr * (pwu * sqrt(sigma_sq))^2)
  sigx7 <- (sigx5/sigx3 - dmusig * dwu * pmustar * sqrt(sigma_sq)/sigx6)
  sigx8 <- (dmusig * muepsi^2 * pmustar/(sigma_sq)^2)
  sigx9 <- (0.5 * (prU * ewv1/sigmastar) + sigmastar)
  sigx10 <- (sigx9 * muvu/ssq^2 + S * (epsilon)/ssq)
  sigx11 <- (0.5 * sigx8 - sigx10 * dmustar * dmusig)
  sigx12 <- (0.5 * (ewu * pwu/sqrt(sigma_sq)) - 0.5 * musq)
  sigx13 <- (sigx11 * ewu/(pwu * sqrt(sigma_sq)) - sigx12 *
    dmusig * pmustar/(pwu * sqrt(sigma_sq))^2)
  sigx14 <- (mu/ssq - (0.5 * (prV * ewu/sigmastar) + sigmastar) *
    muvu/ssq^2)
  sigx15 <- (0.5 * sigx8 + dmustar * dmusig * sigx14)
  sigx16 <- (dmusig * pmustar * pwu/(pwu * sqrt(sigma_sq))^2)
  sigx17 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi)
  sigx18 <- (dmusig * pmustar/(pwu * sqrt(sigma_sq)) - depsi/ewvsr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx2,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx13 *
    (1 - prZ)/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = (sigx15/pwu - 0.5 * sigx16) * (1 - prZ) * ewv1/(sigx2 *
      sqrt(sigma_sq)), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = sigx17 * prZ/(sigx2 * ewvsr2), FUN = "*"), sweep(Zvar,
    MARGIN = 1, STATS = sigx18 * prZ * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf halfnormal-normal distribution
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
chesszisfhalfnormlike_logit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  dmusig <- dnorm(S * (epsilon)/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  depsi <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- (dmusig * ewz * pmusig/(wzdeno * sqrt(sigma_sq)))
  sigx2 <- (prC * dwsr/ewvsr + 2 * sigx1)
  sigx3 <- (S * dmusig * pmusig * (epsilon)/(sigma_sq)^2)
  sigx4 <- (wzdeno * dmusig * pmusig/(wzdeno * sqrt(sigma_sq))^2)
  sigx5 <- (depsi * dmusig * ewu/sigmastar + S * dmusig * pmusig *
    (epsilon))
  sigx6 <- (sigx5 * ewz/(wzdeno * (sigma_sq)^(3/2)))
  sigx7 <- (2 * sigx6 + S * prC * dwsr * (epsilon)/ewvsr^3)
  sigx8 <- (0.5 * ((1 - ewv/(sigma_sq)) * ewu/sigmastar) +
    sigmastar)
  sigx9 <- (wzdeno * sigx5/((wzdeno * sqrt(sigma_sq))^2 * (sigma_sq)))
  sigx10 <- sigx8 * depsi * dmusig * ewu/(ssq)^2
  sigx11 <- (depsi * ewu/sigmastar + S * pmusig * (epsilon))
  sigx12 <- (S * sigx11 * (epsilon)/(sigma_sq) - pmusig)
  sigx13 <- (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
    sigmastar)
  sigx14 <- (1/(ssq) - sigx13 * ewu/(ssq)^2)
  sigx15 <- (dmusig * ewu/ewv + dmusig)
  ddmu <- (0.5 * sigx3 - sigx14 * depsi * dmusig)
  sigx16 <- (S * ddmu * (epsilon)/wzdeno - 0.5 * sigx4)
  sigx17 <- (S * pmusig * (epsilon)/(sigma_sq)^2)
  sigx18 <- (S * dmusig * (S * (0.5 * sigx17 - sigx14 * depsi) *
    (epsilon) - 2 * (pmusig/(sigma_sq))) * (epsilon)/(sigma_sq)^2)
  sigx19 <- ((1/(wzdeno * sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/(wzdeno *
    sqrt(sigma_sq))^2) * dmusig * pmusig)
  sigx20 <- (2 * sigx19 - prC * dwsr/(wzdeno * ewvsr))
  sigx21 <- (ewv * ewz * (S * (sigx10 + 0.5 * sigx3) * (epsilon)/wzdeno -
    0.5 * sigx4)/sqrt(sigma_sq))
  sigx22 <- ((0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 *
    dwsr) * prC/ewvsr + 2 * sigx21)
  sigx23 <- (0.5/sqrt(sigma_sq) - wzdeno^2 * sqrt(sigma_sq)/(wzdeno *
    sqrt(sigma_sq))^2)
  sigx24 <- (sigx23 * ewz + 0.5 * (wzdeno/sqrt(sigma_sq)))
  sigx25 <- (S^2 * (1/(wzdeno * sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/(wzdeno *
    sqrt(sigma_sq))^2) * dmusig * (epsilon)^2/(sigma_sq)^2)
  sigx26 <- (0.5 * sigx25 - sigx24 * dmusig/(wzdeno * sqrt(sigma_sq))^2)
  wzsq <- (wzdeno * sqrt(sigma_sq))^2
  sigx27 <- (wzdeno * (S * ddmu * (epsilon) - wzdeno^2 * dmusig *
    pmusig/wzsq)/wzsq)
  sigx28 <- (1/(wzdeno * sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/wzsq)
  ewvsq <- (1 - ewv/(sigma_sq))
  sigx29 <- (S * (sigx10 + 0.5 * sigx3) * (epsilon)/wzdeno -
    0.5 * sigx4)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (prC * dwsr * (S^2 * (epsilon)^2/ewvsr^2 -
      1)/ewvsr^3 + 2 * ((dmusig * sigx12 + S * depsi *
      sigx15 * ewu * (epsilon)/(ssq)) * ewz/(wzdeno * (sigma_sq)^(3/2))) -
      sigx7^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((sigx14 * depsi *
      dmusig + S * ((0.5 * sigx12 - 0.5 * pmusig) * dmusig/(sigma_sq) -
      S * sigx14 * depsi * sigx15 * (epsilon)) * (epsilon)/(sigma_sq))/wzdeno -
      0.5 * sigx9)/(sigx2 * sqrt(sigma_sq)) - sigx7 * sigx16 *
      sqrt(sigma_sq)/(sigx2 * sqrt(sigma_sq))^2) * ewu *
      ewz), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (2 * (((S * ((0.5 * sigx12 - 0.5 * pmusig) * dmusig/(sigma_sq) +
    S * sigx8 * depsi * sigx15 * ewu * (epsilon)/(ssq)^2) *
    (epsilon)/(sigma_sq) - sigx10)/wzdeno - 0.5 * sigx9) *
    ewv * ewz/sqrt(sigma_sq)) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 -
    2) - 0.5) * prC * dwsr * (epsilon)/ewvsr^3 - sigx22 *
    sigx7/sigx2)/sigx2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (sigx28 * sigx5/(sigma_sq)) -
      (sigx20 * sigx7/sigx2 + S * prC * dwsr * (epsilon)/(wzdeno *
        ewvsr^3))) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((ewu * (S * (0.5 * sigx18 - (0.5 * (S^2 * sigx14 *
    dmusig * (epsilon)^2/(sigma_sq)^2) - (((0.5 * (ewu/(sigma_sq)) +
    1 - 0.5 * (0.5 * (1 - ewu/(sigma_sq)) + ewu/(sigma_sq))) *
    (1 - ewu/(sigma_sq)) * ewv/sigmastar + (2 - 2 * (sigx13^2 *
    ewu * (sigma_sq)/(ssq)^2)) * sigmastar)/(ssq)^2 + S^2 *
    sigx14^2 * ewu * (epsilon)^2/(ssq)) * dmusig) * depsi) *
    (epsilon)/wzdeno - 0.5 * sigx27) + S * ddmu * (epsilon)/wzdeno -
    0.5 * sigx4)/(sigx2 * sqrt(sigma_sq)) - (0.5 * (sigx2/sqrt(sigma_sq)) +
    2 * (ewz * sigx16)) * ewu * sigx16/(sigx2 * sqrt(sigma_sq))^2) *
    ewu * ewz), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * (ewv * (S * (((((0.5 *
      ((1 - ewu/(sigma_sq)) * ewv) - S^2 * sigx8 * sigx14 *
      ewu * (epsilon)^2)/(sigma_sq) + 0.5 * ((ewu/(sigma_sq) -
      1) * ewv/(sigma_sq) + 1 - 0.5 * ((1 - ewu/(sigma_sq)) *
      ewvsq))) * dmusig/sigmastar + 0.5 * (S^2 * sigx8 *
      dmusig * (epsilon)^2/(sigma_sq)^2)) * ewu + sigx8 *
      (1 - 2 * (sigx13 * ewu * ssq/(ssq)^2)) * dmusig) *
      depsi/(ssq)^2 + 0.5 * sigx18) * (epsilon)/wzdeno -
      (0.5 * sigx27 + 0.5 * (sigx29/(sigma_sq))))) - 2 *
      (sigx22 * sigx16/sigx2)) * ewu * ewz/(sigx2 * sqrt(sigma_sq)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx26 * pmusig - S *
      sigx28 * sigx14 * depsi * dmusig * (epsilon)) - 2 *
      (sigx20 * ewz * sigx16/(sigx2 * sqrt(sigma_sq)))) *
      ewu * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (prC * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr *
      (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * (S^2 * dwsr *
      (epsilon)^2/ewvsr^2) - 0.5 * dwsr))/ewvsr + 2 * (ewv *
      ewz * (S * ((((0.5 * (ewv/(sigma_sq)) - 0.5 * (0.5 *
      ewvsq + ewv/(sigma_sq))) * ewvsq + S^2 * sigx8^2 *
      ewu * ewv * (epsilon)^2/((ssq)^2 * (sigma_sq))) *
      dmusig * ewu/sigmastar + ((0.5 * (S^2 * dmusig *
      (epsilon)^2/(sigma_sq)^2) - 2 * (sigx8 * dmusig *
      ssq/(ssq)^2)) * ewv + dmusig) * sigx8) * depsi *
      ewu/(ssq)^2 + S * (0.5 * (ewv * (S * (sigx8 * depsi *
      ewu/(ssq)^2 + 0.5 * sigx17) * (epsilon) - 2 * (pmusig/(sigma_sq)))) +
      0.5 * pmusig) * dmusig * (epsilon)/(sigma_sq)^2) *
      (epsilon)/wzdeno - ((0.5 * (dmusig * pmusig) + 0.5 *
      (ewv * (S * (sigx10 + 0.5 * sigx3) * (epsilon) -
        wzdeno^2 * dmusig * pmusig/wzsq))) * wzdeno/wzsq +
      0.5 * (ewv * sigx29/(sigma_sq))))/sqrt(sigma_sq)) -
      sigx22^2/sigx2)/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (2 * ((sigx26 * pmusig + S * sigx8 * sigx28 * depsi *
      dmusig * ewu * (epsilon)/(ssq)^2) * ewv) - (sigx22 *
      sigx20/sigx2 + (0.5 * (S^2 * dwsr * (epsilon)^2/(wzdeno *
      ewvsr^3)) - 0.5 * (wzdeno * dwsr * ewvsr/(wzdeno *
      ewvsr)^2)) * prC)) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewvsr) +
      ewvsr/(wzdeno * ewvsr)^2) * dwsr - (sigx20^2/sigx2 +
      2 * ((2 - 2 * (wzdeno * (sigma_sq) * ewz/wzsq)) *
        dmusig * pmusig * sqrt(sigma_sq)/wzsq))) * ewz +
      2 * sigx19 - prC * dwsr/(wzdeno * ewvsr)) * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesszisfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewusr <- exp(Wu/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  mu <- 0
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm((mu * ewv - S * ewu * (epsilon))/(ssq))
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  dmu <- dnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  prV <- (1 - ewv/(sigma_sq))
  s2sig <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  dpmu <- (dmusig * muepsi^2 * pmusig/(sigma_sq)^2)
  pssq <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  ewusq <- (ewu * pwu/sqrt(sigma_sq))
  sigx1 <- (mu/(ssq) - s2sig * (mu * ewv - S * ewu * (epsilon))/(ssq)^2)
  sigx2 <- (dmu * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmusig)
  sigx12 <- (pwu * sqrt(sigma_sq))
  sigx3 <- (ewz2 * depsi/ewvsr + ewz1 * dmusig * pmusig/sigx12)
  sigx4 <- (ewz1 * sigx2/pssq + S * ewz2 * depsi * (epsilon)/ewvsr^3)
  sigx5 <- (dmu * dmusig * ewv/sigmastar - dmusig * muepsi *
    pmusig)
  sigx6 <- (ewusr * sigx12^2)
  sigx7 <- (sigx5/pssq - dmusig * dwu * pmusig * sqrt(sigma_sq)/sigx6)
  sigx8 <- (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
    sigmastar)
  sigx9 <- (sigx8 * (mu * ewv - S * ewu * (epsilon))/(ssq)^2 +
    S * (epsilon)/(ssq))
  sigx10 <- (0.5 * ewusq - 0.5 * (mu * dwu * sqrt(sigma_sq)/ewusr))
  sigx11 <- (0.5 * dpmu - sigx9 * dmu * dmusig)
  sigx13 <- (sigx11 * ewu/sigx12 - sigx10 * dmusig * pmusig/sigx12^2)
  sigx14 <- ((0.5 * dpmu + dmu * dmusig * sigx1)/pwu - 0.5 *
    (dmusig * pmusig * pwu/sigx12^2))
  sigx15 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx16 <- (sigx14 * ewz1 * ewv/sqrt(sigma_sq) + ewz2 * sigx15/ewvsr)
  sigx17 <- (dmusig * pmusig/sigx12 - depsi/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((muepsi^2/(sigma_sq) - 1) * pmusig +
      dmu * ewu * muepsi/(ssq)) * dmusig + dmu * (dmusig *
      muepsi - dmusig * (mu * ewv - S * ewu * (epsilon))/ewv) *
      ewu/(ssq)) * ewz1/pssq + ewz2 * depsi * (S^2 * (epsilon)^2/ewvsr^2 -
      1)/ewvsr^3 - sigx4^2/sigx3)/sigx3, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((muepsi^2/(sigma_sq) -
      2) * pmusig + dmu * ewu * muepsi/(ssq)) * dmusig *
      muepsi/(sigma_sq)^2) - (sigx9 * dmusig * muepsi/(sigma_sq) +
      (sigx8 * ewu/(ssq)^2 - (sigx9 * (mu * ewv - S * ewu *
        (epsilon))/ewv + 1/sigmastar)/(sigma_sq)) * dmusig) *
      dmu) * ewu/sigx12 - (sigx13 * sigx4/sigx3 + sigx10 *
      sigx2/((sigma_sq) * sigx12^2))) * ewz1/sigx3, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((((dmusig * muepsi - dmusig * (mu * ewv - S * ewu *
    (epsilon))/ewv) * sigx1/(sigma_sq) - s2sig * dmusig *
    ewu/(ssq)^2) * dmu + 0.5 * (((muepsi^2/(sigma_sq) - 2) *
    pmusig + dmu * ewu * muepsi/(ssq)) * dmusig * muepsi/(sigma_sq)^2))/pwu -
    0.5 * (sigx2 * pwu/((sigma_sq) * sigx12^2))) * ewz1 *
    ewv/sqrt(sigma_sq) + S * ewz2 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 -
    2) - 0.5) * depsi * (epsilon)/ewvsr^3 - sigx16 * sigx4/sigx3)/sigx3,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sigx2/pssq - S * depsi *
      (epsilon)/ewvsr^3)/(pi * sigx3 * ((Wz)^2 + 1)) -
      pi * sigx4 * ((Wz)^2 + 1) * sigx17/(pi * sigx3 *
        ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) * pmusig/(sigma_sq) -
      sigx9 * dmu) * ewu) + 0.5 * pmusig) * dmusig * muepsi^2/(sigma_sq)^2 -
      ((sigx9^2 * ewu * (mu * ewv - S * ewu * (epsilon))/(ssq) +
        ((0.5 * (ewu/(sigma_sq)) - 0.5 * (0.5 * (1 -
          ewu/(sigma_sq)) + ewu/(sigma_sq))) * (1 - ewu/(sigma_sq)) *
          ewv * (mu * ewv - S * ewu * (epsilon))/sigmastar -
          sigx8 * (2 * (sigx8 * (sigma_sq) * (mu * ewv -
          S * ewu * (epsilon)) * sigmastar/(ssq)^2) +
          2 * (S * (epsilon))) * ewu)/(ssq)^2) * dmusig +
        sigx9 * (0.5 * (dmusig * ewu * muepsi^2/(sigma_sq)^2) +
          dmusig)) * dmu)/sigx12 - sigx11 * sigx10/sigx12^2) *
      ewu - ((((0.5 * (((1 - 0.5 * (ewu/(sigma_sq))) *
      pwu - 0.5 * (mu * dwu/ewusr)) * ewu/sqrt(sigma_sq)) -
      0.5 * (mu * ((0.5 * (mu^2/ewusr^2) - 0.5) * sqrt(sigma_sq) +
        0.5 * (ewu/sqrt(sigma_sq))) * dwu/ewusr)) * dmusig +
      0.5 * (sigx10 * dmusig * ewu * muepsi^2/(sigma_sq)^2)) *
      pmusig - (sigx9 * dmu * ewu + 2 * (sigx10 * pmusig *
      pwu * sqrt(sigma_sq)/sigx12^2)) * sigx10 * dmusig)/sigx12^2 +
      sigx13^2 * ewz1/sigx3)) * ewz1/sigx3, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx9 * dmusig * (mu *
      ewv - S * ewu * (epsilon))/sigmastar + 0.5 * (dmusig *
      muepsi^2/(sigma_sq))) * sigx1/(sigma_sq) - ((0.5 *
      ((1 - ewu/(sigma_sq)) * ewv/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) -
      1) * ewv/(sigma_sq) + 1 - 0.5 * ((1 - ewu/(sigma_sq)) *
      prV))) * (mu * ewv - S * ewu * (epsilon))/sigmastar +
      mu * sigx8 - s2sig * (2 * (sigx8 * (sigma_sq) * (mu *
      ewv - S * ewu * (epsilon)) * sigmastar/(ssq)^2) +
      S * (epsilon))) * dmusig/(ssq)^2) * dmu + 0.5 * (((0.5 *
      (muepsi^2/(sigma_sq)) - 2) * pmusig/(sigma_sq) -
      sigx9 * dmu) * dmusig * muepsi^2/(sigma_sq)^2)) *
      ewu + 0.5 * (mu * (0.5 * dpmu + dmu * dmusig * sigx1) *
      dwu/(ewusr * pwu)))/pwu - (0.5 * ((sigx11 * ewu *
      pwu - (0.5 * (mu * dwu/ewusr) + 2 * (sigx10 * pwu^2 *
      sqrt(sigma_sq)/sigx12^2)) * dmusig * pmusig)/sigx12^2) +
      0.5 * (sigx14 * ewu/(sigma_sq)))) * ewv/sqrt(sigma_sq) -
      sigx16 * sigx13/sigx3) * ewz1/sigx3, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx13 * (1/(pi * sigx3 *
      ((Wz)^2 + 1)) - pi * ((Wz)^2 + 1) * ewz1 * sigx17/(pi *
      sigx3 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (dmusig * muepsi^2/(sigma_sq)) -
      dmusig * (mu * ewv - S * ewu * (epsilon)) * sigx1/sigmastar) *
      ewv * sigx1/(sigma_sq) + dmusig * (mu/(ssq) - (((3 *
      (mu) - 2 * (s2sig * (sigma_sq) * (mu * ewv - S *
      ewu * (epsilon)) * sigmastar/(ssq)^2)) * ewv - S *
      ewu * (epsilon)) * s2sig + (0.5 * (ewv/(sigma_sq)) -
      0.5 * (0.5 * prV + ewv/(sigma_sq))) * prV * ewu *
      (mu * ewv - S * ewu * (epsilon))/sigmastar)/(ssq)^2)) *
      dmu + (0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) *
      pmusig/(sigma_sq) + dmu * sigx1) * ewv) + 0.5 * pmusig) *
      dmusig * muepsi^2/(sigma_sq)^2)/pwu - ((0.5 * (((dmu *
      sigx1 - pmusig * pwu^2/sigx12^2) * dmusig + 0.5 *
      dpmu) * ewv) + 0.5 * (dmusig * pmusig)) * pwu/sigx12^2 +
      0.5 * (sigx14 * ewv/(sigma_sq)))) * ewz1 * ewv/sqrt(sigma_sq) +
      ewz2 * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
        1) - 0.25) * depsi * (epsilon)^2/ewvsr^2 - 0.5 *
        sigx15)/ewvsr - sigx16^2/sigx3)/sigx3, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((sigx14 * ewv/sqrt(sigma_sq) - sigx15/ewvsr)/(pi * sigx3 *
      ((Wz)^2 + 1)) - pi * sigx16 * ((Wz)^2 + 1) * sigx17/(pi *
      sigx3 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (pi * Wz * sigx3) +
      dmusig * pmusig/sigx12 - depsi/ewvsr) * sigx17/(pi *
      sigx3 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesszisfhalfnormlike_probit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewusr <- exp(Wu/2)
  mu <- 0
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm((mu * ewv - S * ewu * (epsilon))/(ssq))
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  dmu <- dnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  prV <- (1 - ewv/(sigma_sq))
  s2sig <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  dpmu <- (dmusig * muepsi^2 * pmusig/(sigma_sq)^2)
  pssq <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  ewusq <- (ewu * pwu/sqrt(sigma_sq))
  sigx1 <- (mu/(ssq) - s2sig * (mu * ewv - S * ewu * (epsilon))/(ssq)^2)
  sigx2 <- (dmu * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmusig)
  sigx12 <- (pwu * sqrt(sigma_sq))
  sigx3 <- ((1 - pwZ) * depsi/ewvsr + dmusig * pmusig * pwZ/sigx12)
  sigx4 <- (sigx2 * pwZ/pssq + S * (1 - pwZ) * depsi * (epsilon)/ewvsr^3)
  sigx5 <- (dmu * dmusig * ewv/sigmastar - dmusig * muepsi *
    pmusig)
  sigx6 <- (ewusr * sigx12^2)
  sigx7 <- (sigx5/pssq - dmusig * dwu * pmusig * sqrt(sigma_sq)/sigx6)
  sigx8 <- (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
    sigmastar)
  sigx9 <- (sigx8 * (mu * ewv - S * ewu * (epsilon))/(ssq)^2 +
    S * (epsilon)/(ssq))
  sigx10 <- (0.5 * ewusq - 0.5 * (mu * dwu * sqrt(sigma_sq)/ewusr))
  sigx11 <- (0.5 * dpmu - sigx9 * dmu * dmusig)
  sigx13 <- (sigx11 * ewu/sigx12 - sigx10 * dmusig * pmusig/sigx12^2)
  sigx14 <- ((0.5 * dpmu + dmu * dmusig * sigx1)/pwu - 0.5 *
    (dmusig * pmusig * pwu/sigx12^2))
  sigx15 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx16 <- (sigx14 * ewv * pwZ/sqrt(sigma_sq) + sigx15 * (1 -
    pwZ)/ewvsr)
  sigx17 <- (dmusig * pmusig/sigx12 - depsi/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((muepsi^2/(sigma_sq) - 1) * pmusig +
      dmu * ewu * muepsi/(ssq)) * dmusig + dmu * (dmusig *
      muepsi - dmusig * (mu * ewv - S * ewu * (epsilon))/ewv) *
      ewu/(ssq)) * pwZ/pssq + (1 - pwZ) * depsi * (S^2 *
      (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 - sigx4^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((muepsi^2/(sigma_sq) -
      2) * pmusig + dmu * ewu * muepsi/(ssq)) * dmusig *
      muepsi/(sigma_sq)^2) - (sigx9 * dmusig * muepsi/(sigma_sq) +
      (sigx8 * ewu/(ssq)^2 - (sigx9 * (mu * ewv - S * ewu *
        (epsilon))/ewv + 1/sigmastar)/(sigma_sq)) * dmusig) *
      dmu) * ewu/sigx12 - (sigx13 * sigx4/sigx3 + sigx10 *
      sigx2/((sigma_sq) * sigx12^2))) * pwZ/sigx3, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((((dmusig * muepsi - dmusig * (mu * ewv - S * ewu *
    (epsilon))/ewv) * sigx1/(sigma_sq) - s2sig * dmusig *
    ewu/(ssq)^2) * dmu + 0.5 * (((muepsi^2/(sigma_sq) - 2) *
    pmusig + dmu * ewu * muepsi/(ssq)) * dmusig * muepsi/(sigma_sq)^2))/pwu -
    0.5 * (sigx2 * pwu/((sigma_sq) * sigx12^2))) * ewv *
    pwZ/sqrt(sigma_sq) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 -
    2) - 0.5) * (1 - pwZ) * depsi * (epsilon)/ewvsr^3 - sigx16 *
    sigx4/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx2/pssq - (sigx4 *
      sigx17/sigx3 + S * depsi * (epsilon)/ewvsr^3)) *
      dwZ/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) * pmusig/(sigma_sq) -
      sigx9 * dmu) * ewu) + 0.5 * pmusig) * dmusig * muepsi^2/(sigma_sq)^2 -
      ((sigx9^2 * ewu * (mu * ewv - S * ewu * (epsilon))/(ssq) +
        ((0.5 * (ewu/(sigma_sq)) - 0.5 * (0.5 * (1 -
          ewu/(sigma_sq)) + ewu/(sigma_sq))) * (1 - ewu/(sigma_sq)) *
          ewv * (mu * ewv - S * ewu * (epsilon))/sigmastar -
          sigx8 * (2 * (sigx8 * (sigma_sq) * (mu * ewv -
          S * ewu * (epsilon)) * sigmastar/(ssq)^2) +
          2 * (S * (epsilon))) * ewu)/(ssq)^2) * dmusig +
        sigx9 * (0.5 * (dmusig * ewu * muepsi^2/(sigma_sq)^2) +
          dmusig)) * dmu)/sigx12 - sigx11 * sigx10/sigx12^2) *
      ewu - ((((0.5 * (((1 - 0.5 * (ewu/(sigma_sq))) *
      pwu - 0.5 * (mu * dwu/ewusr)) * ewu/sqrt(sigma_sq)) -
      0.5 * (mu * ((0.5 * (mu^2/ewusr^2) - 0.5) * sqrt(sigma_sq) +
        0.5 * (ewu/sqrt(sigma_sq))) * dwu/ewusr)) * dmusig +
      0.5 * (sigx10 * dmusig * ewu * muepsi^2/(sigma_sq)^2)) *
      pmusig - (sigx9 * dmu * ewu + 2 * (sigx10 * pmusig *
      pwu * sqrt(sigma_sq)/sigx12^2)) * sigx10 * dmusig)/sigx12^2 +
      sigx13^2 * pwZ/sigx3)) * pwZ/sigx3, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx9 * dmusig * (mu *
      ewv - S * ewu * (epsilon))/sigmastar + 0.5 * (dmusig *
      muepsi^2/(sigma_sq))) * sigx1/(sigma_sq) - ((0.5 *
      ((1 - ewu/(sigma_sq)) * ewv/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) -
      1) * ewv/(sigma_sq) + 1 - 0.5 * ((1 - ewu/(sigma_sq)) *
      prV))) * (mu * ewv - S * ewu * (epsilon))/sigmastar +
      mu * sigx8 - s2sig * (2 * (sigx8 * (sigma_sq) * (mu *
      ewv - S * ewu * (epsilon)) * sigmastar/(ssq)^2) +
      S * (epsilon))) * dmusig/(ssq)^2) * dmu + 0.5 * (((0.5 *
      (muepsi^2/(sigma_sq)) - 2) * pmusig/(sigma_sq) -
      sigx9 * dmu) * dmusig * muepsi^2/(sigma_sq)^2)) *
      ewu + 0.5 * (mu * (0.5 * dpmu + dmu * dmusig * sigx1) *
      dwu/(ewusr * pwu)))/pwu - (0.5 * ((sigx11 * ewu *
      pwu - (0.5 * (mu * dwu/ewusr) + 2 * (sigx10 * pwu^2 *
      sqrt(sigma_sq)/sigx12^2)) * dmusig * pmusig)/sigx12^2) +
      0.5 * (sigx14 * ewu/(sigma_sq)))) * ewv/sqrt(sigma_sq) -
      sigx16 * sigx13/sigx3) * pwZ/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx13 * (1 - sigx17 * pwZ/sigx3) *
      dwZ/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (dmusig * muepsi^2/(sigma_sq)) -
      dmusig * (mu * ewv - S * ewu * (epsilon)) * sigx1/sigmastar) *
      ewv * sigx1/(sigma_sq) + dmusig * (mu/(ssq) - (((3 *
      (mu) - 2 * (s2sig * (sigma_sq) * (mu * ewv - S *
      ewu * (epsilon)) * sigmastar/(ssq)^2)) * ewv - S *
      ewu * (epsilon)) * s2sig + (0.5 * (ewv/(sigma_sq)) -
      0.5 * (0.5 * prV + ewv/(sigma_sq))) * prV * ewu *
      (mu * ewv - S * ewu * (epsilon))/sigmastar)/(ssq)^2)) *
      dmu + (0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) *
      pmusig/(sigma_sq) + dmu * sigx1) * ewv) + 0.5 * pmusig) *
      dmusig * muepsi^2/(sigma_sq)^2)/pwu - ((0.5 * (((dmu *
      sigx1 - pmusig * pwu^2/sigx12^2) * dmusig + 0.5 *
      dpmu) * ewv) + 0.5 * (dmusig * pmusig)) * pwu/sigx12^2 +
      0.5 * (sigx14 * ewv/(sigma_sq)))) * ewv * pwZ/sqrt(sigma_sq) +
      (1 - pwZ) * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
        1) - 0.25) * depsi * (epsilon)^2/ewvsr^2 - 0.5 *
        sigx15)/ewvsr - sigx16^2/sigx3)/sigx3, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (sigx14 * ewv/sqrt(sigma_sq) - (sigx16 * sigx17/sigx3 +
      sigx15/ewvsr)) * dwZ/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx17 * dwZ/sigx3 + Wz) *
      sigx17 * dwZ/sigx3), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesszisfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewusr <- exp(Wu/2)
  mu <- 0
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm((mu * ewv - S * ewu * (epsilon))/(ssq))
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  dmu <- dnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  prV <- (1 - ewv/(sigma_sq))
  s2sig <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  dpmu <- (dmusig * muepsi^2 * pmusig/(sigma_sq)^2)
  pssq <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  ewusq <- (ewu * pwu/sqrt(sigma_sq))
  sigx1 <- (mu/(ssq) - s2sig * (mu * ewv - S * ewu * (epsilon))/(ssq)^2)
  sigx2 <- (dmu * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmusig)
  sigx12 <- (pwu * sqrt(sigma_sq))
  sigx3 <- ((1 - prZ) * dmusig * pmusig/sigx12 + depsi * prZ/ewvsr)
  sigx4 <- ((1 - prZ) * sigx2/pssq + S * depsi * prZ * (epsilon)/ewvsr^3)
  sigx5 <- (dmu * dmusig * ewv/sigmastar - dmusig * muepsi *
    pmusig)
  sigx6 <- (ewusr * sigx12^2)
  sigx7 <- (sigx5/pssq - dmusig * dwu * pmusig * sqrt(sigma_sq)/sigx6)
  sigx8 <- (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
    sigmastar)
  sigx9 <- (sigx8 * (mu * ewv - S * ewu * (epsilon))/(ssq)^2 +
    S * (epsilon)/(ssq))
  sigx10 <- (0.5 * ewusq - 0.5 * (mu * dwu * sqrt(sigma_sq)/ewusr))
  sigx11 <- (0.5 * dpmu - sigx9 * dmu * dmusig)
  sigx13 <- (sigx11 * ewu/sigx12 - sigx10 * dmusig * pmusig/sigx12^2)
  sigx14 <- ((0.5 * dpmu + dmu * dmusig * sigx1)/pwu - 0.5 *
    (dmusig * pmusig * pwu/sigx12^2))
  sigx15 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx16 <- (sigx14 * (1 - prZ) * ewv/sqrt(sigma_sq) + sigx15 *
    prZ/ewvsr)
  sigx17 <- (dmusig * pmusig/sigx12 - depsi/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((muepsi^2/(sigma_sq) - 1) * pmusig +
      dmu * ewu * muepsi/(ssq)) * dmusig + dmu * (dmusig *
      muepsi - dmusig * (mu * ewv - S * ewu * (epsilon))/ewv) *
      ewu/(ssq)) * (1 - prZ)/pssq + depsi * prZ * (S^2 *
      (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 - sigx4^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((muepsi^2/(sigma_sq) -
      2) * pmusig + dmu * ewu * muepsi/(ssq)) * dmusig *
      muepsi/(sigma_sq)^2) - (sigx9 * dmusig * muepsi/(sigma_sq) +
      (sigx8 * ewu/(ssq)^2 - (sigx9 * (mu * ewv - S * ewu *
        (epsilon))/ewv + 1/sigmastar)/(sigma_sq)) * dmusig) *
      dmu) * ewu/sigx12 - (sigx13 * sigx4/sigx3 + sigx10 *
      sigx2/((sigma_sq) * sigx12^2))) * (1 - prZ)/sigx3,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((((dmusig * muepsi - dmusig * (mu * ewv - S * ewu *
    (epsilon))/ewv) * sigx1/(sigma_sq) - s2sig * dmusig *
    ewu/(ssq)^2) * dmu + 0.5 * (((muepsi^2/(sigma_sq) - 2) *
    pmusig + dmu * ewu * muepsi/(ssq)) * dmusig * muepsi/(sigma_sq)^2))/pwu -
    0.5 * (sigx2 * pwu/((sigma_sq) * sigx12^2))) * (1 - prZ) *
    ewv/sqrt(sigma_sq) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 -
    2) - 0.5) * depsi * prZ * (epsilon)/ewvsr^3 - sigx16 *
    sigx4/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx2/pssq - (sigx4 *
      sigx17/sigx3 + S * depsi * (epsilon)/ewvsr^3)) *
      prZ * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) * pmusig/(sigma_sq) -
      sigx9 * dmu) * ewu) + 0.5 * pmusig) * dmusig * muepsi^2/(sigma_sq)^2 -
      ((sigx9^2 * ewu * (mu * ewv - S * ewu * (epsilon))/(ssq) +
        ((0.5 * (ewu/(sigma_sq)) - 0.5 * (0.5 * (1 -
          ewu/(sigma_sq)) + ewu/(sigma_sq))) * (1 - ewu/(sigma_sq)) *
          ewv * (mu * ewv - S * ewu * (epsilon))/sigmastar -
          sigx8 * (2 * (sigx8 * (sigma_sq) * (mu * ewv -
          S * ewu * (epsilon)) * sigmastar/(ssq)^2) +
          2 * (S * (epsilon))) * ewu)/(ssq)^2) * dmusig +
        sigx9 * (0.5 * (dmusig * ewu * muepsi^2/(sigma_sq)^2) +
          dmusig)) * dmu)/sigx12 - sigx11 * sigx10/sigx12^2) *
      ewu - ((((0.5 * (((1 - 0.5 * (ewu/(sigma_sq))) *
      pwu - 0.5 * (mu * dwu/ewusr)) * ewu/sqrt(sigma_sq)) -
      0.5 * (mu * ((0.5 * (mu^2/ewusr^2) - 0.5) * sqrt(sigma_sq) +
        0.5 * (ewu/sqrt(sigma_sq))) * dwu/ewusr)) * dmusig +
      0.5 * (sigx10 * dmusig * ewu * muepsi^2/(sigma_sq)^2)) *
      pmusig - (sigx9 * dmu * ewu + 2 * (sigx10 * pmusig *
      pwu * sqrt(sigma_sq)/sigx12^2)) * sigx10 * dmusig)/sigx12^2 +
      sigx13^2 * (1 - prZ)/sigx3)) * (1 - prZ)/sigx3, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx9 * dmusig * (mu *
      ewv - S * ewu * (epsilon))/sigmastar + 0.5 * (dmusig *
      muepsi^2/(sigma_sq))) * sigx1/(sigma_sq) - ((0.5 *
      ((1 - ewu/(sigma_sq)) * ewv/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) -
      1) * ewv/(sigma_sq) + 1 - 0.5 * ((1 - ewu/(sigma_sq)) *
      prV))) * (mu * ewv - S * ewu * (epsilon))/sigmastar +
      mu * sigx8 - s2sig * (2 * (sigx8 * (sigma_sq) * (mu *
      ewv - S * ewu * (epsilon)) * sigmastar/(ssq)^2) +
      S * (epsilon))) * dmusig/(ssq)^2) * dmu + 0.5 * (((0.5 *
      (muepsi^2/(sigma_sq)) - 2) * pmusig/(sigma_sq) -
      sigx9 * dmu) * dmusig * muepsi^2/(sigma_sq)^2)) *
      ewu + 0.5 * (mu * (0.5 * dpmu + dmu * dmusig * sigx1) *
      dwu/(ewusr * pwu)))/pwu - (0.5 * ((sigx11 * ewu *
      pwu - (0.5 * (mu * dwu/ewusr) + 2 * (sigx10 * pwu^2 *
      sqrt(sigma_sq)/sigx12^2)) * dmusig * pmusig)/sigx12^2) +
      0.5 * (sigx14 * ewu/(sigma_sq)))) * ewv/sqrt(sigma_sq) -
      sigx16 * sigx13/sigx3) * (1 - prZ)/sigx3, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx13 * (1 - (1 - prZ) *
      sigx17/sigx3) * prZ * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (dmusig * muepsi^2/(sigma_sq)) -
      dmusig * (mu * ewv - S * ewu * (epsilon)) * sigx1/sigmastar) *
      ewv * sigx1/(sigma_sq) + dmusig * (mu/(ssq) - (((3 *
      (mu) - 2 * (s2sig * (sigma_sq) * (mu * ewv - S *
      ewu * (epsilon)) * sigmastar/(ssq)^2)) * ewv - S *
      ewu * (epsilon)) * s2sig + (0.5 * (ewv/(sigma_sq)) -
      0.5 * (0.5 * prV + ewv/(sigma_sq))) * prV * ewu *
      (mu * ewv - S * ewu * (epsilon))/sigmastar)/(ssq)^2)) *
      dmu + (0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) *
      pmusig/(sigma_sq) + dmu * sigx1) * ewv) + 0.5 * pmusig) *
      dmusig * muepsi^2/(sigma_sq)^2)/pwu - ((0.5 * (((dmu *
      sigx1 - pmusig * pwu^2/sigx12^2) * dmusig + 0.5 *
      dpmu) * ewv) + 0.5 * (dmusig * pmusig)) * pwu/sigx12^2 +
      0.5 * (sigx14 * ewv/(sigma_sq)))) * (1 - prZ) * ewv/sqrt(sigma_sq) +
      prZ * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
        1) - 0.25) * depsi * (epsilon)^2/ewvsr^2 - 0.5 *
        sigx15)/ewvsr - sigx16^2/sigx3)/sigx3, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (sigx14 * ewv/sqrt(sigma_sq) - (sigx16 * sigx17/sigx3 +
      sigx15/ewvsr)) * prZ * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (sigx17 * prZ/sigx3 +
      1) * ewz) * sigx17 * prZ * ewz/sigx3, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Different sigma_v

## logit specification class membership
chessmnsfhalfnormlike_logit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewv1 <- exp(Wv1)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  ssq <- ((sigma_sq) * sigmastar)
  dmusig <- dnorm(-(S * ewu * (epsilon)/ssq))
  pmusig <- pnorm(-(S * ewu * (epsilon)/ssq))
  depsi <- dnorm(S * (epsilon)/sqrt(sigma_sq))
  dwsr <- dnorm(S * (epsilon)/ewv2_h)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- (wzdeno * depsi * pmusig/(wzdeno * sqrt(sigma_sq))^2)
  sigx2 <- (S * depsi * pmusig * (epsilon)/(sigma_sq)^2)
  wsqsq <- (wzdeno * (sigma_sq) * sqrt(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx3 <- (0.5 * (prU * ewv1/sigmastar) + sigmastar)
  sigx4 <- (0.5 * sigx2 - (1/ssq - sigx3 * ewu/ssq^2) * dmusig *
    depsi)
  sigx5 <- (S * sigx4 * (epsilon)/wzdeno - 0.5 * sigx1)
  sigx6 <- (depsi * ewz * pmusig/(wzdeno * sqrt(sigma_sq)))
  sigx7 <- (prC * dwsr/ewv2_h + 2 * sigx6)
  sigx8 <- (sigx7 * sqrt(sigma_sq))
  prV <- (1 - ewv1/(sigma_sq))
  sigx9 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx10 <- (sigx9 * dmusig * depsi * ewu/ssq^2 + 0.5 * sigx2)
  sigx11 <- (S * sigx10 * (epsilon)/wzdeno - 0.5 * sigx1)
  sigx12 <- (S^2 * dwsr * (epsilon)^2/ewv2_h^2)
  sigx13 <- (0.5 * sigx12 - 0.5 * dwsr)
  sigx14 <- (1/(wzdeno * sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/(wzdeno *
    sqrt(sigma_sq))^2)
  sigx15 <- (2 * (sigx14 * depsi * pmusig) - prC * dwsr/(wzdeno *
    ewv2_h))
  sigx16 <- (dmusig * depsi * ewu/sigmastar + S * depsi * pmusig *
    (epsilon))
  sigx17 <- (sigx16 * ewz/wsqsq)
  sigx18 <- (2 * sigx17 + S * prC * dwsr * (epsilon)/ewv2_h^3)
  sigx19 <- (dmusig * ewu/sigmastar + S * pmusig * (epsilon))
  sigx20 <- (S * sigx19 * (epsilon)/(sigma_sq) - pmusig)
  sigx21 <- (1/ssq - sigx3 * ewu/ssq^2)
  sigx22 <- (depsi * ewu/ewv1 + depsi)
  sigx23 <- (0.5 * sigx20 - 0.5 * pmusig)
  sigx24 <- ((wzdeno * sqrt(sigma_sq))^2 * (sigma_sq))
  sigx25 <- (wzdeno * sigx16/sigx24)
  sigx26 <- (S * pmusig * (epsilon)/(sigma_sq)^2)
  sigx27 <- (S * depsi * (S * (0.5 * sigx26 - sigx21 * dmusig) *
    (epsilon) - 2 * (pmusig/(sigma_sq))) * (epsilon)/(sigma_sq)^2)
  sigx28 <- (S * sigx4 * (epsilon) - wzdeno^2 * depsi * pmusig/(wzdeno *
    sqrt(sigma_sq))^2)
  sigx29 <- (wzdeno * sigx28/(wzdeno * sqrt(sigma_sq))^2)
  sigx30 <- (0.5 * (sigx7/sqrt(sigma_sq)) + 2 * (ewz * sigx5))
  sigx31 <- (S^2 * sigx14 * depsi * (epsilon)^2/(sigma_sq)^2)
  sigx32 <- (0.5/sqrt(sigma_sq) - wzdeno^2 * sqrt(sigma_sq)/(wzdeno *
    sqrt(sigma_sq))^2)
  sigx33 <- (sigx32 * ewz + 0.5 * (wzdeno/sqrt(sigma_sq)))
  sigx34 <- (0.5 * sigx31 - sigx33 * depsi/(wzdeno * sqrt(sigma_sq))^2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (prC * dwsr * (S^2 * (epsilon)^2/ewv2_h^2 -
      1)/ewv2_h^3 + 2 * ((depsi * sigx20 + S * dmusig *
      sigx22 * ewu * (epsilon)/ssq) * ewz/wsqsq) - sigx18^2/sigx7)/sigx7,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((sigx21 * dmusig *
      depsi + S * (sigx23 * depsi/(sigma_sq) - S * sigx21 *
      dmusig * sigx22 * (epsilon)) * (epsilon)/(sigma_sq))/wzdeno -
      0.5 * sigx25)/sigx8 - sigx18 * sigx5 * sqrt(sigma_sq)/sigx8^2) *
      ewu * ewz), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * (((S * (sigx23 * depsi/(sigma_sq) + S * sigx9 *
    dmusig * sigx22 * ewu * (epsilon)/ssq^2) * (epsilon)/(sigma_sq) -
    sigx9 * dmusig * depsi * ewu/ssq^2)/wzdeno - 0.5 * sigx25)/sigx8 -
    sigx18 * sigx11 * sqrt(sigma_sq)/sigx8^2) * ewv1 * ewz),
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * prC * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * dwsr * (epsilon)/(sigx7 * ewv2_h^3) -
      sigx13 * sigx18 * ewv2_h/(sigx7 * ewv2_h)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (sigx14 * sigx16/(sigma_sq)) -
      (sigx15 * sigx18/sigx7 + S * prC * dwsr * (epsilon)/(wzdeno *
        ewv2_h^3))) * ewz/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((ewu * (S * (0.5 * sigx27 - (0.5 * (S^2 * sigx21 *
    depsi * (epsilon)^2/(sigma_sq)^2) - (((0.5 * (ewu/(sigma_sq)) +
    1 - 0.5 * (0.5 * prU + ewu/(sigma_sq))) * prU * ewv1/sigmastar +
    (2 - 2 * (sigx3^2 * ewu * (sigma_sq)/ssq^2)) * sigmastar)/ssq^2 +
    S^2 * sigx21^2 * ewu * (epsilon)^2/ssq) * depsi) * dmusig) *
    (epsilon)/wzdeno - 0.5 * sigx29) + S * sigx4 * (epsilon)/wzdeno -
    0.5 * sigx1)/sigx8 - sigx30 * ewu * sigx5/sigx8^2) *
    ewu * ewz), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * (((((0.5 * (prU *
      ewv1) - S^2 * sigx9 * sigx21 * ewu * (epsilon)^2)/(sigma_sq) +
      0.5 * ((ewu/(sigma_sq) - 1) * ewv1/(sigma_sq) + 1 -
        0.5 * (prU * prV))) * depsi/sigmastar + 0.5 *
      (S^2 * sigx9 * depsi * (epsilon)^2/(sigma_sq)^2)) *
      ewu + sigx9 * (1 - 2 * (sigx3 * ewu * (sigma_sq) *
      sigmastar/ssq^2)) * depsi) * dmusig/ssq^2 + 0.5 *
      sigx27) * (epsilon)/wzdeno - 0.5 * sigx29)/sigx8 -
      sigx30 * sigx11/sigx8^2) * ewu * ewv1 * ewz), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (sigx13 * prC * ewu *
      ewv2_h * ewz * sigx5/((sigx7 * ewv2_h)^2 * sqrt(sigma_sq)))),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx34 * pmusig - S *
      sigx14 * sigx21 * dmusig * depsi * (epsilon)) - 2 *
      (sigx15 * ewz * sigx5/sigx8)) * ewu * ewz/sigx7,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * ((((0.5 * (ewv1/(sigma_sq)) -
      0.5 * (0.5 * prV + ewv1/(sigma_sq))) * prV + S^2 *
      sigx9^2 * ewu * ewv1 * (epsilon)^2/(ssq^2 * (sigma_sq))) *
      depsi * ewu/sigmastar + ((0.5 * (S^2 * depsi * (epsilon)^2/(sigma_sq)^2) -
      2 * (sigx9 * depsi * (sigma_sq) * sigmastar/ssq^2)) *
      ewv1 + depsi) * sigx9) * dmusig * ewu/ssq^2 + S *
      (0.5 * (ewv1 * (S * (sigx9 * dmusig * ewu/ssq^2 +
        0.5 * sigx26) * (epsilon) - 2 * (pmusig/(sigma_sq)))) +
        0.5 * pmusig) * depsi * (epsilon)/(sigma_sq)^2) *
      (epsilon)/wzdeno - (0.5 * (depsi * pmusig) + 0.5 *
      (ewv1 * (S * sigx10 * (epsilon) - wzdeno^2 * depsi *
        pmusig/(wzdeno * sqrt(sigma_sq))^2))) * wzdeno/(wzdeno *
      sqrt(sigma_sq))^2)/sigx8 - (0.5 * (sigx7/sqrt(sigma_sq)) +
      2 * (ewz * sigx11)) * ewv1 * sigx11/sigx8^2) * ewv1 *
      ewz), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (2 * (sigx13 * prC * ewv1 * ewv2_h * ewz * sigx11/((sigx7 *
      ewv2_h)^2 * sqrt(sigma_sq)))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx34 * pmusig + S *
      sigx9 * sigx14 * dmusig * depsi * ewu * (epsilon)/ssq^2) -
      2 * (sigx15 * ewz * sigx11/sigx8)) * ewv1 * ewz/sigx7,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * dwsr * (epsilon)^2/(sigx7 * ewv2_h^3) -
      (sigx13 * prC + 0.5 * (sigx7 * ewv2_h)) * sigx13/(sigx7 *
        ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx13 * sigx15/sigx7 +
      0.5 * (S^2 * dwsr * (epsilon)^2/(wzdeno * ewv2_h^2)))/ewv2_h -
      0.5 * (wzdeno * dwsr * ewv2_h/(wzdeno * ewv2_h)^2)) *
      prC * ewz/sigx7), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv2_h) +
      ewv2_h/(wzdeno * ewv2_h)^2) * dwsr - (sigx15^2/sigx7 +
      2 * ((2 - 2 * (wzdeno * (sigma_sq) * ewz/(wzdeno *
        sqrt(sigma_sq))^2)) * depsi * pmusig * sqrt(sigma_sq)/(wzdeno *
        sqrt(sigma_sq))^2))) * ewz + 2 * (sigx14 * depsi *
      pmusig) - prC * dwsr/(wzdeno * ewv2_h)) * ewz/sigx7,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmnsfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  mu <- 0
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  ssq <- ((sigma_sq) * sigmastar)
  muvu <- (mu * ewv1 - S * ewu * (epsilon))
  mustar <- muvu/ssq
  dmustar <- dnorm(mustar, 0, 1)
  pmustar <- pnorm(mustar)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  depsi <- dnorm(S * (epsilon)/ewvsr2)
  musq <- (mu * dwu * sqrt(sigma_sq)/ewusr)
  prV <- (1 - ewv1/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx1 <- (dmustar * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmustar)
  sigx2 <- (ewz2 * depsi/ewvsr2 + ewz1 * dmusig * pmustar/(pwu *
    sqrt(sigma_sq)))
  sigx3 <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  sigx4 <- (ewz1 * sigx1/sigx3 + S * ewz2 * depsi * (epsilon)/ewvsr2^3)
  sigx5 <- (dmustar * dmusig * ewv1/sigmastar - dmusig * muepsi *
    pmustar)
  sigx6 <- (ewusr * (pwu * sqrt(sigma_sq))^2)
  sigx7 <- (sigx5/sigx3 - dmusig * dwu * pmustar * sqrt(sigma_sq)/sigx6)
  sigx8 <- (dmusig * muepsi^2 * pmustar/(sigma_sq)^2)
  sigx9 <- (0.5 * (prU * ewv1/sigmastar) + sigmastar)
  sigx10 <- (sigx9 * muvu/ssq^2 + S * (epsilon)/ssq)
  sigx11 <- (0.5 * sigx8 - sigx10 * dmustar * dmusig)
  sigx12 <- (0.5 * (ewu * pwu/sqrt(sigma_sq)) - 0.5 * musq)
  sigx13 <- (sigx11 * ewu/(pwu * sqrt(sigma_sq)) - sigx12 *
    dmusig * pmustar/(pwu * sqrt(sigma_sq))^2)
  sigx14 <- (mu/ssq - (0.5 * (prV * ewu/sigmastar) + sigmastar) *
    muvu/ssq^2)
  sigx15 <- (0.5 * sigx8 + dmustar * dmusig * sigx14)
  sigx16 <- (dmusig * pmustar * pwu/(pwu * sqrt(sigma_sq))^2)
  sigx17 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi)
  sigx18 <- (dmusig * pmustar/(pwu * sqrt(sigma_sq)) - depsi/ewvsr2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((muepsi^2/(sigma_sq) - 1) * pmustar +
      dmustar * ewu * muepsi/ssq) * dmusig + dmustar *
      (dmusig * muepsi - dmusig * muvu/ewv1) * ewu/ssq) *
      ewz1/sigx3 + ewz2 * depsi * (S^2 * (epsilon)^2/ewvsr2^2 -
      1)/ewvsr2^3 - sigx4^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((muepsi^2/(sigma_sq) -
      2) * pmustar + dmustar * ewu * muepsi/ssq) * dmusig *
      muepsi/(sigma_sq)^2) - (sigx10 * dmusig * muepsi/(sigma_sq) +
      (sigx9 * ewu/ssq^2 - (sigx10 * muvu/ewv1 + 1/sigmastar)/(sigma_sq)) *
        dmusig) * dmustar) * ewu/(pwu * sqrt(sigma_sq)) -
      (sigx13 * sigx4/sigx2 + sigx12 * sigx1/((sigma_sq) *
        (pwu * sqrt(sigma_sq))^2))) * ewz1/sigx2, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((((dmusig * muepsi - dmusig * muvu/ewv1) * sigx14/(sigma_sq) -
    (0.5 * (prV * ewu/sigmastar) + sigmastar) * dmusig *
      ewu/ssq^2) * dmustar + 0.5 * (((muepsi^2/(sigma_sq) -
    2) * pmustar + dmustar * ewu * muepsi/ssq) * dmusig *
    muepsi/(sigma_sq)^2))/pwu - 0.5 * (sigx1 * pwu/((sigma_sq) *
    (pwu * sqrt(sigma_sq))^2)))/(sigx2 * sqrt(sigma_sq)) -
    (sigx15/pwu - 0.5 * sigx16) * sigx4 * sqrt(sigma_sq)/(sigx2 *
      sqrt(sigma_sq))^2) * ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ewz2 * (S * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx2 * ewvsr2^3) -
      sigx4 * sigx17 * ewvsr2/(sigx2 * ewvsr2)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sigx1/sigx3 - S * depsi *
      (epsilon)/ewvsr2^3)/(pi * sigx2 * ((Wz)^2 + 1)) -
      pi * sigx4 * ((Wz)^2 + 1) * sigx18/(pi * sigx2 *
        ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) * pmustar/(sigma_sq) -
      sigx10 * dmustar) * ewu) + 0.5 * pmustar) * dmusig *
      muepsi^2/(sigma_sq)^2 - ((sigx10^2 * ewu * muvu/ssq +
      ((0.5 * (ewu/(sigma_sq)) - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv1 * muvu/sigmastar - sigx9 * (2 * (sigx9 *
        (sigma_sq) * muvu * sigmastar/ssq^2) + 2 * (S *
        (epsilon))) * ewu)/ssq^2) * dmusig + sigx10 *
      (0.5 * (dmusig * ewu * muepsi^2/(sigma_sq)^2) + dmusig)) *
      dmustar)/(pwu * sqrt(sigma_sq)) - sigx11 * sigx12/(pwu *
      sqrt(sigma_sq))^2) * ewu - ((((0.5 * (((1 - 0.5 *
      (ewu/(sigma_sq))) * pwu - 0.5 * (mu * dwu/ewusr)) *
      ewu/sqrt(sigma_sq)) - 0.5 * (mu * ((0.5 * (mu^2/ewusr^2) -
      0.5) * sqrt(sigma_sq) + 0.5 * (ewu/sqrt(sigma_sq))) *
      dwu/ewusr)) * dmusig + 0.5 * (sigx12 * dmusig * ewu *
      muepsi^2/(sigma_sq)^2)) * pmustar - (sigx10 * dmustar *
      ewu + 2 * (sigx12 * pmustar * pwu * sqrt(sigma_sq)/(pwu *
      sqrt(sigma_sq))^2)) * sigx12 * dmusig)/(pwu * sqrt(sigma_sq))^2 +
      sigx13^2 * ewz1/sigx2)) * ewz1/sigx2, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx10 * dmusig * muvu/sigmastar +
      0.5 * (dmusig * muepsi^2/(sigma_sq))) * sigx14/(sigma_sq) -
      ((0.5 * (prU * ewv1/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) -
        1) * ewv1/(sigma_sq) + 1 - 0.5 * (prU * prV))) *
        muvu/sigmastar + mu * sigx9 - (0.5 * (prV * ewu/sigmastar) +
        sigmastar) * (2 * (sigx9 * (sigma_sq) * muvu *
        sigmastar/ssq^2) + S * (epsilon))) * dmusig/ssq^2) *
      dmustar + 0.5 * (((0.5 * (muepsi^2/(sigma_sq)) -
      2) * pmustar/(sigma_sq) - sigx10 * dmustar) * dmusig *
      muepsi^2/(sigma_sq)^2)) * ewu + 0.5 * (mu * sigx15 *
      dwu/(ewusr * pwu)))/pwu - 0.5 * ((sigx11 * ewu *
      pwu - (0.5 * (mu * dwu/ewusr) + 2 * (sigx12 * pwu^2 *
      sqrt(sigma_sq)/(pwu * sqrt(sigma_sq))^2)) * dmusig *
      pmustar)/(pwu * sqrt(sigma_sq))^2))/(sigx2 * sqrt(sigma_sq)) -
      (sigx13 * ewz1 * sqrt(sigma_sq) + 0.5 * (sigx2 *
        ewu/sqrt(sigma_sq))) * (sigx15/pwu - 0.5 * sigx16)/(sigx2 *
        sqrt(sigma_sq))^2) * ewz1 * ewv1, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx13 * ewz2 * sigx17 *
      ewz1 * ewvsr2/(sigx2 * ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx13 * (1/(pi * sigx2 *
      ((Wz)^2 + 1)) - pi * ((Wz)^2 + 1) * ewz1 * sigx18/(pi *
      sigx2 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (dmusig * muepsi^2/(sigma_sq)) -
      dmusig * muvu * sigx14/sigmastar) * ewv1 * sigx14/(sigma_sq) +
      dmusig * (mu/ssq - (((3 * (mu) - 2 * ((0.5 * (prV *
        ewu/sigmastar) + sigmastar) * (sigma_sq) * muvu *
        sigmastar/ssq^2)) * ewv1 - S * ewu * (epsilon)) *
        (0.5 * (prV * ewu/sigmastar) + sigmastar) + (0.5 *
        (ewv1/(sigma_sq)) - 0.5 * (0.5 * prV + ewv1/(sigma_sq))) *
        prV * ewu * muvu/sigmastar)/ssq^2)) * dmustar +
      (0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) * pmustar/(sigma_sq) +
        dmustar * sigx14) * ewv1) + 0.5 * pmustar) *
        dmusig * muepsi^2/(sigma_sq)^2)/pwu - (0.5 *
      (((dmustar * sigx14 - pmustar * pwu^2/(pwu * sqrt(sigma_sq))^2) *
        dmusig + 0.5 * sigx8) * ewv1) + 0.5 * (dmusig *
      pmustar)) * pwu/(pwu * sqrt(sigma_sq))^2)/(sigx2 *
      sqrt(sigma_sq)) - ((sigx15/pwu - 0.5 * sigx16) *
      ewz1 + 0.5 * (sigx2/sqrt(sigma_sq))) * (sigx15/pwu -
      0.5 * sigx16) * ewv1/(sigx2 * sqrt(sigma_sq))^2) *
      ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx15/pwu - 0.5 * sigx16) * ewz2 * sigx17 * ewz1 *
      ewv1 * ewvsr2/((sigx2 * ewvsr2)^2 * sqrt(sigma_sq))),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx15/pwu - 0.5 * sigx16) *
      (1/(pi * sigx2 * ((Wz)^2 + 1)) - pi * ((Wz)^2 + 1) *
        ewz1 * sigx18/(pi * sigx2 * ((Wz)^2 + 1))^2) *
      ewv1/sqrt(sigma_sq), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ewz2 * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx2 * ewvsr2^3) -
      (ewz2 * sigx17 + 0.5 * (sigx2 * ewvsr2)) * sigx17/(sigx2 *
        ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx17 * (1/(pi * sigx2 *
      ((Wz)^2 + 1)) + pi * ((Wz)^2 + 1) * ewz2 * sigx18/(pi *
      sigx2 * ((Wz)^2 + 1))^2)/ewvsr2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (pi * Wz * sigx2) +
      dmusig * pmustar/(pwu * sqrt(sigma_sq)) - depsi/ewvsr2) *
      sigx18/(pi * sigx2 * ((Wz)^2 + 1))^2), FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmnsfhalfnormlike_probit <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  mu <- 0
  sigma_sq <- ewu + ewv1
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  ssq <- ((sigma_sq) * sigmastar)
  muvu <- (mu * ewv1 - S * ewu * (epsilon))
  mustar <- muvu/ssq
  dmustar <- dnorm(mustar, 0, 1)
  pmustar <- pnorm(mustar)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  depsi <- dnorm(S * (epsilon)/ewvsr2)
  musq <- (mu * dwu * sqrt(sigma_sq)/ewusr)
  prV <- (1 - ewv1/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx1 <- (dmustar * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmustar)
  sigx2 <- ((1 - pwZ) * depsi/ewvsr2 + dmusig * pmustar * pwZ/(pwu *
    sqrt(sigma_sq)))
  sigx3 <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  sigx4 <- (sigx1 * pwZ/sigx3 + S * (1 - pwZ) * depsi * (epsilon)/ewvsr2^3)
  sigx5 <- (dmustar * dmusig * ewv1/sigmastar - dmusig * muepsi *
    pmustar)
  sigx6 <- (ewusr * (pwu * sqrt(sigma_sq))^2)
  sigx7 <- (sigx5/sigx3 - dmusig * dwu * pmustar * sqrt(sigma_sq)/sigx6)
  sigx8 <- (dmusig * muepsi^2 * pmustar/(sigma_sq)^2)
  sigx9 <- (0.5 * (prU * ewv1/sigmastar) + sigmastar)
  sigx10 <- (sigx9 * muvu/ssq^2 + S * (epsilon)/ssq)
  sigx11 <- (0.5 * sigx8 - sigx10 * dmustar * dmusig)
  sigx12 <- (0.5 * (ewu * pwu/sqrt(sigma_sq)) - 0.5 * musq)
  sigx13 <- (sigx11 * ewu/(pwu * sqrt(sigma_sq)) - sigx12 *
    dmusig * pmustar/(pwu * sqrt(sigma_sq))^2)
  sigx14 <- (mu/ssq - (0.5 * (prV * ewu/sigmastar) + sigmastar) *
    muvu/ssq^2)
  sigx15 <- (0.5 * sigx8 + dmustar * dmusig * sigx14)
  sigx16 <- (dmusig * pmustar * pwu/(pwu * sqrt(sigma_sq))^2)
  sigx17 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi)
  sigx18 <- (dmusig * pmustar/(pwu * sqrt(sigma_sq)) - depsi/ewvsr2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((muepsi^2/(sigma_sq) - 1) * pmustar +
      dmustar * ewu * muepsi/ssq) * dmusig + dmustar *
      (dmusig * muepsi - dmusig * muvu/ewv1) * ewu/ssq) *
      pwZ/sigx3 + (1 - pwZ) * depsi * (S^2 * (epsilon)^2/ewvsr2^2 -
      1)/ewvsr2^3 - sigx4^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((muepsi^2/(sigma_sq) -
      2) * pmustar + dmustar * ewu * muepsi/ssq) * dmusig *
      muepsi/(sigma_sq)^2) - (sigx10 * dmusig * muepsi/(sigma_sq) +
      (sigx9 * ewu/ssq^2 - (sigx10 * muvu/ewv1 + 1/sigmastar)/(sigma_sq)) *
        dmusig) * dmustar) * ewu/(pwu * sqrt(sigma_sq)) -
      (sigx13 * sigx4/sigx2 + sigx12 * sigx1/((sigma_sq) *
        (pwu * sqrt(sigma_sq))^2))) * pwZ/sigx2, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((((dmusig * muepsi - dmusig * muvu/ewv1) * sigx14/(sigma_sq) -
    (0.5 * (prV * ewu/sigmastar) + sigmastar) * dmusig *
      ewu/ssq^2) * dmustar + 0.5 * (((muepsi^2/(sigma_sq) -
    2) * pmustar + dmustar * ewu * muepsi/ssq) * dmusig *
    muepsi/(sigma_sq)^2))/pwu - 0.5 * (sigx1 * pwu/((sigma_sq) *
    (pwu * sqrt(sigma_sq))^2)))/(sigx2 * sqrt(sigma_sq)) -
    (sigx15/pwu - 0.5 * sigx16) * sigx4 * sqrt(sigma_sq)/(sigx2 *
      sqrt(sigma_sq))^2) * ewv1 * pwZ, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (1 - pwZ) * (S * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx2 * ewvsr2^3) -
      sigx4 * sigx17 * ewvsr2/(sigx2 * ewvsr2)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1/sigx3 - (sigx4 *
      sigx18/sigx2 + S * depsi * (epsilon)/ewvsr2^3)) *
      dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) * pmustar/(sigma_sq) -
      sigx10 * dmustar) * ewu) + 0.5 * pmustar) * dmusig *
      muepsi^2/(sigma_sq)^2 - ((sigx10^2 * ewu * muvu/ssq +
      ((0.5 * (ewu/(sigma_sq)) - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv1 * muvu/sigmastar - sigx9 * (2 * (sigx9 *
        (sigma_sq) * muvu * sigmastar/ssq^2) + 2 * (S *
        (epsilon))) * ewu)/ssq^2) * dmusig + sigx10 *
      (0.5 * (dmusig * ewu * muepsi^2/(sigma_sq)^2) + dmusig)) *
      dmustar)/(pwu * sqrt(sigma_sq)) - sigx11 * sigx12/(pwu *
      sqrt(sigma_sq))^2) * ewu - ((((0.5 * (((1 - 0.5 *
      (ewu/(sigma_sq))) * pwu - 0.5 * (mu * dwu/ewusr)) *
      ewu/sqrt(sigma_sq)) - 0.5 * (mu * ((0.5 * (mu^2/ewusr^2) -
      0.5) * sqrt(sigma_sq) + 0.5 * (ewu/sqrt(sigma_sq))) *
      dwu/ewusr)) * dmusig + 0.5 * (sigx12 * dmusig * ewu *
      muepsi^2/(sigma_sq)^2)) * pmustar - (sigx10 * dmustar *
      ewu + 2 * (sigx12 * pmustar * pwu * sqrt(sigma_sq)/(pwu *
      sqrt(sigma_sq))^2)) * sigx12 * dmusig)/(pwu * sqrt(sigma_sq))^2 +
      sigx13^2 * pwZ/sigx2)) * pwZ/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx10 * dmusig * muvu/sigmastar +
      0.5 * (dmusig * muepsi^2/(sigma_sq))) * sigx14/(sigma_sq) -
      ((0.5 * (prU * ewv1/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) -
        1) * ewv1/(sigma_sq) + 1 - 0.5 * (prU * prV))) *
        muvu/sigmastar + mu * sigx9 - (0.5 * (prV * ewu/sigmastar) +
        sigmastar) * (2 * (sigx9 * (sigma_sq) * muvu *
        sigmastar/ssq^2) + S * (epsilon))) * dmusig/ssq^2) *
      dmustar + 0.5 * (((0.5 * (muepsi^2/(sigma_sq)) -
      2) * pmustar/(sigma_sq) - sigx10 * dmustar) * dmusig *
      muepsi^2/(sigma_sq)^2)) * ewu + 0.5 * (mu * sigx15 *
      dwu/(ewusr * pwu)))/pwu - 0.5 * ((sigx11 * ewu *
      pwu - (0.5 * (mu * dwu/ewusr) + 2 * (sigx12 * pwu^2 *
      sqrt(sigma_sq)/(pwu * sqrt(sigma_sq))^2)) * dmusig *
      pmustar)/(pwu * sqrt(sigma_sq))^2))/(sigx2 * sqrt(sigma_sq)) -
      (sigx13 * pwZ * sqrt(sigma_sq) + 0.5 * (sigx2 * ewu/sqrt(sigma_sq))) *
        (sigx15/pwu - 0.5 * sigx16)/(sigx2 * sqrt(sigma_sq))^2) *
      ewv1 * pwZ, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx13 * sigx17 * (1 -
      pwZ) * ewvsr2 * pwZ/(sigx2 * ewvsr2)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx13 * (1 - sigx18 * pwZ/sigx2) *
      dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (dmusig * muepsi^2/(sigma_sq)) -
      dmusig * muvu * sigx14/sigmastar) * ewv1 * sigx14/(sigma_sq) +
      dmusig * (mu/ssq - (((3 * (mu) - 2 * ((0.5 * (prV *
        ewu/sigmastar) + sigmastar) * (sigma_sq) * muvu *
        sigmastar/ssq^2)) * ewv1 - S * ewu * (epsilon)) *
        (0.5 * (prV * ewu/sigmastar) + sigmastar) + (0.5 *
        (ewv1/(sigma_sq)) - 0.5 * (0.5 * prV + ewv1/(sigma_sq))) *
        prV * ewu * muvu/sigmastar)/ssq^2)) * dmustar +
      (0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) * pmustar/(sigma_sq) +
        dmustar * sigx14) * ewv1) + 0.5 * pmustar) *
        dmusig * muepsi^2/(sigma_sq)^2)/pwu - (0.5 *
      (((dmustar * sigx14 - pmustar * pwu^2/(pwu * sqrt(sigma_sq))^2) *
        dmusig + 0.5 * sigx8) * ewv1) + 0.5 * (dmusig *
      pmustar)) * pwu/(pwu * sqrt(sigma_sq))^2)/(sigx2 *
      sqrt(sigma_sq)) - ((sigx15/pwu - 0.5 * sigx16) *
      pwZ + 0.5 * (sigx2/sqrt(sigma_sq))) * (sigx15/pwu -
      0.5 * sigx16) * ewv1/(sigx2 * sqrt(sigma_sq))^2) *
      ewv1 * pwZ, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx15/pwu - 0.5 * sigx16) * sigx17 * (1 - pwZ) * ewv1 *
      ewvsr2 * pwZ/((sigx2 * ewvsr2)^2 * sqrt(sigma_sq))),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx15/pwu - 0.5 * sigx16) *
      (1 - sigx18 * pwZ/sigx2) * dwZ * ewv1/(sigx2 * sqrt(sigma_sq)),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (1 - pwZ) * (S^2 * (0.5 * (0.5 * (S^2 *
      (epsilon)^2/ewvsr2^2) - 1) - 0.25) * depsi * (epsilon)^2/(sigx2 *
      ewvsr2^3) - (sigx17 * (1 - pwZ) + 0.5 * (sigx2 *
      ewvsr2)) * sigx17/(sigx2 * ewvsr2)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((1 - pwZ) * sigx18/sigx2 +
      1) * sigx17 * dwZ/(sigx2 * ewvsr2)), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx18 * dwZ/sigx2 + Wz) *
      sigx18 * dwZ/sigx2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmnsfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar,
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
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq <- ewu + ewv1
  mu <- 0
  sigmastar <- sqrt(ewu * ewv1/(sigma_sq))
  ssq <- ((sigma_sq) * sigmastar)
  muvu <- (mu * ewv1 - S * ewu * (epsilon))
  mustar <- muvu/ssq
  dmustar <- dnorm(mustar, 0, 1)
  pmustar <- pnorm(mustar)
  pwu <- pnorm(mu/ewusr)
  dwu <- dnorm(mu/ewusr, 0, 1)
  muepsi <- (mu + S * (epsilon))
  dmusig <- dnorm(muepsi/sqrt(sigma_sq))
  depsi <- dnorm(S * (epsilon)/ewvsr2)
  musq <- (mu * dwu * sqrt(sigma_sq)/ewusr)
  prV <- (1 - ewv1/(sigma_sq))
  prU <- (1 - ewu/(sigma_sq))
  sigx1 <- (dmustar * dmusig * ewu/sigmastar + dmusig * muepsi *
    pmustar)
  sigx2 <- ((1 - prZ) * dmusig * pmustar/(pwu * sqrt(sigma_sq)) +
    depsi * prZ/ewvsr2)
  sigx3 <- ((sigma_sq) * pwu * sqrt(sigma_sq))
  sigx4 <- ((1 - prZ) * sigx1/sigx3 + S * depsi * prZ * (epsilon)/ewvsr2^3)
  sigx5 <- (dmustar * dmusig * ewv1/sigmastar - dmusig * muepsi *
    pmustar)
  sigx6 <- (ewusr * (pwu * sqrt(sigma_sq))^2)
  sigx7 <- (sigx5/sigx3 - dmusig * dwu * pmustar * sqrt(sigma_sq)/sigx6)
  sigx8 <- (dmusig * muepsi^2 * pmustar/(sigma_sq)^2)
  sigx9 <- (0.5 * (prU * ewv1/sigmastar) + sigmastar)
  sigx10 <- (sigx9 * muvu/ssq^2 + S * (epsilon)/ssq)
  sigx11 <- (0.5 * sigx8 - sigx10 * dmustar * dmusig)
  sigx12 <- (0.5 * (ewu * pwu/sqrt(sigma_sq)) - 0.5 * musq)
  sigx13 <- (sigx11 * ewu/(pwu * sqrt(sigma_sq)) - sigx12 *
    dmusig * pmustar/(pwu * sqrt(sigma_sq))^2)
  sigx14 <- (mu/ssq - (0.5 * (prV * ewu/sigmastar) + sigmastar) *
    muvu/ssq^2)
  sigx15 <- (0.5 * sigx8 + dmustar * dmusig * sigx14)
  sigx16 <- (dmusig * pmustar * pwu/(pwu * sqrt(sigma_sq))^2)
  sigx17 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr2^2) - 0.5 *
    depsi)
  sigx18 <- (dmusig * pmustar/(pwu * sqrt(sigma_sq)) - depsi/ewvsr2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((muepsi^2/(sigma_sq) - 1) * pmustar +
      dmustar * ewu * muepsi/ssq) * dmusig + dmustar *
      (dmusig * muepsi - dmusig * muvu/ewv1) * ewu/ssq) *
      (1 - prZ)/sigx3 + depsi * prZ * (S^2 * (epsilon)^2/ewvsr2^2 -
      1)/ewvsr2^3 - sigx4^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((muepsi^2/(sigma_sq) -
      2) * pmustar + dmustar * ewu * muepsi/ssq) * dmusig *
      muepsi/(sigma_sq)^2) - (sigx10 * dmusig * muepsi/(sigma_sq) +
      (sigx9 * ewu/ssq^2 - (sigx10 * muvu/ewv1 + 1/sigmastar)/(sigma_sq)) *
        dmusig) * dmustar) * ewu/(pwu * sqrt(sigma_sq)) -
      (sigx13 * sigx4/sigx2 + sigx12 * sigx1/((sigma_sq) *
        (pwu * sqrt(sigma_sq))^2))) * (1 - prZ)/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((((dmusig * muepsi - dmusig * muvu/ewv1) * sigx14/(sigma_sq) -
    (0.5 * (prV * ewu/sigmastar) + sigmastar) * dmusig *
      ewu/ssq^2) * dmustar + 0.5 * (((muepsi^2/(sigma_sq) -
    2) * pmustar + dmustar * ewu * muepsi/ssq) * dmusig *
    muepsi/(sigma_sq)^2))/pwu - 0.5 * (sigx1 * pwu/((sigma_sq) *
    (pwu * sqrt(sigma_sq))^2)))/(sigx2 * sqrt(sigma_sq)) -
    (sigx15/pwu - 0.5 * sigx16) * sigx4 * sqrt(sigma_sq)/(sigx2 *
      sqrt(sigma_sq))^2) * (1 - prZ) * ewv1, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * prZ * (S * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2 -
      2) - 0.5) * depsi * (epsilon)/(sigx2 * ewvsr2^3) -
      sigx4 * sigx17 * ewvsr2/(sigx2 * ewvsr2)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1/sigx3 - (sigx4 *
      sigx18/sigx2 + S * depsi * (epsilon)/ewvsr2^3)) *
      prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) * pmustar/(sigma_sq) -
      sigx10 * dmustar) * ewu) + 0.5 * pmustar) * dmusig *
      muepsi^2/(sigma_sq)^2 - ((sigx10^2 * ewu * muvu/ssq +
      ((0.5 * (ewu/(sigma_sq)) - 0.5 * (0.5 * prU + ewu/(sigma_sq))) *
        prU * ewv1 * muvu/sigmastar - sigx9 * (2 * (sigx9 *
        (sigma_sq) * muvu * sigmastar/ssq^2) + 2 * (S *
        (epsilon))) * ewu)/ssq^2) * dmusig + sigx10 *
      (0.5 * (dmusig * ewu * muepsi^2/(sigma_sq)^2) + dmusig)) *
      dmustar)/(pwu * sqrt(sigma_sq)) - sigx11 * sigx12/(pwu *
      sqrt(sigma_sq))^2) * ewu - ((((0.5 * (((1 - 0.5 *
      (ewu/(sigma_sq))) * pwu - 0.5 * (mu * dwu/ewusr)) *
      ewu/sqrt(sigma_sq)) - 0.5 * (mu * ((0.5 * (mu^2/ewusr^2) -
      0.5) * sqrt(sigma_sq) + 0.5 * (ewu/sqrt(sigma_sq))) *
      dwu/ewusr)) * dmusig + 0.5 * (sigx12 * dmusig * ewu *
      muepsi^2/(sigma_sq)^2)) * pmustar - (sigx10 * dmustar *
      ewu + 2 * (sigx12 * pmustar * pwu * sqrt(sigma_sq)/(pwu *
      sqrt(sigma_sq))^2)) * sigx12 * dmusig)/(pwu * sqrt(sigma_sq))^2 +
      sigx13^2 * (1 - prZ)/sigx2)) * (1 - prZ)/sigx2, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx10 * dmusig * muvu/sigmastar +
      0.5 * (dmusig * muepsi^2/(sigma_sq))) * sigx14/(sigma_sq) -
      ((0.5 * (prU * ewv1/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) -
        1) * ewv1/(sigma_sq) + 1 - 0.5 * (prU * prV))) *
        muvu/sigmastar + mu * sigx9 - (0.5 * (prV * ewu/sigmastar) +
        sigmastar) * (2 * (sigx9 * (sigma_sq) * muvu *
        sigmastar/ssq^2) + S * (epsilon))) * dmusig/ssq^2) *
      dmustar + 0.5 * (((0.5 * (muepsi^2/(sigma_sq)) -
      2) * pmustar/(sigma_sq) - sigx10 * dmustar) * dmusig *
      muepsi^2/(sigma_sq)^2)) * ewu + 0.5 * (mu * sigx15 *
      dwu/(ewusr * pwu)))/pwu - 0.5 * ((sigx11 * ewu *
      pwu - (0.5 * (mu * dwu/ewusr) + 2 * (sigx12 * pwu^2 *
      sqrt(sigma_sq)/(pwu * sqrt(sigma_sq))^2)) * dmusig *
      pmustar)/(pwu * sqrt(sigma_sq))^2))/(sigx2 * sqrt(sigma_sq)) -
      (sigx13 * (1 - prZ) * sqrt(sigma_sq) + 0.5 * (sigx2 *
        ewu/sqrt(sigma_sq))) * (sigx15/pwu - 0.5 * sigx16)/(sigx2 *
        sqrt(sigma_sq))^2) * (1 - prZ) * ewv1, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx13 * sigx17 * (1 -
      prZ) * prZ * ewvsr2/(sigx2 * ewvsr2)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx13 * (1 - (1 - prZ) *
      sigx18/sigx2) * prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (dmusig * muepsi^2/(sigma_sq)) -
      dmusig * muvu * sigx14/sigmastar) * ewv1 * sigx14/(sigma_sq) +
      dmusig * (mu/ssq - (((3 * (mu) - 2 * ((0.5 * (prV *
        ewu/sigmastar) + sigmastar) * (sigma_sq) * muvu *
        sigmastar/ssq^2)) * ewv1 - S * ewu * (epsilon)) *
        (0.5 * (prV * ewu/sigmastar) + sigmastar) + (0.5 *
        (ewv1/(sigma_sq)) - 0.5 * (0.5 * prV + ewv1/(sigma_sq))) *
        prV * ewu * muvu/sigmastar)/ssq^2)) * dmustar +
      (0.5 * (((0.5 * (muepsi^2/(sigma_sq)) - 2) * pmustar/(sigma_sq) +
        dmustar * sigx14) * ewv1) + 0.5 * pmustar) *
        dmusig * muepsi^2/(sigma_sq)^2)/pwu - (0.5 *
      (((dmustar * sigx14 - pmustar * pwu^2/(pwu * sqrt(sigma_sq))^2) *
        dmusig + 0.5 * sigx8) * ewv1) + 0.5 * (dmusig *
      pmustar)) * pwu/(pwu * sqrt(sigma_sq))^2)/(sigx2 *
      sqrt(sigma_sq)) - ((sigx15/pwu - 0.5 * sigx16) *
      (1 - prZ) + 0.5 * (sigx2/sqrt(sigma_sq))) * (sigx15/pwu -
      0.5 * sigx16) * ewv1/(sigx2 * sqrt(sigma_sq))^2) *
      (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx15/pwu - 0.5 * sigx16) * sigx17 * (1 - prZ) * prZ *
      ewv1 * ewvsr2/((sigx2 * ewvsr2)^2 * sqrt(sigma_sq))),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx15/pwu - 0.5 * sigx16) *
      (1 - (1 - prZ) * sigx18/sigx2) * prZ * ewv1 * ewz/(sigx2 *
      sqrt(sigma_sq)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prZ * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr2^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx2 * ewvsr2^3) -
      (sigx17 * prZ + 0.5 * (sigx2 * ewvsr2)) * sigx17/(sigx2 *
        ewvsr2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx18 * prZ/sigx2 + 1) *
      sigx17 * prZ * ewz/(sigx2 * ewvsr2)), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (sigx18 * prZ/sigx2 +
      1) * ewz) * sigx18 * prZ * ewz/sigx2, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf halfnormal-normal distribution
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
zisfhalfnormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfhalfnormlike_logit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(czisfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfhalfnormlike_logit,
      grad = cgradzisfhalfnormlike_logit, hess = chesszisfhalfnormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfhalfnormlike_logit(mleObj$par,
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
      mleObj$hessian <- chesszisfhalfnormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfhalfnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfhalfnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfhalfnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

## cauchit specification class membership
zisfhalfnormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfhalfnormlike_cauchit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(czisfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfhalfnormlike_cauchit,
      grad = cgradzisfhalfnormlike_cauchit, hess = chesszisfhalfnormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfhalfnormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- chesszisfhalfnormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfhalfnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfhalfnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfhalfnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

## probit specification class membership
zisfhalfnormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfhalfnormlike_probit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(czisfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfhalfnormlike_probit,
      grad = cgradzisfhalfnormlike_probit, hess = chesszisfhalfnormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfhalfnormlike_probit(mleObj$par,
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
      mleObj$hessian <- chesszisfhalfnormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfhalfnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfhalfnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfhalfnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

## cloglog specification class membership
zisfhalfnormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfhalfnormlike_cloglog(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(czisfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfhalfnormlike_cloglog,
      grad = cgradzisfhalfnormlike_cloglog, hess = chesszisfhalfnormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfhalfnormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- chesszisfhalfnormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfhalfnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfhalfnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfhalfnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

# Different sigma_v

## logit specification class membership
mnsfhalfnormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfhalfnormlike_logit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmnsfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfhalfnormlike_logit,
      grad = cgradmnsfhalfnormlike_logit, hess = chessmnsfhalfnormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfhalfnormlike_logit(mleObj$par,
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
      mleObj$hessian <- chessmnsfhalfnormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfhalfnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfhalfnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfhalfnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

## cauchit specification class membership
mnsfhalfnormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfhalfnormlike_cauchit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmnsfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfhalfnormlike_cauchit,
      grad = cgradmnsfhalfnormlike_cauchit, hess = chessmnsfhalfnormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfhalfnormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- chessmnsfhalfnormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfhalfnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfhalfnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfhalfnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

## probit specification class membership
mnsfhalfnormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfhalfnormlike_probit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmnsfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfhalfnormlike_probit,
      grad = cgradmnsfhalfnormlike_probit, hess = chessmnsfhalfnormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfhalfnormlike_probit(mleObj$par,
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
      mleObj$hessian <- chessmnsfhalfnormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfhalfnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfhalfnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfhalfnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

## cloglog specification class membership
mnsfhalfnormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfhalfnormlike_cloglog(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmnsfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfhalfnormlike_cloglog,
      grad = cgradmnsfhalfnormlike_cloglog, hess = chessmnsfhalfnormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfhalfnormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- chessmnsfhalfnormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfhalfnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfhalfnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfhalfnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

# Conditional efficiencies estimation ----------
#' efficiencies for zisf halfnormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v

## logit specification class membership
czisfhalfnormeff_logit <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cauchit specification class membership
czisfhalfnormeff_cauchit <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## probit specification class membership
czisfhalfnormeff_probit <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cloglog specification class membership
czisfhalfnormeff_cloglog <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

# Different sigma_v

## logit specification class membership
cmnsfhalfnormeff_logit <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cauchit specification class membership
cmnsfhalfnormeff_cauchit <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## probit specification class membership
cmnsfhalfnormeff_probit <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cloglog specification class membership
cmnsfhalfnormeff_cloglog <- function(object, level) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for zisf halfnormal-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v

## logit specification class membership
czisfmarghalfnorm_Eu_logit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarghalfnorm_Vu_logit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
czisfmarghalfnorm_Eu_cauchit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarghalfnorm_Vu_cauchit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## probit specification class membership
czisfmarghalfnorm_Eu_probit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarghalfnorm_Vu_probit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
czisfmarghalfnorm_Eu_cloglog <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarghalfnorm_Vu_cloglog <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
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
cmnsfmarghalfnorm_Eu_logit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarghalfnorm_Vu_logit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmnsfmarghalfnorm_Eu_cauchit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarghalfnorm_Vu_cauchit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmnsfmarghalfnorm_Eu_probit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarghalfnorm_Vu_probit <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmnsfmarghalfnorm_Eu_cloglog <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarghalfnorm_Vu_cloglog <- function(object) {
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}
