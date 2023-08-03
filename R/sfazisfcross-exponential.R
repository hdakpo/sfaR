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
# Convolution: exponential - normal                                            #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf exponential-normal distribution
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
czisfexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) * pnorm(-S *
    epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cauchit specification class membership
czisfexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) * pnorm(-S *
    epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## probit specification class membership
czisfexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) * pnorm(-S *
    epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

## cloglog specification class membership
czisfexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) * pnorm(-S *
    epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

# Different sigma_v

## logit specification class membership
cmnsfexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) * pnorm(-S *
    epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
cmnsfexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) * pnorm(-S *
    epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
cmnsfexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) * pnorm(-S *
    epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
cmnsfexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) * pnorm(-S *
    epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf exponential-normal distribution
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
cstzisfexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- cstexponorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initExpo <- NULL
  } else {
    cat("Initialization: SFA + exponential - normal distributions...\n")
    initExpo <- maxLik::maxLik(logLik = cexponormlike, start = cstexponorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradexponormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initExpo$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("ZISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initExpo = initExpo))
}

# Different sigma_v
cstmnsfexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- cstexponorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initExpo <- NULL
  } else {
    cat("Initialization: SFA + exponential - normal distributions...\n")
    initExpo <- maxLik::maxLik(logLik = cexponormlike, start = cstexponorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradexponormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initExpo$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar +
    2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("ZISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initExpo = initExpo))
}

# Gradient of the likelihood function ----------
#' gradient for zisf exponential-normal distribution
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
cgradzisfexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewusr <- exp(Wu/2)
  ewz <- exp(Wz)
  mustar <- -(ewvsr/ewusr + S * (epsilon)/ewvsr)
  pmustar <- pnorm(mustar)
  dmustar <- dnorm(mustar)
  epsi <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  wzdeno <- (1 - ewz/(1 + ewz))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  sigx1 <- (wzdeno * dwsr/ewvsr + epsi * ewz * pmustar/((1 + ewz) * ewusr))
  dpmu <- (dmustar/ewvsr - pmustar/ewusr)
  sigx2 <- (dpmu * epsi * ewz/((1 + ewz) * ewusr) + S * wzdeno * dwsr * (epsilon)/ewvsr^3)
  deemu <- (dmustar * ewvsr/ewusr)
  sigx3 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx4 <- ((1 + ewz) * ewusr * pmustar/((1 + ewz) * ewusr)^2)
  sigx5 <- ((0.5 * deemu - sigx3 * pmustar)/((1 + ewz) * ewusr) - 0.5 * sigx4)
  sigx6 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx7 <- (ewv * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr)
  sigx9 <- (sigx8 * wzdeno/ewvsr + epsi * sigx7 * ewz/((1 + ewz) * ewusr))
  sigx10 <- (1/((1 + ewz) * ewusr) - ewusr * ewz/((1 + ewz) * ewusr)^2)
  sigx11 <- (sigx10 * epsi * pmustar - wzdeno * dwsr/((1 + ewz) * ewvsr))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx1, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * epsi * ewz/sigx1, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx9/sigx1, FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx11 * ewz/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradzisfexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  mustar <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  sigx1 <- (dmustar/ewvsr - pmustar/ewusr)
  sigx2 <- (ewz1 * sigx1 * epsi/ewusr + S * ewz2 * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (ewz2 * dwsr/ewvsr + ewz1 * epsi * pmustar/ewusr)
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar * ewvsr/ewusr) - sigx4 * pmustar)
  sigx6 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr)
  sigx7 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx8 <- (ewz2 * sigx6/ewvsr + ewz1 * epsi * (ewv * pmustar/(2 * ewu) - sigx7 *
    dmustar)/ewusr)
  sigx9 <- (epsi * pmustar/ewusr - dwsr/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * ewz1 * epsi/(sigx3 * ewusr), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx8/sigx3, FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx9/(pi * sigx3 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradzisfexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  mustar <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  sigx1 <- (dmustar/ewvsr - pmustar/ewusr)
  sigx2 <- (sigx1 * epsi * pwZ/ewusr + S * (1 - pwZ) * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - pwZ) * dwsr/ewvsr + epsi * pmustar * pwZ/ewusr)
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar * ewvsr/ewusr) - sigx4 * pmustar)
  sigx6 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr)
  sigx7 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx8 <- (sigx6 * (1 - pwZ)/ewvsr + epsi * (ewv * pmustar/(2 * ewu) - sigx7 *
    dmustar) * pwZ/ewusr)
  sigx9 <- (epsi * pmustar/ewusr - dwsr/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * epsi * pwZ/(sigx3 * ewusr), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx8/sigx3, FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = dwZ * sigx9/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradzisfexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  mustar <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  sigx1 <- (dmustar/ewvsr - pmustar/ewusr)
  sigx2 <- ((1 - prZ) * sigx1 * epsi/ewusr + S * dwsr * prZ * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - prZ) * epsi * pmustar/ewusr + dwsr * prZ/ewvsr)
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar * ewvsr/ewusr) - sigx4 * pmustar)
  sigx6 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr)
  sigx7 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx8 <- (sigx6 * prZ/ewvsr + (1 - prZ) * epsi * (ewv * pmustar/(2 * ewu) - sigx7 *
    dmustar)/ewusr)
  sigx9 <- (epsi * pmustar/ewusr - dwsr/ewvsr)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * (1 - prZ) * epsi/(sigx3 * ewusr),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx8/sigx3, FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = prZ * sigx9 * ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v

## logit specification class membership
cgradmnsfexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  mustar <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  dwsr <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  dpmu <- (dmustar/ewv1_h - pmustar/ewu_h)
  deemu <- (dmustar * ewv1_h/ewu_h)
  sigx1 <- (dpmu * epsi * ewz/(wzdeno * ewu_h) + S * prC * dwsr * (epsilon)/ewv2_h^3)
  sigx2 <- (prC * dwsr/ewv2_h + epsi * ewz * pmustar/(wzdeno * ewu_h))
  sigx3 <- (0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx4 <- (wzdeno * ewu_h * pmustar/(wzdeno * ewu_h)^2)
  sigx5 <- ((0.5 * deemu - sigx3 * pmustar)/(wzdeno * ewu_h) - 0.5 * sigx4)
  sigx6 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx7 <- (ewv1 * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewv2_h^2) - 0.5 * dwsr)
  sigx9 <- (1/(wzdeno * ewu_h) - ewu_h * ewz/(wzdeno * ewu_h)^2)
  sigx10 <- (sigx9 * epsi * pmustar - prC * dwsr/(wzdeno * ewv2_h))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * epsi * ewz/sigx2, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = epsi * sigx7 * ewz/(sigx2 * wzdeno * ewu_h), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx8 * prC/(sigx2 * ewv2_h), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx10 * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmnsfexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  mustar <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  dwsr <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  dpmu <- (dmustar/ewv1_h - pmustar/ewu_h)
  deemu <- (dmustar * ewv1_h/ewu_h)
  sigx1 <- (ewz1 * dpmu * epsi/ewu_h + S * ewz2 * dwsr * (epsilon)/ewv2_h^3)
  sigx2 <- (ewz2 * dwsr/ewv2_h + ewz1 * epsi * pmustar/ewu_h)
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * deemu - sigx4 * pmustar)
  sigx6 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx7 <- (ewv1 * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewv2_h^2) - 0.5 * dwsr)
  sigx9 <- (epsi * pmustar/ewu_h - dwsr/ewv2_h)
  sigx10 <- (pi * sigx2 * ((Wz)^2 + 1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * ewz1 * epsi/(sigx2 * ewu_h), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = ewz1 * epsi * sigx7/(sigx2 * ewu_h), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = ewz2 * sigx8/(sigx2 * ewv2_h), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx9/sigx10, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmnsfexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  mustar <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  dwsr <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  dpmu <- (dmustar/ewv1_h - pmustar/ewu_h)
  deemu <- (dmustar * ewv1_h/ewu_h)
  sigx1 <- (dpmu * epsi * pwZ/ewu_h + S * (1 - pwZ) * dwsr * (epsilon)/ewv2_h^3)
  sigx2 <- ((1 - pwZ) * dwsr/ewv2_h + epsi * pmustar * pwZ/ewu_h)
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * deemu - sigx4 * pmustar)
  sigx6 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx7 <- (ewv1 * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewv2_h^2) - 0.5 * dwsr)
  sigx9 <- (epsi * pmustar/ewu_h - dwsr/ewv2_h)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * epsi * pwZ/(sigx2 * ewu_h), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = epsi * sigx7 * pwZ/(sigx2 * ewu_h), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx8 * (1 - pwZ)/(sigx2 * ewv2_h), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = dwZ * sigx9/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmnsfexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  mustar <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  dwsr <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  dpmu <- (dmustar/ewv1_h - pmustar/ewu_h)
  deemu <- (dmustar * ewv1_h/ewu_h)
  sigx1 <- ((1 - prZ) * dpmu * epsi/ewu_h + S * dwsr * prZ * (epsilon)/ewv2_h^3)
  sigx2 <- ((1 - prZ) * epsi * pmustar/ewu_h + dwsr * prZ/ewv2_h)
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * deemu - sigx4 * pmustar)
  sigx6 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx7 <- (ewv1 * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewv2_h^2) - 0.5 * dwsr)
  sigx9 <- (epsi * pmustar/ewu_h - dwsr/ewv2_h)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * (1 - prZ) * epsi/(sigx2 * ewu_h),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (1 - prZ) * epsi * sigx7/(sigx2 *
      ewu_h), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx8 * prZ/(sigx2 *
      ewv2_h), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = prZ * sigx9 * ewz/sigx2,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf exponential-normal distribution
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
chesszisfexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewusr <- exp(Wu/2)
  ewz <- exp(Wz)
  mustar <- -(ewvsr/ewusr + S * (epsilon)/ewvsr)
  pmustar <- pnorm(mustar)
  dmustar <- dnorm(mustar)
  epsi <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  wzdeno <- (1 - ewz/(1 + ewz))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  sigx1 <- (wzdeno * dwsr/ewvsr + epsi * ewz * pmustar/((1 + ewz) * ewusr))
  dpmu <- (dmustar/ewvsr - pmustar/ewusr)
  sigx2 <- (dpmu * epsi * ewz/((1 + ewz) * ewusr) + S * wzdeno * dwsr * (epsilon)/ewvsr^3)
  deemu <- (dmustar * ewvsr/ewusr)
  sigx3 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx4 <- ((1 + ewz) * ewusr * pmustar/((1 + ewz) * ewusr)^2)
  sigx5 <- ((0.5 * deemu - sigx3 * pmustar)/((1 + ewz) * ewusr) - 0.5 * sigx4)
  sigx6 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx7 <- (ewv * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr)
  sigx9 <- (sigx8 * wzdeno/ewvsr + epsi * sigx7 * ewz/((1 + ewz) * ewusr))
  sigx10 <- (1/((1 + ewz) * ewusr) - ewusr * ewz/((1 + ewz) * ewusr)^2)
  sigx11 <- (sigx10 * epsi * pmustar - wzdeno * dwsr/((1 + ewz) * ewvsr))
  sigx12 <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (((sigx12/ewvsr - 1/ewusr) * dmustar/ewvsr - dpmu/ewusr) * epsi * ewz/((1 +
    ewz) * ewusr) + wzdeno * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 -
    sigx2^2/sigx1)/sigx1, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 *
      ewu)^2)) * pmustar - 0.5 * deemu)/ewusr + (0.5 * (sigx12/ewusr) - sigx3/ewvsr) *
      dmustar)/(1 + ewz) + 0.5 * sigx4)/ewusr - (sigx5 * sigx2/sigx1 + 0.5 *
      ((1 + ewz) * dmustar * ewusr/(((1 + ewz) * ewusr)^2 * ewvsr)))) * epsi *
      ewz/sigx1, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmustar * (ewv/(2 * ewu) - (sigx6 * sigx12 +
      0.5))/ewvsr - sigx7/ewusr) * epsi * ewz/((1 + ewz) * ewusr) + S * (0.5 *
      (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * wzdeno * dwsr * (epsilon)/ewvsr^3 -
      sigx9 * sigx2/sigx1)/sigx1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx10 *
    dpmu * epsi - (sigx11 * sigx2/sigx1 + S * wzdeno * dwsr * (epsilon)/((1 +
    ewz) * ewvsr^3))) * ewz/sigx1, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr * sigx12/ewusr) - 0.5) -
      0.5 * sigx3) * dmustar * ewvsr/ewusr - (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) *
      ewu * ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pmustar)/((1 +
      ewz) * ewusr) - ((sigx5 * epsi * ewz/sigx1 + 0.5 * (S * (epsilon)/ewusr) +
      2 * (ewu * ewv/(2 * ewu)^2)) * sigx5 + (0.5 * ((0.5 - (1 + ewz)^2 * ewusr^2/((1 +
      ewz) * ewusr)^2) * ewusr * pmustar + 0.5 * (dmustar * ewvsr)) + 0.5 *
      ((0.5 * deemu - sigx3 * pmustar) * ewusr)) * (1 + ewz)/((1 + ewz) * ewusr)^2)) *
      epsi * ewz/sigx1, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((dmustar *
    ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pmustar/(2 * ewu)^2)) * ewv - ((0.5 *
    (sigx6 * sigx12) - 0.25) * dmustar * ewvsr/ewusr + sigx3 * sigx7))/((1 +
    ewz) * ewusr) - (sigx5 * sigx9/sigx1 + 0.5 * ((1 + ewz) * ewusr * sigx7/((1 +
    ewz) * ewusr)^2))) * epsi * ewz/sigx1, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (0.5 * (sigx10 * dmustar * ewvsr/ewusr) - ((((0.5 - (1 + ewz)^2 * ewusr^2/((1 +
      ewz) * ewusr)^2) * ewz + 0.5 * (1 + ewz)) * ewusr/((1 + ewz) * ewusr)^2 +
      sigx3 * sigx10) * pmustar + sigx5 * sigx11 * ewz/sigx1)) * epsi * ewz/sigx1,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((sigx7/2 + (pmustar - sigx6 * dmustar)/2) * ewv/ewu - (0.25 * (ewvsr/ewusr) +
      0.25 * (S * (epsilon)/ewvsr) - sigx6^2 * sigx12) * dmustar) * epsi *
      ewz/((1 + ewz) * ewusr) + wzdeno * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
      1) - 0.25) * dwsr * (epsilon)^2/ewvsr^2 - 0.5 * sigx8)/ewvsr - sigx9^2/sigx1)/sigx1,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx10 * epsi * sigx7 - (sigx9 * sigx11/sigx1 +
      (0.5 * (S^2 * dwsr * (epsilon)^2/((1 + ewz) * ewvsr^3)) - 0.5 * ((1 +
        ewz) * dwsr * ewvsr/((1 + ewz) * ewvsr)^2)) * wzdeno)) * ewz/sigx1,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((wzdeno * (1/((1 + ewz)^2 * ewvsr) + ewvsr/((1 +
      ewz) * ewvsr)^2) * dwsr - (sigx11^2/sigx1 + (2 - 2 * ((1 + ewz) * ewusr^2 *
      ewz/((1 + ewz) * ewusr)^2)) * ewusr * epsi * pmustar/((1 + ewz) * ewusr)^2)) *
      ewz + sigx10 * epsi * pmustar - wzdeno * dwsr/((1 + ewz) * ewvsr)) *
      ewz/sigx1, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesszisfexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  mustar <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  sigx1 <- (dmustar/ewvsr - pmustar/ewusr)
  sigx2 <- (ewz1 * sigx1 * epsi/ewusr + S * ewz2 * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (ewz2 * dwsr/ewvsr + ewz1 * epsi * pmustar/ewusr)
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar * ewvsr/ewusr) - sigx4 * pmustar)
  sigx6 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr)
  sigx7 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx8 <- (ewz2 * sigx6/ewvsr + ewz1 * epsi * (ewv * pmustar/(2 * ewu) - sigx7 *
    dmustar)/ewusr)
  sigx9 <- (epsi * pmustar/ewusr - dwsr/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (((mustar/ewvsr - 1/ewusr) * dmustar/ewvsr - sigx1/ewusr) * ewz1 *
    epsi/ewusr + ewz2 * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((0.5 * (S * (epsilon)/ewusr) + 1 + 2 * (ewu * ewv/(2 *
      ewu)^2)) * pmustar - 0.5 * (dmustar * ewvsr/ewusr))/ewusr + (0.5 * (mustar/ewusr) -
      sigx4/ewvsr) * dmustar)/(sigx3 * ewusr) - sigx2 * sigx5 * ewusr/(sigx3 *
      ewusr)^2) * ewz1 * epsi, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (ewz1 * (dmustar * (ewv/(2 * ewu) - (sigx7 *
      mustar + 0.5))/ewvsr - (ewv * pmustar/(2 * ewu) - sigx7 * dmustar)/ewusr) *
      epsi/ewusr + S * ewz2 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) *
      dwsr * (epsilon)/ewvsr^3 - sigx8 * sigx2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((sigx1 *
    epsi/ewusr - S * dwsr * (epsilon)/ewvsr^3)/(pi * sigx3 * ((Wz)^2 + 1)) -
    pi * sigx2 * ((Wz)^2 + 1) * sigx9/(pi * sigx3 * ((Wz)^2 + 1))^2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr * mustar/ewusr) - 0.5) -
      0.5 * sigx4) * dmustar * ewvsr/ewusr - (sigx5 * (0.5 * (S * (epsilon)/ewusr) +
      2 * (ewu * ewv/(2 * ewu)^2)) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) *
      ewu * ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pmustar))/(sigx3 *
      ewusr) - (sigx5 * ewz1 * epsi + 0.5 * (sigx3 * ewusr)) * sigx5/(sigx3 *
      ewusr)^2) * ewz1 * epsi, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((dmustar *
    ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pmustar/(2 * ewu)^2)) * ewv - (sigx8 *
    sigx5/sigx3 + (0.5 * (sigx7 * mustar) - 0.25) * dmustar * ewvsr/ewusr + sigx4 *
    (ewv * pmustar/(2 * ewu) - sigx7 * dmustar))) * ewz1 * epsi/(sigx3 * ewusr),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx5 * (1/(pi * sigx3 * ((Wz)^2 + 1)) - pi * ((Wz)^2 + 1) * ewz1 * sigx9/(pi *
    sigx3 * ((Wz)^2 + 1))^2) * epsi/ewusr, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((ewv * pmustar/(2 * ewu) - sigx7 * dmustar)/2 + (pmustar - sigx7 * dmustar)/2) *
      ewv/ewu - (0.25 * (ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) - sigx7^2 *
      mustar) * dmustar) * ewz1 * epsi/ewusr + ewz2 * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr * (epsilon)^2/ewvsr^2 -
      0.5 * sigx6)/ewvsr - sigx8^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((epsi * (ewv * pmustar/(2 * ewu) - sigx7 * dmustar)/ewusr -
      sigx6/ewvsr)/(pi * sigx3 * ((Wz)^2 + 1)) - pi * sigx8 * ((Wz)^2 + 1) *
      sigx9/(pi * sigx3 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (pi * Wz * sigx3) + epsi * pmustar/ewusr -
      dwsr/ewvsr) * sigx9/(pi * sigx3 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesszisfexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  mustar <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  sigx1 <- (dmustar/ewvsr - pmustar/ewusr)
  sigx2 <- (sigx1 * epsi * pwZ/ewusr + S * (1 - pwZ) * dwsr * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - pwZ) * dwsr/ewvsr + epsi * pmustar * pwZ/ewusr)
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar * ewvsr/ewusr) - sigx4 * pmustar)
  sigx6 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr)
  sigx7 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx8 <- (sigx6 * (1 - pwZ)/ewvsr + epsi * (ewv * pmustar/(2 * ewu) - sigx7 *
    dmustar) * pwZ/ewusr)
  sigx9 <- (epsi * pmustar/ewusr - dwsr/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (((mustar/ewvsr - 1/ewusr) * dmustar/ewvsr - sigx1/ewusr) * epsi * pwZ/ewusr +
      (1 - pwZ) * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((0.5 * (S * (epsilon)/ewusr) + 1 + 2 * (ewu * ewv/(2 *
      ewu)^2)) * pmustar - 0.5 * (dmustar * ewvsr/ewusr))/ewusr + (0.5 * (mustar/ewusr) -
      sigx4/ewvsr) * dmustar)/(sigx3 * ewusr) - sigx2 * sigx5 * ewusr/(sigx3 *
      ewusr)^2) * epsi * pwZ, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmustar * (ewv/(2 * ewu) - (sigx7 * mustar +
      0.5))/ewvsr - (ewv * pmustar/(2 * ewu) - sigx7 * dmustar)/ewusr) * epsi *
      pwZ/ewusr + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * (1 -
      pwZ) * dwsr * (epsilon)/ewvsr^3 - sigx8 * sigx2/sigx3)/sigx3, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx1 *
    epsi/ewusr - (sigx2 * sigx9/sigx3 + S * dwsr * (epsilon)/ewvsr^3)) * dwZ/sigx3,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr * mustar/ewusr) - 0.5) -
      0.5 * sigx4) * dmustar * ewvsr/ewusr - (sigx5 * (0.5 * (S * (epsilon)/ewusr) +
      2 * (ewu * ewv/(2 * ewu)^2)) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) *
      ewu * ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pmustar))/(sigx3 *
      ewusr) - (sigx5 * epsi * pwZ + 0.5 * (sigx3 * ewusr)) * sigx5/(sigx3 *
      ewusr)^2) * epsi * pwZ, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((dmustar *
    ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pmustar/(2 * ewu)^2)) * ewv - (sigx8 *
    sigx5/sigx3 + (0.5 * (sigx7 * mustar) - 0.25) * dmustar * ewvsr/ewusr + sigx4 *
    (ewv * pmustar/(2 * ewu) - sigx7 * dmustar))) * epsi * pwZ/(sigx3 * ewusr),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx5 * (1 - sigx9 * pwZ/sigx3) * dwZ * epsi/(sigx3 * ewusr), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((ewv * pmustar/(2 * ewu) - sigx7 * dmustar)/2 + (pmustar - sigx7 * dmustar)/2) *
      ewv/ewu - (0.25 * (ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) - sigx7^2 *
      mustar) * dmustar) * epsi * pwZ/ewusr + (1 - pwZ) * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr * (epsilon)^2/ewvsr^2 -
      0.5 * sigx6)/ewvsr - sigx8^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * dwZ * (epsi * (ewv * pmustar/(2 * ewu) - sigx7 *
      dmustar)/ewusr - (sigx8 * sigx9/sigx3 + sigx6/ewvsr))/sigx3, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (dwZ * (dwZ * sigx9/sigx3 + Wz) * sigx9/sigx3),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesszisfexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  mustar <- (ewvsr/ewusr + S * (epsilon)/ewvsr)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  sigx1 <- (dmustar/ewvsr - pmustar/ewusr)
  sigx2 <- ((1 - prZ) * sigx1 * epsi/ewusr + S * dwsr * prZ * (epsilon)/ewvsr^3)
  sigx3 <- ((1 - prZ) * epsi * pmustar/ewusr + dwsr * prZ/ewvsr)
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmustar * ewvsr/ewusr) - sigx4 * pmustar)
  sigx6 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr)
  sigx7 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx8 <- (sigx6 * prZ/ewvsr + (1 - prZ) * epsi * (ewv * pmustar/(2 * ewu) - sigx7 *
    dmustar)/ewusr)
  sigx9 <- (epsi * pmustar/ewusr - dwsr/ewvsr)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar, ncol = nXvar + nuZUvar +
    nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (((mustar/ewvsr - 1/ewusr) * dmustar/ewvsr - sigx1/ewusr) * (1 - prZ) * epsi/ewusr +
      dwsr * prZ * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((0.5 * (S * (epsilon)/ewusr) + 1 + 2 * (ewu * ewv/(2 *
      ewu)^2)) * pmustar - 0.5 * (dmustar * ewvsr/ewusr))/ewusr + (0.5 * (mustar/ewusr) -
      sigx4/ewvsr) * dmustar)/(sigx3 * ewusr) - sigx2 * sigx5 * ewusr/(sigx3 *
      ewusr)^2) * (1 - prZ) * epsi, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((1 - prZ) * (dmustar * (ewv/(2 * ewu) -
      (sigx7 * mustar + 0.5))/ewvsr - (ewv * pmustar/(2 * ewu) - sigx7 * dmustar)/ewusr) *
      epsi/ewusr + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * dwsr *
      prZ * (epsilon)/ewvsr^3 - sigx8 * sigx2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx1 *
    epsi/ewusr - (sigx2 * sigx9/sigx3 + S * dwsr * (epsilon)/ewvsr^3)) * prZ *
    ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr * mustar/ewusr) - 0.5) -
      0.5 * sigx4) * dmustar * ewvsr/ewusr - (sigx5 * (0.5 * (S * (epsilon)/ewusr) +
      2 * (ewu * ewv/(2 * ewu)^2)) + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) *
      ewu * ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pmustar))/(sigx3 *
      ewusr) - (sigx5 * (1 - prZ) * epsi + 0.5 * (sigx3 * ewusr)) * sigx5/(sigx3 *
      ewusr)^2) * (1 - prZ) * epsi, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((dmustar *
    ewvsr/(4 * (ewu * ewusr)) - 2 * (ewu * pmustar/(2 * ewu)^2)) * ewv - (sigx8 *
    sigx5/sigx3 + (0.5 * (sigx7 * mustar) - 0.25) * dmustar * ewvsr/ewusr + sigx4 *
    (ewv * pmustar/(2 * ewu) - sigx7 * dmustar))) * (1 - prZ) * epsi/(sigx3 *
    ewusr), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx5 * (1 - (1 - prZ) * sigx9/sigx3) * prZ * epsi * ewz/(sigx3 * ewusr),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((ewv * pmustar/(2 * ewu) - sigx7 * dmustar)/2 + (pmustar - sigx7 * dmustar)/2) *
      ewv/ewu - (0.25 * (ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) - sigx7^2 *
      mustar) * dmustar) * (1 - prZ) * epsi/ewusr + prZ * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr * (epsilon)^2/ewvsr^2 -
      0.5 * sigx6)/ewvsr - sigx8^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prZ * (epsi * (ewv * pmustar/(2 * ewu) - sigx7 *
      dmustar)/ewusr - (sigx8 * sigx9/sigx3 + sigx6/ewvsr)) * ewz/sigx3, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (1 + prZ * sigx9/sigx3) * ewz) * prZ * sigx9 *
      ewz/sigx3, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Different sigma_v

## logit specification class membership
chessmnsfexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  mustar <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  dwsr <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  dpmu <- (dmustar/ewv1_h - pmustar/ewu_h)
  deemu <- (dmustar * ewv1_h/ewu_h)
  sigx1 <- (dpmu * epsi * ewz/(wzdeno * ewu_h) + S * prC * dwsr * (epsilon)/ewv2_h^3)
  sigx2 <- (prC * dwsr/ewv2_h + epsi * ewz * pmustar/(wzdeno * ewu_h))
  sigx3 <- (0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx4 <- (wzdeno * ewu_h * pmustar/(wzdeno * ewu_h)^2)
  sigx5 <- ((0.5 * deemu - sigx3 * pmustar)/(wzdeno * ewu_h) - 0.5 * sigx4)
  sigx6 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx7 <- (ewv1 * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewv2_h^2) - 0.5 * dwsr)
  sigx9 <- (1/(wzdeno * ewu_h) - ewu_h * ewz/(wzdeno * ewu_h)^2)
  sigx10 <- (sigx9 * epsi * pmustar - prC * dwsr/(wzdeno * ewv2_h))
  sigx11 <- (sigx2 * wzdeno * ewu_h)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (((mustar/ewv1_h - 1/ewu_h) * dmustar/ewv1_h - dpmu/ewu_h) * epsi *
    ewz/(wzdeno * ewu_h) + prC * dwsr * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 -
    sigx1^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((((0.5 + 0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 *
      ewu)^2)) * pmustar - 0.5 * deemu)/ewu_h + (0.5 * (mustar/ewu_h) - sigx3/ewv1_h) *
      dmustar)/wzdeno + 0.5 * sigx4)/ewu_h - (sigx5 * sigx1/sigx2 + 0.5 * (wzdeno *
      dmustar * ewu_h/((wzdeno * ewu_h)^2 * ewv1_h)))) * epsi * ewz/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmustar * (ewv1/(2 * ewu) - (sigx6 * mustar +
      0.5))/ewv1_h - sigx7/ewu_h)/sigx11 - sigx1 * wzdeno * ewu_h * sigx7/sigx11^2) *
      epsi * ewz, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * prC * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * dwsr * (epsilon)/(sigx2 * ewv2_h^3) - sigx1 * sigx8 * ewv2_h/(sigx2 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx9 *
    dpmu * epsi - (sigx10 * sigx1/sigx2 + S * prC * dwsr * (epsilon)/(wzdeno *
    ewv2_h^3))) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewv1_h * mustar/ewu_h) - 0.5) -
      0.5 * sigx3) * dmustar * ewv1_h/ewu_h - (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) *
      ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewu_h)) * pmustar)/(wzdeno *
      ewu_h) - ((sigx5 * epsi * ewz/sigx2 + 0.5 * (S * (epsilon)/ewu_h) + 2 *
      (ewu * ewv1/(2 * ewu)^2)) * sigx5 + (0.5 * ((0.5 - wzdeno^2 * ewu_h^2/(wzdeno *
      ewu_h)^2) * ewu_h * pmustar + 0.5 * (dmustar * ewv1_h)) + 0.5 * ((0.5 *
      deemu - sigx3 * pmustar) * ewu_h)) * wzdeno/(wzdeno * ewu_h)^2)) * epsi *
      ewz/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((dmustar *
    ewv1_h/(4 * (ewu * ewu_h)) - 2 * (ewu * pmustar/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx6 * mustar) - 0.25) * dmustar * ewv1_h/ewu_h + sigx3 * sigx7))/sigx11 -
    (sigx5 * epsi * ewz + 0.5 * sigx2) * wzdeno * ewu_h * sigx7/sigx11^2) * epsi *
    ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx5 * sigx8 * prC * epsi * ewv2_h * ewz/(sigx2 * ewv2_h)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (0.5 * (sigx9 * dmustar * ewv1_h/ewu_h) - ((((0.5 - wzdeno^2 * ewu_h^2/(wzdeno *
      ewu_h)^2) * ewz + 0.5 * wzdeno) * ewu_h/(wzdeno * ewu_h)^2 + sigx3 *
      sigx9) * pmustar + sigx5 * sigx10 * ewz/sigx2)) * epsi * ewz/sigx2, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((sigx7/2 + (pmustar - sigx6 * dmustar)/2) * ewv1/ewu - (0.25 * (ewv1_h/ewu_h) +
      0.25 * (S * (epsilon)/ewv1_h) - sigx6^2 * mustar) * dmustar)/sigx11 -
      epsi * sigx7^2 * ewz/sigx11^2) * epsi * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx8 * prC * epsi * sigx7 * ewv2_h * ewz/((sigx2 * ewv2_h)^2 *
      wzdeno * ewu_h)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1/(wzdeno * ewu_h) - (sigx10/sigx11 + ewu_h/(wzdeno *
      ewu_h)^2) * ewz) * epsi * sigx7 * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * dwsr * (epsilon)^2/(sigx2 * ewv2_h^3) - (sigx8 * prC + 0.5 *
      (sigx2 * ewv2_h)) * sigx8/(sigx2 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx10 * sigx8/sigx2 + 0.5 * (S^2 * dwsr *
      (epsilon)^2/(wzdeno * ewv2_h^2)))/ewv2_h - 0.5 * (wzdeno * dwsr * ewv2_h/(wzdeno *
      ewv2_h)^2)) * prC * ewz/sigx2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv2_h) + ewv2_h/(wzdeno *
      ewv2_h)^2) * dwsr - (sigx10^2/sigx2 + (2 - 2 * (wzdeno * ewu_h^2 * ewz/(wzdeno *
      ewu_h)^2)) * ewu_h * epsi * pmustar/(wzdeno * ewu_h)^2)) * ewz + sigx9 *
      epsi * pmustar - prC * dwsr/(wzdeno * ewv2_h)) * ewz/sigx2, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmnsfexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  mustar <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  dwsr <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  dpmu <- (dmustar/ewv1_h - pmustar/ewu_h)
  deemu <- (dmustar * ewv1_h/ewu_h)
  sigx1 <- (ewz1 * dpmu * epsi/ewu_h + S * ewz2 * dwsr * (epsilon)/ewv2_h^3)
  sigx2 <- (ewz2 * dwsr/ewv2_h + ewz1 * epsi * pmustar/ewu_h)
  sigx3 <- (0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * deemu - sigx4 * pmustar)
  sigx6 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx7 <- (ewv1 * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewv2_h^2) - 0.5 * dwsr)
  sigx9 <- (epsi * pmustar/ewu_h - dwsr/ewv2_h)
  sigx10 <- (pi * sigx2 * ((Wz)^2 + 1))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (((mustar/ewv1_h - 1/ewu_h) * dmustar/ewv1_h - dpmu/ewu_h) * ewz1 * epsi/ewu_h +
      ewz2 * dwsr * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 - sigx1^2/sigx2)/sigx2,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((0.5 * (S * (epsilon)/ewu_h) + 1 + 2 * (ewu * ewv1/(2 *
      ewu)^2)) * pmustar - 0.5 * deemu)/ewu_h + (0.5 * (mustar/ewu_h) - sigx4/ewv1_h) *
      dmustar)/(sigx2 * ewu_h) - sigx1 * sigx5 * ewu_h/(sigx2 * ewu_h)^2) *
      ewz1 * epsi, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmustar * (ewv1/(2 * ewu) - (sigx6 * mustar +
      0.5))/ewv1_h - sigx7/ewu_h)/(sigx2 * ewu_h) - sigx1 * ewu_h * sigx7/(sigx2 *
      ewu_h)^2) * ewz1 * epsi, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ewz2 * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * dwsr * (epsilon)/(sigx2 * ewv2_h^3) - sigx1 * sigx8 * ewv2_h/(sigx2 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((dpmu *
    epsi/ewu_h - S * dwsr * (epsilon)/ewv2_h^3)/sigx10 - pi * sigx1 * ((Wz)^2 +
    1) * sigx9/sigx10^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewv1_h * mustar/ewu_h) - 0.5) -
      0.5 * sigx4) * dmustar * ewv1_h/ewu_h - (sigx5 * sigx3 + (2 * ((1 - 8 *
      (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewu_h)) *
      pmustar))/(sigx2 * ewu_h) - (sigx5 * ewz1 * epsi + 0.5 * (sigx2 * ewu_h)) *
      sigx5/(sigx2 * ewu_h)^2) * ewz1 * epsi, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((dmustar *
    ewv1_h/(4 * (ewu * ewu_h)) - 2 * (ewu * pmustar/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx6 * mustar) - 0.25) * dmustar * ewv1_h/ewu_h + sigx3 * sigx7))/(sigx2 *
    ewu_h) - (sigx5 * ewz1 * epsi + 0.5 * (sigx2 * ewu_h)) * sigx7/(sigx2 * ewu_h)^2) *
    ewz1 * epsi, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (ewz2 * sigx5 * sigx8 * ewz1 * epsi * ewv2_h/((sigx2 * ewv2_h)^2 * ewu_h)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx5 * (1/sigx10 - pi * ((Wz)^2 + 1) * ewz1 * sigx9/sigx10^2) * epsi/ewu_h,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((sigx7/2 + (pmustar - sigx6 * dmustar)/2) * ewv1/ewu - (0.25 * (ewv1_h/ewu_h) +
      0.25 * (S * (epsilon)/ewv1_h) - sigx6^2 * mustar) * dmustar)/(sigx2 *
      ewu_h) - ewz1 * epsi * sigx7^2/(sigx2 * ewu_h)^2) * ewz1 * epsi, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (ewz2 * sigx8 * ewz1 * epsi * sigx7 * ewv2_h/((sigx2 * ewv2_h)^2 *
      ewu_h)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1/sigx10 - pi * ((Wz)^2 + 1) * ewz1 * sigx9/sigx10^2) *
      epsi * sigx7/ewu_h, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ewz2 * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * dwsr * (epsilon)^2/(sigx2 * ewv2_h^3) - (ewz2 * sigx8 +
      0.5 * (sigx2 * ewv2_h)) * sigx8/(sigx2 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx8 * (1/sigx10 + pi * ((Wz)^2 + 1) * ewz2 *
      sigx9/sigx10^2)/ewv2_h), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (pi * Wz * sigx2) + epsi * pmustar/ewu_h -
      dwsr/ewv2_h) * sigx9/sigx10^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmnsfexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  mustar <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  dwsr <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  dpmu <- (dmustar/ewv1_h - pmustar/ewu_h)
  deemu <- (dmustar * ewv1_h/ewu_h)
  sigx1 <- (dpmu * epsi * pwZ/ewu_h + S * (1 - pwZ) * dwsr * (epsilon)/ewv2_h^3)
  sigx2 <- ((1 - pwZ) * dwsr/ewv2_h + epsi * pmustar * pwZ/ewu_h)
  sigx3 <- (0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * deemu - sigx4 * pmustar)
  sigx6 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx7 <- (ewv1 * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewv2_h^2) - 0.5 * dwsr)
  sigx9 <- (epsi * pmustar/ewu_h - dwsr/ewv2_h)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (((mustar/ewv1_h - 1/ewu_h) * dmustar/ewv1_h - dpmu/ewu_h) * epsi * pwZ/ewu_h +
      (1 - pwZ) * dwsr * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 - sigx1^2/sigx2)/sigx2,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((0.5 * (S * (epsilon)/ewu_h) + 1 + 2 * (ewu * ewv1/(2 *
      ewu)^2)) * pmustar - 0.5 * deemu)/ewu_h + (0.5 * (mustar/ewu_h) - sigx4/ewv1_h) *
      dmustar)/(sigx2 * ewu_h) - sigx1 * sigx5 * ewu_h/(sigx2 * ewu_h)^2) *
      epsi * pwZ, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmustar * (ewv1/(2 * ewu) - (sigx6 * mustar +
      0.5))/ewv1_h - sigx7/ewu_h)/(sigx2 * ewu_h) - sigx1 * ewu_h * sigx7/(sigx2 *
      ewu_h)^2) * epsi * pwZ, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (1 - pwZ) * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * dwsr * (epsilon)/(sigx2 * ewv2_h^3) - sigx1 * sigx8 * ewv2_h/(sigx2 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (dpmu *
    epsi/ewu_h - (sigx1 * sigx9/sigx2 + S * dwsr * (epsilon)/ewv2_h^3)) * dwZ/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewv1_h * mustar/ewu_h) - 0.5) -
      0.5 * sigx4) * dmustar * ewv1_h/ewu_h - (sigx5 * sigx3 + (2 * ((1 - 8 *
      (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewu_h)) *
      pmustar))/(sigx2 * ewu_h) - (sigx5 * epsi * pwZ + 0.5 * (sigx2 * ewu_h)) *
      sigx5/(sigx2 * ewu_h)^2) * epsi * pwZ, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((dmustar *
    ewv1_h/(4 * (ewu * ewu_h)) - 2 * (ewu * pmustar/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx6 * mustar) - 0.25) * dmustar * ewv1_h/ewu_h + sigx3 * sigx7))/(sigx2 *
    ewu_h) - (sigx5 * epsi * pwZ + 0.5 * (sigx2 * ewu_h)) * sigx7/(sigx2 * ewu_h)^2) *
    epsi * pwZ, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx5 * sigx8 * (1 - pwZ) * epsi * ewv2_h * pwZ/((sigx2 * ewv2_h)^2 * ewu_h)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx5 * (1 - sigx9 * pwZ/sigx2) * dwZ * epsi/(sigx2 * ewu_h), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((sigx7/2 + (pmustar - sigx6 * dmustar)/2) * ewv1/ewu - (0.25 * (ewv1_h/ewu_h) +
      0.25 * (S * (epsilon)/ewv1_h) - sigx6^2 * mustar) * dmustar)/(sigx2 *
      ewu_h) - epsi * sigx7^2 * pwZ/(sigx2 * ewu_h)^2) * epsi * pwZ, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx8 * (1 - pwZ) * epsi * sigx7 * ewv2_h * pwZ/((sigx2 *
      ewv2_h)^2 * ewu_h)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx9 * pwZ/sigx2) * dwZ * epsi * sigx7/(sigx2 *
      ewu_h), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - pwZ) * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * dwsr * (epsilon)^2/(sigx2 * ewv2_h^3) - (sigx8 * (1 - pwZ) +
      0.5 * (sigx2 * ewv2_h)) * sigx8/(sigx2 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((1 - pwZ) * sigx9/sigx2 + 1) * sigx8 * dwZ/(sigx2 *
      ewv2_h)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (dwZ * (dwZ * sigx9/sigx2 + Wz) * sigx9/sigx2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmnsfexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  mustar <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmustar <- pnorm(-mustar)
  dmustar <- dnorm(-mustar, 0, 1)
  epsi <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  dwsr <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  dpmu <- (dmustar/ewv1_h - pmustar/ewu_h)
  deemu <- (dmustar * ewv1_h/ewu_h)
  sigx1 <- ((1 - prZ) * dpmu * epsi/ewu_h + S * dwsr * prZ * (epsilon)/ewv2_h^3)
  sigx2 <- ((1 - prZ) * epsi * pmustar/ewu_h + dwsr * prZ/ewv2_h)
  sigx3 <- (0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx4 <- (0.5 + 0.5 * (S * (epsilon)/ewu_h) + 2 * (ewu * ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * deemu - sigx4 * pmustar)
  sigx6 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx7 <- (ewv1 * pmustar/(2 * ewu) - sigx6 * dmustar)
  sigx8 <- (0.5 * (S^2 * dwsr * (epsilon)^2/ewv2_h^2) - 0.5 * dwsr)
  sigx9 <- (epsi * pmustar/ewu_h - dwsr/ewv2_h)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (((mustar/ewv1_h - 1/ewu_h) * dmustar/ewv1_h - dpmu/ewu_h) * (1 - prZ) *
      epsi/ewu_h + dwsr * prZ * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3 -
      sigx1^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((0.5 * (S * (epsilon)/ewu_h) + 1 + 2 * (ewu * ewv1/(2 *
      ewu)^2)) * pmustar - 0.5 * deemu)/ewu_h + (0.5 * (mustar/ewu_h) - sigx4/ewv1_h) *
      dmustar)/(sigx2 * ewu_h) - sigx1 * sigx5 * ewu_h/(sigx2 * ewu_h)^2) *
      (1 - prZ) * epsi, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmustar * (ewv1/(2 * ewu) - (sigx6 * mustar +
      0.5))/ewv1_h - sigx7/ewu_h)/(sigx2 * ewu_h) - sigx1 * ewu_h * sigx7/(sigx2 *
      ewu_h)^2) * (1 - prZ) * epsi, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * prZ * (S * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 -
      2) - 0.5) * dwsr * (epsilon)/(sigx2 * ewv2_h^3) - sigx1 * sigx8 * ewv2_h/(sigx2 *
      ewv2_h)^2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (dpmu *
    epsi/ewu_h - (sigx1 * sigx9/sigx2 + S * dwsr * (epsilon)/ewv2_h^3)) * prZ *
    ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewv1_h * mustar/ewu_h) - 0.5) -
      0.5 * sigx4) * dmustar * ewv1_h/ewu_h - (sigx5 * sigx3 + (2 * ((1 - 8 *
      (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewu_h)) *
      pmustar))/(sigx2 * ewu_h) - (sigx5 * (1 - prZ) * epsi + 0.5 * (sigx2 *
      ewu_h)) * sigx5/(sigx2 * ewu_h)^2) * (1 - prZ) * epsi, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((dmustar *
    ewv1_h/(4 * (ewu * ewu_h)) - 2 * (ewu * pmustar/(2 * ewu)^2)) * ewv1 - ((0.5 *
    (sigx6 * mustar) - 0.25) * dmustar * ewv1_h/ewu_h + sigx3 * sigx7))/(sigx2 *
    ewu_h) - (sigx5 * (1 - prZ) * epsi + 0.5 * (sigx2 * ewu_h)) * sigx7/(sigx2 *
    ewu_h)^2) * (1 - prZ) * epsi, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx5 * sigx8 * (1 - prZ) * prZ * epsi * ewv2_h/((sigx2 * ewv2_h)^2 * ewu_h)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx5 * (1 - (1 - prZ) * sigx9/sigx2) * prZ * epsi * ewz/(sigx2 * ewu_h),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((sigx7/2 + (pmustar - sigx6 * dmustar)/2) * ewv1/ewu - (0.25 * (ewv1_h/ewu_h) +
      0.25 * (S * (epsilon)/ewv1_h) - sigx6^2 * mustar) * dmustar)/(sigx2 *
      ewu_h) - (1 - prZ) * epsi * sigx7^2/(sigx2 * ewu_h)^2) * (1 - prZ) *
    epsi, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx8 * (1 - prZ) * prZ * epsi * sigx7 * ewv2_h/((sigx2 *
      ewv2_h)^2 * ewu_h)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - (1 - prZ) * sigx9/sigx2) * prZ * epsi *
      sigx7 * ewz/(sigx2 * ewu_h), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prZ * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * dwsr * (epsilon)^2/(sigx2 * ewv2_h^3) - (sigx8 * prZ + 0.5 *
      (sigx2 * ewv2_h)) * sigx8/(sigx2 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar), (nXvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx8 * (1 + prZ * sigx9/sigx2) * prZ * ewz/(sigx2 *
      ewv2_h)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (1 + prZ * sigx9/sigx2) * ewz) * prZ * sigx9 *
      ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf exponential-normal distribution
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
zisfexponormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfexponormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfexponormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfexponormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfexponormlike_logit,
    grad = cgradzisfexponormlike_logit, hess = chesszisfexponormlike_logit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfexponormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfexponormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfexponormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfexponormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfexponormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfexponormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initExpo = initExpo))
}

## cauchit specification class membership
zisfexponormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfexponormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfexponormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfexponormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfexponormlike_cauchit,
    grad = cgradzisfexponormlike_cauchit, hess = chesszisfexponormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfexponormlike_cauchit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfexponormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfexponormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfexponormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfexponormlike_cauchit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfexponormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initExpo = initExpo))
}

## probit specification class membership
zisfexponormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfexponormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfexponormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfexponormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfexponormlike_probit,
    grad = cgradzisfexponormlike_probit, hess = chesszisfexponormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfexponormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfexponormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfexponormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfexponormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfexponormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfexponormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initExpo = initExpo))
}

## cloglog specification class membership
zisfexponormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstzisfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(czisfexponormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(czisfexponormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfexponormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisfexponormlike_cloglog,
    grad = cgradzisfexponormlike_cloglog, hess = chesszisfexponormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(czisfexponormlike_cloglog(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(czisfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(czisfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfexponormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chesszisfexponormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfexponormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfexponormlike_cloglog(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfexponormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initExpo = initExpo))
}

# Different sigma_v

## logit specification class membership
mnsfexponormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfexponormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfexponormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfexponormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfexponormlike_logit,
    grad = cgradmnsfexponormlike_logit, hess = chessmnsfexponormlike_logit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfexponormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfexponormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfexponormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfexponormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfexponormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfexponormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initExpo = initExpo))
}

## cauchit specification class membership
mnsfexponormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfexponormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfexponormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfexponormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfexponormlike_cauchit,
    grad = cgradmnsfexponormlike_cauchit, hess = chessmnsfexponormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfexponormlike_cauchit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfexponormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfexponormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfexponormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfexponormlike_cauchit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfexponormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initExpo = initExpo))
}

## probit specification class membership
mnsfexponormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfexponormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfexponormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfexponormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfexponormlike_probit,
    grad = cgradmnsfexponormlike_probit, hess = chessmnsfexponormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfexponormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfexponormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfexponormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfexponormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfexponormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfexponormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initExpo = initExpo))
}

## cloglog specification class membership
mnsfexponormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmnsfexponormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmnsfexponormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmnsfexponormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsfexponormlike_cloglog,
    grad = cgradmnsfexponormlike_cloglog, hess = chessmnsfexponormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmnsfexponormlike_cloglog(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmnsfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmnsfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmnsfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfexponormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmnsfexponormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfexponormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfexponormlike_cloglog(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfexponormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initExpo = initExpo))
}

# Conditional efficiencies estimation ----------
#' efficiencies for zisf exponential-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v

## logit specification class membership
czisfexponormeff_logit <- function(object, level) {
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
  mustar <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sqrt(exp(Wv)) * dnorm(mustar/sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * exp(Wv)) * pnorm(mustar/sqrt(exp(Wv)) - sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * exp(Wv)) * pnorm(mustar/sqrt(exp(Wv)) +
      sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
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
czisfexponormeff_cauchit <- function(object, level) {
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
  mustar <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sqrt(exp(Wv)) * dnorm(mustar/sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * exp(Wv)) * pnorm(mustar/sqrt(exp(Wv)) - sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * exp(Wv)) * pnorm(mustar/sqrt(exp(Wv)) +
      sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
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
czisfexponormeff_probit <- function(object, level) {
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
  mustar <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sqrt(exp(Wv)) * dnorm(mustar/sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * exp(Wv)) * pnorm(mustar/sqrt(exp(Wv)) - sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * exp(Wv)) * pnorm(mustar/sqrt(exp(Wv)) +
      sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
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
czisfexponormeff_cloglog <- function(object, level) {
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
  mustar <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sqrt(exp(Wv)) * dnorm(mustar/sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * exp(Wv)) * pnorm(mustar/sqrt(exp(Wv)) - sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * exp(Wv)) * pnorm(mustar/sqrt(exp(Wv)) +
      sqrt(exp(Wv)))/pnorm(mustar/sqrt(exp(Wv)))
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
cmnsfexponormeff_logit <- function(object, level) {
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
  mustar <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sqrt(exp(Wv1)) * dnorm(mustar/sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * exp(Wv1)) * pnorm(mustar/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * exp(Wv1)) * pnorm(mustar/sqrt(exp(Wv1)) +
      sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
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
cmnsfexponormeff_cauchit <- function(object, level) {
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
  mustar <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sqrt(exp(Wv1)) * dnorm(mustar/sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * exp(Wv1)) * pnorm(mustar/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * exp(Wv1)) * pnorm(mustar/sqrt(exp(Wv1)) +
      sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
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
cmnsfexponormeff_probit <- function(object, level) {
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
  mustar <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sqrt(exp(Wv1)) * dnorm(mustar/sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * exp(Wv1)) * pnorm(mustar/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * exp(Wv1)) * pnorm(mustar/sqrt(exp(Wv1)) +
      sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
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
cmnsfexponormeff_cloglog <- function(object, level) {
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
  mustar <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sqrt(exp(Wv1)) * dnorm(mustar/sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * exp(Wv1)) * pnorm(mustar/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * exp(Wv1)) * pnorm(mustar/sqrt(exp(Wv1)) +
      sqrt(exp(Wv1)))/pnorm(mustar/sqrt(exp(Wv1)))
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
#' marginal impact on efficiencies for zisf exponential-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v

## logit specification class membership
czisfmargexponorm_Eu_logit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargexponorm_Vu_logit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
czisfmargexponorm_Eu_cauchit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargexponorm_Vu_cauchit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
czisfmargexponorm_Eu_probit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargexponorm_Vu_probit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
czisfmargexponorm_Eu_cloglog <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

czisfmargexponorm_Vu_cloglog <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu),
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
cmnsfmargexponorm_Eu_logit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargexponorm_Vu_logit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmnsfmargexponorm_Eu_cauchit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargexponorm_Vu_cauchit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmnsfmargexponorm_Eu_probit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargexponorm_Vu_probit <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmnsfmargexponorm_Eu_cloglog <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu/2),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmnsfmargexponorm_Vu_cloglog <- function(object) {
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
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
