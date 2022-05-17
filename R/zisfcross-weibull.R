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
# Convolution: weibull - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf weibull-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
# Same sigma_v

## logit specification class membership
czisfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  if (k < 0)
    return(NA)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cauchit specification class membership
czisfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  if (k < 0)
    return(NA)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## probit specification class membership
czisfweibullnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  if (k < 0)
    return(NA)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cloglog specification class membership
czisfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  if (k < 0)
    return(NA)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# Different sigma_v

## logit specification class membership
cmnsfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  k <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  if (k < 0)
    return(NA)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cauchit specification class membership
cmnsfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  k <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  if (k < 0)
    return(NA)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## probit specification class membership
cmnsfweibullnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  k <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  if (k < 0)
    return(NA)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cloglog specification class membership
cmnsfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  k <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  if (k < 0)
    return(NA)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf weibull-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
# Same sigma_v
cstzisfweibullnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + weibull - normal distributions...\n")
  initWeibull <- maxLik(logLik = cweibullnormlike, start = cstweibullnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradweibullnormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat)
  Esti <- initWeibull$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1, rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "k",
    paste0("ZI_", colnames(Zvar)))
  names(initWeibull$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "k")
  return(list(StartVal = StartVal, initWeibull = initWeibull))
}

# Different sigma_v
cstmnsfweibullnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + weibull - normal distributions...\n")
  initWeibull <- maxLik(logLik = cweibullnormlike, start = cstweibullnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradweibullnormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat)
  Esti <- initWeibull$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1, rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), "k", paste0("ZI_", colnames(Zvar)))
  names(initWeibull$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "k")
  return(list(StartVal = StartVal, initWeibull = initWeibull))
}

# Gradient of the likelihood function ----------
#' gradient for zisf weibull-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
# Same sigma_v

## logit specification class membership
cgradzisfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewu_h <- exp(Wu/2)
  ewv_h <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  depsi <- dnorm(S * (epsilon)/ewv_h, 0, 1)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewu_h, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewu_h/(k^2 * ewv_h^3),
    FUN = "*")
  Fx1 <- sweep(F1^2, MARGIN = 1, STATS = 1/ewv_h^2, FUN = "*") -
    1
  Fx2 <- sweep(Fx1 * F2, MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  F5 <- sweep(F2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*")
  sDiv <- apply(F5, 1, sum)
  F6 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewu_h/ewv_h^3, FUN = "*")
  vzdeno <- (prC * depsi/ewv_h + ewz * sDiv/(Q * wzdeno))
  sigx1 <- (Q * vzdeno * wzdeno)
  sigx2 <- prC * depsi * (epsilon)/ewv_h^3
  F7 <- sweep(sweep(F2 * F1^2, MARGIN = 1, STATS = 0.5/ewv_h^2,
    FUN = "*") - 0.5 * F2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*")
  XF3 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF3[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz/(Q * wzdeno * vzdeno)
  }
  gx <- XF3 + sweep(Xvar, MARGIN = 1, STATS = S^2 * sigx2/vzdeno,
    FUN = "*")
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz/sigx1
  }
  VF7 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF7[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz/(Q * wzdeno * vzdeno)
  }
  gv <- VF7 + sweep(vHvar, MARGIN = 1, STATS = (0.5 * (S^2 *
    depsi * (epsilon)^2/ewv_h^2) - 0.5 * depsi) * prC/ewv_h/vzdeno,
    FUN = "*")
  Qw <- (1/(Q * wzdeno) - Q * ewz/(Q * wzdeno)^2)
  gradll <- gradll <- cbind(gx, gu, gv, apply(F4, 1, sum) *
    ewz/sigx1, sweep(Zvar, MARGIN = 1, STATS = (Qw * sDiv -
    prC * depsi/(wzdeno * ewv_h)) * ewz/vzdeno, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradzisfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  depsi <- dnorm(S * (epsilon)/ewvsr, 0, 1)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewusr, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewvsr^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewusr/(k^2 * ewvsr^3),
    FUN = "*")
  F5 <- sweep(F2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  sDiv <- apply(F5, 1, sum)
  F6 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr^3, FUN = "*")
  F7 <- sweep(sweep(F2 * F1^2, MARGIN = 1, STATS = 0.5/ewvsr^2,
    FUN = "*") - 0.5 * F2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  sigx1 <- (ewz2 * depsi/ewvsr + ewz1 * sDiv/Q)
  sigx2 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx3 <- (sDiv/Q - depsi/ewvsr)
  sigx4 <- (pi * sigx1 * ((Wz)^2 + 1))
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz1/(Q * sigx1)
  }
  gx <- XF1 + sweep(Xvar, MARGIN = 1, STATS = (S^2 * ewz2 *
    depsi * (epsilon)/ewvsr^3)/sigx1, FUN = "*")
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz1/(Q *
      sigx1)
  }
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz1/(Q * sigx1)
  }
  gv <- VF1 + sweep(vHvar, MARGIN = 1, STATS = (ewz2 * sigx2/ewvsr)/sigx1,
    FUN = "*")
  gradll <- gradll <- cbind(gx, gu, gv, apply(F4, 1, sum) *
    ewz1/(Q * sigx1), sweep(Zvar, MARGIN = 1, STATS = sigx3/sigx4,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradzisfweibullnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr, 0, 1)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewusr, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewvsr^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewusr/(k^2 * ewvsr^3),
    FUN = "*")
  F5 <- sweep(F2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  sDiv <- apply(F5, 1, sum)
  F6 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr^3, FUN = "*")
  F7 <- sweep(sweep(F2 * F1^2, MARGIN = 1, STATS = 0.5/ewvsr^2,
    FUN = "*") - 0.5 * F2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  sigx1 <- ((1 - pwZ) * depsi/ewvsr + pwZ * sDiv/Q)
  sigx2 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx3 <- (sDiv/Q - depsi/ewvsr)
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * pwZ/(Q * sigx1)
  }
  gx <- XF1 + sweep(Xvar, MARGIN = 1, STATS = (1 - pwZ) * depsi *
    (epsilon)/ewvsr^3/sigx1, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * pwZ/(Q *
      sigx1)
  }
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * pwZ/(Q * sigx1)
  }
  gv <- VF1 + sweep(vHvar, MARGIN = 1, STATS = sigx2 * (1 -
    pwZ)/ewvsr/sigx1, FUN = "*")
  gradll <- gradll <- cbind(gx, UF1, gv, apply(F4, 1, sum) *
    pwZ/(Q * sigx1), sweep(Zvar, MARGIN = 1, STATS = dwZ *
    sigx3/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradzisfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  depsi <- dnorm(S * (epsilon)/ewvsr, 0, 1)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewusr, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewvsr^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewusr/(k^2 * ewvsr^3),
    FUN = "*")
  F5 <- sweep(F2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  sDiv <- apply(F5, 1, sum)
  F6 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr^3, FUN = "*")
  F7 <- sweep(sweep(F2 * F1^2, MARGIN = 1, STATS = 0.5/ewvsr^2,
    FUN = "*") - 0.5 * F2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  sigx1 <- ((1 - prZ) * sDiv/Q + depsi * prZ/ewvsr)
  sigx2 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx3 <- (sDiv/Q - depsi/ewvsr)
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/(Q * sigx1)
  }
  gx <- XF1 + sweep(Xvar, MARGIN = 1, STATS = depsi * prZ *
    (epsilon)/ewvsr^3/sigx1, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * (1 - prZ)/(Q *
      sigx1)
  }
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/(Q * sigx1)
  }
  gv <- VF1 + sweep(vHvar, MARGIN = 1, STATS = sigx2 * prZ/ewvsr/sigx1,
    FUN = "*")
  gradll <- cbind(gx, UF1, gv, apply(F4, 1, sum) * (1 - prZ)/(Q *
    sigx1), sweep(Zvar, MARGIN = 1, STATS = prZ * ewz * sigx3/sigx1,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v

## logit specification class membership
cgradmnsfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  k <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewu_h, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewv1_h^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewu_h/(k^2 * ewv1_h^3),
    FUN = "*")
  sDiv <- apply(sweep(F2, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"),
    1, sum)
  F5 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewu_h/ewv1_h^3, FUN = "*")
  sigx1 <- (prC * depsi/ewv2_h + ewz * sDiv/(Q * wzdeno))
  sigx2 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 *
    depsi)
  sigx3 <- ((1/(Q * wzdeno) - Q * ewz/(Q * wzdeno)^2) * sDiv -
    prC * depsi/(wzdeno * ewv2_h))
  F6 <- sweep(F2 * F1^2, MARGIN = 1, STATS = 1/ewv1_h^2, FUN = "*")
  F7 <- sweep((0.5 * F6 - 0.5 * F2), MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*")
  XF3 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF3[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz/(Q * wzdeno * sigx1)
  }
  gx <- XF3 + sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * depsi *
    (epsilon)/ewv2_h^3/sigx1, FUN = "*")
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(-(0.5 * (S * F5)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz/(Q *
      sigx1 * wzdeno)
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz/(Q * sigx1 * wzdeno)
  }
  gradll <- cbind(gx, gu, gv, sweep(vHvar, MARGIN = 1, STATS = sigx2 *
    prC/(sigx1 * ewv2_h), FUN = "*"), apply(F4, 1, sum) *
    ewz/(Q * sigx1 * wzdeno), sweep(Zvar, MARGIN = 1, STATS = sigx3 *
    ewz/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmnsfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  k <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  dv2 <- dnorm(S * epsilon/ewvsr2, 0, 1)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewusr, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewvsr1^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewusr/(k^2 * ewvsr1^3),
    FUN = "*")
  F5 <- sweep(F2, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*")
  sDiv <- apply(F5, 1, sum)
  F6 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr1^3, FUN = "*")
  F7 <- sweep(sweep(F2 * F1^2, MARGIN = 1, STATS = 0.5/ewvsr1^2,
    FUN = "*") - 0.5 * F2, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*")
  sigx1 <- (ewz2 * dv2/ewvsr2 + ewz1 * sDiv/Q)
  sigx2 <- (pi * sigx1 * ((Wz)^2 + 1))
  sigx3 <- (sDiv/Q - dv2/ewvsr2)
  sigx18 <- (0.5 * (S^2 * dv2 * (epsilon)^2/ewvsr2^2) - 0.5 *
    dv2)
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz1/(Q * sigx1)
  }
  gx <- XF1 + sweep(Xvar, MARGIN = 1, STATS = ewz2 * dv2 *
    (epsilon)/ewvsr2^3/sigx1, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz1/(Q *
      sigx1)
  }
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz1/(Q * sigx1)
  }
  gradll <- cbind(gx, UF1, VF1, sweep(vHvar, MARGIN = 1, STATS = ewz2 *
    sigx18/(sigx1 * ewvsr2), FUN = "*"), apply(F4, 1, sum) *
    ewz1/(Q * sigx1), sweep(Zvar, MARGIN = 1, STATS = sigx3/sigx2,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmnsfweibullnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  k <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  dv2 <- dnorm(S * epsilon/ewvsr2, 0, 1)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewusr, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewvsr1^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewusr/(k^2 * ewvsr1^3),
    FUN = "*")
  F5 <- sweep(F2, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*")
  sDiv <- apply(F5, 1, sum)
  F6 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr1^3, FUN = "*")
  F7 <- sweep(sweep(F2 * F1^2, MARGIN = 1, STATS = 0.5/ewvsr1^2,
    FUN = "*") - 0.5 * F2, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*")
  sigx1 <- ((1 - pwZ) * dv2/ewvsr2 + pwZ * sDiv/Q)
  sigx3 <- (sDiv/Q - dv2/ewvsr2)
  sigx18 <- (0.5 * (S^2 * dv2 * (epsilon)^2/ewvsr2^2) - 0.5 *
    dv2)
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * pwZ/(Q * sigx1)
  }
  gx <- XF1 + sweep(Xvar, MARGIN = 1, STATS = (1 - pwZ) * dv2 *
    (epsilon)/ewvsr2^3/sigx1, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * pwZ/(Q *
      sigx1)
  }
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * pwZ/(Q * sigx1)
  }
  gradll <- cbind(gx, UF1, VF1, sweep(vHvar, MARGIN = 1, STATS = sigx18 *
    (1 - pwZ)/(sigx1 * ewvsr2), FUN = "*"), apply(F4, 1,
    sum) * pwZ/(Q * sigx1), sweep(Zvar, MARGIN = 1, STATS = dwZ *
    sigx3/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmnsfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  k <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  dv2 <- dnorm(S * epsilon/ewvsr2, 0, 1)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewusr, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewvsr1^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewusr/(k^2 * ewvsr1^3),
    FUN = "*")
  F5 <- sweep(F2, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*")
  sDiv <- apply(F5, 1, sum)
  F6 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr1^3, FUN = "*")
  F7 <- sweep(sweep(F2 * F1^2, MARGIN = 1, STATS = 0.5/ewvsr1^2,
    FUN = "*") - 0.5 * F2, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*")
  sigx1 <- ((1 - prZ) * sDiv/Q + dv2 * prZ/ewvsr2)
  sigx3 <- (sDiv/Q - dv2/ewvsr2)
  sigx18 <- (0.5 * (S^2 * dv2 * (epsilon)^2/ewvsr2^2) - 0.5 *
    dv2)
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/(Q * sigx1)
  }
  gx <- XF1 + sweep(Xvar, MARGIN = 1, STATS = dv2 * prZ * (epsilon)/ewvsr2^3/sigx1,
    FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * (1 - prZ)/(Q *
      sigx1)
  }
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/(Q * sigx1)
  }
  gradll <- cbind(gx, UF1, VF1, sweep(vHvar, MARGIN = 1, STATS = sigx18 *
    prZ/(sigx1 * ewvsr2), FUN = "*"), apply(F4, 1, sum) *
    (1 - prZ)/(Q * sigx1), sweep(Zvar, MARGIN = 1, STATS = prZ *
    ewz * sigx3/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf weibull-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
# Same sigma_v

## logit specification class membership
chesszisfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewu_h <- exp(Wu/2)
  ewv_h <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  depsi <- dnorm(S * (epsilon)/ewv_h, 0, 1)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewu_h, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewu_h/(k^2 * ewv_h^3),
    FUN = "*")
  Fx1 <- sweep(F1^2, MARGIN = 1, STATS = 1/ewv_h^2, FUN = "*") -
    1
  Fx2 <- sweep(Fx1 * F2, MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  F5 <- sweep(F2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*")
  sDiv <- apply(F5, 1, sum)
  F6 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewu_h/ewv_h^3, FUN = "*")
  vzdeno <- (prC * depsi/ewv_h + ewz * sDiv/(Q * wzdeno))
  sigx1 <- (Q * vzdeno * wzdeno)
  sigx2 <- prC * depsi * (epsilon)/ewv_h^3
  F7 <- sweep(sweep(F2 * F1^2, MARGIN = 1, STATS = 0.5/ewv_h^2,
    FUN = "*") - 0.5 * F2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*")
  F8 <- sweep(F1^2, MARGIN = 1, STATS = 1/ewv_h^2, FUN = "*")
  F9 <- sweep(((-log(1 - FiMat))^(1/k) * F8 - (-log(1 - FiMat))^(1/k)) *
    F2, MARGIN = 1, STATS = ewu_h/ewv_h^3, FUN = "*")
  F10 <- sweep(((-log(1 - FiMat))^(1/k) * F8 - (-log(1 - FiMat))^(1/k)) *
    F2 * log(-log(1 - FiMat)), MARGIN = 1, STATS = ewu_h/(k^2 *
    ewv_h^3), FUN = "*")
  F11 <- sweep((S * (-log(1 - FiMat))^(2/k) * F1), MARGIN = 1,
    STATS = ewu_h/ewv_h^2, FUN = "*")
  F12 <- sweep((S * (-log(1 - FiMat))^(2/k)), MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  Qw <- (1/(Q * wzdeno) - Q * ewz/(Q * wzdeno)^2)
  sigx3 <- (Qw * sDiv - prC * depsi/(wzdeno * ewv_h))
  sigx4 <- (1/(Q * wzdeno) - (sigx3/sigx1 + Q/(Q * wzdeno)^2) *
    ewz)
  F13 <- sweep(((0.5 * (-log(1 - FiMat))^(1/k) - 0.5 * F11) *
    F1 + 0.5 * F12) * F2, MARGIN = 1, STATS = ewu_h/ewv_h^3,
    FUN = "*")
  F14 <- sweep((0.25 * (-log(1 - FiMat))^(1/k) + 0.5 * ((-log(1 -
    FiMat))^(1/k) - 0.5 * ((-log(1 - FiMat))^(1/k) * F8))) *
    F2 * F1, MARGIN = 1, STATS = ewu_h/ewv_h^3, FUN = "*")
  F15 <- sweep(((0.5 * (-log(1 - FiMat))^(1/k) - 0.5 * F11) *
    F1 + 0.5 * F12) * F2 * log(-log(1 - FiMat)), MARGIN = 1,
    STATS = ewu_h/(k^2 * ewv_h^3), FUN = "*")
  F16 <- sweep(((0.5 * (0.5 * (F8) - 1) - 0.25) * F2 * F8 -
    0.5 * (0.5 * (F2 * F8) - 0.5 * F2)), MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*")
  F17 <- sweep((k^2 * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewv_h^3/(k^2 * ewv_h^3)^2, FUN = "*")
  F18 <- sweep(((-log(1 - FiMat))^(1/k) * F1^2), MARGIN = 1,
    STATS = 1/(k^2 * ewv_h^5), FUN = "*")
  F19 <- sweep((0.5 * F18 - 1.5 * F17) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewu_h, FUN = "*")
  F20 <- sweep((-log(1 - FiMat))^(2/k) * F1, MARGIN = 1, STATS = ewu_h/ewv_h^2,
    FUN = "*")
  F21 <- sweep((-log(1 - FiMat))^(2/k), MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  F22 <- sweep((F1 * (S * F20 - (-log(1 - FiMat))^(1/k)) -
    S * F21) * log(-log(1 - FiMat)), MARGIN = 1, STATS = 1/(k^4 *
    ewv_h^3), FUN = "*")
  F23 <- sweep((k * (-log(1 - FiMat))^(1/k) * F1), MARGIN = 1,
    STATS = ewv_h^3/(k^2 * ewv_h^3)^2, FUN = "*")
  F24 <- sweep(S * (F22 - 2 * F23) * F2 * log(-log(1 - FiMat)),
    MARGIN = 1, STATS = ewu_h, FUN = "*")
  XF3 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF3[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz/(Q * wzdeno)
  }
  gx <- XF3 + sweep(Xvar, MARGIN = 1, STATS = S^2 * sigx2,
    FUN = "*")
  XF3_2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF3_2[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
  }
  UF6 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF6[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum)
  }
  UF15 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF15[, k] <- apply(sweep(S * F15, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)
  }
  VF7 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF7[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz/(Q * wzdeno)
  }
  gv <- VF7 + sweep(vHvar, MARGIN = 1, STATS = (0.5 * (S^2 *
    depsi * (epsilon)^2/ewv_h^2) - 0.5 * depsi) * prC/ewv_h,
    FUN = "*")
  VF7_2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF7_2[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  VF19 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF19[, k] <- apply(sweep(S * F19, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  XF10 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF10[, k] <- apply(sweep(S * F10, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz/sigx1
  }
  sF4 <- apply(F4, 1, sum)
  HX1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in 1:Q) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      Fx2[, r] * ewz/(Q * wzdeno)/vzdeno, FUN = "*"), Xvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * F9))[, r] * ewz/sigx1, FUN = "*"), uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      ((0.5 * (F8 - 2) - 0.5) * F3)[, r] * ewz/(Q * wzdeno)/vzdeno,
      FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * F13))[, r] * ewz/sigx1, FUN = "*"), uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      S * F14[, r] * ewz/sigx1, FUN = "*"), vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      F16[, r] * ewz/(Q * wzdeno)/vzdeno, FUN = "*"), vHvar)
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar +
    1, ncol = nXvar + nuZUvar + nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) + crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S^2 * prC * depsi * (S^2 *
      (epsilon)^2/ewv_h^2 - 1)/ewv_h^3/vzdeno, FUN = "*"),
    Xvar) - crossprod(sweep(gx, MARGIN = 1, STATS = wHvar/vzdeno^2,
    FUN = "*"), gx)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- Reduce("+",
    HXU1) - crossprod(sweep(gx, MARGIN = 1, STATS = wHvar *
    ewz * Q * wzdeno/sigx1^2, FUN = "*"), UF6)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- Reduce("+", HXV1) + crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S^2 * (0.5 * (S^2 * (epsilon)^2/ewv_h^2 -
      2) - 0.5) * sigx2/vzdeno, FUN = "*"), vHvar) - crossprod(sweep(gx,
    MARGIN = 1, STATS = wHvar/vzdeno^2, FUN = "*"), gv)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(XF10,
    MARGIN = 1, STATS = wHvar, FUN = "*") - sweep(gx, MARGIN = 1,
    STATS = wHvar * ewz * Q * wzdeno * sF4/sigx1^2, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(XF3_2,
    MARGIN = 1, STATS = wHvar * Qw * ewz/vzdeno, FUN = "*") -
    (sweep(gx, MARGIN = 1, STATS = wHvar * sigx3/vzdeno *
      ewz/vzdeno, FUN = "*") + sweep(Xvar, MARGIN = 1,
      STATS = wHvar * S^2 * prC * depsi * (epsilon)/(wzdeno *
        ewv_h^3) * ewz/vzdeno, FUN = "*")), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- Reduce("+", HU1) - crossprod(sweep(UF6,
    MARGIN = 1, STATS = wHvar * ewz^2/sigx1^2, FUN = "*"),
    UF6)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+", HUV1) -
    crossprod(sweep(UF6, MARGIN = 1, STATS = wHvar * ewz/(vzdeno *
      sigx1), FUN = "*"), gv)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- colSums(sweep(UF15, MARGIN = 1, STATS = wHvar *
    ewz/sigx1, FUN = "*") - sweep(UF6, MARGIN = 1, STATS = wHvar *
    ewz^2 * sF4/sigx1^2, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(UF6,
    MARGIN = 1, STATS = wHvar * sigx4 * ewz/vzdeno, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HV1) + crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv_h^2) -
    1) - 0.25) * depsi * (epsilon)^2/ewv_h^2 - 0.5 * (0.5 *
    (S^2 * depsi * (epsilon)^2/ewv_h^2) - 0.5 * depsi))/ewv_h/vzdeno,
    FUN = "*"), vHvar) - crossprod(sweep(gv, MARGIN = 1,
    STATS = wHvar/vzdeno^2, FUN = "*"), gv)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(VF19,
    MARGIN = 1, STATS = wHvar * ewz/sigx1, FUN = "*") - sweep(gv,
    MARGIN = 1, STATS = wHvar * ewz * Q * wzdeno * sF4/sigx1^2,
    FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar +
      nZHvar + 1)] <- crossprod(sweep(VF7_2, MARGIN = 1,
    STATS = wHvar * Qw * ewz/vzdeno, FUN = "*") - (sweep(gv,
    MARGIN = 1, STATS = wHvar * sigx3/vzdeno * ewz/vzdeno,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (0.5 * (S^2 * depsi * (epsilon)^2/(wzdeno * ewv_h^3)) -
      0.5 * (wzdeno * depsi * ewv_h/(wzdeno * ewv_h)^2)) *
    prC * ewz/vzdeno, FUN = "*")), Zvar)
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(apply(F24, 1, sum) * wHvar * ewz/sigx1 - wHvar *
    ewz^2 * sF4^2/sigx1^2)
  hessll[nXvar + nuZUvar + nvZVvar + 1, (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * sigx4 * ewz * sF4/vzdeno,
    FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv_h) +
      ewv_h/(wzdeno * ewv_h)^2) * depsi - (sigx3^2/vzdeno +
      Q * (2 - 2 * (Q^2 * wzdeno * ewz/(Q * wzdeno)^2)) *
        sDiv/(Q * wzdeno)^2)) * ewz + Qw * sDiv - prC *
      depsi/(wzdeno * ewv_h)) * ewz/vzdeno, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Different sigma_v

## logit specification class membership
chessmnsfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  k <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  depsi <- dnorm(S * (epsilon)/ewv2_h)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewu_h, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"))
  F3 <- sweep(F2 * F1, MARGIN = 1, STATS = 1/ewv1_h^3, FUN = "*")
  F4 <- sweep(S * (-log(1 - FiMat))^(1/k) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewu_h/(k^2 * ewv1_h^3),
    FUN = "*")
  sDiv <- apply(sweep(F2, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"),
    1, sum)
  F5 <- sweep((-log(1 - FiMat))^(1/k) * F2 * F1, MARGIN = 1,
    STATS = ewu_h/ewv1_h^3, FUN = "*")
  sigx1 <- (prC * depsi/ewv2_h + ewz * sDiv/(Q * wzdeno))
  sigx2 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewv2_h^2) - 0.5 *
    depsi)
  sigx3 <- ((1/(Q * wzdeno) - Q * ewz/(Q * wzdeno)^2) * sDiv -
    prC * depsi/(wzdeno * ewv2_h))
  F6 <- sweep(F2 * F1^2, MARGIN = 1, STATS = 1/ewv1_h^2, FUN = "*")
  F7 <- sweep((0.5 * F6 - 0.5 * F2), MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*")
  sigx4 <- prC * depsi * (S^2 * (epsilon)^2/ewv2_h^2 - 1)/ewv2_h^3
  sigx5 <- (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 - 2) - 0.5)
  sigx6 <- (1/(Q * wzdeno) - Q * ewz/(Q * wzdeno)^2)
  sigx7 <- (1/(Q * wzdeno) - (sigx3/(Q * sigx1 * wzdeno) +
    Q/(Q * wzdeno)^2) * ewz)
  sigx8 <- (Q * (sigx1 * ewv2_h)^2 * wzdeno)
  sigx9 <- (ewv2_h * (Q * sigx1 * wzdeno)^2)
  sigx10 <- (1/(wzdeno^2 * ewv2_h) + ewv2_h/(wzdeno * ewv2_h)^2)
  sDiv2 <- apply(F4, 1, sum)
  F8 <- sweep(F1^2, MARGIN = 1, STATS = 1/ewv1_h^2, FUN = "*")
  F9 <- sweep((F8 - 1) * F2, MARGIN = 1, STATS = 1/ewv1_h^3,
    FUN = "*")
  F10 <- sweep(((-log(1 - FiMat))^(1/k) * F8 - (-log(1 - FiMat))^(1/k)) *
    F2, MARGIN = 1, STATS = ewu_h/ewv1_h^3, FUN = "*")
  F11 <- sweep(((-log(1 - FiMat))^(1/k) * F8 - (-log(1 - FiMat))^(1/k)) *
    F2 * log(-log(1 - FiMat)), MARGIN = 1, STATS = ewu_h/(k^2 *
    ewv1_h^3), FUN = "*")
  F12 <- sweep((S * (-log(1 - FiMat))^(2/k)) * F1, MARGIN = 1,
    STATS = ewu_h/ewv1_h^2, FUN = "*")
  F13 <- sweep((S * (-log(1 - FiMat))^(2/k)), MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  F14 <- sweep(((0.5 * (-log(1 - FiMat))^(1/k) - 0.5 * F12) *
    F1 + 0.5 * F13) * F2, MARGIN = 1, STATS = ewu_h/ewv1_h^3,
    FUN = "*")
  F15 <- sweep((0.25 * (-log(1 - FiMat))^(1/k) + 0.5 * ((-log(1 -
    FiMat))^(1/k) - 0.5 * ((-log(1 - FiMat))^(1/k) * F8))) *
    F2 * F1, MARGIN = 1, STATS = ewu_h/ewv1_h^3, FUN = "*")
  F16 <- sweep(((0.5 * (-log(1 - FiMat))^(1/k) - 0.5 * F12) *
    F1 + 0.5 * F13) * F2 * log(-log(1 - FiMat)), MARGIN = 1,
    STATS = ewu_h/(k^2 * ewv1_h^3), FUN = "*")
  F17 <- sweep(((0.5 * (0.5 * (F8) - 1) - 0.25) * F6 - 0.5 *
    (0.5 * (F6) - 0.5 * F2)), MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*")
  F18 <- sweep(((-log(1 - FiMat))^(1/k) * F1^2), MARGIN = 1,
    STATS = 1/(k^2 * ewv1_h^5), FUN = "*")
  F19 <- sweep((k^2 * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewv1_h^3/(k^2 * ewv1_h^3)^2, FUN = "*")
  F20 <- sweep((0.5 * F18 - 1.5 * F19) * F2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = ewu_h, FUN = "*")
  F21 <- sweep((-log(1 - FiMat))^(2/k) * F1, MARGIN = 1, STATS = ewu_h/ewv1_h^2,
    FUN = "*")
  F22 <- sweep((-log(1 - FiMat))^(2/k), MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  F23 <- sweep((F1 * (S * F21 - (-log(1 - FiMat))^(1/k)) -
    S * F22) * log(-log(1 - FiMat)), MARGIN = 1, STATS = 1/(k^4 *
    ewv1_h^3), FUN = "*")
  F24 <- sweep((k * (-log(1 - FiMat))^(1/k) * F1), MARGIN = 1,
    STATS = ewv1_h^3/(k^2 * ewv1_h^3)^2, FUN = "*")
  F25 <- sweep(S * (F23 - 2 * F24) * F2 * log(-log(1 - FiMat)),
    MARGIN = 1, STATS = ewu_h, FUN = "*")
  sDiv3 <- apply(F25, 1, sum)
  XF1 <- matrix(nrow = N, ncol = nXvar)
  XF4 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
    XF4[, k] <- apply(sweep(S * F11, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = ewz/(Q * wzdeno), FUN = "*")
  XF3 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * depsi *
    (epsilon)/ewv2_h^3, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  UF2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(-(0.5 * (S * F5)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum)
    UF2[, k] <- apply(sweep(S * F16, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)
  }
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  VF2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
    VF2[, k] <- apply(sweep(S * F20, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  HX1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in 1:Q) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      F9[, r] * ewz/(Q * sigx1 * wzdeno), FUN = "*"), Xvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * F10))[, r] * ewz/(Q * sigx1 * wzdeno),
      FUN = "*"), uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      ((0.5 * (F8 - 2) - 0.5) * F3)[, r] * ewz/(Q * sigx1 *
      wzdeno), FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * F14))[, r] * ewz/(Q * sigx1 * wzdeno),
      FUN = "*"), uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      S * F15[, r] * ewz/(Q * sigx1 * wzdeno), FUN = "*"),
      vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      F17[, r] * ewz/(Q * sigx1 * wzdeno), FUN = "*"),
      vHvar)
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1, ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) + crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S^2 * sigx4/sigx1, FUN = "*"),
    Xvar) - crossprod(sweep(XF2 + XF3, MARGIN = 1, STATS = wHvar/sigx1^2,
    FUN = "*"), XF2 + XF3)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- Reduce("+",
    HXU1) - crossprod(sweep((XF2 + XF3), MARGIN = 1, STATS = wHvar *
    ewz * Q * wzdeno/(Q * sigx1 * wzdeno)^2, FUN = "*"),
    UF1)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- Reduce("+", HXV1) - crossprod(sweep((XF2 +
    XF3), MARGIN = 1, STATS = wHvar * ewz * Q * wzdeno/(Q *
    sigx1 * wzdeno)^2, FUN = "*"), VF1)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * sigx5 * depsi * (epsilon)/(sigx1 *
      ewv2_h^3) * prC, FUN = "*") - sweep((XF2 + XF3),
    MARGIN = 1, STATS = wHvar * sigx2 * ewv2_h/(sigx1 * ewv2_h)^2 *
      prC, FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + 2 * nvZVvar + 1] <- colSums(sweep(XF4,
    MARGIN = 1, STATS = wHvar * ewz/(Q * sigx1 * wzdeno),
    FUN = "*") - sweep((XF2 + XF3), MARGIN = 1, STATS = wHvar *
    ewz * Q * wzdeno * sDiv2/(Q * sigx1 * wzdeno)^2, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(XF1,
    MARGIN = 1, STATS = wHvar * sigx6 * ewz/sigx1, FUN = "*") -
    (sweep((XF2 + XF3), MARGIN = 1, STATS = wHvar * sigx3/sigx1 *
      ewz/sigx1, FUN = "*") + sweep(Xvar, MARGIN = 1, STATS = wHvar *
      S^2 * prC * depsi * (epsilon)/(wzdeno * ewv2_h^3) *
      ewz/sigx1, FUN = "*")), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- Reduce("+", HU1) - crossprod(sweep(UF1,
    MARGIN = 1, STATS = wHvar * ewz^2/(Q * sigx1 * wzdeno)^2,
    FUN = "*"), UF1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+", HUV1) -
    crossprod(sweep(UF1, MARGIN = 1, STATS = wHvar * ewz^2/(Q *
      sigx1 * wzdeno)^2, FUN = "*"), VF1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(UF1,
    MARGIN = 1, STATS = -wHvar * (sigx2 * prC * ewv2_h *
      ewz/sigx8), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + 2 *
    nvZVvar + 1] <- colSums(sweep(UF2, MARGIN = 1, STATS = wHvar *
    ewz/(Q * sigx1 * wzdeno), FUN = "*") - sweep(UF1, MARGIN = 1,
    STATS = wHvar * ewz^2 * sDiv2/(Q * sigx1 * wzdeno)^2,
    FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- crossprod(sweep(UF1, MARGIN = 1, STATS = wHvar *
    sigx7 * ewz/sigx1, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HV1) - crossprod(sweep(VF1, MARGIN = 1, STATS = wHvar *
    ewz^2/(Q * sigx1 * wzdeno)^2, FUN = "*"), VF1)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(VF1, MARGIN = 1, STATS = -wHvar *
    (sigx2 * prC * ewv2_h * ewz/sigx8), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + 2 * nvZVvar + 1] <- colSums(sweep(VF2,
    MARGIN = 1, STATS = wHvar * ewz/(Q * sigx1 * wzdeno),
    FUN = "*") - sweep(VF1, MARGIN = 1, STATS = wHvar * ewz^2 *
    sDiv2/(Q * sigx1 * wzdeno)^2, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(VF1,
    MARGIN = 1, STATS = wHvar * sigx7 * ewz/sigx1, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * depsi * (epsilon)^2/(sigx1 * ewv2_h^3) -
      (sigx2 * prC + 0.5 * (sigx1 * ewv2_h)) * sigx2/(sigx1 *
        ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), nXvar + nuZUvar + 2 * nvZVvar + 1] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * -(Q * sigx2 * prC * wzdeno *
      ewz * sDiv2/sigx9), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx3 * sigx2/sigx1 +
      0.5 * (S^2 * depsi * (epsilon)^2/(wzdeno * ewv2_h^2)))/ewv2_h -
      0.5 * (wzdeno * depsi * ewv2_h/(wzdeno * ewv2_h)^2)) *
      prC * ewz/sigx1), FUN = "*"), Zvar)
  hessll[nXvar + nuZUvar + 2 * nvZVvar + 1, nXvar + nuZUvar +
    2 * nvZVvar + 1] <- sum(wHvar * (ewz * (sDiv3/(Q * sigx1 *
    wzdeno) - ewz * sDiv2^2/(Q * sigx1 * wzdeno)^2)))
  hessll[nXvar + nuZUvar + 2 * nvZVvar + 1, (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- colSums(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    sigx7 * ewz * sDiv2/sigx1, FUN = "*"))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * sigx10 * depsi -
      (sigx3^2/sigx1 + Q * (2 - 2 * (Q^2 * wzdeno * ewz/(Q *
        wzdeno)^2)) * sDiv/(Q * wzdeno)^2)) * ewz + sigx6 *
      sDiv - prC * depsi/(wzdeno * ewv2_h)) * ewz/sigx1,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf weibull-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
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
zisfweibullnormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfweibullnormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfweibullnormlike_logit,
      grad = cgradzisfweibullnormlike_logit, hess = chesszisfweibullnormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = function(parm) as(-chesszisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hess = function(parm) -chesszisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfweibullnormlike_logit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- chesszisfweibullnormlike_logit(mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfweibullnormlike_logit(mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfweibullnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfweibullnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cauchit specification class membership
zisfweibullnormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfweibullnormlike_cauchit(startVal,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfweibullnormlike_cauchit,
      grad = cgradzisfweibullnormlike_cauchit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradzisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfweibullnormlike_cauchit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- czisfweibullnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfweibullnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## probit specification class membership
zisfweibullnormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfweibullnormlike_probit(startVal,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfweibullnormlike_probit,
      grad = cgradzisfweibullnormlike_probit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradzisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfweibullnormlike_probit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- czisfweibullnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfweibullnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cloglog specification class membership
zisfweibullnormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfweibullnormlike_cloglog(startVal,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfweibullnormlike_cloglog,
      grad = cgradzisfweibullnormlike_cloglog, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradzisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfweibullnormlike_cloglog(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- czisfweibullnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfweibullnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

# Different sigma_v

## logit specification class membership
mnsfweibullnormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfweibullnormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfweibullnormlike_logit,
      grad = cgradmnsfweibullnormlike_logit, hess = chessmnsfweibullnormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = function(parm) as(-chessmnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hess = function(parm) -chessmnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfweibullnormlike_logit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- chessmnsfweibullnormlike_logit(mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfweibullnormlike_logit(mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfweibullnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfweibullnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cauchit specification class membership
mnsfweibullnormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfweibullnormlike_cauchit(startVal,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfweibullnormlike_cauchit,
      grad = cgradmnsfweibullnormlike_cauchit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfweibullnormlike_cauchit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmnsfweibullnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfweibullnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## probit specification class membership
mnsfweibullnormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfweibullnormlike_probit(startVal,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfweibullnormlike_probit,
      grad = cgradmnsfweibullnormlike_probit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfweibullnormlike_probit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmnsfweibullnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfweibullnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cloglog specification class membership
mnsfweibullnormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfweibullnormlike_cloglog(startVal,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfweibullnormlike_cloglog,
      grad = cgradmnsfweibullnormlike_cloglog, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfweibullnormlike_cloglog(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmnsfweibullnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfweibullnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

# Conditional efficiencies estimation ----------
#' efficiencies for zisf weibull-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v

## logit specification class membership
czisfweibullnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
  }
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
    }

    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
czisfweibullnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
  }
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
    }

    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
czisfweibullnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
  }
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
    }

    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
czisfweibullnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
  }
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
    }

    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmnsfweibullnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
  }
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
    }

    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmnsfweibullnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
  }
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
    }

    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmnsfweibullnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
  }
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
    }

    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmnsfweibullnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
  }
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
    }

    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
#' marginal impact on efficiencies for zisf weibull-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v

## logit specification class membership
czisfmargweibullnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmargweibullnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
czisfmargweibullnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmargweibullnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## probit specification class membership
czisfmargweibullnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmargweibullnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
czisfmargweibullnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmargweibullnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# Different sigma_v

## logit specification class membership
cmnsfmargweibullnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmargweibullnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmnsfmargweibullnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmargweibullnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmnsfmargweibullnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmargweibullnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmnsfmargweibullnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmargweibullnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Pi1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}
