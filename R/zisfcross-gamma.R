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
# Convolution: gamma - normal                                                  #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf gamma-normal distribution
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
czisfgammanormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cauchit specification class membership
czisfgammanormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## probit specification class membership
czisfgammanormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cloglog specification class membership
czisfgammanormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# Different sigma_v

## logit specification class membership
cmnsfgammanormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cauchit specification class membership
cmnsfgammanormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## probit specification class membership
cmnsfgammanormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

## cloglog specification class membership
cmnsfgammanormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf gamma-normal distribution
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
cstzisfgammanorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + gamma - normal distributions...\n")
  initGamma <- maxLik(logLik = cgammanormlike, start = cstgammanorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradgammanormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat)
  Esti <- initGamma$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "P",
    paste0("ZI_", colnames(Zvar)))
  names(initGamma$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "P")
  return(list(StartVal = StartVal, initGamma = initGamma))
}

# Different sigma_v
cstmnsfgammanorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + gamma - normal distributions...\n")
  initGamma <- maxLik(logLik = cgammanormlike, start = cstgammanorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradgammanormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat)
  Esti <- initGamma$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), "P", paste0("ZI_", colnames(Zvar)))
  names(initGamma$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "P")
  return(list(StartVal = StartVal, initGamma = initGamma))
}

# Gradient of the likelihood function ----------
#' gradient for zisf gamma-normal distribution
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
cgradzisfgammanormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewv <- exp(Wv)
  ewv_h <- exp(Wv/2)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewu_p <- exp(-(Wu * P/2))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- exp(ewv/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx2 <- (ewv/ewu_h + S * (epsilon))
  pmusig <- pnorm(-(ewv_h/ewu_h + S * (epsilon)/ewv_h))
  dmusig <- dnorm(sigx2/ewv_h, 0, 1)
  pepsi <- pnorm(sigx2/ewv_h)
  depsi <- dnorm(-(ewv_h/ewu_h + S * (epsilon)/ewv_h), 0, 1)
  depsiv <- dnorm(S * (epsilon)/ewv_h)
  sigx3 <- (depsi/ewv_h - pmusig/ewu_h)
  sigx4 <- (Q * wzdeno * gamma(P))
  sigx5 <- (depsi * ewv_h/ewu_h)
  sigx6 <- (0.5 * (S * (epsilon)/ewu_h) + 0.5 * P + 2 * (ewu *
    ewv/(2 * ewu)^2))
  sigx7 <- (0.5 * sigx5 - sigx6 * pmusig)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig, FUN = "*")
  F3 <- sweep(sweep((qFi), MARGIN = 1, STATS = ewv_h, FUN = "*"),
    MARGIN = 1, STATS = sigx2, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sigx12 <- (prC * depsiv/ewv_h + ewu_p * sigx1 * ewz * pmusig *
    sDiv/sigx4)
  sigx8 <- (Q * sigx12 * wzdeno * gamma(P))
  sigx9 <- (ewv * pmusig/(2 * ewu) - (0.5 * (ewv_h/ewu_h) -
    0.5 * (S * (epsilon)/ewv_h)) * depsi)
  sigx10 <- (0.5 * (S^2 * depsiv * (epsilon)^2/ewv_h^2) - 0.5 *
    depsiv)
  sigx11 <- (1/sigx4 - Q * ewz * gamma(P)/sigx4^2)
  sigx13 <- (ewv/ewu_h - 0.5 * sigx2)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx14 <- (sF3 - 0.5 * (Wu * sDiv))
  sigx15 <- (sigx14/sigx4 - Q * wzdeno * digamma(P) * gamma(P) *
    sDiv/sigx4^2)
  sigx16 <- (sigx11 * ewu_p * sigx1 * pmusig * sDiv - prC *
    depsiv/(wzdeno * ewv_h))
  F4 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig *
    sigx13, FUN = "*")
  F5 <- sweep(sweep(qFi, MARGIN = 1, STATS = 0.5 * (ewv_h),
    FUN = "*"), MARGIN = 1, STATS = ewv/ewu_h, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx3 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx4,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * depsiv *
    (epsilon)/ewv_h^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx12, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep((0.5 - 0.5 * (F2)) * F3^(P -
      2) * (P - 1), MARGIN = 1, STATS = ewv/ewu_h * uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx7 * sDiv, FUN = "*")
  gu <- sweep(UF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx8,
    FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F4 + F5) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv <- sweep(VF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx4/sigx12,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = sigx10 *
    prC/ewv_h/sigx12, FUN = "*")
  gradll <- cbind(gx, gu, gv, sigx15 * ewu_p * sigx1 * ewz *
    pmusig/sigx12, sweep(Zvar, MARGIN = 1, STATS = sigx16 *
    ewz/sigx12, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradzisfgammanormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewusrp <- exp(-(Wu * P/2))
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigx1 <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2 <- (ewv/ewusr + S * (epsilon))
  pmusig <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  dmusig <- dnorm(sigx2/ewvsr, 0, 1)
  pepsi <- pnorm(sigx2/ewvsr)
  depsi <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr), 0, 1)
  depsiv <- dnorm(S * (epsilon)/ewvsr)
  sigx3 <- (depsi/ewvsr - pmusig/ewusr)
  sigx4 <- (Q * gamma(P))
  sigx5 <- (depsi * ewvsr/ewusr)
  sigx6 <- (0.5 * (S * (epsilon)/ewusr) + 0.5 * P + 2 * (ewu *
    ewv/(2 * ewu)^2))
  sigx7 <- (0.5 * sigx5 - sigx6 * pmusig)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig, FUN = "*")
  F3 <- sweep(sweep((qFi), MARGIN = 1, STATS = ewvsr, FUN = "*"),
    MARGIN = 1, STATS = sigx2, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sigx8 <- (ewz2 * depsiv/ewvsr + ewz1 * ewusrp * sigx1 * pmusig *
    sDiv/sigx4)
  sigx9 <- (ewv * pmusig/(2 * ewu) - (0.5 * (ewvsr/ewusr) -
    0.5 * (S * (epsilon)/ewvsr)) * depsi)
  sigx10 <- (0.5 * (S^2 * depsiv * (epsilon)^2/ewvsr^2) - 0.5 *
    depsiv)
  sigx11 <- (ewusrp * sigx1 * pmusig * sDiv/sigx4 - depsiv/ewvsr)
  sigx12 <- (Q * sigx8 * gamma(P))
  sigx13 <- (ewv/ewusr - 0.5 * sigx2)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx14 <- (sF3 - 0.5 * (Wu * sDiv))
  sigx15 <- (sigx14/sigx4 - Q * digamma(P) * gamma(P) * sDiv/sigx4^2)
  sigx16 <- (pi * sigx8 * ((Wz)^2 + 1))
  F4 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig *
    sigx13, FUN = "*")
  F5 <- sweep(sweep(qFi, MARGIN = 1, STATS = 0.5 * (ewvsr),
    FUN = "*"), MARGIN = 1, STATS = ewv/ewusr, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx3 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewz1 * ewusrp * sigx1/sigx4,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * ewz2 * depsiv *
    (epsilon)/ewvsr^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx8, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep((0.5 - 0.5 * (F2)) * F3^(P -
      2) * (P - 1), MARGIN = 1, STATS = ewv/ewusr * uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx7 * sDiv, FUN = "*")
  gu <- sweep(UF2, MARGIN = 1, STATS = ewz1 * ewusrp * sigx1/sigx12,
    FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F4 + F5) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv <- sweep(VF2, MARGIN = 1, STATS = ewz1 * ewusrp * sigx1/sigx4/sigx8,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = ewz2 *
    sigx10/ewvsr/sigx8, FUN = "*")
  gradll <- cbind(gx, gu, gv, sigx15 * ewz1 * ewusrp * sigx1 *
    pmusig/sigx8, sweep(Zvar, MARGIN = 1, STATS = sigx11/sigx16,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradzisfgammanormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewusrp <- exp(-(Wu * P/2))
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigx1 <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2 <- (ewv/ewusr + S * (epsilon))
  pmusig <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  dmusig <- dnorm(sigx2/ewvsr, 0, 1)
  pepsi <- pnorm(sigx2/ewvsr)
  depsi <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr), 0, 1)
  depsiv <- dnorm(S * (epsilon)/ewvsr)
  sigx3 <- (depsi/ewvsr - pmusig/ewusr)
  sigx4 <- (Q * gamma(P))
  sigx5 <- (depsi * ewvsr/ewusr)
  sigx6 <- (0.5 * (S * (epsilon)/ewusr) + 0.5 * P + 2 * (ewu *
    ewv/(2 * ewu)^2))
  sigx7 <- (0.5 * sigx5 - sigx6 * pmusig)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig, FUN = "*")
  F3 <- sweep(sweep((qFi), MARGIN = 1, STATS = ewvsr, FUN = "*"),
    MARGIN = 1, STATS = sigx2, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sigx8 <- ((1 - pwZ) * depsiv/ewvsr + ewusrp * sigx1 * pmusig *
    pwZ * sDiv/sigx4)
  sigx9 <- (ewv * pmusig/(2 * ewu) - (0.5 * (ewvsr/ewusr) -
    0.5 * (S * (epsilon)/ewvsr)) * depsi)
  sigx10 <- (0.5 * (S^2 * depsiv * (epsilon)^2/ewvsr^2) - 0.5 *
    depsiv)
  sigx11 <- (ewusrp * sigx1 * pmusig * sDiv/sigx4 - depsiv/ewvsr)
  sigx12 <- (Q * sigx8 * gamma(P))
  sigx13 <- (ewv/ewusr - 0.5 * sigx2)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx14 <- (sF3 - 0.5 * (Wu * sDiv))
  sigx15 <- (sigx14/sigx4 - Q * digamma(P) * gamma(P) * sDiv/sigx4^2)
  F4 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig *
    sigx13, FUN = "*")
  F5 <- sweep(sweep(qFi, MARGIN = 1, STATS = 0.5 * (ewvsr),
    FUN = "*"), MARGIN = 1, STATS = ewv/ewusr, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx3 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewusrp * sigx1 * pwZ/sigx4,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * (1 - pwZ) *
    depsiv * (epsilon)/ewvsr^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx8, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep((0.5 - 0.5 * (F2)) * F3^(P -
      2) * (P - 1), MARGIN = 1, STATS = ewv/ewusr * uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx7 * sDiv, FUN = "*")
  gu <- sweep(UF2, MARGIN = 1, STATS = ewusrp * sigx1 * pwZ/sigx12,
    FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F4 + F5) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv <- sweep(VF2, MARGIN = 1, STATS = ewusrp * sigx1 * pwZ/sigx4/sigx8,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = sigx10 *
    (1 - pwZ)/ewvsr/sigx8, FUN = "*")
  gradll <- cbind(gx, gu, gv, sigx15 * ewusrp * sigx1 * pmusig *
    pwZ/sigx8, sweep(Zvar, MARGIN = 1, STATS = dwZ * sigx11/sigx8,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradzisfgammanormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewusrp <- exp(-(Wu * P/2))
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigx1 <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2 <- (ewv/ewusr + S * (epsilon))
  pmusig <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  dmusig <- dnorm(sigx2/ewvsr, 0, 1)
  pepsi <- pnorm(sigx2/ewvsr)
  depsi <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr), 0, 1)
  depsiv <- dnorm(S * (epsilon)/ewvsr)
  sigx3 <- (depsi/ewvsr - pmusig/ewusr)
  sigx4 <- (Q * gamma(P))
  sigx5 <- (depsi * ewvsr/ewusr)
  sigx6 <- (0.5 * (S * (epsilon)/ewusr) + 0.5 * P + 2 * (ewu *
    ewv/(2 * ewu)^2))
  sigx7 <- (0.5 * sigx5 - sigx6 * pmusig)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig, FUN = "*")
  F3 <- sweep(sweep((qFi), MARGIN = 1, STATS = ewvsr, FUN = "*"),
    MARGIN = 1, STATS = sigx2, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sigx8 <- ((1 - prZ) * ewusrp * sigx1 * pmusig * sDiv/sigx4 +
    depsiv * prZ/ewvsr)
  sigx9 <- (ewv * pmusig/(2 * ewu) - (0.5 * (ewvsr/ewusr) -
    0.5 * (S * (epsilon)/ewvsr)) * depsi)
  sigx10 <- (0.5 * (S^2 * depsiv * (epsilon)^2/ewvsr^2) - 0.5 *
    depsiv)
  sigx11 <- (ewusrp * sigx1 * pmusig * sDiv/sigx4 - depsiv/ewvsr)
  sigx12 <- (Q * sigx8 * gamma(P))
  sigx13 <- (ewv/ewusr - 0.5 * sigx2)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx14 <- (sF3 - 0.5 * (Wu * sDiv))
  sigx15 <- (sigx14/sigx4 - Q * digamma(P) * gamma(P) * sDiv/sigx4^2)
  F4 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig *
    sigx13, FUN = "*")
  F5 <- sweep(sweep(qFi, MARGIN = 1, STATS = 0.5 * (ewvsr),
    FUN = "*"), MARGIN = 1, STATS = ewv/ewusr, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx3 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = (1 - prZ) * ewusrp *
    sigx1/sigx4, FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * depsiv * prZ *
    (epsilon)/ewvsr^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx8, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep((0.5 - 0.5 * (F2)) * F3^(P -
      2) * (P - 1), MARGIN = 1, STATS = ewv/ewusr * uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx7 * sDiv, FUN = "*")
  gu <- sweep(UF2, MARGIN = 1, STATS = (1 - prZ) * ewusrp *
    sigx1/sigx12, FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F4 + F5) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv <- sweep(VF2, MARGIN = 1, STATS = (1 - prZ) * ewusrp *
    sigx1/sigx4/sigx8, FUN = "*") + sweep(vHvar, MARGIN = 1,
    STATS = sigx10 * prZ/ewvsr/sigx8, FUN = "*")
  gradll <- cbind(gx, gu, gv, sigx15 * (1 - prZ) * ewusrp *
    sigx1 * pmusig/sigx8, sweep(Zvar, MARGIN = 1, STATS = sigx11 *
    prZ * ewz/sigx8, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v

## logit specification class membership
cgradmnsfgammanormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewu_h <- exp(Wu/2)
  ewu <- exp(Wu)
  ewz <- exp(Wz)
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewuP <- exp(-(Wu * P/2))
  ewuv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  ewvuepsi <- (ewv1/ewu_h + S * (epsilon))
  depsi <- dnorm(ewvuepsi/ewv1_h, 0, 1)
  pepsi <- pnorm(ewvuepsi/ewv1_h)
  dvsi <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  ev1epsi <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmuv <- pnorm(-ev1epsi)
  dmuv <- dnorm(-ev1epsi, 0, 1)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi, FUN = "*")
  F3 <- sweep(sweep(qFi, MARGIN = 1, STATS = ewv1_h, FUN = "*"),
    MARGIN = 1, STATS = ewvuepsi, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx1 <- (dmuv/ewv1_h - pmuv/ewu_h)
  sigx2 <- (Q * wzdeno * gamma(P))
  sigx3 <- (prC * dvsi/ewv2_h + ewuP * ewuv * ewz * pmuv *
    sDiv/sigx2)
  sigx4 <- (0.5 * (S * (epsilon)/ewu_h) + 0.5 * P + 2 * (ewu *
    ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmuv * ewv1_h/ewu_h) - sigx4 * pmuv)
  sigx6 <- (Q * sigx3 * wzdeno * gamma(P))
  sigx7 <- (ewv1/ewu_h - 0.5 * ewvuepsi)
  sigx8 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx9 <- (ewv1 * pmuv/(2 * ewu) - sigx8 * dmuv)
  sigx10 <- (0.5 * (S^2 * dvsi * (epsilon)^2/ewv2_h^2) - 0.5 *
    dvsi)
  sigx11 <- (1/sigx2 - Q * ewz * gamma(P)/sigx2^2)
  sigx12 <- (sigx11 * ewuP * ewuv * pmuv * sDiv - prC * dvsi/(wzdeno *
    ewv2_h))
  sigx13 <- ((sF3 - 0.5 * (Wu * sDiv))/sigx2 - Q * wzdeno *
    digamma(P) * gamma(P) * sDiv/sigx2^2)
  F4 <- sweep((0.5 - 0.5 * (F2)) * (F3)^(P - 2) * (P - 1),
    MARGIN = 1, STATS = ewv1/ewu_h, FUN = "*")
  F5 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi * sigx7,
    FUN = "*")
  F6 <- sweep(qFi, MARGIN = 1, STATS = ewv1_h, FUN = "*")
  F7 <- sweep(0.5 * F6, MARGIN = 1, STATS = ewv1/ewu_h, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * (F3)^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx1 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewuP * ewuv * ewz/sigx2,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * dvsi *
    (epsilon)/ewv2_h^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx3, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(F4, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * sDiv, FUN = "*")
  gu <- sweep(UF2, MARGIN = 1, STATS = ewuP * ewuv * ewz/sigx6,
    FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F5 + F7) * (F3)^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv1 <- sweep(VF2, MARGIN = 1, STATS = ewuP * ewuv * ewz/sigx6,
    FUN = "*")
  gradll <- cbind(gx, gu, gv1, sweep(vHvar, MARGIN = 1, STATS = sigx10 *
    prC/(sigx3 * ewv2_h), FUN = "*"), sigx13 * ewuP * ewuv *
    ewz * pmuv/sigx3, sweep(Zvar, MARGIN = 1, STATS = sigx12 *
    ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmnsfgammanormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewusr <- exp(Wu/2)
  ewu <- exp(Wu)
  ewv1 <- exp(Wv1)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  ewusrp <- exp(-(Wu * P/2))
  ewvuv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  ewvuepsi <- (ewv1/ewusr + S * (epsilon))
  depsi <- dnorm(ewvuepsi/ewvsr1, 0, 1)
  pepsi <- pnorm(ewvuepsi/ewvsr1)
  dvsi <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  ev1epsi <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  pmuv <- pnorm(-ev1epsi)
  dmuv <- dnorm(-ev1epsi, 0, 1)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi, FUN = "*")
  F3 <- sweep(sweep(qFi, MARGIN = 1, STATS = ewvsr1, FUN = "*"),
    MARGIN = 1, STATS = ewvuepsi, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx1 <- (dmuv/ewvsr1 - pmuv/ewusr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 0.5 * P + 2 * (ewu *
    ewv1/(2 * ewu)^2))
  sigx2 <- (ewz2 * dvsi/ewvsr2 + ewz1 * ewusrp * ewvuv * pmuv *
    sDiv/sigx4)
  sigx3 <- (pi * sigx2 * ((Wz)^2 + 1))
  sigx5 <- (0.5 * (dmuv * ewvsr1/ewusr) - sigx4 * pmuv)
  sigx6 <- (Q * sigx2 * gamma(P))
  sigx7 <- (ewv1/ewusr - 0.5 * ewvuepsi)
  sigx8 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx9 <- (ewv1 * pmuv/(2 * ewu) - sigx8 * dmuv)
  sigx10 <- (0.5 * (S^2 * dvsi * (epsilon)^2/ewvsr2^2) - 0.5 *
    dvsi)
  sigx11 <- ((sF3 - 0.5 * (Wu * sDiv))/sigx4 - Q * digamma(P) *
    gamma(P) * sDiv/sigx4^2)
  sigx12 <- (ewusrp * ewvuv * pmuv * sDiv/sigx4 - dvsi/ewvsr2)
  sigx13 <- (Q * gamma(P))
  sigx14 <- (ewz2 * dvsi/ewvsr2 + ewz1 * ewusrp * ewvuv * pmuv *
    sDiv/sigx13)
  sigx15 <- (Q * sigx14 * gamma(P))
  F4 <- sweep((0.5 - 0.5 * (F2)) * (F3)^(P - 2) * (P - 1),
    MARGIN = 1, STATS = ewv1/ewusr, FUN = "*")
  F5 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi * sigx7,
    FUN = "*")
  F6 <- sweep(qFi, MARGIN = 1, STATS = ewvsr1, FUN = "*")
  F7 <- sweep(0.5 * F6, MARGIN = 1, STATS = ewv1/ewusr, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * (F3)^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx1 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewz1 * ewusrp * ewvuv/sigx13,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * ewz2 * dvsi *
    (epsilon)/ewvsr2^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx14, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(F4, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * sDiv, FUN = "*")
  gu <- sweep(UF2, MARGIN = 1, STATS = ewz1 * ewusrp * ewvuv/sigx15,
    FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F5 + F7) * (F3)^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv1 <- sweep(VF2, MARGIN = 1, STATS = ewz1 * ewusrp * ewvuv/sigx15,
    FUN = "*")
  gradll <- cbind(gx, gu, gv1, sweep(vHvar, MARGIN = 1, STATS = ewz2 *
    sigx10/(sigx14 * ewvsr2), FUN = "*"), ((sF3 - 0.5 * (Wu *
    sDiv))/sigx13 - Q * digamma(P) * gamma(P) * sDiv/sigx13^2) *
    ewz1 * ewusrp * ewvuv * pmuv/sigx14, sweep(Zvar, MARGIN = 1,
    STATS = (ewusrp * ewvuv * pmuv * sDiv/sigx13 - dvsi/ewvsr2)/(pi *
      sigx14 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmnsfgammanormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewusr <- exp(Wu/2)
  ewu <- exp(Wu)
  ewv1 <- exp(Wv1)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  ewusrp <- exp(-(Wu * P/2))
  ewvuv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  ewvuepsi <- (ewv1/ewusr + S * (epsilon))
  depsi <- dnorm(ewvuepsi/ewvsr1, 0, 1)
  pepsi <- pnorm(ewvuepsi/ewvsr1)
  dvsi <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  ev1epsi <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  pmuv <- pnorm(-ev1epsi)
  dmuv <- dnorm(-ev1epsi, 0, 1)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi, FUN = "*")
  F3 <- sweep(sweep(qFi, MARGIN = 1, STATS = ewvsr1, FUN = "*"),
    MARGIN = 1, STATS = ewvuepsi, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx1 <- (dmuv/ewvsr1 - pmuv/ewusr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 0.5 * P + 2 * (ewu *
    ewv1/(2 * ewu)^2))
  sigx2 <- (sF3 - 0.5 * (Wu * sDiv))
  sigx5 <- (0.5 * (dmuv * ewvsr1/ewusr) - sigx4 * pmuv)
  sigx7 <- (ewv1/ewusr - 0.5 * ewvuepsi)
  sigx8 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx9 <- (ewv1 * pmuv/(2 * ewu) - sigx8 * dmuv)
  sigx10 <- (0.5 * (S^2 * dvsi * (epsilon)^2/ewvsr2^2) - 0.5 *
    dvsi)
  sigx11 <- (sigx2/sigx4 - Q * digamma(P) * gamma(P) * sDiv/sigx4^2)
  sigx12 <- (ewusrp * ewvuv * pmuv * sDiv/sigx4 - dvsi/ewvsr2)
  sigx13 <- (Q * gamma(P))
  sigx14 <- ((1 - pwZ) * dvsi/ewvsr2 + ewusrp * ewvuv * pmuv *
    pwZ * sDiv/sigx13)
  sigx15 <- (Q * sigx14 * gamma(P))
  F4 <- sweep((0.5 - 0.5 * (F2)) * (F3)^(P - 2) * (P - 1),
    MARGIN = 1, STATS = ewv1/ewusr, FUN = "*")
  F5 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi * sigx7,
    FUN = "*")
  F6 <- sweep(qFi, MARGIN = 1, STATS = ewvsr1, FUN = "*")
  F7 <- sweep(0.5 * F6, MARGIN = 1, STATS = ewv1/ewusr, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * (F3)^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx1 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewusrp * ewvuv * pwZ/sigx13,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * (1 - pwZ) *
    dvsi * (epsilon)/ewvsr2^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx14, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(F4, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * sDiv, FUN = "*")
  gu <- sweep(UF2, MARGIN = 1, STATS = ewusrp * ewvuv * pwZ/sigx15,
    FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F5 + F7) * (F3)^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv1 <- sweep(VF2, MARGIN = 1, STATS = ewusrp * ewvuv * pwZ/sigx15,
    FUN = "*")
  gradll <- cbind(gx, gu, gv1, sweep(vHvar, MARGIN = 1, STATS = sigx10 *
    (1 - pwZ)/(sigx14 * ewvsr2), FUN = "*"), (sigx2/sigx13 -
    Q * digamma(P) * gamma(P) * sDiv/sigx13^2) * ewusrp *
    ewvuv * pmuv * pwZ/sigx14, sweep(Zvar, MARGIN = 1, STATS = dwZ *
    (ewusrp * ewvuv * pmuv * sDiv/sigx13 - dvsi/ewvsr2)/sigx14,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmnsfgammanormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewusr <- exp(Wu/2)
  ewu <- exp(Wu)
  ewv1 <- exp(Wv1)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  ewusrp <- exp(-(Wu * P/2))
  ewvuv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  ewvuepsi <- (ewv1/ewusr + S * (epsilon))
  depsi <- dnorm(ewvuepsi/ewvsr1, 0, 1)
  pepsi <- pnorm(ewvuepsi/ewvsr1)
  dvsi <- dnorm(S * (epsilon)/ewvsr2, 0, 1)
  ev1epsi <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  pmuv <- pnorm(-ev1epsi)
  dmuv <- dnorm(-ev1epsi, 0, 1)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi, FUN = "*")
  F3 <- sweep(sweep(qFi, MARGIN = 1, STATS = ewvsr1, FUN = "*"),
    MARGIN = 1, STATS = ewvuepsi, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx1 <- (dmuv/ewvsr1 - pmuv/ewusr)
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 0.5 * P + 2 * (ewu *
    ewv1/(2 * ewu)^2))
  sigx2 <- (sF3 - 0.5 * (Wu * sDiv))
  sigx5 <- (0.5 * (dmuv * ewvsr1/ewusr) - sigx4 * pmuv)
  sigx7 <- (ewv1/ewusr - 0.5 * ewvuepsi)
  sigx8 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx9 <- (ewv1 * pmuv/(2 * ewu) - sigx8 * dmuv)
  sigx10 <- (0.5 * (S^2 * dvsi * (epsilon)^2/ewvsr2^2) - 0.5 *
    dvsi)
  sigx11 <- (sigx2/sigx4 - Q * digamma(P) * gamma(P) * sDiv/sigx4^2)
  sigx12 <- (ewusrp * ewvuv * pmuv * sDiv/sigx4 - dvsi/ewvsr2)
  sigx13 <- (Q * gamma(P))
  sigx14 <- ((1 - prZ) * ewusrp * ewvuv * pmuv * sDiv/sigx13 +
    dvsi * prZ/ewvsr2)
  sigx15 <- (Q * sigx14 * gamma(P))
  F4 <- sweep((0.5 - 0.5 * (F2)) * (F3)^(P - 2) * (P - 1),
    MARGIN = 1, STATS = ewv1/ewusr, FUN = "*")
  F5 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi * sigx7,
    FUN = "*")
  F6 <- sweep(qFi, MARGIN = 1, STATS = ewvsr1, FUN = "*")
  F7 <- sweep(0.5 * F6, MARGIN = 1, STATS = ewv1/ewusr, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * (F3)^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx1 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = (1 - prZ) * ewusrp *
    ewvuv/sigx13, FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * dvsi * prZ *
    (epsilon)/ewvsr2^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx14, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(F4, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * sDiv, FUN = "*")
  gu <- sweep(UF2, MARGIN = 1, STATS = (1 - prZ) * ewusrp *
    ewvuv/sigx15, FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F5 + F7) * (F3)^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv1 <- sweep(VF2, MARGIN = 1, STATS = (1 - prZ) * ewusrp *
    ewvuv/sigx15, FUN = "*")
  gradll <- cbind(gx, gu, gv1, sweep(vHvar, MARGIN = 1, STATS = sigx10 *
    prZ/(sigx14 * ewvsr2), FUN = "*"), (sigx2/sigx13 - Q *
    digamma(P) * gamma(P) * sDiv/sigx13^2) * (1 - prZ) *
    ewusrp * ewvuv * pmuv/sigx14, sweep(Zvar, MARGIN = 1,
    STATS = (ewusrp * ewvuv * pmuv * sDiv/sigx13 - dvsi/ewvsr2) *
      prZ * ewz/sigx14, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf gamma-normal distribution
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
chesszisfgammanormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewv <- exp(Wv)
  ewv_h <- exp(Wv/2)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewu_p <- exp(-(Wu * P/2))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- exp(ewv/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx2 <- (ewv/ewu_h + S * (epsilon))
  musig <- (ewv_h/ewu_h + S * (epsilon)/ewv_h)
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(sigx2/ewv_h, 0, 1)
  pepsi <- pnorm(sigx2/ewv_h)
  depsi <- dnorm(-musig, 0, 1)
  depsiv <- dnorm(S * (epsilon)/ewv_h)
  sigx3 <- (depsi/ewv_h - pmusig/ewu_h)
  sigx4 <- (Q * wzdeno * gamma(P))
  sigx5 <- (depsi * ewv_h/ewu_h)
  sigx6 <- (0.5 * (S * (epsilon)/ewu_h) + 0.5 * P + 2 * (ewu *
    ewv/(2 * ewu)^2))
  sigx7 <- (0.5 * sigx5 - sigx6 * pmusig)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig, FUN = "*")
  F3 <- sweep(sweep((qFi), MARGIN = 1, STATS = ewv_h, FUN = "*"),
    MARGIN = 1, STATS = sigx2, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sigx12 <- (prC * depsiv/ewv_h + ewu_p * sigx1 * ewz * pmusig *
    sDiv/sigx4)
  sigx8 <- (Q * sigx12 * wzdeno * gamma(P))
  sigx9 <- (ewv * pmusig/(2 * ewu) - (0.5 * (ewv_h/ewu_h) -
    0.5 * (S * (epsilon)/ewv_h)) * depsi)
  sigx10 <- (0.5 * (S^2 * depsiv * (epsilon)^2/ewv_h^2) - 0.5 *
    depsiv)
  sigx11 <- (1/sigx4 - Q * ewz * gamma(P)/sigx4^2)
  sigx13 <- (ewv/ewu_h - 0.5 * sigx2)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx14 <- (sF3 - 0.5 * (Wu * sDiv))
  sigx15 <- (sigx14/sigx4 - Q * wzdeno * digamma(P) * gamma(P) *
    sDiv/sigx4^2)
  sigx16 <- (sigx11 * ewu_p * sigx1 * pmusig * sDiv - prC *
    depsiv/(wzdeno * ewv_h))
  F4 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig *
    sigx13, FUN = "*")
  F5 <- sweep(sweep(qFi, MARGIN = 1, STATS = 0.5 * (ewv_h),
    FUN = "*"), MARGIN = 1, STATS = ewv/ewu_h, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx3 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx4,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * depsiv *
    (epsilon)/ewv_h^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx12, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep((0.5 - 0.5 * (F2)) * F3^(P -
      2) * (P - 1), MARGIN = 1, STATS = ewv/ewu_h * uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx7 * sDiv, FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F4 + F5) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv <- sweep(VF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx4/sigx12,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = sigx10 *
    prC/ewv_h/sigx12, FUN = "*")
  F6 <- sweep(-sweep((1 - FiMat) * F1 * qFi/F1^2, MARGIN = 1,
    STATS = dmusig, FUN = "*"), MARGIN = 1, STATS = sigx2/ewv_h,
    FUN = "+")
  F7 <- sweep(F6 * (1 - FiMat) * F3^(P - 2)/F1, MARGIN = 1,
    STATS = dmusig/ewv_h, FUN = "*")
  XF5 <- sweep(XF1, MARGIN = 1, STATS = sigx3, FUN = "*")
  XF6 <- sweep(Xvar, MARGIN = 1, STATS = S * (musig/ewv_h -
    1/ewu_h) * sDiv, FUN = "*")
  XF7 <- sweep((XF6 + XF1), MARGIN = 1, STATS = depsi/ewv_h,
    FUN = "*")
  XF8 <- sweep(XF2, MARGIN = 1, STATS = 1/ewu_h, FUN = "*")
  F8 <- sweep(((0.5 - 0.5 * (F2)) * (1 - F2) * F3^(P - 3) *
    (P - 2) - 0.5 * (F7)) * (P - 1), MARGIN = 1, STATS = ewv/ewu_h,
    FUN = "*")
  UF3 <- sweep(UF1, MARGIN = 1, STATS = depsi/ewv_h, FUN = "*") -
    sweep(UF2, MARGIN = 1, STATS = 1/ewu_h, FUN = "*")
  F9 <- sweep(F6 * (1 - FiMat) * F3^(P - 2)/F1, MARGIN = 1,
    STATS = dmusig * sigx13/ewv_h, FUN = "*")
  XF9 <- sweep(Xvar, MARGIN = 1, STATS = S * depsi/ewv_h, FUN = "*")
  XF10 <- sweep((XF3 + XF4), MARGIN = 1, STATS = pmusig/sigx12,
    FUN = "*")
  XF11 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF11[, k] <- apply(sweep(S * (F3^(P - 2) + F3^(P - 2) *
      log(F3) * (P - 1)) * (1 - F2), MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
  }
  XF12 <- sweep(XF1, MARGIN = 1, STATS = (Wu), FUN = "*")
  F10 <- sweep(((1 - FiMat) * F1 * qFi/F1^2), MARGIN = 1, STATS = dmusig,
    FUN = "*")
  F11 <- sweep((-0.5 * F10), MARGIN = 1, STATS = 0.5 * (sigx2/ewv_h),
    FUN = "+")
  F12 <- sweep((F11 * (1 - FiMat) * F3^(P - 2)/F1), MARGIN = 1,
    STATS = dmusig/ewv_h, FUN = "*")
  F13 <- sweep(((0.5 - 0.5 * (F2))^2 * F3^(P - 3) * (P - 2) -
    0.5 * F12), MARGIN = 1, STATS = ewv/ewu_h, FUN = "*")
  F14 <- sweep((F13 - 0.5 * ((0.5 - 0.5 * (F2)) * F3^(P - 2))),
    MARGIN = 1, STATS = ewv * (P - 1)/ewu_h, FUN = "*")
  F15 <- sweep(F10, MARGIN = 1, STATS = sigx2/ewv_h, FUN = "-")
  F16 <- sweep(sweep(F15, MARGIN = 1, STATS = sigx13^2/ewv_h,
    FUN = "*"), MARGIN = 1, STATS = 0.5 * (ewv/ewu_h), FUN = "+")
  F17 <- 0.5 * (sweep(0.5 * (qFi), MARGIN = 1, STATS = ewv_h,
    FUN = "*") + F4)
  F18 <- sweep((F16 * F2 + F17), MARGIN = 1, STATS = ewv/ewu_h,
    FUN = "-")
  F19 <- sweep(F11, MARGIN = 1, STATS = sigx13/ewv_h, FUN = "*")
  F20 <- sweep((((F19 - 0.5) * F2 + 0.5) * (F3)^(P - 2) + (F4 +
    F5) * (0.5 - 0.5 * F2) * (F3)^(P - 3) * (P - 2)), MARGIN = 1,
    STATS = ewv * (P - 1)/ewu_h, FUN = "*")
  HX1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in 1:Q) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      (S^2 * ((1 - F2)^2 * F3^(P - 3) * (P - 2) - F7) *
        (P - 1))[, r] * pmusig * ewu_p * sigx1 * ewz/sigx4/sigx12,
      FUN = "*"), Xvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      (S * F8)[, r] * pmusig/sigx8 * ewu_p * sigx1 * ewz,
      FUN = "*"), uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      (S * ((F4 + F5) * (1 - F2) * F3^(P - 3) * (P - 2) +
        F9) * (P - 1))[, r] * pmusig * ewu_p * sigx1 *
      ewz/sigx4/sigx12, FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      (F14)[, r] * pmusig/sigx8 * ewu_p * sigx1 * ewz,
      FUN = "*"), uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      F20[, r] * pmusig * ewu_p * sigx1 * ewz/sigx8, FUN = "*"),
      vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      ((F18 * (F3)^(P - 2) + (F4 + F5)^2 * (F3)^(P - 3) *
        (P - 2)) * (P - 1))[, r] * pmusig * ewu_p * sigx1 *
      ewz/sigx4/sigx12, FUN = "*"), vHvar)
  }
  UF4 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF4[, k] <- apply(sweep(((F3)^(P - 2) + (F3)^(P - 2) *
      log(F3) * (P - 1)) * (0.5 - 0.5 * F2), MARGIN = 1,
      STATS = uHvar[, k] * ewv/ewu_h, FUN = "*"), 1, sum)
  }
  VF3 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF3[, k] <- apply(sweep((F4 + F5) * ((F3)^(P - 2) + (F3)^(P -
      2) * log(F3) * (P - 1)), MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar +
    1, ncol = nXvar + nuZUvar + nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) + crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * ewu_p * sigx1 * ewz/sigx4/sigx12,
    FUN = "*"), S * (XF5 + XF7 - XF8)) + crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S^2 * prC * depsiv * (S^2 *
      (epsilon)^2/ewv_h^2 - 1)/ewv_h^3/sigx12, FUN = "*"),
    Xvar) - crossprod(sweep(XF3 + XF4, MARGIN = 1, STATS = wHvar/sigx12^2,
    FUN = "*"), XF3 + XF4)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- Reduce("+",
    HXU1) + crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ewu_p * sigx1 * ewz/sigx8, FUN = "*"), S * UF3) + crossprod(sweep(XF1,
    MARGIN = 1, STATS = wHvar * sigx7/sigx8 * ewu_p * sigx1 *
      ewz, FUN = "*") + sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((0.5 * (musig/ewu_h) - sigx6/ewv_h) * depsi + 0.5 *
    (pmusig/ewu_h)) * sDiv/sigx8 * ewu_p * sigx1 * ewz, FUN = "*"),
    uHvar) - crossprod(sweep((XF3 + XF4), MARGIN = 1, STATS = wHvar *
    Q * wzdeno * gamma(P)/sigx8^2 * ewu_p * sigx1 * ewz,
    FUN = "*"), UF2)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- Reduce("+", HXV1) + crossprod(S * Xvar,
    sweep(VF1, MARGIN = 1, STATS = wHvar * depsi/ewv_h *
      ewu_p * sigx1 * ewz/sigx4/sigx12, FUN = "*") - sweep(VF2,
      MARGIN = 1, STATS = wHvar * ewu_p * sigx1 * ewz/sigx4/ewu_h/sigx12,
      FUN = "*")) + crossprod(sweep(XF1, MARGIN = 1, STATS = wHvar *
    sigx9 * ewu_p * sigx1 * ewz/sigx4/sigx12, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = wHvar * S * depsi * (ewv/(2 *
      ewu) - ((0.5 * (ewv_h/ewu_h) - 0.5 * (S * (epsilon)/ewv_h)) *
      musig + 0.5)) * sDiv/ewv_h * ewu_p * sigx1 * ewz/sigx4/sigx12,
      FUN = "*"), vHvar) + crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (0.5 * (S^2 * (epsilon)^2/ewv_h^2 -
      2) - 0.5) * prC * depsiv * (epsilon)/ewv_h^3/sigx12,
    FUN = "*"), vHvar) - crossprod(sweep((XF3 + XF4), MARGIN = 1,
    STATS = wHvar/sigx12, FUN = "*"), gv)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep((XF9 -
    XF10), MARGIN = 1, STATS = wHvar * sigx15 * ewu_p * sigx1 *
    ewz/sigx12, FUN = "*") + (sweep((XF11 - 0.5 * XF12),
    MARGIN = 1, STATS = wHvar * pmusig/sigx4 * ewu_p * sigx1 *
      ewz/sigx12, FUN = "*") - (sweep(XF1, MARGIN = 1,
    STATS = wHvar * Q * wzdeno * digamma(P) * gamma(P)/sigx4^2 *
      pmusig * ewu_p * sigx1 * ewz/sigx12, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = wHvar * S * sigx15/ewu_h *
      pmusig * ewu_p * sigx1 * ewz/sigx12, FUN = "*"))))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(XF2,
    MARGIN = 1, STATS = wHvar * sigx11 * ewu_p * sigx1 *
      ewz/sigx12, FUN = "*") - (sweep(gx, MARGIN = 1, STATS = wHvar *
    sigx16 * ewz/sigx12, FUN = "*") + sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * prC * depsiv * (epsilon)/(wzdeno *
      ewv_h^3) * ewz/sigx12, FUN = "*")), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- Reduce("+", HU1) + crossprod(uHvar, sweep(UF1,
    MARGIN = 1, STATS = wHvar * (depsi * ewv_h/ewu_h - sigx6 *
      pmusig)/sigx8 * ewu_p * sigx1 * ewz, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = wHvar * ((0.5 * (0.5 *
      (ewv_h * musig/ewu_h) - 0.5) - 0.5 * sigx6) * depsi *
      ewv_h/ewu_h - (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) *
      ewu * ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewu_h)) *
      pmusig) * sDiv/sigx8 * ewu_p * sigx1 * ewz, FUN = "*") -
    sweep(UF2, MARGIN = 1, STATS = wHvar * sigx6/sigx8 *
      ewu_p * sigx1 * ewz, FUN = "*")) - crossprod(sweep(UF2,
    MARGIN = 1, STATS = wHvar * ewu_p * sigx1 * ewz/sigx8^2 *
      ewu_p * sigx1 * ewz, FUN = "*"), UF2)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+", HUV1) +
    crossprod(uHvar, sweep(VF1, MARGIN = 1, STATS = wHvar *
      0.5 * (depsi * ewv_h/ewu_h) * ewu_p * sigx1 * ewz/sigx8,
      FUN = "*") - sweep(VF2, MARGIN = 1, STATS = wHvar *
      sigx6 * ewu_p * sigx1 * ewz/sigx8, FUN = "*")) +
    crossprod(sweep(UF1, MARGIN = 1, STATS = wHvar * sigx9 *
      ewu_p * sigx1 * ewz/sigx8, FUN = "*") + sweep(uHvar,
      MARGIN = 1, STATS = wHvar * ((depsi * ewv_h/(4 *
        (ewu * ewu_h)) - 2 * (ewu * pmusig/(2 * ewu)^2)) *
        ewv - (0.5 * ((0.5 * (ewv_h/ewu_h) - 0.5 * (S *
        (epsilon)/ewv_h)) * musig) - 0.25) * depsi *
        ewv_h/ewu_h) * sDiv * ewu_p * sigx1 * ewz/sigx8,
      FUN = "*"), vHvar) - crossprod(UF2, sweep(VF2, MARGIN = 1,
    STATS = wHvar * ewu_p * sigx1 * ewz/sigx4/sigx12 * ewu_p *
      sigx1 * ewz/sigx8, FUN = "*") + sweep(vHvar, MARGIN = 1,
    STATS = wHvar * sigx10 * prC/ewv_h/sigx12 * ewu_p * sigx1 *
      ewz/sigx8, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- colSums(sweep(UF4, MARGIN = 1, STATS = wHvar *
    pmusig * ewu_p * sigx1 * ewz/sigx8, FUN = "*") + sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx7 * sF3 * ewu_p * sigx1 *
      ewz/sigx8, FUN = "*") - (sweep(UF2, MARGIN = 1, STATS = wHvar *
    0.5 * Wu * ewu_p * sigx1 * ewz/sigx8, FUN = "*") + sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 0.5 * (pmusig * sDiv) * ewu_p *
      sigx1 * ewz/sigx8, FUN = "*")) - sweep(UF2, MARGIN = 1,
    STATS = wHvar * Q * (sigx12 * digamma(P) + sigx15 * ewu_p *
      sigx1 * ewz * pmusig) * wzdeno * gamma(P)/sigx8^2 *
      ewu_p * sigx1 * ewz, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(UF2,
    sweep(Zvar, MARGIN = 1, STATS = wHvar * (1/sigx4 - (sigx16/sigx8 +
      Q * gamma(P)/sigx4^2) * ewz) * ewu_p * sigx1 * ewz/sigx12,
      FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HV1) + crossprod(vHvar, sweep(VF1, MARGIN = 1, STATS = wHvar *
    (ewv * pmusig/(2 * ewu) - 2 * ((0.5 * (ewv_h/ewu_h) -
      0.5 * (S * (epsilon)/ewv_h)) * depsi)) * ewu_p *
    sigx1 * ewz/sigx4/sigx12, FUN = "*") + sweep(VF2, MARGIN = 1,
    STATS = wHvar * ewv/(2 * ewu) * ewu_p * sigx1 * ewz/sigx4/sigx12,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (ewv * (pmusig - (0.5 * (ewv_h/ewu_h) - 0.5 * (S * (epsilon)/ewv_h)) *
      depsi)/(2 * ewu) - (0.25 * (ewv_h/ewu_h) + 0.25 *
      (S * (epsilon)/ewv_h) - (0.5 * (ewv_h/ewu_h) - 0.5 *
      (S * (epsilon)/ewv_h))^2 * musig) * depsi) * sDiv *
    ewu_p * sigx1 * ewz/sigx4/sigx12, FUN = "*")) + crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon)^2/ewv_h^2) - 1) - 0.25) * depsiv *
      (epsilon)^2/ewv_h^2 - 0.5 * sigx10)/ewv_h/sigx12,
    FUN = "*"), vHvar) - crossprod(sweep(VF2, MARGIN = 1,
    STATS = wHvar * ewu_p * sigx1 * ewz/sigx4/sigx12, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = wHvar * sigx10 * prC/ewv_h/sigx12,
      FUN = "*"), sweep(VF2, MARGIN = 1, STATS = ewu_p *
    sigx1 * ewz/sigx4/sigx12, FUN = "*") + sweep(vHvar, MARGIN = 1,
    STATS = sigx10 * prC/ewv_h/sigx12, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(VF3,
    MARGIN = 1, STATS = wHvar * pmusig/sigx4 * ewu_p * sigx1 *
      ewz/sigx12, FUN = "*") + sweep(vHvar, MARGIN = 1,
    STATS = wHvar * sigx9 * sF3/sigx4 * ewu_p * sigx1 * ewz/sigx12,
    FUN = "*") - sweep(0.5 * (VF2), MARGIN = 1, STATS = wHvar *
    Wu/sigx4 * ewu_p * sigx1 * ewz/sigx12, FUN = "*") - (sweep(VF2,
    MARGIN = 1, STATS = wHvar * ewu_p * sigx1 * ewz/sigx4 *
      sigx15 * pmusig/sigx12 * ewu_p * sigx1 * ewz/sigx12,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = wHvar *
    sigx10 * prC/ewv_h * sigx15 * pmusig/sigx12 * ewu_p *
    sigx1 * ewz/sigx12, FUN = "*") + sweep(VF2, MARGIN = 1,
    STATS = wHvar * Q * wzdeno * digamma(P) * gamma(P)/sigx4^2 *
      ewu_p * sigx1 * ewz/sigx12, FUN = "*")))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar +
      nZHvar + 1)] <- crossprod(sweep(VF2, MARGIN = 1,
    STATS = wHvar * sigx11 * ewu_p * sigx1 * ewz/sigx12,
    FUN = "*") - (sweep(gv, MARGIN = 1, STATS = wHvar * sigx16 *
    ewz/sigx12, FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (0.5 * (S^2 * depsiv * (epsilon)^2/(wzdeno * ewv_h^3)) -
      0.5 * (wzdeno * depsiv * ewv_h/(wzdeno * ewv_h)^2)) *
    prC * ewz/sigx12, FUN = "*")), Zvar)
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(((apply((F3)^(P - 1) * log(F3)^2, 1, sum) -
    0.5 * (Wu * sF3)) * ewu_p * sigx1 * ewz * pmusig/sigx12/sigx4 -
    ((sigx15 * ewu_p * sigx1 * ewz * pmusig/sigx12 + 0.5 *
      (Wu)) * sigx15 + Q * ((2 * sF3 - (0.5 * (Wu) + 2 *
      (Q^2 * wzdeno^2 * digamma(P) * gamma(P)^2/sigx4^2)) *
      sDiv) * digamma(P) + (digamma(P)^2 + trigamma(P)) *
      sDiv) * wzdeno * gamma(P)/sigx4^2) * ewu_p * sigx1 *
      ewz * pmusig/sigx12) * wHvar)
  hessll[nXvar + nuZUvar + nvZVvar + 1, (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (sigx11 * sF3 - (sigx16 *
      sigx15 * ewz/sigx12 + (0.5 * (Wu * sigx11) + Q *
      ((2 - 2 * (Q^2 * wzdeno^2 * gamma(P)^2/sigx4^2)) *
        ewz + 1) * digamma(P) * gamma(P)/sigx4^2) * sDiv)) *
      ewu_p * sigx1 * ewz * pmusig/sigx12, FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv_h) +
      ewv_h/(wzdeno * ewv_h)^2) * depsiv - (sigx16^2/sigx12 +
      Q * (2 - 2 * (Q^2 * wzdeno * ewz * gamma(P)^2/sigx4^2)) *
        ewu_p * sigx1 * gamma(P) * pmusig * sDiv/sigx4^2)) *
      ewz + sigx11 * ewu_p * sigx1 * pmusig * sDiv - prC *
      depsiv/(wzdeno * ewv_h)) * ewz/sigx12, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Different sigma_v

## logit specification class membership
chessmnsfgammanormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewu_h <- exp(Wu/2)
  ewu <- exp(Wu)
  ewz <- exp(Wz)
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewuP <- exp(-(Wu * P/2))
  ewuv <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  ewvuepsi <- (ewv1/ewu_h + S * (epsilon))
  depsi <- dnorm(ewvuepsi/ewv1_h, 0, 1)
  pepsi <- pnorm(ewvuepsi/ewv1_h)
  dvsi <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  ev1epsi <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  pmuv <- pnorm(-ev1epsi)
  dmuv <- dnorm(-ev1epsi, 0, 1)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi, FUN = "*")
  F3 <- sweep(sweep(qFi, MARGIN = 1, STATS = ewv1_h, FUN = "*"),
    MARGIN = 1, STATS = ewvuepsi, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sF3lo <- apply(F3^(P - 1) * log(F3)^2, 1, sum)
  sigx1 <- (dmuv/ewv1_h - pmuv/ewu_h)
  sigx2 <- (Q * wzdeno * gamma(P))
  sigx3 <- (prC * dvsi/ewv2_h + ewuP * ewuv * ewz * pmuv *
    sDiv/sigx2)
  sigx4 <- (0.5 * (S * (epsilon)/ewu_h) + 0.5 * P + 2 * (ewu *
    ewv1/(2 * ewu)^2))
  sigx5 <- (0.5 * (dmuv * ewv1_h/ewu_h) - sigx4 * pmuv)
  sigx6 <- (Q * sigx3 * wzdeno * gamma(P))
  sigx7 <- (ewv1/ewu_h - 0.5 * ewvuepsi)
  sigx8 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx9 <- (ewv1 * pmuv/(2 * ewu) - sigx8 * dmuv)
  sigx10 <- (0.5 * (S^2 * dvsi * (epsilon)^2/ewv2_h^2) - 0.5 *
    dvsi)
  sigx11 <- (1/sigx2 - Q * ewz * gamma(P)/sigx2^2)
  sigx12 <- (sigx11 * ewuP * ewuv * pmuv * sDiv - prC * dvsi/(wzdeno *
    ewv2_h))
  sigx13 <- ((sF3 - 0.5 * (Wu * sDiv))/sigx2 - Q * wzdeno *
    digamma(P) * gamma(P) * sDiv/sigx2^2)
  F4 <- sweep((0.5 - 0.5 * (F2)) * (F3)^(P - 2) * (P - 1),
    MARGIN = 1, STATS = ewv1/ewu_h, FUN = "*")
  F5 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = depsi * sigx7,
    FUN = "*")
  F6 <- sweep(qFi, MARGIN = 1, STATS = ewv1_h, FUN = "*")
  F7 <- sweep(0.5 * F6, MARGIN = 1, STATS = ewv1/ewu_h, FUN = "-")
  sigx14 <- (ev1epsi/ewv1_h - 1/ewu_h)
  sigx15 <- (S^2 * (epsilon)^2/ewv2_h^2 - 1)
  sigx16 <- ((0.5 * (ev1epsi/ewu_h) - sigx4/ewv1_h) * dmuv +
    0.5 * (pmuv/ewu_h))
  sigx17 <- (ewv1/(2 * ewu) - (sigx8 * ev1epsi + 0.5))
  sigx18 <- (0.5 * (S^2 * (epsilon)^2/ewv2_h^2 - 2) - 0.5)
  sigx19 <- (dmuv * ewv1_h/ewu_h - sigx4 * pmuv)
  sigx20 <- (0.5 * (0.5 * (ewv1_h * ev1epsi/ewu_h) - 0.5) -
    0.5 * sigx4)
  sigx21 <- ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 *
    ewu)^2)
  sigx22 <- (sigx20 * dmuv * ewv1_h/ewu_h - (2 * sigx21 - 0.25 *
    (S * (epsilon)/ewu_h)) * pmuv)
  sigx23 <- (dmuv * ewv1_h/(4 * (ewu * ewu_h)) - 2 * (ewu *
    pmuv/(2 * ewu)^2))
  sigx24 <- (sigx23 * ewv1 - (0.5 * (sigx8 * ev1epsi) - 0.25) *
    dmuv * ewv1_h/ewu_h)
  sigx25 <- (Q * (sigx3 * ewv2_h)^2 * wzdeno * gamma(P))
  sigx26 <- (1/sigx2 - (sigx12/sigx6 + Q * gamma(P)/sigx2^2) *
    ewz)
  sigx27 <- ((2 - 2 * (Q^2 * wzdeno^2 * gamma(P)^2/sigx2^2)) *
    ewz + 1)
  sigx28 <- (sigx11 * sF3 - (sigx12 * sigx13 * ewz/sigx3 +
    (0.5 * (Wu * sigx11) + Q * sigx27 * digamma(P) * gamma(P)/sigx2^2) *
      sDiv))
  sigx29 <- (sigx13 * ewuP * ewuv * ewz * pmuv/sigx3 + 0.5 *
    (Wu))
  sigx30 <- (2 * sF3 - (0.5 * (Wu) + 2 * (Q^2 * wzdeno^2 *
    digamma(P) * gamma(P)^2/sigx2^2)) * sDiv)
  F8 <- sweep((1 - FiMat) * F1 * qFi/F1^2, MARGIN = 1, STATS = depsi,
    FUN = "*")
  F9 <- sweep((-F8), MARGIN = 1, STATS = ewvuepsi/ewv1_h, FUN = "+")
  F10 <- sweep(F9 * (1 - FiMat) * (F3)^(P - 2)/(F1), MARGIN = 1,
    STATS = depsi/ewv1_h, FUN = "*")
  F11 <- sweep(((0.5 - 0.5 * (F2)) * (1 - F2) * (F3)^(P - 3) *
    (P - 2) - 0.5 * (F10)) * (P - 1), MARGIN = 1, STATS = ewv1/ewu_h,
    FUN = "*")
  F12 <- sweep(F9 * (1 - FiMat) * (F3)^(P - 2)/(F1), MARGIN = 1,
    STATS = depsi * sigx7/ewv1_h, FUN = "*")
  F13 <- sweep((-0.5 * (F8)), MARGIN = 1, STATS = 0.5 * (ewvuepsi/ewv1_h),
    FUN = "+")
  F14 <- sweep((F13 * (1 - FiMat) * (F3)^(P - 2)/(F1)), MARGIN = 1,
    STATS = depsi/ewv1_h, FUN = "*")
  F15 <- sweep(((0.5 - 0.5 * (F2))^2 * (F3)^(P - 3) * (P -
    2) - 0.5 * F14), MARGIN = 1, STATS = ewv1/ewu_h, FUN = "*")
  F16 <- sweep((F15 - 0.5 * ((0.5 - 0.5 * (F2)) * (F3)^(P -
    2))) * (P - 1), MARGIN = 1, STATS = ewv1/ewu_h, FUN = "*")
  F17 <- sweep(((0.5 - 0.5 * (F2)) * (1 - F2) * F3^(P - 3) *
    (P - 2) - F14) * (P - 1), MARGIN = 1, STATS = ewv1/ewu_h,
    FUN = "*")
  F18 <- sweep(F9 * (1 - FiMat) * F3^(P - 2)/(F1), MARGIN = 1,
    STATS = depsi/ewv1_h, FUN = "*")
  F19 <- sweep(F13, MARGIN = 1, STATS = sigx7/ewv1_h, FUN = "*")
  F20 <- sweep((((F19 - 0.5) * F2 + 0.5) * F3^(P - 2) + (F5 +
    F7) * (0.5 - 0.5 * (F2)) * F3^(P - 3) * (P - 2)) * (P -
    1), MARGIN = 1, STATS = ewv1/ewu_h, FUN = "*")
  F21 <- sweep((0.5 - 0.5 * (F2)), MARGIN = 1, STATS = ewv1/ewu_h,
    FUN = "*")
  F22 <- sweep(F8, MARGIN = 1, STATS = ewvuepsi/ewv1_h, FUN = "-")
  F23 <- sweep(F22, MARGIN = 1, STATS = sigx7^2/ewv1_h, FUN = "*")
  F24 <- sweep(F23, MARGIN = 1, STATS = 0.5 * (ewv1/ewu_h),
    FUN = "+")
  F25 <- sweep(F24 * F2 + 0.5 * (F5 + 0.5 * F6), MARGIN = 1,
    STATS = ewv1/ewu_h, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  XF7 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * (F3)^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum)
    XF7[, k] <- apply(sweep(S * (F3^(P - 2) + F3^(P - 2) *
      log(F3) * (P - 1)) * (1 - F2), MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx1 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewuP * ewuv * ewz/sigx2,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * dvsi *
    (epsilon)/ewv2_h^3, FUN = "*")
  XF5 <- sweep(Xvar, MARGIN = 1, STATS = S * sigx14 * sDiv,
    FUN = "*")
  XF6 <- sweep(XF1, MARGIN = 1, STATS = sigx1, FUN = "*") +
    sweep(XF5 + XF1, MARGIN = 1, STATS = dmuv/ewv1_h, FUN = "*") -
    sweep(XF2, MARGIN = 1, STATS = 1/ewu_h, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  UF4 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(F4, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)
    UF4[, k] <- apply(sweep((F3^(P - 2) + F3^(P - 2) * log(F3) *
      (P - 1)) * F21, MARGIN = 1, STATS = uHvar[, k], FUN = "*"),
      1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx5 * sDiv, FUN = "*")
  UF3 <- sweep(UF1, MARGIN = 1, STATS = dmuv/ewv1_h, FUN = "*") -
    sweep(UF2, MARGIN = 1, STATS = 1/ewu_h, FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  VF4 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F5 + F7) * (F3)^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
    VF4[, k] <- apply(sweep((F5 + F7) * (F3^(P - 2) + F3^(P -
      2) * log(F3) * (P - 1)), MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmuv, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  VF3 <- sweep(VF1, MARGIN = 1, STATS = dmuv/ewv1_h, FUN = "*") -
    sweep(VF2, MARGIN = 1, STATS = 1/ewu_h, FUN = "*")
  HX1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in 1:Q) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      (((1 - F2)^2 * F3^(P - 3) * (P - 2) - F18) * (P -
        1))[, r] * pmuv * ewuP * ewuv * ewz/(sigx2 *
      sigx3), FUN = "*"), Xvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      S * F11[, r] * pmuv/sigx6 * ewuP * ewuv * ewz, FUN = "*"),
      uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      S * (((F5 + F7) * (1 - F2) * F3^(P - 3) * (P - 2) +
      F12) * (P - 1))[, r] * pmuv * ewuP * ewuv * ewz/sigx6,
      FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      F16[, r] * pmuv/sigx6 * ewuP * ewuv * ewz, FUN = "*"),
      uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      F20[, r] * pmuv/sigx6 * ewuP * ewuv * ewz, FUN = "*"),
      vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      ((F25 * F3^(P - 2) + (F5 + F7)^2 * F3^(P - 3) * (P -
        2)) * (P - 1))[, r] * pmuv/sigx6 * ewuP * ewuv *
      ewz, FUN = "*"), vHvar)
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1, ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) + crossprod(sweep(S *
    Xvar, MARGIN = 1, STATS = wHvar * ewuP * ewuv * ewz/(sigx2 *
    sigx3), FUN = "*"), XF6) + crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * prC * dvsi * sigx15/(ewv2_h^3 * sigx3),
    FUN = "*"), Xvar) - crossprod(sweep(XF3 + XF4, MARGIN = 1,
    STATS = wHvar/sigx3^2, FUN = "*"), XF3 + XF4)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- Reduce("+",
    HXU1) + crossprod(sweep(S * Xvar, MARGIN = 1, STATS = wHvar/sigx6 *
    ewuP * ewuv * ewz, FUN = "*"), UF3) + crossprod(sweep(XF1,
    MARGIN = 1, STATS = wHvar * sigx5/sigx6 * ewuP * ewuv *
      ewz, FUN = "*") + sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * sigx16 * sDiv/sigx6 * ewuP * ewuv * ewz, FUN = "*"),
    uHvar) - crossprod(sweep((XF3 + XF4), MARGIN = 1, STATS = wHvar *
    Q * wzdeno * gamma(P)/sigx6^2 * ewuP * ewuv * ewz, FUN = "*"),
    UF2)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- Reduce("+", HXV1) + crossprod(sweep(S *
    Xvar, MARGIN = 1, STATS = wHvar * ewuP * ewuv * ewz/sigx6,
    FUN = "*"), VF3) + crossprod(sweep(XF1, MARGIN = 1, STATS = wHvar *
    sigx9 * ewuP * ewuv * ewz/sigx6, FUN = "*") + sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * dmuv * sigx17 * sDiv/ewv1_h *
      ewuP * ewuv * ewz/sigx6, FUN = "*"), vHvar) - crossprod(sweep((XF3 +
    XF4), MARGIN = 1, STATS = wHvar * Q * wzdeno * gamma(P)/sigx6^2 *
    ewuP * ewuv * ewz, FUN = "*"), VF2)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * sigx18 * dvsi * (epsilon)/(sigx3 *
      ewv2_h^3) * prC, FUN = "*") - sweep((XF3 + XF4),
    MARGIN = 1, STATS = wHvar * sigx10 * ewv2_h/(sigx3 *
      ewv2_h)^2 * prC, FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + 2 * nvZVvar + 1] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * sigx13 * S * dmuv/ewv1_h *
      ewuP * ewuv * ewz/sigx3, FUN = "*") - sweep((XF3 +
    XF4), MARGIN = 1, STATS = wHvar * sigx13 * pmuv/sigx3 *
    ewuP * ewuv * ewz/sigx3, FUN = "*") + (sweep(XF7, MARGIN = 1,
    STATS = wHvar/sigx2 * pmuv * ewuP * ewuv * ewz/sigx3,
    FUN = "*") - sweep(XF1, MARGIN = 1, STATS = wHvar * 0.5 *
    (Wu)/sigx2 * pmuv * ewuP * ewuv * ewz/sigx3, FUN = "*") -
    (sweep(XF1, MARGIN = 1, STATS = wHvar * Q * wzdeno *
      digamma(P) * gamma(P)/sigx2^2 * pmuv * ewuP * ewuv *
      ewz/sigx3, FUN = "*") + sweep(Xvar, MARGIN = 1, STATS = wHvar *
      S * sigx13/ewu_h * pmuv * ewuP * ewuv * ewz/sigx3,
      FUN = "*"))))
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(XF2,
    MARGIN = 1, STATS = wHvar * sigx11 * ewuP * ewuv * ewz/sigx3,
    FUN = "*") - (sweep((XF3 + XF4), MARGIN = 1, STATS = wHvar *
    sigx12/sigx3 * ewz/sigx3, FUN = "*") + sweep(Xvar, MARGIN = 1,
    STATS = wHvar * prC * dvsi * (epsilon)/(wzdeno * ewv2_h^3) *
      ewz/sigx3, FUN = "*")), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- Reduce("+", HU1) + crossprod(sweep(UF1,
    MARGIN = 1, STATS = wHvar * sigx19/sigx6 * ewuP * ewuv *
      ewz, FUN = "*") + sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx22 * sDiv/sigx6 * ewuP * ewuv * ewz, FUN = "*") -
    sweep(UF2, MARGIN = 1, STATS = wHvar * sigx4/sigx6 *
      ewuP * ewuv * ewz, FUN = "*"), uHvar) - crossprod(sweep(UF2,
    MARGIN = 1, STATS = wHvar * ewuP * ewuv * ewz/sigx6^2 *
      ewuP * ewuv * ewz, FUN = "*"), UF2)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+", HUV1) +
    crossprod(uHvar, sweep(VF1, MARGIN = 1, STATS = wHvar *
      0.5 * (dmuv * ewv1_h/ewu_h)/sigx6 * ewuP * ewuv *
      ewz, FUN = "*") - sweep(VF2, MARGIN = 1, STATS = wHvar *
      sigx4/sigx6 * ewuP * ewuv * ewz, FUN = "*")) + crossprod(sweep(UF1,
    MARGIN = 1, STATS = wHvar * sigx9/sigx6 * ewuP * ewuv *
      ewz, FUN = "*") + sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx24 * sDiv/sigx6 * ewuP * ewuv * ewz, FUN = "*"),
    vHvar) - crossprod(sweep(UF2, MARGIN = 1, STATS = wHvar *
    ewuP * ewuv * ewz/sigx6^2 * ewuP * ewuv * ewz, FUN = "*"),
    VF2)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(UF2,
    MARGIN = 1, STATS = -wHvar * (sigx10 * prC * ewuP * ewuv *
      ewv2_h * ewz/sigx25), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + 2 *
    nvZVvar + 1] <- colSums(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    0.5 * (dmuv * ewv1_h/ewu_h) * sigx13 * ewuP * ewuv *
    ewz/sigx3, FUN = "*") - sweep(UF2, MARGIN = 1, STATS = wHvar *
    ewuP * ewuv * ewz * pmuv/sigx6 * sigx13 * ewuP * ewuv *
    ewz/sigx3, FUN = "*") + sweep(UF4, MARGIN = 1, STATS = wHvar/sigx2 *
    pmuv * ewuP * ewuv * ewz/sigx3, FUN = "*") - 0.5 * (sweep(UF1,
    MARGIN = 1, STATS = wHvar * Wu/sigx2 * pmuv * ewuP *
      ewuv * ewz/sigx3, FUN = "*") + sweep(uHvar, MARGIN = 1,
    STATS = wHvar * sDiv/sigx2 * pmuv * ewuP * ewuv * ewz/sigx3,
    FUN = "*")) - (sweep(UF1, MARGIN = 1, STATS = wHvar *
    Q * wzdeno * digamma(P) * gamma(P)/sigx2^2 * pmuv * ewuP *
    ewuv * ewz/sigx3, FUN = "*") + sweep(uHvar, MARGIN = 1,
    STATS = wHvar * sigx13 * sigx4 * pmuv * ewuP * ewuv *
      ewz/sigx3, FUN = "*")))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- crossprod(sweep(UF2, MARGIN = 1, STATS = wHvar *
    sigx26 * ewuP * ewuv * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HV1) + crossprod(sweep(VF1, MARGIN = 1, STATS = wHvar *
    (ewv1 * pmuv/(2 * ewu) - 2 * (sigx8 * dmuv))/sigx6 *
    ewuP * ewuv * ewz, FUN = "*") + sweep(VF2, MARGIN = 1,
    STATS = wHvar * ewv1/(2 * ewu)/sigx6 * ewuP * ewuv *
      ewz, FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (ewv1 * (pmuv - sigx8 * dmuv)/(2 * ewu) - (0.25 * (ewv1_h/ewu_h) +
      0.25 * (S * (epsilon)/ewv1_h) - sigx8^2 * ev1epsi) *
      dmuv) * sDiv/sigx6 * ewuP * ewuv * ewz, FUN = "*"),
    vHvar) - crossprod(sweep(VF2, MARGIN = 1, STATS = wHvar *
    ewuP * ewuv * ewz/sigx6^2 * ewuP * ewuv * ewz, FUN = "*"),
    VF2)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(VF2, MARGIN = 1, STATS = -wHvar *
    (sigx10 * prC * ewuP * ewuv * ewv2_h * ewz/sigx25), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + 2 * nvZVvar + 1] <- colSums(sweep(VF4,
    MARGIN = 1, STATS = wHvar/sigx2 * pmuv * ewuP * ewuv *
      ewz/sigx3, FUN = "*") - sweep(VF1, MARGIN = 1, STATS = wHvar *
    0.5 * (Wu)/sigx2 * pmuv * ewuP * ewuv * ewz/sigx3, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = wHvar * sigx13 * ewv1/(2 *
      ewu) * pmuv * ewuP * ewuv * ewz/sigx3, FUN = "*") -
    sweep(VF1, MARGIN = 1, STATS = wHvar * Q * wzdeno * digamma(P) *
      gamma(P)/sigx2^2 * pmuv * ewuP * ewuv * ewz/sigx3,
      FUN = "*") - (sweep(VF2, MARGIN = 1, STATS = wHvar *
    sigx13 * ewuP * ewuv * ewz * pmuv/sigx6 * ewuP * ewuv *
    ewz/sigx3, FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = wHvar *
    sigx13 * sigx8 * dmuv * ewuP * ewuv * ewz/sigx3, FUN = "*")))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(VF2,
    MARGIN = 1, STATS = wHvar * sigx26 * ewuP * ewuv * ewz/sigx3,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) -
      1) - 0.25) * dvsi * (epsilon)^2/(sigx3 * ewv2_h^3) -
      (sigx10 * prC + 0.5 * (sigx3 * ewv2_h)) * sigx10/(sigx3 *
        ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), nXvar + nuZUvar + 2 * nvZVvar + 1] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx13 * sigx10 * prC *
      ewuP * ewuv * ewz * pmuv/(sigx3^2 * ewv2_h)), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx12 * sigx10/sigx3 +
      0.5 * (S^2 * dvsi * (epsilon)^2/(wzdeno * ewv2_h^2)))/ewv2_h -
      0.5 * (wzdeno * dvsi * ewv2_h/(wzdeno * ewv2_h)^2)) *
      prC * ewz/sigx3), FUN = "*"), Zvar)
  hessll[nXvar + nuZUvar + 2 * nvZVvar + 1, nXvar + nuZUvar +
    2 * nvZVvar + 1] <- sum(wHvar * (((sF3lo - 0.5 * (Wu *
    sF3))/sigx2 - (sigx29 * sigx13 + Q * (sigx30 * digamma(P) +
    (digamma(P)^2 + trigamma(P)) * sDiv) * wzdeno * gamma(P)/sigx2^2)) *
    ewuP * ewuv * ewz * pmuv/sigx3))
  hessll[nXvar + nuZUvar + 2 * nvZVvar + 1, (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- colSums(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    sigx28 * ewuP * ewuv * ewz * pmuv/sigx3, FUN = "*"))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv2_h) +
      ewv2_h/(wzdeno * ewv2_h)^2) * dvsi - (sigx12^2/sigx3 +
      Q * (2 - 2 * (Q^2 * wzdeno * ewz * gamma(P)^2/sigx2^2)) *
        ewuP * ewuv * gamma(P) * pmuv * sDiv/sigx2^2)) *
      ewz + sigx11 * ewuP * ewuv * pmuv * sDiv - prC *
      dvsi/(wzdeno * ewv2_h)) * ewz/sigx3, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf gamma-normal distribution
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
zisfgammanormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfgammanormlike_logit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(czisfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfgammanormlike_logit,
      grad = cgradzisfgammanormlike_logit, hess = chesszisfgammanormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = function(parm) as(-chesszisfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hess = function(parm) -chesszisfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgammanormlike_logit(mleObj$par,
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
      mleObj$hessian <- chesszisfgammanormlike_logit(mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfgammanormlike_logit(mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfgammanormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgammanormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## cauchit specification class membership
zisfgammanormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfgammanormlike_cauchit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(czisfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfgammanormlike_cauchit,
      grad = cgradzisfgammanormlike_cauchit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradzisfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgammanormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- czisfgammanormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgammanormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## probit specification class membership
zisfgammanormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfgammanormlike_probit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(czisfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfgammanormlike_probit,
      grad = cgradzisfgammanormlike_probit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradzisfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgammanormlike_probit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- czisfgammanormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgammanormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## cloglog specification class membership
zisfgammanormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfgammanormlike_cloglog(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(czisfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfgammanormlike_cloglog,
      grad = cgradzisfgammanormlike_cloglog, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradzisfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgammanormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradzisfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- czisfgammanormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgammanormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

# Different sigma_v

## logit specification class membership
mnsfgammanormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfgammanormlike_logit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfgammanormlike_logit,
      grad = cgradmnsfgammanormlike_logit, hess = chessmnsfgammanormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = function(parm) as(-chessmnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hess = function(parm) -chessmnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfgammanormlike_logit(mleObj$par,
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
      mleObj$hessian <- chessmnsfgammanormlike_logit(mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfgammanormlike_logit(mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfgammanormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfgammanormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## cauchit specification class membership
mnsfgammanormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfgammanormlike_cauchit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfgammanormlike_cauchit,
      grad = cgradmnsfgammanormlike_cauchit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfgammanormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmnsfgammanormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfgammanormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## probit specification class membership
mnsfgammanormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfgammanormlike_probit(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfgammanormlike_probit,
      grad = cgradmnsfgammanormlike_probit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfgammanormlike_probit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmnsfgammanormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfgammanormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## cloglog specification class membership
mnsfgammanormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfgammanormlike_cloglog(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfgammanormlike_cloglog,
      grad = cgradmnsfgammanormlike_cloglog, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfgammanormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmnsfgammanormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfgammanormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

# Conditional efficiencies estimation ----------
#' efficiencies for zisf gamma-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v

## logit specification class membership
czisfgammanormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1/Hi2
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) -
      exp(Wv)
    mui_Ki <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) +
      exp(Wv)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv)/exp(Wu/2) + object$S * epsilon +
      exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) - object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)) * Gi/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv)/exp(Wu/2) - object$S *
      epsilon + exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) -
      object$S * epsilon/exp(Wv/2) + exp(Wv/2)) * Ki/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
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
czisfgammanormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1/Hi2
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) -
      exp(Wv)
    mui_Ki <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) +
      exp(Wv)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv)/exp(Wu/2) + object$S * epsilon +
      exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) - object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)) * Gi/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv)/exp(Wu/2) - object$S *
      epsilon + exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) -
      object$S * epsilon/exp(Wv/2) + exp(Wv/2)) * Ki/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
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
czisfgammanormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1/Hi2
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) -
      exp(Wv)
    mui_Ki <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) +
      exp(Wv)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv)/exp(Wu/2) + object$S * epsilon +
      exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) - object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)) * Gi/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv)/exp(Wu/2) - object$S *
      epsilon + exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) -
      object$S * epsilon/exp(Wv/2) + exp(Wv/2)) * Ki/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
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
czisfgammanormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1/Hi2
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) -
      exp(Wv)
    mui_Ki <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) +
      exp(Wv)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv)/exp(Wu/2) + object$S * epsilon +
      exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) - object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)) * Gi/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv)/exp(Wu/2) - object$S *
      epsilon + exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) -
      object$S * epsilon/exp(Wv/2) + exp(Wv/2)) * Ki/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
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
cmnsfgammanormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1/Hi2
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) -
      exp(Wv1)
    mui_Ki <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) +
      exp(Wv1)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2)
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
cmnsfgammanormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1/Hi2
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) -
      exp(Wv1)
    mui_Ki <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) +
      exp(Wv1)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2)
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
cmnsfgammanormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1/Hi2
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) -
      exp(Wv1)
    mui_Ki <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) +
      exp(Wv1)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2)
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
cmnsfgammanormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1/Hi2
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) -
      exp(Wv1)
    mui_Ki <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) +
      exp(Wv1)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2)
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
#' marginal impact on efficiencies for zisf gamma-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v

## logit specification class membership
czisfmarggammanorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarggammanorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
czisfmarggammanorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarggammanorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## probit specification class membership
czisfmarggammanorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarggammanorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
czisfmarggammanorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarggammanorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
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
cmnsfmarggammanorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarggammanorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmnsfmarggammanorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarggammanorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmnsfmarggammanorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarggammanorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmnsfmarggammanorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarggammanorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

